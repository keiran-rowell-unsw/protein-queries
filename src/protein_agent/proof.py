"""The no-memory-only-answers enforcement layer.

Two mechanisms, both required:

1. `tool_choice="required"` on every model request — the model cannot reply with prose
   until it has asked for data (verified against the live gateway).
2. This module — after the loop ends, check that the answer's verifiable tokens actually
   appear in recorded tool results, and that enough tool calls ran. An answer that cites
   an identifier no tool returned is blocked rather than shown.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

# UniProt canonical accession grammar, documented as: 8 or 10 characters, first is A,
# O or Q, second is a digit, and "-" separates isoform suffixes. The previous naive
# [A-Z0-9]{5,6} flagged ordinary uppercase words (BLAST, KKKDDD) as accessions.
UNIPROT_ACCESSION = re.compile(
    r"\b(?:UPI-[0-9A-F]{8}|[ABOPQR][0-9][A-Z0-9]{3,4}|[ABOPQR][0-9][A-Z0-9]{6,8})"
    r"(?:-[0-9]{1,2})?\b"
)
PDB_ID = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")
GO_TERM = re.compile(r"\bGO:[0-9]{7}\b")
# Only words that claim an actual rendered IMAGE. Deliberately excludes "profile" and
# "curve": a hydropathy profile is legitimate tabular data returned by
# compute_sequence_properties, and treating the word as a plot claim blocks correct answers.
PLOT_LANGUAGE = re.compile(
    r"\b(?:the|this|a|your)?\s*(?:plot|figure|graph|chart)\b|"
    r"\b(?:plotted|rendered|depicted|as shown)\b|"
    r"\b(?:x-axis|y-axis|axes)\b",
    re.IGNORECASE,
)
RESIDUE_POS = re.compile(r"\b(?:res(?:idue)?|aa)\s*#?\s*[0-9]+", re.IGNORECASE)

# Identifiers that look factual but may legitimately appear without a tool call.
GENE_NAMES = re.compile(r"\b[A-Z0-9]{2,10}\b")


class AnswerProofError(RuntimeError):
    """Raised when an answer cannot be tied back to tool evidence."""


@dataclass(slots=True)
class Evidence:
    """Everything the tools returned during one run."""

    calls: list[dict[str, Any]] = field(default_factory=list)  # name + arguments
    raw_results: list[str] = field(default_factory=list)  # tool message payloads

    @property
    def call_count(self) -> int:
        return len(self.calls)

    def tool_names(self) -> list[str]:
        return [str(c.get("name")) for c in self.calls]

    def used_tool(self, name: str) -> bool:
        return any(str(c.get("name")).endswith(name) for c in self.calls)

    def subjects(self) -> list[str]:
        """Identifiers the tools were actually asked about, in call order."""
        out: list[str] = []
        for call in self.calls:
            args = call.get("arguments") or {}
            for key in ("accession", "pdb_id", "query", "sequence", "label"):
                value = args.get(key)
                if isinstance(value, str) and value.strip():
                    out.append(value.strip())
                    break
        return out

    def corpus(self) -> str:
        return "\n".join(self.raw_results)

    def add_call(self, name: str, arguments: dict[str, Any]) -> None:
        self.calls.append({"name": name, "arguments": arguments})

    def add_result(self, payload: str) -> None:
        self.raw_results.append(payload)


def extract_claims(text: str) -> dict[str, set[str]]:
    """Verifiable tokens pulled out of a draft answer."""
    return {
        "accessions": set(m.upper() for m in UNIPROT_ACCESSION.findall(text)),
        "pdb_ids": set(m.upper() for m in PDB_ID.findall(text)),
        "go_terms": set(GO_TERM.findall(text)),
    }


def subject_mismatch(
    question: str, evidence: Evidence, *, ignore_accessions: frozenset[str] = frozenset()
) -> list[str]:
    """Detect the specific failure where the model retrieves *something* to satisfy the
    tool requirement but changes the subject.

    A sequence-pasted question is about the sequence the user supplied. If the run then
    fetched an accession that the user never named and the assistant never verified,
    that is a substitution: the answer describes a different protein than was asked about.
    """
    problems: list[str] = []
    asked = question.strip().upper()
    supplied_sequence = bool(re.search(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", asked.replace(" ", "")))
    if not supplied_sequence:
        return problems
    named = set(UNIPROT_ACCESSION.findall(question))
    for subject in evidence.subjects():
        for candidate in UNIPROT_ACCESSION.findall(subject):
            if candidate in named or candidate in ignore_accessions:
                continue
            problems.append(
                f"{candidate} was fetched, but the user asked about a pasted sequence and "
                f"never named {candidate}; results describe a different protein"
            )
    return problems


def unsupported_claims(
    answer: str,
    evidence: Evidence,
    *,
    min_calls: int = 1,
    question: str = "",
    question_subject_check: list[str] | None = None,
) -> list[str]:
    """Return human-readable reasons the answer is not backed by tool output."""
    problems: list[str] = []
    if evidence.call_count < min_calls:
        problems.append(
            f"no tool call was made ({evidence.call_count} of >= {min_calls} required); "
            "the answer would be model memory"
        )
    if not answer.strip():
        problems.append("the model returned no content")
        return problems

    problems.extend(question_subject_check or [])
    corpus = evidence.corpus().upper()
    claims = extract_claims(answer)
    for kind, tokens in claims.items():
        for token in sorted(tokens):
            if token.upper() in corpus:
                continue
            # A run of residues quoted from the user's own pasted sequence is not a
            # database identifier. Only treat it as a claim about the record when it is
            # absent from the question as well as the tool output.
            if kind == "accessions" and token.upper().replace("-01", "") in question.upper():
                continue
            problems.append(
                f"{kind[:-1]} {token} appears in the answer but was never returned by a tool"
            )
    # Sequence/residue claims need a sequence-bearing call: fetch_fasta, or an entry
    # payload carrying the sequence/length that the claim depends on.
    if RESIDUE_POS.search(answer) and not any(
        "sequence" in r.lower() or "fasta" in r.lower() for r in evidence.raw_results
    ):
        problems.append("residue-level claim with no sequence retrieval in this run")
    if RESIDUE_POS.search(answer) and not (
        evidence.used_tool("fetch_fasta") or evidence.used_tool("fetch_uniprot_entry")
    ):
        problems.append(
            "residue-level claim without fetch_fasta or fetch_uniprot_entry in this run"
        )
    # A figure may only be described by the tool that drew it. Without this, a model that
    # has the underlying numbers will happily "describe" a plot it never rendered.
    if PLOT_LANGUAGE.search(answer) and not evidence.used_tool("plot_hydropathy"):
        problems.append(
            "describes a plot/figure but plot_hydropathy was never called in this run"
        )
    if evidence.used_tool("plot_hydropathy"):
        rendered = [
            r for r in evidence.raw_results if '"image"' in r or "hydropathy" in r.lower()
        ]
        if not rendered:
            problems.append("plot_hydropathy result is not on the evidence record")
    return problems


def prove_answer(
    answer: str,
    evidence: Evidence,
    *,
    question: str = "",
    min_calls: int = 1,
    on_unproven: str = "block",
    trusted_accessions: frozenset[str] = frozenset(),
) -> str:
    """Return the answer if it is grounded; otherwise raise or replace it."""
    mismatch = subject_mismatch(question, evidence, ignore_accessions=trusted_accessions)
    problems = unsupported_claims(
        answer,
        evidence,
        min_calls=min_calls,
        question=question,
        question_subject_check=mismatch,
    )
    if not problems:
        return answer
    if on_unproven == "warn":
        return answer + "\n\n> Unverified claims flagged: " + "; ".join(problems)
    raise AnswerProofError(
        "Answer blocked: not every claim is backed by a tool result.\n- "
        + "\n- ".join(problems)
        + "\nRe-run with retrieval, or report the gap to the user instead of answering."
    )


__all__ = [
    "AnswerProofError",
    "Evidence",
    "extract_claims",
    "prove_answer",
    "unsupported_claims",
]
