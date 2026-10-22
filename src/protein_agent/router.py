"""Which retrieval capabilities a question needs — decided by the model, not hardcoded.

This module only chooses *tool schemas to advertise*. The model still chooses which to
call, with what arguments, and how many times; nothing here fetches data or writes an
answer. Pre-selecting UniProt vs PDBe keeps prompts small and latency bounded without
ever letting the model answer from memory.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from .tools import (
    FETCH_ENTRY,
    FETCH_FASTA,
    PDBE_SEARCH,
    SEARCH_UNIPROT,
    UNIPROT_TOOLS,
)

PDB_ID_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")
ACCESSION_RE = re.compile(r"\b[A-Z0-9]{5,6}(?:-\d{1,2})?\b")
STRUCTURE_HINTS = re.compile(
    r"\b(pdb|structure|structural|crystal|crystallograph|nmr|cryo-?em|resolution|"
    r"chain|assembly|assemblies|ligand|complex|domain structure|conformer)\w*\b",
    re.IGNORECASE,
)
SEQUENCE_HINTS = re.compile(
    r"\b(sequence|residue|amino acid|length|motif|isoform|variant|mutation|peptide|"
    r"mass|mw|start codon)\w*\b",
    re.IGNORECASE,
)
ANNOTATION_HINTS = re.compile(
    r"\b(function|locat\w+|subcellular|pathway|disease|interact\w+|annotation|activity|"
    r"ptm|modification|binding site|go term|taxonomy|organism|gene name|reviewed)\w*\b",
    re.IGNORECASE,
)


@dataclass(slots=True)
class Plan:
    services: tuple[str, ...]
    tool_names: tuple[str, ...]
    reason: str

    @property
    def tools(self) -> list[dict[str, Any]]:
        wanted = set(self.tool_names)
        return [t for t in UNIPROT_TOOLS if t["function"]["name"] in wanted]


def build_plan(question: str, *, pdbe_tools: list[dict[str, Any]] | None = None) -> Plan:
    """Return the tools worth advertising for this question."""
    names = {FETCH_ENTRY}  # the minimum: no claim about a named protein without an entry
    services = ["uniprot"]
    reasons = ["resolve and verify the entry before summarising"]

    if ACCESSION_RE.search(question):
        names.add(FETCH_ENTRY)
        reasons.append("explicit accession supplied")
    elif not PDB_ID_RE.search(question):
        names.add(SEARCH_UNIPROT)
        reasons.append("name/organism given, so identity must be resolved first")

    if SEQUENCE_HINTS.search(question):
        names.add(FETCH_FASTA)
        reasons.append("sequence-derived detail requested")

    wants_structure = bool(STRUCTURE_HINTS.search(question)) or bool(PDB_ID_RE.search(question))
    if wants_structure:
        services.append("pdbe")
        names.add(SEARCH_UNIPROT)  # an unambiguous accession strengthens structure search
        reasons.append("experimental structure evidence requested")
        # PDBe tools are discovered from the live MCP server; never guessed here.
        if pdbe_tools is None:
            reasons.append(
                "PDBe MCP not configured on this node: report structural questions as "
                "unverified rather than answering them from knowledge"
            )

    if ANNOTATION_HINTS.search(question):
        reasons.append("annotation fields must be read from the entry payload")

    return Plan(tuple(services), tuple(sorted(names)), "; ".join(reasons))


def pdbe_tool_names(pdbe_tools: list[dict[str, Any]] | None) -> list[str]:
    """Names advertised by a connected PDBe MCP server, if any."""
    out: list[str] = []
    for tool in pdbe_tools or []:
        name = (tool.get("function") or {}).get("name")
        if name:
            out.append(str(name))
    return out or [PDBE_SEARCH]


__all__ = ["Plan", "build_plan", "pdbe_tool_names"]
