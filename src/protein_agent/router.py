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
    CHROMOPHORE,
    COMPUTE_PROPERTIES,
    FETCH_ENTRY,
    FETCH_FASTA,
    PLOT_HYDROPATHY,
    SEARCH_UNIPROT,
    UNIPROT_TOOLS,
)

PLOT_HINTS = re.compile(
    r"\b(plot|graph|figure|chart|hydropathy|hydrophobicity|gravy|transmembrane|"
    r"tm helix|topolog\w*|profile)\w*\b",
    re.IGNORECASE,
)
PDB_ID_RE = re.compile(r"\b[0-9][A-Za-z0-9]{3}\b")
ACCESSION_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:-[0-9]{1,2})?\b"
)
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
    extra_tool_names: tuple[str, ...] = ()

    @property
    def tools(self) -> list[dict[str, Any]]:
        from .tools import SEQUENCE_TOOLS

        wanted = set(self.tool_names) | set(self.extra_tool_names)
        return [
            t
            for t in (*UNIPROT_TOOLS, *SEQUENCE_TOOLS)
            if t["function"]["name"] in wanted
        ]


def build_plan(question: str, *, pdbe_tools: list[dict[str, Any]] | None = None) -> Plan:
    """Return the tools worth advertising for this question."""
    names: set[str] = set()
    services: list[str] = []
    reasons: list[str] = []

    is_sequence_query = bool(re.search(r"[ACDEFGHIKLMNPQRSTVWY]{20,}", question.upper()))
    if is_sequence_query and not ACCESSION_RE.search(question):
        # The subject is the sequence itself. It must not be attached to an accession, so
        # the entry/fasta tools are withheld rather than merely discouraged.
        names.update({COMPUTE_PROPERTIES, CHROMOPHORE})
        services.append("local-compute")
        reasons.append(
            "sequence pasted: identity is unresolved from sequence alone, so answer only "
            "from computed properties; UniProt text search cannot match a sequence"
        )
        if PLOT_HINTS.search(question):
            names.add(PLOT_HYDROPATHY)
            reasons.append(
                "plot requested: only plot_hydropathy may describe a curve, so it is "
                "advertised and an answer about a plot without calling it is blocked"
            )
        return Plan(tuple(services), tuple(sorted(names)), "; ".join(reasons))

    names.add(FETCH_ENTRY)  # the minimum: no claim about a named protein without an entry
    services.append("uniprot")
    reasons.append("resolve and verify the entry before summarising")
    if ACCESSION_RE.search(question):
        names.add(FETCH_ENTRY)
        reasons.append("explicit accession supplied")
    elif is_sequence_query:
        names.update({COMPUTE_PROPERTIES, CHROMOPHORE})
        reasons.append(
            "sequence pasted: answer from computed properties; UniProt text search "
            "cannot identify a protein from its sequence, so it is not advertised"
        )
    elif not PDB_ID_RE.search(question):
        names.add(SEARCH_UNIPROT)
        reasons.append("name/organism given, so identity must be resolved first")

    if SEQUENCE_HINTS.search(question):
        names.add(FETCH_FASTA)
        reasons.append("sequence-derived detail requested")

    wants_structure = bool(STRUCTURE_HINTS.search(question)) or bool(PDB_ID_RE.search(question))
    direct_pdb_only = bool(PDB_ID_RE.search(question)) and not bool(ACCESSION_RE.search(question))
    if direct_pdb_only:
        names.difference_update({FETCH_ENTRY, FETCH_FASTA, SEARCH_UNIPROT})
        services = [service for service in services if service != "uniprot"]
        reasons.append(
            "explicit PDB ID supplied: query PDBe directly without inventing a UniProt entry"
        )
    if wants_structure:
        services.append("pdbe")
        if not direct_pdb_only and not ACCESSION_RE.search(question):
            names.add(SEARCH_UNIPROT)  # resolve a named protein before structure search
        reasons.append("experimental structure evidence requested")
        # PDBe tools are discovered from the live MCP server; never guessed here.
        if pdbe_tools is None:
            reasons.append("PDBe MCP tools must be discovered before structural retrieval")
        else:
            reasons.append("PDBe API/Search MCP tools were discovered live")

    if ANNOTATION_HINTS.search(question):
        reasons.append("annotation fields must be read from the entry payload")

    return Plan(
        tuple(services),
        tuple(sorted(names)),
        "; ".join(reasons),
        tuple(pdbe_tool_names(pdbe_tools)) if wants_structure and pdbe_tools else (),
    )


def pdbe_tool_names(pdbe_tools: list[dict[str, Any]] | None) -> list[str]:
    """Names advertised by a connected PDBe MCP server, if any."""
    out: list[str] = []
    for tool in pdbe_tools or []:
        name = (tool.get("function") or {}).get("name")
        if name:
            out.append(str(name))
    return out


__all__ = ["Plan", "build_plan", "pdbe_tool_names"]
