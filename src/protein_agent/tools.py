"""Tool schemas advertised to the model, plus name constants.

The model may only call these; it has no other way to learn a fact. Names match the
"Expose ... as tools" bullets in `uniprot/AGENTS.md` and the MCP-discovered names in
`pdbe-kb/AGENTS.md`.
"""

from __future__ import annotations

from typing import Any

SEARCH_UNIPROT = "search_uniprotkb"
FETCH_ENTRY = "fetch_uniprot_entry"
FETCH_FASTA = "fetch_fasta"
COMPUTE_PROPERTIES = "compute_sequence_properties"
CHROMOPHORE = "find_chromophore_candidates"
PDBE_API = "pdbe_api"  # placeholder: replace with names from the connected MCP server
PDBE_SEARCH = "run_pdbe_search_query"
PDBE_SCHEMA = "get_pdbe_search_schema"

_STRING = {"type": "string"}


def _fn(
    name: str,
    description: str,
    properties: dict[str, Any],
    required: list[str],
) -> dict[str, Any]:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": required,
                "additionalProperties": False,
            },
        },
    }


UNIPROT_TOOLS: list[dict[str, Any]] = [
    _fn(
        SEARCH_UNIPROT,
        "Search UniProtKB by name, gene, organism, or accession and return candidates. "
        "REQUIRED before choosing an accession from a gene or protein name. Never answer "
        "from remembered accessions. Cannot identify a protein from its sequence.",
        {
            "query": {
                **_STRING,
                "description": (
                    "UniProtKB query, e.g. 'gene:TP53 AND organism_id:9606 AND "
                    "reviewed:true'. Use a name, gene, accession, or organism. Do NOT "
                    "paste an amino-acid sequence here: full-text search does not match "
                    "sequences and will return nothing. Sequence-based identification "
                    "needs a similarity search, which this tool is not."
                ),
            },
            "fields": {
                "type": "array",
                "items": _STRING,
                "description": "Fields to return, e.g. accession, protein_name, "
                "gene_names, organism_name, length, reviewed",
            },
            "size": {"type": "integer", "minimum": 1, "maximum": 25},
        },
        ["query"],
    ),
    _fn(
        FETCH_ENTRY,
        "Fetch one full UniProtKB entry as JSON by primary accession. REQUIRED for any "
        "claim about function, location, domains, PTMs, disease, pathway, or "
        "cross-references.",
        {"accession": {**_STRING, "description": "Primary accession, e.g. P04637"}},
        ["accession"],
    ),
    _fn(
        FETCH_FASTA,
        "Fetch the canonical FASTA sequence for an accession. REQUIRED before stating a "
        "sequence length, residue position, or sequence-derived detail.",
        {"accession": _STRING},
        ["accession"],
    ),
]

SEQUENCE_TOOLS: list[dict[str, Any]] = [
    _fn(
        COMPUTE_PROPERTIES,
        "Compute biophysicochemical properties from an amino-acid sequence that the user "
        "pasted or uploaded (FASTA). REQUIRED before reporting composition, MW, pI, charge, "
        "extinction coefficient, GRAVY, aromaticity, instability index, or a hydropathy "
        "profile. These values must never be estimated from the sequence by eye or recalled.",
        {
            "sequence": {
                **_STRING,
                "description": "Amino-acid sequence (FASTA header line is tolerated; gaps "
                "and non-standard residues are reported as warnings, not silently removed)",
            },
            "label": {**_STRING, "description": "Optional name for the record"},
            "hydropathy_window": {"type": "integer", "minimum": 3, "maximum": 41},
        },
        ["sequence"],
    ),
    _fn(
        CHROMOPHORE,
        "Nominate candidate chromophore-forming tripeptide cores in a sequence. Returns "
        "CANDIDATES ONLY — it does not establish photon absorption or that the protein is "
        "photo-active. Always pair with fetch_uniprot_entry to read the annotated "
        "Modified residue features and absorption maxima before saying anything about a "
        "chromophore.",
        {"sequence": _STRING},
        ["sequence"],
    ),
]

# PDBe tools come from the connected MCP server's own schemas. Do not hard-code guessed
# names: load them with pdbe_mcp.load_pdbe_tools() and fall back to these placeholders.
PDBE_TOOL_NOTE = (
    "Register the PDBe API and Search server tools discovered from the MCP connection; "
    "see pdbe-kb/AGENTS.md."
)


def all_tools(
    *,
    include_pdbe: bool = True,
    pdbe_tools: list[dict[str, Any]] | None = None,
    include_sequence: bool = True,
) -> list[dict[str, Any]]:
    tools = list(UNIPROT_TOOLS)
    if include_sequence:
        tools.extend(SEQUENCE_TOOLS)
    if include_pdbe:
        tools.extend(pdbe_tools or [])
    return tools


__all__ = [
    "FETCH_ENTRY",
    "FETCH_FASTA",
    "PDBE_API",
    "PDBE_SCHEMA",
    "PDBE_SEARCH",
    "PDBE_TOOL_NOTE",
    "SEARCH_UNIPROT",
    "SEQUENCE_TOOLS",
    "UNIPROT_TOOLS",
    "all_tools",
]
