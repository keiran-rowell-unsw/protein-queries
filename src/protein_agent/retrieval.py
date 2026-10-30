"""Retrieval clients: UniProt REST (implemented) and PDBe MCP (stubs).

These are the only sanctioned sources of protein facts. Every result carries `_provenance`
so the proof layer can hold the answer to it.
"""

from __future__ import annotations

import json
import re
from datetime import UTC, datetime
from typing import Any

import httpx

from .config import Settings
from .pdbe_mcp import PdbeMcpError
from .tools import (
    CHROMOPHORE,
    COMPUTE_PROPERTIES,
    FETCH_ENTRY,
    FETCH_FASTA,
    PLOT_HYDROPATHY,
    SEARCH_UNIPROT,
)

PROTEIN_MIN_RUN = 8
# Isoform suffixes are optional and 1-2 digits. Written as a single alternation inside one
# group on purpose: the previous `^[A-Z0-9]{6,10}$|^[A-Z0-9]{5,10}-\d{2}$` was ordered and
# greedy, so it silently rejected the common single-digit case (P04637-1) while accepting
# P04637-12.
ACCESSION_RE = re.compile(r"(?:[A-Z0-9]{6,10}|[A-Z0-9]{5,10}-[0-9]{1,2})\Z")
DEFAULT_FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "length",
    "reviewed",
    "entry_type",
]


PTM_FEATURES = ("Modified residue", "Glycosylation", "Signal", "Disulfide bond", "Cross-link")


def extract_photochemistry(entry: dict[str, Any]) -> dict[str, Any]:
    """Pull absorption and PTM/Processing facts out of a fetched UniProt entry.

    Returns only what the entry states. An absent block is reported as *not annotated*
    for this entry — never as *no absorption* or *not photo-active*.
    """
    out: dict[str, Any] = {"absorption": None, "ptm_processing": [], "signal_peptide": None}
    for comment in entry.get("comments", []) or []:
        if comment.get("commentType") == "BIOPHYSICOCHEMICAL PROPERTIES":
            absb = comment.get("absorption")
            if absb:
                out["absorption"] = {
                    "max_nm": absb.get("max"),
                    "approximate": absb.get("approximate"),
                    "note": _texts(absb.get("note")),
                    "evidence": _evidence(absb),
                }
            for key, label in (
                ("thermodynamicStability", "thermodynamic_stability"),
                ("kineticStability", "kinetic_stability"),
            ):
                if comment.get(key):
                    out[label] = comment[key]
        if comment.get("commentType") == "PTM" and comment.get("similarity") is None:
            out.setdefault("ptm_comments", []).append(
                {"texts": _texts(comment), "evidence": _evidence(comment)}
            )
    for feature in entry.get("features", []) or []:
        ftype = feature.get("type")
        if ftype in PTM_FEATURES or ftype == "Signal":
            location = feature.get("location", {}) or {}
            out["ptm_processing"].append(
                {
                    "type": ftype,
                    "start": (location.get("start") or {}).get("value"),
                    "end": (location.get("end") or {}).get("value"),
                    "description": feature.get("description"),
                    "evidence": _evidence(feature),
                }
            )
            if ftype == "Signal":
                out["signal_peptide"] = out["ptm_processing"][-1]
    if out["absorption"] is None:
        out["absorption_note"] = (
            "no absorption data annotated in this entry; this is not evidence that the "
            "protein fails to absorb photons"
        )
    return out


def _texts(node: Any) -> list[str]:
    if not isinstance(node, dict):
        return []
    collected: list[str] = []
    for entry in node.get("texts", []) or []:
        if isinstance(entry, dict) and entry.get("value"):
            collected.append(entry["value"])
    return collected


def _evidence(node: Any) -> list[dict[str, Any]]:
    items = node.get("evidences", []) if isinstance(node, dict) else []
    return [
        {"code": e.get("evidenceCode"), "source": e.get("source"), "id": e.get("id")}
        for e in items or []
        if isinstance(e, dict)
    ]


class RetrievalError(RuntimeError):
    """A retrieval call failed; report it as a limitation, never paper over it."""


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _provenance(source: str, endpoint: str, subject: str, **extra: Any) -> dict[str, Any]:
    record = {
        "source": source,
        "endpoint": endpoint,
        "subject": subject,
        "retrieved_at": _now(),
    }
    record.update(extra)
    return record


PROTEIN_LIKE_RE = re.compile(
    "^[ACDEFGHIKLMNPQRSTVWY]{" + str(PROTEIN_MIN_RUN) + ",}$"
)


def looks_like_amino_acid_sequence(query: str) -> bool:
    """True when a search string is a raw peptide, which UniProtKB search cannot match.

    Established empirically on this node: searching the exact sequence of P42212 returns
    0 results (total null), quoted or not. Guarding here stops the agent from burning its
    tool budget re-issuing a query that is structurally incapable of succeeding.
    """
    candidate = query.strip().strip('"').upper().replace(" ", "")
    # Any run of >=8 unbroken residues is treated as sequence input, so a model that
    # shortens or fragments the query to dodge the guard still cannot reach the endpoint.
    if not PROTEIN_LIKE_RE.search(candidate):
        return False
    if not set(candidate) <= set("ACDEFGHIKLMNPQRSTVWY"):
        return False
    # a real query carries operators or organism hints; a pasted sequence does not
    return not any(tok in candidate for tok in ("AND", "OR", ":", "GENE:", "ORGANISM"))


class UniProtClient:
    def __init__(self, settings: Settings, *, client: httpx.Client | None = None) -> None:
        self.settings = settings
        self._owns = client is None
        self._client = client or httpx.Client(
            base_url=settings.uniprot_base_url,
            timeout=30.0,
            headers={"Accept": "application/json"},
        )

    def close(self) -> None:
        if self._owns:
            self._client.close()

    def search(
        self, query: str, *, fields: list[str] | None = None, size: int = 10
    ) -> dict[str, Any]:
        if looks_like_amino_acid_sequence(query):
            return {
                "results": [],
                "query_rejected": "amino_acid_sequence",
                "_provenance": _provenance(
                    "uniprot", "GET /uniprotkb/search", query, status="not_supported"
                ),
                "explanation": (
                    "UniProtKB text search does not match raw sequences; this query "
                    "cannot return a hit and was not sent. Identifying a protein from its "
                    "sequence requires a similarity search (BLAST), which is unavailable "
                    "on this node."
                ),
                "do_not": (
                    "Do not retry with this sequence, a fragment of it, reordered fields, "
                    "or a different size — none of those can succeed. Report identity as "
                    "unresolved and answer from the computed properties on record."
                ),
            }
        params = {
            "query": query,  # httpx URL-encodes; never concatenate raw user text
            "format": "json",
            "size": max(1, min(int(size), 25)),
            "fields": ",".join(fields or DEFAULT_FIELDS),
        }
        response = self._client.get("/uniprotkb/search", params=params)
        if response.status_code == 404 or not response.text.strip():
            return {
                "results": [],
                "_provenance": _provenance(
                    "uniprot", "GET /uniprotkb/search", query, status="no_match"
                ),
            }
        response.raise_for_status()
        payload = response.json()
        return {
            "results": payload.get("results", []),
            "_provenance": _provenance(
                "uniprot",
                "GET /uniprotkb/search",
                query,
                total=payload.get("totalResults"),
                next_link=response.headers.get("link"),
            ),
        }

    def entry(self, accession: str) -> dict[str, Any]:
        acc = accession.strip().upper()
        if not ACCESSION_RE.match(acc):
            raise RetrievalError(f"{accession!r} is not a well-formed UniProt accession.")
        response = self._client.get(f"/uniprotkb/{acc}.json")
        if response.status_code in {301, 302, 303, 307, 308}:
            raise RetrievalError(
                f"UniProt accession {acc} is inactive or redirected; "
                "verify its current primary accession explicitly."
            )
        if response.status_code == 404:
            raise RetrievalError(
                f"UniProt has no entry {acc} (404). Treat the identity as unresolved; do not guess."
            )
        response.raise_for_status()
        entry = response.json()
        sequence = entry.get("sequence", {}) or {}
        organism = entry.get("organism", {}) or {}
        return {
            # A curated, bounded full-entry view: enough for identity, function, sequence,
            # features, evidence and common links without feeding the model the raw record.
            "primaryAccession": entry.get("primaryAccession", acc),
            "secondaryAccessions": entry.get("secondaryAccessions", []),
            "uniProtkbId": entry.get("uniProtkbId"),
            "entryType": entry.get("entryType"),
            "annotationScore": entry.get("annotationScore"),
            "proteinExistence": entry.get("proteinExistence"),
            "organism": organism.get("scientificName"),
            "organismDetails": {
                "scientificName": organism.get("scientificName"),
                "commonName": organism.get("commonName"),
                "taxonId": organism.get("taxonId"),
                "lineage": organism.get("lineage", []),
            },
            "proteinDescription": entry.get("proteinDescription"),
            "genes": entry.get("genes", []),
            "sequence": {
                "value": sequence.get("value", ""),
                "length": sequence.get("length", len(sequence.get("value", ""))),
                "molWeight": sequence.get("molWeight"),
                "crc64": sequence.get("crc64"),
                "md5": sequence.get("md5"),
            },
            "comments": entry.get("comments", []),
            "keywords": entry.get("keywords", []),
            "proteinInteractions": entry.get("proteinInteractions"),
            "ftFeatures": entry.get("features", []),
            "references": [
                {
                    "citation": ref.get("citation"),
                    "referencePositions": ref.get("referencePositions", []),
                    "referenceComments": ref.get("referenceComments", []),
                }
                for ref in entry.get("references") or []
            ],
            "uniProtKBCrossReferences": [
                {
                    "database": x.get("database"),
                    "id": x.get("id"),
                    "properties": x.get("properties", []),
                }
                for x in entry.get("uniProtKBCrossReferences") or []
                if x.get("database") in {
                    "PDB", "GO", "Ensembl", "RefSeq", "AlphaFoldDB", "InterPro",
                    "Pfam", "Reactome", "ChEBI", "IntAct", "BioGRID", "MINT", "DIP",
                    "STRING", "ComplexPortal", "CORUM", "AlphaMissense",
                }
            ],
            "entryAudit": {
                "firstPublicDate": (entry.get("entryAudit") or {}).get("firstPublicDate"),
                "lastAnnotationUpdateDate": (entry.get("entryAudit") or {}).get(
                    "lastAnnotationUpdateDate"
                ),
                "lastSequenceUpdateDate": (entry.get("entryAudit") or {}).get(
                    "lastSequenceUpdateDate"
                ),
            },
            "entryVersion": entry.get("entryVersion"),
            "photochemistry": extract_photochemistry(entry),
            "_provenance": _provenance("uniprot", f"GET /uniprotkb/{acc}.json", acc),
        }

    def fasta(self, accession: str) -> dict[str, Any]:
        acc = accession.strip().upper()
        if not ACCESSION_RE.match(acc):
            raise RetrievalError(f"{accession!r} is not a well-formed UniProt accession.")
        response = self._client.get(f"/uniprotkb/{acc}.fasta")
        response.raise_for_status()
        text = response.text
        header, *seq_lines = text.splitlines()
        return {
            "accession": acc,
            "header": header,
            "sequence": "".join(s for s in seq_lines if not s.startswith(">")),
            "_provenance": _provenance("uniprot", f"GET /uniprotkb/{acc}.fasta", acc),
        }


class PdbeClient:
    """MCP-backed PDBe retrieval with live API/Search tool discovery."""

    def __init__(self, settings: Settings, *, mcp: Any | None = None) -> None:
        from .pdbe_mcp import PdbeMcpClient

        self.settings = settings
        self._mcp = mcp or PdbeMcpClient(settings)

    def tool_schemas(self) -> list[dict[str, Any]]:
        """Return OpenAI-compatible schemas generated from live MCP discovery."""
        return self._mcp.openai_tools()

    def supports(self, name: str) -> bool:
        return self._mcp.supports(name)

    def call(self, name: str, arguments: dict[str, Any]) -> dict[str, Any]:
        return self._mcp.call(name, arguments)

    def close(self) -> None:
        self._mcp.close()


def dispatch(
    name: str,
    arguments: dict[str, Any],
    *,
    uniprot: UniProtClient,
    pdbe: PdbeClient | None = None,
) -> dict[str, Any]:
    """Execute one tool call and return a JSON-serialisable result (or a limitation note)."""
    try:
        if name == SEARCH_UNIPROT:
            return uniprot.search(
                arguments["query"],
                fields=arguments.get("fields"),
                size=arguments.get("size", 10),
            )
        if name == FETCH_ENTRY:
            return uniprot.entry(arguments["accession"])
        if name == FETCH_FASTA:
            return uniprot.fasta(arguments["accession"])
        if name == PLOT_HYDROPATHY:
            from .plots import plot_hydropathy

            out = plot_hydropathy(
                arguments["sequence"],
                label=arguments.get("label") or "(user-supplied sequence)",
                window=int(arguments.get("window") or 19),
                out_path=arguments.get("out_path") or "out/hydropathy.png",
            )
            out["_provenance"] = _provenance(
                "local-compute", "plot_hydropathy", "sequence", scale="Kyte-Doolittle"
            )
            return out
        if name in {COMPUTE_PROPERTIES, CHROMOPHORE}:
            from . import properties as local

            sequence = arguments["sequence"]
            if name == CHROMOPHORE:
                return {
                    **local.chromophore_candidates(sequence),
                    "_provenance": _provenance(
                        "local-compute", "chromophore_candidates", "sequence"
                    ),
                }
            cleaned, warnings = local.clean_fasta_sequence(sequence)
            result = local.compute_properties(
                cleaned,
                header=arguments.get("label") or "(user-supplied sequence)",
            ).as_dict()
            result["warnings"] = warnings + result.get("warnings", [])
            result["input_length"] = len(cleaned)
            result["_provenance"] = _provenance("local-compute", "compute_properties", "sequence")
            return result
        if pdbe is not None and pdbe.supports(name):
            return pdbe.call(name, arguments)
    except (RetrievalError, PdbeMcpError, httpx.HTTPError, KeyError, json.JSONDecodeError) as exc:
        return {
            "error": str(exc),
            "_provenance": _provenance("tool-failure", name, json.dumps(arguments, sort_keys=True)),
        }
    raise RetrievalError(f"Unknown tool requested by the model: {name}")


__all__ = ["PdbeClient", "RetrievalError", "UniProtClient", "dispatch"]
