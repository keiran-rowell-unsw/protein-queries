"""Deterministic PDBe structure summary for a sequence-verified UniProt accession."""

from __future__ import annotations

import json
import re
from typing import Any

from .config import Settings
from .pdbe_mcp import PdbeMcpError
from .retrieval import PdbeClient

PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}\Z")
MAX_STRUCTURES = 10
_BEST_STRUCTURE_DESCRIPTION = "Provides a ranked list of PDB structures mapped"
_SEARCH_SCHEMA_DESCRIPTION = "Retrieves the Solr search schema"
_SEARCH_QUERY_DESCRIPTION = "Executes a search query against the PDBe Solr search service"
_SEARCH_FIELDS = (
    "pdb_id,chain_id,entity_id,experimental_method,resolution,uniprot_coverage,"
    "prefered_assembly_id,overall_quality"
)


def _search_documents(result: dict[str, Any]) -> list[dict[str, Any]]:
    """Parse the bounded, human-readable document list returned by PDBe Search MCP."""
    texts = [
        item.get("text") for item in result.get("content", [])
        if isinstance(item, dict) and item.get("type") == "text" and item.get("text")
    ]
    if not texts:
        return []
    documents: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in texts[0].splitlines():
        if re.fullmatch(r"Document \d+:", line):
            current = {}
            documents.append(current)
            continue
        if current is None or not line.startswith("  ") or ": " not in line:
            continue
        key, raw = line.strip().split(": ", 1)
        value: Any = raw
        if re.fullmatch(r"-?\d+", raw):
            value = int(raw)
        elif re.fullmatch(r"-?(?:\d+\.\d*|\d*\.\d+)", raw):
            value = float(raw)
        current[key] = value
    return documents


def _fetch_bounded_search_fallback(
    accession: str,
    *,
    lookup: PdbeClient,
    schemas: list[dict[str, Any]],
    limit: int,
) -> dict[str, Any]:
    """Retry oversized ranked mappings as a field- and row-bounded PDBe search."""
    schema_tools = [
        item for item in schemas
        if (item.get("function") or {}).get("description", "").lstrip().startswith(
            _SEARCH_SCHEMA_DESCRIPTION
        )
    ]
    query_tools = [
        item for item in schemas
        if (item.get("function") or {}).get("description", "").lstrip().startswith(
            _SEARCH_QUERY_DESCRIPTION
        )
    ]
    if len(schema_tools) != 1 or len(query_tools) != 1:
        raise PdbeMcpError("PDBe Search did not expose an unambiguous bounded fallback")
    lookup.call(schema_tools[0]["function"]["name"], {})
    rows = max(1, min(int(limit), MAX_STRUCTURES))
    result = lookup.call(query_tools[0]["function"]["name"], {
        "query": f"uniprot_accession:{accession}",
        "fl": _SEARCH_FIELDS,
        "sort": "overall_quality desc",
        "rows": rows,
    })
    search_text = "\n".join(
        str(item.get("text") or "") for item in result.get("content", [])
        if isinstance(item, dict)
    )
    if result.get("isError") is True or search_text.startswith("Search query failed:"):
        raise PdbeMcpError("PDBe Search returned an error for the bounded fallback")
    records = _search_documents(result)
    structures: list[dict[str, Any]] = []
    for record in records[:rows]:
        pdb_id = str(record.get("pdb_id") or "").upper()
        if not PDB_ID.fullmatch(pdb_id):
            continue
        structures.append({
            "pdb_id": pdb_id,
            "chain_id": record.get("chain_id"),
            "entity_id": record.get("entity_id"),
            "experimental_method": record.get("experimental_method"),
            "resolution_angstrom": record.get("resolution"),
            "coverage_fraction": record.get("uniprot_coverage"),
            "uniprot_start": None,
            "uniprot_end": None,
            "observed_regions": [],
            "preferred_assembly_id": record.get("prefered_assembly_id"),
            "source_url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id.lower()}",
        })
    if not structures:
        return {
            "status": "not_found",
            "reason": "No experimental structures were returned under this accession.",
        }
    return {
        "status": "verified",
        "accession": accession,
        "structures": structures,
        "display_limit": MAX_STRUCTURES,
        "ranking": "PDBe Search overall quality, after the ranked mapping exceeded its cap",
        "_provenance": result.get("_provenance"),
    }


def fetch_pdbe_structures(
    accession: str,
    *,
    client: PdbeClient | None = None,
    limit: int = MAX_STRUCTURES,
) -> dict[str, Any]:
    """Fetch a bounded ranked structure list using a semantically discovered API tool."""
    owned = client is None
    lookup = client or PdbeClient(Settings.from_env(require_key=False))
    try:
        schemas = lookup.tool_schemas()
        candidates = [
            item for item in schemas
            if (item.get("function") or {}).get("description", "").startswith(
                _BEST_STRUCTURE_DESCRIPTION
            )
        ]
        if len(candidates) != 1:
            return {
                "status": "unavailable",
                "reason": (
                    "The live PDBe API tools did not expose one unambiguous "
                    "ranked-structure lookup."
                ),
            }
        function = candidates[0]["function"]
        parameters = function.get("parameters") or {}
        required = parameters.get("required") or []
        if len(required) != 1:
            return {
                "status": "unavailable",
                "reason": "The ranked-structure tool advertised an unsupported input schema.",
            }
        try:
            result = lookup.call(function["name"], {required[0]: accession})
        except PdbeMcpError as exc:
            if "above the configured" not in str(exc):
                raise
            return _fetch_bounded_search_fallback(
                accession, lookup=lookup, schemas=schemas, limit=limit
            )
        if result.get("isError") is True:
            return {
                "status": "unavailable",
                "reason": "PDBe returned an error while retrieving mapped structures.",
            }
        text_items = [
            item.get("text") for item in result.get("content", [])
            if isinstance(item, dict) and item.get("type") == "text" and item.get("text")
        ]
        if not text_items:
            return {"status": "not_found", "reason": "PDBe returned no mapped structures."}
        response_text = text_items[0]
        if response_text.startswith("API request failed:"):
            if "404" in response_text or "does not contain any data" in response_text:
                return {
                    "status": "not_found",
                    "reason": (
                        "No experimental structures were found in PDBe under this "
                        "UniProt accession."
                    ),
                }
            status = next(
                (line for line in response_text.splitlines() if line.startswith("Status Code:")),
                "status unavailable",
            )
            return {
                "status": "unavailable",
                "reason": "The PDBe API request failed (" + status.lower() + ").",
            }
        if response_text.startswith("Input validation error:"):
            return {
                "status": "unavailable",
                "reason": "The live PDBe tool rejected the structure lookup request.",
            }
        payload = json.loads(response_text)
        records = payload.get(accession) or payload.get(accession.upper()) or []
        if not isinstance(records, list) or not records:
            return {
                "status": "not_found",
                "reason": "No experimental structures were returned under this accession.",
            }
        structures: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()
        for record in records:
            if not isinstance(record, dict):
                continue
            pdb_id = str(record.get("pdb_id") or "").upper()
            chain = str(record.get("chain_id") or "")
            if not PDB_ID.fullmatch(pdb_id) or (pdb_id, chain) in seen:
                continue
            seen.add((pdb_id, chain))
            structures.append({
                "pdb_id": pdb_id,
                "chain_id": chain or None,
                "entity_id": record.get("entity_id"),
                "experimental_method": record.get("experimental_method"),
                "resolution_angstrom": record.get("resolution"),
                "coverage_fraction": record.get("coverage"),
                "uniprot_start": record.get("unp_start"),
                "uniprot_end": record.get("unp_end"),
                "observed_regions": record.get("observed_regions") or [],
                "preferred_assembly_id": record.get("preferred_assembly_id"),
                "source_url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id.lower()}",
            })
            if len(structures) >= max(1, min(int(limit), MAX_STRUCTURES)):
                break
        if not structures:
            return {
                "status": "not_found",
                "reason": "PDBe returned no valid mapped PDB records under this accession.",
            }
        return {
            "status": "verified",
            "accession": accession,
            "structures": structures,
            "display_limit": MAX_STRUCTURES,
            "ranking": "PDBe ranked structures mapped to the verified UniProt accession",
            "_provenance": result.get("_provenance"),
        }
    except PdbeMcpError as exc:
        return {
            "status": "unavailable",
            "reason": "PDBe structure retrieval was unavailable: " + str(exc),
        }
    except (json.JSONDecodeError, KeyError, TypeError, ValueError):
        return {
            "status": "unavailable",
            "reason": "PDBe returned an unexpected structure-response format.",
        }
    finally:
        if owned:
            lookup.close()


__all__ = ["fetch_pdbe_structures"]
