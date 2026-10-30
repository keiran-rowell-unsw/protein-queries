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
        result = lookup.call(function["name"], {required[0]: accession})
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
