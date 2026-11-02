"""Prepare a verified PDBe-mapped chain for Foldseek's structure-only search.

A FASTA is not a 3D structure. Never submit one to Foldseek as though it were a CIF.
"""

from __future__ import annotations

import io
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import httpx
from Bio.PDB import MMCIFIO, MMCIFParser, Select
from Bio.PDB.PDBExceptions import PDBException

PDB_ID = re.compile(r"[0-9][A-Za-z0-9]{3}\Z")
MAX_CIF_BYTES = 8_000_000


class _ChainOnly(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_model(self, model: Any) -> bool:
        return model.id == 0

    def accept_chain(self, chain: Any) -> bool:
        return chain.id == self.chain_id


def prepare_query(
    pdbe: dict[str, Any], destination: Path, *, client: httpx.Client | None = None
) -> dict[str, Any]:
    """Save one mapped experimental chain as mmCIF; return provenance or a limitation."""
    if pdbe.get("status") == "unavailable":
        return {
            "status": "unavailable",
            "reason": "PDBe mapping was unavailable; Foldseek could not obtain a query structure.",
        }
    if pdbe.get("status") != "verified":
        return {
            "status": "skipped",
            "reason": "No verified PDBe-mapped experimental structure was available for Foldseek.",
        }
    candidates = [
        item for item in pdbe.get("structures") or []
        if PDB_ID.fullmatch(str(item.get("pdb_id") or ""))
        and isinstance(item.get("chain_id"), str)
        and item["chain_id"].strip()
    ]
    if not candidates:
        return {
            "status": "skipped",
            "reason": "PDBe returned no mapped structure with a chain ID for Foldseek.",
        }
    structure = candidates[0]  # PDBe already ranked these mappings for the verified accession.
    pdb_id = str(structure["pdb_id"]).lower()
    chain_id = structure["chain_id"].strip()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
    owned = client is None
    http = client or httpx.Client(timeout=15.0, follow_redirects=True)
    try:
        with http.stream("GET", url) as response:
            response.raise_for_status()
            content = bytearray()
            for chunk in response.iter_bytes():
                content.extend(chunk)
                if len(content) > MAX_CIF_BYTES:
                    return {
                        "status": "skipped",
                        "reason": "Mapped structure exceeds the 8 MB query limit.",
                    }
        text = content.decode("utf-8")
        if not text.lstrip().startswith("data_"):
            raise ValueError("PDBe did not return an mmCIF structure")
        parsed = MMCIFParser(QUIET=True).get_structure(pdb_id, io.StringIO(text))
        first_model = next(parsed.get_models())
        if chain_id not in first_model or not list(first_model[chain_id].get_atoms()):
            return {
                "status": "skipped",
                "reason": "The PDBe-mapped chain was not present in the downloaded structure.",
            }
        writer = MMCIFIO()
        writer.set_structure(parsed)
        writer.save(str(destination), select=_ChainOnly(chain_id))
        if destination.stat().st_size > 2_000_000:
            destination.unlink()
            return {
                "status": "skipped",
                "reason": "Mapped chain exceeds the 2 MB Foldseek query limit.",
            }
        return {
            "status": "ready", "pdb_id": pdb_id.upper(), "chain_id": chain_id,
            "source_url": url,
            "entry_url": structure.get("source_url"),
            "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "mapped_coverage": structure.get("coverage_fraction"),
            "mapping_source": pdbe.get("_provenance"),
        }
    except (
        httpx.HTTPError, UnicodeError, ValueError, StopIteration, OSError, KeyError,
        IndexError, TypeError, PDBException,
    ):
        destination.unlink(missing_ok=True)
        return {
            "status": "unavailable",
            "reason": "The mapped PDBe structure could not be fetched or prepared for Foldseek.",
        }
    finally:
        if owned:
            http.close()
