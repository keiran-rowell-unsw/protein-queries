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
from .tools import FETCH_ENTRY, FETCH_FASTA, SEARCH_UNIPROT

ACCESSION_RE = re.compile(r"^[A-Z0-9]{6,10}$|^[A-Z0-9]{5,10}-\d{2}$")
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
        if response.status_code == 404:
            raise RetrievalError(
                f"UniProt has no entry {acc} (404). Treat the identity as unresolved; do not guess."
            )
        response.raise_for_status()
        entry = response.json()
        sequence = entry.get("sequence", {}) or {}
        return {
            # Trimmed to what answers need; keeps prompts and latency bounded (uniprot/AGENTS.md).
            "primaryAccession": entry.get("primaryAccession", acc),
            "entryType": entry.get("entryType"),
            "organism": entry.get("organism", {}).get("scientificName"),
            "proteinDescription": entry.get("proteinDescription"),
            "genes": entry.get("genes", []),
            "sequence": {
                "value": sequence.get("value", ""),
                "length": sequence.get("length", len(sequence.get("value", ""))),
            },
            "comments": entry.get("comments", []),
            "proteinInteractions": entry.get("proteinInteractions"),
            "ftFeatures": entry.get("features", []),
            "uniProtKBCrossReferences": [
                {
                    "database": x.get("database"),
                    "id": x.get("id"),
                    "properties": x.get("properties", []),
                }
                for x in entry.get("uniProtKBCrossReferences", [])
                if x.get("database") in {"PDB", "GO", "PubMed", "Ensembl", "RefSeq"}
            ],
            "entryVersion": entry.get("entryVersion"),
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
    """MCP-backed structural retrieval. Not implemented and not installed on this node.

    Connect with the stdio config in `pdbe-kb/AGENTS.md` (`uvx pdbe-mcp-server --server-type
    pdbe_api_server` and `pdbe_search_server`), discover tool names from the server, and
    implement the same `_provenance` discipline as `UniProtClient`.
    """

    def __init__(self, settings: Settings) -> None:
        self.settings = settings

    def search_structures(self, uniprot_accession: str, *, rows: int = 10) -> dict[str, Any]:
        raise NotImplementedError(
            "PDBe MCP is not configured here. See pdbe-kb/AGENTS.md for the server "
            "registration and schema-discovery rules before implementing."
        )

    def entry(self, pdb_id: str) -> dict[str, Any]:
        raise NotImplementedError("See pdbe-kb/AGENTS.md: fetch via the PDBe API MCP server.")


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
        if pdbe is not None and name.startswith("pdbe"):
            return pdbe.entry(arguments.get("pdb_id", ""))
    except (RetrievalError, httpx.HTTPError, KeyError, json.JSONDecodeError) as exc:
        return {
            "error": str(exc),
            "_provenance": _provenance("tool-failure", name, json.dumps(arguments, sort_keys=True)),
        }
    raise RetrievalError(f"Unknown tool requested by the model: {name}")


__all__ = ["PdbeClient", "RetrievalError", "UniProtClient", "dispatch"]
