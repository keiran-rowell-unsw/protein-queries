"""Small client for the public Foldseek search-server API."""

from __future__ import annotations

import json
import re
import time
from datetime import UTC, datetime
from typing import Any

import httpx

from .config import Settings
from .retrieval import RetrievalError


def _provenance(source: str, endpoint: str, subject: str, **extra: Any) -> dict[str, Any]:
    record = {
        "source": source,
        "endpoint": endpoint,
        "subject": subject,
        "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    record.update(extra)
    return record

_DEFAULT_DATABASE = "pdb100"
_DATABASE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,80}$")
_MODE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,40}$")
_MAX_RESULT_BYTES = 1_000_000


class FoldseekClient:
    """Submit a query, poll its ticket, and return the bounded result payload."""

    def __init__(self, settings: Settings, *, client: httpx.Client | None = None) -> None:
        self.settings = settings
        self._owns = client is None
        self._client = client or httpx.Client(
            base_url=settings.foldseek_base_url,
            timeout=settings.request_timeout,
            headers={"Accept": "application/json"},
        )

    def close(self) -> None:
        if self._owns:
            self._client.close()

    def search(
        self,
        query: str,
        *,
        database: str = _DEFAULT_DATABASE,
        mode: str | None = None,
    ) -> dict[str, Any]:
        query = query.strip()
        database = database.strip()
        mode = mode.strip() if mode else "3diaa"  # required by the live Foldseek API
        if not query:
            raise RetrievalError("Foldseek query is empty.")
        # Foldseek's structural-search endpoint accepts PDB/mmCIF coordinates, not
        # amino-acid FASTA. Reject sequence input rather than report a false search.
        if not query.lstrip().startswith("data_") and not any(
            line.startswith(("ATOM  ", "HETATM")) for line in query.splitlines()
        ):
            raise RetrievalError(
                "Foldseek requires PDB or mmCIF coordinates; FASTA is not a structure."
            )
        if len(query) > 2_000_000:
            raise RetrievalError("Foldseek query exceeds the 2 MB limit.")
        if not _DATABASE_RE.fullmatch(database):
            raise RetrievalError("Foldseek database contains unsupported characters.")
        if mode and not _MODE_RE.fullmatch(mode):
            raise RetrievalError("Foldseek mode contains unsupported characters.")

        data = {"database[]": database, "mode": mode}
        try:
            response = self._client.post(
                "/api/ticket",
                files={
                    "q": (
                        "query.cif" if query.lstrip().startswith("data_") else "query.pdb",
                        query.encode("utf-8"),
                        (
                            "chemical/x-mmcif" if query.lstrip().startswith("data_")
                            else "chemical/x-pdb"
                        ),
                    )
                },
                data=data,
            )
            response.raise_for_status()
            ticket_payload = response.json()
        except (httpx.HTTPError, ValueError) as exc:
            raise RetrievalError(f"Foldseek submission failed: {exc}") from exc

        ticket = str(ticket_payload.get("id") or ticket_payload.get("ticket") or "")
        if not ticket:
            return self._result(ticket_payload, database, mode, status="submission_error")

        deadline = time.monotonic() + min(max(self.settings.request_timeout, 10.0), 180.0)
        status_payload = ticket_payload
        while True:
            status = str(status_payload.get("status") or "").upper()
            if status in {"COMPLETE", "FAILED", "ERROR", "UNKNOWN", "RATELIMIT", "MAINTENANCE"}:
                break
            if time.monotonic() >= deadline:
                return self._result(
                    status_payload, database, mode, status="timeout", ticket=ticket
                )
            time.sleep(1.0)
            try:
                status_response = self._client.get(f"/api/ticket/{ticket}")
                status_response.raise_for_status()
                status_payload = status_response.json()
            except (httpx.HTTPError, ValueError) as exc:
                raise RetrievalError(f"Foldseek status query failed: {exc}") from exc

        if status != "COMPLETE":
            return self._result(
                status_payload, database, mode, status=status.lower(), ticket=ticket
            )
        try:
            result_response = self._client.get(
                f"/api/result/{ticket}/0",
                params={"format": "brief", "database": database},
            )
            result_response.raise_for_status()
            result = self._summarise_hits(result_response.json(), database)
        except (httpx.HTTPError, ValueError) as exc:
            raise RetrievalError(f"Foldseek result query failed: {exc}") from exc
        return self._result(result, database, mode, status="complete", ticket=ticket)

    @staticmethod
    def _summarise_hits(payload: Any, database: str) -> dict[str, Any]:
        """Preserve bounded server-ranked candidates, not raw alignment/coordinate dumps."""
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            raise RetrievalError("Foldseek returned an unexpected hit format.")
        if len(payload.get("queries") or []) != 1:
            raise RetrievalError("Foldseek did not return exactly one structural query.")
        fields = (
            "target", "seqId", "alnLength", "qStartPos", "qEndPos",
            "dbStartPos", "dbEndPos", "prob", "eval", "score", "qLen", "dbLen",
        )
        datasets = []
        for group in payload["results"]:
            if not isinstance(group, dict) or group.get("db") != database:
                continue
            nested = group.get("alignments") or []
            if not isinstance(nested, list):
                raise RetrievalError("Foldseek alignments had an unexpected format.")
            hits = [hit for partition in nested if isinstance(partition, list)
                    for hit in partition if isinstance(hit, dict)]
            datasets.append({
                "database": database,
                "total_hits_reported": len(hits),
                "hits": [{key: hit[key] for key in fields if key in hit} for hit in hits[:5]],
            })
        if len(datasets) != 1:
            raise RetrievalError("Foldseek did not return one matching database result.")
        return {"mode": payload.get("mode"), "query_count": 1, **datasets[0]}

    @staticmethod
    def _result(
        payload: Any,
        database: str,
        mode: str | None,
        *,
        status: str,
        ticket: str | None = None,
    ) -> dict[str, Any]:
        bounded = payload
        try:
            raw = json.dumps(payload, ensure_ascii=False)
            if len(raw.encode("utf-8")) > _MAX_RESULT_BYTES:
                bounded = {"truncated": True, "original_bytes": len(raw.encode("utf-8"))}
                if isinstance(payload, dict):
                    bounded.update(
                        {key: payload[key] for key in ("status", "id", "ticket") if key in payload}
                    )
        except (TypeError, ValueError):
            bounded = {"unserialisable_payload": True}
        return {
            "status": status,
            "ticket": ticket,
            "database": database,
            "mode": mode,
            "result": bounded,
            "_provenance": _provenance(
                "foldseek",
                "POST /api/ticket -> GET /api/result/{ticket}",
                database,
                ticket=ticket,
                status=status,
                retrieved_at=datetime.now(UTC).isoformat(timespec="seconds"),
            ),
        }


__all__ = ["FoldseekClient"]
