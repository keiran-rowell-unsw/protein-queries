"""Small client for the public Foldseek search-server API."""

from __future__ import annotations

import json
import re
import time
from collections.abc import Sequence
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


_DATABASE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,80}$")
_MODE_RE = re.compile(r"^[A-Za-z0-9_.-]{1,40}$")
_MAX_DATABASES = 50
_MAX_RESULT_BYTES = 1_000_000


class FoldseekClient:
    """Discover databases, submit a query, poll its ticket, and bound its results."""

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

    def _discover_databases(self) -> list[str]:
        """Return every live standard protein-search database offered by the server UI."""
        try:
            response = self._client.get("/api/databases/all")
            response.raise_for_status()
            payload = response.json()
        except (httpx.HTTPError, ValueError) as exc:
            raise RetrievalError(f"Foldseek database discovery failed: {exc}") from exc
        records = payload.get("databases") if isinstance(payload, dict) else None
        if not isinstance(records, list):
            raise RetrievalError("Foldseek returned an unexpected database list.")
        databases = [
            str(record["path"])
            for record in records
            if isinstance(record, dict)
            and record.get("status") == "COMPLETE"
            and record.get("default") is True
            and not record.get("interface")
            and not record.get("motif")
            and not record.get("rna")
            and _DATABASE_RE.fullmatch(str(record.get("path") or ""))
        ]
        if not databases:
            raise RetrievalError("Foldseek reported no available protein-search databases.")
        return databases[:_MAX_DATABASES]

    def search(
        self,
        query: str,
        *,
        database: str | Sequence[str] | None = None,
        mode: str | None = None,
    ) -> dict[str, Any]:
        query = query.strip()
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
        if mode and not _MODE_RE.fullmatch(mode):
            raise RetrievalError("Foldseek mode contains unsupported characters.")

        if database is None:
            databases = self._discover_databases()
        elif isinstance(database, str):
            databases = [database.strip()]
        else:
            databases = [str(item).strip() for item in database]
        databases = list(dict.fromkeys(databases))
        if not databases or len(databases) > _MAX_DATABASES:
            raise RetrievalError("Foldseek requires between 1 and 50 databases.")
        if not all(_DATABASE_RE.fullmatch(item) for item in databases):
            raise RetrievalError("Foldseek database contains unsupported characters.")

        data = {"database[]": databases, "mode": mode}
        try:
            response = self._client.post(
                "/api/ticket",
                files={
                    "q": (
                        "query.cif" if query.lstrip().startswith("data_") else "query.pdb",
                        query.encode("utf-8"),
                        (
                            "chemical/x-mmcif"
                            if query.lstrip().startswith("data_")
                            else "chemical/x-pdb"
                        ),
                    )
                },
                data=data,
            )
            if response.status_code == 429:
                try:
                    rate_payload = response.json()
                except ValueError:
                    rate_payload = {
                        "status": "RATELIMIT",
                        "reason": "The Foldseek service rate-limited this search.",
                    }
                return self._result(rate_payload, databases, mode, status="ratelimit")
            response.raise_for_status()
            ticket_payload = response.json()
        except httpx.HTTPStatusError as exc:
            raise RetrievalError(
                f"Foldseek submission failed: HTTP {exc.response.status_code}"
            ) from exc
        except (httpx.HTTPError, ValueError) as exc:
            raise RetrievalError(f"Foldseek submission failed: {exc}") from exc

        if not isinstance(ticket_payload, dict):
            raise RetrievalError("Foldseek submission returned an unexpected response.")
        ticket = str(ticket_payload.get("id") or ticket_payload.get("ticket") or "")
        if not ticket:
            status = str(ticket_payload.get("status") or "submission_error").lower()
            return self._result(ticket_payload, databases, mode, status=status)

        poll_timeout = min(max(self.settings.foldseek_poll_timeout, 10.0), 600.0)
        deadline = time.monotonic() + poll_timeout
        status_payload = ticket_payload
        poll_interval = 1.0
        while True:
            status = str(status_payload.get("status") or "").upper()
            if status in {
                "COMPLETE",
                "FAILED",
                "ERROR",
                "UNKNOWN",
                "RATELIMIT",
                "MAINTENANCE",
            }:
                break
            if time.monotonic() >= deadline:
                return self._result(
                    status_payload, databases, mode, status="timeout", ticket=ticket
                )
            time.sleep(poll_interval)
            poll_interval = min(5.0, poll_interval * 1.5)
            try:
                status_response = self._client.get(f"/api/ticket/{ticket}")
                status_response.raise_for_status()
                status_payload = status_response.json()
            except (httpx.HTTPError, ValueError) as exc:
                raise RetrievalError(f"Foldseek status query failed: {exc}") from exc

        if status != "COMPLETE":
            return self._result(
                status_payload, databases, mode, status=status.lower(), ticket=ticket
            )
        try:
            result_response = self._client.get(f"/api/result/{ticket}/0", params={"format": "json"})
            result_response.raise_for_status()
            result = self._summarise_hits(result_response.json(), databases)
        except (httpx.HTTPError, ValueError) as exc:
            raise RetrievalError(f"Foldseek result query failed: {exc}") from exc
        return self._result(result, databases, mode, status="complete", ticket=ticket)

    @staticmethod
    def _summarise_hits(payload: Any, databases: Sequence[str]) -> dict[str, Any]:
        """Merge bounded candidates across databases, ranked by reported E-value."""
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            raise RetrievalError("Foldseek returned an unexpected hit format.")
        if len(payload.get("queries") or []) != 1:
            raise RetrievalError("Foldseek did not return exactly one structural query.")
        fields = (
            "target",
            "seqId",
            "alnLength",
            "qStartPos",
            "qEndPos",
            "dbStartPos",
            "dbEndPos",
            "prob",
            "eval",
            "score",
            "qLen",
            "dbLen",
        )
        requested = set(databases)
        hits: list[dict[str, Any]] = []
        totals: dict[str, int] = {}
        for group in payload["results"]:
            if not isinstance(group, dict) or group.get("db") not in requested:
                continue
            db = str(group["db"])
            nested = group.get("alignments") or []
            if not isinstance(nested, list):
                raise RetrievalError("Foldseek alignments had an unexpected format.")
            group_hits = [
                hit
                for partition in nested
                if isinstance(partition, list)
                for hit in partition
                if isinstance(hit, dict)
            ]
            totals[db] = len(group_hits)
            hits.extend(
                {"database": db, **{key: hit[key] for key in fields if key in hit}}
                for hit in group_hits
            )
        if not totals:
            raise RetrievalError("Foldseek returned no matching database results.")

        def rank(hit: dict[str, Any]) -> tuple[float, float, float]:
            try:
                evalue = float(hit.get("eval", float("inf")))
            except (TypeError, ValueError):
                evalue = float("inf")
            try:
                probability = float(hit.get("prob", 0.0))
            except (TypeError, ValueError):
                probability = 0.0
            try:
                score = float(hit.get("score", 0.0))
            except (TypeError, ValueError):
                score = 0.0
            return evalue, -probability, -score

        hits.sort(key=rank)
        return {
            "mode": payload.get("mode"),
            "query_count": 1,
            "databases": list(databases),
            "database_hit_counts": totals,
            "total_hits_reported": sum(totals.values()),
            "hits": hits[:5],
        }

    @staticmethod
    def _result(
        payload: Any,
        databases: Sequence[str],
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
        database_list = list(databases)
        result = {
            "status": status,
            "ticket": ticket,
            "database": database_list[0] if len(database_list) == 1 else None,
            "databases": database_list,
            "mode": mode,
            "result": bounded,
            "_provenance": _provenance(
                "foldseek",
                "POST /api/ticket -> GET /api/result/{ticket}",
                ",".join(database_list),
                ticket=ticket,
                status=status,
                retrieved_at=datetime.now(UTC).isoformat(timespec="seconds"),
            ),
        }
        if isinstance(payload, dict) and payload.get("reason"):
            result["reason"] = str(payload["reason"])
        elif status == "timeout":
            result["reason"] = "Foldseek did not finish before the polling deadline."
        elif status not in {"complete", "pending", "running"}:
            result["reason"] = f"Foldseek ended with status {status}."
        return result


__all__ = ["FoldseekClient"]
