"""Optional UniProt annotations for an exact canonical-sequence match only.

An accession is supplied by the researcher; no sequence-to-identity guess is made.
User-facing annotations are drawn exclusively from the fetched UniProt entry in this run.
"""

from __future__ import annotations

import re
import time
from typing import Any

import httpx

from .config import Settings
from .retrieval import RetrievalError, UniProtClient

ACCESSION = re.compile(r"[A-Z0-9]{6,10}\Z")
ANNOTATION_TYPES = {"Modified residue", "Glycosylation", "Signal", "Cross-link"}


def fetch_annotations(
    sequence: str,
    accession: str,
    *,
    client: UniProtClient | None = None,
) -> dict[str, Any]:
    """Return verified entry facts or an explicit mismatch/failure; never partial facts."""
    acc = accession.strip().upper()
    if not ACCESSION.fullmatch(acc):
        return {"status": "invalid", "reason": "Enter a valid UniProt accession."}
    owned = client is None
    lookup = client or UniProtClient(Settings())
    try:
        for attempt in range(2):
            try:
                entry = lookup.entry(acc)
                break
            except httpx.HTTPStatusError as exc:
                if attempt == 0 and exc.response.status_code in {429, 500, 502, 503, 504}:
                    time.sleep(0.5)
                    continue
                raise
        primary = entry.get("primaryAccession") or acc
        if primary != acc:
            return {
                "status": "mismatch",
                "reason": (
                    "This accession resolves to a different primary accession. "
                    "Confirm the current UniProt entry before attaching annotations."
                ),
            }
        canonical = (entry.get("sequence") or {}).get("value") or ""
        if not canonical or canonical != sequence:
            return {
                "status": "mismatch",
                "reason": (
                    "The submitted sequence does not exactly match this entry's canonical "
                    "sequence. UniProt annotations were not attached; check isoforms, "
                    "constructs, tags, and variants."
                ),
            }
        features = []
        for feature in entry.get("ftFeatures") or []:
            if feature.get("type") not in ANNOTATION_TYPES:
                continue
            location = feature.get("location") or {}
            features.append({
                "type": feature["type"],
                "start": (location.get("start") or {}).get("value"),
                "end": (location.get("end") or {}).get("value"),
                "description": feature.get("description"),
                "evidence": [
                    {"code": item.get("evidenceCode"), "source": item.get("source"),
                     "id": item.get("id")}
                    for item in feature.get("evidences") or [] if isinstance(item, dict)
                ],
            })
        absorption = (entry.get("photochemistry") or {}).get("absorption")
        return {
            "status": "verified",
            "accession": entry.get("primaryAccession") or acc,
            "organism": entry.get("organism"),
            "entry_type": entry.get("entryType"),
            "entry_version": entry.get("entryVersion"),
            "retrieved_at": (entry.get("_provenance") or {}).get("retrieved_at"),
            "canonical_sequence_match": True,
            "features": features,
            "absorption": absorption,
            "source_url": (
                "https://www.uniprot.org/uniprotkb/"
                + (entry.get("primaryAccession") or acc) + "/entry"
            ),
        }
    except (RetrievalError, httpx.HTTPError, ValueError, KeyError) as exc:
        if isinstance(exc, RetrievalError) and "404" in str(exc):
            return {"status": "not_found", "reason": "No UniProt entry found for this accession."}
        return {"status": "unavailable", "reason": "UniProt retrieval is temporarily unavailable."}
    finally:
        if owned:
            lookup.close()
