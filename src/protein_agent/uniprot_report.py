"""Optional UniProt annotations for an exact canonical-sequence match only.

An accession is supplied by the researcher; no sequence-to-identity guess is made.
User-facing annotations are drawn exclusively from the fetched UniProt entry in this run.
"""

from __future__ import annotations

import re
import time
from datetime import UTC, datetime
from typing import Any

import httpx
from Bio.SeqUtils.CheckSum import crc64

from .config import Settings
from .retrieval import RetrievalError, UniProtClient

# Two validators with two different jobs; one regex cannot serve both honestly.
# USER_ACCESSION accepts what a person may type, including isoform suffixes.
# ACCESSION filters the bare accessions embedded in UniParc `CRC.accession` strings,
# which are split on "." upstream and must never carry a hyphen or isoform suffix.
USER_ACCESSION = re.compile(r"(?:[A-Z0-9]{6,10}|[A-Z0-9]{5,10}-[0-9]{1,2})\Z")
ACCESSION = re.compile(r"[A-Z0-9]{6,10}\Z")
ANNOTATION_TYPES = {"Modified residue", "Glycosylation", "Signal", "Cross-link"}
MAX_CANDIDATES = 4


def fetch_annotations(
    sequence: str,
    accession: str,
    *,
    client: UniProtClient | None = None,
    max_attempts: int = 2,
) -> dict[str, Any]:
    """Return verified entry facts or an explicit mismatch/failure; never partial facts."""
    acc = accession.strip().upper()
    if not USER_ACCESSION.fullmatch(acc):
        return {"status": "invalid", "reason": "Enter a valid UniProt accession."}
    owned = client is None
    lookup = client or UniProtClient(Settings.from_env(require_key=False))
    try:
        for attempt in range(max_attempts):
            try:
                entry = lookup.entry(acc)
                break
            except httpx.HTTPStatusError as exc:
                if (
                    attempt + 1 < max_attempts
                    and exc.response.status_code in {429, 500, 502, 503, 504}
                ):
                    time.sleep(0.5)
                    continue
                raise
        primary = entry.get("primaryAccession") or acc
        # UniProt's `primaryAccession` is not a stable identity signal for isoforms. Live
        # checks: Q15465-1 reports 'Q15465' (stripped) because that isoform is identical to
        # canonical at 462 aa, while P04637-2 and P12110-2 echo '-2' since they genuinely
        # differ. A strict != therefore rejects legitimate isoform submissions depending on
        # whether the isoform happens to match canonical.
        # The exact-sequence comparison below is the authoritative gate, so accept the
        # accession when it matches either as submitted or modulo an isoform suffix, and let
        # the sequence check catch any real mismatch.
        if primary not in {acc, acc.split("-", 1)[0]}:
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
        if isinstance(exc, RetrievalError):
            if "404" in str(exc):
                return {
                    "status": "not_found",
                    "reason": "No UniProt entry found for this accession.",
                }
            if "inactive or redirected" in str(exc):
                return {
                    "status": "inactive",
                    "reason": (
                        "This accession is inactive or redirected; "
                        "confirm its current primary entry."
                    ),
                }
        return {"status": "unavailable", "reason": "UniProt retrieval is temporarily unavailable."}
    finally:
        if owned:
            lookup.close()


def find_sequence_annotations(
    sequence: str, *, settings: Settings | None = None
) -> dict[str, Any]:
    """Search UniParc by Biopython CRC64; verify actual residues and each KB entry.

    No fuzzy matching or first-hit selection. Multiple plausible entries remain candidates
    for the user to resolve; failures are not misreported as absence of a match.
    """
    checksum = crc64(sequence).removeprefix("CRC-")
    settings = settings or Settings.from_env(require_key=False)
    try:
        with httpx.Client(base_url=settings.uniprot_base_url, timeout=10.0) as http:
            response = http.get(
                "/uniparc/search",
                params={"query": f"checksum:{checksum}", "format": "json", "size": 5},
            )
            response.raise_for_status()
            payload = response.json()
            if int(response.headers.get("x-total-results", "0")) > 5:
                return {
                    "status": "ambiguous",
                    "reason": (
                        "Too many checksum hits to verify safely; "
                        "enter a known UniProt accession."
                    ),
                }
            matches = [
                item for item in payload.get("results", [])
                if (item.get("sequence") or {}).get("value") == sequence
                and (item.get("sequence") or {}).get("length") == len(sequence)
            ]
            if not matches:
                return {
                    "status": "not_found",
                    "reason": (
                        "No exact sequence found in UniParc under these criteria; "
                        "UniProt identity is unresolved."
                    ),
                }
            candidates = sorted({
                acc.split(".", 1)[0]
                for item in matches
                for acc in item.get("uniProtKBAccessions") or []
                if ACCESSION.fullmatch(acc.split(".", 1)[0])
            })
            if not candidates:
                return {
                    "status": "not_found",
                    "reason": (
                        "An exact sequence was found in UniParc, but no UniProtKB "
                        "accession was listed."
                    ),
                }
            if len(candidates) > MAX_CANDIDATES:
                return {
                    "status": "ambiguous", "candidates": candidates[:8],
                    "reason": (
                        "Several UniProtKB accessions are linked to this exact sequence. "
                        "Enter the accession you want to annotate."
                    ),
                }
            # Bound entry lookups to four calls of at most eight seconds each.
            with httpx.Client(base_url=settings.uniprot_base_url, timeout=8.0) as kb_http:
                kb = UniProtClient(settings, client=kb_http)
                results = [
                    fetch_annotations(sequence, acc, client=kb, max_attempts=1)
                    for acc in candidates
                ]
            verified = [item for item in results if item["status"] == "verified"]
            if any(item["status"] == "unavailable" for item in results):
                return {
                    "status": "unavailable",
                    "reason": "Entry verification was unavailable; no accession was selected.",
                }
            if len(verified) == 1:
                entry = verified[0]
                entry["lookup"] = {
                    "method": (
                        "UniParc CRC64 checksum search; full sequence and "
                        "canonical entry verified"
                    ),
                    "checksum": checksum,
                    "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
                }
                return entry
            if len(verified) > 1:
                return {
                    "status": "ambiguous",
                    "candidates": [item["accession"] for item in verified],
                    "reason": (
                        "More than one canonical UniProtKB entry matches this exact "
                        "sequence. Enter an accession to choose which entry to annotate."
                    ),
                }
            return {
                "status": "not_found", "candidates": candidates,
                "reason": (
                    "This sequence was found in UniParc, but no linked current UniProtKB "
                    "canonical sequence matched exactly. Try an isoform-specific lookup "
                    "or enter a known accession."
                ),
            }
    except (httpx.HTTPError, ValueError, KeyError, TypeError, AttributeError) as exc:
        return {
            "status": "unavailable",
            "reason": (
                f"Sequence lookup was unavailable ({type(exc).__name__}); "
                "annotations were not attached."
            ),
        }
