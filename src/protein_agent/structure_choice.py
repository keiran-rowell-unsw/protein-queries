"""Choosing which structure to draw, and refusing to draw a misleading one.

The report renders one image, so the choice is a claim. Ranking PDBe mappings by resolution
alone is the obvious rule and the wrong one: high-resolution entries are frequently
fabodies, nanobodies, Fc fusions, or isolated domains, and their excellent geometry belongs
to a fragment of the protein rather than to the protein. A live check on APP1 P05067 (770
residues) returned 8T82 at 1.1 Å with 0.01 coverage and 5ONR at 1.39 Å with 0.0 coverage.
Rendering either as "the structure" would show a designer a piece of the protein while the
label implied the whole of it — precisely the failure this module exists to prevent.

Resolution is therefore a tiebreaker applied *within* a coverage tier, never the primary
sort. The tiers are decided against the sequence the user is actually designing, whose
length is a required argument, and every returned candidate carries the reason it was
chosen so the report can show it next to the image.
"""

from __future__ import annotations

from typing import Any

# Fraction of the analysed sequence that must be covered before an entry is treated as
# depicting it. Below this the entry is a construct, and the report says so.
FULL_DISPLAY_THRESHOLD = 0.85
CONSTRUCT_DISPLAY_THRESHOLD = 0.30
MIN_VALID_RESOLUTION = 0.1  # Å; cryo-EM below this and NMR report resolution as null
MAX_PLAUSIBLE_RESOLUTION = 999.0  # Å; a nonsense resolution must not win a ranking

STATUS_LABELS: dict[str, str] = {
    "full": "covers the analysed sequence",
    "partial": "construct/domain coverage — not the full protein",
    "fragment_only": "no mapping covers a usable fraction; nothing is rendered",
    "none": "no PDBe mapping is available for this accession",
}


def _resolution(value: Any) -> float | None:
    """Parse a resolution in Å, rejecting nonsense without touching low coverage values.

    The guard applies to resolution only. The same bound must never be used for coverage:
    a 0.01-coverage fragment is precisely the case this module exists to catch, and
    filtering it out would let the fragment fall through as "coverage unreported".
    """
    number = _as_float(value)
    if number is None or number < MIN_VALID_RESOLUTION:
        return None
    return number if number <= MAX_PLAUSIBLE_RESOLUTION else None


def _coverage(value: Any) -> float | None:
    """Parse a coverage fraction; low values are meaningful and are kept."""
    number = _as_float(value)
    if number is None:
        return None
    return min(max(number, 0.0), 1.0)


def _as_float(value: Any) -> float | None:
    """Parse a finite non-negative number, applying no domain bounds.

    Deliberately value-blind: the validity limits for resolution and coverage are
    different and belong in their own parsers. A 0.01 coverage fraction is a real and
    decisive measurement, whereas a 0.01 Å resolution is not a measurement at all.
    """
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number or number in (float("inf"), float("-inf")) or number < 0:
        return None
    return number


def _coverage_from_regions(structure: dict[str, Any], length: int) -> float | None:
    """Derive coverage from observed UniProt segments when the API's fraction is absent.

    ``observed_regions`` is the residue-level SIFTS mapping, so it states what is actually
    modelled. It is preferred over the reported fraction because a construct can carry a
    generous ``coverage`` value while leaving the residues the user cares about unmodelled.
    """
    regions = structure.get("observed_regions") or []
    covered: set[int] = set()
    for region in regions:
        if not isinstance(region, dict):
            continue
        start, end = region.get("unp_start"), region.get("unp_end")
        if isinstance(start, int) and isinstance(end, int) and end >= start >= 1:
            covered.update(range(start, min(end, length) + 1))
    if not covered:
        return None
    return len(covered) / length if length else None


def effective_coverage(structure: dict[str, Any], length: int) -> float | None:
    """Best available estimate of the analysed sequence this entry depicts."""
    derived = _coverage_from_regions(structure, length)
    reported = _coverage(structure.get("coverage_fraction"))
    if derived is not None:
        # Keep the tighter of the two: a reported fraction that exceeds the observed
        # segments is describing the construct, not the modelled residues.
        if reported is None:
            return derived
        return min(derived, reported)
    if reported is not None:
        return reported
    # No fraction and no SIFTS segments: fall back to the declared mapping range so a
    # range-only record is tiered on its actual span instead of escaping as "unknown".
    return _range_coverage(structure, length)


def coverage_basis(structure: dict[str, Any], length: int) -> str:
    if structure.get("observed_regions"):
        return "PDBe SIFTS observed segments"
    if _coverage(structure.get("coverage_fraction")) is not None:
        return "PDBe reported coverage fraction"
    start, end = structure.get("uniprot_start"), structure.get("uniprot_end")
    if isinstance(start, int) and isinstance(end, int) and length:
        return "PDBe UniProt mapping range"
    return "no coverage field returned"


def _range_coverage(structure: dict[str, Any], length: int) -> float | None:
    """Coverage implied by the UniProt mapping range when no fraction or segments exist."""
    start, end = structure.get("uniprot_start"), structure.get("uniprot_end")
    if not (isinstance(start, int) and isinstance(end, int)) or length < 1:
        return None
    lo, hi = max(1, start), min(end, length)
    return (hi - lo + 1) / length if hi >= lo else None


def is_crystallographic(structure: dict[str, Any]) -> bool:
    method = str(structure.get("experimental_method") or "").lower()
    return "x-ray" in method or "diffraction" in method


def _sort_key(structure: dict[str, Any], length: int) -> tuple[float, float, str]:
    coverage = effective_coverage(structure, length)
    resolution = _resolution(structure.get("resolution_angstrom"))
    return (
        -(coverage if coverage is not None else -1.0),  # most coverage first, unknown last
        resolution if resolution is not None else MAX_PLAUSIBLE_RESOLUTION,
        str(structure.get("pdb_id") or ""),
    )


def choose_display_structure(
    structures: list[dict[str, Any]],
    *,
    sequence_length: int,
    preferred_pdb_id: str | None = None,
) -> dict[str, Any]:
    """Rank PDBe mappings for display and state why the winner won.

    ``sequence_length`` is the length of the sequence the user submitted, not of any
    PDB chain: the threshold question is how much of *their* protein an entry depicts.
    A missing coverage estimate disqualifies an entry from the ``full`` tier rather than
    letting it win by default.
    """
    if sequence_length < 1:
        raise ValueError("sequence_length must be a positive residue count")
    candidates = [item for item in structures if isinstance(item, dict) and item.get("pdb_id")]
    if not candidates:
        return {
            "status": "none",
            "selected": None,
            "alternative_count": 0,
            "reason": "PDBe returned no structures mapped to this accession.",
            "thresholds": {
                "full": FULL_DISPLAY_THRESHOLD,
                "construct": CONSTRUCT_DISPLAY_THRESHOLD,
            },
        }
    ranked = sorted(candidates, key=lambda item: _sort_key(item, sequence_length))

    if preferred_pdb_id:
        wanted = preferred_pdb_id.strip().upper()
        chosen = next((item for item in ranked if str(item.get("pdb_id")).upper() == wanted), None)
        if chosen is not None:
            return _verdict(
                chosen,
                chosen,
                ranked,
                sequence_length,
                basis="user-requested entry",
                rank_note="highest-ranked available mapping",
            )

    top = ranked[0]
    return _verdict(top, top, ranked, sequence_length, rank_note="highest-ranked available mapping")


def _verdict(
    chosen: dict[str, Any],
    ranked_top: dict[str, Any],
    ranked: list[dict[str, Any]],
    sequence_length: int,
    *,
    basis: str | None = None,
    rank_note: str = "",
) -> dict[str, Any]:
    """One guard, one place: the fragment threshold applies whether or not a caller named
    the entry, so a requested PDB ID cannot bypass it by being requested."""
    coverage = effective_coverage(chosen, sequence_length)
    thresholds = {"full": FULL_DISPLAY_THRESHOLD, "construct": CONSTRUCT_DISPLAY_THRESHOLD}
    prefix = f"{basis}; " if basis else ""
    if coverage is None:
        return {
            "status": "partial",
            "selected": chosen,
            "alternative_count": len(ranked) - 1,
            "coverage_estimate": None,
            "coverage_basis": coverage_basis(chosen, sequence_length),
            "selection_reason": (
                prefix + "PDBe returned no coverage estimate for this mapping, so coverage "
                "could not be verified; shown with coverage unconfirmed rather than asserted "
                "as a full model"
            ),
            "thresholds": thresholds,
        }
    tier = _tier_for(coverage)
    start, end = chosen.get("uniprot_start"), chosen.get("uniprot_end")
    span_text = (
        f" UniProt residues {start}-{end}."
        if isinstance(start, int) and isinstance(end, int)
        else ""
    )
    if tier == "full":
        reason = (
            prefix + f"{rank_note}: highest-resolution mapping covering at least "
            f"{FULL_DISPLAY_THRESHOLD:.0%} of the analysed {sequence_length}-residue sequence"
            f" (est. {coverage:.0%} coverage).{span_text}"
        )
    elif tier == "partial":
        reason = (
            prefix + f"{rank_note}: best available coverage is {coverage:.0%} of the analysed "
            f"{sequence_length}-residue sequence; this is a construct or domain, not the "
            f"full protein.{span_text}"
        )
    else:
        best = ", ".join(
            f"{item.get('pdb_id')} {_pct(effective_coverage(item, sequence_length))}"
            for item in ranked[:3]
        )
        reason = (
            prefix + f"every PDBe mapping covers under "
            f"{CONSTRUCT_DISPLAY_THRESHOLD:.0%} of the analysed sequence (best: {best}); "
            "no structure is rendered rather than present a fragment as the protein"
        )
    return {
        "status": tier,
        "selected": chosen if tier != "fragment_only" else None,
        "alternative_count": len(ranked) - 1,
        "coverage_estimate": round(coverage, 4),
        "coverage_basis": coverage_basis(chosen, sequence_length),
        "selection_reason": reason,
        "thresholds": thresholds,
    }


def _tier_for(coverage: float | None) -> str:
    if coverage is None:
        return "partial"
    if coverage >= FULL_DISPLAY_THRESHOLD:
        return "full"
    if coverage >= CONSTRUCT_DISPLAY_THRESHOLD:
        return "partial"
    return "fragment_only"


def _pct(coverage: float | None) -> str:
    return "coverage unreported" if coverage is None else f"{coverage:.0%}"


__all__ = [
    "CONSTRUCT_DISPLAY_THRESHOLD",
    "FULL_DISPLAY_THRESHOLD",
    "STATUS_LABELS",
    "choose_display_structure",
    "coverage_basis",
    "effective_coverage",
    "is_crystallographic",
]
