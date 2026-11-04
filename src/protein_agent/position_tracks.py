"""Track grouping for residue position maps — matplotlib-free by design.

Everything here is evidence grouping: coordinates and labels passed in by the caller,
sorted into tracks whose provenance is encoded in geometry at render time (filled =
annotated in a retrieved entry, hollow = motif prediction). No property is computed and
no fact is asserted, and nothing in this module may import plotting libraries, so the
Passenger front end can build predicted tracks without a scientific environment.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

# Track order is the priority a designer reads in: glycosylation first, then other PTMs,
# then the features that bound a construct.
TRACK_ORDER: tuple[str, ...] = (
    "Glycosylation (annotated)",
    "N-glycan sequon (predicted)",
    "Other PTM",
    "Disulfide bond",
    "Active / binding site",
    "Transmembrane / topology",
    "Domain",
    "Deamidation risk",
    "Variants",
)

TRACK_COLORS: dict[str, str] = {
    "Glycosylation (annotated)": "#0e7490",
    "N-glycan sequon (predicted)": "#0e7490",
    "Other PTM": "#7c3aed",
    "Disulfide bond": "#b45309",
    "Active / binding site": "#dc2626",
    "Transmembrane / topology": "#2563eb",
    "Domain": "#475569",
    "Deamidation risk": "#ca8a04",
    "Variants": "#16a34a",
}

SPAN_TRACKS = {"Transmembrane / topology", "Domain", "Variants"}
MAX_LABELLED_POSITIONS = 24
MAX_RENDERED_PER_TRACK = 400


def _positions(features: Sequence[dict[str, Any]], types: set[str]) -> list[tuple[int, int, str]]:
    """Normalise fetched/derived features to (start, end, detail) tuples."""
    out: list[tuple[int, int, str]] = []
    for feature in features:
        if feature.get("type") not in types:
            continue
        start, end = feature.get("start"), feature.get("end")
        if not isinstance(start, int):
            continue
        stop = end if isinstance(end, int) else start
        detail = str(feature.get("description") or feature.get("type") or "")
        out.append((start, max(start, stop), detail))
    return out


def build_tracks(
    *,
    uniprot_features: Sequence[dict[str, Any]] | None,
    sequons: Sequence[dict[str, Any]] | None,
    deamidation: Sequence[dict[str, Any]] | None,
) -> dict[str, list[tuple[int, int, str]]]:
    """Group the caller's evidence into plottable tracks, keeping provenance separate.

    ``uniprot_features`` must come from a UniProt response in this run. Predicted entries
    are only ever drawn on their own tracks, never merged onto an annotated one.
    """
    features = list(uniprot_features or [])
    annotated_glycan = {
        site["position"] for site in (sequons or []) if site.get("annotation_state") == "confirmed"
    }
    tracks: dict[str, list[tuple[int, int, str]]] = {name: [] for name in TRACK_ORDER}

    for start, end, detail in _positions(features, {"Glycosylation"}):
        tracks["Glycosylation (annotated)"].append((start, end, detail))

    for site in sequons or []:
        position = site.get("position")
        if not isinstance(position, int):
            continue
        # A confirmed sequon is already drawn as an annotated site; plotting it twice would
        # read as two independent pieces of evidence, so predicted tracks stay predicted-only.
        if position in annotated_glycan or site.get("annotation_state") != "predicted_only":
            continue
        tracks["N-glycan sequon (predicted)"].append(
            (position, position, f"{site.get('sequon', '')} sequon; occupancy not annotated")
        )

    tracks["Other PTM"] = _positions(
        features,
        {"Modified residue", "Lipidation", "Cross-link", "Non-standard residue"},
    )
    tracks["Disulfide bond"] = _positions(features, {"Disulfide bond"})
    tracks["Active / binding site"] = _positions(
        features,
        {
            "Active site",
            "Binding site",
            "Site",
            "Metal binding",
            "Calcium binding",
            "Nucleotide binding",
            "DNA binding",
        },
    )
    tracks["Transmembrane / topology"] = _positions(
        features,
        {"Transmembrane", "Intramembrane", "Topological domain", "Signal", "Transit peptide"},
    )
    tracks["Domain"] = _positions(
        features, {"Domain", "Repeat", "Zinc finger", "Coiled coil", "Region", "Motif"}
    )
    tracks["Deamidation risk"] = [
        (item["position"], item["position"], str(item.get("sequon", "")))
        for item in (deamidation or [])
        if isinstance(item.get("position"), int)
    ]
    tracks["Variants"] = _positions(
        features,
        {"Natural variant", "Mutagenesis", "Sequence conflict"},
    )
    return {name: items[:MAX_RENDERED_PER_TRACK] for name, items in tracks.items()}


__all__ = [
    "MAX_LABELLED_POSITIONS",
    "MAX_RENDERED_PER_TRACK",
    "SPAN_TRACKS",
    "TRACK_COLORS",
    "TRACK_ORDER",
    "build_tracks",
]
