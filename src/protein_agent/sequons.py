"""Predicted sequence motifs for design, including N-glycosylation sequons.

This module computes nothing that a database asserted. Its entire output class is
``motif``: pattern matches on the submitted sequence, with no occupancy, no evidence, and
no claim that the modification occurs. It exists because a designed or novel sequence has
no UniProt entry to annotate it, so annotation-only reporting is silent exactly where a
designer needs the most information.

Every predicted site is reconciled against the annotated feature list in
``compare_to_annotations`` so the two classes can never be conflated: a sequon that UniProt
also annotates is ``confirmed``, an unannotated one is ``predicted_only``, and — the case
that most deserves surfacing — an annotated site whose sequon is absent from this exact
sequence is ``annotated_without_sequon``, which usually means isoform, construct, or
numbering mismatch rather than a mysterious modification mechanism.

Cys-Gly (thioester) cores are deliberately *excluded*: they are not an N-glycan sequon and
would contaminate the track with a different chemistry.
"""

from __future__ import annotations

import re
from typing import Any

CLASSIFICATION = "heuristic / computational — never report as curated or measured"

# The canonical N-linked sequon: Asn, any residue except Pro, then Ser/Thr.
# The X position may itself be N, C, G, S, or T; only Pro is excluded, because it
# cannot adopt the backbone conformation the oligotransferase requires.
N_SEQUON = re.compile(r"(?=N[^P][ST])")
N_GLYCAN_RESIDUES = "NST"

# Deamidation is a liability, not a modification site: N-G is the fast, unstructured-case
# motif, N-S and N-T are the slower forms that nonetheless matter in a stored construct.
DEAMIDATION = re.compile(r"(?=NG)")

# Cleavage frequency alone is the liability signal; whether a site is surface-accessible
# needs a structure, which this module never claims to have.
TRYPSIN = re.compile(r"(?=[KR])(?!P)")
CHYMOTRYPSIN = re.compile(r"(?=[FYW])(?!P)")

MAX_FEATURE_LABELS = 60


def _residue_type(sequence: str, position: int) -> str:
    """Residue identity at a 1-based position; '?' when the position is off-sequence."""
    if 1 <= position <= len(sequence):
        return sequence[position - 1]
    return "?"


def _match_context(sequence: str, start: int, motif: str) -> str:
    """Human-readable sequon spelling, e.g. ``N542-T-S`` for a 1-based match start."""
    if motif == "N-X-[ST]":
        return (
            f"{_residue_type(sequence, start)}{start}-"
            f"{_residue_type(sequence, start + 1)}-"
            f"{_residue_type(sequence, start + 2)}"
        )
    return f"{_residue_type(sequence, start)}{start}-{_residue_type(sequence, start + 1)}"


def n_glycan_sequons(sequence: str) -> list[dict[str, Any]]:
    """Every N-X-[ST] sequon in the exact submitted sequence, 1-based at the Asn."""
    return [
        {
            "position": match.start() + 1,
            "motif": "N-X-[ST]",
            "sequon": _match_context(sequence, match.start() + 1, "N-X-[ST]"),
            "classification": "motif",
        }
        for match in N_SEQUON.finditer(sequence)
    ]


def deamidation_sites(sequence: str) -> list[dict[str, Any]]:
    """N-G (fast) and N-S/N-T (slow) deamidation liability positions.

    The slow class is matched on N-[ST] specifically, so the motif label always names the
    residue that actually satisfied the rule. Labelling an N-N neighbour as a slow site
    would advertise a deamidation liability the pattern does not support.
    """
    fast = [
        {
            "position": match.start() + 1,
            "motif": "N-G",
            "sequon": _match_context(sequence, match.start() + 1, "N-G"),
            "rate_class": "fast",
            "classification": "motif",
        }
        for match in DEAMIDATION.finditer(sequence)
    ]
    slow = [
        {
            "position": asparagine + 1,
            "motif": f"N-{sequence[acceptor]}",
            "sequon": f"N{asparagine + 1}-{sequence[acceptor]}",
            "rate_class": "slow",
            "classification": "motif",
        }
        for asparagine, acceptor in _slow_deamidation_pairs(sequence)
    ]
    return sorted([*fast, *slow], key=lambda item: (item["position"], item["rate_class"]))


def _slow_deamidation_pairs(sequence: str) -> list[tuple[int, int]]:
    """0-based (Asn, acceptor) index pairs where the Asn sits at +1 of a Ser or Thr."""
    return [
        (index - 1, index)
        for index in range(1, len(sequence))
        if sequence[index] in "ST" and sequence[index - 1] == "N"
    ]


def cleavage_sites(sequence: str, enzyme: str = "trypsin") -> dict[str, Any]:
    """Likely cleavage positions; frequency only, with no claim of accessibility."""
    pattern = {"trypsin": TRYPSIN, "chymotrypsin": CHYMOTRYPSIN}.get(enzyme)
    if pattern is None:
        raise ValueError(f"unsupported protease {enzyme!r}; use trypsin or chymotrypsin")
    positions = sorted({match.start() + 1 for match in pattern.finditer(sequence)})
    return {
        "enzyme": enzyme,
        "specificity": (
            "K/R not followed by P" if enzyme == "trypsin" else "F/Y/W not followed by P"
        ),
        "positions": positions,
        "count": len(positions),
        "classification": "motif",
        "caveat": (
            "sequence-level specificity only; cleavage also requires the site to be "
            "surface-accessible and unliganded, which cannot be established from sequence"
        ),
    }


def termini_flags(sequence: str) -> dict[str, Any]:
    """Terminal residues, the practical input to tag-site choice."""
    return {
        "n_term": _residue_type(sequence, 1),
        "c_term": _residue_type(sequence, len(sequence)) if sequence else "?",
        "n_term_proline": _residue_type(sequence, 1) == "P",
        "classification": "computed",
        "caveat": (
            "terminal residues only; tag tolerance at either end is not inferable from "
            "sequence and needs an experimental or predicted structure to assess"
        ),
    }


def compare_to_annotations(
    predicted: list[dict[str, Any]],
    annotated_glycosylation: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Label each predicted sequon by whether the fetched entry also annotates it.

    ``annotated_glycosylation`` must already come from a UniProt response in this run; this
    function reads its positions and never supplies new ones.
    """
    annotated_at: dict[int, list[str]] = {}
    for feature in annotated_glycosylation:
        start = feature.get("start")
        if isinstance(start, int):
            kinds = [str(feature.get("type"))]
            if feature.get("description"):
                kinds.append(str(feature["description"]))
            annotated_at.setdefault(start, []).extend(kinds)
    reconciled: list[dict[str, Any]] = []
    for site in predicted:
        position = site["position"]
        if position in annotated_at:
            state = "confirmed"
            detail = "annotated in this UniProt entry"
        elif any(position <= known <= position + 2 for known in annotated_at):
            state = "confirmed"
            detail = "annotated within this sequon span"
        else:
            state = "predicted_only"
            detail = "sequon present; occupancy not annotated in this entry"
        reconciled.append({**site, "annotation_state": state, "annotation_detail": detail})
    covered = {site["position"] for site in predicted} | {
        site["position"] + 2 for site in predicted
    }
    for position in sorted(annotated_at):
        if position not in covered and not any(
            position - 2 <= site["position"] <= position for site in predicted
        ):
            reconciled.append(
                {
                    "position": position,
                    "motif": "annotated",
                    "sequon": f"{_residue_type('', position)}{position}",
                    "classification": "curated",
                    "annotation_state": "annotated_without_sequon",
                    "annotation_detail": (
                        "UniProt annotates this glycosylation site but no N-X-[ST] sequon was "
                        "found at it in the analysed sequence; check isoform, construct, "
                        "tagging, and numbering before treating either as fact"
                    ),
                }
            )
    return sorted(reconciled, key=lambda item: (item["position"], item["classification"]))


def predicted_motifs(
    sequence: str, annotated_glycosylation: list[dict] | None = None
) -> dict[str, Any]:
    """Full motif inventory for one sequence, self-describing about what it cannot say."""
    sequons = n_glycan_sequons(sequence)
    annotated = annotated_glycosylation or []
    reconciled = compare_to_annotations(sequons, annotated)
    counts = {
        state: sum(item["annotation_state"] == state for item in reconciled)
        for state in ("confirmed", "predicted_only", "annotated_without_sequon")
    }
    if sequons:
        verdict = (
            f"{len(sequons)} N-glycan sequon(s) present at position(s) "
            f"{', '.join(str(item['position']) for item in sequons[:MAX_FEATURE_LABELS])}"
            + (
                f" (+{len(sequons) - MAX_FEATURE_LABELS} more)"
                if len(sequons) > MAX_FEATURE_LABELS
                else ""
            )
            + f". {counts['confirmed']} confirmed by annotation, "
            f"{counts['predicted_only']} predicted only. Presence of a sequon does not "
            "establish that the site is glycosylated."
        )
    else:
        verdict = (
            "No N-X-[ST] sequon was found in the analysed sequence. Absence of a sequon is "
            "NOT evidence that the protein carries no N-glycan: C-linked and non-canonical "
            "occupancy exist, and this scan covers neither."
        )
    return {
        "status": "computed" if sequence else "skipped",
        "classification": "motif",
        "verdict": verdict,
        "n_glycan_sequons": reconciled,
        "raw_sequon_positions": [item["position"] for item in sequons],
        "counts": counts,
        "deamidation": deamidation_sites(sequence),
        "cleavage": {
            enzyme: cleavage_sites(sequence, enzyme) for enzyme in ("trypsin", "chymotrypsin")
        },
        "termini": termini_flags(sequence),
        "method": {
            "n_glycan_sequon": "N-X-[ST], X != P (Barb & Kroc / Green-Yedney core sequon)",
            "deamidation": "N-G fast, N-S/N-T slow (Wright 1991 liability classes); "
            "the fast class is Gly-specific, not Ser-inclusive",
            "cleavage": "K/R and F/Y/W not followed by P",
        },
        "known_limitations": [
            "sequon occupancy depends on local flexibility and surface exposure, neither of "
            "which a linear pattern can see",
            "the Asn1 residue must be in the correct backbone conformation; a sequon in a "
            "buried position is unlikely to be occupied",
            "O-linked sites have no reliable linear consensus and are therefore not predicted "
            "here at all",
            "positions are relative to the analysed sequence; if that is an isoform or "
            "construct they will not match UniProt numbering",
        ],
        "meaning": (
            "motif matches nominate sites for engineering or for experimental confirmation; "
            "they do not establish that any modification occurs"
        ),
        "classification_note": CLASSIFICATION,
    }


__all__ = [
    "CLASSIFICATION",
    "cleavage_sites",
    "compare_to_annotations",
    "deamidation_sites",
    "n_glycan_sequons",
    "predicted_motifs",
    "termini_flags",
]
