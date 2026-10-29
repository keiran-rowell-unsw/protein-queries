"""Optional metapredict V3 sequence-based disorder scores and plot.

Predicted disorder is not an experimental structure or a pLDDT measurement. Missing or
broken optional dependencies are reported explicitly; no default scores are fabricated.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

MAX_PREDICTION_LENGTH = 2000  # keep a synchronous OOD request within the SSH time limit
MIN_PREDICTION_LENGTH = 25  # V3 domain smoothing window is 2 * minimum_IDR_size = 24


def predict_disorder(sequence: str, out_path: Path) -> dict[str, Any]:
    """Run installed metapredict V3 once; include per-residue scores and IDR boundaries."""
    if len(sequence) < MIN_PREDICTION_LENGTH:
        return {
            "status": "skipped",
            "reason": (
                f"Sequence is shorter than {MIN_PREDICTION_LENGTH} residues; "
                "too short for this disorder-region prediction."
            ),
        }
    if len(sequence) > MAX_PREDICTION_LENGTH:
        return {
            "status": "skipped",
            "reason": (
                f"Sequence exceeds {MAX_PREDICTION_LENGTH} residues; "
                "no disorder prediction was run."
            ),
        }
    try:
        import metapredict
    except (ImportError, OSError):
        return {
            "status": "unavailable",
            "reason": "Disorder predictor is not installed on this node.",
        }

    try:
        result = metapredict.predict_disorder(
            sequence, version="V3", round_values=False, return_domains=True
        )
        values = [float(value) for value in result.disorder]
        boundaries = [
            {"start": int(start) + 1, "end": int(end)}
            for start, end in result.disordered_domain_boundaries
        ]
        if (
            len(values) != len(sequence)
            or any(not math.isfinite(value) or not 0 <= value <= 1 for value in values)
            or any(not 1 <= region["start"] <= region["end"] <= len(sequence)
                   for region in boundaries)
        ):
            raise ValueError("prediction shape, scores, or boundaries are invalid")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(9.0, 3.4))
        ax.plot(range(1, len(values) + 1), values, color="#4c1688", lw=1.2)
        for region in boundaries:
            ax.axvspan(region["start"], region["end"], color="#4c1688", alpha=0.12)
        ax.axhline(0.5, color="#555", lw=0.9, ls="--", label="V3 threshold (0.5)")
        ax.set(xlim=(1, max(2, len(values))), ylim=(0, 1),
               xlabel="residue position (1-based)", ylabel="predicted disorder score")
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    except (RuntimeError, ValueError, ImportError, OSError, TypeError, AttributeError) as exc:
        return {
            "status": "unavailable",
            "reason": f"Disorder prediction failed ({type(exc).__name__}).",
        }

    return {
        "status": "predicted",
        "predictor": "metapredict V3",
        "threshold": 0.5,
        "disordered_regions": boundaries,
        "mean_score": round(sum(values) / len(values), 3),
        "per_residue_scores": [round(value, 4) for value in values],
        "image": out_path.name,
        "caveat": (
            "Sequence-based prediction only; regions are model calls, not experimentally "
            "observed disorder or a structural confidence score."
        ),
    }
