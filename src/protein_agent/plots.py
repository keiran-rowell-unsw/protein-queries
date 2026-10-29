"""Hydropathy plot rendering.

Purely a visualisation of `properties.hydropathy_profile`: this module computes nothing
new. Every value drawn comes from the Kyte-Doolittle centred window in `properties.py`,
so the plot inherits the same provenance and must not be described as a structure or
disorder prediction.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless node: no display, write files only
import matplotlib.pyplot as plt

from .properties import clean_fasta_sequence, hydropathy_profile


def plot_hydropathy(
    sequence: str,
    *,
    label: str = "user sequence",
    window: int = 19,
    out_path: str | Path = "hydropathy.png",
    dpi: int = 150,
) -> dict[str, Any]:
    """Render a Kyte-Doolittle hydropathy trace to PNG; returns plot provenance.

    Never raises on odd input it cannot handle: unusable sequence is reported through
    `SequenceError` from `clean_fasta_sequence`, and the caller surfaces it.
    """
    cleaned, warnings = clean_fasta_sequence(sequence)
    profile = hydropathy_profile(cleaned, window=window)
    positions = [p for p, _ in profile]
    values = [v for _, v in profile]

    fig, ax = plt.subplots(figsize=(9.0, 3.4))
    ax.axhspan(0, max(values) + 0.2, color="#f4f2ec", zorder=0)
    ax.axhline(0, color="#54585E", lw=0.8, ls="--", zorder=1)
    ax.plot(positions, values, color="#3538CD", lw=1.6, zorder=3)
    ax.fill_between(
        positions, values, 0, where=[v >= 0 for v in values],
        color="#3538CD", alpha=0.18, zorder=2,
    )
    ax.fill_between(
        positions, values, 0, where=[v < 0 for v in values],
        color="#0057B8", alpha=0.14, zorder=2,
    )

    peak = max(profile, key=lambda t: t[1])
    trough = min(profile, key=lambda t: t[1])
    ax.annotate(
        f"most hydrophobic\n{peak[0]} ({peak[1]:+.2f})",
        xy=peak,
        xytext=(peak[0], min(max(peak[1] + 0.35, 0.4), max(values) + 0.7)),
        fontsize=7.5,
        ha="center",
        color="#0F0F0F",
    )
    ax.annotate(
        f"most hydrophilic\n{trough[0]} ({trough[1]:+.2f})",
        xy=trough,
        xytext=(trough[0], max(trough[1] - 0.55, min(values) - 0.9)),
        fontsize=7.5,
        ha="center",
        color="#0F0F0F",
    )

    ax.set_xlabel("residue position (1-based)", fontsize=9)
    ax.set_ylabel("hydropathy index", fontsize=9)
    ax.set_title(
        f"{label} — Kyte-Doolittle, {window}-residue centred window",
        fontsize=10.5,
        loc="left",
    )
    ax.set_xlim(min(positions), max(positions))
    ax.set_ylim(min(values) - 1.0, max(values) + 0.9)
    ax.tick_params(labelsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.text(
        1.0,
        -0.32,
        "computed from composition; not a structure or disorder prediction",
        transform=ax.transAxes,
        ha="right",
        fontsize=7,
        color="#54585E",
        style="italic",
    )
    fig.tight_layout()

    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return {
        "image": str(destination.resolve()),
        "length_residues": len(cleaned),
        "window": window,
        "scale": "Kyte-Doolittle (1982)",
        "most_hydrophobic": {"position": peak[0], "value": round(peak[1], 4)},
        "most_hydrophilic": {"position": trough[0], "value": round(trough[1], 4)},
        "y_range": [round(min(values), 4), round(max(values), 4)],
        "warnings": warnings,
        "interpretation": (
            "positive values are hydrophobic and negative values hydrophilic; a sustained "
            "positive stretch of ~19+ residues is consistent with a transmembrane or "
            "buried-hydrophobic segment, which is a hypothesis to test, not a conclusion"
        ),
        "caveats": [
            "the plot visualises computed values and adds no information of its own",
            "no experimental or predicted structure is implied",
        ],
    }


__all__ = ["plot_hydropathy"]
