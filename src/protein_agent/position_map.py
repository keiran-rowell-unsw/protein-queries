"""At-a-glance residue position map: the report's centrepiece.

A single 1-based axis carrying one track per annotation class, so locations read as
*where* rather than as a list of rows. This module draws and annotates only; it computes
no property and asserts no fact. Every coordinate it plots is passed in by the caller
from either a fetched UniProt feature (``curated``) or a pattern match on the analysed
sequence (``motif``).

Provenance is encoded in geometry, not just colour, so the distinction survives greyscale
printing and colour-vision deficiency: **filled markers are annotated in the retrieved
entry**, **hollow markers are motif predictions with no annotation**. The legend states it
in words as well, and predicted O-linked sites are never drawn because the sequon scan
does not attempt them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless node: no display, write files only
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Track grouping lives in a matplotlib-free module so the Passenger front end can build
# predicted tracks without a scientific environment; this module re-exports for callers
# that only know the plotting module.
from .position_tracks import (
    MAX_LABELLED_POSITIONS,
    MAX_RENDERED_PER_TRACK,
    SPAN_TRACKS,
    TRACK_COLORS,
    TRACK_ORDER,
    build_tracks,
)

__all__ = [
    "MAX_LABELLED_POSITIONS",
    "MAX_RENDERED_PER_TRACK",
    "SPAN_TRACKS",
    "TRACK_COLORS",
    "TRACK_ORDER",
    "build_tracks",
    "plot_position_map",
]

def plot_position_map(
    length: int,
    tracks: dict[str, list[tuple[int, int, str]]],
    *,
    label: str = "sequence",
    out_path: str | Path,
    dpi: int = 150,
) -> dict[str, Any]:
    """Render the position map; returns what was drawn plus its provenance caveats."""
    if length < 1:
        return {"status": "skipped", "reason": "No residues to plot."}
    active = [name for name in TRACK_ORDER if tracks.get(name)]
    if not active:
        return {
            "status": "empty",
            "reason": (
                "No annotated features and no predicted sequons were available for this "
                "sequence; the map is intentionally blank rather than assumed."
            ),
        }

    height = max(2.6, 0.52 * len(active) + 1.5)
    fig, ax = plt.subplots(figsize=(9.6, height))
    for row, name in enumerate(active):
        y = len(active) - row
        color = TRACK_COLORS.get(name, "#334155")
        predicted = "predicted" in name
        for start, end, _detail in tracks[name]:
            if name == "Disulfide bond":
                # UniProt encodes a disulfide as the two bonded cysteines in start/end,
                # NOT as a covered range. Drawing that as a bar would imply a domain
                # spanning e.g. 149 residues, which is the wrong claim entirely.
                for position in {start, end}:
                    ax.plot(
                        [position],
                        [y],
                        marker="|",
                        markersize=13,
                        color=color,
                        markeredgewidth=2.2,
                        linestyle="none",
                        zorder=3,
                    )
                if start != end:
                    ax.plot(
                        [start, end],
                        [y, y],
                        color=color,
                        lw=0.7,
                        alpha=0.45,
                        linestyle=":",
                        zorder=1,
                    )
                continue
            if start == end or name in {"Deamidation risk", "N-glycan sequon (predicted)"}:
                ax.plot(
                    [start],
                    [y],
                    marker="o",
                    markersize=7.5,
                    markerfacecolor="white" if predicted else color,
                    markeredgecolor=color,
                    markeredgewidth=1.6,
                    linestyle="none",
                    zorder=3,
                )
            elif end - start < 4:
                # A 1-3 residue feature is a marker, not a band: at 1/770 of the axis a
                # rectangle is an invisible hairline that reads as "nothing annotated".
                for position in {start, end}:
                    ax.plot(
                        [position],
                        [y],
                        marker="o",
                        markersize=7.5,
                        markerfacecolor="white" if predicted else color,
                        markeredgecolor=color,
                        markeredgewidth=1.6,
                        linestyle="none",
                        zorder=3,
                    )
            else:
                ax.add_patch(
                    Rectangle(
                        (start, y - 0.2),
                        max(1, end - start),
                        0.4,
                        facecolor=color,
                        edgecolor=color,
                        alpha=0.35,
                        zorder=2,
                    )
                )
        ax.text(
            -length * 0.012,
            y,
            name,
            ha="right",
            va="center",
            fontsize=7.8,
            color="#1f2937",
            clip_on=False,
        )

    ax.axhline(0.45, color="#94a3b8", lw=0.8, zorder=1)
    ax.set_xlim(0, length)
    ax.set_ylim(0.2, len(active) + 0.8)
    ax.set_yticks([])
    #ax.set_xlabel("residue position", fontsize=9)
    ax.tick_params(labelsize=8)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        f"{label} — feature positions ({length} residues)",
        fontsize=10.5,
        loc="left",
        pad=10,
    )

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=7.5,
            markerfacecolor="#0e7490",
            markeredgecolor="#0e7490",
            label="annotated in the retrieved UniProt entry",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=7.5,
            markerfacecolor="white",
            markeredgecolor="#0e7490",
            label="motif prediction — occupancy not established",
        ),
        Rectangle(
            (0, 0),
            1,
            1,
            facecolor="#475569",
            edgecolor="#475569",
            alpha=0.35,
            label="span (domain / topological region)",
        ),
        Line2D(
            [],
            [],
            marker="|",
            linestyle="none",
            markersize=11,
            color="#b45309",
            markeredgewidth=2.2,
            label="disulfide: paired tick marks join two bonded cysteines",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.16 / max(1, height / 3)),
        fontsize=7.5,
        frameon=False,
        ncol=2,
        columnspacing=1.4,
        handletextpad=0.5,
    )
    # The filled/hollow distinction is carried by the legend entries themselves. A second
    # free-text footnote at the same sub-axes offset collided with the last legend row,
    # so the caveat lives in the returned payload and the page instead of on the image.
    fig.tight_layout()

    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    counts = {name: len(tracks[name]) for name in active}
    return {
        "status": "rendered",
        "image": destination.name,
        "length_residues": length,
        "tracks": counts,
        "labelled_positions": {
            name: sorted({start for start, _end, _ in tracks[name]})[:MAX_LABELLED_POSITIONS]
            for name in active
        },
        "provenance": {
            "annotated_tracks": "UniProt feature list retrieved in this run",
            "predicted_tracks": "sequon/deamidation pattern match on the analysed sequence",
        },
        "caveats": [
            "hollow markers are sequence-pattern predictions and do not establish that any "
            "modification occurs",
            "positions refer to the analysed sequence, which may be an isoform or construct "
            "rather than the canonical entry",
            "absence of a marker is not evidence that a modification is absent",
        ],
    }
