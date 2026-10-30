"""Deterministic structure images for the report: CA-trace cartoon and ligand map.

Chosen over a PyMOL-backed renderer because this deployment must produce an image on
every run. The node's glibc is 2.28 and the ``pymol-open-source`` wheels require 2.35, so
PyMOL cannot load here; an image that never appears is worse than a plainer image that
does. This module depends only on Biopython and matplotlib, both already required.

What the image is: a C-alpha trace from the deposited coordinates, drawn in two
orientations. It is not a raytraced cartoon, has no shaded surface, no depth cues, and no
side chains. It is not a prediction — every atom plotted comes from the entry PDBe served
for that PDB ID in this run.

What the image is not: a depiction of the analysed sequence. The deposited model carries
its own construct boundaries, engineered mutations, and truncations, and residues absent
from the coordinates are absent from the image rather than interpolated. Coverage is
computed here from *observed* CA atoms and must be shown beside the image, because the
deposited residue span is not the observed span.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless node: no display, write files only
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

MAX_COORDINATE_BYTES = 12_000_000  # bound a synchronous OOD request
MAX_PLOTTED_RESIDUES = 6000  # above this a CA trace is visual noise, not information
WATER_AND_ADDUCTS = {"HOH", "DOD", "WAT"}
BINARY_BOUNDARIES = (0, 28, 14, 29, 15, 30, 16, 31, 17, 32, 18, 33, 19, 34, 20, 35)
PROVENANCE = (
    "Biopython PDB/MMCIF parser + matplotlib 2D projection of C-alpha atoms; "
    "backbone trace only, not a raytraced cartoon"
)


def parse_coordinates(content: bytes, *, preferred: str = "pdb") -> dict[str, Any]:
    """Parse served coordinates into chains of observed CA trace plus ligand instances.

    Returns raw structural facts only; no interpretation, and no claim about the
    relationship between these atoms and any query sequence.
    """
    from Bio.PDB import MMCIFParser, PDBParser

    stripped = content.lstrip()[:80].decode("utf-8", "replace").upper()
    if preferred == "pdb" and stripped.startswith(("ATOM", "HETATM", "HEADER")):
        source, parser = "pdb", PDBParser(QUIET=True)
    elif stripped.startswith(("DATA_", "LOOP_", "_PDBDFLINE", "HEADER")) or preferred == "cif":
        source, parser = "cif", MMCIFParser(QUIET=True)
    else:
        source, parser = "pdb", PDBParser(QUIET=True)
    if source == "pdb" and _accepts_header_only(parser):
        # Waters carry no design information and inflate the residue count.
        structure = parser.get_structure("entry", io.StringIO(_as_text(content)), header_only=False)
    else:
        structure = parser.get_structure("entry", io.StringIO(_as_text(content)))
    model = next(iter(structure))

    chains: list[dict[str, Any]] = []
    ligands: dict[str, dict[str, Any]] = {}
    for chain in model:
        trace: list[tuple[int, float, float, float]] = []
        for residue in chain:
            hetflag, identifier, insertion = residue.get_id()
            if not isinstance(identifier, int):
                continue
            atom = residue.child_dict.get("CA")
            if hetflag == " " and insertion.strip() == "" and atom is not None:
                # Standard polymer residue: no hetflag, no insertion code, CA modelled.
                x, y, z = atom.get_coord()
                trace.append((identifier, float(x), float(y), float(z)))
            elif hetflag.strip() != "":
                # Bound ligand: a HETATM group that is not water. Classified from the
                # record itself, never from a residue-name lookup table.
                name = residue.get_resname().strip() or "UNLABELLED"
                if name.upper() in WATER_AND_ADDUCTS:
                    continue
                bucket = ligands.setdefault(name, {"name": name, "instances": 0, "chains": set()})
                bucket["instances"] += 1
                bucket["chains"].add(chain.id)
        if trace:
            chains.append(
                {
                    "chain_id": chain.id,
                    "ca_positions": [position for position, *_ in trace],
                    "coordinates": [(x, y, z) for _, x, y, z in trace],
                }
            )
    for bucket in ligands.values():
        bucket["chains"] = sorted(bucket["chains"])
    return {
        "parsed_as": source,
        "chains": chains,
        "ligands": sorted(ligands.values(), key=lambda item: item["name"]),
    }


def _as_text(content: bytes) -> str:
    """Decode served coordinates; Biopython's header parser wants str, not bytes."""
    if isinstance(content, str):
        return content
    return content.decode("utf-8", errors="replace")


def _accepts_header_only(parser: Any) -> bool:
    import inspect

    try:
        return "header_only" in inspect.signature(parser.get_structure).parameters
    except (TypeError, ValueError):
        return False


def observed_coverage(parsed: dict[str, Any]) -> dict[str, Any]:
    """Contiguous observed segments per chain, from CA atoms actually in the file."""
    summary: dict[str, Any] = {"chains": {}, "total_ca_atoms": 0}
    for chain in parsed["chains"]:
        positions = sorted(set(chain["ca_positions"]))
        segments: list[list[int]] = []
        start = previous = positions[0]
        for position in positions[1:]:
            if position != previous + 1:
                segments.append([start, previous])
                start = position
            previous = position
        segments.append([start, previous])
        summary["chains"][chain["chain_id"]] = {
            "observed_ca_count": len(positions),
            "deposited_span": [positions[0], positions[-1]],
            "observed_segments": segments,
            "fragmented": len(segments) > 1,
        }
        summary["total_ca_atoms"] += len(positions)
    return summary


def _autoscale(points: list[tuple[float, float]]) -> tuple[float, float]:
    """Square data bounds around a point cloud so equal aspect is never degenerate."""
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
    centre_x, centre_y = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    half = span / 2 + span * 0.05
    return (centre_x - half, centre_x + half), (centre_y - half, centre_y + half)


def _project(coordinates: list[tuple[float, float, float]], view: str) -> list[tuple[float, float]]:
    """Drop the viewing axis; orthonormal projections so relative distances stay true."""
    if view == "xy":
        return [(x, y) for x, y, _z in coordinates]
    if view == "xz":
        return [(x, z) for x, _y, z in coordinates]
    return [(y, z) for _x, y, z in coordinates]


def _breaks_at_gaps(points: list[tuple[float, float]], positions: list[int]) -> list[Any]:
    """Split the polyline wherever residue numbering jumps, so links are not fabricated."""
    segments: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = [points[0]]
    for index in range(1, len(points)):
        if positions[index] != positions[index - 1] + 1:
            if len(current) > 1:
                segments.append(current)
            current = [points[index]]
        else:
            current.append(points[index])
    if len(current) > 1:
        segments.append(current)
    return segments


def _label_chain(axis: Any, point: tuple[float, float], name: str, colour: str) -> None:
    """Place a chain label clear of the data, flipping into the frame at the edges.

    Equal-aspect autoscaling can push the first CA atom against the top of the axes, where
    an outward text offset is clipped by bbox_inches='tight'. A clipped chain label reads as
    an anonymous trace, so the offset is inverted whenever the point sits near a limit.
    """
    bounds = axis.get_xlim(), axis.get_ylim()
    x_lo, x_hi = bounds[0]
    y_lo, y_hi = bounds[1]
    right = point[0] > x_lo + (x_hi - x_lo) * 0.85
    below = point[1] > y_lo + (y_hi - y_lo) * 0.88
    axis.annotate(
        name,
        point,
        fontsize=7.5,
        color=colour,
        fontweight="bold",
        zorder=5,
        ha="right" if right else "left",
        xytext=(-4 if right else 2, -8 if below else 5),
        textcoords="offset points",
    )


def plot_cartoon(
    parsed: dict[str, Any],
    out_path: Path | str,
    *,
    views: tuple[str, ...] = ("xy", "xz"),
    title: str = "structure",
    dpi: int = 150,
    width: float = 9.6,
) -> dict[str, Any]:
    """Render CA traces side by side in the requested projections, labelled by chain."""
    if not parsed["chains"]:
        return {
            "status": "unavailable",
            "reason": "The parsed coordinate file contained no protein CA atoms to draw.",
        }
    total = sum(len(chain["ca_positions"]) for chain in parsed["chains"])
    if total > MAX_PLOTTED_RESIDUES:
        return {
            "status": "skipped",
            "reason": (
                f"{total} CA atoms exceed the {MAX_PLOTTED_RESIDUES}-atom plotting limit; "
                "a CA trace at that density is noise rather than information."
            ),
        }

    palette = ["#0e7490", "#b45309", "#7c3aed", "#166534", "#be123c", "#1d4ed8"]
    fig, axes = plt.subplots(1, len(views), figsize=(width, width * 0.44), squeeze=False)
    chain_summaries = observed_coverage(parsed)
    labels: list[str] = []
    for axis, view in zip(axes[0], views):
        all_points: list[tuple[float, float]] = []
        for index, chain in enumerate(parsed["chains"]):
            colour = palette[index % len(palette)]
            points = _project(chain["coordinates"], view)
            positions = chain["ca_positions"]
            all_points.extend(points)
            segments = _breaks_at_gaps(points, positions)
            if segments:
                axis.add_collection(
                    LineCollection(segments, colors=colour, linewidths=0.9, alpha=0.9, zorder=2)
                )
            axis.scatter(
                [point[0] for point in points],
                [point[1] for point in points],
                # CA markers carry the fold legibly; connecting lines alone read as
                # ambiguous spaghetti once a chain is compact and globular.
                s=5.5,
                color=colour,
                zorder=3,
                linewidths=0,
            )
            _label_chain(axis, points[0], chain["chain_id"], colour)
            stats = chain_summaries["chains"][chain["chain_id"]]
            labels.append(
                f"chain {chain['chain_id']}: {stats['observed_ca_count']} CA, observed "
                f"{stats['deposited_span'][0]}-{stats['deposited_span'][1]}"
                + (" (fragmented)" if stats["fragmented"] else "")
            )
        if not all_points:
            continue
        # Square, equal-aspect frame in Angstroms, deliberately unticked: numeric axis
        # labels invite misreading coordinate extents as residue numbering.
        (x_lo, x_hi), (y_lo, y_hi) = _autoscale(all_points)
        axis.set(xlim=(x_lo, x_hi), ylim=(y_lo, y_hi))
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(f"{view.upper()} projection", fontsize=9, loc="left")
        axis.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
        for spine in ("top", "right", "left", "bottom"):
            axis.spines[spine].set_visible(False)
        axis.grid(alpha=0.15, lw=0.5)
    fig.suptitle(title, fontsize=10.5, x=0.012, ha="left", y=0.99)
    # Captions are drawn inside the figure in space reserved by the axes rect, not placed
    # below y=0: text outside the canvas is clipped by savefig, and the earlier
    # bbox_inches='tight' workaround resized the two equal-aspect panels unevenly.
    caption_lines = []
    if labels:
        caption_lines.append(
            ("modelled: " + " | ".join(dict.fromkeys(labels)), "#334155", "normal")
        )
    caption_lines.append(
        (
            f"CA trace from deposited coordinates \u00b7 {total} residues modelled \u00b7 "
            "axes are in \u00c5, unlabelled by design \u00b7 not a raytrace",
            "#54585E",
            "italic",
        )
    )
    caption_lines.append(
        (
            "\u03b1/\u03b2 cartoons omitted: DSSP is not installed on this node, and secondary "
            "structure is not inferred from sequence here.",
            "#94a3b8",
            "italic",
        )
    )
    reserved = 0.055 + 0.075 * len(caption_lines)
    fig.tight_layout(rect=(0, reserved, 1, 0.95))
    for line_index, (line, colour, style) in enumerate(caption_lines):
        fig.text(
            0.012,
            reserved - 0.052 * (line_index + 1),
            line,
            fontsize=7.2,
            color=colour,
            style=style,
        )
    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=dpi)
    plt.close(fig)
    return {
        "status": "rendered",
        "image": destination.name,
        "views": list(views),
        "ca_atoms_plotted": total,
        "chains_drawn": sorted(chain_summaries["chains"]),
        "renderer": PROVENANCE,
        "caveats": [
            "CA trace only: no side chains, no surface, no depth shading",
            "modelled ligands and ions are counted and reported, but are not drawn into "
            "this backbone image; see the bound-molecules figure for what is present",
            "residues unmodelled in the deposited entry are absent from the image, not "
            "interpolated; polylines are broken at numbering gaps rather than bridged",
            "the entry may carry engineered mutations, tags, or truncations relative to the "
            "analysed sequence",
            "axes are Angstrom coordinates with a fixed equal aspect ratio, so relative "
            "distances are true but the ticks are intentionally suppressed",
            "no secondary-structure cartoons are drawn: DSSP is unavailable on this node and "
            "\u03b1/\u03b2 assignment is not inferred from sequence",
        ],
    }


def plot_ligand_map(
    parsed: dict[str, Any],
    coverage: dict[str, Any],
    out_path: Path | str,
    *,
    length: int,
    title: str = "bound molecules",
    dpi: int = 150,
) -> dict[str, Any]:
    """Ligand identity and multiplicity beside the observed range of each protein chain.

    Ligand instances are counted from the coordinates; UniProt cross-reference chains are
    deliberately not shown here, because that table covers the whole entry and can overstate
    what any single chain binds.
    """
    if not parsed["ligands"]:
        return {
            "status": "empty",
            "reason": "No non-water ligand is modelled in these coordinates.",
        }
    destination = Path(out_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    names = parsed["ligands"]
    fig, axis = plt.subplots(figsize=(9.6, max(1.5, 0.34 * len(names) + 1.2)))
    axis.barh(
        [item["name"] for item in names],
        [item["instances"] for item in names],
        color="#b45309",
        height=0.6,
    )
    axis.set_xlabel("modelled instances in this entry")
    axis.set_title(title, fontsize=10.5, loc="left", pad=8)
    axis.tick_params(labelsize=8)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    for axis_index, item in enumerate(names):
        axis.text(
            item["instances"] + 0.15,
            axis_index,
            "chain " + ",".join(item["chains"][:4]) + ("…" if len(item["chains"]) > 4 else ""),
            fontsize=7,
            color="#54585E",
            va="center",
        )
    observed_text = "; ".join(
        f"{chain}: {stats['observed_ca_count']} CA, "
        f"{stats['deposited_span'][0]}-{stats['deposited_span'][1]}"
        + (" (fragmented)" if stats["fragmented"] else "")
        for chain, stats in sorted(coverage["chains"].items())
    )
    fig.text(
        0.995,
        0.02,
        f"observed in coordinates: {observed_text}" if observed_text else "",
        ha="right",
        fontsize=6.8,
        color="#54585E",
        style="italic",
    )
    fig.tight_layout()
    fig.savefig(destination, dpi=dpi)
    plt.close(fig)
    return {
        "status": "rendered",
        "image": destination.name,
        "ligands": [
            {"name": item["name"], "instances": item["instances"], "chains": item["chains"]}
            for item in names
        ],
        "renderer": PROVENANCE,
        "caveats": [
            "counts are modelled instances in this one entry; stoichiometry is not implied",
            "waters are excluded",
            "alternative conformations and partial occupancies are counted as modelled",
        ],
    }


__all__ = [
    "PROVENANCE",
    "observed_coverage",
    "parse_coordinates",
    "plot_cartoon",
    "plot_ligand_map",
]
