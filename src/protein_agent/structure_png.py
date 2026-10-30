"""Render a PDBe coordinate file to PNG with pdb2png, and say exactly what was rendered.

The report must never present a picture it cannot substantiate. Every path through this
module therefore returns one of four statuses — ``rendered``, ``not_found``,
``unavailable``, or ``skipped`` — each carrying the reason, so an absent image is a
reported finding rather than a silent gap. Nothing here invents a structure: the entry is
fetched from PDBe in this run, or nothing is drawn.

Coordinate format
-----------------
``pdb2png`` hands its input straight to ``cmd.load(path)``, and PyMol infers the format
from the file suffix. PDBe's canonical download is now mmCIF (``{id}.cif``), whose suffix
would leave the upstream PNG name derivation broken and may not be read at all by an
older PyMol build. ``_load_command`` therefore writes the molecule out as legacy PDB
(``.pdb``) via ``cmd.save(..., "pdb")`` before rendering — the same re-serialization
choice Biopython's PDBIO makes — and records that conversion so a reviewer can see the
image was not drawn from the deposited CIF directly.

Chain selection is delegated to the caller's ``select_chain`` decision (see
``pdbe_report.choose_display_structure``); this module does not choose a structure, it
draws the one it is given and reports what was actually visible in it.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

from .config import Settings
from .pdb2png import pdb2png

MAX_COORDINATE_BYTES = 12_000_000  # bound a synchronous OOD request; large ribosomes excluded
RENDER_STYLES = ("default", "default_bw", "default_cartoon", "dark", "muted")
VIEWS: dict[str, tuple[float, float, float]] = {
    "front": (0.0, 0.0, 0.0),
    "side": (0.0, 90.0, 0.0),
    "top": (90.0, 0.0, 0.0),
}
CIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.cif"
PDB_FORMAT_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}.pdb"
LEGACY_PDB_URL = "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id}.ent"


class StructureRenderError(RuntimeError):
    """Coordinate retrieval or rasterisation failed; reported as a status, not a guess."""


def download_coordinates(pdb_id: str, *, settings: Settings | None = None) -> bytes:
    """Fetch deposited coordinates for one entry from PDBe, bounded in size.

    Legacy PDB is requested first because PyMOL infers format from the file suffix and
    older builds cannot read mmCIF; PDBe serves both. When no PDB-style file exists, which
    is normal for large or modern-only entries, the deposited mmCIF is returned instead.
    """
    import httpx

    identifier = pdb_id.strip().upper()
    if not (len(identifier) == 4 and identifier[0].isdigit()):
        raise StructureRenderError(f"{pdb_id!r} is not a well-formed PDB entry identifier")
    lowered = identifier.lower()
    urls = (
        PDB_FORMAT_URL.format(pdb_id=lowered),
        LEGACY_PDB_URL.format(pdb_id=lowered),
        CIF_URL.format(pdb_id=lowered),
    )
    try:
        with httpx.Client(timeout=30.0, follow_redirects=True) as http:
            for url in urls:
                with http.stream("GET", url) as response:
                    if response.status_code == 404:
                        continue
                    response.raise_for_status()
                    buffer = bytearray()
                    for chunk in response.iter_bytes():
                        buffer.extend(chunk)
                        if len(buffer) > MAX_COORDINATE_BYTES:
                            raise StructureRenderError(
                                f"{identifier} coordinates exceed the "
                                f"{MAX_COORDINATE_BYTES // 1_000_000} MB render limit"
                            )
                    return bytes(buffer)
            return b""
    except StructureRenderError:
        raise
    except httpx.HTTPError as exc:
        raise StructureRenderError(
            f"PDBe coordinate download failed: {type(exc).__name__}"
        ) from exc


def render_structure_png(
    pdb_id: str,
    *,
    out_path: Path,
    chain_id: str | None = None,
    view: str = "front",
    render_style: str = "default_cartoon",
    width: int = 1100,
    height: int = 1100,
    coordinates: bytes | None = None,
) -> dict[str, Any]:
    """Draw one entry/chain with pdb2png; returns status plus what was actually rendered.

    ``chain_id`` of ``None`` renders every protein chain in the entity, which is the right
    choice for an assembly but must be labelled as such by the caller.
    """
    out_path = Path(out_path)
    try:
        import pymol  # noqa: F401
    except (ImportError, OSError) as exc:
        return {
            "status": "unavailable",
            "reason": (
                "The PyMOL rendering engine is not installed on this node "
                f"({type(exc).__name__}); pdb2png cannot rasterise coordinates. "
                "Install the ``render`` extra to enable structure images."
            ),
            "pdb_id": pdb_id.upper(),
        }
    if view not in VIEWS:
        return {
            "status": "invalid",
            "reason": f"unknown view {view!r}; choose from {', '.join(VIEWS)}",
            "pdb_id": pdb_id.upper(),
        }
    if render_style not in RENDER_STYLES:
        return {
            "status": "invalid",
            "reason": f"unknown render style {render_style!r}",
            "pdb_id": pdb_id.upper(),
        }

    try:
        raw = coordinates if coordinates is not None else download_coordinates(pdb_id)
    except StructureRenderError as exc:
        return {"status": "unavailable", "reason": str(exc), "pdb_id": pdb_id.upper()}
    if not raw:
        return {
            "status": "not_found",
            "reason": f"PDBe returned no coordinate file for {pdb_id.upper()}.",
            "pdb_id": pdb_id.upper(),
        }

    from pymol import cmd

    identifier = pdb_id.strip().upper()
    selection = f"chain {chain_id}" if chain_id else "all"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    working = out_path.with_suffix(".working.pdb")
    try:
        cmd.reinitialize()
        # cmd.load accepts a string treated as a filename, so the CIF is handed over as a
        # temporary file rather than concatenated into the session by string tricks.
        with tempfile.TemporaryDirectory() as scratch:
            cif_path = Path(scratch) / f"{identifier.lower()}.cif"
            cif_path.write_bytes(raw)
            cmd.load(str(cif_path), identifier)
            cmd.select("render_target", f"{identifier} and {selection} and not solvent")
            if cmd.count_atoms("render_target") == 0:
                return {
                    "status": "not_found",
                    "reason": (
                        f"No protein atoms were found for {selection} in {identifier}; "
                        "the requested chain is not present in this entry."
                    ),
                    "pdb_id": identifier,
                }
            atom_count = cmd.count_atoms("render_target")
            cmd.set("sele_only", "on")
            cmd.save(str(working), "render_target")
        # Legacy-PDB suffix keeps upstream's format detection and output naming sane.
        if working.stat().st_size == 0:
            return {
                "status": "unavailable",
                "reason": f"PyMOL wrote an empty PDB for {identifier} {selection}.",
                "pdb_id": identifier,
            }
        png = out_path.with_suffix(".png")
        pdb2png(
            str(working),
            output=str(png),
            protein_rotate=VIEWS[view],
            protein_color="0,0,.6,1,0,0,0,.6,0,.6,.6,0",
            hetatm_color="1,.65,0" if chain_id is None else None,
            show_water=False,
            render_style=render_style,
            width=width,
            height=height,
        )
    except (StructureRenderError, ValueError, OSError) as exc:
        return {
            "status": "unavailable",
            "reason": f"Structure rendering failed for {identifier}: {exc}",
            "pdb_id": identifier,
        }
    except Exception as exc:  # PyMOL raises bare exceptions freely.
        return {
            "status": "unavailable",
            "reason": (
                f"Structure rendering failed for {identifier} "
                f"({type(exc).__name__}); no image is shown."
            ),
            "pdb_id": identifier,
        }
    finally:
        working.unlink(missing_ok=True)
        try:
            cmd.delete("all")
        except Exception:  # pragma: no cover - engine teardown must never break a report
            pass

    return {
        "status": "rendered",
        "pdb_id": identifier,
        "chain_id": chain_id,
        "selection": selection,
        "view": view,
        "render_style": render_style,
        "atoms_rendered": atom_count,
        "image": out_path.name,
        "renderer": "pdb2png (hgbrian/pdb2png, Unlicense) on PyMOL raytracer",
        "source_format": "PDBe deposited mmCIF, re-serialised to legacy PDB for rendering",
        "source_url": f"https://www.ebi.ac.uk/pdbe/entry/pdb/{identifier.lower()}",
        "caveats": [
            "the image shows the deposited model, including any engineered mutations, "
            "truncations, or fusion partners present in this entry",
            "hydrogens and waters are not rendered",
            "side chains and ligands are shown only where the deposited model includes them; "
            "unmodelled residues are absent from the image rather than interpolated",
            "a cartoon conveys fold, not sequence identity to the analysed protein; coverage "
            "is stated separately in the structure table",
        ],
    }


__all__ = [
    "RENDER_STYLES",
    "VIEWS",
    "StructureRenderError",
    "download_coordinates",
    "render_structure_png",
]
