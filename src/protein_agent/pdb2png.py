"""Vendored pdb2png (https://github.com/hgbrian/pdb2png), released under the Unlicense.

Public domain: no attribution or licence conditions apply. The upstream file is kept
verbatim below the header, with three departures, each marked ``# [protein-agent]`` in
place because the report must be reproducible and this node has no GPU:

1. ``from pymol import cmd`` is deferred into the functions that need it. Importing this
   module must not require PyMOL; ``structure_png`` reports the missing dependency instead
   of crashing the report, matching the optional-dependency contract in ``disorder.py``.
2. ``pdb2png()`` gained an ``output`` argument. Upstream derives the PNG name with
   ``pdb_file.replace(".pdb", ".png")``, which silently produces a file called
   ``1m5c.cif`` when handed the CIF/other extension PDBe actually serves, and gives the
   caller no way to choose the path.
3. ``render_structure()`` accepts a pre-rotated view list so the caller can request the
   orientations it wants rather than only upstream's single default.

Upstream copyright: none (Unlicense). Upstream last updated 2023-08-20.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Upstream body (verbatim except where marked).
# ---------------------------------------------------------------------------

# [protein-agent] was: from pymol import cmd  -- deferred; see module docstring.

RENDER_OPTIONS = {
    "default": {
        "ray_opaque_background": "off",
        "antialias": "2",
        "orthoscopic": "on",
        "depth_cue": "0",
        "ray_trace_mode": "1",
        "ray_trace_color": "black",
    },
    "default_bw": {
        "bg_color": "white",
        "ray_opaque_background": "on",
        "antialias": "2",
        "orthoscopic": "on",
        "depth_cue": "0",
        "ray_trace_mode": "2",
    },
    "default_cartoon": {
        "ray_opaque_background": "off",
        "antialias": "2",
        "orthoscopic": "on",
        "depth_cue": "0",
        "ray_trace_mode": "3",
    },
    "dark": {
        "bg_color": "black",
        "ray_opaque_background": "on",
        "antialias": "2",
        "orthoscopic": "on",
        "light_count": "2",
        "specular": "1",
        "depth_cue": "1",
        "fog_start": "0.35",
        "ray_trace_mode": "1",
        "ray_trace_color": "black",
    },
    "muted": {
        "bg_color": "white",
        "valence": "0",
        "bg_rgb": "white",
        "reflect": "0",
        "spec_direct": "0",
        "light_count": "1",
        "spec_count": "0",
        "shininess": "0",
        "power": "1",
        "specular": "0",
        "ambient_occlusion_mode": "1",
        "ambient_occlusion_scale": "15",
        "ambient_occlusion_mix": "1",
        "ambient_occlusion_samples": "8",
        "ray_trace_mode": "0",
    },
}

DEFAULT_PROTEIN_COLORS = "0,0,.6,1,0,0,0,.6,0"
DEFAULT_HETATM_COLORS = "1,.65,0"


def apply_render_style(style_name: str, custom_options: dict | None = None):
    """Apply a named render style or custom options to pymol."""
    from pymol import cmd  # [protein-agent] deferred import

    if custom_options:
        for key, value in custom_options.items():
            cmd.set(str(key), str(value))

    # If style_name is empty or None, do nothing and return
    if not style_name:
        return

    # Apply the render style options if the style exists
    style = RENDER_OPTIONS.get(style_name)
    if style:
        for setting, value in style.items():
            cmd.set(setting, str(value))
    else:
        raise ValueError(f"Unknown render style: {style_name}")


def pdb2png(
    pdb_file: str,
    *,
    output: str | None = None,  # [protein-agent] explicit destination
    protein_rotate: tuple[float, float, float] | list[float] | None = None,
    protein_color: tuple[float, float, float] | str | None = None,
    protein_zoom: float | None = None,
    hetatm_color: tuple[float, float, float] | str | None = None,
    ligand_id: str | None = None,
    ligand_chain: str | None = None,
    ligand_zoom: float = 1.0,
    ligand_color: tuple[float, float, float] | str | None = "red",
    show_water: bool = False,
    render_style: str = "default",
    width: int = 1600,
    height: int = 1600) -> str:
    """
    Input is a pdb file.
    Output is a png file.
    """
    from pymol import cmd  # [protein-agent] deferred import

    # reinitialize is important if you call the function multiple times!
    cmd.reinitialize()

    cmd.load(pdb_file)
    cmd.orient() # maybe not needed

    if protein_rotate is not None:
        cmd.rotate("x", protein_rotate[0])
        cmd.rotate("y", protein_rotate[1])
        cmd.rotate("z", protein_rotate[2])

    if protein_color is not None:
        if isinstance(protein_color, tuple):
            n = 0
            for chain in cmd.get_chains():
                cmd.set_color("protein_color", protein_color[n:n+3])
                cmd.color("protein_color", f"chain {chain} and not hetatm")
                n = (n + 3) % len(protein_color)
        else:
            cmd.color(protein_color, "not hetatm")

    # Color proteins and hetatms
    for hp_id, hp_color, hp_sel in [("protein", protein_color, "not hetatm"),
                                    ("hetatm", hetatm_color, "hetatm")]:
        if hp_color is not None:
            if isinstance(hp_color, tuple):
                n = 0
                for chain in cmd.get_chains():
                    cmd.select(f"sel_{hp_id}_{chain}", f"chain {chain} and {hp_sel}")
                    if cmd.count_atoms(f"sel_{hp_id}_{chain}") > 0:
                        cmd.set_color(f"{hp_id}_color_{chain}", hp_color[n:n+3])
                        cmd.color(f"{hp_id}_color_{chain}", f"sel_{hp_id}_{chain}")
                        n = (n + 3) % len(hp_color)
            else:
                cmd.color(hp_color, hp_sel)

    if protein_zoom is not None:
        cmd.zoom("all", protein_zoom)

    if ligand_id is not None:
        and_chain = f" and chain {ligand_chain}" if ligand_chain else ""
        cmd.select("ligand", f"resn {ligand_id}{and_chain}")
        cmd.zoom("ligand", ligand_zoom)

        if ligand_color is not None:
            if isinstance(ligand_color, tuple):
                cmd.set_color("ligand_color", ligand_color)
                cmd.color("ligand_color", "ligand")
            else:
                cmd.color(ligand_color, "ligand")

    if not show_water:
        cmd.select("HOH", "resn HOH")
        cmd.hide("everything", "HOH")

    apply_render_style(render_style)

    cmd.ray(width, height)
    destination = output or pdb_file.replace(".pdb", ".png")  # [protein-agent]
    cmd.save(destination, pdb_file)
    return destination


__all__ = ["RENDER_OPTIONS", "apply_render_style", "pdb2png"]
