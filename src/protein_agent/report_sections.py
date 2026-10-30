"""Report sections that combine retrieved evidence with derived predictions.

This module is the wiring between the single-sequence report and the analysis modules
that were built ahead of it. It owns three things the base report lacked:

* a post-translational-modification position map (annotated UniProt features on their
  own tracks, sequence-pattern predictions kept visibly separate),
* a predicted-motif inventory (N-glycan sequons, deamidation liabilities, cleavage),
* structure display chosen from the PDBe mappings the report already fetched.

Nothing here retrieves new facts on its own initiative. Callers pass in the UniProt
feature list and the PDBe structure list returned by tool calls earlier in this same
run; this module only ranks, plots, and renders them. The one exception is the
coordinate download for the entry that ranking selected, which is reported as such.
An empty or missing input produces a "not available" section, never an inferred one.
"""

# Template strings are intentionally readable rather than wrapped to the line limit.
# ruff: noqa: E501
from __future__ import annotations

import html
import logging
import re
from pathlib import Path
from typing import Any

from .cartoon import observed_coverage, parse_coordinates, plot_cartoon, plot_ligand_map
from .position_map import build_tracks, plot_position_map
from .sequons import predicted_motifs
from .structure_choice import choose_display_structure
from .structure_png import download_coordinates, render_structure_png

LOGGER = logging.getLogger(__name__)

# A PDBe entry identifier is four alphanumeric characters. Checked before rendering so a
# malformed value cannot reach an href or an alt attribute as though it were a real entry.
_PDB_ID = re.compile(r"[0-9A-Za-z]{4}\Z")

# Words that carry no information when deciding whether a reason merely restates the lead.
_DETAIL_STOPWORDS = frozenset(
    {"no", "not", "the", "a", "an", "was", "were", "been", "is", "are", "of", "in", "to",
     "and", "for", "this", "that", "under", "by", "on", "we", "it", "its"}
)


def _restates_lead(reason: str, lead: str) -> bool:
    """True when the reason adds nothing beyond vocabulary the lead already used.

    Any content word absent from the lead -- most importantly a number or an accession --
    makes the reason worth printing.
    """
    lead_words = _content_words(lead)
    reason_words = _content_words(reason)
    if not reason_words:
        return True
    return reason_words <= lead_words


def _content_words(text: str) -> frozenset[str]:
    return frozenset(
        w for w in re.findall(r"[a-z0-9]+", text.lower()) if w not in _DETAIL_STOPWORDS
    )

# A structure is only worth downloading and rasterising below this size; larger entries
# have stalled the report before and render as an unreadable crowd of chains.
MAX_RENDER_RESIDUES = 5_000
MAX_LIGAND_LABELS = 12

# Tracks that can only be populated from a retrieved UniProt entry. Without a verified
# accession these lists are absent, which is not the same as UniProt reporting none.
ANNOTATED_TRACKS = frozenset(
    {
        "Glycosylation (annotated)",
        "Other PTM",
        "Disulfide bond",
        "Active / binding site",
        "Transmembrane / topology",
        "Domain",
        "Variants",
    }
)


def _glycan_annotations(uniprot: dict[str, Any]) -> list[dict[str, Any]]:
    """Return only the glycosylation features already present in the fetched entry."""
    features = uniprot.get("features") or []
    return [item for item in features if item.get("type") == "Glycosylation"]


def build_motifs(sequence: str, uniprot: dict[str, Any]) -> dict[str, Any]:
    """Predict motifs, then state what the annotation could and could not confirm."""
    verified = uniprot.get("status") == "verified"
    payload = predicted_motifs(sequence, _glycan_annotations(uniprot) if verified else None)
    confirmed = (payload.get("counts") or {}).get("confirmed", 0)
    if not verified:
        payload["annotation_reconciliation"] = (
            "No verified UniProt accession was resolved for this sequence, so no predicted "
            "site could be reconciled against annotation. Every site listed is a sequence-pattern "
            "prediction whose occupancy is unknown."
        )
    elif confirmed:
        payload["annotation_reconciliation"] = (
            f"{confirmed} predicted site(s) coincide with a glycosylation annotation in the "
            "retrieved entry. Coincidence supports, but does not prove, occupancy."
        )
    else:
        payload["annotation_reconciliation"] = (
            "No predicted sequon coincides with a glycosylation annotation in the retrieved "
            "entry. Absence of an annotation is not absence of a glycan."
        )
    return payload


def build_position_map(
    sequence: str, uniprot: dict[str, Any], motifs: dict[str, Any], out_path: Path
) -> dict[str, Any]:
    """Plot annotated features and predicted sequons, keeping the two provenances apart."""
    verified = uniprot.get("status") == "verified"
    tracks = build_tracks(
        uniprot_features=(uniprot.get("features") or []) if verified else None,
        sequons=motifs.get("n_glycan_sequons"),
        deamidation=motifs.get("deamidation"),
    )
    if not verified:
        # An absent feature list is indistinguishable from "UniProt reports none", so the
        # annotated tracks are cleared rather than left empty-by-omission. The map then
        # shows predictions only, and says so.
        for name in ANNOTATED_TRACKS:
            tracks[name] = []
    result = plot_position_map(
        len(sequence), tracks, label="Analysed sequence", out_path=out_path
    )
    if result.get("status") == "rendered" and not verified:
        result["caveats"] = [
            *list(result.get("caveats", [])),
            "no UniProt accession was verified, so this map shows sequence-pattern "
            "predictions only",
        ]
    return result


def _rows(pairs: list[tuple[str, str]]) -> str:
    return "".join(
        f"<tr><th scope='row'>{html.escape(name)}</th><td>{html.escape(value)}</td></tr>"
        for name, value in pairs
    )


def render_motif_section(data: dict[str, Any]) -> str:
    """Summarise predicted PTM motifs, holding the prediction/annotation line explicitly."""
    motifs = data.get("predicted_motifs") or {}
    if motifs.get("status") != "computed":
        reason = motifs.get("reason") or "no analysed sequence was available"
        return f"<h2>Predicted PTM motifs</h2><p>Not computed: {html.escape(str(reason))}.</p>"

    sequons = motifs.get("n_glycan_sequons") or []
    deamidation = motifs.get("deamidation") or []
    cleavage = motifs.get("cleavage") or {}
    counts = motifs.get("counts") or {}
    method = motifs.get("method") or {}

    positions = ", ".join(str(item["position"]) for item in sequons[:MAX_LIGAND_LABELS]) or "none"
    if len(sequons) > MAX_LIGAND_LABELS:
        positions += f" (+{len(sequons) - MAX_LIGAND_LABELS} more)"
    sites = "".join(
        "<li>position {pos} · sequon {seq} · {state}{detail}</li>".format(
            pos=html.escape(str(item["position"])),
            seq=html.escape(str(item.get("sequon", ""))),
            state=html.escape(str(item.get("annotation_state", ""))),
            detail=(
                " ({})".format(html.escape(str(item["annotation_detail"])))
                if item.get("annotation_detail")
                else ""
            ),
        )
        for item in sequons[:MAX_LIGAND_LABELS]
    )
    list_items = f"<ul>{sites}</ul>" if sites else ""
    overflow = (
        f"<p>{len(sequons) - MAX_LIGAND_LABELS} further sequon(s) are listed in report.json.</p>"
        if len(sequons) > MAX_LIGAND_LABELS
        else ""
    )
    summary = _rows(
        [
            (
                "N-glycan sequons",
                f"{len(sequons)} ({counts.get('confirmed', 0)} annotated in the entry, "
                f"{counts.get('predicted_only', 0)} predicted only)",
            ),
            ("Deamidation liabilities", str(len(deamidation))),
            ("Predicted tryptic peptides", str(cleavage.get("peptide_count", "not returned"))),
        ]
    )
    return (
        "<h2>Predicted PTM motifs</h2>"
        f"<table>{summary}</table>"
        f"<p><strong>N-glycan sequon positions:</strong> {html.escape(positions)}.</p>"
        f"<p>{html.escape(str(motifs.get('verdict') or ''))}</p>"
        f"<p>{html.escape(str(motifs.get('annotation_reconciliation') or ''))}</p>"
        f"{list_items}{overflow}"
        "<p><em>Method:</em> "
        f"{html.escape(str(method.get('n_glycan_sequon', '')))}; "
        f"{html.escape(str(method.get('deamidation', '')))}. "
        "These are sequence-pattern predictions, not measurements of occupancy.</p>"
    )


def render_position_map_section(data: dict[str, Any]) -> str:
    """Embed the PTM position map, or state plainly why it is absent."""
    record = data.get("position_map") or {}
    if record.get("status") != "rendered":
        reason = record.get("reason") or "no position map was produced"
        return f"<h2>Modifications along the sequence</h2><p>Not shown: {html.escape(str(reason))}.</p>"

    drawn = ", ".join(
        f"{name} ({count})" for name, count in (record.get("tracks") or {}).items() if count
    )
    caveats = "".join(f"<li>{html.escape(str(item))}</li>" for item in record.get("caveats", []))
    return (
        "<h2>Modifications along the sequence</h2>"
        f"<p>Tracks drawn: {html.escape(drawn) or 'none'}.</p>"
        '<img src="positions.png" alt="Feature positions along the analysed sequence">'
        '<p><a href="positions.png" download>Download the position map</a></p>'
        f"<ul>{caveats}</ul>"
    )


def _cartoon(coordinates: bytes, length: int, output_dir: Path, title: str) -> dict[str, Any]:
    """Draw the CA trace and ligand map from coordinates already in hand."""
    try:
        parsed = parse_coordinates(coordinates)
    except Exception as exc:  # pragma: no cover - malformed deposits are upstream's bug
        LOGGER.warning("Coordinate parse failed: %s", exc)
        return {
            "cartoon": {
                "status": "unavailable",
                "reason": f"Coordinates could not be parsed ({type(exc).__name__}).",
            }
        }
    cartoon = plot_cartoon(parsed, output_dir / "structure.png", title=title)
    result: dict[str, Any] = {
        "cartoon": cartoon,
        "chains": len(parsed.get("chains", [])),
        "parsed_as": parsed.get("parsed_as"),
    }
    ligands = plot_ligand_map(
        parsed,
        observed_coverage(parsed),
        output_dir / "ligands.png",
        length=length,
        title=title,
    )
    result["ligand_map"] = ligands
    if ligands.get("status") == "rendered":
        result["bound_molecules"] = ligands.get("ligands", [])
    return result


def build_structure_display(data: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    """Choose and render one experimental structure for the report.

    Ranking reads only the PDBe mappings already fetched. The single new request is the
    coordinate download for the entry that ranking named, and its absence is reported
    rather than papered over.
    """
    pdbe = data.get("pdbe") or {}
    if pdbe.get("status") != "verified":
        # Propagate the upstream state: "PDBe has no structure for this accession" and "we
        # never asked because identity was unverified" are different findings, and flattening
        # both to one message would let a lookup gap read as an absence of experimental data.
        return {
            "status": pdbe.get("status") or "not_requested",
            "reason": pdbe.get("reason")
            or "A verified UniProt accession with PDBe mappings is required before display.",
        }
    length = int(data.get("length_residues") or 0)
    verdict = choose_display_structure(pdbe.get("structures") or [], sequence_length=length)
    selected = verdict.get("selected")
    if not selected:
        # Fragment-only rankings return no selected entry while still carrying the reason
        # the ranking declined it. Dropping that would report "nothing was found" when the
        # truth is "structures exist but none covers enough of this sequence to draw".
        return {
            "status": verdict.get("status") or "none",
            "reason": verdict.get("selection_reason")
            or verdict.get("reason")
            or "No usable PDBe mapping was returned.",
            "alternative_count": verdict.get("alternative_count", 0),
            "coverage_estimate": verdict.get("coverage_estimate"),
            "coverage_basis": verdict.get("coverage_basis"),
            "thresholds": verdict.get("thresholds"),
            "verdict": verdict,
        }

    pdb_id = str(selected["pdb_id"]).upper()
    display: dict[str, Any] = {
        "status": verdict.get("status"),
        "pdb_id": pdb_id,
        "selection_reason": verdict.get("selection_reason"),
        "coverage_estimate": verdict.get("coverage_estimate"),
        "coverage_basis": verdict.get("coverage_basis"),
        "alternative_count": verdict.get("alternative_count", 0),
        "chain_id": selected.get("chain_id"),
        "source_url": selected.get("source_url"),
        "experimental_method": selected.get("experimental_method"),
        "resolution_angstrom": selected.get("resolution_angstrom"),
    }

    if length > MAX_RENDER_RESIDUES:
        display["render"] = {
            "status": "too_large",
            "reason": (
                f"{length} residues exceeds the {MAX_RENDER_RESIDUES}-residue display limit."
            ),
        }
        return display

    # PyMOL is an optional extra. The C-alpha cartoon is the always-available rendering and
    # is drawn first; the PyMOL raster replaces it only when the engine is actually present,
    # so a node without PyMOL still gets a structure image rather than none.
    try:
        coordinates = download_coordinates(pdb_id)
    except Exception as exc:  # pragma: no cover - network failure is a reported state
        LOGGER.warning("Coordinate download failed for %s: %s", pdb_id, exc)
        display["render"] = {
            "status": "unavailable",
            "reason": f"Coordinates could not be retrieved ({type(exc).__name__}).",
        }
        return display

    display.update(_cartoon(coordinates, length, output_dir, pdb_id))
    raster = render_structure_png(
        pdb_id, out_path=output_dir / "structure.png", coordinates=coordinates
    )
    display["pymol_render"] = raster
    cartoon = display.get("cartoon") or {}
    if raster.get("status") == "rendered":
        display["render"] = raster
    else:
        # PyMOL declined (normally: it is not installed on this node), so what sits on disk
        # at structure.png is the C-alpha cartoon. Its caveats must survive the swap: the
        # cartoon states plainly that a residue absent from the model is not a residue
        # absent from the protein, and losing that would read as "the model is complete".
        display["render"] = {
            **cartoon,
            "pymol_status": raster.get("status"),
            "pymol_reason": raster.get("reason"),
            "caveats": [
                *list(cartoon.get("caveats", [])),
                *list(raster.get("caveats", [])),
            ],
        }
        display["cartoon_is_fallback"] = True
    return display


def _coverage_qualifier(
    status: object, coverage: object, basis: object
) -> str:
    """State what is known about coverage, wording it to match the evidence actually held.

    "PDBe gave a coverage number we could not corroborate" and "PDBe gave no number at all"
    are different findings and must not share a sentence.
    """
    if status == "fragment_only":
        shown = f"{coverage:.0%} " if isinstance(coverage, (int, float)) else ""
        return (
            "<p><strong>Partial model:</strong> this entry covers "
            f"{shown}of the analysed sequence. It depicts a construct, domain, or fragment, "
            "and is not a picture of the full-length protein.</p>"
        )
    if status == "partial":
        if isinstance(coverage, (int, float)):
            return (
                "<p><strong>Coverage approximate:</strong> this entry covers "
                f"{coverage:.0%} of the analysed sequence, estimated from {html.escape(str(basis or 'a PDBe coverage field'))}. "
                "It may be a construct or domain rather than the full-length protein.</p>"
            )
        return (
            "<p><strong>Coverage unconfirmed:</strong> PDBe returned no coverage field for "
            "this mapping, so completeness could not be estimated and is not asserted.</p>"
        )
    return ""


# The PDBe lookup states that justify a different sentence in the report. "not_found" is an
# assertion about this accession's coverage of the database; "unavailable" is a failure of the
# lookup itself and must never be phrased as though no structures exist.
# Whitelist of ranking verdicts that mean "an entry was actually selected and may be drawn".
# Deliberately a whitelist, not a blacklist: an unrecognised status must fall closed to "not
# shown", because this renderer's whole job is never to assert a structure without evidence.
_RENDERABLE_STATUSES = frozenset({"full", "partial", "fragment_only"})

_STRUCTURE_ABSENT_STATES = frozenset(
    {None, "", "not_requested", "not_found", "unavailable", "invalid", "none", "too_large"}
)


def _structure_absent_heading(status: object, reason: object) -> str:
    """Name the gap accurately, without implying absence of data from a failed lookup."""
    if status == "not_found":
        lead = "No experimental structure has been deposited in PDBe under this accession."
    elif status == "unavailable":
        lead = (
            "PDBe could not be queried for structures, so the absence of an image here says "
            "nothing about whether structures exist."
        )
    elif status == "not_requested":
        lead = (
            "No structure was requested: UniProt identity was not verified, and an unverified "
            "sequence must not be illustrated with a named entry."
        )
    elif status == "too_large":
        # Structures were found and one was chosen; only the image was withheld. Saying
        # "no mapping was returned" here would be factually wrong.
        lead = (
            "Structures were retrieved for this accession but no image is drawn, because the "
            "entry is too large to render legibly. The selected entry is still named below."
        )
    else:
        lead = "No usable PDBe mapping was returned for this sequence."
    detail = ""
    if reason:
        text = str(reason).strip().rstrip(".")
        # Suppress only a pure restatement of the lead. Any reason carrying a figure or an
        # identifier that the lead does not mention is kept, since that is the information
        # a reader needs in order to act on the gap.
        if text and not _restates_lead(text, lead):
            detail = f" {html.escape(text)}."
    return f"<h2>Representative structure</h2><p>Not shown: {lead}{detail}</p>"


def render_structure_section(data: dict[str, Any]) -> str:
    """Render the chosen structure, its provenance, and every caveat that qualifies it."""
    display = data.get("structure_display") or {}
    status = display.get("status")
    # A whitelist, not a blacklist: every recognised absent state, every missing or
    # non-string status, and every status this module does not recognise all land here
    # rather than in the success template, which would draw a structure box as though an
    # entry had actually been selected. The specific status is passed through so a refusal
    # keeps its own wording instead of collapsing to the generic one.
    if not isinstance(status, str) or status.strip() not in _RENDERABLE_STATUSES:
        normalised = status.strip() if isinstance(status, str) else None
        return _structure_absent_heading(normalised, display.get("reason"))

    # Never coerce: str(1234) == "1234" would satisfy the four-character check and let a
    # non-string value be attributed to a report as if it were a real deposited entry.
    declared_id = display.get("pdb_id")
    raw_id = declared_id.strip() if isinstance(declared_id, str) else ""
    if not raw_id:
        return _structure_absent_heading(
            status, display.get("reason") or "no PDBe entry identifier was returned"
        )
    if not _PDB_ID.fullmatch(raw_id):
        return _structure_absent_heading(
            status, f"{raw_id!r} is not a well-formed four-character PDBe entry identifier"
        )
    pdb_id = html.escape(raw_id.upper(), quote=True)
    url = html.escape(str(display.get("source_url") or ""), quote=True)
    facts = [str(display.get("experimental_method") or "method not returned")]
    resolution = display.get("resolution_angstrom")
    if resolution is not None:
        facts.append(f"{resolution} Å")
    coverage = display.get("coverage_estimate")
    if isinstance(coverage, (int, float)):
        facts.append(f"{coverage:.0%} of the analysed sequence")
    if display.get("chain_id"):
        facts.append(f"chain {display['chain_id']}")

    qualifier = _coverage_qualifier(status, coverage, display.get("coverage_basis"))

    render = display.get("render") or {}
    render_status = render.get("status")
    images = ""
    if render_status == "rendered":
        caveats = "".join(
            f"<li>{html.escape(str(item))}</li>" for item in (render.get("caveats") or [])
        )
        images += (
            f'<img src="structure.png" alt="Structure of {pdb_id}">'
            '<p><a href="structure.png" download>Download the structure image</a> · '
            + (
                "C-alpha trace rendered from deposited coordinates; the PyMOL engine is not "
                "installed on this node, so side chains, surface, and depth shading are absent."
                if display.get("cartoon_is_fallback")
                else "PyMOL cartoon rendering of the deposited coordinates."
            )
            + "</p>"
            + (f"<ul>{caveats}</ul>" if caveats else "")
        )
    else:
        images += (
            "<p>Structure image not shown: "
            f"{html.escape(str(render.get('reason') or 'rendering was unavailable'))}.</p>"
        )

    ligands = display.get("bound_molecules") or []
    if ligands:
        names = ", ".join(
            f"{item['name']} ({item.get('instances', 1)} modelled)"
            for item in ligands[:MAX_LIGAND_LABELS]
        )
        images += (
            '<img src="ligands.png" alt="Bound molecules modelled in this entry">'
            '<p><a href="ligands.png" download>Download the ligand map</a> · Modelled non-water '
            f"molecules: {html.escape(names)}.</p>"
        )
    elif (display.get("ligand_map") or {}).get("status") == "empty":
        images += "<p>No non-water ligand is modelled in these coordinates.</p>"

    alternatives = display.get("alternative_count") or 0
    alternative_text = (
        f"<p>{alternatives} further mapped entr{'y' if alternatives == 1 else 'ies'} exist for "
        "this accession; the full ranked list is in the experimental structures table.</p>"
        if alternatives
        else ""
    )
    return (
        "<h2>Representative structure</h2>"
        f'<p>Selected <a href="{url}"><strong>{pdb_id}</strong></a> from the PDBe mappings for '
        f"this accession · {html.escape(' · '.join(facts))}.</p>"
        "<p><em>Why this entry:</em> "
        f"{html.escape(str(display.get('selection_reason') or 'not returned'))}.</p>"
        f"{qualifier}{images}{alternative_text}"
        "<p>One structure represents this accession. It is not necessarily the best, the newest, "
        "or the closest to the analysed isoform, and deposited coordinates describe that deposit "
        "alone.</p>"
    )


__all__ = [
    "build_motifs",
    "build_position_map",
    "build_structure_display",
    "render_motif_section",
    "render_position_map_section",
    "render_structure_section",
]
