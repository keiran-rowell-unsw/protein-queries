"""Offline, deterministic single-sequence report for the local web app.

No UniProt identity, structure, expression, or experimental claims are inferred from a
pasted sequence. This entry point never invokes the inference gateway or requires a key.
"""

# HTML template is intentionally readable rather than wrapped to Python's line limit.
# ruff: noqa: E501
from __future__ import annotations

import argparse
import html
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from .disorder import predict_disorder
from .plots import plot_hydropathy
from .properties import STANDARD_AA, SequenceError, chromophore_candidates, compute_properties
from .uniprot_report import fetch_annotations, find_sequence_annotations

MAX_INPUT_BYTES = 100_000
MAX_RESIDUES = 10_000


def parse_sequence(text: str) -> tuple[str, str]:
    """Accept one FASTA record or bare sequence; reject ambiguous/lossy inputs."""
    if len(text.encode("utf-8")) > MAX_INPUT_BYTES:
        raise SequenceError("input exceeds the 100 KB limit")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise SequenceError("no sequence provided")
    if any(line.startswith(">") for line in lines[1:]):
        raise SequenceError("provide exactly one FASTA record")
    header = lines.pop(0)[1:].strip() if lines[0].startswith(">") else "User sequence"
    if header != "User sequence" and not header:
        raise SequenceError("FASTA header is empty")
    if not lines:
        raise SequenceError("no residues found")
    sequence = "".join("".join(line.split()) for line in lines).upper()
    if len(sequence) > MAX_RESIDUES:
        raise SequenceError("sequence exceeds the 10,000-residue limit")
    invalid = sorted(set(sequence) - STANDARD_AA)
    if invalid:
        raise SequenceError(
            "unsupported characters: " + ", ".join(repr(char) for char in invalid)
            + "; use the 20 standard amino-acid letters (no gaps or ambiguity codes)"
        )
    return header[:200], sequence


def _row(name: str, value: object) -> str:
    return f"<tr><th scope='row'>{html.escape(name)}</th><td>{html.escape(str(value))}</td></tr>"


def _uniprot_section(record: dict) -> str:
    """Display only fetched fields for an exact canonical-sequence match."""
    if record["status"] == "not_requested":
        return (
            "<h2>UniProt annotations</h2><p>No accession was supplied. "
            "To add annotations, enter a UniProt accession with your sequence.</p>"
        )
    if record["status"] != "verified":
        candidates = record.get("candidates") or []
        suggestion = (
            " Candidate accessions: " + ", ".join(html.escape(acc) for acc in candidates) + "."
            if candidates else ""
        )
        return (
            "<h2>UniProt annotations</h2><p>Not included: "
            + html.escape(record["reason"]) + suggestion + "</p>"
        )
    rows = []
    for feature in record["features"]:
        start, end = feature["start"], feature["end"]
        position = str(start) if start == end else f"{start}-{end}"
        evidence = ", ".join(
            " ".join(str(item[key]) for key in ("code", "source", "id") if item.get(key))
            for item in feature["evidence"]
        )
        text = f"Residue(s) {position}: {feature['description'] or 'No description annotated'}"
        if evidence:
            text += f" (evidence: {evidence})"
        rows.append(_row(feature["type"], text))
    annotations = (
        f"<table>{''.join(rows)}</table>" if rows else
        "<p>No matching processing or modification features are annotated in this entry.</p>"
    )
    absorption = record["absorption"]
    if absorption:
        maximum = absorption.get("max_nm")
        notes = "; ".join(absorption.get("note") or [])
        absorption_text = html.escape(
            f"Annotated absorption maximum: {maximum} nm" if maximum is not None
            else "Absorption annotated (no maximum reported)"
        )
        if notes:
            absorption_text += ". " + html.escape(notes)
        evidence = ", ".join(
            " ".join(str(item[key]) for key in ("code", "source", "id") if item.get(key))
            for item in absorption.get("evidence") or []
        )
        if evidence:
            absorption_text += " (evidence: " + html.escape(evidence) + ")"
        annotations += f"<p>{absorption_text}</p>"
    else:
        annotations += "<p>Absorption is not annotated in this entry.</p>"
    url = html.escape(record["source_url"], quote=True)
    identity = html.escape(
        f"{record['accession']} · {record['organism'] or 'organism not provided'} "
        f"· {record['entry_type'] or 'review status not provided'}"
    )
    date = html.escape(record["retrieved_at"] or "retrieval date not available")
    return (
        "<h2>UniProt annotations</h2>"
        f"<p>Canonical sequence matched exactly: <a href=\"{url}\">{identity}</a> "
        f"(retrieved {date}). Positions refer to that canonical sequence.</p>"
        + annotations
    )


def render_html(data: dict, *, generated_at: str) -> str:
    """Render only escaped input/data; never execute arbitrary HTML from a FASTA header."""
    ext = data["molar_extinction_coefficient_M-1cm-1"]
    rows = [
        _row("Length", f"{data['length_residues']} residues"),
        _row("Molecular weight", f"{data['molecular_weight_da']} Da"),
        _row("Isoelectric point (pI)", data["isoelectric_point"]),
        _row("Net charge at pH 7", data["net_charge_pH_7"]),
        _row("Average hydropathy (GRAVY)", data["gravy_kyte_doolittle"]),
        _row("Aromatic amino-acid fraction", f"{data['aromaticity_iron_dos_santos']:.1%}"),
        _row("Instability index (sequence estimate)", data["instability_index_guruprasad"]),
        _row("A280 extinction coefficient (reduced)", f"{ext['reduced']} M⁻¹ cm⁻¹"),
        _row("A280 extinction coefficient (oxidised)", f"{ext['oxidised']} M⁻¹ cm⁻¹"),
    ]
    composition = "".join(
        _row(aa, f"{item['number']} residues · {item['mole_percent']}% · "
             f"{item['dayhoff_stat']}x reference frequency")
        for aa, item in sorted(data["amino_acid_dayhoff_stats"].items())
        if aa in STANDARD_AA
    )
    charge = "".join(
        _row(f"pH {ph}", value)
        for ph, value in data["charge_vs_ph"]["ph_to_net_charge"].items()
    )
    classes = "".join(
        _row(name, f"{values['mole_percent']}%")
        for name, values in data["physico_chemical_classes"].items()
    )
    disorder = data["disorder_prediction"]
    if disorder["status"] == "predicted":
        regions = disorder["disordered_regions"]
        range_text = ", ".join(f"{item['start']}-{item['end']}" for item in regions)
        region_text = f"Predicted regions: {range_text}." if regions else "No regions called by this predictor."
        disorder_section = (
            "<h2>Predicted disorder</h2>"
            f"<p>{html.escape(region_text)} This is a sequence-based prediction, not a "
            "measured structure.</p>"
            '<img src="disorder.png" alt="Predicted disorder scores along the sequence">'
            '<p><a href="disorder.png" download>Download the disorder plot</a></p>'
        )
    else:
        disorder_section = (
            "<h2>Predicted disorder</h2>"
            f"<p>Not available: {html.escape(disorder['reason'])}</p>"
        )
    motifs = data["chromophore_candidates"]
    motif_section = ""
    if motifs["status"] == "candidates_found":
        positions = sorted({pos for group in motifs["candidates"].values() for pos in group})
        motif_section = (
            "<details><summary>Possible chromophore motif matches</summary>"
            f"<p>Pattern matches start at residue(s) {', '.join(map(str, positions))}. "
            "A match alone does not establish that the protein absorbs light or forms "
            "a chromophore. Verify against an identified, annotated protein entry.</p></details>"
        )
    ext_mg = data["a280_extinction_per_mg_ml_M-1cm-1"]
    estimates = "".join([
        _row("Average residue weight", f"{data['average_residue_weight_da']} Da"),
        _row("A280 per 1 mg/mL (reduced / oxidised)",
             f"{ext_mg['reduced']} / {ext_mg['oxidised']}"),
    ])
    structure = data["secondary_structure_composition_estimate"]
    structure_rows = "".join(
        _row(label.title(), f"{value:.1%}") for label, value in structure.items()
    )
    expression = data["inclusion_body_propensity"]
    expression_flag = "higher" if expression["class"] == "likely_insoluble" else "lower"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sequence analysis report</title><style>
body{{font:1rem/1.5 system-ui,sans-serif;max-width:850px;margin:2rem auto;padding:0 1rem;color:#17243a}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;margin-bottom:2rem}}
th,td{{padding:.45rem;border-bottom:1px solid #ccd;text-align:left;vertical-align:top}}
th{{width:48%}}img{{max-width:100%;height:auto}}.note{{background:#eef2f7;padding:1rem}}
details{{margin:1rem 0;padding:.5rem;border:1px solid #ccd}}summary{{cursor:pointer;font-weight:600}}
a{{color:#0057b8}}
</style></head><body><h1>Sequence analysis report</h1>
<p><strong>Your sequence:</strong> {html.escape(data['header'])} · {data['length_residues']} residues</p>
<p>Report created {html.escape(generated_at)} UTC. <a href="../../">Analyse another sequence</a></p>
<h2>At a glance</h2><table>{''.join(rows)}</table>
<h2>Hydropathy along the sequence</h2>
<p>Positive values are more hydrophobic; negative values are more hydrophilic.</p>
<img src="hydropathy.png" alt="Hydropathy profile across the submitted sequence">
<p><a href="hydropathy.png" download>Download the hydropathy image</a></p>
{disorder_section}
<details><summary>Amino-acid composition</summary>
<p>Counts, share of the sequence, and enrichment relative to a reference composition.</p>
<table>{composition}</table></details>
<details><summary>Charge at different pH values</summary><table>{charge}</table></details>
<details><summary>Residue groups</summary><table>{classes}</table>
<p>Groups overlap, so percentages do not add to 100%.</p></details>
<details><summary>Other sequence estimates</summary><table>{estimates}</table>
<p>Secondary-structure composition estimate (not a structure prediction):</p>
<table>{structure_rows}</table>
<p>Exploratory E. coli inclusion-body composition flag: {expression_flag}. This simple
sequence-only flag does not establish actual expression or solubility.</p></details>
{motif_section}
{_uniprot_section(data['uniprot'])}
<h2>Download all calculated results</h2>
<p><a href="report.json" download>Download full results (JSON)</a> — includes additional
composition statistics, disorder scores (if available), and calculation details.</p>
<div class="note"><strong>What this report cannot tell you:</strong> These are sequence-derived
estimates, not measurements. A UniProt entry is shown only when its canonical sequence
matches exactly; a match does not prove an experimental structure or expression outcome.
A missing annotation or prediction is not evidence that a property is absent.</div></body></html>
"""


def generate_report(input_path: Path, output_dir: Path, *, accession: str = "") -> Path:
    if input_path.stat().st_size > MAX_INPUT_BYTES:
        raise SequenceError("input exceeds the 100 KB limit")
    header, sequence = parse_sequence(input_path.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    data = compute_properties(sequence, header=header).as_dict()
    data["input_source"] = "user-supplied sequence; UniProt identity is assessed separately"
    data["generated_at_utc"] = datetime.now(UTC).isoformat(timespec="seconds")
    plot_hydropathy(sequence, label=header, out_path=output_dir / "hydropathy.png")
    data["disorder_prediction"] = predict_disorder(sequence, output_dir / "disorder.png")
    data["chromophore_candidates"] = chromophore_candidates(sequence)
    data["uniprot"] = (
        fetch_annotations(sequence, accession) if accession.strip()
        else find_sequence_annotations(sequence)
    )
    (output_dir / "report.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    report = output_dir / "report.html"
    report.write_text(render_html(data, generated_at=data["generated_at_utc"]), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Create an offline report from one protein sequence")
    parser.add_argument("--input", type=Path, required=True, help="single FASTA record or bare sequence")
    parser.add_argument("--output", type=Path, required=True, help="report directory")
    parser.add_argument("--accession", default="", help="optional UniProt accession; otherwise search exact sequence in UniParc")
    args = parser.parse_args()
    try:
        report = generate_report(args.input, args.output, accession=args.accession)
    except (SequenceError, OSError, UnicodeError, ValueError) as exc:
        print(f"Report not generated: {exc}", file=sys.stderr)
        return 2
    print(f"Report generated: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
