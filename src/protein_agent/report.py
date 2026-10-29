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

from .plots import plot_hydropathy
from .properties import STANDARD_AA, SequenceError, compute_properties

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
        _row(aa, f"{fraction:.2%}") for aa, fraction in sorted(data["amino_acid_percent"].items())
    )
    charge = "".join(
        _row(f"pH {ph}", value)
        for ph, value in data["charge_vs_ph"]["ph_to_net_charge"].items()
    )
    classes = "".join(
        _row(name, f"{values['mole_percent']}%")
        for name, values in data["physico_chemical_classes"].items()
    )
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
<details><summary>Amino-acid composition</summary><table>{composition}</table></details>
<details><summary>Charge at different pH values</summary><table>{charge}</table></details>
<details><summary>Residue groups</summary><table>{classes}</table>
<p>Groups overlap, so percentages do not add to 100%.</p></details>
<h2>Download all calculated results</h2>
<p><a href="report.json" download>Download full results (JSON)</a> — includes additional
composition statistics, extinction estimates, and calculation details.</p>
<div class="note"><strong>What this report cannot tell you:</strong> These are calculations
from your sequence, not measurements. The sequence has not been identified in a database.
No structure, disorder, or expression outcome was determined.</div></body></html>
"""


def generate_report(input_path: Path, output_dir: Path) -> Path:
    if input_path.stat().st_size > MAX_INPUT_BYTES:
        raise SequenceError("input exceeds the 100 KB limit")
    header, sequence = parse_sequence(input_path.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    data = compute_properties(sequence, header=header).as_dict()
    data["input_source"] = "user-supplied sequence (not database-identified)"
    data["generated_at_utc"] = datetime.now(UTC).isoformat(timespec="seconds")
    plot_hydropathy(sequence, label=header, out_path=output_dir / "hydropathy.png")
    (output_dir / "report.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    report = output_dir / "report.html"
    report.write_text(render_html(data, generated_at=data["generated_at_utc"]), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description="Create an offline report from one protein sequence")
    parser.add_argument("--input", type=Path, required=True, help="single FASTA record or bare sequence")
    parser.add_argument("--output", type=Path, required=True, help="report directory")
    args = parser.parse_args()
    try:
        report = generate_report(args.input, args.output)
    except (SequenceError, OSError, UnicodeError, ValueError) as exc:
        print(f"Report not generated: {exc}", file=sys.stderr)
        return 2
    print(f"Report generated: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
