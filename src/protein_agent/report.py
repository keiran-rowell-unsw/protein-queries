"""Offline, deterministic single-sequence report for the local web app.

No UniProt identity, structure, expression, or experimental claims are inferred from a
pasted sequence. This entry point never invokes the inference gateway or requires a key.
"""

# HTML template is intentionally readable rather than wrapped to Python's line limit.
# ruff: noqa: E501, RUF001
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
        _row("GRAVY (Kyte–Doolittle)", data["gravy_kyte_doolittle"]),
        _row("Aromaticity", data["aromaticity_iron_dos_santos"]),
        _row("Instability index", data["instability_index_guruprasad"]),
        _row("Extinction coefficient (reduced)", f"{ext['reduced']} M⁻¹ cm⁻¹"),
        _row("Extinction coefficient (oxidised)", f"{ext['oxidised']} M⁻¹ cm⁻¹"),
        _row("Average residue weight", f"{data['average_residue_weight_da']} Da"),
    ]
    composition = "".join(
        _row(aa, f"{fraction:.2%}") for aa, fraction in sorted(data["amino_acid_percent"].items())
    )
    methods = "".join(_row(key.replace("_", " "), value) for key, value in data["method"].items())
    caveats = "".join(f"<li>{html.escape(note)}</li>" for note in data["caveats"])
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Sequence analysis report</title><style>
body{{font:1rem/1.5 system-ui,sans-serif;max-width:850px;margin:2rem auto;padding:0 1rem;color:#17243a}}
h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;margin-bottom:2rem}}
th,td{{padding:.45rem;border-bottom:1px solid #ccd;text-align:left;vertical-align:top}}
th{{width:48%}}img{{max-width:100%;height:auto}}.note{{background:#eef2f7;padding:1rem}}
</style></head><body><h1>Sequence analysis report</h1>
<p><strong>Input:</strong> {html.escape(data['header'])}</p>
<p>Generated: {html.escape(generated_at)} (UTC). Input: user supplied; identity unresolved.</p>
<h2>Computed properties</h2><table>{''.join(rows)}</table>
<h2>Amino-acid composition</h2><table>{composition}</table>
<h2>Hydropathy profile</h2><p>Kyte–Doolittle, 19-residue centred moving average; values are composition-derived.</p>
<img src="hydropathy.png" alt="Hydropathy profile across the submitted sequence">
<h2>Methods</h2><table>{methods}</table>
<h2>Limitations</h2><div class="note"><ul>{caveats}
<li>No UniProt identity, experimental structure, disorder prediction, or expression result was obtained.</li>
<li>Local calculations are not experimental measurements. See <a href="report.json">report.json</a> for machine-readable results, including charge versus pH and pepstats-equivalent metrics.</li>
</ul></div></body></html>
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
