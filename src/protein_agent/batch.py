"""Batch designability screen: many candidate sequences, one comparable table.

The single-sequence report answers "what is true of this protein". This answers the
selection question a designer actually asks first: "of these candidates, which is worth
building" — glycosylation liabilities, deamidation sites, cysteine pairing, aggregation
and expression signals, compared side by side instead of one webform at a time.

Designability, not identity
--------------------------
Every ranked column is derived from the sequence alone and therefore exists for a novel
binder with no database record. Identity is *reported*, never *assumed*: an accession is
only attached after exact canonical-sequence verification, and unverified rows say so. A
ranked table of unverifiable candidates is still useful; a confident table of invented
annotations is not.

Provenance is preserved through the rollup, not flattened into it. Motif counts are
labelled predicted, annotated glycans are labelled curated, and the two are never summed
into a single "glycosylation" number — that is the distinction the whole report rests on.

No inference gateway, and no model call. The ranking is a documented sum of computed
values, so it can be audited, disagreed with, and reproduced.
"""

from __future__ import annotations

import argparse
import html
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .properties import SequenceError, clean_fasta_sequence, compute_properties
from .sequons import predicted_motifs

MAX_RECORDS = 200
MAX_RESIDUES_PER_RECORD = 10_000
MAX_INPUT_BYTES = 5_000_000

# Ranking weights are a stated opinion, not a measurement. They order candidates; they do
# not classify anything as buildable. Deliberately coarse: a 0.1-point weight cannot be
# defended, and a designer should override the whole column with their own judgement.
RANK_WEIGHTS: dict[str, float] = {
    "instability_index_guruprasad": 0.30,  # higher = worse, normalised 0-100
    "gravy_abs": 0.20,  # extreme hydrophobicity in either direction
    "n_glycan_sequons": 0.20,  # predicted; occupancy unknown
    "deamidation_sites": 0.15,  # N-G counts double: the fast class
    "free_cysteines": 0.10,  # unpaired Cys, from sequence when no disulfides annotated
    "long_hydrophobic_stretch": 0.05,
}
FAST_DEAMIDATION_MULTIPLIER = 2  # N-G deamidation is the fast class; weighted heavier
LONG_STRETCH_REPORT_LIMIT = 40


@dataclass(slots=True)
class Candidate:
    """One row of the screen, with its evidence state attached."""

    label: str
    sequence: str
    properties: dict[str, Any]
    motifs: dict[str, Any]
    uniprot: dict[str, Any]
    warnings: list[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.warnings is None:
            self.warnings = []

    @property
    def accession(self) -> str:
        return str(self.uniprot.get("accession") or "")

    @property
    def verified(self) -> bool:
        return self.uniprot.get("status") == "verified"


def read_multi_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse every record in a FASTA file; blank and malformed records are rejected."""
    if path.stat().st_size > MAX_INPUT_BYTES:
        raise SequenceError(f"input exceeds the {MAX_INPUT_BYTES // 1_000_000} MB limit")
    text = path.read_text(encoding="utf-8")
    records: list[tuple[str, str]] = []
    header: str | None = None
    body: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header is not None:
                records.append((header, "".join(body)))
            header, body = stripped[1:].strip(), []
        elif header is None:
            header, body = "(no header)", [stripped]
        else:
            body.append(stripped)
    if header is not None:
        records.append((header, "".join(body)))
    if not records:
        raise SequenceError("no FASTA records found")
    return [(label or "(unlabelled)", seq) for label, seq in records]


def _longest_hydrophobic_stretch(sequence: str, threshold: float = 0.0) -> int:
    """Longest run whose centred-window hydropathy stays above `threshold`.

    Local import: this keeps the batch screen importable for testing without pulling
    matplotlib in on nodes where the plotting backend is absent.
    """
    from .properties import hydropathy_profile

    profile = hydropathy_profile(sequence)
    longest = current = 0
    for _position, value in profile:
        if value > threshold:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def _counts(candidate: Candidate) -> dict[str, int]:
    motifs = candidate.motifs
    deamidation = motifs.get("deamidation") or []
    fast = sum(1 for item in deamidation if item.get("rate_class") == "fast")
    sequons = motifs.get("raw_sequon_positions") or []
    annotated = [
        item for item in (candidate.uniprot.get("features") or [])
        if item.get("type") == "Glycosylation"
    ] if candidate.verified else []
    disulfide_pairs = {
        frozenset((item.get("start"), item.get("end")))
        for item in (candidate.uniprot.get("features") or [])
        if item.get("type") == "Disulfide bond" and item.get("start") and item.get("end")
    } if candidate.verified else set()
    bonded = {position for pair in disulfide_pairs for position in pair if isinstance(position, int)}
    cysteine_count = sum(1 for residue in candidate.sequence if residue == "C")
    # Fast (N-G) deamidation is weighted heavier, but the reported figure must stay an
    # integer count: a screening table showing "2.0" sites implies a precision it does not
    # have, and breaks any consumer that treats the column as an int.
    deamidation_score = int(fast * FAST_DEAMIDATION_MULTIPLIER + (len(deamidation) - fast))
    return {
        "n_glycan_sequons": len(sequons),
        "annotated_glycans": len(annotated),
        "deamidation_sites": deamidation_score,
        "free_cysteines": max(0, cysteine_count - len(bonded)),
        "cysteines": cysteine_count,
        "long_hydrophobic_stretch": _longest_hydrophobic_stretch(candidate.sequence),
    }


def analyse_record(
    label: str,
    sequence: str,
    *,
    lookup_identity: bool = True,
    supplied_uniprot: dict[str, Any] | None = None,
) -> Candidate:
    """Compute the designability row for one sequence, verifying identity if permitted.

    Sequence-only by default. ``lookup_identity`` attaches the UniProt annotation block,
    which requires a network call and only ever applies to an exact canonical match.
    """
    cleaned, warnings = clean_fasta_sequence(sequence)
    if not cleaned:
        raise SequenceError(f"{label}: no standard amino acids to analyse")
    if len(cleaned) > MAX_RESIDUES_PER_RECORD:
        raise SequenceError(f"{label}: exceeds the {MAX_RESIDUES_PER_RECORD:,} residue limit")
    properties = compute_properties(cleaned, header=label).as_dict()
    uniprot: dict[str, Any] = {"status": "not_requested"}
    annotated_glycans: list[dict] = []
    if supplied_uniprot is not None:
        # The accession was named by the user and fetched with accession precedence, so it
        # is authoritative. Re-deriving it by checksum would second-guess a deliberate
        # choice into "ambiguous" and discard the annotation the user asked for.
        uniprot = supplied_uniprot
    elif lookup_identity:
        from .uniprot_report import find_sequence_annotations

        uniprot = find_sequence_annotations(cleaned)
        if uniprot.get("status") == "verified":
            annotated_glycans = [
                item for item in (uniprot.get("features") or [])
                if item.get("type") == "Glycosylation"
            ]
    motifs = predicted_motifs(cleaned, annotated_glycans)
    return Candidate(
        label=label,
        sequence=cleaned,
        properties=properties,
        motifs=motifs,
        uniprot=uniprot,
        warnings=list(warnings),
    )


def _normalised(candidate: Candidate) -> dict[str, float]:
    props = candidate.properties
    counts = _counts(candidate)
    gravy = float(props["gravy_kyte_doolittle"])
    return {
        "instability_index_guruprasad": max(0.0, min(1.0, float(props["instability_index_guruprasad"]) / 100.0)),
        "gravy_abs": max(0.0, min(1.0, abs(gravy) / 2.5)),
        "n_glycan_sequons": min(1.0, counts["n_glycan_sequons"] / 10.0),
        "deamidation_sites": min(1.0, counts["deamidation_sites"] / 10.0),
        "free_cysteines": min(1.0, counts["free_cysteines"] / 4.0),
        "long_hydrophobic_stretch": min(1.0, counts["long_hydrophobic_stretch"] / 30.0),
    }


def rank_key(candidate: Candidate) -> float:
    """Weighted designability risk; higher means more likely to be troublesome.

    Every term is bounded to 0-1 before weighting so a single extreme cannot dominate, and
    the weights are public in RANK_WEIGHTS so a designer can substitute their own.
    """
    normalised = _normalised(candidate)
    return round(
        sum(RANK_WEIGHTS[column] * value for column, value in normalised.items()), 4
    )


def screen(
    records: list[tuple[str, str]],
    *,
    lookup_identity: bool = True,
    workers: int = 4,
    supplied: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Analyse and rank every record; failures are reported per row, never fatal."""
    if len(records) > MAX_RECORDS:
        raise SequenceError(f"{len(records)} records exceed the {MAX_RECORDS}-record batch limit")
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []

    supplied = dict(supplied or {})

    def work(item: tuple[str, str]) -> Candidate:
        label, sequence = item
        return analyse_record(
            label,
            sequence,
            lookup_identity=lookup_identity and label not in supplied,
            supplied_uniprot=supplied.get(label),
        )

    with ThreadPoolExecutor(max_workers=max(1, min(workers, 8))) as pool:
        futures = [(label, sequence, pool.submit(work, (label, sequence)))
                   for label, sequence in records]
        for label, _sequence, future in futures:
            try:
                candidate = future.result()
            except (SequenceError, ValueError, OSError) as exc:
                failures.append({"label": label, "reason": str(exc)})
                continue
            props = candidate.properties
            counts = _counts(candidate)
            rows.append({
                "label": candidate.label,
                "accession": candidate.accession or None,
                "identity": candidate.uniprot.get("status"),
                "length": props["length_residues"],
                "rank_score": rank_key(candidate),
                "rank_components": _normalised(candidate),
                "counts": counts,
                "pI": props["isoelectric_point"],
                "net_charge_pH_7": props["net_charge_pH_7"],
                "gravy": props["gravy_kyte_doolittle"],
                "instability_index": props["instability_index_guruprasad"],
                "molecular_weight_da": props["molecular_weight_da"],
                "sequon_positions": candidate.motifs["raw_sequon_positions"],
                "warnings": list(candidate.warnings),
                "deamidation": [
                    {"position": item["position"], "motif": item["motif"]}
                    for item in candidate.motifs.get("deamidation") or []
                ],
                "annotations": sorted({
                    str(item.get("type")) for item in (candidate.uniprot.get("features") or [])
                }) if candidate.verified else [],
            })
    rows.sort(key=lambda row: row["rank_score"])
    for position, row in enumerate(rows, start=1):
        row["rank"] = position
    return {
        "rows": rows,
        "failures": failures,
        "ranking": {
            "weights": RANK_WEIGHTS,
            "direction": "ascending; rank 1 is the lowest predicted design risk",
            "fast_deamidation_multiplier": FAST_DEAMIDATION_MULTIPLIER,
            "note": (
                "a screening order over sequence-derived estimates, not a prediction that any "
                "candidate will or will not express, fold, or be modified"
            ),
        },
        "identity_lookup": lookup_identity,
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }


IDENTITY_LABELS: dict[str, str] = {
    "verified": "matched",
    "ambiguous": "multiple",
    "not_found": "none",
    "unavailable": "lookup failed",
    "mismatch": "mismatch",
    "inactive": "inactive",
    "invalid": "invalid acc",
    "not_requested": "not looked up",
}


def _identity_display(row: dict[str, Any]) -> str:
    """Render identity so a failed lookup never reads as a biological negative.

    "none" means UniParc held no exact match. "lookup failed" means the service did not
    answer, which is not a statement about the protein at all — collapsing the two would
    tell a designer to discard a candidate for the wrong reason.
    """
    status = str(row.get("identity") or "not_requested")
    label = IDENTITY_LABELS.get(status, status)
    if row.get("counts", {}).get("annotated_glycans"):
        return f"{label} \u00b7 {row['counts']['annotated_glycans']} annotated glycan(s)"
    if status == "verified":
        return f"{label} \u00b7 none annotated"
    if status in {"unavailable", "ambiguous"}:
        return f"{label} \u00b7 glycans unknown"
    return label


def _cell(value: Any) -> str:
    return "<td>" + html.escape("" if value is None else str(value)) + "</td>"


COLUMNS: tuple[tuple[str, str], ...] = (
    ("rank", "#"),
    ("label", "Candidate"),
    ("accession", "UniProt"),
    ("identity", "Identity"),
    ("length", "Len"),
    ("counts.n_glycan_sequons", "N-seq*"),
    ("counts.annotated_glycans", "Glyc†"),
    ("counts.deamidation_sites", "Deam"),
    ("counts.free_cysteines", "Cys‡"),
    ("pI", "pI"),
    ("net_charge_pH_7", "q@7"),
    ("gravy", "GRAVY"),
    ("instability_index", "Instab"),
    ("rank_score", "Risk"),
)


def _dig(row: dict[str, Any], path: str) -> Any:
    current: Any = row
    for part in path.split("."):
        if not isinstance(current, dict):
            return None
        current = current.get(part)
    return current


def _risk_breakdown(row: dict[str, Any]) -> str:
    """Show the weighted terms behind the Risk column so the ordering can be audited.

    Without this the score is a black box on the page: the components exist in
    batch.json but a designer reading the table cannot see which term put a candidate
    at the bottom, which is the point of ranking rather than listing.
    """
    components = row.get("rank_components") or {}
    if not components:
        return ""
    parts = sorted(
        (components.get(name, 0.0) * weight, name) for name, weight in RANK_WEIGHTS.items()
    )
    shown = "; ".join(f"{name} {value:.3f}" for value, name in reversed(parts) if value > 0)
    if not shown:
        return ""
    return " &nbsp;·&nbsp; <span class='k'>Risk from</span> " + html.escape(shown)


def accessions_to_records(
    accessions: list[str], *, skip_identity: bool = False
) -> list[tuple[str, str, dict[str, Any]]]:
    """Resolve accessions to canonical sequences so a named protein can be screened.

    The sequence comes from the fetched entry, never from recall. An accession that does
    not resolve raises with its reason rather than being dropped, so a silent shortfall in
    the candidate list is not mistaken for a screened-and-rejected one.
    """
    from .uniprot_report import fetch_annotations

    records: list[tuple[str, str]] = []
    for raw in accessions:
        acc = raw.strip().upper()
        if not acc:
            continue
        if skip_identity:
            raise SequenceError(f"{acc}: cannot resolve an accession with identity lookup off")
        result = fetch_annotations("", acc, accession_precedence=True)
        sequence = result.get("canonical_sequence")
        if result.get("status") != "verified" or not sequence:
            raise SequenceError(
                f"{acc}: {result.get('reason') or result.get('status')}; "
                "no canonical sequence returned"
            )
        organism = (result.get("organism_details") or {}).get("scientificName") or "organism not returned"
        records.append((f"{acc} ({organism})", str(sequence), result))
    return records


def _sites(row: dict[str, Any]) -> str:
    """Per-candidate site positions, grouped by provenance so the classes never merge."""
    sequons = ", ".join(str(item) for item in (row.get("sequon_positions") or [])[:40]) or "none"
    deamidation = ", ".join(
        f"{item['motif']}@{item['position']}"
        for item in (row.get("deamidation") or [])[:40]
    ) or "none"
    annotations = ", ".join(row.get("annotations") or []) or (
        "no verified accession; no annotation attached"
    )
    return (
        "<span class='k'>N-glycan sequons (predicted)</span> " + html.escape(sequons)
        + " &nbsp;·&nbsp; <span class='k'>Deamidation (motif)</span> " + html.escape(deamidation)
        + " &nbsp;·&nbsp; <span class='k'>Annotated feature types</span> "
        + html.escape(annotations)
        + _risk_breakdown(row)
        + (
            " &nbsp;·&nbsp; <span class='k warn'>Input warnings</span> "
            + html.escape("; ".join(row.get("warnings") or []))
            if row.get("warnings")
            else ""
        )
    )


def render_batch_html(data: dict[str, Any], *, source: str) -> str:
    """Render the screen table; only escaped computed values reach the page."""
    head = "".join(f"<th>{html.escape(title)}</th>" for _key, title in COLUMNS)
    body = []
    for row in data["rows"]:
        cells = []
        for key, _title in COLUMNS:
            value = _identity_display(row) if key == "identity" else _dig(row, key)
            cells.append(_cell(value))
        body.append("<tr>" + "".join(cells) + "</tr>")
        body.append(f"<tr class='detail'><td></td><td colspan='{len(COLUMNS) - 1}'>{_sites(row)}</td></tr>")
    failures = "".join(
        f"<li><strong>{html.escape(item['label'])}</strong>: {html.escape(item['reason'])}</li>"
        for item in data["failures"]
    )
    weights = ", ".join(f"{key} {value:+.2f}" for key, value in RANK_WEIGHTS.items())
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Batch designability screen</title><style>
:root{{--ink:#172033;--muted:#64748b;--line:#dbe4ee;--brand:#0e7490}}*{{box-sizing:border-box}}body{{margin:0;font:1rem/1.55 system-ui;color:var(--ink);background:#f4f7f9}}.shell{{width:min(1180px,calc(100% - 2rem));margin:auto;padding-block:2rem 4rem}}h1{{font-size:clamp(1.7rem,5vw,2.6rem);letter-spacing:-.03em;margin:0 0 .4rem}}p{{max-width:90ch;color:var(--muted)}}table{{display:table;width:100%;border-collapse:separate;border-spacing:0;border:1px solid var(--line);border-radius:14px;background:#fff;box-shadow:0 7px 24px rgba(15,23,42,.045);font-variant-numeric:tabular-nums}}th,td{{padding:.5rem .6rem;border-bottom:1px solid var(--line);text-align:left;font-size:.9rem}}th{{background:#f1f5f9;font-weight:700;position:sticky;top:0}}tr:last-child td{{border-bottom:0}}td:first-child{{color:var(--brand);font-weight:700}}tr.detail td{{background:#f8fafc;color:var(--muted);font-size:.82rem;padding-top:.1rem}}tr.detail .k{{color:var(--ink);font-weight:650}}tr.detail .warn{{color:#9f1239}}.legend{{margin-top:1rem;font-size:.85rem;color:var(--muted)}}.legend code{{background:#e2e8f0;padding:.05rem .3rem;border-radius:5px}}ul{{font-size:.9rem}}
</style></head><body><div class="shell"><h1>Batch designability screen</h1>
<p>{len(data['rows'])} candidate(s) analysed from {html.escape(source)} at {html.escape(data['generated_at_utc'])} UTC. Rank 1 carries the lowest weighted sequence-derived risk; the order is a screening aid, not a prediction that any candidate will express, fold, or carry a modification.</p>
<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>
<p class="legend"><strong>*</strong> predicted N-X-[ST] sequon — occupancy not established, no evidence behind it.
<strong>†</strong> glycosylation annotated in the verified UniProt entry — only populated for an exact canonical-sequence match.
<strong>‡</strong> cysteines not accounted for by an annotated disulfide, or all cysteines when no identity was verified.
<br>Risk weights: <code>{html.escape(weights)}</code>. Fast (N-G) deamidation counts double.</p>
{f"<h2>Not analysed</h2><ul>{failures}</ul>" if failures else ""}
<div class="legend"><p><strong>What this screen cannot tell you:</strong> every column is computed from sequence. None is a measurement. A candidate with no verified accession has no annotation attached at all, and absence of a glycan annotation is not absence of a glycan.</p></div>
</div></body></html>
"""


def run(
    input_path: Path | None,
    output_dir: Path,
    *,
    no_identity: bool = False,
    accession: list[str] | None = None,
) -> Path:
    """Screen FASTA records and/or accessions into one HTML table plus JSON."""
    lookup_identity = not no_identity
    records = read_multi_fasta(input_path) if input_path is not None else []
    resolved = accessions_to_records(list(accession or []), skip_identity=no_identity)
    records = records + [(lab, seq) for lab, seq, _rec in resolved]
    supplied = {lab: rec for lab, _seq, rec in resolved}
    if not records:
        raise SequenceError("provide either --input FASTA records or --accession identifiers")
    data = screen(records, lookup_identity=lookup_identity, supplied=supplied)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "batch.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    report = output_dir / "batch.html"
    report.write_text(render_batch_html(data, source=str(input_path)), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Screen many candidate sequences for design-relevant properties"
    )
    parser.add_argument("--input", type=Path, help="multi-record FASTA")
    parser.add_argument("--output", type=Path, required=True, help="output directory")
    parser.add_argument(
        "--no-identity",
        action="store_true",
        help="skip UniProt lookup; pure sequence-only screening, faster and offline-safe",
    )
    parser.add_argument(
        "--accession",
        action="append",
        default=[],
        help="UniProt accession to screen; repeatable, and combined with --input records",
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    try:
        report = run(
            args.input, args.output, no_identity=args.no_identity, accession=args.accession
        )
    except (SequenceError, OSError, UnicodeError, ValueError) as exc:
        print(f"Batch screen not generated: {exc}", file=sys.stderr)
        return 2
    print(f"Batch screen generated: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
