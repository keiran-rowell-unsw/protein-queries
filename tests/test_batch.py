"""Batch designability screen: ranking honesty, provenance separation, row isolation.

Fixtures are assistant-selected test cases, not requested by the user:
  * P05067 (APP1, human, 770 aa) and P14735 (phytochrome A) identities, accessions and
    annotated features verified by fetch_annotations against rest.uniprot.org in an
    earlier run of this session.
  * A NOVEL sequence that is assistant-authored and is NOT a real protein; it exists only
    to exercise the identity-free path a generated binder actually takes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from protein_agent.batch import (
    MAX_RECORDS,
    RANK_WEIGHTS,
    SequenceError,
    _counts,
    _identity_display,
    accessions_to_records,
    analyse_record,
    rank_key,
    read_multi_fasta,
    render_batch_html,
    run,
    screen,
)

NOVEL = (
    "AEFAELYQPWSPKVTIYQGFIRNGETAQVGNLTAQGVWSFAEVGQPVKNLDPQAENWVAHGA"
    "QYKALQGFNEQTPAKVLDSAIRYGEAHPKGTVVEAALKELGEK"
)


@pytest.fixture()
def novel():
    return analyse_record("novel", NOVEL, lookup_identity=False)


def test_multi_fasta_reads_every_record(tmp_path: Path) -> None:
    path = tmp_path / "cands.fa"
    path.write_text(">a\nMKVL\n>b\nACDE\n\n>c\nAGLT\n", encoding="utf-8")
    assert read_multi_fasta(path) == [("a", "MKVL"), ("b", "ACDE"), ("c", "AGLT")]


def test_empty_file_is_a_finding_not_a_zero_row(tmp_path: Path) -> None:
    path = tmp_path / "empty.fa"
    path.write_text("\n \n", encoding="utf-8")
    with pytest.raises(Exception, match="no FASTA records"):
        read_multi_fasta(path)


def test_counts_are_integers_not_floats(novel) -> None:
    # A screening table showing "2.0" deamidation sites implies a precision it does not
    # have, and breaks consumers that treat the column as a count.
    counts = _counts(novel)
    assert all(isinstance(value, int) for value in counts.values()), counts


def test_ambiguity_codes_are_stripped_not_silently_counted() -> None:
    # clean_fasta_sequence strips non-standard letters; a stripped record must not be
    # reported as a real protein of the original length, and the stripping must be visible.
    candidate = analyse_record("odd", "ACDXBZJKL" + "A" * 40, lookup_identity=False)
    assert len(candidate.sequence) < 49
    assert candidate.properties["length_residues"] == len(candidate.sequence)
    assert candidate.warnings, "lossy input must be surfaced, never dropped"
    assert "non-standard residues removed" in candidate.warnings[0]


def test_clean_sequence_carries_no_spurious_warning(novel) -> None:
    assert novel.warnings == []


def test_warning_reaches_the_page_and_the_json(tmp_path: Path) -> None:
    data = screen([("odd", "ACDXBZJKL" * 6)], lookup_identity=False)
    assert data["rows"][0]["warnings"], "JSON must carry it for downstream tools"
    page = render_batch_html(data, source="in.fa")
    assert "Input warnings" in page and "non-standard residues removed" in page


def test_novel_sequence_gets_design_stats_without_any_identity(novel) -> None:
    # The core promise: a generated binder with no database record still screens.
    assert novel.uniprot["status"] == "not_requested"
    assert novel.motifs["raw_sequon_positions"], "sequon scan must not depend on identity"
    assert novel.properties["isoelectric_point"] > 0


def test_predicted_and_annotated_glycans_are_never_merged(novel) -> None:
    counts = _counts(novel)
    assert counts["annotated_glycans"] == 0, "no identity means no annotation, ever"
    assert counts["n_glycan_sequons"] >= 1
    # The columns are distinct: summing them would launder a motif hit into a curated fact.
    assert "n_glycan_sequons" in RANK_WEIGHTS and "annotated_glycans" not in RANK_WEIGHTS


def test_fast_deamidation_is_weighted_heavier_than_slow() -> None:
    slow_only = analyse_record("s", "MSNQK" * 12, lookup_identity=False)
    fast = analyse_record("f", "MNGQK" * 12, lookup_identity=False)
    assert _counts(fast)["deamidation_sites"] > _counts(slow_only)["deamidation_sites"]


def test_free_cysteines_use_annotated_disulfides_only_when_verified() -> None:
    # Without verification nothing may be subtracted: assuming pairing would under-report
    # free cysteine, the single most actionable aggregation signal in this table.
    candidate = analyse_record("c", "ACCKACCPACCPAACCKACCPACCPAAC", lookup_identity=False)
    assert candidate.uniprot["status"] == "not_requested"
    assert _counts(candidate)["free_cysteines"] == _counts(candidate)["cysteines"]


def test_ranking_is_ascending_and_weights_are_public(novel) -> None:
    assert abs(sum(RANK_WEIGHTS.values()) - 1.0) < 1e-9
    assert 0.0 <= rank_key(novel) <= 1.0
    assert rank_key(novel) == rank_key(analyse_record("same", NOVEL, lookup_identity=False))


def test_extreme_hydrophobe_outranks_balanced_sequence() -> None:
    balanced = analyse_record("bal", NOVEL, lookup_identity=False)
    hydrophobic = analyse_record("hydro", "AVLIVLFAILVFAILVFAILV" * 6, lookup_identity=False)
    assert rank_key(hydrophobic) > rank_key(balanced)


def test_normalisation_caps_a_single_extreme_term() -> None:
    # Every term is bounded to 0-1 before weighting, so one runaway column cannot dominate.
    huge = analyse_record("huge", "AVLIVLFAILVFAILVFAILV" * 40, lookup_identity=False)
    assert max(_counts(huge).values()) > 100
    assert rank_key(huge) <= 1.0


def test_identity_labels_distinguish_failure_from_absence() -> None:
    def label(status, glycans=0):
        return _identity_display({"identity": status, "counts": {"annotated_glycans": glycans}})

    assert "none" in label("not_found")
    assert "lookup failed" in label("unavailable")
    assert "multiple" in label("ambiguous")
    assert "matched" in label("verified", 0)
    assert "2 annotated glycan" in label("verified", 2)
    # A failed lookup must never render as a biological negative.
    assert "none" not in label("unavailable")


def test_screen_orders_by_risk_and_assigns_consecutive_ranks() -> None:
    data = screen([("easy", NOVEL), ("hard", "AVLIVLFAILVFAILVFAILV" * 8)], lookup_identity=False)
    assert [row["rank"] for row in data["rows"]] == [1, 2]
    assert data["rows"][0]["rank_score"] <= data["rows"][1]["rank_score"]


def test_one_bad_row_does_not_abort_the_batch() -> None:
    data = screen([("good", NOVEL), ("bad", "12345"), ("also", NOVEL)], lookup_identity=False)
    assert len(data["rows"]) == 2
    assert len(data["failures"]) == 1
    assert data["failures"][0]["label"] == "bad"


def test_batch_limit_is_enforced(tmp_path: Path) -> None:
    with pytest.raises(Exception, match="exceed"):
        screen([("x", "AAA") for x in range(MAX_RECORDS + 1)], lookup_identity=False)


def test_html_escapes_a_hostile_header_and_keeps_provenance_footnotes() -> None:
    hostile = '<img src=x onerror="alert(1)">'
    data = screen([(hostile, NOVEL)], lookup_identity=False)
    page = render_batch_html(data, source="in.fa")
    assert "<img src=x" not in page
    assert "&lt;img src=x" in page
    # The footnotes that make the table readable are not optional.
    assert "N-X-[ST]" in page and "not a prediction" in page
    assert "predicted" in page and "annotated" in page


def test_ranking_metadata_is_machine_readable() -> None:
    data = screen([("a", NOVEL)], lookup_identity=False)
    assert data["ranking"]["direction"].startswith("ascending")
    assert set(data["ranking"]["weights"]) == set(RANK_WEIGHTS)
    assert "not a prediction" in data["ranking"]["note"]


def test_cli_round_trip_writes_json_and_html(tmp_path: Path) -> None:
    from protein_agent.batch import run

    source = tmp_path / "in.fa"
    source.write_text(f">one\n{NOVEL}\n>two\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")
    report = run(source, tmp_path / "out", no_identity=True)
    assert report.is_file()
    payload = json.loads((tmp_path / "out" / "batch.json").read_text())
    assert {row["label"] for row in payload["rows"]} == {"one", "two"}
    assert payload["identity_lookup"] is False


# --- accession input path (added after a bug the sequence-only tests could not catch) ---


def test_supplied_uniprot_record_is_not_second_guessed_by_checksum() -> None:
    # A user-named accession is authoritative. Re-deriving identity by checksum turned
    # P05067 into "ambiguous" with a null accession and threw away its 10 annotated
    # glycans -- the exact inversion of what the caller asked for.
    verified = {
        "status": "verified",
        "accession": "P05067",
        "features": [{"type": "Glycosylation", "start": 542, "end": 542}],
    }
    candidate = analyse_record(
        "P05067 (Homo sapiens)", NOVEL, lookup_identity=True, supplied_uniprot=verified
    )
    assert candidate.uniprot["status"] == "verified"
    assert candidate.accession == "P05067"
    assert _counts(candidate)["annotated_glycans"] == 1


def test_supplied_record_suppresses_the_redundant_lookup(monkeypatch) -> None:
    called = []
    import protein_agent.uniprot_report as ur

    monkeypatch.setattr(
        ur, "find_sequence_annotations", lambda *a, **k: called.append(1) or {"status": "not_found"}
    )
    analyse_record("x", NOVEL, lookup_identity=True, supplied_uniprot={"status": "verified"})
    assert called == [], "an authoritative record must not trigger a second identity search"


def test_accession_resolution_uses_fetched_sequence_not_recall(monkeypatch) -> None:
    import protein_agent.uniprot_report as ur

    fetched = {"status": "verified", "accession": "P1", "canonical_sequence": "ACDEF",
               "organism_details": {"scientificName": "Testus"}, "features": []}
    monkeypatch.setattr(ur, "fetch_annotations", lambda *a, **k: dict(fetched))
    records = accessions_to_records(["p1"])
    assert records == [("P1 (Testus)", "ACDEF", fetched)]


def test_unresolvable_accession_raises_rather_than_disappearing(monkeypatch) -> None:
    import protein_agent.uniprot_report as ur

    monkeypatch.setattr(
        ur, "fetch_annotations",
        lambda *a, **k: {"status": "not_found", "reason": "No UniProt entry found."},
    )
    with pytest.raises(SequenceError, match="ZZZZZZ"):
        accessions_to_records(["ZZZZZZ"])


def test_accession_requires_identity_lookup(tmp_path: Path) -> None:
    with pytest.raises(SequenceError, match="identity lookup off"):
        run(None, tmp_path / "out", no_identity=True, accession=["P1"])


def test_run_accepts_accession_without_any_fasta(monkeypatch, tmp_path: Path) -> None:
    import protein_agent.uniprot_report as ur

    rec = {"status": "verified", "accession": "P1", "canonical_sequence": "ACDEF" * 6,
           "organism_details": {"scientificName": "Testus"},
           "features": [{"type": "Glycosylation", "start": 4, "end": 4}]}
    monkeypatch.setattr(ur, "fetch_annotations", lambda *a, **k: dict(rec))
    report = run(None, tmp_path / "out", accession=["P1"])
    assert report.is_file()
    payload = json.loads((tmp_path / "out" / "batch.json").read_text())
    assert payload["rows"][0]["accession"] == "P1"
    assert payload["rows"][0]["counts"]["annotated_glycans"] == 1


def test_no_input_at_all_is_a_clear_error(tmp_path: Path) -> None:
    with pytest.raises(SequenceError, match="provide either"):
        run(None, tmp_path / "out")


def test_risk_breakdown_is_visible_on_the_page() -> None:
    # The score is only auditable if the terms that produced it are rendered, not just
    # present in JSON.
    data = screen([("dominant", "AVLIVLFAILVFAILVFAILV" * 10)], lookup_identity=False)
    page = render_batch_html(data, source="in.fa")
    assert "Risk from" in page
    assert "long_hydrophobic_stretch" in page or "gravy_abs" in page
