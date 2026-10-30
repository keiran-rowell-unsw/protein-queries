"""Which structure gets drawn, and the refusal to draw a misleading one.

These tests exist because ranking by resolution — the obvious implementation — is wrong.
The fixtures encode a live observation on APP1 P05067 (770 aa, verified via
fetch_pdbe_structures in this session): the ranked fallback returned 8T82 at 1.1 Å with
0.01 coverage and 5ONR at 1.39 Å with 0.0 coverage. A resolution-first sort renders those;
a coverage-first sort declines to render anything.
"""

from __future__ import annotations

import pytest

from protein_agent.structure_choice import (
    CONSTRUCT_DISPLAY_THRESHOLD,
    FULL_DISPLAY_THRESHOLD,
    choose_display_structure,
    coverage_basis,
    effective_coverage,
    is_crystallographic,
)

# The real APP1 fallback set: best resolution belongs to the smallest fragment.
APP1_FRAGMENTS = [
    {"pdb_id": "8T82", "chain_id": "A", "resolution_angstrom": 1.1, "coverage_fraction": 0.01},
    {"pdb_id": "3PZZ", "chain_id": "B, A", "resolution_angstrom": 1.29, "coverage_fraction": 0.01},
    {"pdb_id": "5ONR", "chain_id": "B", "resolution_angstrom": 1.39, "coverage_fraction": 0.0},
    {"pdb_id": "2FMA", "chain_id": "A", "resolution_angstrom": 0.85, "coverage_fraction": 0.08},
]


def test_fragment_only_set_refuses_to_render() -> None:
    result = choose_display_structure(APP1_FRAGMENTS, sequence_length=770)
    assert result["status"] == "fragment_only"
    assert result["selected"] is None
    assert "no structure is rendered" in result["selection_reason"]
    assert "2FMA 8%" in result["selection_reason"]  # names the near-miss it declined


def test_resolution_is_a_tiebreaker_not_the_primary_sort() -> None:
    full = {
        "pdb_id": "1AAA",
        "resolution_angstrom": 2.9,
        "coverage_fraction": 0.97,
        "observed_regions": [{"unp_start": 11, "unp_end": 768}],
    }
    fragment = {"pdb_id": "1BBB", "resolution_angstrom": 1.1, "coverage_fraction": 0.01}
    result = choose_display_structure([fragment, full], sequence_length=770)
    assert result["status"] == "full"
    assert result["selected"]["pdb_id"] == "1AAA"  # the 2.9 A model, not the 1.1 A fragment
    assert "highest-resolution mapping covering" in result["selection_reason"]


def test_resolution_ordering_within_a_coverage_tier() -> None:
    a = {"pdb_id": "1CCC", "resolution_angstrom": 2.5, "coverage_fraction": 0.95}
    b = {"pdb_id": "1DDD", "resolution_angstrom": 1.6, "coverage_fraction": 0.95}
    assert choose_display_structure([a, b], sequence_length=200)["selected"]["pdb_id"] == "1DDD"


def test_construct_is_labelled_rather_than_rejected() -> None:
    only = {
        "pdb_id": "1EEE",
        "resolution_angstrom": 1.8,
        "coverage_fraction": 0.5,
        "uniprot_start": 100,
        "uniprot_end": 485,
        "observed_regions": [{"unp_start": 100, "unprot_end": 485}],
    }
    result = choose_display_structure([only], sequence_length=770)
    assert result["status"] == "partial"
    assert "construct or domain" in result["selection_reason"]
    assert "not the full protein" in result["selection_reason"]
    assert result["selected"]["pdb_id"] == "1EEE"


def test_observed_segments_beat_a_generous_reported_fraction() -> None:
    # A construct claiming 90% coverage that only models 20 residues must not be drawn as
    # the full protein; the residue-level mapping is the tighter, honest number.
    structure = {
        "pdb_id": "1FFF",
        "resolution_angstrom": 2.0,
        "coverage_fraction": 0.9,
        "observed_regions": [{"unp_start": 10, "unp_end": 29}],
    }
    assert effective_coverage(structure, 770) == pytest.approx(20 / 770)
    assert coverage_basis(structure, 770) == "PDBe SIFTS observed segments"
    assert choose_display_structure([structure], sequence_length=770)["status"] == "fragment_only"


def test_missing_coverage_never_wins_by_default() -> None:
    unknown = {"pdb_id": "1GGG", "resolution_angstrom": 1.0}
    known = {"pdb_id": "1HHH", "resolution_angstrom": 3.4, "coverage_fraction": 0.92}
    result = choose_display_structure([unknown, known], sequence_length=300)
    assert result["selected"]["pdb_id"] == "1HHH"  # unverified does not outrank verified
    assert result["status"] == "full"


def test_mapping_range_is_used_when_no_fraction_or_segments() -> None:
    structure = {
        "pdb_id": "1HHH",
        "uniprot_start": 1,
        "uniprot_end": 290,
        "resolution_angstrom": 2.0,
    }
    assert effective_coverage(structure, 300) == pytest.approx(290 / 300)
    assert coverage_basis(structure, 300) == "PDBe UniProt mapping range"


def test_bad_resolution_values_are_discarded_without_touching_coverage() -> None:
    # The guard must apply to resolution alone: 0.01 coverage is the decisive measurement
    # this module exists to act on, and must survive parsing.
    assert effective_coverage({"coverage_fraction": 0.01}, 770) == 0.01
    assert effective_coverage({"coverage_fraction": -3}, 770) is None
    assert effective_coverage({"coverage_fraction": True}, 770) is None
    assert effective_coverage({"coverage_fraction": 4.2}, 770) == 1.0  # clamped, not trusted
    assert is_crystallographic({"experimental_method": "X-ray diffraction"}) is True
    assert is_crystallographic({"experimental_method": "solution NMR"}) is False


def test_alternative_structures_are_counted_not_dropped() -> None:
    result = choose_display_structure(APP1_FRAGMENTS, sequence_length=770, preferred_pdb_id="2FMA")
    assert result["selected"] is None  # user-named entry cannot override the fragment guard
    assert result["status"] == "fragment_only"
    # Declining a named entry must still say it was the named entry, so the user learns
    # their request was seen and why it was refused, rather than seeing a generic failure.
    assert "user-requested entry" in result["selection_reason"]
    chosen = choose_display_structure(APP1_FRAGMENTS, sequence_length=770, preferred_pdb_id="8T82")
    assert chosen["alternative_count"] == 3
    assert set(result["thresholds"]) == {"full", "construct"}


def test_thresholds_are_explicit_in_the_result() -> None:
    result = choose_display_structure(
        [{"pdb_id": "1JJJ", "coverage_fraction": 0.95}], sequence_length=100
    )
    assert result["thresholds"] == {
        "full": FULL_DISPLAY_THRESHOLD,
        "construct": CONSTRUCT_DISPLAY_THRESHOLD,
    }


def test_nonsense_inputs_are_handled_as_states_not_crashes() -> None:
    assert choose_display_structure([], sequence_length=100)["status"] == "none"
    assert choose_display_structure([{"chain_id": "A"}], sequence_length=100)["status"] == "none"
    assert choose_display_structure([None, 3], sequence_length=100)["status"] == "none"  # type: ignore[list-item]
    with pytest.raises(ValueError, match="positive residue count"):
        choose_display_structure(APP1_FRAGMENTS, sequence_length=0)


def test_selection_is_deterministic_for_equal_candidates() -> None:
    tied = [
        {"pdb_id": "1ZZZ", "coverage_fraction": 0.9},
        {"pdb_id": "1AAA", "coverage_fraction": 0.9},
    ]
    first = choose_display_structure(tied, sequence_length=100)["selected"]["pdb_id"]
    second = choose_display_structure(list(reversed(tied)), sequence_length=100)["selected"][
        "pdb_id"
    ]
    assert first == second == "1AAA"
