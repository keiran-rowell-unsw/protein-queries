"""The report must never show a structure it did not earn.

These cover the wiring between the single-sequence report and the analysis modules:
what gets drawn when PDBe yields a good entry, and what is withheld — and how it is
worded — when it does not. The distinction that matters is between "no structure exists"
and "we could not look", which a reader cannot see unless the wording carries it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from protein_agent import report_sections as rs
from protein_agent.cartoon import parse_coordinates

# A committed subset of PDBe entry 1M5C (retinoid-binding receptor, RBP4): C-alpha atoms
# plus the modelled BRH ligand, waters dropped. Biopython reads only CA for the trace and
# only non-water HETATMs for the ligand census, so this parses identically to the deposit
# while staying small enough to commit. It exists so the suite needs no network and no
# gitignored scratch directory.
COORDINATES = Path(__file__).parent / "data" / "1m5c_ca.pdb"

FULL_MAPPING = {
    "pdb_id": "1M5C",
    "chain_id": "A",
    "resolution_angstrom": 1.9,
    "experimental_method": "X-ray diffraction",
    "uniprot_start": 1,
    "uniprot_end": 236,
    "coverage_fraction": 0.97,
    "source_url": "https://www.ebi.ac.uk/pdbe/entry/pdb/1m5c",
}

# The observed APP1 fallback set: best resolution belongs to the smallest fragment.
FRAGMENT_MAPPINGS = [
    {"pdb_id": "8T82", "chain_id": "A", "resolution_angstrom": 1.1, "coverage_fraction": 0.01},
    {"pdb_id": "2FMA", "chain_id": "A", "resolution_angstrom": 0.85, "coverage_fraction": 0.08},
]

SEQUENCE = (
    "MVPSKDKDINKIYGEENVTCQMNASFGELSPRTTFGTTVKNPPDHTYTSPHAFHISPGNKDLKPENLLLK"
    "KSSGGGLTDPVGGAADGGNKKCQQWQNVTSKPSRSPSKPKDAPQQNQNSSPRTVSQINEVLSKKGSDSS"
    "HPYRQVDYYDPDPKKQDKEPWSYKTEPSSSTNKQPTPQSTLQKQVFGNPSQGLPSKWIPDQTHQGLSSM"
)


def test_the_committed_coordinate_fixture_parses() -> None:
    """If this fails, the structure tests below would silently stop exercising the renderer."""
    assert COORDINATES.is_file(), f"missing committed fixture {COORDINATES}"
    parsed = parse_coordinates(COORDINATES.read_bytes())
    assert parsed["chains"] and parsed["ligands"], "fixture lost its CA trace or its ligand"


@pytest.fixture
def no_network(monkeypatch):
    """A report must not reach the network from these tests; coordinates come from disk."""
    def fake(pdb_id: str) -> bytes:
        return COORDINATES.read_bytes()

    monkeypatch.setattr(rs, "download_coordinates", fake)
    return fake


def _data(structures: list[dict], length: int) -> dict:
    return {"length_residues": length, "pdbe": {"status": "verified", "structures": structures}}


# ---------------------------------------------------------------- choice and rendering


def test_a_covering_entry_is_chosen_and_rendered(tmp_path, no_network) -> None:
    display = rs.build_structure_display(_data([FULL_MAPPING], 236), tmp_path)
    assert display["status"] == "full"
    assert display["pdb_id"] == "1M5C"
    assert display["render"]["status"] == "rendered"
    assert (tmp_path / "structure.png").is_file()
    assert (tmp_path / "ligands.png").is_file()


def test_the_cartoon_survives_a_missing_pymol(tmp_path, no_network) -> None:
    """PyMOL is optional; the CA trace is what makes the section useful without it."""
    display = rs.build_structure_display(_data([FULL_MAPPING], 236), tmp_path)
    assert display["cartoon_is_fallback"] is True
    assert display["pymol_render"]["status"] == "unavailable"
    # The caveats that matter belong to the renderer that actually produced the image.
    caveats = " ".join(display["render"]["caveats"])
    assert "CA trace only" in caveats
    assert "unmodelled" in caveats


def test_why_an_entry_was_picked_is_stated(tmp_path, no_network) -> None:
    display = rs.build_structure_display(_data([FULL_MAPPING], 236), tmp_path)
    assert "coverage" in display["selection_reason"].lower()
    html = rs.render_structure_section({"structure_display": display})
    assert "Why this entry" in html
    assert "1M5C" in html


def test_resolution_alone_does_not_win_a_slot_in_the_report(tmp_path, no_network) -> None:
    """A 0.85 Å fragment must not displace a 1.9 Å near-complete model."""
    display = rs.build_structure_display(
        _data([FULL_MAPPING, *FRAGMENT_MAPPINGS], 236), tmp_path
    )
    assert display["pdb_id"] == "1M5C"


def test_all_fragments_decline_to_render_but_say_why(tmp_path, no_network) -> None:
    """Structures exist and were declined; that is not the same as having none.

    The wording is the substance of this test. "No structure was found" would tell a
    reader the protein is undiscovered, when in fact ten entries exist for it.
    """
    display = rs.build_structure_display(_data(FRAGMENT_MAPPINGS, 770), tmp_path)
    assert display["status"] == "fragment_only"
    assert display.get("pdb_id") is None
    assert not (tmp_path / "structure.png").exists()

    html = rs.render_structure_section({"structure_display": display})
    assert "<img" not in html
    assert "2FMA" in html and "8T82" in html  # the near-misses are named
    assert "no structure is rendered" in html
    assert "under 30%" in html


@pytest.mark.parametrize(
    ("pdbe_status", "must_not_say"),
    [
        # A failed lookup must never read as an absence of experimental data.
        ("unavailable", "has been deposited"),
        # Nor may "we did not ask" read as "there is nothing to find".
        ("not_requested", "No experimental structure has been deposited"),
    ],
)
def test_a_lookup_gap_is_not_reported_as_absence_of_data(
    pdbe_status: str, must_not_say: str
) -> None:
    display = {
        "status": pdbe_status,
        "reason": "the recorded reason from the upstream call",
        "pdb_id": "1M5C",
    }
    html = rs.render_structure_section({"structure_display": display})
    assert "<img" not in html
    assert must_not_say not in html
    assert "the recorded reason from the upstream call" in html


def test_not_found_says_the_database_has_nothing_under_this_accession() -> None:
    html = rs.render_structure_section(
        {"structure_display": {"status": "not_found", "reason": "none returned"}}
    )
    assert "deposited in PDBe under this accession" in html
    assert "<img" not in html


@pytest.mark.parametrize("status", [None, "", "   ", "weird", 7, ["full"], "not_a_state"])
def test_an_unrecognised_status_cannot_borrow_the_success_template(status: object) -> None:
    """The dangerous default. Anything unknown must fail closed to "not shown".

    A whitelist matters more here than a blacklist: if the ranking ever grows a status,
    the report should withhold an image rather than draw a confident empty structure box.
    """
    html = rs.render_structure_section(
        {"structure_display": {
            "status": status, "pdb_id": "1M5C", "render": {"status": "rendered"},
        }}
    )
    assert "<img" not in html
    assert "not necessarily the best" not in html


@pytest.mark.parametrize("identifier", ["", "   ", "1M5", "1M5CC", "<script>", None, 1234])
def test_a_malformed_identifier_is_never_rendered_as_an_entry(identifier: object) -> None:
    html = rs.render_structure_section(
        {"structure_display": {
            "status": "full", "pdb_id": identifier, "render": {"status": "rendered"},
        }}
    )
    assert "<img" not in html
    assert "<script>" not in html


def test_oversized_entries_are_named_but_not_drawn(tmp_path, no_network) -> None:
    display = rs.build_structure_display(_data([FULL_MAPPING], 9_000), tmp_path)
    assert display["pdb_id"] == "1M5C"  # the entry is still identified
    assert display["render"]["status"] == "too_large"
    assert not (tmp_path / "structure.png").exists()

    html = rs.render_structure_section({"structure_display": display})
    assert "1M5C" in html
    # Naming an entry while withholding its image must not also claim nothing was found.
    assert "No usable PDBe mapping was returned" not in html


def test_coordinate_failure_is_reported_rather_than_blanked(tmp_path, monkeypatch) -> None:
    def boom(pdb_id: str) -> bytes:
        raise ConnectionError("PDBe refused")

    monkeypatch.setattr(rs, "download_coordinates", boom)
    display = rs.build_structure_display(_data([FULL_MAPPING], 236), tmp_path)
    assert display["render"]["status"] == "unavailable"
    assert "ConnectionError" in display["render"]["reason"]
    assert "1M5C" in rs.render_structure_section({"structure_display": display})


# ------------------------------------------------------------------ motifs / position map


def test_motifs_are_reconciled_against_fetched_annotation_only() -> None:
    verified = {"status": "verified", "features": [
        {"type": "Glycosylation", "start": 13, "end": 13, "description": "N-linked"}]}
    motifs = rs.build_motifs(SEQUENCE, verified)
    assert motifs["status"] == "computed"
    # Position 13 sits in an N-Q sequon: reconciliation reads the position from the payload.
    assert any(item["position"] == 13 for item in motifs["n_glycan_sequons"])
    assert "retrieved entry" in motifs["annotation_reconciliation"] or "annotation" in motifs[
        "annotation_reconciliation"
    ]


def test_without_identity_no_predicted_site_is_called_confirmed() -> None:
    """An unverified sequence gets no annotation, and must not be shown any."""
    motifs = rs.build_motifs(SEQUENCE, {"status": "not_found", "reason": "no match"})
    assert motifs["counts"]["confirmed"] == 0
    assert "No verified UniProt accession" in motifs["annotation_reconciliation"]


def test_an_absent_feature_list_does_not_become_an_empty_annotation_track(tmp_path) -> None:
    """UniProt returning nothing and UniProt not being asked look identical on a plot.

    The annotated tracks are cleared and the map says it shows predictions only.
    """
    uniprot = {"status": "unavailable", "reason": "UniProt could not be reached"}
    motifs = rs.build_motifs(SEQUENCE, uniprot)
    record = rs.build_position_map(SEQUENCE, uniprot, motifs, tmp_path / "positions.png")
    assert record["status"] == "rendered"
    for name in rs.ANNOTATED_TRACKS:
        assert not record["tracks"].get(name), f"{name} must be empty without a verified entry"
    assert any("predictions only" in caveat for caveat in record["caveats"])


def test_annotated_tracks_appear_when_the_entry_was_retrieved(tmp_path) -> None:
    uniprot = {"status": "verified", "features": [
        {"type": "Glycosylation", "start": 13, "end": 13, "description": "N-linked (GlcNAc)"},
        {"type": "Domain", "start": 30, "end": 90, "description": "Ig-like"}]}
    motifs = rs.build_motifs(SEQUENCE, uniprot)
    record = rs.build_position_map(SEQUENCE, uniprot, motifs, tmp_path / "positions.png")
    assert record["tracks"]["Glycosylation (annotated)"] == 1
    assert record["tracks"]["Domain"] == 1
    assert (tmp_path / "positions.png").is_file()


def test_section_renderers_survive_a_partial_record() -> None:
    """Every section is rendered from a report that may be missing any one of them."""
    assert "Not computed" in rs.render_motif_section({})
    assert "Not shown" in rs.render_position_map_section({})
    assert "Not shown" in rs.render_structure_section({})


def test_rendered_sections_carry_their_qualifiers(tmp_path) -> None:
    uniprot = {"status": "verified", "features": [
        {"type": "Glycosylation", "start": 13, "end": 13, "description": "N-linked"}]}
    motifs = rs.build_motifs(SEQUENCE, uniprot)
    position_map = rs.build_position_map(SEQUENCE, uniprot, motifs, tmp_path / "positions.png")
    data = {"predicted_motifs": motifs, "position_map": position_map}
    motif_html = rs.render_motif_section(data)
    assert "sequence-pattern predictions" in motif_html
    assert "occupancy" in motif_html
    map_html = rs.render_position_map_section(data)
    assert "positions.png" in map_html


# --------------------------------------------------------------------- serialisation


def test_added_records_are_json_serialisable_for_report_json(tmp_path, no_network) -> None:
    """report.json is written from the same dict, so nothing here may hold a Path or a model."""
    uniprot = {"status": "verified", "features": [
        {"type": "Glycosylation", "start": 13, "end": 13, "description": "N-linked"}]}
    motifs = rs.build_motifs(SEQUENCE, uniprot)
    payload = {
        "predicted_motifs": motifs,
        "position_map": rs.build_position_map(
            SEQUENCE, uniprot, motifs, tmp_path / "positions.png"
        ),
        "structure_display": rs.build_structure_display(_data([FULL_MAPPING], 236), tmp_path),
    }
    json.dumps(payload)  # raises TypeError on anything unserialisable
