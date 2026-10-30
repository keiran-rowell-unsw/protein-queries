"""Predicted-motif tests: sequon truth, provenance separation, and claim limits.

Fixtures are sequence-derived only, so no network or inference call is needed: a pattern
match is a fact about the string, and these tests hold the module to that. UniProt
positions quoted in comments come from fetch_annotations against rest.uniprot.org in an
earlier run of this session (APP1 P05067, 770 aa).
"""

from __future__ import annotations

from protein_agent.sequons import (
    CLASSIFICATION,
    cleavage_sites,
    compare_to_annotations,
    deamidation_sites,
    n_glycan_sequons,
    predicted_motifs,
    termini_flags,
)

# Assistant-selected test string (not requested by the user): built to contain, in order,
# a canonical sequon, a Pro-excluded sequon, an N-G deamidation, and a K/R-before-P site.
DESIGN = "MSNALLKRGHELONGTAVCFNGTLHPQNPSGTKRY"  # assistant-selected test string, not a real protein


def test_sequon_requires_non_proline_x_position() -> None:
    # Positions: 1 N-A-T (sequon) | 5 N-B-C is not | 8 N-P-S excluded because Pro cannot
    # adopt the conformation the oligotransferase requires.
    assert [item["position"] for item in n_glycan_sequons("NATXXNBCNXXNPST")] == [1]
    # X may itself be N, C, G, S, or T; only Pro is excluded.
    assert [
        item["position"] for item in n_glycan_sequons("NNT NCT NGT NST NTT".replace(" ", ""))
    ] == [1, 4, 7, 10, 13]


def test_proline_blocked_sequon_is_absent_not_flagged() -> None:
    # N-P-S would be a sequon if Pro were permitted; it must simply not appear.
    assert n_glycan_sequons("NPST") == []
    assert n_glycan_sequons("NAST") and n_glycan_sequons("NAST")[0]["position"] == 1


def test_sequon_spelling_reports_actual_residues() -> None:
    site = n_glycan_sequons("NAT")[0]
    assert site["sequon"] == "N1-A-T"
    assert site["classification"] == "motif"


def test_deamidation_classes_separate_fast_and_slow() -> None:
    items = deamidation_sites("NNGNST")
    assert [(item["position"], item["motif"]) for item in items] == [(2, "N-G"), (4, "N-S")]
    assert [item["rate_class"] for item in items] == ["fast", "slow"]


def test_slow_motif_label_names_the_residue_that_matched() -> None:
    # The slow class requires Ser/Thr immediately after the Asn. A mislabelled liability is
    # worse than a missing one: a designer would mutate a residue that was never at risk.
    # NGNAT: N1-G2 is fast; N3 is followed by A, so it is neither class and must not appear.
    assert [(i["position"], i["motif"]) for i in deamidation_sites("NGNAT")] == [(1, "N-G")]

    # Ser/Thr at +1 is the slow class; Gly at +1 is fast; anything else at +1 is no site.
    def triples(seq):
        return [(i["position"], i["motif"], i["rate_class"]) for i in deamidation_sites(seq)]

    assert triples("NST") == [(1, "N-S", "slow")]
    assert triples("NGS") == [(1, "N-G", "fast")]
    assert triples("NAS") == []  # Ala at +1 blocks both classes
    assert triples("NAT") == []  # Ser/Thr must be at +1, not merely present at +2
    # NNST: residue 2 is Asn followed by Ser. Ser is the slow class only; the fast class is
    # Gly-specific, so exactly one entry is reported. Residue 1 sees N at +1: no class.
    assert [(i["position"], i["motif"], i["rate_class"]) for i in deamidation_sites("NNST")] == [
        (2, "N-S", "slow")
    ]
    assert all(item["motif"] in {"N-G", "N-S", "N-T"} for item in deamidation_sites("NNGNSTNAT"))


def test_cleavage_sites_respect_proline_exclusion() -> None:
    trypsin = cleavage_sites("KRAPKRG", "trypsin")
    # K1 and R2 are followed by A, so both cleave; A3-P4 blocks; K5-R6-G7 both cleave.
    assert trypsin["positions"] == [1, 2, 5, 6]
    # K immediately before P must never appear.
    assert 3 not in trypsin["positions"]
    # A terminal K/R has no following residue, so the not-before-P rule cannot exclude it.
    assert cleavage_sites("KPG", "trypsin")["positions"] == [1]
    # R4 precedes P5? No: R2 is followed by A3 and cleaves; the K5-R6 pair is blocked only
    # where the next residue is Pro, so R6 before P is excluded while K5 is kept.
    assert cleavage_sites("KRP", "trypsin")["positions"] == [1, 2]
    assert trypsin["classification"] == "motif"
    assert "surface-accessible" in trypsin["caveat"]


def test_chymotrypsin_tracks_aromatics_not_basic_residues() -> None:
    sites = cleavage_sites("KFWAPY", "chymotrypsin")
    # F2 and W3 cleave; K is not a chymotrypsin site; Y6 before end-of-chain cleaves.
    assert sites["positions"] == [2, 3, 6]
    assert sites["specificity"] == "F/Y/W not followed by P"


def test_unsupported_protease_is_rejected_not_ignored() -> None:
    try:
        cleavage_sites("ACDE", "proteinase_k")
    except ValueError as exc:
        assert "proteinase_k" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("unknown protease must raise")


def test_reconciliation_keeps_annotated_and_predicted_apart() -> None:
    sequons = n_glycan_sequons("NATXNNST")  # positions 1 and 5-ish; second is also sequonic
    annotated = [{"type": "Glycosylation", "start": 1, "end": 1, "description": "N-linked"}]
    reconciled = compare_to_annotations(sequons, annotated)
    states = {item["position"]: item["annotation_state"] for item in reconciled}
    assert states[1] == "confirmed"
    motif_only = [item for item in reconciled if item["classification"] == "motif"]
    assert any(item["annotation_state"] == "predicted_only" for item in motif_only)
    # Predicted sites never arrive labelled as curated.
    assert all(item["classification"] != "curated" for item in motif_only)


def test_annotated_site_without_sequon_is_surfaced_as_suspicious() -> None:
    # A fetched entry asserting glycosylation at a residue whose sequon is absent from the
    # analysed sequence is an identity/numbering flag, not a prediction to be averaged away.
    result = compare_to_annotations([], [{"type": "Glycosylation", "start": 50, "end": 50}])
    assert len(result) == 1
    assert result[0]["annotation_state"] == "annotated_without_sequon"
    assert result[0]["classification"] == "curated"
    assert "isoform" in result[0]["annotation_detail"]


def test_predicted_motifs_carries_claim_limits() -> None:
    payload = predicted_motifs(DESIGN, [])
    assert payload["status"] == "computed"
    assert payload["classification"] == "motif"
    assert payload["classification_note"] == CLASSIFICATION
    assert any("O-linked" in note for note in payload["known_limitations"])
    # Absence of a sequon must not be phrased as absence of glycosylation.
    empty = predicted_motifs("ACDEFGHIKLM", [])
    assert empty["raw_sequon_positions"] == []
    assert "NOT" in empty["verdict"]


def test_termini_flags_do_not_claim_tag_tolerance() -> None:
    flags = termini_flags("PMQ")
    assert flags == {
        "n_term": "P",
        "c_term": "Q",
        "n_term_proline": True,
        "classification": "computed",
        "caveat": flags["caveat"],
    }
    assert "tag tolerance" in flags["caveat"]


def test_empty_sequence_is_skipped_not_invented() -> None:
    assert predicted_motifs("")["status"] == "skipped"
