"""Local computation + the chromophore claim limits, anchored on a retrieved fixture.

Fixture: green fluorescent protein P42212 (Aequorea victoria). Selected by the assistant
as a test case, not requested by the user; accession and annotations verified by
UniProtClient against rest.uniprot.org in the test run. Kept honest by the facts that
fixture revealed: wtGFP carries Ser65, so TYG/G[FYL]G cores miss the real chromophore.
"""

from __future__ import annotations

import json

import httpx
import pytest
import respx

from protein_agent.config import Settings
from protein_agent.properties import (
    INVALIDATED_MOTIFS,
    SequenceError,
    chromophore_candidates,
    clean_fasta_sequence,
    compute_properties,
)
from protein_agent.retrieval import UniProtClient


def test_composition_values_are_sane() -> None:
    props = compute_properties("ACDEFGHIKLMNPQRSTVWY", "std20").as_dict()
    assert props["length_residues"] == 20
    total = round(sum(props["amino_acid_percent"].values()), 1)
    assert total == pytest.approx(100.0, abs=0.5), "composition must be percent summing to 100"
    assert props["molecular_weight_da"] > 2000
    # extinction is reported as a reduced/oxidised pair, never a single figure
    ext = props["molar_extinction_coefficient_M-1cm-1"]
    assert set(ext) == {"reduced", "oxidised"}
    for value in props["method"].values():
        assert "ProteinAnalysis" in value or "Kyte" in value, "every number names its method"


def test_every_number_names_its_scale() -> None:
    props = compute_properties("MTMTRYSTALK", "x").as_dict()
    for key in ("molecular_weight", "isoelectric_point", "gravy", "instability_index"):
        assert any(key.split("_")[0] in v for v in props["method"].values())
    assert props["caveats"], "computed values must ship with their limitations"


def test_hydropathy_is_a_profile_not_a_claim() -> None:
    props = compute_properties("IIIIIVVVVVKKKKK", "hp").as_dict()
    lo_pos, lo = props["hydropathy_extremes"]["min"]
    hi_pos, hi = props["hydropathy_extremes"]["max"]
    assert lo < hi and 1 <= lo_pos <= 16 and 1 <= hi_pos <= 16


def test_input_is_never_silently_repaired() -> None:
    seq, warnings = clean_fasta_sequence("MTMTRY-STAL\nKQ")
    assert "-" not in seq
    assert any("gap" in w for w in warnings), "gaps must be reported, not swallowed"
    _, w2 = clean_fasta_sequence("MTMXRJZQ")
    assert any("non-standard" in w for w in w2)
    with pytest.raises(SequenceError):
        clean_fasta_sequence("12345")


def test_chromophore_detector_only_nominates() -> None:
    result = chromophore_candidates("MTSYGVQAPLAGYQ")
    assert result["status"] == "candidates_found"
    assert result["complete_result"] is True
    assert "does not establish" in result["meaning"].lower()
    assert result["classification"].startswith("heuristic")
    assert any("uniprot" in step.lower() for step in result["to_confirm"])


def test_no_candidates_is_a_terminal_result_not_a_retry_prompt() -> None:
    """An empty scan must stop the loop, not invite a mutated re-call."""
    result = chromophore_candidates("MMMMMMKKKKKKDDDDDD")
    assert result["status"] == "no_candidates"
    assert result["complete_result"] is True
    assert "do NOT retry" in result["verdict"]
    # absence of the motif is not absence of photo-activity
    assert "NOT" in result["verdict"] and "fails to absorb" in result["verdict"]


def test_invalidated_motifs_stay_invalidated() -> None:
    """The cores that failed against the retrieved fixture must never come back."""
    assert "TYG" in INVALIDATED_MOTIFS and "G[FYL]G" in INVALIDATED_MOTIFS
    result = chromophore_candidates("MTSYGVQAPLAGYQ")
    assert not any("TYG" in k for k in result["candidates"])


@respx.mock
def test_absorption_is_read_from_the_entry_not_inferred() -> None:
    entry = {
        "sequence": {"value": "MTMTRY", "length": 6},
        "comments": [
            {
                "commentType": "BIOPHYSICOCHEMICAL PROPERTIES",
                "absorption": {
                    "max": 395,
                    "approximate": False,
                    "texts": [{"value": "peaks at 395 nm"}],
                    "evidences": [
                        {
                            "evidenceCode": "ECO:0000269",
                            "source": "PubMed",
                            "id": "8448132",
                        }
                    ],
                },
            }
        ],
        "features": [
            {
                "type": "Cross-link",
                "location": {"start": {"value": 65}, "end": {"value": 67}},
                "description": "5-imidazolinone (Ser-Gly)",
                "evidences": [
                    {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "8448132"}
                ],
            }
        ],
    }
    respx.get("https://uniprot.test/uniprotkb/P42212.json").mock(
        return_value=httpx.Response(200, json=entry)
    )
    client = UniProtClient(Settings(uniprot_base_url="https://uniprot.test", api_key="k"))
    result = client.entry("P42212")
    photo = result["photochemistry"]
    assert photo["absorption"]["max_nm"] == 395
    assert photo["absorption"]["evidence"][0]["code"] == "ECO:0000269"
    assert photo["ptm_processing"][0]["type"] == "Cross-link"


@respx.mock
def test_missing_absorption_is_not_annotated_never_absent() -> None:
    respx.get("https://uniprot.test/uniprotkb/P12345.json").mock(
        return_value=httpx.Response(200, json={"sequence": {"value": "AAA", "length": 3}})
    )
    client = UniProtClient(Settings(uniprot_base_url="https://uniprot.test", api_key="k"))
    result = client.entry("P12345")
    photo = result["photochemistry"]
    assert photo["absorption"] is None
    assert "not evidence" in photo["absorption_note"]


# --- pepstats-equivalence (constants read from EMBOSS source, not recalled) ---------


def test_pepstats_equivalent_fields_present() -> None:
    props = compute_properties(
        "MTMTRYSTALKSTALEGCNHNSQITWQFYVPKAGDLLLLLLLLIVVVVVIKKKDDDIIIIIIIIIIIIILLLLLL",
        "x",
    ).as_dict()
    for key in (
        "average_residue_weight_da",
        "charge_vs_ph",
        "amino_acid_dayhoff_stats",
        "physico_chemical_classes",
        "inclusion_body_propensity",
        "a280_extinction_per_mg_ml_M-1cm-1",
    ):
        assert props.get(key), f"missing pepstats-equivalent field: {key}"


def test_charge_table_spans_pH_and_brackets_the_pI() -> None:
    from protein_agent.pepstats import charge_table

    table = charge_table("MTMTRYSTALKSTALEGCNHNSQITWQFYVPKAGDLLLLLLLL")
    charges = table["ph_to_net_charge"]
    assert "3.0" in charges and "11.0" in charges
    # charge must decrease monotonically as pH rises, and cross zero near the pI
    values = [charges[k] for k in sorted(charges, key=float)]
    assert values == sorted(values, reverse=True), "charge must fall as pH rises"
    assert any(v <= 0 for v in values) and any(v >= 0 for v in values)


def test_dayhoff_stat_is_mole_percent_over_reference_frequency() -> None:
    from protein_agent.pepstats import DAYHOFF_FREQUENCY, dayhoff_stats

    seq = "AAAAA" + "C" * 0  # 5 Ala out of 5 residues = 100 mol%
    stats = dayhoff_stats(seq)
    assert stats["A"]["number"] == 5
    assert stats["A"]["mole_percent"] == pytest.approx(100.0)
    assert stats["A"]["dayhoff_stat"] == pytest.approx(100.0 / DAYHOFF_FREQUENCY["A"], abs=1e-3)
    assert stats["W"]["number"] == 0 and stats["W"]["dayhoff_stat"] == 0.0


def test_inclusion_body_probability_is_reported_not_quoted_when_invalid() -> None:
    """The EMBOSS quadratic is a parabola: outside cv in [0, 2*CVDISC] the value is
    meaningless as a probability and a negative one can occur. Direction must come from
    `class`, which is what EMBOSS itself branches on."""
    from protein_agent.pepstats import inclusion_body_propensity

    extreme = inclusion_body_propensity("NNNNGGGGPPPPSSSS" * 3 + "KKKDDD")
    assert extreme["probability_inclusion_bodies"] < 0.0  # impossible as a probability
    assert extreme["probability_is_meaningful"] is False
    assert extreme["read_the_class_not_the_probability"] is True
    assert extreme["class"] == "likely_insoluble"  # direction still correct

    mild = inclusion_body_propensity("MTMTRYSTALKSTALEGCNHNSQITWQFYVPKAGDL")
    assert 0.0 <= mild["probability_inclusion_bodies"] <= 1.0
    assert mild["probability_is_meaningful"] is True


def test_pepstats_caveats_are_carried_through() -> None:
    props = compute_properties("MTMTRYSTALK", "x").as_dict()
    text = json.dumps(props)
    assert "EMBOSS pepstats.c" in text or "LAM1" in text, "provenance of the constants"
    assert props["inclusion_body_propensity"]["caveats"]
