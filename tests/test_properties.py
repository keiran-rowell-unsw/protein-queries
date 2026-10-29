"""Local computation + the chromophore claim limits, anchored on a retrieved fixture.

Fixture: green fluorescent protein P42212 (Aequorea victoria). Selected by the assistant
as a test case, not requested by the user; accession and annotations verified by
UniProtClient against rest.uniprot.org in the test run. Kept honest by the facts that
fixture revealed: wtGFP carries Ser65, so TYG/G[FYL]G cores miss the real chromophore.
"""

from __future__ import annotations

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
