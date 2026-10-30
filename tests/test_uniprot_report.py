"""UniProt fields reach a report only after an exact canonical sequence match."""

import httpx
import respx
from Bio.SeqUtils.CheckSum import crc64

from protein_agent.config import Settings
from protein_agent.report import _uniprot_section
from protein_agent.retrieval import UniProtClient
from protein_agent.uniprot_report import fetch_annotations, find_sequence_annotations

SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 2
ENTRY = {
    "primaryAccession": "Q12345",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "organism": {"scientificName": "Example organism"},
    "entryVersion": 12,
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Example protein"}}},
    "genes": [{"geneName": {"value": "EXAMPLE1"}}],
    "sequence": {"value": SEQUENCE, "length": len(SEQUENCE)},
    "features": [
        {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 5}},
         "description": "Signal peptide", "evidences": [{"evidenceCode": "ECO:0000256"}]},
        {"type": "Glycosylation", "location": {"start": {"value": 12}, "end": {"value": 12}},
         "description": "<img src=x onerror=alert(1)>"},
        {"type": "Modified residue", "location": {"start": {"value": 30}, "end": {"value": 30}},
         "description": "Modified site"},
    ],
    "comments": [
        {"commentType": "BIOPHYSICOCHEMICAL PROPERTIES",
         "absorption": {"max": 405, "note": {"texts": [{"value": "Example note"}]}}},
    ],
}


def make_client():
    return UniProtClient(Settings(uniprot_base_url="https://uniprot.test"))


@respx.mock
def test_exact_match_displays_entry_annotations_and_escapes_text():
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(200, json=ENTRY)
    )
    record = fetch_annotations(SEQUENCE, "q12345", client=make_client())
    assert record["status"] == "verified"
    assert record["canonical_sequence_match"] is True
    assert record["accession"] == "Q12345"
    assert len(record["features"]) == 3
    assert record["absorption"]["max_nm"] == 405
    page = _uniprot_section(record)
    assert "Glycosylation" in page and "Residue(s) 12" in page
    assert "Signal" in page and "Modified residue" in page
    assert "&lt;img" in page and "<img src=x" not in page
    assert "https://www.uniprot.org/uniprotkb/Q12345/entry" in page
    assert "ECO:0000256" in page
    assert "Example protein" in page and "EXAMPLE1" in page
    assert record["protein_description"] == ENTRY["proteinDescription"]


@respx.mock
def test_mismatched_sequence_never_receives_annotations():
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(200, json=ENTRY)
    )
    result = fetch_annotations("A" * len(SEQUENCE), "Q12345", client=make_client())
    assert result["status"] == "mismatch"
    assert "features" not in result and "absorption" not in result
    assert "Glycosylation" not in _uniprot_section(result)


@respx.mock
def test_not_found_is_not_a_negative_biological_claim():
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(404)
    )
    result = fetch_annotations(SEQUENCE, "Q12345", client=make_client())
    assert result["status"] == "not_found"
    assert "features" not in result


def test_invalid_accession_never_queries():
    assert fetch_annotations(SEQUENCE, "../etc/passwd")["status"] == "invalid"


@respx.mock
def test_detailed_annotation_sections_render_only_returned_uniprot_fields():
    detailed = {
        **ENTRY,
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "features": [
            *ENTRY["features"],
            {"type": "Domain", "location": {"start": {"value": 4}, "end": {"value": 20}},
             "description": "Catalytic domain"},
            {"type": "Active site", "location": {"start": {"value": 18},
                                                    "end": {"value": 18}},
             "description": "Proton donor", "evidences": [{"evidenceCode": "ECO:0000269",
                                                               "source": "PubMed", "id": "1"}]},
            {"type": "Natural variant", "location": {"start": {"value": 21},
                                                        "end": {"value": 21}},
             "featureId": "VAR_000001", "description": "Example variant",
             "alternativeSequence": {"originalSequence": "A",
                                     "alternativeSequences": ["V"]}},
        ],
        "comments": [
            *ENTRY["comments"],
            {"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [
                {"location": {"value": "Cytoplasm", "evidences": [
                    {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "3"}
                ]}, "topology": {"value": "Peripheral"}}
            ]},
            {"commentType": "MASS SPECTROMETRY", "method": "Electrospray",
             "molWeight": 1234.5, "range": {"start": {"value": 2},
                                              "end": {"value": 39}}},
            {"commentType": "PTM", "texts": [{"value": "Phosphorylated after activation."}]},
            {"commentType": "SUBUNIT", "texts": [{"value": "Forms a homodimer."}]},
            {"commentType": "FUNCTION", "texts": [{"value": "Cleaves an example substrate."}]},
            {"commentType": "CATALYTIC ACTIVITY", "reaction": {
                "name": "A + H2O = B.", "ecNumber": "1.2.3.4",
                "evidences": [{"evidenceCode": "ECO:0000269", "source": "PubMed",
                               "id": "2"}],
            }},
            {"commentType": "INTERACTION", "interactions": [{
                "interactantOne": {"uniProtKBAccession": "Q12345"},
                "interactantTwo": {"uniProtKBAccession": "Q54321", "geneName": "PARTNER",
                                    "intActId": "EBI-1"},
                "numberOfExperiments": 3,
            }]},
        ],
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "1ABC", "properties": [
                {"key": "Method", "value": "X-ray"},
                {"key": "Resolution", "value": "1.5 A"},
                {"key": "Chains", "value": "A=1-40"},
            ]},
            {"database": "IntAct", "id": "Q12345", "properties": [
                {"key": "Interactions", "value": "1"}
            ]},
            {"database": "InterPro", "id": "IPR000001", "properties": []},
        ],
    }
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(200, json=detailed)
    )
    record = fetch_annotations(SEQUENCE, "Q12345", client=make_client())
    page = _uniprot_section(record)
    for expected in (
        "Experimental structures", "1ABC", "X-ray", "Mass-spectrometry protocol",
        "Electrospray", "Localisation", "Cytoplasm; Peripheral; evidence: ECO:0000269",
        "Glycosylation sites", "Post-translational modifications",
        "Phosphorylated after activation", "Stoichiometry / subunit organisation",
        "Forms a homodimer", "Domain boundaries", "Catalytic domain",
        "Functional / active sites", "Proton donor", "Catalytic activity", "EC 1.2.3.4",
        "Protein-protein interactions", "PARTNER", "Variants", "VAR_000001",
        "AlphaMissense classifications were not returned",
    ):
        assert expected in page
    assert len(record["features"]) == len(detailed["features"])
    assert record["organism_details"]["taxonId"] == 9606


@respx.mock
def test_sequence_checksum_search_verifies_candidate_and_entry():
    checksum = crc64(SEQUENCE).removeprefix("CRC-")
    endpoint = respx.get("https://uniprot.test/uniparc/search").mock(
        return_value=httpx.Response(200, json={"results": [{
            "sequence": {"value": SEQUENCE, "length": len(SEQUENCE), "crc64": checksum},
            "uniProtKBAccessions": ["Q12345"],
        }]})
    )
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(200, json=ENTRY)
    )
    result = find_sequence_annotations(SEQUENCE, settings=Settings(uniprot_base_url="https://uniprot.test"))
    assert result["status"] == "verified"
    assert result["accession"] == "Q12345"
    assert result["lookup"]["checksum"] == checksum
    assert endpoint.called
    assert endpoint.calls[0].request.url.params["query"] == f"checksum:{checksum}"
    assert SEQUENCE not in str(endpoint.calls[0].request.url)


@respx.mock
def test_checksum_collision_must_verify_full_sequence():
    respx.get("https://uniprot.test/uniparc/search").mock(
        return_value=httpx.Response(200, json={"results": [{
            "sequence": {"value": "A" * len(SEQUENCE), "length": len(SEQUENCE)},
            "uniProtKBAccessions": ["Q12345"],
        }]})
    )
    result = find_sequence_annotations(SEQUENCE, settings=Settings(uniprot_base_url="https://uniprot.test"))
    assert result["status"] == "not_found"
    assert "accession" not in result


@respx.mock
def test_inactive_accession_does_not_block_unique_current_match():
    respx.get("https://uniprot.test/uniparc/search").mock(
        return_value=httpx.Response(200, json={"results": [{
            "sequence": {"value": SEQUENCE, "length": len(SEQUENCE)},
            "uniProtKBAccessions": ["Q12345", "Q54321.1"],
        }]})
    )
    respx.get("https://uniprot.test/uniprotkb/Q12345.json").mock(
        return_value=httpx.Response(200, json=ENTRY)
    )
    respx.get("https://uniprot.test/uniprotkb/Q54321.json").mock(
        return_value=httpx.Response(303, json={"entryType": "Inactive"})
    )
    result = find_sequence_annotations(SEQUENCE, settings=Settings(uniprot_base_url="https://uniprot.test"))
    assert result["status"] == "verified"
    assert result["accession"] == "Q12345"
    assert fetch_annotations(SEQUENCE, "Q54321", client=make_client())["status"] == "inactive"


@respx.mock
def test_multiple_exact_accessions_are_ambiguous_not_first_hit():
    respx.get("https://uniprot.test/uniparc/search").mock(
        return_value=httpx.Response(200, json={"results": [{
            "sequence": {"value": SEQUENCE, "length": len(SEQUENCE)},
            "uniProtKBAccessions": ["Q12345", "Q54321"],
        }]})
    )
    for acc in ("Q12345", "Q54321"):
        respx.get(f"https://uniprot.test/uniprotkb/{acc}.json").mock(
            return_value=httpx.Response(200, json={**ENTRY, "primaryAccession": acc})
        )
    result = find_sequence_annotations(SEQUENCE, settings=Settings(uniprot_base_url="https://uniprot.test"))
    assert result["status"] == "ambiguous"
    assert set(result["candidates"]) == {"Q12345", "Q54321"}
    assert "features" not in result
    assert "Q12345" in _uniprot_section(result)


@respx.mock
def test_no_checksum_hit_is_not_a_biological_absence():
    respx.get("https://uniprot.test/uniparc/search").mock(
        return_value=httpx.Response(200, json={"results": []})
    )
    result = find_sequence_annotations(SEQUENCE, settings=Settings(uniprot_base_url="https://uniprot.test"))
    assert result["status"] == "not_found"
    assert "identity is unresolved" in result["reason"]
