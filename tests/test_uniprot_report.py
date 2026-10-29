"""UniProt fields reach a report only after an exact canonical sequence match."""

import httpx
import respx

from protein_agent.config import Settings
from protein_agent.report import _uniprot_section
from protein_agent.retrieval import UniProtClient
from protein_agent.uniprot_report import fetch_annotations

SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 2
ENTRY = {
    "primaryAccession": "Q12345",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "organism": {"scientificName": "Example organism"},
    "entryVersion": 12,
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
