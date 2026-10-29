"""The report shows only bounded PDBe mappings for a verified UniProt accession."""

import json

from protein_agent.pdbe_report import fetch_pdbe_structures
from protein_agent.report import _pdbe_section


class FakePdbe:
    def tool_schemas(self):
        return [{
            "type": "function",
            "function": {
                "name": "generated_live_name",
                "description": "Provides a ranked list of PDB structures mapped to UniProt.",
                "parameters": {
                    "type": "object",
                    "properties": {"live_accession_field": {}},
                    "required": ["live_accession_field"],
                },
            },
        }]

    def call(self, name, arguments):
        assert name == "generated_live_name"
        assert arguments == {"live_accession_field": "Q12345"}
        payload = {"Q12345": [{
            "pdb_id": "1abc",
            "chain_id": "A",
            "entity_id": 1,
            "experimental_method": "X-ray diffraction",
            "resolution": 1.8,
            "coverage": 0.75,
            "unp_start": 10,
            "unp_end": 39,
            "observed_regions": [{"unp_start": 10, "unp_end": 39}],
            "preferred_assembly_id": 1,
        }]}
        return {
            "content": [{"type": "text", "text": json.dumps(payload)}],
            "isError": False,
            "_provenance": {"source": "PDBe MCP", "server": "pdbe_api_server"},
        }

    def close(self):
        raise AssertionError("injected clients are not owned")


def test_verified_pdbe_mapping_is_bounded_and_rendered() -> None:
    record = fetch_pdbe_structures("Q12345", client=FakePdbe())
    assert record["status"] == "verified"
    assert record["structures"][0]["pdb_id"] == "1ABC"
    page = _pdbe_section(record)
    assert "1ABC" in page
    assert "1.8 Å" in page
    assert "75% coverage" in page
    assert "https://www.ebi.ac.uk/pdbe/entry/pdb/1abc" in page


def test_pdbe_empty_mapping_is_not_reported_as_global_absence() -> None:
    class Empty(FakePdbe):
        def call(self, name, arguments):
            return {"content": [{"type": "text", "text": '{"Q12345": []}'}]}

    record = fetch_pdbe_structures("Q12345", client=Empty())
    assert record["status"] == "not_found"
    assert "under this accession" in record["reason"]
