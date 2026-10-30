"""The report shows only bounded PDBe mappings for a verified UniProt accession."""

import json

from protein_agent.pdbe_mcp import PdbeMcpError
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


def test_pdbe_failure_reason_is_actionable() -> None:
    class Broken(FakePdbe):
        def tool_schemas(self):
            raise PdbeMcpError("server executable was not found")

    record = fetch_pdbe_structures("Q12345", client=Broken())
    assert record["status"] == "unavailable"
    assert "server executable was not found" in record["reason"]


def test_oversized_ranked_mapping_uses_bounded_search_fallback() -> None:
    class Oversized(FakePdbe):
        def tool_schemas(self):
            return [
                *super().tool_schemas(),
                {
                    "type": "function",
                    "function": {
                        "name": "live_schema_name",
                        "description": "Retrieves the Solr search schema for PDBe.",
                        "parameters": {"type": "object", "properties": {}},
                    },
                },
                {
                    "type": "function",
                    "function": {
                        "name": "live_query_name",
                        "description": (
                            "Executes a search query against the PDBe Solr search service."
                        ),
                        "parameters": {
                            "type": "object",
                            "properties": {"query": {}},
                            "required": ["query"],
                        },
                    },
                },
            ]

        def call(self, name, arguments):
            if name == "generated_live_name":
                raise PdbeMcpError(
                    "PDBe MCP result was 724704 bytes, above the configured 262144-byte limit"
                )
            if name == "live_schema_name":
                assert arguments == {}
                return {"content": [{"type": "text", "text": "schema"}]}
            assert name == "live_query_name"
            assert arguments["rows"] == 10
            assert arguments["query"] == "uniprot_accession:Q12345"
            assert "pdb_id" in arguments["fl"]
            return {
                "content": [{
                    "type": "text",
                    "text": (
                        "Results metadata:\nNumber of documents found: 1\n"
                        "Start index: 0\nDocuments:\nDocument 1:\n"
                        "  chain_id: A\n  entity_id: 1\n"
                        "  experimental_method: X-ray diffraction\n"
                        "  overall_quality: 82.0\n  pdb_id: 8t82\n"
                        "  prefered_assembly_id: 1\n  resolution: 1.1\n"
                        "  uniprot_coverage: 0.25"
                    ),
                }],
                "_provenance": {"source": "PDBe MCP", "server": "pdbe_search_server"},
            }

    record = fetch_pdbe_structures("Q12345", client=Oversized())
    assert record["status"] == "verified"
    assert record["structures"][0]["pdb_id"] == "8T82"
    assert record["structures"][0]["coverage_fraction"] == 0.25
    assert "overall quality" in record["ranking"]


def test_pdbe_wrapped_404_is_a_bounded_no_match() -> None:
    class NoData(FakePdbe):
        def call(self, name, arguments):
            return {
                "content": [{
                    "type": "text",
                    "text": (
                        "API request failed: 404 Client Error\nStatus Code: 404\n"
                        'Response: {"message":"Requested endpoint does not contain any data"}'
                    ),
                }],
                "isError": False,
            }

    record = fetch_pdbe_structures("Q12345", client=NoData())
    assert record["status"] == "not_found"
    assert "under this UniProt accession" in record["reason"]


def test_pdbe_empty_mapping_is_not_reported_as_global_absence() -> None:
    class Empty(FakePdbe):
        def call(self, name, arguments):
            return {"content": [{"type": "text", "text": '{"Q12345": []}'}]}

    record = fetch_pdbe_structures("Q12345", client=Empty())
    assert record["status"] == "not_found"
    assert "under this accession" in record["reason"]
