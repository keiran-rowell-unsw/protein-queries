"""PDBe tools are discovered from MCP and invoked with bounded, provenanced calls."""

from __future__ import annotations

from typing import Any

import pytest

from protein_agent.config import Settings
from protein_agent.pdbe_mcp import SEARCH_SERVER, PdbeMcpClient, PdbeMcpError
from protein_agent.proof import Evidence, unsupported_claims
from protein_agent.retrieval import PdbeClient, UniProtClient, dispatch
from protein_agent.router import build_plan


def _tool(name: str, description: str, properties: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": name,
        "description": description,
        "inputSchema": {
            "type": "object",
            "properties": properties,
            "required": ["query"] if "query" in properties else [],
            "additionalProperties": False,
        },
    }


class FakeRunner:
    def __init__(self) -> None:
        self.calls: list[tuple[str, str, str | None, dict[str, Any] | None]] = []

    def __call__(
        self,
        server: str,
        operation: str,
        name: str | None,
        arguments: dict[str, Any] | None,
    ) -> Any:
        self.calls.append((server, operation, name, arguments))
        if operation == "list":
            if server == SEARCH_SERVER:
                return [
                    _tool("schema_live", "search schema", {}),
                    _tool(
                        "search_live",
                        "bounded PDBe structure search",
                        {"query": {"type": "string"}, "rows": {"type": "integer"}},
                    ),
                ]
            return [_tool("entry_live_hash", "PDB entry summary", {"pdb_id": {}})]
        return {"content": [{"type": "text", "text": "PDB 1ABC"}], "isError": False}


def test_discovery_preserves_server_schemas_and_routes_calls() -> None:
    runner = FakeRunner()
    mcp = PdbeMcpClient(Settings(), runner=runner)
    schemas = mcp.openai_tools()
    assert {item["function"]["name"] for item in schemas} == {
        "schema_live", "search_live", "entry_live_hash"
    }
    assert len([call for call in runner.calls if call[1] == "list"]) == 2

    result = mcp.call("search_live", {"query": "uniprot_accession:P12345", "rows": 500})
    assert runner.calls[-1] == (
        SEARCH_SERVER,
        "call",
        "search_live",
        {"query": "uniprot_accession:P12345", "rows": 25},
    )
    assert result["_provenance"]["server"] == SEARCH_SERVER
    assert result["_provenance"]["tool"] == "search_live"


def test_fielded_search_options_require_schema_call_first() -> None:
    runner = FakeRunner()
    mcp = PdbeMcpClient(Settings(), runner=runner)
    mcp.discover()
    with pytest.raises(PdbeMcpError, match="schema tool"):
        mcp.call("search_live", {"query": "*:*", "fl": ["pdb_id"]})
    mcp.call("schema_live", {})
    mcp.call("search_live", {"query": "*:*", "fl": ["pdb_id"]})


def test_unknown_or_duplicate_tools_are_rejected() -> None:
    mcp = PdbeMcpClient(Settings(), runner=FakeRunner())
    with pytest.raises(PdbeMcpError, match="did not advertise"):
        mcp.call("remembered_tool_name", {})

    def duplicate(server: str, operation: str, name: str | None, args: dict | None):
        return [_tool("same_name", server, {})]

    with pytest.raises(PdbeMcpError, match="duplicated"):
        PdbeMcpClient(Settings(), runner=duplicate).discover()


def test_dispatch_uses_discovered_name_not_pdbe_prefix() -> None:
    runner = FakeRunner()
    pdbe = PdbeClient(Settings(), mcp=PdbeMcpClient(Settings(), runner=runner))
    uniprot = UniProtClient(Settings())
    try:
        result = dispatch(
            "entry_live_hash", {"pdb_id": "1ABC"}, uniprot=uniprot, pdbe=pdbe
        )
    finally:
        uniprot.close()
    assert result["_provenance"]["tool"] == "entry_live_hash"


def test_failed_result_cannot_ground_an_identifier() -> None:
    evidence = Evidence()
    evidence.add_call("entry_live_hash", {"pdb_id": "1ABC"})
    evidence.add_result(
        '{"isError":true,"content":"No PDB 1ABC",'
        '"_provenance":{"source":"PDBe MCP","server":"pdbe_api_server"}}'
    )
    problems = unsupported_claims("PDB 1ABC has a structure.", evidence)
    assert any("successful tool" in item for item in problems)
    assert any("never returned" in item for item in problems)


def test_structure_details_need_api_confirmation_not_search_only() -> None:
    evidence = Evidence()
    evidence.add_call("search_live", {"query": "pdb_id:1ABC"})
    evidence.add_result(
        '{"content":"PDB 1ABC, resolution 2.0",'
        '"_provenance":{"source":"PDBe MCP","server":"pdbe_search_server"}}'
    )
    problems = unsupported_claims(
        "PDB 1ABC has a resolution of 2.0 Å and contains a ligand.", evidence
    )
    assert any("PDBe API MCP" in item for item in problems)

    evidence.add_call("entry_live_hash", {"pdb_id": "1ABC"})
    evidence.add_result(
        '{"content":"PDB 1ABC, resolution 2.0, ligand", "isError":false,'
        '"_provenance":{"source":"PDBe MCP","server":"pdbe_api_server"}}'
    )
    assert not unsupported_claims(
        "PDB 1ABC has a resolution of 2.0 Å and contains a ligand.", evidence
    )


def test_router_advertises_all_live_pdbe_tools_only_for_structure_questions() -> None:
    live = PdbeMcpClient(Settings(), runner=FakeRunner()).openai_tools()
    structure = build_plan("What structures and ligands exist for P12345?", pdbe_tools=live)
    assert "pdbe" in structure.services
    assert set(structure.extra_tool_names) == {
        "schema_live", "search_live", "entry_live_hash"
    }

    annotation = build_plan("What is the function of P12345?", pdbe_tools=live)
    assert "pdbe" not in annotation.services
    assert annotation.extra_tool_names == ()

    direct = build_plan("What is the resolution of PDB 1ABC?", pdbe_tools=live)
    assert direct.services == ("pdbe",)
    assert "fetch_uniprot_entry" not in direct.tool_names
    assert "search_uniprotkb" not in direct.tool_names
