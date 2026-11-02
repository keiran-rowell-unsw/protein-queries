"""Foldseek auto-reporting only searches a retrieved, PDBe-mapped structure."""

import io

import httpx
import numpy as np
import pytest
from Bio.PDB import MMCIFIO, Atom, Chain, MMCIFParser, Model, Residue, Structure

from protein_agent import report as report_module
from protein_agent.config import Settings
from protein_agent.foldseek import FoldseekClient
from protein_agent.foldseek_query import prepare_query
from protein_agent.report import generate_report, render_foldseek_section
from protein_agent.retrieval import RetrievalError
from protein_agent.router import build_plan
from protein_agent.tools import FOLDSEEK_SEARCH


def synthetic_cif(tmp_path):
    """Synthetic coordinates, not an actual protein or asserted PDB entry."""
    structure = Structure.Structure("synthetic")
    model = Model.Model(0)
    structure.add(model)
    for chain_id in ("A", "B"):
        chain = Chain.Chain(chain_id)
        model.add(chain)
        residue = Residue.Residue((" ", 1, " "), "ALA", " ")
        chain.add(residue)
        residue.add(Atom.Atom("CA", np.array([1.0, 2.0, 3.0]), 1.0, 1.0,
                              " ", "CA", 1, element="C"))
    destination = tmp_path / "synthetic.cif"
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(destination))
    return destination.read_bytes()


def test_unmapped_or_missing_chain_never_submits(tmp_path):
    assert prepare_query({"status": "not_found"}, tmp_path / "query.cif")["status"] == "skipped"
    assert not (tmp_path / "query.cif").exists()


def test_extracts_only_pdbe_mapped_chain(tmp_path):
    cif = synthetic_cif(tmp_path)
    seen = []

    def respond(request):
        seen.append(request.url)
        return httpx.Response(200, content=cif)

    client = httpx.Client(transport=httpx.MockTransport(respond))
    pdbe = {"status": "verified", "structures": [{
        "pdb_id": "1abc", "chain_id": "A", "coverage_fraction": 0.75,
    }], "_provenance": {"source": "synthetic test mapping"}}
    result = prepare_query(pdbe, tmp_path / "query.cif", client=client)
    assert result["status"] == "ready"
    assert result["chain_id"] == "A"
    assert result["source_url"] == "https://www.ebi.ac.uk/pdbe/entry-files/download/1abc.cif"
    assert len(seen) == 1
    parsed = MMCIFParser(QUIET=True).get_structure("query", io.StringIO(
        (tmp_path / "query.cif").read_text()
    ))
    assert [chain.id for chain in next(parsed.get_models())] == ["A"]
    pdbe["structures"][0]["chain_id"] = "Z"
    assert prepare_query(pdbe, tmp_path / "missing.cif", client=client)["status"] == "skipped"


def test_sequence_only_question_cannot_offer_foldseek_structure_search():
    plan = build_plan("Foldseek this sequence ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY")
    assert FOLDSEEK_SEARCH not in plan.tool_names


def test_foldseek_rejects_fasta_before_http():
    with httpx.Client(
        transport=httpx.MockTransport(lambda _req: pytest.fail("HTTP called"))
    ) as http:
        foldseek = FoldseekClient(Settings(), client=http)
        with pytest.raises(RetrievalError, match="FASTA is not a structure"):
            foldseek.search(">sequence\nACDEFGHIKLMNPQRSTVWY", database="pdb100")


def test_foldseek_submits_cif_and_records_ticket(tmp_path):
    cif = synthetic_cif(tmp_path).decode()
    seen = []

    def respond(request):
        seen.append(request)
        if request.url.path == "/api/ticket" and request.method == "POST":
            return httpx.Response(200, json={"id": "testticket", "status": "COMPLETE"})
        if request.url.path == "/api/result/testticket/0":
            return httpx.Response(200, json={
                "queries": [{"header": "synthetic"}], "mode": "3diaa",
                "results": [{"db": "pdb100", "alignments": [[
                    {"target": "synthetic-hit", "alnLength": 20, "eval": 0.1,
                     "score": 30, "qAln": "unbounded alignment not for report"}
                ]]}],
            })
        raise AssertionError(request.url)

    with httpx.Client(transport=httpx.MockTransport(respond), base_url="https://fold.test") as http:
        result = FoldseekClient(Settings(), client=http).search(cif)
    assert result["status"] == "complete"
    assert result["ticket"] == "testticket"
    assert result["result"]["hits"][0]["target"] == "synthetic-hit"
    assert "qAln" not in result["result"]["hits"][0]
    assert b"query.cif" in seen[0].content
    assert b"name=\"mode\"" in seen[0].content
    assert b"3diaa" in seen[0].content
    assert b"data_" in seen[0].content


def test_report_skips_foldseek_when_no_verified_structure(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "protein_agent.report.find_sequence_annotations",
        lambda _seq: {"status": "not_found", "reason": "No exact record."},
    )
    source = tmp_path / "sequence.fa"
    source.write_text(">synthetic sequence\n" + "ACDEFGHIKLMNPQRSTVWY" * 2 + "\n")
    path = generate_report(source, tmp_path / "results", run_foldseek=True)
    assert path.is_file()
    import json

    data = json.loads((path.parent / "report.json").read_text())
    assert data["foldseek"]["status"] == "skipped"
    assert "structure" in data["foldseek"]["reason"].lower()
    assert "Not completed" in path.read_text()


def test_report_search_uses_prepared_structure_and_cleans_query(monkeypatch, tmp_path):
    prepared = {"status": "ready", "pdb_id": "1ABC", "chain_id": "A"}
    seen = {}

    def fake_prepare(pdbe, destination):
        seen["mapping"] = pdbe
        destination.write_text("data_synthetic\n")
        return prepared

    class FakeFoldseek:
        def __init__(self, settings):
            seen["timeout"] = settings.request_timeout

        def search(self, query, **kwargs):
            seen["query"] = query
            return {"status": "complete", "database": kwargs["database"],
                    "result": {"total_hits_reported": 0, "hits": []}}

        def close(self):
            pass

    monkeypatch.setattr(report_module, "prepare_query", fake_prepare)
    monkeypatch.setattr(report_module, "FoldseekClient", FakeFoldseek)
    mapping = {"status": "verified", "structures": [{"pdb_id": "1ABC", "chain_id": "A"}]}
    result = report_module.run_foldseek_search(mapping, tmp_path)
    assert result["status"] == "complete"
    assert result["query_structure"] == prepared
    assert seen["mapping"] is mapping and seen["query"].startswith("data_")
    assert seen["timeout"] <= 15
    assert not (tmp_path / ".foldseek-query.cif").exists()


def test_foldseek_report_is_user_focused_and_escapes_hits():
    content = render_foldseek_section({"foldseek": {
        "status": "complete", "database": "pdb100",
        "query_structure": {"pdb_id": "1ABC", "chain_id": "A"},
        "result": {"total_hits_reported": 1, "hits": [
            {"target": "<script>synthetic</script>", "alnLength": 20, "eval": 0.1}
        ]},
    }})
    assert "1ABC chain A" in content
    assert "<script>" not in content
    assert "&lt;script&gt;" in content
    assert "<pre>" not in content
    skipped = render_foldseek_section({"foldseek": {
        "status": "skipped", "reason": "No mapped structure was available.",
    }})
    assert "Not completed" in skipped
    assert "similar structures are absent" in skipped
