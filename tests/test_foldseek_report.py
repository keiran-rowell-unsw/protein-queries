"""Foldseek auto-reporting only searches a retrieved, PDBe-mapped structure."""

import io

import httpx
import numpy as np
import pytest
from Bio.PDB import MMCIFIO, Atom, Chain, MMCIFParser, Model, Residue, Structure

from protein_agent import report as report_module
from protein_agent.config import ConfigError, Settings
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
        residue.add(Atom.Atom("CA", np.array([1.0, 2.0, 3.0]), 1.0, 1.0, " ", "CA", 1, element="C"))
        water = Residue.Residue(("W", 2, " "), "HOH", " ")
        chain.add(water)
        water.add(Atom.Atom("O", np.array([2.0, 3.0, 4.0]), 1.0, 1.0, " ", "O", 2, element="O"))
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
    pdbe = {
        "status": "verified",
        "structures": [
            {
                "pdb_id": "1abc",
                "chain_id": "A",
                "coverage_fraction": 0.75,
            }
        ],
        "_provenance": {"source": "synthetic test mapping"},
    }
    result = prepare_query(pdbe, tmp_path / "query.cif", client=client)
    assert result["status"] == "ready"
    assert result["chain_id"] == "A"
    assert result["source_url"] == "https://www.ebi.ac.uk/pdbe/entry-files/download/1abc.cif"
    assert len(seen) == 1
    parsed = MMCIFParser(QUIET=True).get_structure(
        "query", io.StringIO((tmp_path / "query.cif").read_text())
    )
    selected = next(parsed.get_models())
    assert [chain.id for chain in selected] == ["A"]
    assert [residue.resname for residue in selected["A"]] == ["ALA"]
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
            return httpx.Response(
                200,
                json={
                    "queries": [{"header": "synthetic"}],
                    "mode": "3diaa",
                    "results": [
                        {
                            "db": "pdb100",
                            "alignments": [
                                [
                                    {
                                        "target": "synthetic-hit",
                                        "alnLength": 20,
                                        "eval": 0.1,
                                        "score": 30,
                                        "qAln": "unbounded alignment not for report",
                                    }
                                ]
                            ],
                        }
                    ],
                },
            )
        raise AssertionError(request.url)

    with httpx.Client(transport=httpx.MockTransport(respond), base_url="https://fold.test") as http:
        result = FoldseekClient(Settings(), client=http).search(cif, database="pdb100")
    assert result["status"] == "complete"
    assert result["ticket"] == "testticket"
    assert result["result"]["hits"][0]["target"] == "synthetic-hit"
    assert "qAln" not in result["result"]["hits"][0]
    assert b"query.cif" in seen[0].content
    assert b'name="mode"' in seen[0].content
    assert b"3diaa" in seen[0].content
    assert b"data_" in seen[0].content
    assert b"format=json" in seen[-1].url.query


def test_foldseek_discovers_and_merges_all_standard_databases(tmp_path):
    cif = synthetic_cif(tmp_path).decode()
    seen = []

    def respond(request):
        seen.append(request)
        if request.url.path == "/api/databases/all":
            return httpx.Response(
                200,
                json={
                    "databases": [
                        {
                            "path": "pdb100",
                            "status": "COMPLETE",
                            "default": True,
                            "interface": False,
                            "motif": False,
                            "rna": False,
                        },
                        {
                            "path": "afdb50",
                            "status": "COMPLETE",
                            "default": True,
                            "interface": False,
                            "motif": False,
                            "rna": False,
                        },
                        {
                            "path": "motif-db",
                            "status": "COMPLETE",
                            "default": True,
                            "interface": False,
                            "motif": True,
                            "rna": False,
                        },
                        {
                            "path": "offline",
                            "status": "PENDING",
                            "default": True,
                            "interface": False,
                            "motif": False,
                            "rna": False,
                        },
                    ]
                },
            )
        if request.url.path == "/api/ticket":
            return httpx.Response(200, json={"id": "all-ticket", "status": "COMPLETE"})
        return httpx.Response(
            200,
            json={
                "queries": [{"header": "synthetic"}],
                "mode": "3diaa",
                "results": [
                    {
                        "db": "pdb100",
                        "alignments": [[{"target": "experimental", "eval": 1e-5, "score": 20}]],
                    },
                    {
                        "db": "afdb50",
                        "alignments": [[{"target": "predicted", "eval": 1e-20, "score": 50}]],
                    },
                ],
            },
        )

    with httpx.Client(transport=httpx.MockTransport(respond), base_url="https://fold.test") as http:
        result = FoldseekClient(Settings(), client=http).search(cif)
    submission = next(request for request in seen if request.url.path == "/api/ticket")
    assert submission.content.count(b'name="database[]"') == 2
    assert b"pdb100" in submission.content and b"afdb50" in submission.content
    assert b"motif-db" not in submission.content and b"offline" not in submission.content
    assert result["databases"] == ["pdb100", "afdb50"]
    assert result["result"]["total_hits_reported"] == 2
    assert result["result"]["hits"][0]["target"] == "predicted"
    assert result["result"]["hits"][0]["database"] == "afdb50"


def test_foldseek_reports_live_style_rate_limit():
    def respond(_request):
        return httpx.Response(
            429,
            json={
                "status": "RATELIMIT",
                "reason": "The Foldseek server is a shared resource.",
            },
        )

    with httpx.Client(transport=httpx.MockTransport(respond), base_url="https://fold.test") as http:
        result = FoldseekClient(Settings(), client=http).search(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000",
            database="pdb100",
        )
    assert result["status"] == "ratelimit"
    assert "shared resource" in result["reason"]


def test_foldseek_reports_plain_text_rate_limit():
    with httpx.Client(
        transport=httpx.MockTransport(
            lambda _request: httpx.Response(429, text="Too many requests")
        ),
        base_url="https://fold.test",
    ) as http:
        result = FoldseekClient(Settings(), client=http).search(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000",
            database="pdb100",
        )
    assert result["status"] == "ratelimit"
    assert "rate-limited" in result["reason"]


def test_foldseek_poll_timeout_config_is_validated(monkeypatch):
    monkeypatch.setenv("FOLDSEEK_POLL_TIMEOUT", "not-a-number")
    with pytest.raises(ConfigError, match="numeric environment setting"):
        Settings.from_env(require_key=False)
    monkeypatch.setenv("FOLDSEEK_POLL_TIMEOUT", "0")
    with pytest.raises(ConfigError, match="must be > 0"):
        Settings.from_env(require_key=False)


def test_foldseek_queue_deadline_is_not_http_timeout(monkeypatch):
    polls = 0

    def respond(request):
        nonlocal polls
        if request.method == "POST":
            return httpx.Response(200, json={"id": "queued", "status": "PENDING"})
        if request.url.path == "/api/ticket/queued":
            polls += 1
            return httpx.Response(
                200,
                json={
                    "id": "queued",
                    "status": "COMPLETE" if polls == 2 else "RUNNING",
                },
            )
        return httpx.Response(
            200,
            json={
                "queries": [{"header": "synthetic"}],
                "mode": "3diaa",
                "results": [{"db": "pdb100", "alignments": []}],
            },
        )

    monkeypatch.setattr("protein_agent.foldseek.time.sleep", lambda _seconds: None)
    settings = Settings(request_timeout=0.01, foldseek_poll_timeout=30.0)
    with httpx.Client(transport=httpx.MockTransport(respond), base_url="https://fold.test") as http:
        result = FoldseekClient(settings, client=http).search(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000",
            database="pdb100",
        )
    assert result["status"] == "complete"
    assert polls == 2


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
            seen["kwargs"] = kwargs
            return {
                "status": "complete",
                "databases": ["db-one", "db-two"],
                "result": {"total_hits_reported": 0, "hits": []},
            }

        def close(self):
            pass

    monkeypatch.setattr(report_module, "prepare_query", fake_prepare)
    monkeypatch.setattr(report_module, "FoldseekClient", FakeFoldseek)
    mapping = {"status": "verified", "structures": [{"pdb_id": "1ABC", "chain_id": "A"}]}
    result = report_module.run_foldseek_search(mapping, tmp_path)
    assert result["status"] == "complete"
    assert result["query_structure"] == prepared
    assert seen["mapping"] is mapping and seen["query"].startswith("data_")
    assert seen["kwargs"] == {}
    assert seen["timeout"] <= 30
    assert not (tmp_path / ".foldseek-query.cif").exists()


def test_foldseek_report_is_user_focused_and_escapes_hits():
    content = render_foldseek_section(
        {
            "foldseek": {
                "status": "complete",
                "databases": ["pdb100", "afdb50"],
                "query_structure": {"pdb_id": "1ABC", "chain_id": "A"},
                "result": {
                    "total_hits_reported": 1,
                    "hits": [
                        {
                            "target": "<script>synthetic</script>",
                            "database": "afdb50",
                            "alnLength": 20,
                            "eval": 0.1,
                        }
                    ],
                },
            }
        }
    )
    assert "1ABC chain A" in content
    assert "<script>" not in content
    assert "&lt;script&gt;" in content
    assert "Databases searched (2): pdb100, afdb50" in content
    assert "<td>afdb50</td>" in content
    assert "<pre>" not in content
    skipped = render_foldseek_section(
        {
            "foldseek": {
                "status": "skipped",
                "reason": "No mapped structure was available.",
            }
        }
    )
    assert "Not completed" in skipped
    assert "similar structures are absent" in skipped
