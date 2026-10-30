"""Offline report and local web submission regressions."""


import json
from http.server import ThreadingHTTPServer
from threading import Thread
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pytest

from protein_agent.properties import SequenceError
from protein_agent.report import generate_report, parse_sequence
from protein_agent.web import make_handler


@pytest.fixture(autouse=True)
def no_live_uniparc(monkeypatch):
    monkeypatch.setattr(
        "protein_agent.report.find_sequence_annotations",
        lambda _seq: {"status": "not_found", "reason": "No exact match under these criteria."},
    )


def _accession_record(sequence):
    return {
        "status": "verified",
        "accession": "Q12345",
        "organism": "Example organism",
        "organism_details": {"taxonId": 1234},
        "entry_type": "UniProtKB reviewed (Swiss-Prot)",
        "entry_version": 1,
        "retrieved_at": "2026-01-01T00:00:00+00:00",
        "canonical_sequence_match": None,
        "identity_basis": "supplied_accession",
        "canonical_sequence": sequence,
        "protein_description": None,
        "genes": [],
        "comments": [],
        "keywords": [],
        "cross_references": [],
        "protein_interactions": [],
        "entry_audit": {},
        "features": [],
        "absorption": None,
        "source_url": "https://www.uniprot.org/uniprotkb/Q12345/entry",
    }


def test_report_is_complete_and_escapes_header(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">hello <script>alert(1)</script>\nACDEFGHIKLMNPQRSTVWY\n")
    path = generate_report(source, tmp_path / "results")
    page = path.read_text()
    assert "<script>" not in page
    assert "&lt;script&gt;" in page
    assert "<h2>Methods</h2>" not in page
    assert "Charge at different pH values" in page
    assert "Download full results (JSON)" in page
    assert "neither establishes an expression outcome" in page
    assert "Predicted disorder" in page
    assert "Other sequence estimates" in page
    assert "Amino-acid composition" in page
    assert (path.parent / "hydropathy.png").stat().st_size > 0
    data = json.loads((path.parent / "report.json").read_text())
    assert data["length_residues"] == 20
    assert "charge_vs_ph" in data
    assert "method" in data  # provenance remains in the downloadable record
    assert "chromophore_candidates" in data
    assert "inclusion_body_propensity" in data
    assert data["disorder_prediction"]["status"] == "skipped"
    assert "shorter than 25" in data["disorder_prediction"]["reason"]
    assert not (path.parent / "disorder.png").exists()


def test_accession_only_uses_canonical_sequence(monkeypatch, tmp_path):
    canonical = "ACDEFGHIKLMNPQRSTVWY"
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: _accession_record(canonical),
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    path = generate_report(None, tmp_path / "accession-only", accession="q12345")
    data = json.loads((path.parent / "report.json").read_text())
    assert data["length_residues"] == len(canonical)
    assert data["input_source"].startswith("canonical sequence retrieved")
    assert "Supplied accession resolved to" in path.read_text()


def test_accession_precedes_different_submitted_sequence(monkeypatch, tmp_path):
    canonical = "ACDEFGHIKLMNPQRSTVWY"
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: _accession_record(canonical),
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    source = tmp_path / "different.fa"
    source.write_text("A" * 25)
    path = generate_report(source, tmp_path / "precedence", accession="Q12345")
    data = json.loads((path.parent / "report.json").read_text())
    assert data["length_residues"] == len(canonical)
    assert data["uniprot"]["submitted_sequence_match"] is False
    assert "accession precedence was applied" in path.read_text()


@pytest.mark.parametrize("sequence", [">a\nACD\n>b\nEFG", "ACD-X", "ACD123", "A" * 10001])
def test_invalid_input_is_rejected_not_modified(sequence):
    with pytest.raises(SequenceError):
        parse_sequence(sequence)


def test_local_web_accepts_accession_without_sequence(monkeypatch, tmp_path):
    canonical = "ACDEFGHIKLMNPQRSTVWY"
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: _accession_record(canonical),
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            payload = urlencode({"accession": "Q12345"}).encode()
            with urlopen(Request(base + "/run", data=payload)) as response:
                assert "Supplied accession resolved to" in response.read().decode()
        finally:
            server.shutdown()
            thread.join()


def test_local_web_submission(tmp_path):
    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            with urlopen(
                Request(
                    base + "/run", data=urlencode({"sequence": "ACDEFGHIKLMNPQRSTVWY"}).encode()
                )
            ) as response:
                assert "Protein analysis report" in response.read().decode()
                assert response.url.startswith(base + "/results/")
            with pytest.raises(HTTPError) as error:
                urlopen(base + "/results/" + "a" * 32 + "/../../input.fasta")
            assert error.value.code == 404
        finally:
            server.shutdown()
            thread.join()


def test_report_json_carries_the_new_sections(monkeypatch, tmp_path):
    """The three added sections must reach report.json, not just the rendered page."""
    canonical = "ACDEFGHIKLMNPQRSTVWY" * 3
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: _accession_record(canonical),
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    path = generate_report(None, tmp_path / "sections", accession="Q12345")
    data = json.loads((path.parent / "report.json").read_text())
    assert data["predicted_motifs"]["status"] == "computed"
    assert data["position_map"]["status"] in {"rendered", "empty", "skipped"}
    assert data["structure_display"]["status"] == "not_found"
    # A withheld structure must not be phrased as an absence of data.
    assert "not deposited" not in path.read_text().lower() or data["pdbe"]["status"] == "not_found"


def test_local_web_serves_every_image_the_report_references(monkeypatch, tmp_path):
    """An <img> tag whose file the server refuses to serve is a broken report.

    This guards the web.py allowlist: each new figure has to be both written to disk and
    reachable by URL, and the two lists (allowlist, MIME map) have to agree.
    """
    canonical = "ACDEFGHIKLMNPQRSTVWY" * 3
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: _accession_record(canonical),
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    token = "a" * 32
    output_root = tmp_path / "served"
    output = output_root / token
    report = generate_report(None, output, accession="Q12345")
    page = report.read_text()
    # The fixture entry annotates no features and its sequence holds no sequon, so the map
    # is legitimately absent. Asserting it here anyway would test the assertion, not the
    # report; the features-present case below is what proves the figure reaches the page.
    data = json.loads((output / "report.json").read_text())
    assert data["position_map"]["status"] == "empty"
    assert data["structure_display"]["status"] == "not_found"

    referenced = {
        name
        for name in (
            "hydropathy.png", "disorder.png", "positions.png",
            "structure.png", "ligands.png",
        )
        if 'src="' + name + '"' in page or 'href="' + name + '"' in page
    }
    assert referenced == {"hydropathy.png", "disorder.png"}, referenced

    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(output_root)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            for name in sorted(referenced):
                with urlopen(f"{base}/results/{token}/{name}") as response:
                    assert response.status == 200, name
                    assert response.headers["Content-Type"] == "image/png", name
                    assert len(response.read()) > 100, f"{name} served empty"
        finally:
            server.shutdown()
            thread.join()


def test_local_web_serves_the_position_map_when_features_exist(monkeypatch, tmp_path):
    """With annotated features the map renders, and its file must be servable."""
    canonical = "ACDEFGHIKLMNPQRSTVWY" * 3
    record = _accession_record(canonical)
    record["features"] = [
        {"type": "Glycosylation", "start": 4, "end": 4, "description": "N-linked (GlcNAc)"},
    ]
    monkeypatch.setattr(
        "protein_agent.report.fetch_annotations",
        lambda sequence, accession, accession_precedence=False: record,
    )
    monkeypatch.setattr(
        "protein_agent.report.fetch_pdbe_structures",
        lambda accession: {"status": "not_found", "reason": "No mapping."},
    )
    token = "b" * 32
    output_root = tmp_path / "mapped"
    output = output_root / token
    page = generate_report(None, output, accession="Q12345").read_text()
    assert 'src="positions.png"' in page
    assert (output / "positions.png").is_file()

    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(output_root)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            with urlopen(
                f"http://127.0.0.1:{server.server_port}/results/{token}/positions.png"
            ) as response:
                assert response.status == 200
                assert response.headers["Content-Type"] == "image/png"
                assert len(response.read()) > 100
        finally:
            server.shutdown()
            thread.join()


def test_local_web_still_refuses_files_that_are_not_report_assets(tmp_path):
    """Widening the allowlist must not become a file read."""
    output = tmp_path / "strict"
    (output / ("a" * 32)).mkdir(parents=True)
    secret = output / ("a" * 32) / "secret.txt"
    secret.write_text("do not serve me", encoding="utf-8")
    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(output)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            with pytest.raises(HTTPError) as error:
                urlopen(f"{base}/results/{'a' * 32}/secret.txt")
            assert error.value.code == 404
            for name in ("positions.png", "structure.png", "ligands.png"):
                path = output / ("a" * 32) / name
                path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 200)
                with urlopen(f"{base}/results/{'a' * 32}/{name}") as response:
                    assert response.status == 200
                    assert response.headers["Content-Type"] == "image/png"
        finally:
            server.shutdown()
            thread.join()
