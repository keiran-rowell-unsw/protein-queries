"""Offline report and local web submission regressions."""

# ruff: noqa: E501

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
    assert "No\nexperimental structure or expression outcome was determined." in page
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


@pytest.mark.parametrize("sequence", [">a\nACD\n>b\nEFG", "ACD-X", "ACD123", "A" * 10001])
def test_invalid_input_is_rejected_not_modified(sequence):
    with pytest.raises(SequenceError):
        parse_sequence(sequence)


def test_local_web_submission(tmp_path):
    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            with urlopen(Request(base + "/run", data=urlencode({"sequence": "ACDEFGHIKLMNPQRSTVWY"}).encode())) as response:
                assert "Sequence analysis report" in response.read().decode()
                assert response.url.startswith(base + "/results/")
            with pytest.raises(HTTPError) as error:
                urlopen(base + "/results/" + "a" * 32 + "/../../input.fasta")
            assert error.value.code == 404
        finally:
            server.shutdown()
            thread.join()
