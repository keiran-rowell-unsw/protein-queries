"""Passenger WSGI contract without an OOD server or remote SSH."""

import shlex
from io import BytesIO
from types import SimpleNamespace
from urllib.parse import urlencode

import ood_app
import pytest

from protein_agent.report import generate_report


def request(method, path, *, body=b"", cookie="", prefix="/pun/dev/protein-queries"):
    captured = {}

    def start(status, headers):
        captured["status"] = status
        captured["headers"] = dict(headers)

    environ = {
        "REQUEST_METHOD": method,
        "PATH_INFO": path,
        "SCRIPT_NAME": prefix,
        "CONTENT_LENGTH": str(len(body)),
        "CONTENT_TYPE": "application/x-www-form-urlencoded",
        "HTTP_COOKIE": cookie,
        "wsgi.input": BytesIO(body),
    }
    captured["body"] = b"".join(ood_app.application(environ, start))
    return captured


def test_manifest_and_wsgi_entrypoint():
    from pathlib import Path

    assert "role: passenger" in Path("manifest.yml").read_text()
    from passenger_wsgi import application

    assert application is ood_app.application


def test_session_prefix_csrf_and_report(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    monkeypatch.setattr(
        "protein_agent.report.find_sequence_annotations",
        lambda _seq: {"status": "not_found", "reason": "No exact match under these criteria."},
    )

    def compute(text, folder, accession=""):
        (folder / "input.fasta").write_text(text)
        generate_report(folder / "input.fasta", folder)
        (folder / "disorder.png").write_bytes(b"test image")
        (folder / "input.fasta").unlink()

    monkeypatch.setattr(ood_app, "run_report", compute)
    form = request("GET", "/")
    assert form["status"] == "200 OK"
    assert b"/pun/dev/protein-queries/run" in form["body"]
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    params = {"csrf": token, "sequence": ">example\nACDEFGHIKLMNPQRSTVWY"}
    submit = request("POST", "/run", body=urlencode(params).encode(), cookie=cookie)
    assert submit["status"] == "303 See Other"
    url = submit["headers"]["Location"]
    assert url.startswith("/pun/dev/protein-queries/results/")
    report = request("GET", url.removeprefix("/pun/dev/protein-queries"))
    assert report["status"] == "200 OK"
    assert b"Sequence analysis report" in report["body"]
    image = request(
        "GET", url.removeprefix("/pun/dev/protein-queries").replace("report.html", "disorder.png")
    )
    assert image["status"] == "200 OK"
    assert image["headers"]["Content-Type"] == "image/png"
    assert (
        request("GET", "/results/" + "a" * 32 + "/../../input.fasta")["status"] == "404 Not Found"
    )


def test_accession_only_submission_is_accepted(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    captured = {}

    def fake_report(text, folder, accession=""):
        captured.update(text=text, accession=accession)
        (folder / "report.html").write_text("generated")

    monkeypatch.setattr(ood_app, "run_report", fake_report)
    form = request("GET", "/")
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    params = {"csrf": token, "accession": "Q12345"}
    result = request("POST", "/run", body=urlencode(params).encode(), cookie=cookie)
    assert result["status"] == "303 See Other"
    assert captured == {"text": "", "accession": "Q12345"}


def test_csrf_and_invalid_input(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    assert request("POST", "/run", body=b"sequence=ACDE")["status"] == "403 Forbidden"
    form = request("GET", "/")
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    params = {"csrf": token, "sequence": ">a\nACD\n>b\nEFG"}
    result = request("POST", "/run", body=urlencode(params).encode(), cookie=cookie)
    assert result["status"] == "400 Bad Request"
    assert list(tmp_path.iterdir()) == []


def test_remote_run_uses_ssh_without_shell_or_pbs(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    seen = []

    def fake_run(args, **kwargs):
        seen.append((args, kwargs))
        (folder / "report.html").write_text("generated")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(ood_app.subprocess, "run", fake_run)
    ood_app.run_report("ACDE", folder)
    args, kwargs = seen[0]
    assert args[0] == "ssh" and args[-2] == ood_app.HOST
    assert shlex.split(args[-1])[-4:] == [
        "--input",
        str(folder / "input.fasta"),
        "--output",
        str(folder),
    ]
    assert kwargs["timeout"] == 90 and kwargs["check"] is False
    assert not (folder / "input.fasta").exists()
    ood_app.run_report("ACDE", folder, "Q12345")
    assert shlex.split(seen[-1][0][-1])[-2:] == ["--accession", "Q12345"]


@pytest.mark.parametrize("host", ["bad host", "-bad", "host;echo"])
def test_ssh_host_config_rejects_shell_input(monkeypatch, tmp_path, host):
    monkeypatch.setattr(ood_app, "HOST", host)
    with pytest.raises(RuntimeError):
        ood_app.run_report("ACDE", tmp_path)
