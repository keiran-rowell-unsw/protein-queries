"""Passenger WSGI contract: validate and queue; never calculate on the OOD host."""

import json
import re
from http.server import ThreadingHTTPServer
from io import BytesIO
from pathlib import Path
from threading import Thread
from types import SimpleNamespace
from urllib.parse import urlencode
from urllib.request import urlopen

import ood_app
import pytest


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


def form_cookie():
    response = request("GET", "/")
    cookie = response["headers"]["Set-Cookie"]
    return response, cookie, cookie.split(";", 1)[0].split("=", 1)[1]


def test_manifest_wsgi_and_form_describe_queued_reports():
    assert "role: passenger" in Path("manifest.yml").read_text()
    from passenger_wsgi import application

    assert application is ood_app.application
    form, _cookie, _csrf = form_cookie()
    assert form["status"] == "200 OK"
    assert b"Queues a job to fetch results" in form["body"]
    assert b"Queues a full report job" in form["body"]
    assert b"maximum 5,000 residues" in form["body"]
    assert b"Queue report" in form["body"]


@pytest.mark.parametrize(
    ("sequence", "accession"),
    [
        (">submitted\nACDEFGHIKLMNPQRSTVWY", ""),
        ("", "Q12345"),
        ("ACDEFGHIKLMNPQRSTVWY", "Q12345"),
    ],
)
def test_every_submission_queues_without_local_report(
    monkeypatch, tmp_path, sequence, accession
):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    monkeypatch.setattr(ood_app, "PBS_STATUS_COMMAND", "/definitely/missing/qstat")
    observed = {}

    def fake_submit(folder):
        observed["folder"] = folder
        observed["sequence"] = (folder / "input.fasta").read_text()
        observed["request"] = json.loads((folder / "request.json").read_text())
        assert not (folder / "report.html").exists()
        (folder / "job.json").write_text('{"status":"queued","job_id":"123.server"}')
        return {"status": "queued", "job_id": "123.server"}

    monkeypatch.setattr(ood_app, "submit_report_job", fake_submit)
    _form, cookie, csrf = form_cookie()
    response = request(
        "POST", "/run",
        body=urlencode({"csrf": csrf, "sequence": sequence, "accession": accession}).encode(),
        cookie=cookie,
    )
    assert response["status"] == "303 See Other"
    assert observed["sequence"] == sequence
    assert observed["request"] == {"accession": accession}
    pending = request(
        "GET", response["headers"]["Location"].removeprefix("/pun/dev/protein-queries")
    )
    assert pending["status"] == "200 OK"
    assert b"Report queued" in pending["body"]
    assert b"123.server" in pending["body"]
    assert b"OOD Active Jobs" in pending["body"]
    assert b"/pun/sys/dashboard/activejobs" in pending["body"]
    assert b"OOD Active Jobs" in pending["body"]
    assert b"<style>" in pending["body"]


def test_reconcile_pbs_status_marks_disappeared_job_unavailable(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    job = {"status": "running", "job_id": "123.server"}
    monkeypatch.setattr(ood_app, "PBS_STATUS_COMMAND", "/usr/bin/qstat")
    monkeypatch.setattr(ood_app.Path, "is_file", lambda self: str(self) == "/usr/bin/qstat")
    monkeypatch.setattr(ood_app.os, "access", lambda *_args: True)
    monkeypatch.setattr(
        ood_app.subprocess, "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=153, stdout=""),
    )
    result = ood_app.reconcile_pbs_status(folder, job)
    assert result["status"] == "unavailable"
    assert "no longer reports" in result["reason"]


def test_reconcile_pbs_status_maps_scheduler_state(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    monkeypatch.setattr(ood_app, "PBS_STATUS_COMMAND", "/usr/bin/qstat")
    monkeypatch.setattr(ood_app.Path, "is_file", lambda _self: True)
    monkeypatch.setattr(ood_app.os, "access", lambda *_args: True)
    monkeypatch.setattr(
        ood_app.subprocess, "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout="job_state = R\n"),
    )
    result = ood_app.reconcile_pbs_status(folder, {"status": "queued", "job_id": "123.server"})
    assert result["status"] == "running"


def test_submission_retains_optional_email(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    monkeypatch.setattr(ood_app, "submit_report_job", lambda _folder: {"status": "queued"})
    _form, cookie, csrf = form_cookie()
    response = request(
        "POST", "/run",
        body=urlencode({
            "csrf": csrf, "sequence": "ACDEFGHIKLMNPQRSTVWY",
            "email": "user@unsw.edu.au",
        }).encode(),
        cookie=cookie,
    )
    assert response["status"] == "303 See Other"
    folder = next(tmp_path.iterdir())
    assert json.loads((folder / "request.json").read_text()) == {
        "accession": "", "email": "user@unsw.edu.au"
    }


def test_submission_failure_removes_partial_request(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    monkeypatch.setattr(
        ood_app, "submit_report_job",
        lambda _folder: (_ for _ in ()).throw(RuntimeError("qsub unavailable")),
    )
    _form, cookie, csrf = form_cookie()
    response = request(
        "POST", "/run",
        body=urlencode({"csrf": csrf, "sequence": "ACDEFGHIKLMNPQRSTVWY"}).encode(),
        cookie=cookie,
    )
    assert response["status"] == "503 Service Unavailable"
    assert b"could not be queued" in response["body"]
    assert list(tmp_path.iterdir()) == []


def test_csrf_invalid_fasta_and_residue_limit_are_rejected(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    assert request("POST", "/run", body=b"sequence=ACDE")["status"] == "403 Forbidden"
    _form, cookie, csrf = form_cookie()
    for sequence in (">a\nACD\n>b\nEFG", "A" * 5001):
        response = request(
            "POST", "/run",
            body=urlencode({"csrf": csrf, "sequence": sequence}).encode(),
            cookie=cookie,
        )
        assert response["status"] == "400 Bad Request"
    malformed_accession = request(
        "POST", "/run",
        body=urlencode({"csrf": csrf, "accession": "A123456"}).encode(),
        cookie=cookie,
    )
    assert malformed_accession["status"] == "400 Bad Request"
    assert list(tmp_path.iterdir()) == []
    ood_app.validate_sequence(">synthetic\n" + "A" * 2500 + "\n" + "A" * 2500)


def test_report_email_is_optional_and_added_to_pbs_job(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    (folder / "input.fasta").write_text("ACDEFGHIKLMNPQRSTVWY")
    monkeypatch.setattr(ood_app, "EXTENDED_PYTHON", "/compute/project/.venv/bin/python")
    monkeypatch.setattr(ood_app, "REPORT_EMAIL", "requestor@unsw.edu.au")
    monkeypatch.setattr(
        ood_app.subprocess, "run",
        lambda *_a, **_kw: SimpleNamespace(returncode=0, stdout="123.server\n"),
    )
    ood_app.submit_report_job(folder)
    script = (folder / ".report.pbs").read_text()
    assert "--recipient requestor@unsw.edu.au" in script
    assert "-m abe" in script


def test_submit_report_job_writes_one_full_pbs_job(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    (folder / "input.fasta").write_text(">x\nACDEFGHIKLMNPQRSTVWY")
    (folder / "request.json").write_text(json.dumps({"accession": "Q12345"}))
    seen = {}

    def fake_run(args, **kwargs):
        seen.update(args=args, kwargs=kwargs)
        return SimpleNamespace(returncode=0, stdout="123.server\n")

    monkeypatch.setattr(ood_app, "EXTENDED_PYTHON", "/compute/project/.venv/bin/python")
    monkeypatch.setattr(ood_app.subprocess, "run", fake_run)
    result = ood_app.submit_report_job(folder)
    assert result == {"status": "queued", "job_id": "123.server"}
    assert seen["args"] == [
        ood_app.PBS_SUBMIT_COMMAND,
        "-o", str(folder / "job.log"),
        "-j", "oe",
        str(folder / ".report.pbs"),
    ]
    assert seen["kwargs"]["timeout"] == 15 and seen["kwargs"]["check"] is False
    assert (folder / "job-id.txt").read_text().strip() == "123.server"
    script = (folder / ".report.pbs").read_text()
    assert "#PBS -N protein-report" in script
    assert "--foldseek" in script and "--accession Q12345" in script
    assert "--lightweight" not in script and "ssh" not in script
    assert "#PBS -l ncpus=2" in script and "#PBS -l mem=4gb" in script
    assert "#PBS -m abe" in script


def test_default_compute_environment_and_missing_environment(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    (folder / "input.fasta").write_text("ACDEFGHIKLMNPQRSTVWY")
    python = tmp_path / ".compute-venv" / "bin" / "python"
    monkeypatch.setattr(ood_app, "EXTENDED_PYTHON", "")
    monkeypatch.setattr(ood_app, "DEFAULT_COMPUTE_PYTHON", str(python))
    monkeypatch.setattr(
        ood_app.subprocess, "run", lambda *_a, **_kw: pytest.fail("qsub must not run")
    )
    with pytest.raises(RuntimeError, match="bin/setup-compute-env"):
        ood_app.submit_report_job(folder)
    assert not (folder / ".report.pbs").exists()


def test_compute_setup_is_locked_and_separate_from_passenger():
    script = Path("bin/setup-compute-env").read_text()
    assert "UV_PROJECT_ENVIRONMENT=" in script
    assert "sync --locked --no-dev" in script
    assert "import metapredict" in script and "import torch" in script
    assert ".compute-venv/" in Path(".gitignore").read_text()


def test_pending_failure_escapes_reason_and_completed_report_is_served(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    token = "a" * 32
    folder = tmp_path / token
    folder.mkdir()
    (folder / "job.json").write_text(
        json.dumps({"status": "failed", "reason": "<script>failed</script>"})
    )
    failed = request("GET", f"/results/{token}/report.html")
    assert b"Report unavailable" in failed["body"]
    assert b"<script>" not in failed["body"] and b"&lt;script&gt;" in failed["body"]
    (folder / "report.html").write_text("<html><body>complete</body></html>")
    complete = request("GET", f"/results/{token}/report.html")
    assert complete["status"] == "200 OK" and b"complete" in complete["body"]
    assert b"Retry enrichment" not in complete["body"]


def test_health_is_read_only_and_reports_pbs_only(monkeypatch):
    monkeypatch.setattr(
        ood_app.subprocess, "run", lambda *_a, **_kw: pytest.fail("health must not run tools")
    )
    response = request("GET", "/health")
    payload = json.loads(response["body"])
    assert payload["report_execution"] == "pbs_only"
    assert payload["pbs_submit_command_exists"] is True
    assert payload["active_jobs_url"] == "/pun/sys/dashboard/activejobs"
    assert payload["uniprot"] == "deferred_to_compute_job"


def test_every_report_figure_is_allowed_and_path_traversal_is_blocked():
    report_source = Path("src/protein_agent/report.py").read_text()
    sections_source = Path("src/protein_agent/report_sections.py").read_text()
    served = set(re.findall(r'(?:src|href)="([a-z]+\.png)"', report_source + sections_source))
    assert served and served <= set(ood_app.FILES)
    traversal = request("GET", "/results/" + "a" * 32 + "/../../input.fasta")
    assert traversal["status"] == "404 Not Found"


def test_local_prototype_serves_every_report_figure(tmp_path):
    from protein_agent.web import make_handler

    source = Path("src/protein_agent/report.py").read_text()
    source += Path("src/protein_agent/report_sections.py").read_text()
    served = set(re.findall(r'(?:src|href)="([a-z]+\.png)"', source))
    token = "c" * 32
    folder = tmp_path / token
    folder.mkdir()
    for name in served:
        (folder / name).write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 200)
    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            for name in sorted(served):
                url = f"http://127.0.0.1:{server.server_port}/results/{token}/{name}"
                with urlopen(url) as response:
                    assert response.status == 200
        finally:
            server.shutdown()
            thread.join()
