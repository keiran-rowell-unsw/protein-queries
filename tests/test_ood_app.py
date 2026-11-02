"""Passenger WSGI contract without an OOD server or remote SSH."""

import json
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
        generate_report(folder / "input.fasta", folder, lightweight=True)
        (folder / "disorder.png").write_bytes(b"test image")

    queued = []
    monkeypatch.setattr(ood_app, "run_report", compute)
    monkeypatch.setattr(ood_app, "submit_extended_analysis", lambda folder: queued.append(folder))
    form = request("GET", "/")
    assert form["status"] == "200 OK"
    assert b"/pun/dev/protein-queries/run" in form["body"]
    assert b"Supply a UniProt accession" not in form["body"]
    assert b"Turn a protein identifier or sequence" in form["body"]
    assert b"UniProt accession</label><span class=\"hint\">Queues a job to fetch results" in form["body"]
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    params = {"csrf": token, "sequence": ">example\nACDEFGHIKLMNPQRSTVWY"}
    submit = request("POST", "/run", body=urlencode(params).encode(), cookie=cookie)
    assert submit["status"] == "303 See Other"
    url = submit["headers"]["Location"]
    assert url.startswith("/pun/dev/protein-queries/results/")
    report = request("GET", url.removeprefix("/pun/dev/protein-queries"))
    assert report["status"] == "200 OK"
    assert b"Protein analysis report" in report["body"]
    assert b"Preliminary sequence-only report" in report["body"]
    assert len(queued) == 1
    assert (queued[0] / "input.fasta").is_file()
    assert b"Open the report-only assistant" not in report["body"]
    image = request(
        "GET", url.removeprefix("/pun/dev/protein-queries").replace("report.html", "disorder.png")
    )
    assert image["status"] == "200 OK"
    assert image["headers"]["Content-Type"] == "image/png"
    assert (
        request("GET", "/results/" + "a" * 32 + "/../../input.fasta")["status"] == "404 Not Found"
    )


def test_accession_only_submission_is_queued(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    captured = {}

    def fake_queue(folder):
        captured["request"] = json.loads((folder / "request.json").read_text())
        (folder / "job.json").write_text('{"status":"queued","job_id":"123.server"}')
        return {"status": "queued", "job_id": "123.server"}

    monkeypatch.setattr(ood_app, "submit_extended_analysis", fake_queue)
    form = request("GET", "/")
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    params = {"csrf": token, "accession": "Q12345"}
    result = request("POST", "/run", body=urlencode(params).encode(), cookie=cookie)
    assert result["status"] == "303 See Other"
    assert captured == {"request": {"accession": "Q12345"}}


def test_accession_only_failed_job_does_not_show_endless_pending_page(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    token = "a" * 32
    folder = tmp_path / token
    folder.mkdir()
    (folder / "job.json").write_text('{"status":"failed"}')
    result = request("GET", f"/results/{token}/report.html")
    assert result["status"] == "200 OK"
    assert b"Analysis unavailable" in result["body"]
    assert b"refresh" not in result["body"]


def test_sequence_and_accession_reports_immediately_then_queues(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    observed = {}

    def fake_report(sequence, folder, accession=""):
        observed["sequence"] = sequence
        observed["accession"] = accession
        (folder / "report.html").write_text("preliminary")
        (folder / "report.json").write_text('{"uniprot":{"status":"not_requested"}}')

    def fake_queue(folder):
        observed["queued_with_report"] = (folder / "report.html").is_file()
        observed["saved_input"] = (folder / "input.fasta").read_text()
        return {"status": "queued", "job_id": "123.server"}

    monkeypatch.setattr(ood_app, "run_report", fake_report)
    monkeypatch.setattr(ood_app, "submit_extended_analysis", fake_queue)
    form = request("GET", "/")
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    result = request("POST", "/run", body=urlencode({
        "csrf": token, "sequence": ">submitted\nACDEFGHIKLMNPQRSTVWY", "accession": "Q12345",
    }).encode(), cookie=cookie)
    assert result["status"] == "303 See Other"
    assert observed == {
        "sequence": ">submitted\nACDEFGHIKLMNPQRSTVWY",
        "accession": "Q12345", "queued_with_report": True,
        "saved_input": ">submitted\nACDEFGHIKLMNPQRSTVWY",
    }


def test_preliminary_report_refreshes_only_while_enrichment_is_queued(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)
    token = "b" * 32
    folder = tmp_path / token
    folder.mkdir()
    (folder / "report.html").write_text("<html><body>preliminary</body></html>")
    job = folder / "job.json"
    job.write_text('{"status":"queued"}')
    queued = request("GET", f"/results/{token}/report.html")
    assert queued["headers"]["Refresh"] == "10"
    job.write_text('{"status":"complete"}')
    complete = request("GET", f"/results/{token}/report.html")
    assert "Refresh" not in complete["headers"]


def test_queue_failure_keeps_preliminary_report(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "ROOT", tmp_path)

    def fake_report(_sequence, folder, _accession=""):
        (folder / "report.html").write_text("<html><body>preliminary</body></html>")
        (folder / "report.json").write_text("{}")

    monkeypatch.setattr(ood_app, "run_report", fake_report)
    monkeypatch.setattr(ood_app, "submit_extended_analysis", lambda _folder: (_ for _ in ()).throw(
        RuntimeError("compute submission unavailable")
    ))
    form = request("GET", "/")
    cookie = form["headers"]["Set-Cookie"]
    token = cookie.split(";", 1)[0].split("=", 1)[1]
    result = request("POST", "/run", body=urlencode({
        "csrf": token, "sequence": "ACDEFGHIKLMNPQRSTVWY",
    }).encode(), cookie=cookie)
    assert result["status"] == "303 See Other"
    report = request("GET", result["headers"]["Location"].removeprefix(
        "/pun/dev/protein-queries"
    ))
    assert b"preliminary" in report["body"]
    assert b"Retry enrichment" in report["body"]
    folder = next(tmp_path.iterdir())
    assert json.loads((folder / "job.json").read_text())["status"] == "unavailable"
    assert (folder / "input.fasta").is_file()


def test_report_chat_is_temporarily_disabled(tmp_path):
    token = "a" * 32
    folder = tmp_path / token
    folder.mkdir()
    (folder / "report.json").write_text("{}")
    assert request("GET", f"/results/{token}/chat")["status"] == "404 Not Found"
    assert request("POST", f"/results/{token}/chat", body=b"")["status"] == "404 Not Found"


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


def test_legacy_ssh_mode_cannot_run_report_on_login_node(monkeypatch, tmp_path):
    monkeypatch.setattr(ood_app, "EXECUTION_MODE", "ssh")
    monkeypatch.setattr(ood_app.subprocess, "run", lambda *_args, **_kwargs: pytest.fail(
        "must not launch SSH"
    ))
    with pytest.raises(RuntimeError, match="SSH/login nodes are disabled"):
        ood_app.run_report("ACDE", tmp_path)



def test_local_run_avoids_ssh(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    seen = []

    def fake_run(args, **kwargs):
        seen.append((args, kwargs))
        (folder / "report.html").write_text("generated")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(ood_app, "EXECUTION_MODE", "local")
    monkeypatch.setattr(ood_app.subprocess, "run", fake_run)
    ood_app.run_report("ACDE", folder)
    args, kwargs = seen[0]
    assert args[0].endswith("/.venv/bin/python")
    assert "ssh" not in args
    assert args[-1] == "--lightweight"
    assert kwargs["timeout"] == 90
    assert (folder / "input.fasta").is_file()
    ood_app.run_report("ACDE", folder, "Q12345")
    assert seen[-1][0][-2:] == ["--accession", "Q12345"]


def test_extended_analysis_submits_pbs_script(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    (folder / "report.json").write_text("{}")
    (folder / "input.fasta").write_text(">x\nACDEFGHIKLMNPQRSTVWY")
    (folder / "request.json").write_text(json.dumps({"accession": "Q12345"}))
    seen = {}

    def fake_run(args, **kwargs):
        seen["args"] = args
        return SimpleNamespace(returncode=0, stdout="123.server\n")

    monkeypatch.setattr(ood_app, "EXTENDED_PYTHON", "/compute/project/.venv/bin/python")
    monkeypatch.setattr(ood_app.subprocess, "run", fake_run)
    result = ood_app.submit_extended_analysis(folder)
    assert result == {"status": "queued", "job_id": "123.server"}
    assert seen["args"] == ["qsub", str(folder / ".extended.pbs")]
    assert (folder / "job-id.txt").read_text().strip() == "123.server"
    script = (folder / ".extended.pbs").read_text()
    assert "--foldseek" in script and "Q12345" in script
    assert "ssh" not in script


def test_remote_report_chat_uses_stdin_not_shell_arguments(monkeypatch, tmp_path):
    folder = tmp_path / "report"
    folder.mkdir()
    (folder / "report.json").write_text("{}")
    seen = {}

    def fake_run(args, **kwargs):
        seen.update(args=args, kwargs=kwargs)
        return SimpleNamespace(returncode=0, stdout=json.dumps({"answer": "Grounded."}))

    monkeypatch.setattr(ood_app.subprocess, "run", fake_run)
    answer = ood_app.run_report_chat(folder, "question; rm -rf /", [])
    assert answer == "Grounded."
    assert "question; rm -rf /" not in seen["args"][-1]
    assert json.loads(seen["kwargs"]["input"])["question"] == "question; rm -rf /"
    assert "protein_agent.report_chat" in seen["args"][-1]




def test_every_figure_the_report_emits_is_an_allowed_file():
    """The report and the file allowlist are edited together or the page shows broken images.

    Both halves of the OOD app gate a request through a literal filename set, and each
    figure also needs a Content-Type. A name missing from either is a 404 or an inline
    blob rather than an error, so this checks the two lists against the report's own
    output instead of against a copy of the same list.
    """
    import re
    from pathlib import Path

    report_source = Path("src/protein_agent/report.py").read_text()
    sections_source = Path("src/protein_agent/report_sections.py").read_text()
    served = set(re.findall(r'(?:src|href)="([a-z]+\.png)"', report_source + sections_source))
    assert served, "no figure is referenced by the report, so this test has stopped testing"
    assert served <= set(ood_app.FILES), served - set(ood_app.FILES)


def test_the_local_handler_serves_every_figure_the_report_emits(tmp_path):
    """Same check against web.py's separate allowlist and its separate MIME map.

    The two lists are written out twice in that handler; a name present in one and
    absent from the other is a KeyError in the response path rather than a clean 404.
    """
    import re
    from http.server import ThreadingHTTPServer
    from pathlib import Path
    from threading import Thread
    from urllib.request import urlopen

    from protein_agent.web import make_handler

    report_source = Path("src/protein_agent/report.py").read_text()
    sections_source = Path("src/protein_agent/report_sections.py").read_text()
    served = set(re.findall(r'(?:src|href)="([a-z]+\.png)"', report_source + sections_source))
    assert {"hydropathy.png", "positions.png"} <= served

    token = "c" * 32
    folder = tmp_path / token
    folder.mkdir()
    for name in served:
        (folder / name).write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 200)

    with ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base = f"http://127.0.0.1:{server.server_port}"
        try:
            for name in sorted(served):
                with urlopen(f"{base}/results/{token}/{name}") as response:
                    assert response.status == 200, name
                    assert response.headers["Content-Type"] == "image/png", name
        finally:
            server.shutdown()
            thread.join()
