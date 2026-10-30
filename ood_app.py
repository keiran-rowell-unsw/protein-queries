"""Open OnDemand Passenger WSGI front end; computes on katana3 via per-user SSH.

No queue, no network listener on the compute node, no inference key. Inputs and results
live in the user's shared home directory; the OOD PUN provides authentication.
"""

# HTML template and WSGI responses are kept readable as continuous lines.
# ruff: noqa: E501, UP031
import fcntl
import hmac
import html
import json
import os
import re
import secrets
import shlex
import shutil
import subprocess
from pathlib import Path
from urllib.parse import parse_qs

HOST = os.environ.get("PROTEIN_REPORT_SSH_HOST", "katana3")
PROJECT = os.environ.get("PROTEIN_REPORT_PROJECT", "/srv/scratch/sbf/protein_queries")
ROOT = Path.home() / "ondemand" / "data" / "protein-queries" / "reports"
MAX_BODY = 300_000
MAX_CHAT_BODY = 20_000
MAX_CHAT_MESSAGES = 12
TOKEN = re.compile(r"[a-f0-9]{32}\Z")
FILES = {"report.html": "text/html; charset=utf-8", "report.json": "application/json",
         "hydropathy.png": "image/png", "disorder.png": "image/png"}
CSP = "default-src 'none'; img-src 'self'; style-src 'unsafe-inline'; form-action 'self'"


def respond(start_response, code, body, content_type="text/html; charset=utf-8", headers=()):
    if isinstance(body, str):
        body = body.encode("utf-8")
    start_response(code, [
        ("Content-Type", content_type), ("Content-Length", str(len(body))),
        ("Cache-Control", "no-store"), ("X-Content-Type-Options", "nosniff"),
        ("Content-Security-Policy", CSP), *headers,
    ])
    return [body]


def page(prefix, csrf, message=""):
    prefix = html.escape(prefix, quote=True)
    return """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Protein Sequence Report</title>
<style>body{font:1rem/1.5 system-ui;max-width:720px;margin:2rem auto;padding:1rem}
textarea{width:100%%;box-sizing:border-box}button{padding:.6rem 1rem}</style></head><body>
<h1>Protein Sequence Report</h1><p>Supply a UniProt accession, one FASTA record, or a bare
protein sequence (up to 10,000 standard amino acids). An accession takes precedence when
both are supplied. Without one, the sequence is matched exactly to UniProt. A verified
accession is used to retrieve experimental structure mappings from PDBe. No PBS submission
is made.</p>
%s<form method="post" action="%s/run"><input type="hidden" name="csrf" value="%s">
<label for="sequence">Sequence (optional if accession supplied)</label><br><textarea id="sequence" name="sequence" rows="12"
maxlength="100000"></textarea><br>
<label for="accession">UniProt accession (optional; takes precedence)</label><br>
<input id="accession" name="accession" maxlength="13" pattern="[A-Za-z0-9]{6,10}(-[0-9]{1,2})?"><br>
<button type="submit">Generate report</button>
</form></body></html>""" % (message, prefix, html.escape(csrf, quote=True))


def chat_page(prefix, token, csrf, history=(), message=""):
    action = "%s/results/%s/chat" % (html.escape(prefix, quote=True), token)
    turns = "".join(
        "<section><h2>%s</h2><p>%s</p></section>" % (
            "You" if item.get("role") == "user" else "Report assistant",
            html.escape(str(item.get("content") or "")),
        )
        for item in history
        if item.get("role") in {"user", "assistant"}
    )
    return """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Ask about this report</title>
<style>body{font:1rem/1.5 system-ui;max-width:800px;margin:2rem auto;padding:1rem}
textarea{width:100%%;box-sizing:border-box}.assistant{background:#eef2f7;padding:1rem}</style></head><body>
<h1>Ask about this report</h1><p>This assistant can read only this generated report. It cannot
perform fresh UniProt or PDBe retrieval. Every answer re-reads report data through a required tool.</p>
<p><a href="%s/results/%s/report.html">Back to report</a></p>%s%s
<form method="post" action="%s"><input type="hidden" name="csrf" value="%s">
<label for="question">Question</label><br><textarea id="question" name="question" rows="4"
maxlength="4000" required></textarea><br><button type="submit">Ask report assistant</button>
</form></body></html>""" % (
        html.escape(prefix, quote=True), token, message, turns, action,
        html.escape(csrf, quote=True),
    )


def _load_chat(folder):
    path = folder / "chat.json"
    if not path.is_file():
        return []
    try:
        items = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return []
    return items[-MAX_CHAT_MESSAGES:] if isinstance(items, list) else []


def _save_chat(folder, history):
    target = folder / "chat.json"
    temporary = folder / (".chat.%s.tmp" % secrets.token_hex(8))
    try:
        temporary.write_text(
            json.dumps(history[-MAX_CHAT_MESSAGES:], ensure_ascii=False), encoding="utf-8"
        )
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def run_report_chat(folder, question, history):
    """Run report-only inference remotely; user text travels on stdin, never in the shell."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.-]*", HOST) or not PROJECT.startswith("/"):
        raise RuntimeError("invalid report-chat host configuration")
    args = [
        PROJECT + "/.venv/bin/python", "-m", "protein_agent.report_chat",
        "--report", str(folder / "report.json"),
    ]
    remote = "cd " + shlex.quote(PROJECT) + " && " + shlex.join(args)
    payload = json.dumps({"question": question, "history": history}, ensure_ascii=False)
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=yes",
         "-o", "ConnectTimeout=8", HOST, remote],
        input=payload, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        text=True, timeout=150, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError("report assistant unavailable")
    try:
        response = json.loads(result.stdout)
    except (ValueError, TypeError) as exc:
        raise RuntimeError("report assistant returned an invalid response") from exc
    answer = response.get("answer")
    if not isinstance(answer, str) or not answer.strip():
        raise RuntimeError("report assistant could not produce a grounded answer")
    return answer.strip()


def validate_sequence(text):
    """Validate locally so malformed user input is a 400, not an SSH failure."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("no sequence")
    if lines[0].startswith(">"):
        if not lines[0][1:].strip():
            raise ValueError("empty FASTA header")
        lines.pop(0)
    if not lines or any(line.startswith(">") for line in lines):
        raise ValueError("provide exactly one FASTA record")
    sequence = "".join("".join(line.split()) for line in lines).upper()
    if not 1 <= len(sequence) <= 10_000 or re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence) is None:
        raise ValueError("use 1-10,000 standard amino-acid letters with no gaps")


def run_report(sequence, folder, accession=""):
    """SSH as the PUN user. Paths are generated by this app, never supplied by users."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.-]*", HOST):
        raise RuntimeError("invalid SSH host configuration")
    if not PROJECT.startswith("/"):
        raise RuntimeError("project path must be absolute")
    source = folder / "input.fasta"
    source.write_text(sequence, encoding="utf-8")
    try:
        args = [PROJECT + "/.venv/bin/python", "-m", "protein_agent.report",
                "--input", str(source), "--output", str(folder)]
        if accession:
            args.extend(["--accession", accession])
        result = subprocess.run(
            ["ssh", "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=yes",
             "-o", "ConnectTimeout=8", HOST, shlex.join(args)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            timeout=90, check=False,
        )
        if result.returncode != 0 or not (folder / "report.html").is_file():
            raise RuntimeError("analysis failed; check OOD-to-katana3 SSH and the shared environment")
    finally:
        source.unlink(missing_ok=True)


def application(environ, start_response):
    """Small per-user WSGI app; paths work under OOD's SCRIPT_NAME prefix."""
    method = environ.get("REQUEST_METHOD", "GET")
    path = environ.get("PATH_INFO", "/")
    prefix = environ.get("SCRIPT_NAME", "").rstrip("/")
    if method == "GET" and path in ("/", ""):
        csrf = secrets.token_hex(16)
        cookie_path = prefix + "/"
        return respond(start_response, "200 OK", page(prefix, csrf), headers=[
            ("Set-Cookie", "protein_report_csrf=%s; Path=%s; Secure; HttpOnly; SameSite=Strict" %
             (csrf, cookie_path)),
        ])
    chat_parts = path.split("/")
    is_chat = (
        len(chat_parts) == 4
        and chat_parts[1] == "results"
        and TOKEN.fullmatch(chat_parts[2])
        and chat_parts[3] == "chat"
    )
    if method == "GET" and is_chat:
        folder = ROOT / chat_parts[2]
        if not (folder / "report.json").is_file():
            return respond(start_response, "404 Not Found", "Not found")
        csrf = secrets.token_hex(16)
        cookie_path = prefix + "/"
        return respond(
            start_response,
            "200 OK",
            chat_page(prefix, chat_parts[2], csrf, _load_chat(folder)),
            headers=[
                ("Set-Cookie", "protein_report_csrf=%s; Path=%s; Secure; HttpOnly; "
                 "SameSite=Strict" % (csrf, cookie_path)),
            ],
        )
    if method == "POST" and is_chat:
        folder = ROOT / chat_parts[2]
        if not (folder / "report.json").is_file():
            return respond(start_response, "404 Not Found", "Not found")
        try:
            size = int(environ.get("CONTENT_LENGTH") or 0)
            if size < 1 or size > MAX_CHAT_BODY:
                raise ValueError("invalid chat submission size")
            fields = parse_qs(
                environ["wsgi.input"].read(size).decode("utf-8"), keep_blank_values=True
            )
            cookies = dict(
                part.strip().split("=", 1)
                for part in environ.get("HTTP_COOKIE", "").split(";")
                if "=" in part
            )
            stored = cookies.get("protein_report_csrf", "")
            submitted = fields.get("csrf", [])
            if (
                not TOKEN.fullmatch(stored)
                or len(submitted) != 1
                or not hmac.compare_digest(stored, submitted[0])
            ):
                return respond(start_response, "403 Forbidden", "Invalid form token. Reload chat.")
            if set(fields) - {"csrf", "question"} or len(fields.get("question", [])) != 1:
                raise ValueError("submit exactly one question")
            question = fields["question"][0].strip()
            if not 1 <= len(question) <= 4000:
                raise ValueError("question must contain 1-4,000 characters")
            with (folder / ".chat.lock").open("a", encoding="utf-8") as lock:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
                history = _load_chat(folder)
                answer = run_report_chat(folder, question, history)
                history.extend([
                    {"role": "user", "content": question},
                    {"role": "assistant", "content": answer},
                ])
                _save_chat(folder, history)
            return respond(
                start_response,
                "200 OK",
                chat_page(prefix, chat_parts[2], stored, history),
            )
        except (ValueError, UnicodeError):
            return respond(start_response, "400 Bad Request", "Invalid chat submission.")
        except (OSError, RuntimeError, subprocess.TimeoutExpired):
            return respond(
                start_response,
                "503 Service Unavailable",
                "<p>The report assistant is unavailable. The report itself remains valid.</p>",
            )
    if method == "POST" and path == "/run":
        stored = ""
        try:
            size = int(environ.get("CONTENT_LENGTH") or 0)
            if size < 1 or size > MAX_BODY:
                raise ValueError("submission is empty or too large")
            fields = parse_qs(environ["wsgi.input"].read(size).decode("utf-8"), keep_blank_values=True)
            cookies = dict(part.strip().split("=", 1) for part in
                           environ.get("HTTP_COOKIE", "").split(";") if "=" in part)
            submitted = fields.get("csrf", [])
            stored = cookies.get("protein_report_csrf", "")
            if (not TOKEN.fullmatch(stored) or len(submitted) != 1 or
                    not hmac.compare_digest(stored, submitted[0])):
                return respond(start_response, "403 Forbidden", "Invalid form token. Reload the app.")
            if (set(fields) - {"csrf", "sequence", "accession"}
                    or len(fields.get("sequence", [])) > 1
                    or len(fields.get("accession", [])) > 1):
                raise ValueError("submit at most one sequence and one accession")
            text = (fields.get("sequence") or [""])[0]
            accession = (fields.get("accession") or [""])[0].strip().upper()
            if accession and not re.fullmatch(r"[A-Z0-9]{6,10}(?:-[0-9]{1,2})?", accession):
                raise ValueError("invalid UniProt accession")
            if len(text.encode("utf-8")) > 100_000:
                raise ValueError("sequence input exceeds 100 KB")
            if not text.strip() and not accession:
                raise ValueError("provide either a protein sequence or a UniProt accession")
            if text.strip() and not accession:
                validate_sequence(text)
            # The remote parser rechecks the same input before computation.
            token = secrets.token_hex(16)
            ROOT.mkdir(mode=0o700, parents=True, exist_ok=True)
            folder = ROOT / token
            folder.mkdir(mode=0o700)
            try:
                run_report(text, folder, accession)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                shutil.rmtree(folder)
                raise
            return respond(start_response, "303 See Other", b"", headers=[
                ("Location", prefix + "/results/" + token + "/report.html"),
            ])
        except (ValueError, UnicodeError):
            return respond(start_response, "400 Bad Request",
                           page(prefix, stored if TOKEN.fullmatch(stored) else "",
                                "<p role='alert'>Invalid sequence submission. Check the input format and size.</p>"))
        except (OSError, RuntimeError, subprocess.TimeoutExpired):
            return respond(start_response, "503 Service Unavailable",
                           "<p>Report unavailable. Ask the administrator to check SSH access from "
                           "the OOD host to katana3 and the shared project environment.</p>")
    if method == "GET" and path.startswith("/results/"):
        parts = path.split("/")
        if len(parts) == 4 and TOKEN.fullmatch(parts[2]) and parts[3] in FILES:
            result = ROOT / parts[2] / parts[3]
            if result.is_file():
                body = result.read_bytes()
                if parts[3] == "report.html":
                    link = (
                        b'<div class="note"><h2>Ask about this report</h2>'
                        b'<p><a href="chat">Open the report-only assistant</a>. It can explain '
                        b'this report but cannot perform fresh database retrieval.</p></div>'
                    )
                    body = body.replace(b"</body>", link + b"</body>")
                return respond(start_response, "200 OK", body, FILES[parts[3]])
    return respond(start_response, "404 Not Found", "Not found")
