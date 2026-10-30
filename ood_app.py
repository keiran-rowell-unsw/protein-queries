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
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Protein report</title>
<style>
:root{color-scheme:light;--ink:#172033;--muted:#64748b;--line:#dbe4ee;--brand:#155e75;--brand2:#0e7490}
*{box-sizing:border-box}body{margin:0;min-height:100vh;padding:3rem 1rem;font:1rem/1.55 system-ui,-apple-system,sans-serif;color:var(--ink);background:linear-gradient(145deg,#eef7f8,#f8fafc 48%%,#eef2ff)}
.shell{width:min(760px,100%%);margin:auto}.hero{margin:0 0 1.5rem}.eyebrow{margin:0 0 .4rem;color:var(--brand2);font-size:.78rem;font-weight:750;letter-spacing:.12em;text-transform:uppercase}h1{margin:0;font-size:clamp(2rem,7vw,3.4rem);line-height:1.05;letter-spacing:-.04em}.subtitle{margin:.8rem 0 0;color:var(--muted);font-size:1.1rem}
.panel{padding:clamp(1.25rem,4vw,2.2rem);border:1px solid rgba(255,255,255,.8);border-radius:22px;background:rgba(255,255,255,.9);box-shadow:0 24px 70px rgba(30,64,175,.12);backdrop-filter:blur(12px)}form{display:grid;gap:1.25rem}.field{display:grid;gap:.45rem}label{font-weight:700}.hint{color:var(--muted);font-size:.88rem}input,textarea{width:100%%;border:1px solid var(--line);border-radius:12px;background:#fff;color:var(--ink);font:inherit;padding:.8rem .9rem;transition:border-color .15s,box-shadow .15s}textarea{min-height:13rem;resize:vertical;font-family:ui-monospace,SFMono-Regular,monospace;font-size:.9rem}input:focus,textarea:focus{outline:none;border-color:#0891b2;box-shadow:0 0 0 4px rgba(8,145,178,.13)}button{justify-self:start;border:0;border-radius:999px;padding:.78rem 1.25rem;background:linear-gradient(135deg,var(--brand),var(--brand2));color:white;font:inherit;font-weight:750;cursor:pointer;box-shadow:0 8px 20px rgba(14,116,144,.22)}button:hover{filter:brightness(1.06)}[role=alert]{padding:.8rem 1rem;border-radius:10px;background:#fff1f2;color:#9f1239}
</style></head><body><div class="shell"><header class="hero"><p class="eyebrow">Sequence · annotation · structure</p><h1>Protein report</h1><p class="subtitle">Turn a protein identifier or sequence into a clear, downloadable report.</p></header><main class="panel">
%s<form method="post" action="%s/run"><input type="hidden" name="csrf" value="%s">
<div class="field"><label for="accession">UniProt accession</label><span class="hint">Recommended when you know it.</span><input id="accession" name="accession" maxlength="13" autocomplete="off" placeholder="e.g. P05067" pattern="[A-Za-z0-9]{6,10}(-[0-9]{1,2})?"></div>
<div class="field"><label for="sequence">Protein sequence</label><span class="hint">Optional with an accession. FASTA and bare sequences are accepted.</span><textarea id="sequence" name="sequence" rows="12" maxlength="100000" spellcheck="false" placeholder=">protein_name&#10;MKT…"></textarea></div>
<button type="submit">Generate report</button>
</form></main></div></body></html>""" % (message, prefix, html.escape(csrf, quote=True))


def chat_page(prefix, token, csrf, history=(), message=""):
    action = "%s/results/%s/chat" % (html.escape(prefix, quote=True), token)
    turns = "".join(
        "<section class=\"turn %s\"><h2>%s</h2><p>%s</p></section>" % (
            "user" if item.get("role") == "user" else "assistant",
            "You" if item.get("role") == "user" else "Report assistant",
            html.escape(str(item.get("content") or "")),
        )
        for item in history
        if item.get("role") in {"user", "assistant"}
    )
    return """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Ask about this report</title>
<style>:root{--ink:#172033;--muted:#64748b;--line:#dbe4ee;--brand:#0e7490}*{box-sizing:border-box}body{margin:0;min-height:100vh;padding:2.5rem 1rem;font:1rem/1.55 system-ui;color:var(--ink);background:#f1f6f8}.shell{max-width:820px;margin:auto}a{color:var(--brand)}header{margin-bottom:1.5rem}h1{margin:.4rem 0;font-size:clamp(2rem,6vw,3rem);letter-spacing:-.035em}header p{color:var(--muted)}.turn{margin:1rem 0;padding:1rem 1.2rem;border:1px solid var(--line);border-radius:16px;background:white;box-shadow:0 8px 25px rgba(15,23,42,.05)}.turn.user{margin-left:12%%;background:#e6f5f7}.turn h2{margin:0 0 .35rem;font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--brand)}.turn p{margin:0;white-space:pre-wrap}.composer{margin-top:1.5rem;padding:1.2rem;border-radius:18px;background:white;box-shadow:0 12px 35px rgba(15,23,42,.08)}label{display:block;margin-bottom:.4rem;font-weight:700}textarea{width:100%%;box-sizing:border-box;resize:vertical;border:1px solid var(--line);border-radius:12px;padding:.8rem;font:inherit}textarea:focus{outline:none;border-color:var(--brand);box-shadow:0 0 0 4px rgba(14,116,144,.12)}button{margin-top:.75rem;border:0;border-radius:999px;padding:.7rem 1.1rem;background:var(--brand);color:white;font:inherit;font-weight:700;cursor:pointer}</style></head><body><div class="shell">
<header><a href="%s/results/%s/report.html">← Back to report</a><h1>Ask about this report</h1><p>Answers are restricted to facts recorded in this report.</p></header>%s%s
<form class="composer" method="post" action="%s"><input type="hidden" name="csrf" value="%s">
<label for="question">Your question</label><textarea id="question" name="question" rows="4"
maxlength="4000" required placeholder="What are the main findings?"></textarea><button type="submit">Ask report assistant</button>
</form></div></body></html>""" % (
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
                        b'<div class="report-shell"><div class="note"><h2>Ask about this report</h2>'
                        b'<p><a href="chat">Open the report-only assistant</a>. It can explain '
                        b'this report but cannot perform fresh database retrieval.</p></div></div>'
                    )
                    body = body.replace(b"</body>", link + b"</body>")
                return respond(start_response, "200 OK", body, FILES[parts[3]])
    return respond(start_response, "404 Not Found", "Not found")
