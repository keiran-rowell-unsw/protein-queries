"""Open OnDemand Passenger front end for reports executed entirely through PBS.

The PUN validates and stores input, submits one job, and serves status/results. All protein
calculation and retrieval runs on an allocated compute node.
"""

# HTML template and WSGI responses are kept readable as continuous lines.
# ruff: noqa: E501, UP031
import contextlib
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
import sys
from pathlib import Path
from urllib.parse import parse_qs

HOST = os.environ.get("PROTEIN_REPORT_SSH_HOST", "katana3")
PROJECT = os.environ.get("PROTEIN_REPORT_PROJECT", str(Path(__file__).resolve().parent))
PBS_SUBMIT_COMMAND = (
    os.environ.get("PBS_SUBMIT_COMMAND") or shutil.which("qsub") or "/opt/pbs/bin/qsub"
)
PBS_STATUS_COMMAND = (
    os.environ.get("PBS_STATUS_COMMAND") or shutil.which("qstat") or "/opt/pbs/bin/qstat"
)
EXTENDED_PYTHON = os.environ.get("PROTEIN_REPORT_EXTENDED_PYTHON", "")
DEFAULT_COMPUTE_PYTHON = str(Path(PROJECT) / ".compute-venv" / "bin" / "python")
OOD_ACTIVE_JOBS_URL = os.environ.get(
    "OOD_ACTIVE_JOBS_URL", "/pun/sys/dashboard/activejobs"
)
REPORT_EMAIL = os.environ.get("PROTEIN_REPORT_EMAIL", "").strip()
EMAIL_RE = re.compile(r"[^\s@<>]+@[^\s@<>]+\.[^\s@<>]+\Z")
ROOT = Path.home() / "ondemand" / "data" / "protein-queries" / "reports"
MAX_BODY = 300_000
MAX_CHAT_BODY = 20_000
MAX_CHAT_MESSAGES = 12
MAX_SEQUENCE_RESIDUES = 5_000
# Temporarily disabled while the retrieval interface is being expanded.
REPORT_ASSISTANT_ENABLED = False
TOKEN = re.compile(r"[a-f0-9]{32}\Z")
FILES = {"report.html": "text/html; charset=utf-8", "report.json": "application/json",
         "hydropathy.png": "image/png", "disorder.png": "image/png",
         "positions.png": "image/png", "structure.png": "image/png",
         "ligands.png": "image/png"}
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
<div class="field"><label for="accession">UniProt accession</label><span class="hint">Queues a job to fetch results</span><input id="accession" name="accession" maxlength="13" autocomplete="off" placeholder="e.g. P05067" pattern="([A-Za-z0-9]{6}|[A-Za-z0-9]{10})(-[0-9]{1,2})?"></div>
<div class="field"><label for="email">Email report (optional)</label><span class="hint">We will send the completed HTML and JSON report as attachments.</span><input id="email" name="email" type="email" maxlength="254" autocomplete="email" placeholder="you@example.edu"></div>
<div class="field"><label for="sequence">Protein sequence</label><span class="hint">Queues a full report job. FASTA or bare sequence; maximum 5,000 residues (header and whitespace excluded).</span><textarea id="sequence" name="sequence" rows="12" maxlength="100000" spellcheck="false" placeholder=">protein_name&#10;MKT…"></textarea></div>
<button type="submit">Queue report</button>
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
    """Validate locally so malformed user input is rejected before PBS submission."""
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
    if len(sequence) > MAX_SEQUENCE_RESIDUES:
        raise ValueError("sequence exceeds the 5,000-residue input limit")
    if not sequence or re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", sequence) is None:
        raise ValueError("use standard amino-acid letters with no gaps")


def _write_json(path, payload):
    temporary = path.with_name("." + path.name + ".tmp")
    try:
        temporary.write_text(json.dumps(payload), encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def submit_report_job(folder):
    """Submit the complete report to PBS; Passenger performs no report calculations."""
    job_path = folder / "job.json"
    if job_path.is_file():
        try:
            existing = json.loads(job_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            existing = {}
        if existing.get("status") in {"queued", "running", "complete"}:
            return existing
    if not (folder / "input.fasta").is_file():
        raise RuntimeError("the retained input sequence is unavailable")
    if not PROJECT.startswith("/"):
        raise RuntimeError("invalid project configuration")
    pbs_command = Path(PBS_SUBMIT_COMMAND)
    if pbs_command.is_absolute():
        if not pbs_command.is_file() or not os.access(pbs_command, os.X_OK):
            raise RuntimeError("configured PBS submission command is unavailable")
    elif not re.fullmatch(r"[A-Za-z0-9_.-]+", PBS_SUBMIT_COMMAND):
        raise RuntimeError("invalid PBS submission command")
    request = {}
    with contextlib.suppress(OSError, ValueError, TypeError):
        request = json.loads((folder / "request.json").read_text(encoding="utf-8"))
    report_email = str(request.get("email") or REPORT_EMAIL).strip()
    if report_email and not EMAIL_RE.fullmatch(report_email):
        raise RuntimeError("report email is not a valid address")
    python = EXTENDED_PYTHON or DEFAULT_COMPUTE_PYTHON
    if not python.startswith("/"):
        raise RuntimeError("extended Python path must be absolute")
    if not EXTENDED_PYTHON and not Path(python).is_file():
        raise RuntimeError(
            "compute environment is missing; run bin/setup-compute-env on a compute node"
        )
    command = [
        python, "-m", "protein_agent.report", "--input", str(folder / "input.fasta"),
        "--output", str(folder), "--foldseek",
    ]
    accession = request.get("accession")
    if isinstance(accession, str) and accession:
        command.extend(["--accession", accession])
    script = folder / ".report.pbs"
    reservation = folder / ".submission-reserved"
    job_json = str(job_path)
    script.write_text(
        "#!/bin/bash\n"
        "#PBS -N protein-report\n"
        "#PBS -l walltime=00:20:00\n"
        "#PBS -l ncpus=2\n"
        "#PBS -l mem=4gb\n"
        "#PBS -m abe\n"
        "set -u\n"
        f"printf '%s\\n' '{{\"status\":\"running\"}}' > {shlex.quote(job_json)}\n"
        f"if {shlex.join(command)}; then\n"
        + (
            "  " + shlex.join([
                python, "-m", "protein_agent.mail_report", "--report-dir", str(folder),
                "--recipient", report_email,
            ]) + " || true\n"
            if report_email else ""
        )
        + f"  printf '%s\\n' '{{\"status\":\"complete\"}}' > {shlex.quote(job_json)}\n"
        f"  rm -f -- {shlex.quote(str(reservation))}\n"
        "else\n"
        f"  printf '%s\\n' '{{\"status\":\"failed\"}}' > {shlex.quote(job_json)}\n"
        f"  rm -f -- {shlex.quote(str(reservation))}\n"
        "  exit 1\n"
        "fi\n",
        encoding="utf-8",
    )
    script.chmod(0o700)
    try:
        fd = os.open(reservation, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        os.close(fd)
    except FileExistsError:
        return {"status": "queued", "reason": "submission already in progress"}
    try:
        result = subprocess.run(
            [
                PBS_SUBMIT_COMMAND,
                *(('-M', report_email) if report_email else ()),
                "-o", str(folder / "job.log"),
                "-j", "oe",
                str(script),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=15,
            check=False,
        )
        if result.returncode != 0 or not result.stdout.strip():
            raise RuntimeError("report could not be queued")
        job_id = result.stdout.strip().splitlines()[0]
        (folder / "job-id.txt").write_text(job_id + "\n", encoding="utf-8")
        payload = {"status": "queued", "job_id": job_id}
        # A very short job can already have written its final status before qsub returns.
        current = {}
        with contextlib.suppress(OSError, ValueError, TypeError):
            current = json.loads(job_path.read_text(encoding="utf-8"))
        if current.get("status") not in {"running", "complete", "failed"}:
            _write_json(job_path, payload)
        return payload
    except Exception:
        reservation.unlink(missing_ok=True)
        raise


def reconcile_pbs_status(folder, job):
    """Reconcile an unfinished local status with PBS without doing report work."""
    if job.get("status") not in {"queued", "running"}:
        return job
    job_id = str(job.get("job_id") or "").strip()
    if not job_id:
        return job
    status_command = Path(PBS_STATUS_COMMAND)
    if status_command.is_absolute() and (
        not status_command.is_file() or not os.access(status_command, os.X_OK)
    ):
        return job
    try:
        result = subprocess.run(
            [PBS_STATUS_COMMAND, "-f", job_id],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=5,
            check=False,
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return job
    if result.returncode == 0:
        match = re.search(r"^\s*job_state\s*=\s*([A-Z])\s*$", result.stdout, re.MULTILINE)
        if match:
            state = match.group(1)
            if state in {"R", "E", "H", "W", "T"}:
                job["status"] = "running"
            elif state == "Q":
                job["status"] = "queued"
            return job
    if (folder / "report.html").is_file():
        job["status"] = "complete"
    else:
        job.update({
            "status": "unavailable",
            "reason": "PBS no longer reports this job and no report was produced.",
        })
    _write_json(folder / "job.json", job)
    return job


def pending_page(prefix, token, folder):
    job = {}
    with contextlib.suppress(OSError, ValueError, TypeError):
        job = json.loads((folder / "job.json").read_text(encoding="utf-8"))
    job = reconcile_pbs_status(folder, job)
    status = str(job.get("status") or "queued")
    job_id = str(job.get("job_id") or "")
    if not job_id:
        with contextlib.suppress(OSError, UnicodeError):
            job_id = (folder / "job-id.txt").read_text(encoding="utf-8").strip()
    active_jobs = html.escape(OOD_ACTIVE_JOBS_URL, quote=True)
    report_url = html.escape(f"{prefix}/results/{token}/report.html", quote=True)
    if status in {"failed", "unavailable"}:
        reason = html.escape(str(job.get("reason") or "The compute job did not complete."))
        return (
            "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            "<meta name='viewport' content='width=device-width,initial-scale=1'>"
            "<title>Report unavailable</title><style>"
            "body{margin:0;background:#f4f7f9;color:#172033;font:16px system-ui}"
            ".card{max-width:34rem;margin:12vh auto;padding:2rem;background:white;"
            "border:1px solid #dbe4ee;border-radius:16px;box-shadow:0 8px 30px #0f172a12}"
            "a{color:#0369a1}</style></head><body><main class='card'>"
            f"<h1>Report unavailable</h1><p>{reason}</p>"
            f"<a href='{html.escape(prefix, quote=True)}/'>Submit a new report</a>"
            "</main></body></html>"
        )
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        "<meta http-equiv='refresh' content='10'><title>Report queued</title><style>"
        "body{margin:0;background:#f4f7f9;color:#172033;font:16px system-ui}"
        ".card{max-width:34rem;margin:12vh auto;padding:2rem;background:white;"
        "border:1px solid #dbe4ee;border-radius:16px;box-shadow:0 8px 30px #0f172a12}"
        "h1{margin-top:0;font-size:1.5rem}code{padding:.15rem .35rem;background:#eef2f7;"
        "border-radius:5px;font-size:.9em}a{color:#0369a1}.muted{color:#64748b}"
        "</style></head><body><main class='card'>"
        f"<h1>Report queued</h1><p>Status: <strong>{html.escape(status)}</strong>.</p>"
        f"<p class='muted'>{('PBS job: <code>' + html.escape(job_id) + '</code>.' if job_id else 'PBS job submitted.')}</p>"
        f"<p><a href='{active_jobs}' target='_blank' rel='noopener'>Open in OOD Active Jobs</a></p>"
        f"<p><a href='{report_url}'>Open report when ready</a></p>"
        "</main></body></html>"
    )


def application(environ, start_response):
    """Small per-user WSGI app; paths work under OOD's SCRIPT_NAME prefix."""
    method = environ.get("REQUEST_METHOD", "GET")
    path = environ.get("PATH_INFO", "/")
    prefix = environ.get("SCRIPT_NAME", "").rstrip("/")
    if method == "GET" and path in {"/health", "/__health"}:
        compute_python = EXTENDED_PYTHON or DEFAULT_COMPUTE_PYTHON
        payload = {
            "passenger_python": sys.executable,
            "passenger_python_version": sys.version,
            "project": PROJECT,
            "report_execution": "pbs_only",
            "pbs_submit_command": PBS_SUBMIT_COMMAND,
            "pbs_submit_command_exists": (
                Path(PBS_SUBMIT_COMMAND).is_file()
                if Path(PBS_SUBMIT_COMMAND).is_absolute()
                else shutil.which(PBS_SUBMIT_COMMAND) is not None
            ),
            "active_jobs_url": OOD_ACTIVE_JOBS_URL,
            "compute_python": compute_python,
            "compute_python_exists": Path(compute_python).is_file(),
            "uniprot": "deferred_to_compute_job",
        }
        return respond(
            start_response,
            "200 OK",
            json.dumps(payload, indent=2),
            "application/json; charset=utf-8",
        )
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
        and REPORT_ASSISTANT_ENABLED
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
            if (set(fields) - {"csrf", "sequence", "accession", "email"}
                    or len(fields.get("sequence", [])) > 1
                    or len(fields.get("accession", [])) > 1
                    or len(fields.get("email", [])) > 1):
                raise ValueError("submit at most one sequence and one accession")
            text = (fields.get("sequence") or [""])[0]
            accession = (fields.get("accession") or [""])[0].strip().upper()
            email = (fields.get("email") or [""])[0].strip()
            if email and not EMAIL_RE.fullmatch(email):
                raise ValueError("invalid email address")
            if accession and not re.fullmatch(
                r"(?:[A-Z0-9]{6}|[A-Z0-9]{10})(?:-[0-9]{1,2})?", accession
            ):
                raise ValueError("invalid UniProt accession")
            if len(text.encode("utf-8")) > 100_000:
                raise ValueError("sequence input exceeds 100 KB")
            if not text.strip() and not accession:
                raise ValueError("provide either a protein sequence or a UniProt accession")
            if text.strip():
                validate_sequence(text)
            # Passenger validates only; the compute-node parser rechecks the same input.
            token = secrets.token_hex(16)
            ROOT.mkdir(mode=0o700, parents=True, exist_ok=True)
            folder = ROOT / token
            folder.mkdir(mode=0o700)
            try:
                (folder / "input.fasta").write_text(text, encoding="utf-8")
                request = {"accession": accession}
                if email:
                    request["email"] = email
                _write_json(folder / "request.json", request)
                submit_report_job(folder)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                shutil.rmtree(folder)
                raise
            return respond(start_response, "303 See Other", b"", headers=[
                ("Location", prefix + "/results/" + token + "/report.html"),
            ])
        except (ValueError, UnicodeError):
            return respond(start_response, "400 Bad Request",
                           page(prefix, stored if TOKEN.fullmatch(stored) else "",
                                "<p role='alert'>Invalid sequence submission. Use one FASTA "
                                "record or bare sequence, up to 5,000 residues.</p>"))
        except (OSError, RuntimeError, subprocess.TimeoutExpired):
            return respond(
                start_response,
                "503 Service Unavailable",
                "<p>The report could not be queued. Check PBS and compute-environment configuration.</p>",
            )
    if method == "GET" and path.startswith("/results/"):
        parts = path.split("/")
        if len(parts) == 4 and TOKEN.fullmatch(parts[2]) and parts[3] in FILES:
            folder = ROOT / parts[2]
            result = folder / parts[3]
            if parts[3] == "report.html" and not result.is_file() and (folder / "job.json").is_file():
                return respond(start_response, "200 OK", pending_page(prefix, parts[2], folder))
            if result.is_file():
                body = result.read_bytes()
                csrf = secrets.token_hex(16)
                if parts[3] == "report.html" and REPORT_ASSISTANT_ENABLED:
                    link = (
                        b'<div class="report-shell"><div class="note"><h2>Ask about this report</h2>'
                        b'<p><a href="chat">Open the report-only assistant</a>. It can explain '
                        b'this report but cannot perform fresh database retrieval.</p></div></div>'
                    )
                    body = body.replace(b"</body>", link + b"</body>")
                headers = []
                if parts[3] == "report.html":
                    headers.append(("Set-Cookie", "protein_report_csrf=%s; Path=%s; Secure; HttpOnly; SameSite=Strict" % (csrf, prefix + "/")))
                return respond(start_response, "200 OK", body, FILES[parts[3]], headers=headers)
    return respond(start_response, "404 Not Found", "Not found")
