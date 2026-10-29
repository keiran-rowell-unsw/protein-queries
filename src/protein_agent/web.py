"""Local-only sequence report web interface; reverse-proxy behind authenticated OOD.

This is not an authentication server. Never expose it on a public interface; only bind
127.0.0.1 and put an authenticated, site-approved OOD proxy in front of it.
"""

# Keep HTML template and security header intact.
# ruff: noqa: E501
from __future__ import annotations

import argparse
import html
import os
import secrets
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

from .properties import SequenceError
from .report import MAX_INPUT_BYTES, generate_report

PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Protein report</title>
<style>body{{font:1rem/1.5 system-ui;max-width:720px;margin:2rem auto;padding:1rem}}
textarea{{width:100%;box-sizing:border-box}}button{{padding:.6rem 1rem}}</style></head><body>
<h1>Protein sequence report</h1><p>Runs locally on this node. Paste one FASTA record or bare sequence
(up to 10,000 standard amino acids). Exact matches are verified through UniProt; a verified
accession is then used to retrieve ranked experimental structure mappings from PDBe.</p>
{message}<form method="post" action="/run"><label for="sequence">Protein sequence</label><br>
<textarea id="sequence" name="sequence" rows="12" required maxlength="100000"></textarea><br>
<label for="accession">UniProt accession (optional; exact sequence match required)</label><br>
<input id="accession" name="accession" maxlength="13" pattern="[A-Za-z0-9]{{6,10}}(-[0-9]{{1,2}})?" placeholder="e.g. P04637 or P04637-1" inputmode="text" autocomplete="off" spellcheck="false"><br>
<button type="submit">Generate report</button></form></body></html>"""


def make_handler(output_root: Path) -> type[BaseHTTPRequestHandler]:
    class Handler(BaseHTTPRequestHandler):
        def _send(self, status: int, body: bytes, mime: str = "text/html; charset=utf-8") -> None:
            self.send_response(status)
            self.send_header("Content-Type", mime)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy", "default-src 'none'; img-src 'self'; style-src 'unsafe-inline'")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self) -> None:
            path = urlsplit(self.path).path
            if path == "/":
                self._send(200, PAGE.format(message="").encode())
                return
            parts = path.strip("/").split("/")
            if len(parts) != 3 or parts[0] != "results" or len(parts[1]) != 32 or not all(
                c in "0123456789abcdef" for c in parts[1]
            ) or parts[2] not in {"report.html", "report.json", "hydropathy.png", "disorder.png"}:
                self.send_error(404)
                return
            file = output_root / parts[1] / parts[2]
            if not file.is_file():
                self.send_error(404)
                return
            mime = {"report.html": "text/html; charset=utf-8", "report.json": "application/json", "hydropathy.png": "image/png", "disorder.png": "image/png"}[parts[2]]
            self._send(200, file.read_bytes(), mime)

        def do_POST(self) -> None:
            if urlsplit(self.path).path != "/run" or self.headers.get_content_type() != "application/x-www-form-urlencoded":
                self.send_error(400)
                return
            try:
                size = int(self.headers.get("Content-Length", "0"))
                if size < 1 or size > MAX_INPUT_BYTES * 3:
                    raise SequenceError("submission is empty or exceeds the 100 KB input limit")
                from urllib.parse import parse_qs

                raw = self.rfile.read(size).decode("utf-8", errors="strict")
                fields = parse_qs(raw, keep_blank_values=True)
                if (set(fields) - {"sequence", "accession"} or "sequence" not in fields
                        or len(fields["sequence"]) != 1 or len(fields.get("accession", [])) > 1):
                    raise SequenceError("submit exactly one sequence")
                text = fields["sequence"][0]
                accession = (fields.get("accession") or [""])[0].strip().upper()
                token = secrets.token_hex(16)
                folder = output_root / token
                folder.mkdir(mode=0o700, parents=True)
                source = folder / "input.fasta"
                source.write_text(text, encoding="utf-8")
                try:
                    generate_report(source, folder, accession=accession)
                finally:
                    source.unlink(missing_ok=True)
            except (SequenceError, UnicodeError, OSError, ValueError) as exc:
                message = f"<p role='alert'>Report not generated: {html.escape(str(exc))}</p>"
                self._send(400, PAGE.format(message=message).encode())
                return
            self.send_response(303)
            self.send_header("Location", f"/results/{token}/report.html")
            self.send_header("Content-Length", "0")
            self.send_header("Cache-Control", "no-store")
            self.end_headers()

    return Handler


def main() -> None:
    parser = argparse.ArgumentParser(description="Local sequence-report web server (no built-in auth)")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--output", type=Path, default=Path.home() / "protein-reports")
    args = parser.parse_args()
    # Default process umask keeps submissions/results private to the service account.
    os.umask(0o077)
    with ThreadingHTTPServer(("127.0.0.1", args.port), make_handler(args.output)) as server:
        print(f"Serving locally on http://127.0.0.1:{args.port}/", flush=True)
        server.serve_forever()


if __name__ == "__main__":
    main()
