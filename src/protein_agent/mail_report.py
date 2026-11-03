"""Send completed report files through the site's local sendmail relay."""

from __future__ import annotations

import argparse
import mimetypes
import shutil
import subprocess
from email.message import EmailMessage
from pathlib import Path


def send_report(report_dir: Path, recipient: str) -> None:
    sender = shutil.which("sendmail") or "/usr/sbin/sendmail"
    if not Path(sender).is_file():
        raise RuntimeError("sendmail is not available on the compute node")
    message = EmailMessage()
    message["To"] = recipient
    message["Subject"] = "Protein report ready"
    message.set_content(
        "Your protein report is ready. The HTML and JSON reports are attached.\n"
    )
    for name in ("report.html", "report.json"):
        path = report_dir / name
        if not path.is_file():
            raise FileNotFoundError(path)
        data = path.read_bytes()
        mime, _ = mimetypes.guess_type(path.name)
        maintype, subtype = (mime or "application/octet-stream").split("/", 1)
        message.add_attachment(data, maintype=maintype, subtype=subtype, filename=name)
    subprocess.run(
        [sender, "-t", "-oi"], input=message.as_bytes(), check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=30,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--recipient", required=True)
    args = parser.parse_args()
    send_report(args.report_dir, args.recipient)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
