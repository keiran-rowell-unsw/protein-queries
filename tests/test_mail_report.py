from types import SimpleNamespace

from protein_agent import mail_report


def test_send_report_builds_attachments_and_bounds_sendmail(monkeypatch, tmp_path):
    (tmp_path / "report.html").write_text("<html></html>")
    (tmp_path / "report.json").write_text("{}")
    sender = tmp_path / "sendmail"
    sender.write_text("")
    sender.chmod(0o700)
    seen = {}

    def fake_run(args, **kwargs):
        seen["args"] = args
        seen["kwargs"] = kwargs
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(mail_report.shutil, "which", lambda _name: str(sender))
    monkeypatch.setattr(mail_report.subprocess, "run", fake_run)
    mail_report.send_report(tmp_path, "user@example.edu")
    assert seen["args"] == [str(sender), "-t", "-oi"]
    assert seen["kwargs"]["timeout"] == 30
    assert b"report.html" in seen["kwargs"]["input"]
    assert b"report.json" in seen["kwargs"]["input"]
