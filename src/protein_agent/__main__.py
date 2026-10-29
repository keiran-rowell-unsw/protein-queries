"""Allow `python -m protein_agent` to run the CLI."""

from __future__ import annotations

from .cli import app

if __name__ == "__main__":
    app()
