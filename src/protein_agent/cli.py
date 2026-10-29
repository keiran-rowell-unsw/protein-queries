"""CLI: ask a question, or inspect the gateway.

    uv run python -m protein_agent --list-models
    uv run python -m protein_agent "function of human TP53"
    uv run python -m protein_agent --dry-run --explain "any question"
"""

from __future__ import annotations

import json

import typer
from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from rich.table import Table

from .agent import ProteinAgent
from .config import Settings
from .gateway import GatewayError, ResTechGateway
from .proof import AnswerProofError
from .retrieval import RetrievalError

app = typer.Typer(add_completion=False, help="Retrieval-grounded protein data assistant.")
console = Console()


@app.command()
def main(
    question: str = typer.Argument("", help="Natural-language question about a protein."),
    model: str = typer.Option("", "--model", "-m", help="Gateway model id (default from env)."),
    min_calls: int = typer.Option(1, "--min-calls", help="Minimum tool calls per answer."),
    on_unproven: str = typer.Option(
        "block", "--on-unproven", help="block | warn (warn is for review only)."
    ),
    explain: bool = typer.Option(False, "--explain", help="Show the service/tool plan."),
    fasta: str = typer.Option(
        "", "--fasta", "-f", help="FASTA file; computes properties from the first record."
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Stub all network calls; exercises the loop logic only."
    ),
    list_models: bool = typer.Option(False, "--list-models", help="Print the live model roster."),
) -> None:
    settings = Settings.from_env(require_key=not (dry_run or list_models))
    if model:
        settings.model = model

    if dry_run:
        settings.require_tool_calls = False
        console.print("[yellow]DRY RUN:[/yellow] gateway and UniProt are stubbed; "
                      "this proves wiring, not retrieval.")

    if list_models:
        try:
            with ResTechGateway(settings) as gateway:
                ids = gateway.list_models()
        except GatewayError as exc:
            console.print(f"[red]{exc}[/red]")
            raise typer.Exit(2) from None
        table = Table(title=f"Live models @ {settings.base_url}")
        table.add_column("id")
        table.add_column("chat model?")
        for item in ids:
            table.add_row(item, "yes" if item in {settings.model} else "")
        console.print(table)
        return

    if fasta:
        from .properties import SequenceError, compute_properties, read_first_record

        try:
            header, sequence, warnings = read_first_record(fasta)
        except SequenceError as exc:
            console.print(f"[red]Unusable sequence:[/red] {settings.redact(str(exc))}")
            console.print("Nothing was computed; fix the input rather than guessing it.")
            raise typer.Exit(6) from None
        report = compute_properties(sequence, header=header).as_dict()
        report["warnings"] = warnings + report["warnings"]
        for warning in report["warnings"]:
            console.print(f"[yellow]input:[/yellow] {warning}")
        console.print(f"[bold]{header}[/bold] · {report['length_residues']} residues")
        table = Table(box=None, show_header=False)
        table.add_column(justify="right", style="dim")
        table.add_column(justify="left")
        for key in (
            "molecular_weight_da",
            "isoelectric_point",
            "net_charge_pH_7",
            "gravy_kyte_doolittle",
            "aromaticity_iron_dos_santos",
            "instability_index_guruprasad",
        ):
            table.add_row(key.replace("_", " "), str(report[key]))
        table.add_row(
            "molar extinction (red/ox)",
            " / ".join(str(v) for v in report["molar_extinction_coefficient_M-1cm-1"].values()),
        )
        table.add_row(
            "hydropathy min/max",
            f"{report['hydropathy_extremes']['min']} / {report['hydropathy_extremes']['max']}",
        )
        console.print(table)
        console.print("[dim]computed from composition — no structure or expression implied[/dim]")
        console.print_json(json.dumps(report))
        return

    if not question.strip():
        console.print("[red]Give a question, a --fasta file, or use --list-models.[/red]")
        raise typer.Exit(1) from None

    agent = ProteinAgent(settings)
    plan = None
    if explain:
        from .router import build_plan

        plan = build_plan(question)
        console.print(
            Panel(
                f"services: {', '.join(plan.services)}\n"
                f"tools advertised: {', '.join(plan.tool_names)}\n"
                f"why: {plan.reason}",
                title="plan (model still chooses the calls)",
            )
        )

    try:
        result = agent.run(question, min_calls=min_calls, on_unproven=on_unproven)
    except AnswerProofError as exc:
        console.print(f"[red]Not answered:[/red] {exc}")
        console.print("No answer was produced because it could not be tied to tool output.")
        raise typer.Exit(3) from None
    except GatewayError as exc:
        console.print(f"[red]Inference unavailable:[/red] {exc}")
        console.print("Report the outage; do not substitute a remembered answer.")
        raise typer.Exit(4) from None
    except RetrievalError as exc:
        console.print(f"[red]Retrieval failed:[/red] {exc}")
        raise typer.Exit(5) from None

    calls = ", ".join(
        f"{c['name']}({', '.join(f'{k}={v}' for k, v in c['arguments'].items())})"
        for c in result.evidence.calls
    ) or "none"
    prompt_tokens = result.usage.get("prompt_tokens")
    console.print(
        f"[dim]model={result.model} rounds={result.rounds} "
        f"tool_calls={result.evidence.call_count} prompt_tokens={prompt_tokens}[/dim]"
    )
    console.print(f"[green]backed by:[/green] {calls}")
    console.print(Markdown(result.answer))


if __name__ == "__main__":
    app()
