# protein-queries

Retrieval-grounded protein tools for UniProt, PDBe, sequence analysis, and
structure-neighbour searches. Inference uses the UNSW ResTech Local Inference
gateway; it does not use a vendor cloud API.

## Choose a workflow

| Workflow | Use it for | Runs where |
| --- | --- | --- |
| Agent CLI | Questions about a named protein | Local client plus UniProt/PDBe services |
| Local report | Developing or testing a sequence report | The current machine |
| Open OnDemand report | Full reports without login-node computation | PBS compute node |

## Install

Python 3.11 or newer is required. Install the pinned environment with:

```bash
uv sync --extra dev
```

The agent CLI requires an invite-only ResTech key:

```bash
export RESTECH_API_KEY=sk-…       # never commit or print this value
```

## Agent CLI

Ask a retrieval-grounded question or inspect available models:

```bash
uv run python -m protein_agent --model qwen3.8-flash-next \
  "function of human TP53"
uv run python -m protein_agent --list-models
uv run python -m protein_agent --dry-run "…"
```

The first model request must make a retrieval tool call. Later requests may use
`auto` tool choice so the model can finish after retrieval. The proof layer then
checks the answer against recorded tool results and blocks unsupported selected
identifiers and claim patterns. This is a safeguard, not exhaustive semantic
verification; model output is never evidence by itself.

Protein identity, UniProt annotations, and experimental structures must come from
retrieved UniProt or PDBe records. Sequence-derived properties are calculated
locally. Disorder is a prediction, and Foldseek results are structural-neighbour
candidates, not identity claims.

## Local sequence report

Use this mode for development or a local test only. It performs calculations and
external retrieval on the current machine, so do not use it on a login node.

```bash
uv run python -m protein_agent.web --port 8765 --output "$HOME/protein-reports"
# Open http://127.0.0.1:8765/ on the same machine

# Or write a report without starting the web server:
uv run python -m protein_agent.report --input input.fasta --output results
uv run python -m protein_agent.report --output results --accession YOUR_ACCESSION
```

Reports contain HTML, JSON, and available plots. They may include composition,
charge, extinction estimates, residue groups, hydropathy, disorder prediction,
UniProt annotations, PDBe structure mappings, and optional Foldseek candidates.

Input and identity rules:

- Accept one FASTA record or a bare sequence, up to 5,000 residues.
- An accession may be supplied without a sequence and takes precedence when both
  are supplied.
- A pasted sequence is matched to UniProt by exact sequence, not by similarity.
- Ambiguous matches, missing records, network failures, and unavailable optional
  services are reported; they are not filled in from model memory.
- Metapredict reports disorder, not solubility. The report does not run BLAST,
  AlphaFold, structure prediction, SignalP, or other tools not listed above.

## Open OnDemand / PBS report app

The Passenger front end validates input, stores it in the user’s private report
directory, submits one PBS job, and serves the status and completed report. It does
not perform scientific calculations, retrieval, inference, or structure processing.
The compute job runs UniProt verification, metapredict, PDBe retrieval, and the
optional structure-backed Foldseek path.

The app is deployed from the `passenger-fully-queued` branch:

```bash
cd ~/ondemand/dev
git clone --branch passenger-fully-queued \
  git@github.com:keiran-rowell-unsw/protein-queries.git protein-queries-queued
```

Deployment prerequisites:

- Passenger/PUN must be allowed to invoke `qsub`.
- Passenger needs Python 3.11 or newer, but not the scientific dependencies.
- The project clone and report directory must be visible at the same paths on the
  PUN and compute nodes.
- Set `PROTEIN_REPORT_EXTENDED_PYTHON` to a maintained compute Python, or create
  `<project>/.compute-venv` once on an allocated compute node with
  `bin/setup-compute-env`.
- Set `PBS_SUBMIT_COMMAND` only when `qsub` is not discoverable.
- Establish a retention policy; reports are not automatically deleted.

### Email delivery

The form has an optional email field. When supplied, the completed `report.html`
and `report.json` are sent as attachments. `PROTEIN_REPORT_EMAIL` may provide a
default recipient when the form field is blank. Both paths require an approved
`sendmail` relay on the compute node. Attachment delivery is best-effort and does
not change a successful report into a failed report.

PBS also requests its normal start/end/abort notifications (`#PBS -m abe`). An
email address supplied by a user is passed to PBS as the notification recipient;
configure site policy before allowing external addresses. The email address is
kept in the private request directory and is not put in the report URL.

No dependency installation occurs during a web request or report job. Restart the
PUN after changing code or environment configuration. The separate
`protein_agent.web` prototype has no authentication and is not the OOD deployment.

## Report-backed Foldseek path

When enabled, Foldseek runs only after a verified UniProt accession has a usable
PDBe-mapped experimental structure. The report downloads the mapped structure,
extracts the mapped chain, and submits its mmCIF coordinates to Foldseek in
`pdb100` 3Di/AA mode. Raw FASTA is not submitted as a 3D structure.

Up to five candidate hits, with bounded scores and provenance, are shown in the
report and saved in `report.json`. Missing mappings, inaccessible structures,
service failures, and busy searches are reported as skipped, unavailable, or
timed out—not as evidence that no structural neighbours exist.

## Configuration

| Variable | Purpose |
| --- | --- |
| `RESTECH_API_KEY` | Required inference-gateway bearer key; never logged |
| `RESTECH_BASE_URL` | Inference gateway; defaults to the ResTech `/v1` URL |
| `RESTECH_MODEL` | Default chat model; pin a model for reproducibility |
| `RESTECH_TIMEOUT` / `RESTECH_TEMPERATURE` | Optional gateway settings |
| `MAX_TOOL_ROUNDS` | Agent-loop bound; defaults to 6 |
| `UNIPROT_BASE_URL` | UniProt API base; defaults to `https://rest.uniprot.org` |
| `PROTEIN_REPORT_EXTENDED_PYTHON` | Compute-node Python for queued reports |
| `PROTEIN_REPORT_EMAIL` | Optional default report recipient |
| `PBS_SUBMIT_COMMAND` | Optional absolute or named `qsub` command |

See `.env.example`, `pyproject.toml`, and `uv.lock` for the complete environment
and dependency definitions.

## Development checks

```bash
uv run pytest -q
uv run ruff check .
```

Tests cover retrieval grounding, report rendering, OOD validation and submission,
PDBe mapping, Foldseek input selection, and failure handling. Live service access
is not required for the test suite.

## Repository guide

- `src/protein_agent/agent.py` — bounded model/tool loop
- `src/protein_agent/retrieval.py` — UniProt REST retrieval and provenance
- `src/protein_agent/pdbe_mcp.py` — PDBe MCP discovery and calls
- `src/protein_agent/report.py` — sequence report generation
- `ood_app.py` — PBS-only Passenger front end
- `bin/setup-compute-env` — compute-node dependency setup
- `tests/` — automated checks

## Scope and limitations

This project reports database annotations, calculations, predictions, and
candidate matches. It is not a clinical, diagnostic, therapeutic, or experimental
protocol service. Missing annotation means “not annotated in this record,” not
that the biological property is absent. Retain the accession, PDB identifier,
retrieval date, and provenance shown in a report when using its results.
