# protein-queries

Protein reports using UniProt, PDBe, local sequence analysis, and optional
structure-neighbour searches. The Open OnDemand report app is deterministic and
requires no inference service, model, or API key.

## Choose a workflow

| Workflow | Use it for | Runs where |
| --- | --- | --- |
| Local report | Developing or testing a sequence report | The current machine |
| Open OnDemand report | Full reports without login-node computation | PBS compute node |

## Install

Python 3.11 or newer is required. Install the pinned environment with:

```bash
uv sync --extra dev
```

No API credentials are required for report generation. Protein identity, UniProt
annotations, and experimental structures come from retrieved UniProt or PDBe
records. Sequence-derived properties are calculated locally. Disorder is a
prediction, and Foldseek results are structural-neighbour candidates rather than
identity claims.

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
directory, submits one PBS job, and serves the status and completed report. Apart
from the bounded, offline `/protein-prop-summary` calculations described below, it
does not perform scientific calculations, retrieval, inference, or structure
processing.
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
- Set `PROTEIN_REPORT_EXTENDED_PYTHON` to a maintained compute Python. For a central
  `/sys` installation, the administrator should provide this shared environment and
  configure the path; users do not need to run the setup script. Otherwise, create
  `<project>/.compute-venv` once on an allocated compute node with
  `bin/setup-compute-env`. The script provisions and validates the environment; it is
  not run during Passenger requests or on login nodes.
- Set `PBS_SUBMIT_COMMAND` only when `qsub` is not discoverable.
- Configure a 90-day retention cleanup for report directories; automatic deletion is
  an administrator deployment task and is not performed during web requests.

### Immediate preliminary summary (browser retrieval + bounded PUN calculations)

The submission page offers two paths:

- **Immediate preliminary summary.** The user’s browser makes one request to the
  public UniProt API for a typed accession (cross-origin requests are permitted;
  the JSON already contains the canonical sequence, so retrieval is a single
  fetch). The sequence is then posted to the app’s stateless
  `/protein-prop-summary` endpoint, which runs bounded, offline, deterministic
  library calculations (Biopython ProtParam plus the pepstats-equivalent metrics;
  no retrieval, no network, no writes, sub-second for the 5,000-residue limit)
  and returns composition-derived properties, including the CRC64 checksum.
  The response also carries a downsampled Kyte-Doolittle hydropathy profile and
  the charge-versus-pH series, which the page renders as inline SVG curves; the
  numbers come from the same library calculations as the queued report, and the
  full-resolution plots remain part of the PBS output. The summary also renders a
  **preliminary position map** as inline SVG: filled tracks are features annotated
  in the retrieved UniProt entry (held browser-side), hollow tracks are
  library-predicted N-glycan sequons and deamidation motifs computed by the PUN,
  and the two provenances are never mixed on one marker. The queued report renders
  the authoritative, full-resolution map.
  For a pasted sequence the browser then runs one bounded UniParc checksum
  search (`size=5`, total capped via the CORS-exposed `X-Total-Results` header)
  and, when exactly one UniProtKB accession matches after residue-level
  verification, one further entry fetch to show the resolved identity. Zero or
  several matches are reported as unresolved/ambiguous with candidates — never
  guessed. This exact-sequence match mirrors the queued report’s identity logic;
  it is not a similarity search.
- **Queue full report.** Unchanged: one PBS job performs UniProt verification,
  metapredict, PDBe retrieval, and the optional Foldseek path, and its report is
  authoritative.

The Content Security Policy permits exactly two remote destinations for scripts
and connections: the app itself and `https://rest.uniprot.org`. Fetched records
are rendered with text-safe DOM APIs only. The immediate summary requires the
Passenger interpreter to provide Biopython; if it is unavailable the endpoint
reports that the summary is unavailable and queueing still works. The PUN
endpoint never contacts UniProt or any other external service.

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
PUN after changing code or environment configuration. This is a normal queued PBS
application, not a Batch Connect app: it does not request or hold an interactive
session. Each submission is one short-lived PBS job, and the Passenger page serves
its static result after completion. Batch Connect is appropriate for interactive
notebooks, shells, or applications that need a persistent allocation. The separate
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
