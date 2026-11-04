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
directory, submits one PBS job, and serves the status and completed report. It
performs no scientific calculations, retrieval, inference, structure processing,
or plotting.
The compute job runs UniProt verification, metapredict, PDBe retrieval, and the
optional structure-backed Foldseek path.

### System (/sys) deployment runbook

Users never clone anything: they open the app from the OOD dashboard and submit
the form. All git and environment work is administrator-side.

1. **Immutable release.** Copy a validated tree to an administrator-owned directory
   outside any group-writable scratch (for example
   `/srv/ood-apps/protein-queries/releases/<version>/`), with directories `0755` and
   files `0644` and no group write or setgid bit. Strip `.env`, `.compute-venv`,
   `tmp/`, and `__pycache__/` from the release. Every user's PUN and PBS job
   executes this code, so nothing here may be user-writable.
2. **Register the app.** Symlink the release into the OOD system apps location
   (for example `/var/www/ood/apps/sys/protein-queries`). `manifest.yml` already
   places the tile under *Interactive Apps → Biology*.
3. **One central compute environment.** Run `bin/setup-compute-env` once inside an
   administrator PBS allocation, with `UV_PROJECT_ENVIRONMENT` pointed at a shared
   path visible to the compute nodes (for example
   `/srv/ood-apps/protein-queries/compute-env/`). The script provisions and
   validates the environment (metapredict, torch, protein_agent imports).
4. **Configure via OOD's per-app environment** (`/etc/ood/config/apps/<app>/env` on
   standard installs). Recognised settings:
   - `PROTEIN_REPORT_EXTENDED_PYTHON` — the central compute Python above; setting it
     removes any per-user setup step.
   - `PROTEIN_REPORT_WALLTIME` / `PROTEIN_REPORT_MEMORY` / `PROTEIN_REPORT_QUEUE` —
     scheduler requests (defaults `00:10:00`, `2gb`, no queue). Values are validated;
     a 5,000-residue report normally needs a couple of minutes on one core.
   - `PROTEIN_REPORT_EMAIL_DOMAIN` — domain used to build the report recipient from
     the signed-in user (default `unsw.edu.au`).
   - `PBS_SUBMIT_COMMAND` / `PBS_STATUS_COMMAND` — only when `qsub`/`qstat` are not
     discoverable on the web host.
5. **Canary.** Submit one report from the web UI and confirm it runs on a compute
   node, renders its figures, and appears under OOD Active Jobs.
6. **Retention.** Report folders under `~/ondemand/data/protein-queries/reports`
   are persistent user artifacts (a few hundred KB each). Because this is a custom
   Passenger app, not a Batch Connect app, OOD does not manage their lifecycle;
   install a site cron line (for example `find … -mtime +90 -delete`) or accept
   quota growth for a small user base.
7. **Upgrades.** Build a new release directory, validate, then flip a `current`
   symlink; running jobs finish against the old release and restarted PUNs pick up
   the new one.

For development on the current machine, a per-user clone under `~/ondemand/dev`
works exactly as before; the runbook above replaces it for the shared service.

Deployment prerequisites:

- Passenger/PUN must be allowed to invoke `qsub`.
- Passenger needs Python 3.11 or newer, but not the scientific dependencies.
- The release directory and central compute environment must be visible at the
  same paths on the PUN and compute nodes.

### Submission page and identity preview

The form takes an optional UniProt accession and a sequence (one FASTA record or
bare sequence, up to 5,000 residues; at least one of accession or sequence is
required). Report delivery by email is automatic and derived from the sign-in
(see *Email delivery* below); there is nothing else to fill in.

Typing an accession shows an **identity preview**: the browser makes one request
to the public UniProt API (cross-origin requests are permitted) and displays the
entry’s name, organism, length, and reviewed status, so the submitter can confirm
the right protein before queueing. The preview is preliminary and text-only; it
asserts nothing the queued job will not verify.

There is one submission path: **Queue full report** submits one PBS job that
verifies identity (exact-sequence match for a pasted sequence; accession
precedence when both are given), runs disorder prediction, retrieves UniProt
annotations and PDBe structures, optionally searches structural neighbours, and
renders every figure with matplotlib on the compute node. The report it writes is
the authoritative result; the front end performs no calculations.

The Content Security Policy permits exactly two remote destinations for scripts
and connections: the app itself and `https://rest.uniprot.org`. Fetched records
are rendered with text-safe DOM APIs only, and the page draws no charts. The PUN
never contacts any other external service and never imports matplotlib or any
scientific library.

### Email delivery

The form has no email field. The recipient is derived from the authenticated
portal user (`REMOTE_USER`): its local part (a realm suffix is dropped) is
completed with `PROTEIN_REPORT_EMAIL_DOMAIN` (default `unsw.edu.au`), so the
cluster relay can only ever mail the submitter's own address. When no usable
identity is available the report is still queued and served through the web UI.
`PROTEIN_REPORT_EMAIL` remains an administrator override that routes every report
to one mailbox. Delivery uses the compute node's approved `sendmail` relay, is
best-effort, and does not change a successful report into a failed report.

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
