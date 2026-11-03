# protein-queries

Retrieval-grounded protein data assistant. Model generation runs on the UNSW ResTech
**Local Inference** gateway (OpenAI-compatible, on-prem); factual content comes only from
the UniProt REST API and the PDBe-KB MCP servers.

| Path | Purpose |
| --- | --- |
| `AGENTS.md` | Shared agent rules: inference backend, tool-calling mandate, answer format, scope, failure handling |
| `uniprot/AGENTS.md` | UniProt identity, sequence, and annotation retrieval |
| `pdbe-kb/AGENTS.md` | PDBe experimental structure retrieval (MCP API + Search servers) |
| `pyproject.toml` | uv project: Python >=3.11, `httpx`, `pydantic`, `typer`, `rich` |
| `src/protein_agent/` | Gateway client, tool schemas, agent loop, proof layer |
| `tests/test_grounding.py` | Fails if an unbacked answer can reach the user |
| `.env.example` | Template; copy to `.env` (git-ignored) |

## Local sequence-report interface (no PBS queue)

A researcher can paste one protein sequence into a local web form and download a
`report.html`, `report.json`, and `hydropathy.png`. It uses the installed
Biopython/pepstats-equivalent calculations and hydropathy plotter **on this node**, plus
metapredict V3 sequence-based disorder predictions and `disorder.png` when available
(for sequences of 25–2,000 residues). Shorter/longer sequences or predictor failures are
explicitly reported as skipped/unavailable, never filled in. The report includes
composition, charge-versus-pH, extinction estimates, residue groups, a cautious
E. coli inclusion-body composition flag and an optional chromophore motif candidate scan.
Advanced values and method details remain in the downloadable JSON; the HTML prioritises
results and caveats. With **no accession**, the app computes the sequence's UniProt CRC64
checksum (Biopython), searches the UniParc exact-sequence index, verifies that the entire
returned sequence matches, and checks up to four linked UniProtKB canonical entries.
If exactly one matches, the verified accession and its annotations are shown. The report
then uses the dynamically discovered PDBe API MCP ranked-structure tool to show up to ten
experimental PDB mappings with method, resolution when returned, chain and UniProt coverage.
If matches are ambiguous, it lists candidates and asks the researcher to enter a specific
accession; it never selects the first hit. A **UniProt accession may be supplied without a
sequence** and takes precedence when both are present. The report then retrieves and analyses
the accession's canonical sequence. If a submitted sequence differs, that mismatch is stated
and the submitted sequence is not analysed. UniProt sections cover experimental-structure
cross-references, localisation, glycosylation, PTMs, subunit organisation, domains,
functional sites, catalytic activity, PPI records, variants and annotation evidence.
A missing hit, network failure, or ambiguous sequence match is stated, not guessed past.
Only a checksum goes into the search URL; the full sequence is verified against its
returned record. This is an **exact** search, not a BLAST or similarity search: variants,
tags, fragments and isoforms without a canonical match will not be identified.
It does not run BLAST,
AlphaFold confidence analysis, structure prediction, Protein-sol,
SoDoPE, SignalP/Razor, AI summaries, ChimeraX surface MLP or queued programs.
Metapredict predicts disorder, **not solubility**. Identity is unresolved without a
matching accession; experimental structures are reported only from PDBe mappings to a
verified accession. Multiple FASTA records, gaps, ambiguous amino acids, and sequences longer
than 5,000 residues are rejected rather than silently changed.

For a local test only:

```bash
.venv/bin/python -m protein_agent.web --port 8765 --output "$HOME/protein-reports"
# On the same node: visit http://127.0.0.1:8765/
# Or generate a report without a web service:
.venv/bin/python -m protein_agent.report --input input.fasta --output results
# A UniProt accession works without sequence input and takes precedence if both are supplied:
.venv/bin/python -m protein_agent.report --output results --accession YOUR_ACCESSION
```

### Open OnDemand Passenger app: fully queued

This branch provides a **PBS-only** Passenger app. The per-user web process validates a
sequence or UniProt accession, writes the request to the user's private
`~/ondemand/data/protein-queries/reports/<random-id>/` directory, submits one PBS job, and
serves a status page. It performs no Biopython calculations, retrieval, inference, or
structure processing on the OOD host. The compute job produces the complete report,
including UniProt verification, metapredict, PDBe retrieval and optional Foldseek search.
The PBS job ID is retained as `job-id.txt` for scheduler accounting.

```bash
cd ~/ondemand/dev
git clone --branch passenger-fully-queued \
  git@github.com:keiran-rowell-unsw/protein-queries.git protein-queries-queued
```

**Deployment prerequisites:**

- The PUN must support Python Passenger and be allowed to invoke `qsub`.
- `bin/python` must find Python 3.11 or newer; no scientific Python environment is needed
  by Passenger itself.
- The project clone and report directory must be visible at identical paths on the PUN and
  compute nodes.
- Prefer a centrally maintained full compute environment and set
  `PROTEIN_REPORT_EXTENDED_PYTHON` to its absolute Python path. Alternatively, an
  administrator can create `<project>/.compute-venv` once with `bin/setup-compute-env`.
- Set `PBS_SUBMIT_COMMAND` only if `qsub` is not available under that name.
- Establish a report-retention policy; results are not automatically deleted.

No dependency installation occurs during a web request or report job. PBS sends the
requestor a normal start/end/abort notification (`#PBS -m abe`). To attach the completed
`report.html` and `report.json`, set `PROTEIN_REPORT_EMAIL` explicitly in the Passenger
environment; the compute node must provide `/usr/sbin/sendmail` (or `sendmail` on PATH)
and an approved relay. Attachment delivery is best-effort and never changes a successful
report into a failed one. Restart the PUN after changing branches or environment configuration. The separate
`protein_agent.web` prototype has no authentication and is not the OOD deployment.

## Hard rule: no memory-only answers

Every user-visible statement about a specific protein must be backed by a tool call that
ran in that same request. `tool_choice: "required"` enforces it at the API level, and the
client refuses an answer whose facts do not map onto recorded tool results. A model that
"knows" TP53 is not evidence; `GET /uniprotkb/P04637.json` is.

## How the mandate is enforced

Three layers, so a memory-only answer cannot slip through even if one is bypassed:

1. **Prompt** (`prompts.py`) states the rule to the model.
2. **API**: every request in the loop carries `tool_choice: "required"`, so replying with
   prose and no tool call is not an available action.
3. **Proof** (`proof.py`) re-checks the finished answer against recorded tool output and
   raises `AnswerProofError` if any accession, PDB ID, GO term, residue position, or length
   is unbacked, or if no tool call ran at all.

`REQUIRE_TOOL_CALLS=0` exists for tests only and must never be set in a deployment.

## Setup

```bash
uv sync --extra dev          # creates .venv, installs pinned deps + Python >=3.11
export RESTECH_API_KEY=sk-…  # invite-only key from ResTech; never commit it
uv run python -m protein_agent --model qwen3.8-flash-next "function of human TP53"
```

`uv run` uses its own isolated interpreter, so this node's system Python 3.6 is never used.

## Environment

| Variable | Meaning |
| --- | --- |
| `RESTECH_API_KEY` | Bearer key for the inference gateway. Required. Never logged or echoed. |
| `RESTECH_BASE_URL` | Defaults to `https://inference.restech.unsw.edu.au/v1` |
| `RESTECH_MODEL` | Defaults to `qwen3.8-flash-next`; pin it per run for reproducibility |
| `REQUIRE_TOOL_CALLS` | `1` by default; see the enforcement layers above |
| `MAX_TOOL_ROUNDS` | Bounds the agent loop (default 6) |
| `UNIPROT_BASE_URL` | Defaults to `https://rest.uniprot.org` |

## Check

```bash
uv run python -m protein_agent --list-models   # live GET /v1/models
uv run python -m protein_agent --dry-run "…"   # loop logic, network stubbed
uv run pytest -q                               # grounding tests, no network
uv run ruff check .
```

## Status

Verified against the live gateway: `tool_choice: "required"` forces a correct
`fetch_uniprot_entry(accession="P04637")` call on all three chat models, and the loop then
accepts only tool-backed text (an invented accession is blocked). UniProt retrieval is
implemented. Structural questions now lazily start the pinned official PDBe API and Search
MCP stdio servers, discover their live tool names and schemas, and route calls
without hard-coded generated API-tool names. Search fan-out is capped at 25 rows and every
result records its MCP server, tool, arguments, and retrieval time. The OOD report also
shows a bounded ranked PDBe mapping after resolving a verified UniProt accession. Remaining
work includes streaming and structured output in the CLI.

The OOD report-only assistant is temporarily disabled. Foldseek auto-reporting only runs
when a verified accession has a PDBe-mapped experimental structure with an accessible
chain: the report downloads that structure from PDBe, isolates the mapped chain, and
submits its mmCIF coordinates to Foldseek (`pdb100`, 3Di/AA mode). Raw FASTA is **never**
submitted as a 3D structure. Up to five structural-hit candidates and the query PDB/chain
appear in a collapsed section; bounded scores and provenance are saved in `report.json`.
No mapped structure, an inaccessible structure, service failure, or a busy search is
reported as skipped/unavailable/timeout, **not** as absence of structural neighbours.
Foldseek hits are candidates, not protein identity. The CLI uses `--foldseek` to request
this optional structure-backed path.
