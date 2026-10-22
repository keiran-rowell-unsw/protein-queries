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
implemented; PDBe MCP is stubbed with `NotImplementedError` because `pdbe-mcp-server` is not
installed on this node (`uvx` is available). Remaining work: wire the PDBe MCP client, add
streaming and structured output to the CLI, then put both behind a chat UI.
