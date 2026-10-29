#!/bin/bash
# protein-queries agent launcher
#
# Sources the uv virtualenv at /srv/scratch/sbf/protein_queries/.venv so the assistant
# runs on the pinned interpreter, then hands you the grounded protein CLI. The ResTech
# key is never hard-coded here: it is read from the environment, or from a git-ignored
# .env in the project, or you are prompted for it once.
#
#   bash env.sh              # interactive chat with the agent
#   bash env.sh "question"    # one-shot answer, then return
#   bash env.sh -- status     # report venv/model/tool-call state
#   bash env.sh -- tests      # run the grounding suite
#   bash env.sh -- list-models
#   bash env.sh -- any CLI flags the agent accepts

set -euo pipefail

PROJECT_DIR="/srv/scratch/sbf/protein_queries"
VENV_DIR="${PROTEIN_AGENT_VENV:-$PROJECT_DIR/.venv}"
GATEWAY="${RESTECH_BASE_URL:-https://inference.restech.unsw.edu.au/v1}"
MODEL="${RESTECH_MODEL:-qwen3.8-flash-next}"
MIN_TOOL_CALLS="${PROTEIN_AGENT_MIN_CALLS:-1}"

cd "$PROJECT_DIR"

# --- environment -------------------------------------------------------------
if [[ ! -f "$VENV_DIR/bin/activate" ]]; then
  echo "No virtualenv at $VENV_DIR — creating it with uv..." >&2
  command -v uv >/dev/null 2>&1 || {
    echo "uv not found on PATH; cannot build the environment." >&2
    exit 1
  }
  uv venv "$VENV_DIR" >&2
  uv sync --project "$PROJECT_DIR" --extra dev >&2
fi
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

if [[ -f "$PROJECT_DIR/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "$PROJECT_DIR/.env"
  set +a
fi

# The key must exist, but never print it.
if [[ -z "${RESTECH_API_KEY:-}" ]]; then
  if [[ -t 0 ]]; then
    read -r -s -p "ResTech Local Inference key (sk-…, not stored): " RESTECH_API_KEY
    echo >&2
    export RESTECH_API_KEY
  else
    echo "RESTECH_API_KEY is not set and this is not a terminal — cannot run inference." >&2
    exit 1
  fi
fi

export RESTECH_BASE_URL="$GATEWAY"
export RESTECH_MODEL="$MODEL"
export REQUIRE_TOOL_CALLS=1   # hard requirement; see AGENTS.md
export PYTHONUNBUFFERED=1

if [[ "${1:-}" == "--" ]]; then
  shift
  case "${1:-}" in
    status)
      cat <<EOF
python        : $(command -v python) ($(python -V 2>&1))
VIRTUAL_ENV   : ${VIRTUAL_ENV:-<none>}
gateway       : $RESTECH_BASE_URL
model         : $RESTECH_MODEL
tool calls    : required, min $MIN_TOOL_CALLS per answer (REQUIRE_TOOL_CALLS=$REQUIRE_TOOL_CALLS)
api key       : $([[ -n "${RESTECH_API_KEY:-}" ]] && echo "set (not shown)" || echo MISSING)
uniProt       : ${UNIPROT_BASE_URL:-https://rest.uniprot.org}
pdbe MCP      : not configured on this node (structural claims stay unverified)
EOF
      exit 0
      ;;
    tests)
      exec pytest -q "${@:2}"
      ;;
    list-models)
      exec python -m protein_agent --list-models
      ;;
    ""|chat)
      :
      ;;
    *)
      exec protein-agent "$@"
      ;;
  esac
fi

# --- interactive/one-shot chat ------------------------------------------------
if [[ $# -gt 0 ]]; then
  exec python -m protein_agent --model "$MODEL" --min-calls "$MIN_TOOL_CALLS" "$*"
fi

echo "protein-queries agent · $MODEL @ $RESTECH_BASE_URL"
echo "Answers require a tool call; nothing is recalled from model memory."
echo "Type a question ('help' for commands, Ctrl-D to exit)."
echo

while true; do
  printf '\033[35myou\033[0m › '
  read -r line || break
  case "${line%% *}" in
    "") continue ;;
    quit|exit|bye) break ;;
    help)
      cat <<EOF
  <question>        grounded answer (UniProt REST + PDBe MCP tools)
  models            live GET /v1/models roster
  status            interpreter, model, enforcement state
  Ctrl-D / quit     leave
EOF
      continue ;;
    models)
      python -m protein_agent --list-models; continue ;;
    status)
      bash "${BASH_SOURCE[0]}" -- status; continue ;;
  esac
  printf '\033[36magent\033[0m › '
  python -m protein_agent --model "$MODEL" --min-calls "$MIN_TOOL_CALLS" --explain "$line"
done
