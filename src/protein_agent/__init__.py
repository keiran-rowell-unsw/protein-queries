"""Retrieval-grounded protein data assistant.

Inference: UNSW ResTech Local Inference gateway (OpenAI-compatible).
Evidence: UniProt REST + PDBe-KB MCP servers. No memory-only answers.

Convenience exports are lazy: ``import protein_agent`` stays dependency-light so the
Passenger front end can use ``protein_agent.properties`` without pulling httpx/pydantic.
"""

from __future__ import annotations

from typing import Any

__version__ = "0.1.0"

_LAZY: dict[str, tuple[str, str]] = {
    "Settings": (".config", "Settings"),
    "ResTechGateway": (".gateway", "ResTechGateway"),
    "AnswerProofError": (".proof", "AnswerProofError"),
    "prove_answer": (".proof", "prove_answer"),
    "build_plan": (".router", "build_plan"),
}

__all__ = [*sorted(_LAZY), "__version__"]


def __getattr__(name: str) -> Any:
    target = _LAZY.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    module = importlib.import_module(target[0], __name__)
    value = getattr(module, target[1])
    globals()[name] = value
    return value
