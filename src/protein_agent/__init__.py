"""Retrieval-grounded protein data assistant.

Inference: UNSW ResTech Local Inference gateway (OpenAI-compatible).
Evidence: UniProt REST + PDBe-KB MCP servers. No memory-only answers.
"""

from .config import Settings
from .gateway import ResTechGateway
from .proof import AnswerProofError, prove_answer
from .router import build_plan

__all__ = [
    "AnswerProofError",
    "ResTechGateway",
    "Settings",
    "build_plan",
    "prove_answer",
]

__version__ = "0.1.0"
