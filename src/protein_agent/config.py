"""Runtime configuration.

The ResTech key is read from the environment only; it is never written to a log, an
answer, or an exception message.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

DEFAULT_BASE_URL = "https://inference.restech.unsw.edu.au/v1"
DEFAULT_MODEL = "qwen3.8-flash-next"
DEFAULT_UNIPROT_BASE_URL = "https://rest.uniprot.org"
# Live roster verified with GET /v1/models; the endpoint remains the authority.
KNOWN_CHAT_MODELS = ("qwen3.8-flash-next", "qwen3.8-27b", "glm-5.3-flash")


class ConfigError(RuntimeError):
    """Raised when the agent cannot start safely."""


def _env_flag(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off"}


@dataclass(slots=True)
class Settings:
    base_url: str = DEFAULT_BASE_URL
    model: str = DEFAULT_MODEL
    uniprot_base_url: str = DEFAULT_UNIPROT_BASE_URL
    max_tool_rounds: int = 6
    request_timeout: float = 120.0
    temperature: float = 0.0
    # The no-memory-only-answers rule, enforced in code. Disable for tests only.
    require_tool_calls: bool = True
    api_key: str = field(default="", repr=False)

    @classmethod
    def from_env(cls, *, require_key: bool = True) -> Settings:
        try:
            from dotenv import load_dotenv  # optional; uv sync installs it
        except ModuleNotFoundError:
            load_dotenv = None
        else:
            load_dotenv(override=False)

        settings = cls(
            base_url=os.getenv("RESTECH_BASE_URL", DEFAULT_BASE_URL).rstrip("/"),
            model=os.getenv("RESTECH_MODEL", DEFAULT_MODEL),
            uniprot_base_url=os.getenv(
                "UNIPROT_BASE_URL", DEFAULT_UNIPROT_BASE_URL
            ).rstrip("/"),
            max_tool_rounds=int(os.getenv("MAX_TOOL_ROUNDS", "6")),
            request_timeout=float(os.getenv("RESTECH_TIMEOUT", "120")),
            temperature=float(os.getenv("RESTECH_TEMPERATURE", "0")),
            require_tool_calls=_env_flag("REQUIRE_TOOL_CALLS", True),
            api_key=os.getenv("RESTECH_API_KEY", ""),
        )
        if require_key and not settings.api_key:
            raise ConfigError(
                "RESTECH_API_KEY is not set. Export it or put it in .env (git-ignored); "
                "request a key from ResTech Local Inference."
            )
        if settings.max_tool_rounds < 1:
            raise ConfigError("MAX_TOOL_ROUNDS must be >= 1")
        return settings

    def redact(self, text: str) -> str:
        """Strip the key from anything that could reach a terminal, log, or error."""
        if not self.api_key:
            return text
        return text.replace(self.api_key, "[REDACTED]")
