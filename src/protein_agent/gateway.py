"""OpenAI-compatible client for the ResTech Local Inference gateway.

Only this module talks to the model. Retrieval services (UniProt, PDBe MCP) live in
`retrieval.py`; the agent loop that ties them together lives in `agent.py`.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any

import httpx

from .config import Settings


class GatewayError(RuntimeError):
    """Gateway unreachable, refusing, or returning something unusable."""

    def __init__(self, message: str, *, status: int | None = None) -> None:
        super().__init__(message)
        self.status = status


@dataclass(slots=True)
class ChatResponse:
    content: str
    tool_calls: list[dict[str, Any]] = field(default_factory=list)
    reasoning_content: str | None = None  # internal scratchpad; never show to a user
    finish_reason: str = ""
    model: str = ""
    usage: dict[str, Any] = field(default_factory=dict)


class ResTechGateway:
    """Thin, explicit client for POST /v1/chat/completions and GET /v1/models."""

    def __init__(self, settings: Settings, *, client: httpx.Client | None = None) -> None:
        self.settings = settings
        self._owns_client = client is None
        headers = {
            "Authorization": f"Bearer {settings.api_key}",
            "Content-Type": "application/json",
        }
        self._client = client or httpx.Client(
            base_url=settings.base_url,
            headers=headers,
            timeout=settings.request_timeout,
        )

    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def __enter__(self) -> ResTechGateway:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def list_models(self) -> list[str]:
        try:
            response = self._client.get("/models")
        except httpx.HTTPError as exc:
            raise GatewayError(f"Cannot reach inference gateway: {exc}") from exc
        if response.status_code == 401:
            raise GatewayError(
                "Inference gateway rejected the key (401). Check RESTECH_API_KEY.",
                status=401,
            )
        response.raise_for_status()
        return [m["id"] for m in response.json().get("data", [])]

    def chat(
        self,
        messages: list[dict[str, Any]],
        *,
        tools: list[dict[str, Any]] | None = None,
        tool_choice: str | dict[str, Any] | None = None,
        model: str | None = None,
        temperature: float | None = None,
        max_tokens: int | None = None,
    ) -> ChatResponse:
        """One completion request with bounded retries on 429/5xx."""
        payload: dict[str, Any] = {
            "model": model or self.settings.model,
            "messages": messages,
            "temperature": self.settings.temperature
            if temperature is None
            else temperature,
        }
        if tools:
            payload["tools"] = tools
            # "required" makes a memory-only answer structurally impossible, but it also
            # forces a tool call on the turn that should be the summary. Use it for the
            # retrieval rounds only, then "auto" for the final round (see agent.py).
            payload["tool_choice"] = tool_choice or "auto"
        if max_tokens:
            payload["max_tokens"] = max_tokens

        last_error = "inference gateway request failed"
        for attempt in range(3):
            try:
                response = self._client.post("/chat/completions", json=payload)
            except httpx.HTTPError as exc:
                last_error = self.settings.redact(f"Cannot reach inference gateway: {exc}")
            else:
                if response.status_code == 401:
                    raise GatewayError(
                        "Inference gateway rejected the key (401). Check RESTECH_API_KEY.",
                        status=401,
                    )
                if response.status_code in {429, 500, 502, 503, 504}:
                    last_error = f"gateway returned {response.status_code}"
                else:
                    response.raise_for_status()
                    return self._parse(response.json())
            time.sleep(min(2**attempt, 8))  # bounded backoff; no unbounded spinning
        raise GatewayError(f"{last_error} after retries — the service is unavailable.")

    @staticmethod
    def _parse(body: dict[str, Any]) -> ChatResponse:
        choices = body.get("choices") or []
        if not choices:
            raise GatewayError("Inference gateway returned no choices.")
        message = choices[0].get("message") or {}
        return ChatResponse(
            content=message.get("content") or "",
            tool_calls=message.get("tool_calls") or [],
            reasoning_content=message.get("reasoning_content"),
            finish_reason=choices[0].get("finish_reason") or "",
            model=body.get("model", ""),
            usage=body.get("usage") or {},
        )
