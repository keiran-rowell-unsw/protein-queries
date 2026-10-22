"""The agent loop: model request -> tool calls -> retrieval -> proof, repeated.

The loop is bounded, refuses to summarise an empty evidence set, and hands the final text
to `proof.prove_answer`, so no answer reaches the user without tool backing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

from .config import Settings
from .gateway import ChatResponse, ResTechGateway
from .prompts import build_system_prompt
from .proof import AnswerProofError, Evidence, prove_answer
from .retrieval import PdbeClient, UniProtClient, dispatch
from .tools import all_tools


@dataclass(slots=True)
class RunResult:
    answer: str
    evidence: Evidence
    model: str
    rounds: int
    usage: dict[str, Any]
    stop_reason: str


class ProteinAgent:
    def __init__(
        self,
        settings: Settings,
        *,
        gateway: ResTechGateway | None = None,
        uniprot: UniProtClient | None = None,
        pdbe: PdbeClient | None = None,
        pdbe_tools: list[dict[str, Any]] | None = None,
    ) -> None:
        self.settings = settings
        self.gateway = gateway or ResTechGateway(settings)
        self.uniprot = uniprot or UniProtClient(settings)
        # PDBe is opt-in: the MCP server is not installed on every node, and a stub that
        # raises would teach the model to distrust working tools.
        self.pdbe = pdbe
        self.tools = all_tools(pdbe_tools=pdbe_tools)

    def run(
        self,
        question: str,
        *,
        min_calls: int = 1,
        on_unproven: str = "block",
        extra_context: str | None = None,
    ) -> RunResult:
        messages: list[dict[str, Any]] = [
            {"role": "system", "content": build_system_prompt(extra_context)},
            {"role": "user", "content": question},
        ]
        evidence = Evidence()
        usage: dict[str, Any] = {}
        model_used = self.settings.model

        for round_index in range(1, self.settings.max_tool_rounds + 1):
            # tool_choice="required" on every request: the model cannot answer from
            # memory, because replying without a tool call is not an available action.
            response: ChatResponse = self.gateway.chat(
                messages,
                tools=self.tools,
                tool_choice="required" if self.settings.require_tool_calls else "auto",
            )
            model_used = response.model or model_used
            usage = response.usage or usage

            if response.tool_calls:
                messages.append(
                    {
                        "role": "assistant",
                        "content": response.content,
                        "tool_calls": response.tool_calls,
                    }
                )
                for call in response.tool_calls:
                    name, arguments = _parse_call(call)
                    evidence.add_call(name, arguments)
                    result = dispatch(name, arguments, uniprot=self.uniprot, pdbe=self.pdbe)
                    evidence.add_result(json.dumps(result, ensure_ascii=False))
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": call.get("id", name),
                            "content": json.dumps(result, ensure_ascii=False),
                        }
                    )
                continue

            # No tool calls: this is the final answer.
            if evidence.call_count < min_calls:
                raise AnswerProofError(
                    f"Refusing to summarise: {evidence.call_count} tool call(s) made, "
                    f"{min_calls} required. Nothing was retrieved, so nothing can be said."
                )
            answer = prove_answer(
                response.content, evidence, min_calls=min_calls, on_unproven=on_unproven
            )
            return RunResult(
                answer=answer,
                evidence=evidence,
                model=model_used,
                rounds=round_index,
                usage=usage,
                stop_reason=response.finish_reason or "stop",
            )

        raise AnswerProofError(
            f"Tool-call budget exhausted after {self.settings.max_tool_rounds} rounds; "
            "no grounded answer was produced."
        )


def _parse_call(call: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    function = call.get("function") or {}
    name = str(function.get("name", ""))
    raw = function.get("arguments") or "{}"
    if isinstance(raw, dict):
        return name, raw
    try:
        return name, json.loads(raw)
    except json.JSONDecodeError as exc:
        raise AnswerProofError(f"Tool {name!r} received unparseable arguments: {exc}") from exc


__all__ = ["ProteinAgent", "RunResult"]
