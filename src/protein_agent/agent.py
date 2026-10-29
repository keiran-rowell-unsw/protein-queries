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


_COSMETIC_ARGS = {"label", "name", "title", "id", "comment"}


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
        self._pdbe_tools = pdbe_tools
        self.tools = all_tools(pdbe_tools=pdbe_tools)
        # Accessions the assistant itself verified in a prior turn and may reuse. Anything
        # else appearing in a sequence-lookup run is treated as an invented subject.
        self.trusted_accessions: frozenset[str] = frozenset()

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
        # Advertise only what this question can actually use. For a pasted sequence that
        # means the computed-property tools and nothing else: leaving search_uniprotkb
        # reachable lets the model "resolve" an unknown fragment by inventing an organism
        # or accession, which is the substitution this run exists to prevent.
        self.tools = self._tools_for(question)
        evidence = Evidence()
        seen_calls: dict[str, int] = {}
        usage: dict[str, Any] = {}
        model_used = self.settings.model

        for round_index in range(1, self.settings.max_tool_rounds + 1):
            # Required for the retrieval rounds: a memory-only reply is not an available
            # action. Auto once evidence is on record, so the model can actually conclude;
            # proof.prove_answer still blocks anything unbacked.
            require_call = self.settings.require_tool_calls and not evidence.calls
            response: ChatResponse = self.gateway.chat(
                messages,
                tools=self.tools,
                tool_choice="required" if require_call else "auto",
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
                    key = _call_key(name, arguments)
                    repeat_count = seen_calls.get(key, 0)
                    if repeat_count:
                        # A repeat means the model is stuck, not learning more. Rather than
                        # burn the budget, tell it plainly what is already on record.
                        messages.append(
                            {
                                "role": "tool",
                                "tool_call_id": call.get("id", name),
                                "content": json.dumps(
                                    {
                                        "already_returned": True,
                                        "times_called": repeat_count + 1,
                                        "instruction": (
                                            "You already received this exact result in this "
                                            "conversation. Do not call it again, and do not "
                                            "retry it with reordered or slightly changed "
                                            "arguments. Write your answer from the evidence "
                                            "already on record, or state plainly what could "
                                            "not be verified."
                                        ),
                                    }
                                ),
                            }
                        )
                        seen_calls[key] = repeat_count + 1
                        if repeat_count >= 1:
                            # It is asking twice for what it already has: stop offering the
                            # tool so the next turn must be the answer.
                            self.tools = [
                                t for t in self.tools
                                if t["function"]["name"] != name
                            ]
                        continue
                    seen_calls[key] = 1
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
                response.content,
                evidence,
                question=question,
                min_calls=min_calls,
                on_unproven=on_unproven,
                trusted_accessions=self.trusted_accessions,
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


    def _tools_for(self, question: str) -> list[dict[str, Any]]:
        from .router import build_plan

        plan = build_plan(question, pdbe_tools=self._pdbe_tools)
        allowed = set(plan.tool_names) | set(pdbe_tool_names_for(plan))
        candidates = all_tools(pdbe_tools=self._pdbe_tools)
        tools = [t for t in candidates if t["function"]["name"] in allowed]
        return tools or candidates


def pdbe_tool_names_for(plan) -> list[str]:
    return list(plan.extra_tool_names) if hasattr(plan, "extra_tool_names") else []


def _call_key(name: str, arguments: dict[str, Any]) -> str:
    """Identity of a call for repeat detection.

    Argument order, whitespace, and cosmetic labels must not count: a model that has
    finished retrieving but has not noticed tends to re-issue the same computation under a
    new label, which would otherwise look like novel work and silently drain the budget.
    """
    normalised = {
        str(k): (v.strip() if isinstance(v, str) else v)
        for k, v in arguments.items()
        if str(k).lower() not in _COSMETIC_ARGS
    }
    return f"{name}:{json.dumps(normalised, sort_keys=True, separators=(',', ':'))}"


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
