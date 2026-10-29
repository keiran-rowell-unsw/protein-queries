"""The core guarantee: nothing is said without a tool call behind it."""

from __future__ import annotations

import json
from typing import Any

import httpx
import pytest
import respx

from protein_agent.agent import ProteinAgent
from protein_agent.config import Settings
from protein_agent.gateway import ResTechGateway
from protein_agent.proof import AnswerProofError, Evidence, prove_answer, unsupported_claims
from protein_agent.retrieval import UniProtClient

TOOL_RESULT = {
    "primaryAccession": "P04637",
    "entryType": "Reviewed entry (UniProtKB/Swiss-Prot)",
    "sequence": {"length": 393},
}


def _assistant_tool_calls(call_id: str = "call_1") -> dict[str, Any]:
    return {
        "id": "chatcmpl-1",
        "choices": [
            {
                "finish_reason": "tool_calls",
                "message": {
                    "role": "assistant",
                    "content": "",
                    "tool_calls": [
                        {
                            "id": call_id,
                            "type": "function",
                            "function": {
                                "name": "fetch_uniprot_entry",
                                "arguments": json.dumps({"accession": "P04637"}),
                            },
                        }
                    ],
                },
            }
        ],
        "model": "qwen3.8-flash-next",
        "usage": {"prompt_tokens": 10},
    }


def _assistant_text(text: str) -> dict[str, Any]:
    return {
        "id": "chatcmpl-2",
        "choices": [
            {"finish_reason": "stop", "message": {"role": "assistant", "content": text}}
        ],
        "model": "qwen3.8-flash-next",
        "usage": {"prompt_tokens": 20},
    }


def _settings() -> Settings:
    return Settings(
        base_url="https://inference.test/v1",
        uniprot_base_url="https://uniprot.test",
        api_key="sk-test",
        max_tool_rounds=3,
    )


def _agent(settings: Settings) -> ProteinAgent:
    """Agent wired to mocked endpoints; no test below touches a real server."""
    gateway_client = httpx.Client(base_url="https://inference.test/v1")
    uniprot_client = httpx.Client(base_url="https://uniprot.test")
    return ProteinAgent(
        settings,
        gateway=ResTechGateway(settings, client=gateway_client),
        uniprot=UniProtClient(settings, client=uniprot_client),
    )


@pytest.fixture()
def agent_with_mock(respx_mock: respx.MockRouter):
    captured: list[dict[str, Any]] = []
    responses = [
        _assistant_tool_calls(),
        _assistant_text("P04637 is a reviewed tumour suppressor (length 393)."),
    ]

    def chat_handler(request: httpx.Request) -> httpx.Response:
        captured.append(json.loads(request.content))
        return httpx.Response(200, json=responses[min(len(captured) - 1, len(responses) - 1)])

    respx_mock.route(url__startswith="https://inference.test/").side_effect = chat_handler
    respx_mock.get("https://uniprot.test/uniprotkb/P04637.json").mock(
        return_value=httpx.Response(200, json=TOOL_RESULT)
    )
    return _agent(_settings()), captured


def test_happy_path_requires_and_uses_a_tool_call(agent_with_mock) -> None:
    agent, captured = agent_with_mock
    result = agent.run("function of human TP53?")
    assert result.evidence.call_count >= 1
    assert result.evidence.tool_names() == ["fetch_uniprot_entry"]
    assert "P04637" in result.answer
    # First request must forbid a memory-only reply; the follow-up that writes the
    # answer runs with tools offered but not forced, or it could never conclude.
    assert captured[0]["tool_choice"] == "required"
    assert captured[0]["tools"]
    assert all(p["tool_choice"] in {"required", "auto"} for p in captured)


def test_memory_only_answer_is_blocked(agent_with_mock) -> None:
    agent, _ = agent_with_mock
    with pytest.raises(AnswerProofError) as exc:
        agent.run("function of human TP53?", min_calls=5)
    assert "tool call" in str(exc.value)


def test_invented_accession_never_reaches_the_user() -> None:
    evidence = Evidence()
    evidence.add_call("fetch_uniprot_entry", {"accession": "P04637"})
    evidence.add_result(json.dumps(TOOL_RESULT))
    problems = unsupported_claims("See P04637 and also Q12345, plus PDB 6ABC.", evidence)
    assert any("Q12345" in p for p in problems)
    assert any("6ABC" in p for p in problems)
    assert not any("P04637" in p for p in problems)
    with pytest.raises(AnswerProofError):
        prove_answer("See P04637 and also Q12345.", evidence)


def test_annotated_residue_claim_without_entry_fetch_is_flagged() -> None:
    """A position asserted as a database annotation needs an entry behind it."""
    evidence = Evidence()
    evidence.add_call("search_uniprotkb", {"query": "gene:TP53"})
    evidence.add_result(json.dumps({"results": [{"primaryAccession": "P04637"}]}))
    problems = unsupported_claims(
        "The annotated DNA binding residue 120 is phosphorylated.", evidence
    )
    assert any("annotated residue position" in p for p in problems)


def test_computed_residue_position_is_not_flagged() -> None:
    """Positions from the user's own sequence are legitimate; the old rule blocked them."""
    evidence = Evidence()
    evidence.add_call("compute_sequence_properties", {"sequence": "MTM" * 20})
    evidence.add_result(
        json.dumps(
            {
                "hydropathy_extremes": {
                    "min": [14, -0.9158],
                    "max": [66, 4.2789],
                }
            }
        )
    )
    problems = unsupported_claims(
        "The hydropathy profile peaks at residue 66 and is most hydrophilic at residue 14.",
        evidence,
    )
    assert not any("residue" in p.lower() for p in problems), problems


def test_empty_evidence_cannot_summarise() -> None:
    with pytest.raises(AnswerProofError):
        prove_answer("TP53 is a transcription factor.", Evidence())


def test_loop_budget_is_bounded(respx_mock: respx.MockRouter) -> None:
    """A model that keeps calling tools must be stopped, not followed forever."""
    settings = _settings()
    settings.max_tool_rounds = 2
    respx_mock.route(url__startswith="https://inference.test/").mock(
        return_value=httpx.Response(200, json=_assistant_tool_calls())
    )
    respx_mock.get(url__startswith="https://uniprot.test/").mock(
        return_value=httpx.Response(200, json=TOOL_RESULT)
    )
    agent = _agent(settings)
    with pytest.raises(AnswerProofError, match="budget exhausted"):
        agent.run("tell me everything about TP53")
