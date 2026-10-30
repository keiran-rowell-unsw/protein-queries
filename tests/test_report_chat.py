"""Report-only chat must read report data in every request."""

import json

import pytest

from protein_agent.config import Settings
from protein_agent.gateway import ChatResponse
from protein_agent.proof import AnswerProofError
from protein_agent.report_chat import READ_REPORT, ask_report, read_report_section

REPORT = {
    "header": "Example",
    "length_residues": 40,
    "molecular_weight_da": 4321.0,
    "isoelectric_point": 6.5,
    "net_charge_pH_7": -1.2,
    "gravy_kyte_doolittle": -0.3,
    "generated_at_utc": "2026-01-01T00:00:00+00:00",
    "uniprot": {
        "status": "verified",
        "accession": "Q12345",
        "organism": "Example organism",
        "entry_type": "reviewed",
        "features": [{"type": "Active site", "start": 10, "end": 10}],
        "comments": [],
        "cross_references": [],
        "protein_interactions": [],
        "source_url": "https://www.uniprot.org/uniprotkb/Q12345/entry",
    },
    "pdbe": {"status": "not_found", "reason": "No structures under this accession."},
    "disorder_prediction": {"status": "skipped", "reason": "Too short."},
}


class FakeGateway:
    def __init__(self, responses):
        self.responses = iter(responses)
        self.calls = []

    def chat(self, messages, **kwargs):
        self.calls.append((messages, kwargs))
        return next(self.responses)

    def close(self):
        pass


def tool_call(section="summary"):
    return ChatResponse(
        content="",
        finish_reason="tool_calls",
        tool_calls=[{
            "id": "report_1",
            "type": "function",
            "function": {
                "name": READ_REPORT,
                "arguments": json.dumps({"section": section}),
            },
        }],
    )


def test_report_chat_requires_reader_tool_and_answers_from_result():
    gateway = FakeGateway([
        tool_call(),
        ChatResponse(
            content="uniprot.accession: Q12345 [report:summary]",
            finish_reason="stop",
        ),
    ])
    answer = ask_report(
        REPORT,
        "Which accession is in this report?",
        settings=Settings(api_key="test"),
        gateway=gateway,
    )
    assert "Q12345" in answer
    assert gateway.calls[0][1]["tool_choice"] == "required"
    assert gateway.calls[1][1]["tool_choice"] == "auto"
    advertised = gateway.calls[0][1]["tools"]
    assert [item["function"]["name"] for item in advertised] == [READ_REPORT]


@pytest.mark.parametrize(
    "hallucination",
    [
        "The active site is residue 999. [report:summary]",
        "The molecular weight is 999999 Da. [report:summary]",
        "Figure 7 shows the structure. [report:summary]",
        "The report contains 40 active sites. [report:summary]",
        "This protein is a human kinase. [report:summary]",
    ],
)
def test_hallucinated_report_claims_are_blocked(hallucination):
    gateway = FakeGateway([
        tool_call(),
        ChatResponse(content=hallucination, finish_reason="stop"),
    ])
    with pytest.raises(AnswerProofError):
        ask_report(
            REPORT,
            "Walk me through the report.",
            settings=Settings(api_key="test"),
            gateway=gateway,
        )


def test_memory_only_report_answer_is_blocked():
    gateway = FakeGateway([
        ChatResponse(content="Q12345 is the accession.", finish_reason="stop")
    ])
    with pytest.raises(AnswerProofError, match="successful tool"):
        ask_report(
            REPORT,
            "Which accession?",
            settings=Settings(api_key="test"),
            gateway=gateway,
        )


def test_report_sections_are_bounded_and_do_not_expose_sequence():
    report = {
        **REPORT,
        "uniprot": {
            **REPORT["uniprot"],
            "canonical_sequence": "A" * 10_000,
            "features": [{"type": "Site", "start": i} for i in range(250)],
        },
    }
    section = read_report_section(report, "uniprot")
    assert "canonical_sequence" not in section["data"]
    assert len(section["data"]["features"]) == 200
    assert section["data"]["features_truncated"] is True

    report["disorder_prediction"] = {
        "status": "predicted",
        "per_residue_scores": [0.1] * 700,
    }
    predictions = read_report_section(report, "predictions")["disorder_prediction"]
    assert len(predictions["per_residue_scores"]) == 500
    assert predictions["scores_total"] == 700
    assert predictions["scores_truncated"] is True
