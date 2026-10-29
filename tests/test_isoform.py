r"""Isoform accessions must be accepted without weakening the sequence-verification gate.

Regression cover for a real defect: `retrieval.ACCESSION_RE` was written as an ordered
alternation `^[A-Z0-9]{6,10}$|^[A-Z0-9]{5,10}-\d{2}$`, so the greedy first branch won and
the *common* single-digit isoform (P04637-1) was rejected as malformed while the unusual
P04637-12 was accepted.

Live-verified against rest.uniprot.org (no recalled sequences): UniProt's
`primaryAccession` is stripped for isoforms identical to canonical (Q15465-1, 462 aa) and
echoed when they differ (P04637-2, 341 aa), so an equality test on it is not a sound gate.
"""

from __future__ import annotations

import pytest

from protein_agent import retrieval, uniprot_report
from protein_agent.uniprot_report import fetch_annotations

ISOFORMS = ["P04637-1", "P04637-2", "P04637-12", "Q15465-1", "P12110-3", "A0A1B3GSG3"]


@pytest.mark.parametrize("acc", ISOFORMS)
def test_retrieval_accepts_isoform_accessions(acc: str) -> None:
    assert retrieval.ACCESSION_RE.fullmatch(acc), f"{acc} must be well-formed"


@pytest.mark.parametrize(
    "acc",
    ["P04637-100", "P42", "P04637_1", "", "TOOLONGACCESSIONX", "P04637-1-2"],
)
def test_retrieval_still_rejects_malformed(acc: str) -> None:
    assert not retrieval.ACCESSION_RE.fullmatch(acc)


def test_single_digit_isoform_is_not_longer_rejected() -> None:
    """The specific regression: -1 failed while -12 passed under the old ordered regex."""
    assert bool(retrieval.ACCESSION_RE.fullmatch("P04637-1")) is True
    assert bool(retrieval.ACCESSION_RE.fullmatch("P04637-12")) is True


def test_uniparc_gate_stays_hyphen_intolerant() -> None:
    """The UniParc `CRC.accession` filter must never grow an isoform-suffix branch."""
    assert uniprot_report.ACCESSION.fullmatch("P04637")
    assert not uniprot_report.ACCESSION.fullmatch("P04637-1")


def test_user_gate_accepts_isoforms_but_not_overlong_suffixes() -> None:
    gate = uniprot_report.USER_ACCESSION
    assert gate.fullmatch("P04637-1") and gate.fullmatch("Q8WZ42-23")
    assert not gate.fullmatch("P04637-100")


class _StubClient:
    """Minimal UniProtClient stand-in returning a scripted entry."""

    def __init__(self, entry: dict) -> None:
        self._entry = entry
        self.requested: list[str] = []

    def entry(self, accession: str) -> dict:
        self.requested.append(accession)
        return self._entry


def _entry(primary: str, sequence: str) -> dict:
    return {
        "primaryAccession": primary,
        "sequence": {"value": sequence, "length": len(sequence)},
        "features": [],
        "photochemistry": {},
        "organism": {"scientificName": "Homo sapiens"},
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "entryVersion": "1",
    }


def test_isoform_verified_when_primary_is_stripped() -> None:
    """Q15465-1 shape: primary comes back without the suffix; exact sequence matches."""
    seq = "M" * 462
    stub = _StubClient(_entry("Q15465", seq))
    result = fetch_annotations(seq, "Q15465-1", client=stub)
    assert result["status"] == "verified", result


def test_isoform_verified_when_primary_is_echoed() -> None:
    """P04637-2 shape: primary echoes the isoform; exact sequence matches."""
    seq = "M" * 341
    stub = _StubClient(_entry("P04637-2", seq))
    assert fetch_annotations(seq, "P04637-2", client=stub)["status"] == "verified"


def test_wrong_isoform_sequence_is_mismatch_not_silently_accepted() -> None:
    """Relaxing the accession check must NOT let a different isoform's sequence through."""
    stub = _StubClient(_entry("P04637-1", "A" * 393))
    result = fetch_annotations("T" * 341, "P04637-1", client=stub)
    assert result["status"] == "mismatch"
    assert "does not exactly match" in result["reason"]


def test_unrelated_accession_is_not_accepted_by_the_suffix_tolerence() -> None:
    """Asking for P04637-2 but being served something from another gene must fail."""
    stub = _StubClient(_entry("Q15465", "G" * 462))
    assert fetch_annotations("G" * 462, "P04637-2", client=stub)["status"] == "mismatch"
