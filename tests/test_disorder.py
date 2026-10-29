"""Predictions must be labelled, bounded, and unavailable rather than fabricated."""

import sys
from types import SimpleNamespace

from protein_agent.disorder import predict_disorder


def test_disorder_scores_and_boundaries(tmp_path, monkeypatch):
    calls = []

    def fake_predict(sequence, **options):
        calls.append(options)
        return SimpleNamespace(
            disorder=[0.75] * len(sequence),
            disordered_domain_boundaries=[[0, len(sequence)]],
        )

    monkeypatch.setitem(sys.modules, "metapredict", SimpleNamespace(predict_disorder=fake_predict))
    report = predict_disorder("Q" * 30, tmp_path / "disorder.png")
    assert report["status"] == "predicted"
    assert report["disordered_regions"] == [{"start": 1, "end": 30}]
    assert report["threshold"] == 0.5
    assert calls[0]["version"] == "V3" and calls[0]["round_values"] is False
    assert (tmp_path / "disorder.png").stat().st_size > 0


def test_invalid_prediction_is_not_reported(tmp_path, monkeypatch):
    monkeypatch.setitem(
        sys.modules, "metapredict",
        SimpleNamespace(predict_disorder=lambda *_args, **_kwargs: SimpleNamespace(
            disorder=[float("nan")], disordered_domain_boundaries=[]
        )),
    )
    result = predict_disorder("Q" * 30, tmp_path / "disorder.png")
    assert result["status"] == "unavailable"
    assert not (tmp_path / "disorder.png").exists()


def test_short_sequence_skips_without_calling_model(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "metapredict", None)
    result = predict_disorder("Q" * 24, tmp_path / "disorder.png")
    assert result["status"] == "skipped"
    assert "shorter than 25" in result["reason"]


def test_long_sequence_skips_prediction(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "metapredict", None)
    result = predict_disorder("Q" * 2001, tmp_path / "disorder.png")
    assert result["status"] == "skipped"
    assert "2000" in result["reason"]
