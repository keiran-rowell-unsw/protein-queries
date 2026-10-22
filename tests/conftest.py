"""Shared fixtures. Live inference and network egress must be opt-in, never default."""

from __future__ import annotations

import pytest

from protein_agent.config import Settings


@pytest.fixture()
def settings() -> Settings:
    return Settings(
        base_url="https://inference.test/v1",
        model="qwen3.8-flash-next",
        uniprot_base_url="https://uniprot.test",
        api_key="sk-test-not-a-real-key",
        max_tool_rounds=3,
        require_tool_calls=True,
    )


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "live: hits the real ResTech gateway; needs RESTECH_API_KEY"
    )
