"""Shared test fixtures for SaunaFlow."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture()
def project_root() -> Path:
    """Return the project root directory."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture()
def configs_dir(project_root: Path) -> Path:
    """Return the configs directory."""
    return project_root / "configs"


@pytest.fixture()
def sample_case_path(configs_dir: Path) -> Path:
    """Return path to the sample dry sauna case YAML."""
    return configs_dir / "cases" / "dry_sauna_steady.yaml"
