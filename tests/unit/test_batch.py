"""Tests for batch runner."""

from __future__ import annotations

from pathlib import Path

import pytest

from harness.batch import BatchReport, run_batch


class TestRunBatch:
    def test_batch_all_cases(self) -> None:
        yamls = sorted(Path("configs/cases").glob("*.yaml"))
        assert len(yamls) >= 3

        report = run_batch(yamls, max_iter=10000, n_profile=40)
        assert len(report.cases) == len(yamls)

        for c in report.cases:
            assert c.result.upper_layer_temp > 273.15
            assert len(c.kpis) >= 2  # at least K-01 and K-07

    def test_summary_table_format(self) -> None:
        yamls = sorted(Path("configs/cases").glob("*.yaml"))
        report = run_batch(yamls, max_iter=10000, n_profile=40)
        md = report.summary_table()

        assert "# Batch Comparison Report" in md
        assert "T_upper" in md
        assert "K-01" in md
        for y in yamls:
            assert y.stem in md

    def test_empty_batch(self) -> None:
        report = run_batch([], max_iter=100)
        assert len(report.cases) == 0
        md = report.summary_table()
        assert "No cases" in md
