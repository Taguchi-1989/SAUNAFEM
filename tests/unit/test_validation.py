"""Tests for validation module."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from harness.validation import (
    ValidationPoint,
    compare_probes,
    load_experimental_csv,
    time_average,
)


class TestLoadCSV:
    def test_load_sample(self) -> None:
        csv_path = Path("experiments/processed/dry_sauna_sample.csv")
        if not csv_path.exists():
            pytest.skip("Sample CSV not found")
        data = load_experimental_csv(csv_path)
        assert "upper_bench" in data
        assert "time" in data
        assert len(data["time"]) == 6

    def test_load_from_tmp(self, tmp_path: Path) -> None:
        csv = tmp_path / "test.csv"
        csv.write_text(
            "time,probe_a,probe_b\n0,100,200\n1,101,201\n", encoding="utf-8"
        )
        data = load_experimental_csv(csv)
        assert len(data["probe_a"]) == 2
        assert data["probe_b"][1] == 201


class TestTimeAverage:
    def test_full_average(self) -> None:
        times = np.array([0, 1, 2, 3, 4])
        values = np.array([10, 20, 30, 40, 50])
        assert time_average(values, times) == 30.0

    def test_windowed(self) -> None:
        times = np.array([0, 1, 2, 3, 4])
        values = np.array([10, 20, 30, 40, 50])
        avg = time_average(values, times, start=2, end=4)
        assert avg == 40.0  # mean of [30, 40, 50]


class TestCompareProbes:
    def test_all_within_tol(self) -> None:
        sim = {"a": 100.0, "b": 200.0}
        exp = {"a": 102.0, "b": 198.0}
        report = compare_probes(sim, exp, default_tol=5.0)
        assert report.overall_pass is True
        assert len(report.points) == 2

    def test_one_outside_tol(self) -> None:
        sim = {"a": 100.0, "b": 200.0}
        exp = {"a": 110.0, "b": 198.0}
        report = compare_probes(sim, exp, default_tol=5.0)
        assert report.overall_pass is False

    def test_missing_probe(self) -> None:
        sim = {"a": 100.0}
        exp = {"b": 200.0}
        report = compare_probes(sim, exp)
        assert len(report.points) == 0

    def test_custom_tolerance(self) -> None:
        sim = {"a": 100.0}
        exp = {"a": 108.0}
        report = compare_probes(sim, exp, tolerances={"a": 10.0})
        assert report.overall_pass is True
