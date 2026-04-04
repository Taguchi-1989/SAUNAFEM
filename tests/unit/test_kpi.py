"""Tests for KPI calculation module."""

from __future__ import annotations

from harness.kpi import (
    compute_k01,
    compute_k02,
    compute_k03,
    compute_k04,
    compute_k05,
    compute_k06,
    compute_k07,
    evaluate_all_kpis,
    evaluate_phase1_kpis,
)


class TestComputeK01:
    def test_positive_stratification(self) -> None:
        result = compute_k01(upper_temp=353.0, lower_temp=323.0)
        assert result.kpi_id == "K-01"
        assert result.value == 30.0
        assert result.unit == "K"
        assert result.pass_fail == "pass"

    def test_no_stratification(self) -> None:
        result = compute_k01(upper_temp=300.0, lower_temp=300.0)
        assert result.value == 0.0
        assert result.pass_fail == "fail"

    def test_inverted_stratification(self) -> None:
        result = compute_k01(upper_temp=300.0, lower_temp=350.0)
        assert result.value == -50.0
        assert result.pass_fail == "fail"


class TestComputeK07:
    def test_relative_difference(self) -> None:
        result = compute_k07(upper_temp=360.0, lower_temp=320.0)
        assert result.kpi_id == "K-07"
        # (360 - 320) / 340 ≈ 0.1176
        assert abs(result.value - 0.1176) < 0.001
        assert result.unit == "-"
        assert result.pass_fail is None

    def test_equal_temperatures(self) -> None:
        result = compute_k07(upper_temp=300.0, lower_temp=300.0)
        assert result.value == 0.0

    def test_zero_temperatures(self) -> None:
        result = compute_k07(upper_temp=0.0, lower_temp=0.0)
        assert result.value == 0.0


class TestEvaluatePhase1Kpis:
    def test_returns_two_kpis(self) -> None:
        values = {"upper_bench": 358.2, "lower_bench": 327.5}
        results = evaluate_phase1_kpis(values)
        assert len(results) == 2
        assert results[0].kpi_id == "K-01"
        assert results[1].kpi_id == "K-07"

    def test_correct_values(self) -> None:
        values = {"upper_bench": 358.2, "lower_bench": 327.5}
        results = evaluate_phase1_kpis(values)
        assert results[0].value == 30.7  # 358.2 - 327.5
        assert results[0].pass_fail == "pass"

    def test_missing_probes(self) -> None:
        values = {}
        results = evaluate_phase1_kpis(values)
        assert results[0].value == 0.0
        assert results[0].pass_fail == "fail"


class TestKPI_K02:
    def test_peak_rise(self) -> None:
        series = [350.0, 360.0, 370.0, 365.0, 355.0]
        result = compute_k02(series, baseline_temp=350.0)
        assert result.kpi_id == "K-02"
        assert result.value == 20.0  # 370 - 350
        assert result.pass_fail == "pass"

    def test_no_rise(self) -> None:
        result = compute_k02([350.0, 350.0], baseline_temp=350.0)
        assert result.value == 0.0
        assert result.pass_fail == "fail"


class TestKPI_K03:
    def test_peak_humidity(self) -> None:
        series = [0.0, 0.01, 0.05, 0.03, 0.01]
        result = compute_k03(series)
        assert result.value == 50.0  # 0.05 * 1000


class TestKPI_K04:
    def test_arrival_time(self) -> None:
        times = [0, 5, 10, 15, 20]
        temps = [350, 355, 370, 365, 360]
        result = compute_k04(times, temps, event_time=0)
        assert result.value == 10.0  # peak at index 2, time=10


class TestKPI_K05:
    def test_wind_proxy(self) -> None:
        result = compute_k05(0.5)
        assert result.value > 0

    def test_no_aufguss(self) -> None:
        result = compute_k05(0.0)
        assert result.value == 0.0


class TestKPI_K06:
    def test_comfortable(self) -> None:
        result = compute_k06(50.0)
        assert "comfortable" in result.name

    def test_intense(self) -> None:
        result = compute_k06(90.0)
        assert "intense" in result.name


class TestEvaluateAllKPIs:
    def test_basic(self) -> None:
        results = evaluate_all_kpis(
            probe_values={"upper_bench": 370.0, "lower_bench": 340.0},
            perceived_temp_c=95.0,
        )
        ids = [r.kpi_id for r in results]
        assert "K-01" in ids
        assert "K-06" in ids
        assert "K-07" in ids
