"""Tests for batch runner."""

from __future__ import annotations

from pathlib import Path

import pytest

from harness.batch import BatchReport, compare_openfoam_results, parse_openfoam_case, run_batch


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


class TestOpenFOAMComparison:
    def _make_case_dir(self, base: Path, name: str, probe_temps: list[float],
                       heat_fluxes: dict[str, float] | None = None,
                       vol_avg_t: float | None = None) -> Path:
        """Create a mock OpenFOAM case directory with postProcessing data."""
        case_dir = base / name
        # Probe data
        probe_dir = case_dir / "postProcessing" / "probes" / "0"
        probe_dir.mkdir(parents=True)
        lines = ["# Probe 0 (1.5 2.0 1.25)", "# Probe 1 (1.5 0.8 1.25)",
                 "# Time\t0\t1"]
        lines.append(f"300\t{probe_temps[0]}\t{probe_temps[1]}")
        (probe_dir / "T").write_text("\n".join(lines), encoding="utf-8")

        # Wall heat flux
        if heat_fluxes:
            whf_dir = case_dir / "postProcessing" / "wallHeatFlux" / "0"
            whf_dir.mkdir(parents=True)
            for patch, flux in heat_fluxes.items():
                (whf_dir / f"{patch}.dat").write_text(
                    f"# Time flux\n300\t{flux}\n", encoding="utf-8"
                )

        # Vol average T
        if vol_avg_t is not None:
            vat_dir = case_dir / "postProcessing" / "volAverageT" / "0"
            vat_dir.mkdir(parents=True)
            (vat_dir / "volFieldValue.dat").write_text(
                f"# Time volAverage(T)\n300\t{vol_avg_t}\n", encoding="utf-8"
            )

        return case_dir

    def test_parse_openfoam_case(self, tmp_path: Path) -> None:
        case_dir = self._make_case_dir(
            tmp_path, "case_a", [480.15, 353.15],
            {"heater_wall": 18000.0, "floor": -5000.0},
            vol_avg_t=358.2,
        )
        result = parse_openfoam_case(case_dir, ["upper_bench", "lower_bench"])
        assert result.probe_values["upper_bench"] == 480.15
        assert result.heat_balance is not None
        assert result.heat_balance.heater_input_W == 18000.0
        assert result.heat_balance.vol_avg_T == 358.2

    def test_compare_two_cases(self, tmp_path: Path) -> None:
        case_a = self._make_case_dir(
            tmp_path, "surfflux", [480.15, 353.15],
            {"heater_wall": 18000.0, "floor": -10000.0, "ceiling": -8000.0},
            vol_avg_t=358.2,
        )
        case_b = self._make_case_dir(
            tmp_path, "volsource", [440.15, 343.15],
            {"heater_wall": 0.0, "floor": -9000.0, "ceiling": -7000.0},
            vol_avg_t=350.0,
        )

        md = compare_openfoam_results(
            [case_a, case_b],
            ["upper_bench", "lower_bench"],
        )
        assert "A/B Comparison" in md
        assert "surfflux" in md
        assert "volsource" in md
        assert "upper_bench" in md
        assert "Heater input" in md

    def test_missing_postprocessing(self, tmp_path: Path) -> None:
        """Case with no postProcessing doesn't crash."""
        case_dir = tmp_path / "empty_case"
        case_dir.mkdir()
        result = parse_openfoam_case(case_dir, ["upper_bench"])
        assert result.probe_values == {}
        assert result.heat_balance is None
