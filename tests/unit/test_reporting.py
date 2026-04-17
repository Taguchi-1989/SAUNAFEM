"""Tests for reporting module."""

from __future__ import annotations

from pathlib import Path

from harness.heat_balance_parser import HeatBalance
from harness.reporting import heat_balance_to_markdown, report_to_dict, report_to_markdown
from harness.validation import ValidationPoint, ValidationReport


class TestMarkdownReport:
    def test_generates_markdown(self) -> None:
        report = ValidationReport(
            case_name="test_case",
            points=[
                ValidationPoint("probe_a", 100.0, 102.0, -2.0, 0.0196, True),
                ValidationPoint("probe_b", 200.0, 190.0, 10.0, 0.0526, False),
            ],
            overall_pass=False,
            mean_abs_error=6.0,
            max_abs_error=10.0,
        )
        md = report_to_markdown(report)
        assert "test_case" in md
        assert "FAIL" in md
        assert "probe_a" in md
        assert "10.00" in md

    def test_writes_file(self, tmp_path: Path) -> None:
        report = ValidationReport("test", [], True, 0.0, 0.0)
        out = tmp_path / "report.md"
        report_to_markdown(report, out)
        assert out.exists()
        assert "PASS" in out.read_text(encoding="utf-8")


class TestHeatBalanceReport:
    def test_generates_markdown_table(self) -> None:
        hb = HeatBalance(
            heater_input_W=18000.0,
            wall_loss_W=-17500.0,
            vent_loss_W=-200.0,
            vol_avg_T=358.2,
            patch_fluxes={"heater_wall": 18000.0, "floor": -5000.0, "ceiling": -6000.0},
        )
        md = heat_balance_to_markdown(hb)
        assert "Heat Balance Summary" in md
        assert "+18000" in md
        assert "-17500" in md
        assert "-200" in md
        assert "358.2 K" in md
        assert "85.1 C" in md
        assert "heater_wall" in md

    def test_no_vent_row_when_zero(self) -> None:
        hb = HeatBalance(heater_input_W=18000.0, wall_loss_W=-17000.0)
        md = heat_balance_to_markdown(hb)
        assert "Vent" not in md

    def test_no_vol_avg_when_zero(self) -> None:
        hb = HeatBalance(heater_input_W=18000.0, wall_loss_W=-17000.0)
        md = heat_balance_to_markdown(hb)
        assert "Volume-averaged" not in md


class TestDictReport:
    def test_serializable(self) -> None:
        report = ValidationReport(
            case_name="test",
            points=[ValidationPoint("p", 1.0, 1.1, -0.1, 0.091, True)],
            overall_pass=True,
            mean_abs_error=0.1,
            max_abs_error=0.1,
        )
        d = report_to_dict(report)
        assert d["overall_pass"] is True
        assert len(d["points"]) == 1
        import json

        json.dumps(d)  # verify JSON-serializable
