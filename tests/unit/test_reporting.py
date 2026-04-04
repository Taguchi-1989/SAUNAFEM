"""Tests for reporting module."""

from __future__ import annotations

from pathlib import Path

from harness.reporting import report_to_dict, report_to_markdown
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
