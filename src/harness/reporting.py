"""Report generation (plots, tables, pass/fail)."""

from __future__ import annotations

from pathlib import Path

from harness.validation import ValidationReport


def report_to_markdown(
    report: ValidationReport, output_path: Path | None = None
) -> str:
    """Generate a Markdown validation report."""
    lines = [
        f"# Validation Report: {report.case_name}",
        "",
        f"**Overall: {'PASS' if report.overall_pass else 'FAIL'}**",
        "",
        f"- Mean absolute error: {report.mean_abs_error:.2f}",
        f"- Max absolute error: {report.max_abs_error:.2f}",
        "",
        "## Per-probe results",
        "",
        "| Probe | Simulation | Experiment | Error | Rel Error | Status |",
        "|-------|-----------|------------|-------|-----------|--------|",
    ]

    for p in report.points:
        status = "PASS" if p.within_tol else "FAIL"
        lines.append(
            f"| {p.probe_name} | {p.sim_value:.2f} | {p.exp_value:.2f} | "
            f"{p.error:+.2f} | {p.rel_error:.2%} | {status} |"
        )

    text = "\n".join(lines) + "\n"

    if output_path:
        output_path.write_text(text, encoding="utf-8")

    return text


def report_to_dict(report: ValidationReport) -> dict:
    """Convert report to a JSON-serializable dict."""
    return {
        "case_name": report.case_name,
        "overall_pass": report.overall_pass,
        "mean_abs_error": report.mean_abs_error,
        "max_abs_error": report.max_abs_error,
        "points": [
            {
                "probe": p.probe_name,
                "sim": p.sim_value,
                "exp": p.exp_value,
                "error": p.error,
                "rel_error": p.rel_error,
                "within_tol": p.within_tol,
            }
            for p in report.points
        ],
    }
