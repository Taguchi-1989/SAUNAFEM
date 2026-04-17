"""Report generation (plots, tables, pass/fail)."""

from __future__ import annotations

import json
from pathlib import Path

from harness.heat_balance_parser import HeatBalance
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


def report_to_json(
    report: ValidationReport, output_path: Path | None = None, indent: int = 2
) -> str:
    """Generate a JSON validation report."""
    text = json.dumps(report_to_dict(report), indent=indent)
    if output_path:
        output_path.write_text(text + "\n", encoding="utf-8")
    return text


def heat_balance_to_markdown(balance: HeatBalance) -> str:
    """Generate a Markdown heat balance summary table."""
    abs_input = abs(balance.heater_input_W) or 1.0

    lines = [
        "## Heat Balance Summary",
        "",
        "| Component | Value [W] | % of Input |",
        "|-----------|-----------|------------|",
        f"| Heater input | {balance.heater_input_W:+.0f} | 100.0% |",
        f"| Wall losses | {balance.wall_loss_W:+.0f} | {balance.wall_loss_W / abs_input * 100:.1f}% |",
    ]

    if abs(balance.vent_loss_W) > 0.1:
        lines.append(
            f"| Vent losses | {balance.vent_loss_W:+.0f} | {balance.vent_loss_W / abs_input * 100:.1f}% |"
        )

    lines.extend([
        f"| **Imbalance** | **{balance.imbalance_W:+.0f}** | **{balance.imbalance_pct:.1f}%** |",
        "",
    ])

    if balance.vol_avg_T > 0:
        t_c = balance.vol_avg_T - 273.15
        lines.append(f"Volume-averaged T: {balance.vol_avg_T:.1f} K ({t_c:.1f} C)")
        lines.append("")

    if balance.patch_fluxes:
        lines.extend([
            "### Per-patch breakdown",
            "",
            "| Patch | Flux [W] |",
            "|-------|----------|",
        ])
        for patch, flux in sorted(balance.patch_fluxes.items()):
            lines.append(f"| {patch} | {flux:+.1f} |")
        lines.append("")

    return "\n".join(lines)
