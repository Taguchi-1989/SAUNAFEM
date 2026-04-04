"""Batch case runner and comparison report generator."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from harness.kpi import KPIResult, evaluate_all_kpis
from harness.simple_solver import SimpleSolverResult, solve_two_zone


@dataclass
class BatchCaseResult:
    """Result of a single case in a batch run."""
    case_name: str
    case_yaml: Path
    result: SimpleSolverResult
    kpis: list[KPIResult]


@dataclass
class BatchReport:
    """Aggregated results from multiple cases."""
    cases: list[BatchCaseResult]

    def summary_table(self) -> str:
        """Generate a Markdown comparison table."""
        if not self.cases:
            return "No cases to compare.\n"

        lines = [
            "# Batch Comparison Report",
            "",
            f"**Cases:** {len(self.cases)}",
            "",
            "## Temperature Summary",
            "",
            "| Case | T_upper [C] | T_lower [C] | dT [K] | Interface [m] | Wall [C] | Humidity [g/kg] | Perceived [C] |",
            "| ---- | ----------- | ----------- | ------ | ------------- | -------- | --------------- | ------------- |",
        ]

        for c in self.cases:
            r = c.result
            t_up = r.upper_layer_temp - 273.15
            t_low = r.lower_layer_temp - 273.15
            dt = r.upper_layer_temp - r.lower_layer_temp
            wall = r.wall_inner_temp - 273.15
            hum = r.humidity_ratio * 1000
            perc = r.perceived_temp_upper
            lines.append(
                f"| {c.case_name} | {t_up:.1f} | {t_low:.1f} | {dt:.1f} | "
                f"{r.interface_height:.2f} | {wall:.1f} | {hum:.1f} | {perc:.1f} |"
            )

        lines.extend(["", "## KPI Comparison", ""])

        # Collect all KPI IDs across cases
        all_ids: list[str] = []
        for c in self.cases:
            for k in c.kpis:
                if k.kpi_id not in all_ids:
                    all_ids.append(k.kpi_id)

        header = "| KPI | " + " | ".join(c.case_name for c in self.cases) + " |"
        sep = "| --- | " + " | ".join("---" for _ in self.cases) + " |"
        lines.extend([header, sep])

        for kid in all_ids:
            row = f"| {kid} |"
            for c in self.cases:
                kpi = next((k for k in c.kpis if k.kpi_id == kid), None)
                if kpi:
                    status = f" [{kpi.pass_fail}]" if kpi.pass_fail else ""
                    row += f" {kpi.value} {kpi.unit}{status} |"
                else:
                    row += " — |"
            lines.append(row)

        return "\n".join(lines) + "\n"


def run_batch(
    case_yamls: list[Path],
    max_iter: int = 50000,
    n_profile: int = 80,
) -> BatchReport:
    """Run multiple cases and collect results with KPIs.

    Args:
        case_yamls: List of YAML case file paths.
        max_iter: Maximum solver iterations per case.
        n_profile: Number of vertical profile points.

    Returns:
        BatchReport with per-case results and KPIs.
    """
    cases: list[BatchCaseResult] = []

    for yaml_path in case_yamls:
        result = solve_two_zone(yaml_path, n_profile=n_profile, max_iter=max_iter)

        kpis = evaluate_all_kpis(
            probe_values=result.probe_values,
            perceived_temp_c=result.perceived_temp_upper,
            beta_aug=result.beta_aug_applied,
        )

        cases.append(BatchCaseResult(
            case_name=yaml_path.stem,
            case_yaml=yaml_path,
            result=result,
            kpis=kpis,
        ))

    return BatchReport(cases=cases)
