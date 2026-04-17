"""Batch case runner and comparison report generator."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from harness.heat_balance_parser import HeatBalance, compute_heat_balance, parse_vol_average_t, parse_wall_heat_flux
from harness.kpi import KPIResult, evaluate_all_kpis
from harness.probe_parser import parse_probe_file, get_steady_state_values
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
        try:
            result = solve_two_zone(yaml_path, n_profile=n_profile, max_iter=max_iter)
        except Exception as exc:
            import sys
            print(f"WARNING: Skipping {yaml_path.name}: {exc}", file=sys.stderr)
            continue

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


@dataclass
class OpenFOAMCaseResult:
    """Parsed results from one OpenFOAM case directory."""

    case_name: str
    case_dir: Path
    probe_values: dict[str, float] = field(default_factory=dict)
    heat_balance: HeatBalance | None = None


def parse_openfoam_case(
    case_dir: Path,
    probe_names: list[str],
    case_name: str | None = None,
) -> OpenFOAMCaseResult:
    """Parse probe and heat balance data from an OpenFOAM case directory.

    Args:
        case_dir: Path to the OpenFOAM case directory (with postProcessing/).
        probe_names: Ordered list of probe names.
        case_name: Display name (defaults to directory name).
    """
    name = case_name or case_dir.name

    # Parse probes
    probe_values: dict[str, float] = {}
    probe_dirs = sorted(
        (case_dir / "postProcessing" / "probes").glob("[0-9]*"),
    ) if (case_dir / "postProcessing" / "probes").exists() else []

    if probe_dirs:
        t_file = probe_dirs[-1] / "T"
        if t_file.exists():
            probe_data = parse_probe_file(t_file, probe_names)
            probe_values = get_steady_state_values(probe_data)

    # Parse heat balance
    wall_fluxes = parse_wall_heat_flux(case_dir)
    vol_avg = parse_vol_average_t(case_dir)
    heat_balance = compute_heat_balance(wall_fluxes, vol_avg) if wall_fluxes else None

    return OpenFOAMCaseResult(
        case_name=name,
        case_dir=case_dir,
        probe_values=probe_values,
        heat_balance=heat_balance,
    )


def compare_openfoam_results(
    case_dirs: list[Path],
    probe_names: list[str],
    case_names: list[str] | None = None,
) -> str:
    """Parse and compare results from multiple OpenFOAM case directories.

    Returns Markdown comparison report.
    """
    names = case_names or [d.name for d in case_dirs]
    results = [
        parse_openfoam_case(d, probe_names, n)
        for d, n in zip(case_dirs, names)
    ]

    lines = [
        "# OpenFOAM A/B Comparison Report",
        "",
        f"**Cases:** {len(results)}",
        "",
    ]

    # Probe temperature comparison
    lines.extend([
        "## Probe Temperatures (steady-state)",
        "",
        "| Probe | " + " | ".join(r.case_name for r in results) + " |",
        "| ----- | " + " | ".join("---" for _ in results) + " |",
    ])
    for pname in probe_names:
        row = f"| {pname} |"
        for r in results:
            val = r.probe_values.get(pname)
            if val is not None:
                row += f" {val - 273.15:.1f} C |"
            else:
                row += " -- |"
        lines.append(row)

    # Heat balance comparison
    lines.extend(["", "## Heat Balance", ""])
    hb_header = "| Metric | " + " | ".join(r.case_name for r in results) + " |"
    hb_sep = "| ------ | " + " | ".join("---" for _ in results) + " |"
    lines.extend([hb_header, hb_sep])

    for label, getter in [
        ("Heater input [W]", lambda hb: f"{hb.heater_input_W:+.0f}"),
        ("Wall losses [W]", lambda hb: f"{hb.wall_loss_W:+.0f}"),
        ("Vent losses [W]", lambda hb: f"{hb.vent_loss_W:+.0f}"),
        ("Imbalance [W]", lambda hb: f"{hb.imbalance_W:+.0f}"),
        ("Imbalance [%]", lambda hb: f"{hb.imbalance_pct:.1f}%"),
        ("Vol-avg T [C]", lambda hb: f"{hb.vol_avg_T - 273.15:.1f}"),
    ]:
        row = f"| {label} |"
        for r in results:
            if r.heat_balance:
                row += f" {getter(r.heat_balance)} |"
            else:
                row += " -- |"
        lines.append(row)

    lines.append("")
    return "\n".join(lines)
