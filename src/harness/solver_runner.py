"""Solver execution control."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from utils.wsl import wsl_exec

# Regex for OpenFOAM residual lines
_RESIDUAL_RE = re.compile(
    r"Solving for (\w+), Initial residual = ([\d.eE+-]+), "
    r"Final residual = ([\d.eE+-]+), No Iterations (\d+)"
)

# Default convergence thresholds
_DEFAULT_THRESHOLDS: dict[str, float] = {
    "Ux": 1e-4,
    "Uy": 1e-4,
    "Uz": 1e-4,
    "p_rgh": 1e-4,
    "T": 1e-5,
    "k": 1e-4,
    "epsilon": 1e-4,
}


@dataclass
class SolverResult:
    """Result of solver execution."""

    success: bool
    iterations: int
    converged: bool
    final_residuals: dict[str, float] = field(default_factory=dict)
    log_path: Path | None = None


def parse_residuals(log_text: str) -> list[dict[str, float]]:
    """Parse per-iteration final residuals from a solver log.

    Returns a list of dicts (one per time step), each mapping
    field name to final residual value.
    """
    iterations: list[dict[str, float]] = []
    current: dict[str, float] = {}

    for line in log_text.splitlines():
        m = _RESIDUAL_RE.search(line)
        if m:
            field_name = m.group(1)
            final_residual = float(m.group(3))

            # New iteration starts when we see a field we've already recorded
            if field_name in current:
                iterations.append(current)
                current = {}

            current[field_name] = final_residual

    # Don't forget the last iteration
    if current:
        iterations.append(current)

    return iterations


def check_convergence(
    residuals: list[dict[str, float]],
    thresholds: dict[str, float] | None = None,
) -> bool:
    """Check if the final residuals are below convergence thresholds.

    Args:
        residuals: List of per-iteration residual dicts.
        thresholds: Field name -> threshold. Defaults to standard values.

    Returns:
        True if all monitored fields' final residuals are below threshold.
    """
    if not residuals:
        return False

    if thresholds is None:
        thresholds = _DEFAULT_THRESHOLDS

    final = residuals[-1]
    for field_name, threshold in thresholds.items():
        if field_name in final and final[field_name] > threshold:
            return False
    return True


def run_solver(
    case_dir: Path,
    solver_name: str = "buoyantSimpleFoam",
    timeout: int = 3600,
) -> SolverResult:
    """Run the OpenFOAM solver on the case directory.

    Args:
        case_dir: Path to the OpenFOAM case directory.
        solver_name: Name of the solver executable.
        timeout: Maximum execution time in seconds.

    Returns:
        SolverResult with convergence information.
    """
    result = wsl_exec(solver_name, cwd=case_dir, timeout=timeout)
    log_text = result.stdout

    # Save log file
    log_path = case_dir / f"log.{solver_name}"
    log_path.write_text(log_text, encoding="utf-8")

    residuals = parse_residuals(log_text)
    converged = check_convergence(residuals)

    final_residuals = residuals[-1] if residuals else {}

    return SolverResult(
        success=result.returncode == 0,
        iterations=len(residuals),
        converged=converged,
        final_residuals=final_residuals,
        log_path=log_path,
    )
