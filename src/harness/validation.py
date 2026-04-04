"""CFD vs experimental data comparison and validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class ValidationPoint:
    """Single comparison between simulation and experiment."""

    probe_name: str
    sim_value: float
    exp_value: float
    error: float  # sim - exp
    rel_error: float  # |error| / |exp| (or 0 if exp==0)
    within_tol: bool


@dataclass
class ValidationReport:
    """Full validation report for a case."""

    case_name: str
    points: list[ValidationPoint]
    overall_pass: bool
    mean_abs_error: float
    max_abs_error: float


def load_experimental_csv(csv_path: Path) -> dict[str, np.ndarray]:
    """Load experimental data from CSV.

    Expected format: first column is time [s], remaining columns are
    probe measurements. Header row has probe names.

    Returns dict mapping probe_name -> array of values.
    """
    data = np.genfromtxt(csv_path, delimiter=",", names=True, encoding="utf-8")
    if data.dtype.names is None:
        return {}
    data = np.atleast_1d(data)
    result = {}
    for name in data.dtype.names:
        result[name] = np.asarray(data[name], dtype=float)
    return result


def time_average(
    values: np.ndarray,
    times: np.ndarray,
    start: float | None = None,
    end: float | None = None,
) -> float:
    """Compute time-averaged value over a window.

    If start/end are None, averages over the full range.
    """
    if start is None:
        start = times[0]
    if end is None:
        end = times[-1]
    mask = (times >= start) & (times <= end)
    if not np.any(mask):
        return float(np.mean(values))
    return float(np.mean(values[mask]))


def compare_probes(
    sim_values: dict[str, float],
    exp_values: dict[str, float],
    tolerances: dict[str, float] | None = None,
    default_tol: float = 5.0,
    case_name: str = "unknown",
) -> ValidationReport:
    """Compare simulation probe values against experimental data.

    Args:
        sim_values: dict of probe_name -> simulated value (e.g., temperature in K)
        exp_values: dict of probe_name -> experimental value
        tolerances: per-probe absolute tolerance (K or relevant unit)
        default_tol: default tolerance if not specified per-probe

    Returns:
        ValidationReport with per-point and aggregate results.
    """
    if tolerances is None:
        tolerances = {}

    points = []
    for name in sim_values:
        if name not in exp_values:
            continue
        sim = sim_values[name]
        exp = exp_values[name]
        error = sim - exp
        rel = abs(error) / abs(exp) if abs(exp) > 1e-10 else 0.0
        tol = tolerances.get(name, default_tol)
        within = abs(error) <= tol
        points.append(
            ValidationPoint(
                probe_name=name,
                sim_value=sim,
                exp_value=exp,
                error=round(error, 3),
                rel_error=round(rel, 4),
                within_tol=within,
            )
        )

    if not points:
        return ValidationReport(
            case_name=case_name,
            points=[],
            overall_pass=False,
            mean_abs_error=0.0,
            max_abs_error=0.0,
        )

    errors = [abs(p.error) for p in points]
    return ValidationReport(
        case_name=case_name,
        points=points,
        overall_pass=all(p.within_tol for p in points),
        mean_abs_error=round(float(np.mean(errors)), 3),
        max_abs_error=round(float(max(errors)), 3),
    )
