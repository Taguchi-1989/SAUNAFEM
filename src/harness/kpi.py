"""KPI calculation (peak T, humidity, wind speed, arrival time)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class KPIResult:
    """Result of a single KPI computation."""

    kpi_id: str
    name: str
    value: float
    unit: str
    pass_fail: str | None = None  # "pass", "fail", or None


def compute_k01(upper_temp: float, lower_temp: float) -> KPIResult:
    """K-01: Steady-state temperature differential (upper - lower bench).

    A positive value indicates thermal stratification (hot air rises).
    """
    diff = upper_temp - lower_temp
    return KPIResult(
        kpi_id="K-01",
        name="Steady-state temperature differential",
        value=round(diff, 2),
        unit="K",
        pass_fail="pass" if diff > 0 else "fail",
    )


def compute_k07(upper_temp: float, lower_temp: float) -> KPIResult:
    """K-07: Upper/lower relative difference.

    Computed as (upper - lower) / mean_temperature.
    """
    mean = (upper_temp + lower_temp) / 2
    rel_diff = (upper_temp - lower_temp) / mean if mean != 0 else 0.0
    return KPIResult(
        kpi_id="K-07",
        name="Upper/lower relative difference",
        value=round(rel_diff, 4),
        unit="-",
        pass_fail=None,
    )


def compute_k02(
    t_upper_series: list[float] | np.ndarray, baseline_temp: float
) -> KPIResult:
    """K-02: Post-löyly peak temperature rise above baseline [K]."""
    peak = max(t_upper_series) if len(t_upper_series) > 0 else baseline_temp
    rise = peak - baseline_temp
    return KPIResult(
        kpi_id="K-02",
        name="Post-löyly peak temperature rise",
        value=round(rise, 2),
        unit="K",
        pass_fail="pass" if rise > 0.5 else "fail",
    )


def compute_k03(humidity_series: list[float] | np.ndarray) -> KPIResult:
    """K-03: Post-löyly peak absolute humidity [g/kg]."""
    peak_kg = max(humidity_series) if len(humidity_series) > 0 else 0.0
    peak_gkg = peak_kg * 1000  # convert to g/kg for readability
    return KPIResult(
        kpi_id="K-03",
        name="Post-löyly peak humidity",
        value=round(peak_gkg, 2),
        unit="g/kg",
        pass_fail=None,
    )


def compute_k04(
    time_series: list[float] | np.ndarray,
    t_upper_series: list[float] | np.ndarray,
    event_time: float = 0.0,
) -> KPIResult:
    """K-04: Time from event (löyly) to peak temperature [s]."""
    if len(t_upper_series) == 0 or len(time_series) == 0:
        return KPIResult(
            kpi_id="K-04",
            name="Peak arrival time",
            value=0.0,
            unit="s",
            pass_fail=None,
        )
    peak_idx = int(np.argmax(t_upper_series))
    arrival = time_series[peak_idx] - event_time
    return KPIResult(
        kpi_id="K-04",
        name="Peak arrival time",
        value=round(max(arrival, 0.0), 1),
        unit="s",
        pass_fail=None,
    )


def compute_k05(beta_aug: float) -> KPIResult:
    """K-05: Face-level wind speed proxy from Aufguss mixing coefficient.

    Approximate face wind speed from forced mixing rate.
    v ~ beta_aug / (rho * A_face), rough estimate.
    """
    rho = 0.9  # hot air density ~100C [kg/m3]
    a_face = 0.05  # effective face area [m2]
    v_proxy = beta_aug / (rho * a_face) if beta_aug > 0 else 0.0
    return KPIResult(
        kpi_id="K-05",
        name="Face-level wind speed (proxy)",
        value=round(v_proxy, 2),
        unit="m/s",
        pass_fail="pass" if v_proxy > 0.1 else None,
    )


def compute_k06(perceived_temp_c: float) -> KPIResult:
    """K-06: Simplified thermal stress index.

    Based on perceived temperature (includes humidity effect).
    Categories: <60C comfortable, 60-80C moderate, 80-100C intense, >100C extreme.
    """
    if perceived_temp_c < 60:
        category = "comfortable"
    elif perceived_temp_c < 80:
        category = "moderate"
    elif perceived_temp_c < 100:
        category = "intense"
    else:
        category = "extreme"
    return KPIResult(
        kpi_id="K-06",
        name=f"Thermal stress index ({category})",
        value=round(perceived_temp_c, 1),
        unit="C",
        pass_fail=None,
    )


def evaluate_phase1_kpis(probe_values: dict[str, float]) -> list[KPIResult]:
    """Compute all Phase 1 KPIs from probe steady-state values.

    Expects probe_values to contain 'upper_bench' and 'lower_bench' keys.
    """
    upper = probe_values.get("upper_bench", 0.0)
    lower = probe_values.get("lower_bench", 0.0)

    return [
        compute_k01(upper, lower),
        compute_k07(upper, lower),
    ]


def evaluate_all_kpis(
    probe_values: dict[str, float],
    t_upper_series: list[float] | np.ndarray | None = None,
    humidity_series: list[float] | np.ndarray | None = None,
    time_series: list[float] | np.ndarray | None = None,
    baseline_temp: float = 0.0,
    event_time: float = 0.0,
    beta_aug: float = 0.0,
    perceived_temp_c: float = 0.0,
) -> list[KPIResult]:
    """Compute all available KPIs."""
    upper = probe_values.get("upper_bench", 0.0)
    lower = probe_values.get("lower_bench", 0.0)

    results = [
        compute_k01(upper, lower),
        compute_k07(upper, lower),
    ]

    if t_upper_series is not None and len(t_upper_series) > 0:
        results.append(compute_k02(t_upper_series, baseline_temp))
    if humidity_series is not None and len(humidity_series) > 0:
        results.append(compute_k03(humidity_series))
    if time_series is not None and t_upper_series is not None:
        results.append(compute_k04(time_series, t_upper_series, event_time))
    if beta_aug > 0:
        results.append(compute_k05(beta_aug))
    if perceived_temp_c > 0:
        results.append(compute_k06(perceived_temp_c))

    return results
