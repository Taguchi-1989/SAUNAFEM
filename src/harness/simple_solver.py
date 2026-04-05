"""Two-zone plume model for sauna thermal stratification.

Based on the Morton-Taylor-Turner (MTT) entrainment model and
Zukoski's two-layer zone model used in building fire/HVAC engineering.

Physics:
  - Room divided into upper hot layer and lower cool layer
  - Heater drives a buoyant plume that entrains lower-layer air
  - Plume deposits mass and enthalpy into the upper layer
  - Both layers lose heat to surrounding walls
  - Layer interface height descends until steady state

This is NOT a CFD solver. It captures the dominant physics of
thermal stratification (mass conservation, plume entrainment,
layer energy balance) without solving Navier-Stokes.

References:
  - Morton, Taylor & Turner (1956), Proc. R. Soc. A 234:1-23
  - Zukoski (1978), NBS-GCR-78-150
  - Cooper (1982), NBSIR 82-2612
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import math

import numpy as np

from harness.schema import load_yaml


@dataclass
class SimpleSolverResult:
    """Result of two-zone model calculation."""

    y_positions: np.ndarray       # vertical positions [m]
    temperatures: np.ndarray      # temperature at each y [K]
    probe_values: dict[str, float]
    iterations: int
    converged: bool
    residual_history: list[float]
    # Zone model specific outputs
    interface_height: float       # hot/cold layer interface [m]
    upper_layer_temp: float       # upper layer temperature [K]
    lower_layer_temp: float       # lower layer temperature [K]
    plume_mass_flow: float        # plume mass flow at interface [kg/s]
    steam_mass_flow: float = 0.0      # peak evaporation rate [kg/s]
    total_steam_generated: float = 0.0 # cumulative steam mass [kg]
    beta_aug_applied: float = 0.0  # forced mixing coefficient used [kg/s]
    vent_mass_flow: float = 0.0        # ventilation mass flow rate [kg/s]
    wall_inner_temp: float = 293.15    # inner wall surface temperature [K]
    humidity_ratio: float = 0.0        # kg vapor / kg dry air
    relative_humidity: float = 0.0     # 0-1
    perceived_temp_upper: float = 0.0  # perceived temp at upper layer [C]
    perceived_temp_lower: float = 0.0  # perceived temp at lower layer [C]


@dataclass
class TransientResult:
    """Result of transient two-zone simulation with time-series data."""

    time: np.ndarray              # time points [s]
    t_upper_series: np.ndarray    # upper layer temp over time [K]
    t_lower_series: np.ndarray    # lower layer temp over time [K]
    z_int_series: np.ndarray      # interface height over time [m]
    humidity_series: np.ndarray   # humidity ratio over time [kg/kg]
    wall_temp_series: np.ndarray  # wall inner temp over time [K]
    perceived_upper_series: np.ndarray  # perceived temp upper [C]
    # Final state (same as steady)
    steady_result: SimpleSolverResult


def _compute_view_factors(
    room_w: float, room_d: float, room_h: float,
    heater_y: float, heater_h: float, heater_w: float,
) -> dict[str, float]:
    """Compute approximate view factors from heater to room surfaces.

    Uses solid-angle fractions from the heater centroid to each surface.
    A heater near the floor radiates more downward; near the ceiling,
    more upward. The opposite wall receives radiation based on the
    heater's angular extent.

    Returns dict with keys: 'floor', 'lower_walls', 'upper_walls',
    'ceiling', 'body'.  All values sum to ~1.0 (enclosure closure rule).
    """
    h_center = heater_y + heater_h / 2.0

    # Vertical angle fractions from heater centroid
    # Floor: angle subtended looking DOWN from heater
    dist_floor = max(h_center, 0.01)
    # Ceiling: angle subtended looking UP from heater
    dist_ceil = max(room_h - h_center, 0.01)

    # Heater is on the wall (x=0 face, width = room_d direction).
    # heater_w is the heater width along the wall (room_d direction).
    # A wider heater illuminates more of each surface.
    # Scale the characteristic length by heater coverage fraction.
    heater_coverage = min(heater_w / max(room_d, 0.01), 1.0)
    # Wider heater → larger effective radiating solid angle
    # Blend between point source (coverage→0) and full-wall source (coverage→1)
    coverage_factor = 0.3 + 0.7 * heater_coverage  # [0.3, 1.0]

    # Solid-angle fraction ≈ atan(characteristic_length / distance) / pi
    char_len = np.sqrt(room_w**2 + room_d**2) / 2.0

    f_floor_raw = np.arctan(char_len / dist_floor) / np.pi * coverage_factor
    f_ceiling_raw = np.arctan(char_len / dist_ceil) / np.pi * coverage_factor

    # Opposite wall fraction — vertical angle times horizontal coverage
    f_opposite = np.arctan(heater_h / max(room_w, 0.01)) / np.pi * coverage_factor

    # Remaining goes to side walls, split by position relative to heater
    f_sum = f_floor_raw + f_ceiling_raw + f_opposite
    if f_sum > 1.0:
        scale = 1.0 / f_sum
        f_floor_raw *= scale
        f_ceiling_raw *= scale
        f_opposite *= scale

    f_sides = max(0.0, 1.0 - f_floor_raw - f_ceiling_raw - f_opposite)

    # Side walls: split into lower/upper based on the heater's vertical
    # position. Below heater → "lower walls", above → "upper walls".
    # fraction_below = how much wall area is below heater center
    frac_below = h_center / room_h
    frac_above = 1.0 - frac_below

    # Lower zone receives: floor + below-portion of side walls + below-portion of opposite
    f_lower_walls = f_sides * frac_below + f_opposite * frac_below
    f_upper_walls = f_sides * frac_above + f_opposite * frac_above

    # --- Heater-to-body view factor ---
    # Person on upper bench, roughly 1.5 m from heater horizontally.
    # Approximate as small target (human front surface ~0.6 m wide, ~0.5 m high)
    # at distance d from heater.  F = A_person / (2*pi*d^2) (point-to-small-area).
    body_width = 0.6   # m, approximate torso width
    body_height = 0.5  # m, seated torso height
    d_body = max(np.sqrt(room_w**2 + (room_h * 0.8 - h_center)**2), 0.5)
    f_body_raw = (body_width * body_height) / (2.0 * np.pi * d_body**2)
    f_body = float(np.clip(f_body_raw, 0.005, 0.10))

    # Renormalise wall factors so total (walls + body) ≈ 1.0
    f_walls_total = (
        float(np.clip(f_floor_raw, 0.01, 0.6))
        + float(np.clip(f_lower_walls, 0.01, 0.5))
        + float(np.clip(f_upper_walls, 0.01, 0.5))
        + float(np.clip(f_ceiling_raw, 0.01, 0.6))
    )
    scale_walls = (1.0 - f_body) / f_walls_total if f_walls_total > 0 else 1.0

    return {
        "floor": float(np.clip(f_floor_raw, 0.01, 0.6)) * scale_walls,
        "lower_walls": float(np.clip(f_lower_walls, 0.01, 0.5)) * scale_walls,
        "upper_walls": float(np.clip(f_upper_walls, 0.01, 0.5)) * scale_walls,
        "ceiling": float(np.clip(f_ceiling_raw, 0.01, 0.6)) * scale_walls,
        "body": f_body,
    }


def _plume_entrainment(
    q_conv_w: float,
    z: float,
    t_amb: float,
    rho_amb: float = 1.1,
    cp: float = 1005.0,
    g: float = 9.81,
    heater_diameter: float = 0.5,
) -> tuple[float, float]:
    """MTT plume entrainment model.

    Computes plume mass flow rate and temperature at height z above
    the heater center using Zukoski's correlation with Heskestad's
    virtual origin correction for finite-diameter sources.

    Args:
        q_conv_w: Convective heat release rate [W].
        z: Height above heater center [m].
        t_amb: Ambient (lower layer) temperature [K].
        rho_amb: Ambient air density [kg/m3].
        cp: Specific heat [J/(kg*K)].
        g: Gravitational acceleration [m/s2].
        heater_diameter: Heater characteristic diameter [m].

    Returns:
        (m_dot, t_plume): Mass flow rate [kg/s] and plume temperature [K].
    """
    if z <= 0.01 or q_conv_w <= 0:
        return 0.0, t_amb

    # Heskestad virtual origin correction for finite-diameter sources
    # z_0 = -1.02 * D + 0.083 * Q_kw^(2/5)
    # z_0 is typically negative for small heaters, making z_eff > z
    q_kw = q_conv_w / 1000.0
    z_0 = -1.02 * heater_diameter + 0.083 * q_kw ** 0.4
    z_eff = max(0.01, z - z_0)

    # Zukoski plume correlation:
    # m_dot = 0.071 * Q_c^(1/3) * z^(5/3) + 0.0018 * Q_c
    # (Q_c in kW, z in m, m_dot in kg/s)
    m_dot = 0.071 * q_kw ** (1.0 / 3.0) * z_eff ** (5.0 / 3.0) + 0.0018 * q_kw

    # Plume temperature from energy balance: Q = m_dot * cp * (T_plume - T_amb)
    t_plume = t_amb + q_conv_w / (m_dot * cp) if m_dot > 0.001 else t_amb + 100.0

    return m_dot, t_plume


def _evaporation_rate(
    water_mass_kg: float,
    elapsed: float,
    tau_evap: float = 5.0,  # characteristic evaporation time [s]
    stone_temp_k: float = 573.15,  # stone surface temp ~300C
) -> float:
    """Evaporation rate of water on hot stones [kg/s].

    Exponential decay model: m_dot = (water_mass / tau) * exp(-t / tau)
    """
    if water_mass_kg <= 0 or elapsed < 0:
        return 0.0
    return (water_mass_kg / tau_evap) * np.exp(-elapsed / tau_evap)


def _perceived_temperature(t_c: float, rh: float, q_rad_body: float = 0.0) -> float:
    """Skin heat balance equivalent temperature for sauna conditions.

    Replaces the Steadman (1979) outdoor heat-index formula which is only
    valid for 20-50 degC.  This model is physically sound up to ~120 degC.

    Physics:
      q_conv  = h_conv * (T_air - T_skin)          convective heat gain
      q_evap  = Lewis-relation evaporative term     (cooling or condensation heating)
      q_rad   = direct radiative flux from heater   (passed in)
      T_eq    = T_skin + q_total / h_ref            normalised equivalent temperature

    Args:
        t_c: Air dry-bulb temperature [degC].
        rh: Relative humidity [0-1].
        q_rad_body: Direct radiative heat flux on body [W/m2] (from heater).

    Returns:
        Equivalent perceived temperature [degC].
    """
    T_SKIN = 36.0    # mean skin temperature [degC]
    H_CONV = 8.0     # convective HTC [W/(m2*K)]
    H_REF = 10.0     # reference HTC for normalisation [W/(m2*K)]

    # Convective heat gain
    q_conv = H_CONV * (t_c - T_SKIN)

    # Saturation vapour pressures (Magnus / August-Roche-Magnus formula) [Pa]
    p_sat_skin = 610.78 * np.exp(17.27 * T_SKIN / (T_SKIN + 237.3))
    p_sat_air = 610.78 * np.exp(17.27 * t_c / (t_c + 237.3)) if t_c > -40.0 else 610.78
    p_vapor = rh * p_sat_air

    # Evaporative / condensation term via simplified Lewis relation
    # h_e ~ 16.5 * h_conv [W/(m2*kPa)] (Lewis factor for air-water)
    # Skin wettedness w limits actual evaporation (0.4 typical in sauna;
    # full body sweating but not every surface is fully wet).
    W_SKIN = 0.4  # skin wettedness fraction [-]
    Q_EVAP_MAX = 400.0  # physiological max evaporative cooling [W/m2]
    p_sat_skin_kpa = p_sat_skin / 1000.0
    p_vapor_kpa = p_vapor / 1000.0
    if p_vapor > p_sat_skin:
        # Condensation on skin surface → adds heat (wettedness irrelevant)
        q_evap = 16.5 * H_CONV * (p_vapor_kpa - p_sat_skin_kpa)
    else:
        # Evaporative cooling (sweat), limited by skin wettedness and physiology
        q_evap_raw = 16.5 * H_CONV * (p_sat_skin_kpa - p_vapor_kpa)
        q_evap = -min(W_SKIN * q_evap_raw, Q_EVAP_MAX)

    # Total heat flux on body
    q_total = q_conv + q_rad_body + q_evap

    # Equivalent temperature: dry-air temperature at H_REF giving same flux
    return T_SKIN + q_total / H_REF


def _ventilation_flow(
    t_upper: float,
    t_lower: float,
    z_int: float,
    vent_cfg: dict,
    rho_0: float = 1.1,
    g: float = 9.81,
) -> float:
    """Compute natural ventilation mass flow rate via stack effect.

    Stack pressure: dp = rho_amb * g * (z_exhaust - z_supply) * (T_hot - T_amb) / T_amb
    Orifice flow:   m_dot = Cd_eff * A_eff * sqrt(2 * rho * |dp|) * sign(dp)

    Args:
        t_upper: Upper layer temperature [K].
        t_lower: Lower layer temperature [K].
        z_int: Interface height [m].
        vent_cfg: Ventilation config dict with 'supply', 'exhaust', 'T_ambient'.
        rho_0: Reference air density [kg/m3].
        g: Gravitational acceleration [m/s2].

    Returns:
        Mass flow rate [kg/s] (positive = inflow through supply).
    """
    supply = vent_cfg["supply"]
    exhaust = vent_cfg["exhaust"]
    t_ambient = vent_cfg.get("T_ambient", 293.15)

    z_supply = supply["height"]
    z_exhaust = exhaust["height"]
    a_supply = supply["area"]
    a_exhaust = exhaust["area"]
    cd_supply = supply.get("Cd", 0.6)
    cd_exhaust = exhaust.get("Cd", 0.6)

    if z_exhaust <= z_supply:
        return 0.0

    # Temperature at each vent height
    t_at_supply = t_lower if z_supply < z_int else t_upper
    t_at_exhaust = t_upper if z_exhaust > z_int else t_lower

    # Ambient density
    rho_ambient = rho_0 * 300.0 / max(t_ambient, 250.0)

    # Stack effect pressure difference (height-weighted average temperature)
    # dp = rho_amb * g * dz * (T_col - T_amb) / T_col
    # Denominator is T_col (not T_amb) from ideal gas: drho = rho_amb*(1 - T_amb/T_col)
    dz_total = z_exhaust - z_supply
    dz_lower_frac = max(0.0, min(z_int, z_exhaust) - z_supply) / dz_total if dz_total > 0 else 0.5
    dz_upper_frac = 1.0 - dz_lower_frac
    t_avg_col = dz_lower_frac * t_at_supply + dz_upper_frac * t_at_exhaust
    delta_p = rho_ambient * g * dz_total * (t_avg_col - t_ambient) / max(t_avg_col, 250.0)

    # Effective orifice area (balanced flow: limited by smaller vent)
    a_eff = min(cd_supply * a_supply, cd_exhaust * a_exhaust)

    # Upwind density: use the density of the fluid flowing INTO the orifice
    # Inflow (delta_p > 0): ambient air enters → use rho_ambient
    # Outflow (delta_p < 0): indoor air exits → use rho at supply vent
    rho_at_supply = rho_0 * 300.0 / max(t_at_supply, 250.0)
    rho_upwind = rho_ambient if delta_p > 0 else rho_at_supply

    # Orifice mass flow: sgn(dp) * Cd * A * sqrt(2 * rho_upwind * |dp|)
    if abs(delta_p) < 1e-6:
        return 0.0

    m_dot = a_eff * math.sqrt(2.0 * rho_upwind * abs(delta_p))
    if delta_p < 0:
        m_dot = -m_dot

    # Sauna is always hotter than ambient — clamp to non-negative.
    # Reverse flow (cold outdoor air sinking through exhaust) is not modeled.
    return max(m_dot, 0.0)


def _humid_air_properties(
    t_k: float, humidity_ratio: float = 0.0, q_rad_body: float = 0.0,
) -> dict[str, float]:
    """Compute humid air mixture properties.

    Args:
        t_k: Temperature [K].
        humidity_ratio: kg water vapor per kg dry air (absolute humidity).
        q_rad_body: Direct radiative heat flux on body [W/m2] (from heater).

    Returns:
        dict with cp_mix, lambda_mix, h_wall_eff, perceived_temp_c, etc.
    """
    # Dry air properties
    cp_air = 1005.0  # J/(kg*K)
    lambda_air = 0.026  # W/(m*K) at ~350K

    # Water vapor properties
    cp_vapor = 1860.0  # J/(kg*K)
    lambda_vapor = 0.025  # W/(m*K) at ~370K (slightly less than air)

    # Mixture properties (mass-weighted)
    w = max(0.0, humidity_ratio)
    y_vapor = w / (1.0 + w)  # vapor mass fraction
    cp_mix = cp_air * (1.0 - y_vapor) + cp_vapor * y_vapor
    lambda_mix = lambda_air * (1.0 - y_vapor) + lambda_vapor * y_vapor

    # Effective h_wall: humid air has higher thermal capacity → stronger
    # natural convection (Ra ∝ β * ΔT * cp * ρ² / (μ * λ))
    # Simplified: h scales roughly as (cp_mix / cp_air)^0.25 * (lambda_mix / lambda_air)^0.75
    h_ratio = (cp_mix / cp_air) ** 0.25 * (lambda_mix / lambda_air) ** 0.75
    h_wall_eff = 8.0 * h_ratio

    t_c = t_k - 273.15
    # Antoine equation for saturation pressure [Pa] at T [K]
    if t_k > 273.15:
        p_sat = 610.78 * np.exp(17.27 * t_c / (t_c + 237.3))
    else:
        p_sat = 610.78
    # Relative humidity approximation (from humidity ratio)
    p_atm = 101325.0
    p_vapor = w * p_atm / (0.622 + w) if w > 0 else 0.0
    rh = min(p_vapor / p_sat, 1.0) if p_sat > 0 else 0.0

    # Skin heat balance perceived temperature (valid for sauna range 60-120 degC)
    perceived_c = _perceived_temperature(t_c, rh, q_rad_body)

    return {
        "cp_mix": cp_mix,
        "lambda_mix": lambda_mix,
        "h_wall_eff": h_wall_eff,
        "humidity_ratio": w,
        "relative_humidity": rh,
        "perceived_temp_c": perceived_c,
    }


def _q_rad_body(
    power_w: float,
    f_conv: float,
    f_body: float,
    heater_w: float,
    heater_h: float,
    t_wall_inner_k: float,
) -> float:
    """Radiative heat flux from heater directly to a person's body [W/m2].

    Estimates the heater surface temperature from its radiant power output
    and the heater area, then computes the net radiative flux intercepted
    by the body using the heater-to-body view factor.

    Args:
        power_w: Total heater power [W].
        f_conv: Convective fraction of heater output (rest is radiant).
        f_body: View factor from heater to body (from _compute_view_factors).
        heater_w: Heater width [m].
        heater_h: Heater height [m].
        t_wall_inner_k: Inner wall surface temperature [K] (used as ambient
            radiative reference for heater surface temperature estimate).

    Returns:
        Radiative heat flux on the body [W/m2].
    """
    SIGMA = 5.67e-8        # Stefan-Boltzmann constant [W/(m2*K4)]
    EPSILON_HEATER = 0.90  # heater surface emissivity (stone/metal)
    EPSILON_BODY = 0.97    # human skin emissivity
    T_SKIN_K = 36.0 + 273.15

    a_heater = max(heater_w * heater_h, 0.01)
    q_rad_total = power_w * (1.0 - f_conv)  # total radiant output [W]
    q_heater_surface = q_rad_total / a_heater  # [W/m2]

    # Estimate heater surface temperature from heater's own emissivity
    t_heater_k = (q_heater_surface / (EPSILON_HEATER * SIGMA) + t_wall_inner_k**4) ** 0.25

    # Net radiative heat flux intercepted by the body
    return float(
        EPSILON_BODY * SIGMA * f_body * (t_heater_k**4 - T_SKIN_K**4)
    )


def _build_profile_and_probes(
    *,
    data: dict,
    height: float,
    heater_y: float,
    heater_h: float,
    n_profile: int,
    z_int: float,
    t_lower: float,
    t_upper: float,
    t_plume: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Build a smooth vertical profile and sample configured probes."""
    dy = height / n_profile
    y = np.linspace(dy / 2, height - dy / 2, n_profile)
    temperatures = np.zeros(n_profile)

    transition_width = 0.15 * height
    for i in range(n_profile):
        sigma = (y[i] - z_int) / (transition_width / 2.0)
        blend = 1.0 / (1.0 + np.exp(-sigma * 3.0))
        temperatures[i] = t_lower + blend * (t_upper - t_lower)
        if heater_y <= y[i] <= heater_y + heater_h:
            plume_boost = (t_plume - temperatures[i]) * 0.1
            temperatures[i] += max(0, plume_boost)

    probe_values = {}
    for probe in data.get("probes", []):
        py = probe["position"]["y"]
        idx = int(np.clip(py / dy, 0, n_profile - 1))
        probe_values[probe["name"]] = float(temperatures[idx])

    return y, temperatures, probe_values


def _make_simple_result(
    *,
    data: dict,
    n_profile: int,
    height: float,
    heater_y: float,
    heater_h: float,
    z_int: float,
    t_upper: float,
    t_lower: float,
    t_plume: float,
    m_plume: float,
    iterations: int,
    converged: bool,
    residual_history: list[float],
    peak_steam_rate: float,
    total_steam: float,
    beta_aug_applied: float,
    t_wall_inner: float,
    humidity_ratio: float,
    q_rad_body: float = 0.0,
    vent_mass_flow: float = 0.0,
) -> SimpleSolverResult:
    """Assemble a ``SimpleSolverResult`` from the current zone state."""
    y, temperatures, probe_values = _build_profile_and_probes(
        data=data,
        height=height,
        heater_y=heater_y,
        heater_h=heater_h,
        n_profile=n_profile,
        z_int=z_int,
        t_lower=t_lower,
        t_upper=t_upper,
        t_plume=t_plume,
    )
    props_upper = _humid_air_properties(t_upper, humidity_ratio, q_rad_body)
    props_lower = _humid_air_properties(t_lower, humidity_ratio * 0.3)

    return SimpleSolverResult(
        y_positions=y,
        temperatures=temperatures,
        probe_values=probe_values,
        iterations=iterations,
        converged=converged,
        residual_history=residual_history,
        interface_height=float(z_int),
        upper_layer_temp=float(t_upper),
        lower_layer_temp=float(t_lower),
        plume_mass_flow=float(m_plume),
        steam_mass_flow=float(peak_steam_rate),
        total_steam_generated=float(total_steam),
        beta_aug_applied=float(beta_aug_applied),
        vent_mass_flow=float(vent_mass_flow),
        wall_inner_temp=float(t_wall_inner),
        humidity_ratio=float(humidity_ratio),
        relative_humidity=float(props_upper["relative_humidity"]),
        perceived_temp_upper=float(props_upper["perceived_temp_c"]),
        perceived_temp_lower=float(props_lower["perceived_temp_c"]),
    )


def solve_two_zone(
    case_yaml: Path,
    n_profile: int = 80,
    max_iter: int = 10000,
    dt: float = 0.5,
    tol: float = 1e-4,
) -> SimpleSolverResult:
    """Solve sauna thermal stratification using a two-zone plume model.

    The model tracks two state variables:
      - T_upper: temperature of the upper (hot) layer [K]
      - z_int: height of the interface between layers [m]

    The lower layer temperature is determined by wall heat exchange.

    Energy balance for upper layer:
      rho * cp * V_upper * dT_upper/dt = Q_plume_in - Q_wall_upper

    Mass balance (interface movement):
      A_floor * dz_int/dt = -(m_plume / rho_upper) + (m_leak / rho_lower)

    where m_leak represents the slow leakage from upper to lower
    layer through the descending interface (wall-driven return flow).
    """
    data = load_yaml(case_yaml)
    dims = data["geometry"]["dimensions"]
    bc = data["boundary_conditions"]
    heater = bc.get("heater", {})
    walls = bc.get("walls", {})

    height = dims["y"]
    width = dims["x"]
    depth = dims["z"]
    a_floor = width * depth  # floor area [m2]

    power_w = heater.get("power_kw", 9.0) * 1000.0
    t_wall = walls.get("temperature", 293.15)
    heater_y = heater.get("position", {}).get("y", 0.0)
    heater_h = heater.get("height", 0.5)
    heater_center_y = heater_y + heater_h / 2.0

    heater_w = heater.get("width", 0.6)

    # View factor radiation model
    vf = _compute_view_factors(width, depth, height, heater_y, heater_h, heater_w)
    f_rad_lower = vf["floor"] + vf["lower_walls"]
    f_rad_upper = vf["ceiling"] + vf["upper_walls"]
    f_body = vf["body"]

    # Löyly (steam) parameters
    loyly = data.get("loyly")
    if loyly:
        water_kg = loyly.get("water_ml", 0) / 1000.0
        loyly_time = loyly.get("time", 0.0)
        tau_evap = loyly.get("tau_evap", 5.0)
    else:
        water_kg = 0.0
        loyly_time = 0.0
        tau_evap = 5.0

    # Aufguss (forced mixing) parameters
    aufguss = data.get("aufguss")
    if aufguss:
        beta_aug = aufguss.get("beta_aug", 0.0)
        aufguss_start = aufguss.get("start_time", 0.0)
        aufguss_duration = aufguss.get("duration", 30.0)
    else:
        beta_aug = 0.0
        aufguss_start = 0.0
        aufguss_duration = 0.0


    # Ventilation parameters
    vent_cfg = data.get("ventilation")
    vent_enabled = vent_cfg is not None and vent_cfg.get("model", "none") != "none"

    # Convective fraction of heater output (rest is radiant to walls)
    f_conv = 0.7
    q_conv = power_w * f_conv

    # Air properties
    rho_0 = 1.1  # reference density at ~300K [kg/m3]
    cp = 1005.0

    # Wall thermal model
    wall_cfg = walls.get("model", "fixed")  # "fixed" or "lumped"
    wall_thickness = walls.get("thickness", 0.015)  # wood panel [m]
    wall_lambda = walls.get("conductivity", 0.12)  # wood thermal conductivity [W/(m*K)]
    wall_rho_cp = walls.get("rho_cp", 0.5e6)  # wood volumetric heat capacity [J/(m3*K)] (rho=450,cp=1100)
    h_wall_base = 8.0  # base natural convection HTC [W/(m2*K)]

    perimeter = 2 * (width + depth)

    # Initial state
    t_wall_inner = t_wall  # inner wall surface temperature (evolves if lumped model)
    t_upper = t_wall + 1.0
    t_lower = t_wall
    z_int = height * 0.95
    humidity_ratio = 0.0  # kg vapor / kg dry air
    steam_remaining = 0.0  # kg of steam still in the room

    # Total wall area for lumped wall model
    a_wall_total = 2 * (width * height + depth * height) + 2 * a_floor
    wall_mass_cp = wall_rho_cp * wall_thickness * a_wall_total  # [J/K]

    residual_history = []
    converged = False
    total_steam = 0.0
    peak_steam_rate = 0.0
    pseudo_time = 0.0
    steam_applied = False

    for iteration in range(max_iter):
        t_upper_old = t_upper
        t_lower_old = t_lower
        z_int_old = z_int
        t_wall_inner_old = t_wall_inner

        # Layer volumes
        v_upper = a_floor * (height - z_int)
        v_lower = a_floor * z_int

        if v_upper < 0.01 * a_floor * height:
            v_upper = 0.01 * a_floor * height
        if v_lower < 0.01 * a_floor * height:
            v_lower = 0.01 * a_floor * height

        # Density variation (ideal gas approximation)
        rho_upper = rho_0 * t_wall / max(t_upper, 250.0)
        rho_lower = rho_0 * t_wall / max(t_lower, 250.0)

        # Plume at interface height
        z_plume = max(0.01, z_int - heater_center_y)
        m_plume, t_plume = _plume_entrainment(
            q_conv, z_plume, t_lower, heater_diameter=heater_w,
        )

        # --- Humidity-dependent properties ---
        props = _humid_air_properties(t_upper, humidity_ratio)
        h_wall = props["h_wall_eff"]
        cp_eff = props["cp_mix"]

        # --- Ventilation flow ---
        if vent_enabled:
            m_vent = _ventilation_flow(t_upper, t_lower, z_int, vent_cfg, rho_0)
            t_ambient_vent = vent_cfg.get("T_ambient", 293.15)
            w_ambient_vent = vent_cfg.get("w_ambient", 0.005)
        else:
            m_vent = 0.0
            t_ambient_vent = t_wall
            w_ambient_vent = 0.0

        # --- Upper layer energy balance ---
        q_plume_in = m_plume * cp_eff * (t_plume - t_upper)
        # Ventilation: exhausted upper-layer air replaced by lower-layer air
        q_vent_upper = m_vent * cp_eff * (t_lower - t_upper) if vent_enabled else 0.0

        upper_height = height - z_int
        a_wall_upper = perimeter * upper_height + a_floor
        q_wall_upper = h_wall * a_wall_upper * (t_upper - t_wall_inner)

        # Radiative loss from heater to walls (non-convective fraction)
        q_rad_to_walls = power_w * (1.0 - f_conv)

        # Radiation to upper zone (fixed wall: direct to air; lumped: via wall model)
        q_rad_to_upper = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_upper

        m_upper = rho_upper * v_upper
        dt_upper = (q_plume_in - q_wall_upper + q_vent_upper + q_rad_to_upper) / (m_upper * cp_eff) if m_upper > 0.1 else 0.0
        t_upper += dt * dt_upper

        # Steam (löyly) — steady-state treatment
        # Latent heat is drawn FROM the heater/stones, not created.
        # Humidity is a state variable: water_kg of steam is distributed
        # in the current upper layer mass, so humidity_ratio tracks m_upper
        # as it evolves. Ventilation removes moisture, reducing the
        # effective steam remaining in the room.
        if water_kg > 0:
            if not steam_applied:
                steam_remaining = water_kg
                total_steam = water_kg
                peak_steam_rate = water_kg / tau_evap
                steam_applied = True
            # Ventilation removes humid air → reduces steam remaining in room
            if vent_enabled and m_vent > 0 and m_upper > 0.1 and steam_remaining > 0:
                dm_steam_out = m_vent * humidity_ratio / (1.0 + humidity_ratio) * dt
                steam_remaining = max(steam_remaining - dm_steam_out, 0.0)
            # Humidity ratio from current steam mass in current upper layer
            humidity_ratio = steam_remaining / max(m_upper, 0.1)
            v_steam = 0.0
            m_dot_steam = 0.0
        else:
            m_dot_steam = 0.0
            v_steam = 0.0

        # --- Lower layer energy balance ---
        # Lower layer receives heat from: wall radiation, conduction from upper
        # Lower layer loses heat to: plume entrainment, floor, lower walls
        a_wall_lower = perimeter * z_int + a_floor
        q_wall_lower = h_wall * a_wall_lower * (t_lower - t_wall_inner)

        # Conductive exchange across interface (small)
        k_interface = 0.5  # effective interface conductivity [W/(m*K)]
        q_interface = k_interface * a_floor * (t_upper - t_lower) / (height * 0.1)

        # Ventilation: fresh ambient air enters lower layer
        q_vent_lower = m_vent * cp_eff * (t_ambient_vent - t_lower) if vent_enabled else 0.0

        # Heat removed by plume entrainment (air leaves lower layer at t_lower)
        # This is already accounted in the plume model

        m_lower = rho_lower * v_lower
        if m_lower > 0.1:
            # Radiation path depends on wall model:
            # - lumped: all radiation heats walls, then conducts/convects to air
            #   → no direct q_rad to air (handled by wall model)
            # - fixed: radiation goes directly to air (no wall dynamics)
            q_rad_to_lower = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_lower
            dt_lower = (q_rad_to_lower + q_interface - q_wall_lower + q_vent_lower) / (m_lower * cp_eff)
        else:
            dt_lower = 0.0
        t_lower += dt * dt_lower

        v_plume_flow = m_plume / max(rho_upper, 0.5)
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall_inner) / (rho_upper * cp_eff * dt_layers)

        # Aufguss mass exchange: bidirectional mixing moves beta_aug kg/s
        # in each direction. Net volume effect from density difference:
        v_mix = beta_aug * (1.0 / max(rho_lower, 0.5) - 1.0 / max(rho_upper, 0.5))

        # Ventilation: supply adds volume to lower layer (pushes interface up)
        v_vent = m_vent / max(rho_lower, 0.5) if vent_enabled else 0.0
        dz_int = (-v_plume_flow + v_return - v_steam + v_mix + v_vent) / a_floor
        z_int += dt * dz_int

        # Clamp
        z_int = np.clip(z_int, 0.05 * height, 0.95 * height)
        t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
        t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # --- Lumped wall temperature evolution ---
        # Wall absorbs heat from air (convection) and heater (radiation)
        # Wall loses heat to outside through conduction
        if wall_cfg == "lumped" and wall_mass_cp > 0:
            q_to_wall = (q_wall_upper + q_wall_lower)
            # ALL heater radiation eventually heats the walls (upper + lower portions)
            q_rad_wall = q_rad_to_walls
            q_out = wall_lambda / wall_thickness * a_wall_total * (t_wall_inner - t_wall)
            dt_wall = (q_to_wall + q_rad_wall - q_out) / wall_mass_cp
            t_wall_inner += dt * dt_wall
            t_wall_inner = np.clip(t_wall_inner, t_wall, t_wall + 200)

        # Aufguss forced mixing (ROM: transfers heat from upper to lower)
        # In steady-state solver, aufguss is always active (no physical time)
        if beta_aug > 0:
            q_mix = beta_aug * cp_eff * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= dt * q_mix / (m_upper * cp_eff)
            if m_lower > 0.1:
                t_lower += dt * q_mix / (m_lower * cp_eff)
            t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
            t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Convergence check (all state variables)
        res_t_upper = abs(t_upper - t_upper_old)
        res_t_lower = abs(t_lower - t_lower_old)
        res_z = abs(z_int - z_int_old)
        res_wall = abs(t_wall_inner - t_wall_inner_old) if wall_cfg == "lumped" else 0.0
        residual = max(res_t_upper, res_t_lower, res_z * 10, res_wall)
        residual_history.append(residual)

        if residual < tol and iteration > 100:
            converged = True
            break
        pseudo_time += dt

    total_steam = min(total_steam, water_kg) if water_kg > 0 else 0.0

    # Direct radiation from heater to body
    q_rad_body_val = _q_rad_body(
        power_w, f_conv, f_body, heater_w, heater_h, t_wall_inner,
    )

    return _make_simple_result(
        data=data,
        n_profile=n_profile,
        height=height,
        heater_y=heater_y,
        heater_h=heater_h,
        z_int=z_int,
        t_upper=t_upper,
        t_lower=t_lower,
        t_plume=t_plume,
        m_plume=m_plume,
        iterations=iteration + 1,
        converged=converged,
        residual_history=residual_history,
        peak_steam_rate=peak_steam_rate,
        total_steam=total_steam,
        beta_aug_applied=beta_aug,
        t_wall_inner=t_wall_inner,
        humidity_ratio=humidity_ratio,
        q_rad_body=q_rad_body_val,
        vent_mass_flow=m_vent if vent_enabled else 0.0,
    )


def solve_transient(
    case_yaml: Path,
    n_profile: int = 80,
    physical_dt: float = 1.0,
    end_time: float = 300.0,
    record_interval: float = 1.0,
) -> TransientResult:
    """Solve sauna thermal stratification as a real-time transient simulation.

    Uses the same two-zone plume physics as ``solve_two_zone`` but steps
    through real (physical) time from 0 to ``end_time`` and records
    time-series data at every ``record_interval`` seconds.

    Args:
        case_yaml: Path to YAML case definition.
        n_profile: Number of vertical profile points for final steady result.
        physical_dt: Real time step size [s].
        end_time: Total simulation time [s].
        record_interval: How often to record a snapshot [s].

    Returns:
        TransientResult with time-series arrays and a final SimpleSolverResult.
    """
    data = load_yaml(case_yaml)
    dims = data["geometry"]["dimensions"]
    bc = data["boundary_conditions"]
    heater = bc.get("heater", {})
    walls = bc.get("walls", {})

    height = dims["y"]
    width = dims["x"]
    depth = dims["z"]
    a_floor = width * depth

    power_w = heater.get("power_kw", 9.0) * 1000.0
    t_wall = walls.get("temperature", 293.15)
    heater_y = heater.get("position", {}).get("y", 0.0)
    heater_h = heater.get("height", 0.5)
    heater_center_y = heater_y + heater_h / 2.0

    heater_w = heater.get("width", 0.6)

    vf = _compute_view_factors(width, depth, height, heater_y, heater_h, heater_w)
    f_rad_lower = vf["floor"] + vf["lower_walls"]
    f_rad_upper = vf["ceiling"] + vf["upper_walls"]
    f_body = vf["body"]

    # Loyly parameters
    loyly = data.get("loyly")
    if loyly:
        water_kg = loyly.get("water_ml", 0) / 1000.0
        loyly_time = loyly.get("time", 0.0)
        tau_evap = loyly.get("tau_evap", 5.0)
    else:
        water_kg = 0.0
        loyly_time = 0.0
        tau_evap = 5.0

    # Aufguss parameters
    aufguss = data.get("aufguss")
    if aufguss:
        beta_aug = aufguss.get("beta_aug", 0.0)
        aufguss_start = aufguss.get("start_time", 0.0)
        aufguss_duration = aufguss.get("duration", 30.0)
    else:
        beta_aug = 0.0
        aufguss_start = 0.0
        aufguss_duration = 0.0


    # Ventilation parameters
    vent_cfg = data.get("ventilation")
    vent_enabled = vent_cfg is not None and vent_cfg.get("model", "none") != "none"

    f_conv = 0.7
    q_conv = power_w * f_conv

    rho_0 = 1.1
    cp = 1005.0

    L_VAPORIZATION = 2.26e6
    MW_STEAM = 18.015e-3
    R_GAS = 8.314
    P_ATM = 101325.0

    wall_cfg = walls.get("model", "fixed")
    wall_thickness = walls.get("thickness", 0.015)
    wall_lambda = walls.get("conductivity", 0.12)
    wall_rho_cp = walls.get("rho_cp", 0.5e6)

    perimeter = 2 * (width + depth)

    # Initial state
    t_wall_inner = t_wall
    t_upper = t_wall + 1.0
    t_lower = t_wall
    z_int = height * 0.95
    humidity_ratio = 0.0

    a_wall_total = 2 * (width * height + depth * height) + 2 * a_floor
    wall_mass_cp = wall_rho_cp * wall_thickness * a_wall_total

    total_steam = 0.0
    peak_steam_rate = 0.0

    # Time-series recording
    n_records = int(end_time / record_interval) + 1
    time_arr = np.zeros(n_records)
    t_upper_arr = np.zeros(n_records)
    t_lower_arr = np.zeros(n_records)
    z_int_arr = np.zeros(n_records)
    humidity_arr = np.zeros(n_records)
    wall_temp_arr = np.zeros(n_records)
    perceived_upper_arr = np.zeros(n_records)

    record_idx = 0
    next_record_time = 0.0

    # Record initial state
    q_rad_body_val = _q_rad_body(power_w, f_conv, f_body, heater_w, heater_h, t_wall_inner)
    props_init = _humid_air_properties(t_upper, humidity_ratio, q_rad_body_val)
    time_arr[0] = 0.0
    t_upper_arr[0] = t_upper
    t_lower_arr[0] = t_lower
    z_int_arr[0] = z_int
    humidity_arr[0] = humidity_ratio
    wall_temp_arr[0] = t_wall_inner
    perceived_upper_arr[0] = props_init["perceived_temp_c"]
    record_idx = 1
    next_record_time = record_interval

    current_time = 0.0
    t_plume = t_upper
    m_plume = 0.0

    while current_time < end_time - 1e-9:
        # Limit dt_step so we never skip a record point
        dt_step = min(physical_dt, end_time - current_time)
        if record_idx < n_records and current_time + dt_step > next_record_time + 1e-9:
            dt_step = max(next_record_time - current_time, 1e-6)
        current_time += dt_step

        # Layer volumes
        v_upper = a_floor * (height - z_int)
        v_lower = a_floor * z_int

        if v_upper < 0.01 * a_floor * height:
            v_upper = 0.01 * a_floor * height
        if v_lower < 0.01 * a_floor * height:
            v_lower = 0.01 * a_floor * height

        # Density variation
        rho_upper = rho_0 * t_wall / max(t_upper, 250.0)
        rho_lower = rho_0 * t_wall / max(t_lower, 250.0)

        # Plume at interface height
        z_plume = max(0.01, z_int - heater_center_y)
        m_plume, t_plume = _plume_entrainment(
            q_conv, z_plume, t_lower, heater_diameter=heater_w,
        )

        # Humidity-dependent properties
        props = _humid_air_properties(t_upper, humidity_ratio)
        h_wall = props["h_wall_eff"]
        cp_eff = props["cp_mix"]

        # Ventilation flow
        if vent_enabled:
            m_vent = _ventilation_flow(t_upper, t_lower, z_int, vent_cfg, rho_0)
            t_ambient_vent = vent_cfg.get("T_ambient", 293.15)
            w_ambient_vent = vent_cfg.get("w_ambient", 0.005)
        else:
            m_vent = 0.0
            t_ambient_vent = t_wall
            w_ambient_vent = 0.0

        # Upper layer energy balance
        q_plume_in = m_plume * cp_eff * (t_plume - t_upper)
        # Ventilation: exhausted upper-layer air replaced by lower-layer air
        q_vent_upper = m_vent * cp_eff * (t_lower - t_upper) if vent_enabled else 0.0

        upper_height = height - z_int
        a_wall_upper = perimeter * upper_height + a_floor
        q_wall_upper = h_wall * a_wall_upper * (t_upper - t_wall_inner)

        q_rad_to_walls = power_w * (1.0 - f_conv)

        # Radiation to upper zone (fixed wall: direct to air; lumped: via wall model)
        q_rad_to_upper = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_upper

        m_upper = rho_upper * v_upper
        dt_upper = (q_plume_in - q_wall_upper + q_vent_upper + q_rad_to_upper) / (m_upper * cp_eff) if m_upper > 0.1 else 0.0
        t_upper += dt_step * dt_upper

        # Steam injection (loyly) — interval-integrated evaporation
        if water_kg > 0 and current_time >= loyly_time and total_steam < water_kg:
            t0 = max(current_time - dt_step, loyly_time) - loyly_time
            t1 = current_time - loyly_time
            # Integral of (water_kg/tau)*exp(-t/tau) from t0 to t1
            # = water_kg * (exp(-t0/tau) - exp(-t1/tau))
            mass_evap = water_kg * (np.exp(-t0 / tau_evap) - np.exp(-t1 / tau_evap))
            mass_evap = min(mass_evap, water_kg - total_steam)  # cap at remaining
            m_dot_steam = mass_evap / dt_step if dt_step > 1e-9 else 0.0
            peak_steam_rate = max(peak_steam_rate, m_dot_steam)
            total_steam += mass_evap

            # Latent heat is drawn from stones/heater, not added to air.
            # Steam adds humidity but not net thermal energy.
            if m_upper > 0.1:
                humidity_ratio += mass_evap / m_upper

            v_steam = m_dot_steam * R_GAS * t_upper / (P_ATM * MW_STEAM)
        else:
            m_dot_steam = 0.0
            v_steam = 0.0

        # Ventilation humidity decay: fresh air dilutes indoor humidity
        if vent_enabled and m_vent > 0 and m_upper > 0.1:
            # Humidity change: supply brings w_ambient, exhaust removes w_upper
            dw = m_vent * (w_ambient_vent - humidity_ratio) / m_upper
            humidity_ratio += dt_step * dw
            humidity_ratio = max(humidity_ratio, 0.0)

        # Lower layer energy balance
        a_wall_lower = perimeter * z_int + a_floor
        q_wall_lower = h_wall * a_wall_lower * (t_lower - t_wall_inner)

        k_interface = 0.5
        q_interface = k_interface * a_floor * (t_upper - t_lower) / (height * 0.1)

        # Ventilation: fresh ambient air enters lower layer
        q_vent_lower = m_vent * cp_eff * (t_ambient_vent - t_lower) if vent_enabled else 0.0

        m_lower = rho_lower * v_lower
        if m_lower > 0.1:
            q_rad_to_lower = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_lower
            dt_lower = (q_rad_to_lower + q_interface - q_wall_lower + q_vent_lower) / (m_lower * cp_eff)
        else:
            dt_lower = 0.0
        t_lower += dt_step * dt_lower

        # Interface movement
        v_plume_flow = m_plume / max(rho_upper, 0.5)
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall_inner) / (rho_upper * cp_eff * dt_layers)

        # Aufguss mass exchange: bidirectional mixing moves beta_aug kg/s
        # in each direction (only active during Aufguss window).
        if beta_aug > 0 and aufguss_start <= current_time <= aufguss_start + aufguss_duration:
            v_mix = beta_aug * (1.0 / max(rho_lower, 0.5) - 1.0 / max(rho_upper, 0.5))
        else:
            v_mix = 0.0

        # Ventilation: supply adds volume to lower layer (pushes interface up)
        v_vent = m_vent / max(rho_lower, 0.5) if vent_enabled else 0.0
        dz_int = (-v_plume_flow + v_return - v_steam + v_mix + v_vent) / a_floor
        z_int += dt_step * dz_int

        # Clamp
        z_int = np.clip(z_int, 0.05 * height, 0.95 * height)
        t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
        t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Lumped wall temperature evolution
        if wall_cfg == "lumped" and wall_mass_cp > 0:
            q_to_wall = q_wall_upper + q_wall_lower
            q_rad_wall = q_rad_to_walls  # all radiation heats walls
            q_out = wall_lambda / wall_thickness * a_wall_total * (t_wall_inner - t_wall)
            dt_wall = (q_to_wall + q_rad_wall - q_out) / wall_mass_cp
            t_wall_inner += dt_step * dt_wall
            t_wall_inner = np.clip(t_wall_inner, t_wall, t_wall + 200)

        # Aufguss forced mixing
        if beta_aug > 0 and aufguss_start <= current_time <= aufguss_start + aufguss_duration:
            q_mix = beta_aug * cp_eff * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= dt_step * q_mix / (m_upper * cp_eff)
            if m_lower > 0.1:
                t_lower += dt_step * q_mix / (m_lower * cp_eff)
            t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
            t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Record snapshot
        if record_idx < n_records and current_time >= next_record_time - 1e-9:
            q_rad_body_val = _q_rad_body(power_w, f_conv, f_body, heater_w, heater_h, t_wall_inner)
            props_snap = _humid_air_properties(t_upper, humidity_ratio, q_rad_body_val)
            time_arr[record_idx] = current_time
            t_upper_arr[record_idx] = t_upper
            t_lower_arr[record_idx] = t_lower
            z_int_arr[record_idx] = z_int
            humidity_arr[record_idx] = humidity_ratio
            wall_temp_arr[record_idx] = t_wall_inner
            perceived_upper_arr[record_idx] = props_snap["perceived_temp_c"]
            record_idx += 1
            next_record_time += record_interval

    q_rad_body_val = _q_rad_body(power_w, f_conv, f_body, heater_w, heater_h, t_wall_inner)
    if record_idx == 0 or abs(time_arr[record_idx - 1] - end_time) > 1e-9:
        if record_idx >= len(time_arr):
            time_arr = np.append(time_arr, end_time)
            t_upper_arr = np.append(t_upper_arr, t_upper)
            t_lower_arr = np.append(t_lower_arr, t_lower)
            z_int_arr = np.append(z_int_arr, z_int)
            humidity_arr = np.append(humidity_arr, humidity_ratio)
            wall_temp_arr = np.append(wall_temp_arr, t_wall_inner)
            perceived_upper_arr = np.append(
                perceived_upper_arr, _humid_air_properties(t_upper, humidity_ratio, q_rad_body_val)["perceived_temp_c"]
            )
            record_idx += 1
        else:
            props_final = _humid_air_properties(t_upper, humidity_ratio, q_rad_body_val)
            time_arr[record_idx] = end_time
            t_upper_arr[record_idx] = t_upper
            t_lower_arr[record_idx] = t_lower
            z_int_arr[record_idx] = z_int
            humidity_arr[record_idx] = humidity_ratio
            wall_temp_arr[record_idx] = t_wall_inner
            perceived_upper_arr[record_idx] = props_final["perceived_temp_c"]
            record_idx += 1

    # Trim arrays to actual recorded length
    time_arr = time_arr[:record_idx]
    t_upper_arr = t_upper_arr[:record_idx]
    t_lower_arr = t_lower_arr[:record_idx]
    z_int_arr = z_int_arr[:record_idx]
    humidity_arr = humidity_arr[:record_idx]
    wall_temp_arr = wall_temp_arr[:record_idx]
    perceived_upper_arr = perceived_upper_arr[:record_idx]

    total_steam = min(total_steam, water_kg) if water_kg > 0 else 0.0
    # Report beta_aug only if Aufguss window actually overlapped with simulation time
    aufguss_end = aufguss_start + aufguss_duration
    aufguss_was_active = beta_aug > 0 and aufguss_start < end_time and aufguss_end > 0
    beta_aug_applied = beta_aug if aufguss_was_active else 0.0

    # Use time-averaged values from the last 25% of the simulation for summary,
    # rather than instantaneous final values which can be noisy.
    n_avg = max(1, record_idx // 4)
    t_upper_avg = float(t_upper_arr[-n_avg:].mean()) if record_idx > 0 else t_upper
    t_lower_avg = float(t_lower_arr[-n_avg:].mean()) if record_idx > 0 else t_lower
    z_int_avg = float(z_int_arr[-n_avg:].mean()) if record_idx > 0 else z_int

    # Convergence: check if T_upper stabilised in the averaging window
    if record_idx > 1:
        t_upper_std = float(t_upper_arr[-n_avg:].std())
        converged = t_upper_std < 0.5  # < 0.5 K standard deviation
    else:
        converged = False

    steady = _make_simple_result(
        data=data,
        n_profile=n_profile,
        height=height,
        heater_y=heater_y,
        heater_h=heater_h,
        z_int=z_int_avg,
        t_upper=t_upper_avg,
        t_lower=t_lower_avg,
        t_plume=t_plume,
        m_plume=m_plume,
        iterations=int(np.ceil(end_time / physical_dt)) if physical_dt > 0 else 0,
        converged=converged,
        residual_history=[],
        peak_steam_rate=peak_steam_rate,
        total_steam=total_steam,
        beta_aug_applied=beta_aug_applied,
        t_wall_inner=t_wall_inner,
        humidity_ratio=humidity_ratio,
        q_rad_body=q_rad_body_val,
        vent_mass_flow=m_vent if vent_enabled else 0.0,
    )

    return TransientResult(
        time=time_arr,
        t_upper_series=t_upper_arr,
        t_lower_series=t_lower_arr,
        z_int_series=z_int_arr,
        humidity_series=humidity_arr,
        wall_temp_series=wall_temp_arr,
        perceived_upper_series=perceived_upper_arr,
        steady_result=steady,
    )


# Backward compatibility alias
solve_1d_thermal = solve_two_zone
