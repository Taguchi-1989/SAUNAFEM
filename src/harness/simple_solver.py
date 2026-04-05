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

    Returns dict with keys: 'floor', 'lower_walls', 'upper_walls', 'ceiling'.
    All values sum to ~1.0 (enclosure closure rule).
    """
    h_center = heater_y + heater_h / 2.0

    # Vertical angle fractions from heater centroid
    # Floor: angle subtended looking DOWN from heater
    dist_floor = max(h_center, 0.01)
    # Ceiling: angle subtended looking UP from heater
    dist_ceil = max(room_h - h_center, 0.01)

    # Solid-angle fraction ≈ atan(characteristic_length / distance) / pi
    # Characteristic length for floor/ceiling is the room diagonal projected
    char_len = np.sqrt(room_w**2 + room_d**2) / 2.0

    f_floor_raw = np.arctan(char_len / dist_floor) / np.pi
    f_ceiling_raw = np.arctan(char_len / dist_ceil) / np.pi

    # Opposite wall fraction (across the room width)
    f_opposite = np.arctan(heater_h / max(room_w, 0.01)) / np.pi

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

    return {
        "floor": float(np.clip(f_floor_raw, 0.01, 0.6)),
        "lower_walls": float(np.clip(f_lower_walls, 0.01, 0.5)),
        "upper_walls": float(np.clip(f_upper_walls, 0.01, 0.5)),
        "ceiling": float(np.clip(f_ceiling_raw, 0.01, 0.6)),
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


def _humid_air_properties(t_k: float, humidity_ratio: float = 0.0) -> dict[str, float]:
    """Compute humid air mixture properties.

    Args:
        t_k: Temperature [K].
        humidity_ratio: kg water vapor per kg dry air (absolute humidity).

    Returns:
        dict with cp_mix, lambda_mix, h_wall_eff, and perceived_temp.
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

    # Perceived temperature (simplified wet-bulb approximation)
    # Humid air feels hotter because it impedes evaporative cooling from skin
    # WBGT-like index: T_perceived ≈ T_dry + humidity_effect
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

    # Simplified perceived temperature (Steadman 1979 heat index approximation)
    if t_c > 27 and rh > 0.05:
        perceived_c = t_c + 0.33 * (rh * p_sat / 1000.0) - 4.0
    else:
        perceived_c = t_c

    return {
        "cp_mix": cp_mix,
        "lambda_mix": lambda_mix,
        "h_wall_eff": h_wall_eff,
        "humidity_ratio": w,
        "relative_humidity": rh,
        "perceived_temp_c": perceived_c,
    }


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
    props_upper = _humid_air_properties(t_upper, humidity_ratio)
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

    # View factor radiation model
    vf = _compute_view_factors(width, depth, height, heater_y, heater_h,
                                heater.get("width", 0.6))
    f_rad_lower = vf["floor"] + vf["lower_walls"]

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

    # Heater characteristic diameter for virtual origin correction
    heater_width = heater.get("width", 0.6)

    # Convective fraction of heater output (rest is radiant to walls)
    f_conv = 0.7
    q_conv = power_w * f_conv

    # Air properties
    rho_0 = 1.1  # reference density at ~300K [kg/m3]
    cp = 1005.0

    # Steam properties
    L_VAPORIZATION = 2.26e6  # latent heat of water [J/kg]
    MW_STEAM = 18.015e-3     # molecular weight of steam [kg/mol]
    R_GAS = 8.314            # universal gas constant [J/(mol*K)]
    P_ATM = 101325.0         # atmospheric pressure [Pa]

    # Wall thermal model
    wall_cfg = walls.get("model", "fixed")  # "fixed" or "lumped"
    wall_thickness = walls.get("thickness", 0.02)  # wood panel [m]
    wall_lambda = walls.get("conductivity", 0.12)  # wood thermal conductivity [W/(m*K)]
    wall_rho_cp = walls.get("rho_cp", 0.4e6)  # wood volumetric heat capacity [J/(m3*K)]
    h_wall_base = 8.0  # base natural convection HTC [W/(m2*K)]

    perimeter = 2 * (width + depth)

    # Initial state
    t_wall_inner = t_wall  # inner wall surface temperature (evolves if lumped model)
    t_upper = t_wall + 1.0
    t_lower = t_wall
    z_int = height * 0.95
    humidity_ratio = 0.0  # kg vapor / kg dry air

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
            q_conv, z_plume, t_lower, heater_diameter=heater_width,
        )

        # --- Humidity-dependent properties ---
        props = _humid_air_properties(t_upper, humidity_ratio)
        h_wall = props["h_wall_eff"]
        cp_eff = props["cp_mix"]

        # --- Upper layer energy balance ---
        q_plume_in = m_plume * cp_eff * (t_plume - t_upper)

        upper_height = height - z_int
        a_wall_upper = perimeter * upper_height + a_floor
        q_wall_upper = h_wall * a_wall_upper * (t_upper - t_wall_inner)

        # Radiative loss from heater to walls (non-convective fraction)
        q_rad_to_walls = power_w * (1.0 - f_conv)

        m_upper = rho_upper * v_upper
        dt_upper = (q_plume_in - q_wall_upper) / (m_upper * cp_eff) if m_upper > 0.1 else 0.0
        t_upper += dt * dt_upper

        # Steam (löyly) — steady-state treatment
        # Latent heat is drawn FROM the heater/stones, not created.
        # The evaporation reduces the effective convective power.
        # Humidity is set from total water mass (final equilibrium state).
        if water_kg > 0:
            if not steam_applied:
                humidity_ratio = water_kg / max(m_upper, 0.1)
                total_steam = water_kg
                peak_steam_rate = water_kg / tau_evap
                steam_applied = True
                # Note: no t_upper boost — latent heat comes from stones,
                # not created from nothing. The steam adds humidity (which
                # affects h_wall, cp) but does not add net energy to the air.
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

        # Heat removed by plume entrainment (air leaves lower layer at t_lower)
        # This is already accounted in the plume model

        m_lower = rho_lower * v_lower
        if m_lower > 0.1:
            # Radiation path depends on wall model:
            # - lumped: all radiation heats walls, then conducts/convects to air
            #   → no direct q_rad to air (handled by wall model)
            # - fixed: radiation goes directly to air (no wall dynamics)
            q_rad_to_lower = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_lower
            dt_lower = (q_rad_to_lower + q_interface - q_wall_lower) / (m_lower * cp_eff)
        else:
            dt_lower = 0.0
        t_lower += dt * dt_lower

        v_plume_flow = m_plume / max(rho_upper, 0.5)
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall_inner) / (rho_upper * cp_eff * dt_layers)

        # Aufguss mass exchange: bidirectional mixing moves beta_aug kg/s
        # in each direction. Net volume effect from density difference:
        v_mix = beta_aug * (1.0 / max(rho_lower, 0.5) - 1.0 / max(rho_upper, 0.5))

        dz_int = (-v_plume_flow + v_return - v_steam + v_mix) / a_floor
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
            q_mix = beta_aug * cp * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= dt * q_mix / (m_upper * cp)
            if m_lower > 0.1:
                t_lower += dt * q_mix / (m_lower * cp)
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

    vf = _compute_view_factors(width, depth, height, heater_y, heater_h,
                                heater.get("width", 0.6))
    f_rad_lower = vf["floor"] + vf["lower_walls"]

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

    # Heater characteristic diameter for virtual origin correction
    heater_width = heater.get("width", 0.6)

    f_conv = 0.7
    q_conv = power_w * f_conv

    rho_0 = 1.1
    cp = 1005.0

    L_VAPORIZATION = 2.26e6
    MW_STEAM = 18.015e-3
    R_GAS = 8.314
    P_ATM = 101325.0

    wall_cfg = walls.get("model", "fixed")
    wall_thickness = walls.get("thickness", 0.02)
    wall_lambda = walls.get("conductivity", 0.12)
    wall_rho_cp = walls.get("rho_cp", 0.4e6)

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
    props_init = _humid_air_properties(t_upper, humidity_ratio)
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
            q_conv, z_plume, t_lower, heater_diameter=heater_width,
        )

        # Humidity-dependent properties
        props = _humid_air_properties(t_upper, humidity_ratio)
        h_wall = props["h_wall_eff"]
        cp_eff = props["cp_mix"]

        # Upper layer energy balance
        q_plume_in = m_plume * cp_eff * (t_plume - t_upper)

        upper_height = height - z_int
        a_wall_upper = perimeter * upper_height + a_floor
        q_wall_upper = h_wall * a_wall_upper * (t_upper - t_wall_inner)

        q_rad_to_walls = power_w * (1.0 - f_conv)

        m_upper = rho_upper * v_upper
        dt_upper = (q_plume_in - q_wall_upper) / (m_upper * cp_eff) if m_upper > 0.1 else 0.0
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

        # Lower layer energy balance
        a_wall_lower = perimeter * z_int + a_floor
        q_wall_lower = h_wall * a_wall_lower * (t_lower - t_wall_inner)

        k_interface = 0.5
        q_interface = k_interface * a_floor * (t_upper - t_lower) / (height * 0.1)

        m_lower = rho_lower * v_lower
        if m_lower > 0.1:
            q_rad_to_lower = 0.0 if wall_cfg == "lumped" else q_rad_to_walls * f_rad_lower
            dt_lower = (q_rad_to_lower + q_interface - q_wall_lower) / (m_lower * cp_eff)
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

        dz_int = (-v_plume_flow + v_return - v_steam + v_mix) / a_floor
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
            q_mix = beta_aug * cp * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= dt_step * q_mix / (m_upper * cp)
            if m_lower > 0.1:
                t_lower += dt_step * q_mix / (m_lower * cp)
            t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
            t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Record snapshot
        if record_idx < n_records and current_time >= next_record_time - 1e-9:
            props_snap = _humid_air_properties(t_upper, humidity_ratio)
            time_arr[record_idx] = current_time
            t_upper_arr[record_idx] = t_upper
            t_lower_arr[record_idx] = t_lower
            z_int_arr[record_idx] = z_int
            humidity_arr[record_idx] = humidity_ratio
            wall_temp_arr[record_idx] = t_wall_inner
            perceived_upper_arr[record_idx] = props_snap["perceived_temp_c"]
            record_idx += 1
            next_record_time += record_interval

    if record_idx == 0 or abs(time_arr[record_idx - 1] - end_time) > 1e-9:
        if record_idx >= len(time_arr):
            time_arr = np.append(time_arr, end_time)
            t_upper_arr = np.append(t_upper_arr, t_upper)
            t_lower_arr = np.append(t_lower_arr, t_lower)
            z_int_arr = np.append(z_int_arr, z_int)
            humidity_arr = np.append(humidity_arr, humidity_ratio)
            wall_temp_arr = np.append(wall_temp_arr, t_wall_inner)
            perceived_upper_arr = np.append(
                perceived_upper_arr, _humid_air_properties(t_upper, humidity_ratio)["perceived_temp_c"]
            )
            record_idx += 1
        else:
            props_final = _humid_air_properties(t_upper, humidity_ratio)
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
    # Report beta_aug if aufguss was configured (regardless of current time)
    beta_aug_applied = beta_aug if beta_aug > 0 else 0.0
    steady = _make_simple_result(
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
        iterations=int(np.ceil(end_time / physical_dt)) if physical_dt > 0 else 0,
        converged=abs(t_upper_arr[-1] - t_upper_arr[-2]) < 1e-3 if len(t_upper_arr) > 1 else False,
        residual_history=[],
        peak_steam_rate=peak_steam_rate,
        total_steam=total_steam,
        beta_aug_applied=beta_aug_applied,
        t_wall_inner=t_wall_inner,
        humidity_ratio=humidity_ratio,
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
