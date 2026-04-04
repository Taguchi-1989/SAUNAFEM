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

    Uses parallel-plate and perpendicular-plate analytical approximations
    (Hottel & Sarofim, 1967) for a rectangular enclosure.

    Returns dict with keys: 'floor', 'lower_walls', 'upper_walls', 'ceiling'.
    All values sum to ~1.0 (enclosure closure rule).
    """
    # Heater centroid height
    h_center = heater_y + heater_h / 2.0

    # Approximate view factors using solid-angle fractions
    # Floor: heater looks down at the floor
    # Use the angle subtended by the floor from heater centroid
    dist_to_floor = h_center
    floor_half_angle = np.arctan(room_w / max(dist_to_floor, 0.01))
    f_floor = floor_half_angle / np.pi  # fraction of hemisphere

    # Ceiling: heater looks up at ceiling
    dist_to_ceiling = room_h - h_center
    ceiling_half_angle = np.arctan(room_w / max(dist_to_ceiling, 0.01))
    f_ceiling = ceiling_half_angle / np.pi

    # Opposite wall: direct view across room
    f_opposite = np.arctan(heater_h / room_w) / np.pi

    # Side walls and remaining surfaces
    f_remaining = max(0.0, 1.0 - f_floor - f_ceiling - f_opposite)

    # Split into lower and upper wall portions based on heater height
    lower_frac = h_center / room_h
    upper_frac = 1.0 - lower_frac

    return {
        "floor": float(np.clip(f_floor, 0.01, 0.5)),
        "lower_walls": float(np.clip(f_remaining * lower_frac + f_opposite * lower_frac, 0.01, 0.5)),
        "upper_walls": float(np.clip(f_remaining * upper_frac + f_opposite * upper_frac, 0.01, 0.5)),
        "ceiling": float(np.clip(f_ceiling, 0.01, 0.5)),
    }


def _plume_entrainment(
    q_conv_w: float,
    z: float,
    t_amb: float,
    rho_amb: float = 1.1,
    cp: float = 1005.0,
    g: float = 9.81,
) -> tuple[float, float]:
    """MTT plume entrainment model.

    Computes plume mass flow rate and temperature at height z above
    the virtual origin using Zukoski's correlation.

    Args:
        q_conv_w: Convective heat release rate [W].
        z: Height above heater center [m].
        t_amb: Ambient (lower layer) temperature [K].
        rho_amb: Ambient air density [kg/m3].
        cp: Specific heat [J/(kg*K)].
        g: Gravitational acceleration [m/s2].

    Returns:
        (m_dot, t_plume): Mass flow rate [kg/s] and plume temperature [K].
    """
    if z <= 0.01 or q_conv_w <= 0:
        return 0.0, t_amb + 100.0

    # Zukoski plume correlation:
    # m_dot = 0.071 * Q_c^(1/3) * z^(5/3) + 0.0018 * Q_c
    # (Q_c in kW, z in m, m_dot in kg/s)
    q_kw = q_conv_w / 1000.0
    m_dot = 0.071 * q_kw ** (1.0 / 3.0) * z ** (5.0 / 3.0) + 0.0018 * q_kw

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
    w = humidity_ratio
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

    for iteration in range(max_iter):
        t_upper_old = t_upper
        z_int_old = z_int

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
        m_plume, t_plume = _plume_entrainment(q_conv, z_plume, t_lower)

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

        # Steam injection (löyly)
        if water_kg > 0 and pseudo_time >= loyly_time:
            elapsed = pseudo_time - loyly_time
            m_dot_steam = _evaporation_rate(water_kg, elapsed, tau_evap)
            peak_steam_rate = max(peak_steam_rate, m_dot_steam)
            total_steam += m_dot_steam * dt

            # Latent heat adds energy to upper layer
            q_steam = m_dot_steam * L_VAPORIZATION
            t_upper += dt * q_steam / (m_upper * cp_eff) if m_upper > 0.1 else 0.0

            # Update humidity ratio (vapor added to upper layer air mass)
            if m_upper > 0.1:
                humidity_ratio += m_dot_steam * dt / m_upper

            # Steam volume expansion pushes interface down
            v_steam = m_dot_steam * R_GAS * t_upper / (P_ATM * MW_STEAM)
        else:
            m_dot_steam = 0.0
            v_steam = 0.0

        pseudo_time += dt

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
            dt_lower = (q_rad_to_walls * f_rad_lower + q_interface - q_wall_lower) / (m_lower * cp_eff)
        else:
            dt_lower = 0.0
        t_lower += dt * dt_lower

        # --- Interface movement ---
        # Plume deposits mass into upper layer, pulling interface down
        # Return flow (wall-driven) pushes interface back up
        # At steady state these balance
        v_plume_flow = m_plume / max(rho_upper, 0.5)  # volume flow into upper
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall_inner) / (rho_upper * cp_eff * dt_layers)
        dz_int = (-v_plume_flow + v_return - v_steam) / a_floor
        z_int += dt * dz_int

        # Clamp
        z_int = np.clip(z_int, 0.05 * height, 0.95 * height)
        t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
        t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # --- Lumped wall temperature evolution ---
        # Wall absorbs heat from air (convection) and heater (radiation)
        # Wall loses heat to outside through conduction
        if wall_cfg == "lumped" and wall_mass_cp > 0:
            q_to_wall = (q_wall_upper + q_wall_lower)  # heat from air to wall
            q_rad_wall = q_rad_to_walls * (1.0 - f_rad_lower)  # radiation to upper walls
            q_out = wall_lambda / wall_thickness * a_wall_total * (t_wall_inner - t_wall)
            dt_wall = (q_to_wall + q_rad_wall - q_out) / wall_mass_cp
            t_wall_inner += dt * dt_wall
            t_wall_inner = np.clip(t_wall_inner, t_wall, t_wall + 200)

        # Aufguss forced mixing (ROM: transfers heat from upper to lower)
        if beta_aug > 0 and aufguss_start <= pseudo_time <= aufguss_start + aufguss_duration:
            q_mix = beta_aug * cp * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= dt * q_mix / (m_upper * cp)
            if m_lower > 0.1:
                t_lower += dt * q_mix / (m_lower * cp)
            t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
            t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Convergence check
        res_t = abs(t_upper - t_upper_old)
        res_z = abs(z_int - z_int_old)
        residual = max(res_t, res_z * 10)  # scale z residual
        residual_history.append(residual)

        pseudo_time += dt

        if residual < tol and iteration > 100:
            converged = True
            break

    # Build vertical temperature profile from zone model
    dy = height / n_profile
    y = np.linspace(dy / 2, height - dy / 2, n_profile)
    temperatures = np.zeros(n_profile)

    # Transition zone around interface (smooth step)
    transition_width = 0.15 * height  # 15% of room height
    for i in range(n_profile):
        # Sigmoid transition between lower and upper layer
        sigma = (y[i] - z_int) / (transition_width / 2.0)
        blend = 1.0 / (1.0 + np.exp(-sigma * 3.0))  # smooth step
        temperatures[i] = t_lower + blend * (t_upper - t_lower)

        # Add plume signature near heater wall (slight local heating)
        if heater_y <= y[i] <= heater_y + heater_h:
            plume_boost = (t_plume - temperatures[i]) * 0.1
            temperatures[i] += max(0, plume_boost)

    # Extract probe values
    probes = data.get("probes", [])
    probe_values = {}
    for p in probes:
        py = p["position"]["y"]
        idx = int(np.clip(py / dy, 0, n_profile - 1))
        probe_values[p["name"]] = float(temperatures[idx])

    # Compute final perceived temperatures
    props_upper = _humid_air_properties(t_upper, humidity_ratio)
    props_lower = _humid_air_properties(t_lower, humidity_ratio * 0.3)

    return SimpleSolverResult(
        y_positions=y,
        temperatures=temperatures,
        probe_values=probe_values,
        iterations=iteration + 1,
        converged=converged,
        residual_history=residual_history,
        interface_height=float(z_int),
        upper_layer_temp=float(t_upper),
        lower_layer_temp=float(t_lower),
        plume_mass_flow=float(m_plume),
        steam_mass_flow=float(peak_steam_rate),
        total_steam_generated=float(total_steam),
        beta_aug_applied=float(beta_aug) if aufguss else 0.0,
        wall_inner_temp=float(t_wall_inner),
        humidity_ratio=float(humidity_ratio),
        relative_humidity=float(props_upper["relative_humidity"]),
        perceived_temp_upper=float(props_upper["perceived_temp_c"]),
        perceived_temp_lower=float(props_lower["perceived_temp_c"]),
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

    n_steps = int(end_time / physical_dt)
    current_time = 0.0

    for _step in range(n_steps):
        current_time += physical_dt

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
        m_plume, t_plume = _plume_entrainment(q_conv, z_plume, t_lower)

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
        t_upper += physical_dt * dt_upper

        # Steam injection (loyly)
        if water_kg > 0 and current_time >= loyly_time:
            elapsed = current_time - loyly_time
            m_dot_steam = _evaporation_rate(water_kg, elapsed, tau_evap)
            peak_steam_rate = max(peak_steam_rate, m_dot_steam)
            total_steam += m_dot_steam * physical_dt

            q_steam = m_dot_steam * L_VAPORIZATION
            t_upper += physical_dt * q_steam / (m_upper * cp_eff) if m_upper > 0.1 else 0.0

            if m_upper > 0.1:
                humidity_ratio += m_dot_steam * physical_dt / m_upper

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
            dt_lower = (q_rad_to_walls * f_rad_lower + q_interface - q_wall_lower) / (m_lower * cp_eff)
        else:
            dt_lower = 0.0
        t_lower += physical_dt * dt_lower

        # Interface movement
        v_plume_flow = m_plume / max(rho_upper, 0.5)
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall_inner) / (rho_upper * cp_eff * dt_layers)
        dz_int = (-v_plume_flow + v_return - v_steam) / a_floor
        z_int += physical_dt * dz_int

        # Clamp
        z_int = np.clip(z_int, 0.05 * height, 0.95 * height)
        t_upper = np.clip(t_upper, t_wall_inner, t_wall_inner + 200)
        t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Lumped wall temperature evolution
        if wall_cfg == "lumped" and wall_mass_cp > 0:
            q_to_wall = q_wall_upper + q_wall_lower
            q_rad_wall = q_rad_to_walls * (1.0 - f_rad_lower)
            q_out = wall_lambda / wall_thickness * a_wall_total * (t_wall_inner - t_wall)
            dt_wall = (q_to_wall + q_rad_wall - q_out) / wall_mass_cp
            t_wall_inner += physical_dt * dt_wall
            t_wall_inner = np.clip(t_wall_inner, t_wall, t_wall + 200)

        # Aufguss forced mixing
        if beta_aug > 0 and aufguss_start <= current_time <= aufguss_start + aufguss_duration:
            q_mix = beta_aug * cp * (t_upper - t_lower)
            if m_upper > 0.1:
                t_upper -= physical_dt * q_mix / (m_upper * cp)
            if m_lower > 0.1:
                t_lower += physical_dt * q_mix / (m_lower * cp)
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

    # Trim arrays to actual recorded length
    time_arr = time_arr[:record_idx]
    t_upper_arr = t_upper_arr[:record_idx]
    t_lower_arr = t_lower_arr[:record_idx]
    z_int_arr = z_int_arr[:record_idx]
    humidity_arr = humidity_arr[:record_idx]
    wall_temp_arr = wall_temp_arr[:record_idx]
    perceived_upper_arr = perceived_upper_arr[:record_idx]

    # Get final steady result using solve_two_zone
    steady = solve_two_zone(case_yaml, n_profile=n_profile)

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
