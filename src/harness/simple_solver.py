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

    # Löyly (steam) parameters
    loyly = data.get("loyly")
    if loyly:
        water_kg = loyly.get("water_ml", 0) / 1000.0  # mL -> kg
        loyly_time = loyly.get("time", 0.0)  # when water is poured [s]
        tau_evap = loyly.get("tau_evap", 5.0)
    else:
        water_kg = 0.0
        loyly_time = 0.0
        tau_evap = 5.0

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

    # Wall heat transfer
    h_wall = 8.0  # natural convection HTC [W/(m2*K)]
    perimeter = 2 * (width + depth)

    # Initial state
    t_upper = t_wall + 1.0   # upper layer starts near wall temp
    t_lower = t_wall          # lower layer at wall temp
    z_int = height * 0.95     # interface starts near ceiling

    residual_history = []
    converged = False
    total_steam = 0.0
    peak_steam_rate = 0.0
    pseudo_time = 0.0  # accumulated pseudo-time for evaporation model

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

        # --- Upper layer energy balance ---
        # Heat input from plume
        q_plume_in = m_plume * cp * (t_plume - t_upper)

        # Wall heat loss from upper layer
        # Upper layer contacts: ceiling + portion of side walls + portion of front/back
        upper_height = height - z_int
        a_wall_upper = perimeter * upper_height + a_floor  # sides + ceiling
        q_wall_upper = h_wall * a_wall_upper * (t_upper - t_wall)

        # Radiative loss from heater to walls (non-convective fraction)
        q_rad_to_walls = power_w * (1.0 - f_conv)

        # Upper layer temperature change
        m_upper = rho_upper * v_upper
        dt_upper = (q_plume_in - q_wall_upper) / (m_upper * cp) if m_upper > 0.1 else 0.0
        t_upper += dt * dt_upper

        # Steam injection (löyly)
        if water_kg > 0 and pseudo_time >= loyly_time:
            elapsed = pseudo_time - loyly_time
            m_dot_steam = _evaporation_rate(water_kg, elapsed, tau_evap)
            peak_steam_rate = max(peak_steam_rate, m_dot_steam)
            total_steam += m_dot_steam * dt

            # Latent heat adds energy to upper layer
            q_steam = m_dot_steam * L_VAPORIZATION
            t_upper += dt * q_steam / (m_upper * cp) if m_upper > 0.1 else 0.0

            # Steam volume expansion pushes interface down
            v_steam = m_dot_steam * R_GAS * t_upper / (P_ATM * MW_STEAM)
        else:
            m_dot_steam = 0.0
            v_steam = 0.0

        pseudo_time += dt

        # --- Lower layer energy balance ---
        # Lower layer receives heat from: wall radiation, conduction from upper
        # Lower layer loses heat to: plume entrainment, floor, lower walls
        a_wall_lower = perimeter * z_int + a_floor  # sides + floor
        q_wall_lower = h_wall * a_wall_lower * (t_lower - t_wall)

        # Conductive exchange across interface (small)
        k_interface = 0.5  # effective interface conductivity [W/(m*K)]
        q_interface = k_interface * a_floor * (t_upper - t_lower) / (height * 0.1)

        # Heat removed by plume entrainment (air leaves lower layer at t_lower)
        # This is already accounted in the plume model

        m_lower = rho_lower * v_lower
        if m_lower > 0.1:
            dt_lower = (q_rad_to_walls * 0.3 + q_interface - q_wall_lower) / (m_lower * cp)
        else:
            dt_lower = 0.0
        t_lower += dt * dt_lower

        # --- Interface movement ---
        # Plume deposits mass into upper layer, pulling interface down
        # Return flow (wall-driven) pushes interface back up
        # At steady state these balance
        v_plume_flow = m_plume / max(rho_upper, 0.5)  # volume flow into upper
        dt_layers = max(t_upper - t_lower, 1.0)
        v_return = h_wall * a_wall_upper * (t_upper - t_wall) / (rho_upper * cp * dt_layers)
        dz_int = (-v_plume_flow + v_return - v_steam) / a_floor
        z_int += dt * dz_int

        # Clamp
        z_int = np.clip(z_int, 0.05 * height, 0.95 * height)
        t_upper = np.clip(t_upper, t_wall, t_wall + 200)
        t_lower = np.clip(t_lower, t_wall - 1, t_upper)

        # Convergence check
        res_t = abs(t_upper - t_upper_old)
        res_z = abs(z_int - z_int_old)
        residual = max(res_t, res_z * 10)  # scale z residual
        residual_history.append(residual)

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
    )


# Backward compatibility alias
solve_1d_thermal = solve_two_zone
