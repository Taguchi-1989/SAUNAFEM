"""Standalone overall heat transfer coefficient calculator for human skin.

Computes the local heat balance on bare skin under sauna-like conditions
as a pure function of (T_air, RH, V_local, T_mrt). Independent of the
existing CFD/zone solver — intended as a standalone analytical tool to
characterise "perceived heat" U(T, RH, V).

Heat balance (sign convention: q > 0 = heat INTO skin):

    q_total = q_conv + q_rad + q_evap

    q_conv = h_c * (T_air - T_skin)                    [W/m²]
    q_rad  = h_r * (T_mrt - T_skin)                    [W/m²]
    q_evap = h_e * w_skin * (P_air - P_sat,skin)       [W/m²]

Coefficients:

    h_c(V, body_part)  natural+forced convection (Mitchell + de Dear)
    h_r(T_mrt, T_skin) linearised gray-body radiation
    h_e(h_c)           Lewis relation, h_e = 16.5 * h_c   [W/(m²·kPa)]

Overall HTCs reported:
    U_sens     = h_c + h_r                       [sensible only]
    U_overall  = q_total / (T_air - T_skin)      [total, dry-bulb basis]

References:
    Mitchell, D. (1974). Convective heat transfer from man and other animals.
    de Dear, R. et al. (1997). Convective and radiative heat transfer
        coefficients for individual human body segments. Int J Biomet 40.
    ASHRAE Fundamentals (2021), Chapter 9 — Thermal Comfort.
    ISO 9920:2007 — Estimation of thermal insulation and water vapour
        resistance of a clothing ensemble.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# de Dear (1997) forced-convection coefficients for nude body segments:
#   h_forc = a * V^b  [W/(m²·K)],  V in m/s, valid 0.1 <= V <= 4 m/s
_FORCED_COEFFS: dict[str, tuple[float, float]] = {
    "face":    (14.0, 0.61),
    "chest":   (10.4, 0.56),
    "back":    ( 9.1, 0.61),
    "arm":     (13.4, 0.55),
    "thigh":   (10.6, 0.59),
    "calf":    (11.6, 0.59),
    "average": (10.4, 0.60),
}

# Body surface area weights (Du Bois fractions, normalised to sum 1.0
# over the parts modelled here). Used for whole-body aggregation.
_AREA_WEIGHTS: dict[str, float] = {
    "face":    0.05,
    "chest":   0.18,
    "back":    0.18,
    "arm":     0.16,  # both arms combined
    "thigh":   0.20,  # both thighs combined
    "calf":    0.13,  # both calves combined
    # remaining 0.10 goes to head/feet/hands lumped into "average"
}

# Stefan-Boltzmann
_SIGMA = 5.67e-8     # W/(m²·K⁴)
_EPSILON_BODY = 0.97  # human skin emissivity
_LEWIS_FACTOR = 16.5  # K/kPa, Lewis relation for air-water

# Saturation pressure at typical mean skin temperature (used for evap)
# computed at runtime from skin temp via Magnus.

# Skin physiology
_DEFAULT_T_SKIN_C = 36.0
_DEFAULT_W_SKIN = 0.4    # skin wettedness in sauna [-]
_Q_EVAP_MAX = 400.0      # physiological max sweat-evap cooling [W/m²]

# Body parts available for "all parts at once" aggregation
BODY_PARTS = ("face", "chest", "back", "arm", "thigh", "calf", "average")


def saturation_pressure_kpa(t_c: float) -> float:
    """Magnus formula for water-vapour saturation pressure [kPa].

    Valid roughly -40°C to 150°C, accurate to ~0.3% for sauna conditions.
    """
    return 0.61078 * np.exp(17.27 * t_c / (t_c + 237.3))


# Atmospheric pressure used for humidity-ratio / RH conversions [kPa]
_P_ATM_KPA = 101.325
# Specific gas constant of dry air [J/(kg·K)]
_R_DRY_AIR = 287.05


def humidity_ratio_from_rh(rh: float, t_air_c: float) -> float:
    """Humidity ratio [kg vapor / kg dry air] from RH at given air temperature."""
    if rh <= 0:
        return 0.0
    p_sat = saturation_pressure_kpa(t_air_c)
    p_vapor = min(rh, 1.0) * p_sat
    return 0.622 * p_vapor / max(_P_ATM_KPA - p_vapor, 1e-6)


def rh_from_humidity_ratio(humidity_ratio: float, t_air_c: float) -> float:
    """Relative humidity [0-1] from humidity ratio at given air temperature."""
    if humidity_ratio <= 0:
        return 0.0
    p_vapor = humidity_ratio * _P_ATM_KPA / (0.622 + humidity_ratio)
    p_sat = saturation_pressure_kpa(t_air_c)
    return min(p_vapor / p_sat, 1.0) if p_sat > 0 else 0.0


def humidity_ratio_from_water_addition(
    water_kg: float,
    sauna_volume_m3: float,
    t_air_c: float,
    initial_humidity_ratio: float = 0.0,
) -> float:
    """Well-mixed humidity ratio after evaporating ``water_kg`` of water into a
    sauna of volume ``sauna_volume_m3`` at ``t_air_c``.

    Assumes the room is fully mixed and atmospheric pressure stays constant
    (i.e. the small extra vapor mass leaks out via vents on the timescale
    we care about; we just track moisture content).

    1 standard ladle ≈ 100 mL ≈ 0.1 kg water.
    """
    if sauna_volume_m3 <= 0:
        return initial_humidity_ratio
    t_k = t_air_c + 273.15
    rho_dry = (_P_ATM_KPA * 1000.0) / (_R_DRY_AIR * t_k)  # kg/m³
    m_dry = rho_dry * sauna_volume_m3
    if m_dry <= 1e-9:
        return initial_humidity_ratio
    return initial_humidity_ratio + max(water_kg, 0.0) / m_dry


def h_convective(
    v_local: float,
    t_air_c: float,
    t_skin_c: float = _DEFAULT_T_SKIN_C,
    body_part: str = "average",
) -> float:
    """Local convective heat transfer coefficient [W/(m²·K)].

    Combines natural convection (Mitchell 1974) and forced convection
    (de Dear 1997 segment correlations) via Churchill blending (n=3).
    """
    if body_part not in _FORCED_COEFFS:
        raise ValueError(
            f"Unknown body_part {body_part!r}. Use one of {list(_FORCED_COEFFS)}."
        )
    # Natural convection — driven by skin/air temperature difference
    dt = abs(t_skin_c - t_air_c)
    h_nat = 2.38 * dt ** 0.25 if dt > 0.01 else 0.5

    # Forced convection — de Dear's per-segment power law
    a, b = _FORCED_COEFFS[body_part]
    v_eff = max(v_local, 0.05)  # 5 cm/s floor (residual room air motion)
    h_forc = a * v_eff ** b

    # Churchill blend, n=3
    return float((h_nat ** 3 + h_forc ** 3) ** (1.0 / 3.0))


def h_radiative(
    t_mrt_c: float,
    t_skin_c: float = _DEFAULT_T_SKIN_C,
    epsilon: float = _EPSILON_BODY,
) -> float:
    """Linearised radiative heat transfer coefficient [W/(m²·K)].

    Definition: q_rad = h_r * (T_mrt - T_skin)

    Derivation: q_rad = εσ(T_mrt^4 - T_skin^4)
                      = εσ(T_mrt² + T_skin²)(T_mrt + T_skin)(T_mrt - T_skin)
    so h_r = εσ(T_mrt² + T_skin²)(T_mrt + T_skin) with T in Kelvin.

    Returns ~5 W/(m²·K) at 25°C ambient, ~9 W/(m²·K) at 100°C sauna.
    """
    t_mrt_k = t_mrt_c + 273.15
    t_skin_k = t_skin_c + 273.15
    return float(
        epsilon * _SIGMA * (t_mrt_k ** 2 + t_skin_k ** 2) * (t_mrt_k + t_skin_k)
    )


def h_evaporative(h_c: float) -> float:
    """Evaporative heat transfer coefficient [W/(m²·kPa)] via Lewis relation."""
    return _LEWIS_FACTOR * h_c


@dataclass
class BodyPartHTC:
    """Heat balance for a single body segment."""

    name: str
    h_conv: float       # [W/(m²·K)]
    h_rad: float        # [W/(m²·K)]
    h_evap: float       # [W/(m²·kPa)]
    q_conv: float       # [W/m²]
    q_rad: float        # [W/m²]
    q_evap: float       # [W/m²] (positive = into skin = condensation gain;
                        #         negative = sweat cooling)
    q_total: float      # [W/m²]
    u_sensible: float   # h_conv + h_rad   [W/(m²·K)]
    u_overall: float    # q_total / (T_air - T_skin) [W/(m²·K)]
    t_operative_c: float
    t_equivalent_c: float


@dataclass
class SkinHeatBalance:
    """Per-part and whole-body heat balance under given environment."""

    t_air_c: float
    rh: float
    v_local: float
    t_mrt_c: float
    t_skin_c: float
    skin_wettedness: float
    p_vapor_kpa: float
    p_sat_skin_kpa: float
    parts: dict[str, BodyPartHTC] = field(default_factory=dict)
    whole_body: BodyPartHTC | None = None


def _balance_for_part(
    part: str,
    *,
    t_air_c: float,
    rh: float,
    v_local: float,
    t_mrt_c: float,
    t_skin_c: float,
    w_skin: float,
    p_vapor_kpa: float,
    p_sat_skin_kpa: float,
) -> BodyPartHTC:
    h_c = h_convective(v_local, t_air_c, t_skin_c, part)
    h_r = h_radiative(t_mrt_c, t_skin_c)
    h_e = h_evaporative(h_c)

    q_conv = h_c * (t_air_c - t_skin_c)
    q_rad = h_r * (t_mrt_c - t_skin_c)

    # Evap: condensation when P_air > P_sat,skin (no wettedness limit, all
    # incoming vapor condenses); otherwise sweat-evap, limited by wettedness
    # and physiological maximum.
    if p_vapor_kpa > p_sat_skin_kpa:
        q_evap = h_e * (p_vapor_kpa - p_sat_skin_kpa)
    else:
        q_evap_raw = h_e * w_skin * (p_vapor_kpa - p_sat_skin_kpa)  # negative
        q_evap = max(q_evap_raw, -_Q_EVAP_MAX)

    q_total = q_conv + q_rad + q_evap

    u_sens = h_c + h_r
    # Operative temperature (sensible-weighted reference air temp)
    if u_sens > 1e-9:
        t_op = (h_c * t_air_c + h_r * t_mrt_c) / u_sens
        t_eq = t_skin_c + q_total / u_sens
    else:
        t_op = t_air_c
        t_eq = t_skin_c

    # Overall U referred to dry-bulb gradient. Undefined at T_air = T_skin;
    # use NaN as a flag in that edge case so callers can skip it.
    dt_dry = t_air_c - t_skin_c
    u_overall = q_total / dt_dry if abs(dt_dry) > 1e-6 else float("nan")

    return BodyPartHTC(
        name=part,
        h_conv=h_c,
        h_rad=h_r,
        h_evap=h_e,
        q_conv=q_conv,
        q_rad=q_rad,
        q_evap=q_evap,
        q_total=q_total,
        u_sensible=u_sens,
        u_overall=u_overall,
        t_operative_c=t_op,
        t_equivalent_c=t_eq,
    )


def compute_skin_balance(
    t_air_c: float,
    rh: float,
    v_local: float,
    t_mrt_c: float | None = None,
    t_skin_c: float = _DEFAULT_T_SKIN_C,
    w_skin: float = _DEFAULT_W_SKIN,
    parts: tuple[str, ...] = BODY_PARTS,
) -> SkinHeatBalance:
    """Compute the full skin heat balance for all requested body parts.

    Args:
        t_air_c: Dry-bulb air temperature [°C].
        rh: Relative humidity [0-1].
        v_local: Local air velocity at the body surface [m/s].
        t_mrt_c: Mean radiant temperature [°C]. Defaults to t_air_c.
        t_skin_c: Mean skin surface temperature [°C].
        w_skin: Skin wettedness fraction [0-1].
        parts: Body parts to evaluate.

    Returns:
        SkinHeatBalance with per-part results and a whole-body weighted average.
    """
    if not 0.0 <= rh <= 1.0:
        raise ValueError(f"rh must be in [0, 1], got {rh}")
    if v_local < 0:
        raise ValueError(f"v_local must be >= 0, got {v_local}")
    if not 0.0 <= w_skin <= 1.0:
        raise ValueError(f"w_skin must be in [0, 1], got {w_skin}")

    t_mrt_c = t_air_c if t_mrt_c is None else t_mrt_c

    p_sat_air = saturation_pressure_kpa(t_air_c)
    p_sat_skin = saturation_pressure_kpa(t_skin_c)
    p_vapor = rh * p_sat_air

    part_results = {
        part: _balance_for_part(
            part,
            t_air_c=t_air_c,
            rh=rh,
            v_local=v_local,
            t_mrt_c=t_mrt_c,
            t_skin_c=t_skin_c,
            w_skin=w_skin,
            p_vapor_kpa=p_vapor,
            p_sat_skin_kpa=p_sat_skin,
        )
        for part in parts
    }

    # Whole-body weighted aggregate. Use the area weights for parts present;
    # renormalise to whatever subset was requested. "average" is excluded
    # from the area-weighted aggregate because it duplicates other parts.
    weighted_parts = {p: w for p, w in _AREA_WEIGHTS.items() if p in part_results}
    total_w = sum(weighted_parts.values())
    if total_w > 0:
        agg = {
            "h_conv": 0.0, "h_rad": 0.0, "h_evap": 0.0,
            "q_conv": 0.0, "q_rad": 0.0, "q_evap": 0.0,
            "q_total": 0.0, "u_sensible": 0.0,
        }
        for p, w in weighted_parts.items():
            r = part_results[p]
            frac = w / total_w
            for key in agg:
                agg[key] += frac * getattr(r, key)
        # h_rad is part-independent; this still gives the same value
        u_sens_wb = agg["u_sensible"]
        dt_dry = t_air_c - t_skin_c
        u_overall_wb = agg["q_total"] / dt_dry if abs(dt_dry) > 1e-6 else float("nan")
        if u_sens_wb > 1e-9:
            t_op_wb = (agg["h_conv"] * t_air_c + agg["h_rad"] * t_mrt_c) / u_sens_wb
            t_eq_wb = t_skin_c + agg["q_total"] / u_sens_wb
        else:
            t_op_wb = t_air_c
            t_eq_wb = t_skin_c
        whole_body = BodyPartHTC(
            name="whole_body",
            h_conv=agg["h_conv"],
            h_rad=agg["h_rad"],
            h_evap=agg["h_evap"],
            q_conv=agg["q_conv"],
            q_rad=agg["q_rad"],
            q_evap=agg["q_evap"],
            q_total=agg["q_total"],
            u_sensible=u_sens_wb,
            u_overall=u_overall_wb,
            t_operative_c=t_op_wb,
            t_equivalent_c=t_eq_wb,
        )
    else:
        whole_body = None

    return SkinHeatBalance(
        t_air_c=t_air_c,
        rh=rh,
        v_local=v_local,
        t_mrt_c=t_mrt_c,
        t_skin_c=t_skin_c,
        skin_wettedness=w_skin,
        p_vapor_kpa=p_vapor,
        p_sat_skin_kpa=p_sat_skin,
        parts=part_results,
        whole_body=whole_body,
    )


def sweep_grid(
    t_air_c_values: list[float],
    rh_values: list[float],
    v_values: list[float],
    t_mrt_c: float | None = None,
    t_skin_c: float = _DEFAULT_T_SKIN_C,
    w_skin: float = _DEFAULT_W_SKIN,
    parts: tuple[str, ...] = BODY_PARTS,
) -> list[dict]:
    """Return a flat list of records for a (T, RH, V) parameter sweep.

    Each record contains the input conditions plus the whole-body and
    per-part outputs as flat keys, suitable for csv.DictWriter or pandas.
    """
    rows = []
    for t_air in t_air_c_values:
        for rh in rh_values:
            for v in v_values:
                bal = compute_skin_balance(
                    t_air_c=t_air,
                    rh=rh,
                    v_local=v,
                    t_mrt_c=t_mrt_c,
                    t_skin_c=t_skin_c,
                    w_skin=w_skin,
                    parts=parts,
                )
                base = {
                    "t_air_c": t_air,
                    "rh": rh,
                    "v_local": v,
                    "t_mrt_c": bal.t_mrt_c,
                    "p_vapor_kpa": bal.p_vapor_kpa,
                    "p_sat_skin_kpa": bal.p_sat_skin_kpa,
                }
                if bal.whole_body is not None:
                    wb = bal.whole_body
                    base.update({
                        "wb_h_conv": wb.h_conv,
                        "wb_h_rad": wb.h_rad,
                        "wb_q_conv": wb.q_conv,
                        "wb_q_rad": wb.q_rad,
                        "wb_q_evap": wb.q_evap,
                        "wb_q_total": wb.q_total,
                        "wb_u_sensible": wb.u_sensible,
                        "wb_u_overall": wb.u_overall,
                        "wb_t_operative_c": wb.t_operative_c,
                        "wb_t_equivalent_c": wb.t_equivalent_c,
                    })
                for part_name, part in bal.parts.items():
                    base.update({
                        f"{part_name}_h_conv": part.h_conv,
                        f"{part_name}_q_total": part.q_total,
                        f"{part_name}_u_overall": part.u_overall,
                    })
                rows.append(base)
    return rows
