"""Unit tests for skin overall heat transfer coefficient module."""

from __future__ import annotations

import math

import pytest

from harness.skin_htc import (
    BODY_PARTS,
    compute_skin_balance,
    h_convective,
    h_evaporative,
    h_radiative,
    humidity_ratio_from_rh,
    humidity_ratio_from_water_addition,
    rh_from_humidity_ratio,
    saturation_pressure_kpa,
    sweep_grid,
)

# ---------------------------------------------------------------------------
# Saturation pressure
# ---------------------------------------------------------------------------


def test_saturation_pressure_known_values():
    # Reference values from steam tables (Magnus is ~0.3% accurate)
    assert saturation_pressure_kpa(20.0) == pytest.approx(2.34, rel=0.02)
    assert saturation_pressure_kpa(36.0) == pytest.approx(5.94, rel=0.02)
    assert saturation_pressure_kpa(100.0) == pytest.approx(101.3, rel=0.05)


# ---------------------------------------------------------------------------
# Convective HTC
# ---------------------------------------------------------------------------


def test_h_conv_natural_only_at_low_velocity():
    # At V=0.05 (the floor), natural convection should dominate
    h = h_convective(v_local=0.0, t_air_c=80.0, t_skin_c=36.0, body_part="average")
    # Mitchell: 2.38 * 44^0.25 ~= 2.38 * 2.58 ~= 6.13
    # de Dear forced at V=0.05: 10.4 * 0.05^0.6 ~= 10.4 * 0.16 ~= 1.7
    # Churchill blend (n=3) ~= 6.2
    assert 5.5 < h < 7.0


def test_h_conv_forced_dominates_at_high_velocity():
    h_low = h_convective(v_local=0.1, t_air_c=80.0, t_skin_c=36.0, body_part="face")
    h_high = h_convective(v_local=2.0, t_air_c=80.0, t_skin_c=36.0, body_part="face")
    assert h_high > h_low
    # Aufguss-level: should reach 15-25 W/(m2K) on face
    assert 15.0 < h_high < 30.0


def test_h_conv_face_higher_than_back():
    v = 1.0
    h_face = h_convective(v, 80.0, 36.0, "face")
    h_back = h_convective(v, 80.0, 36.0, "back")
    assert h_face > h_back


def test_h_conv_monotonic_in_velocity():
    velocities = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0]
    h_values = [h_convective(v, 80.0, 36.0, "face") for v in velocities]
    assert all(h_values[i] < h_values[i + 1] for i in range(len(h_values) - 1))


def test_h_conv_unknown_body_part_raises():
    with pytest.raises(ValueError):
        h_convective(1.0, 80.0, 36.0, body_part="nonsense")


# ---------------------------------------------------------------------------
# Radiative HTC
# ---------------------------------------------------------------------------


def test_h_rad_room_temperature():
    # Room: T_mrt=25C, T_skin=33C → ~6 W/(m²·K) with ε=0.97
    h = h_radiative(25.0, t_skin_c=33.0)
    assert 5.5 < h < 6.5


def test_h_rad_sauna_temperature():
    # Sauna: T_mrt = 90C, T_skin = 36C → ~7-9 W/(m²·K)
    h = h_radiative(90.0, t_skin_c=36.0)
    assert 7.0 < h < 10.0


def test_h_rad_increases_with_t_mrt():
    h_low = h_radiative(40.0, 36.0)
    h_high = h_radiative(120.0, 36.0)
    assert h_high > h_low


# ---------------------------------------------------------------------------
# Evaporative HTC (Lewis)
# ---------------------------------------------------------------------------


def test_h_evap_lewis_relation():
    h_c = 8.0
    h_e = h_evaporative(h_c)
    # h_e = 16.5 * h_c
    assert h_e == pytest.approx(132.0, rel=1e-6)


# ---------------------------------------------------------------------------
# Skin balance integration
# ---------------------------------------------------------------------------


def test_balance_typical_dry_sauna():
    # 90 C, low humidity, no wind
    bal = compute_skin_balance(t_air_c=90.0, rh=0.05, v_local=0.1)
    wb = bal.whole_body
    assert wb is not None
    assert wb.q_conv > 0  # heat into skin
    assert wb.q_rad > 0   # T_mrt = T_air > T_skin
    assert wb.q_evap < 0  # sweat evap dominates at low RH
    # Net should be positive (sauna feels hot)
    assert wb.q_total > 0
    # All body parts present
    for part in BODY_PARTS:
        assert part in bal.parts


def test_balance_humid_sauna_has_higher_q_total_than_dry():
    dry = compute_skin_balance(t_air_c=80.0, rh=0.05, v_local=0.5)
    humid = compute_skin_balance(t_air_c=80.0, rh=0.95, v_local=0.5)
    # High humidity → less sweat evap cooling, possibly condensation gain
    assert humid.whole_body.q_total > dry.whole_body.q_total


def test_balance_aufguss_increases_q_conv():
    calm = compute_skin_balance(t_air_c=90.0, rh=0.3, v_local=0.1)
    aufguss = compute_skin_balance(t_air_c=90.0, rh=0.3, v_local=2.0)
    # Wind boosts convective transfer dramatically
    assert aufguss.whole_body.q_conv > 1.5 * calm.whole_body.q_conv
    # Face is most exposed
    assert aufguss.parts["face"].h_conv > calm.parts["face"].h_conv * 2.0


def test_balance_condensation_when_humid_enough():
    # Steam-saturated air at 80 C → P_air >> P_sat_skin (5.94 kPa)
    bal = compute_skin_balance(t_air_c=80.0, rh=0.5, v_local=0.5)
    # P_sat(80) ~= 47.4 kPa, * 0.5 = 23.7 kPa >> 5.94 → condensation
    assert bal.parts["face"].q_evap > 0
    # And heating, not cooling
    assert bal.whole_body.q_evap > 0


def test_balance_evap_cooling_capped_at_physiological_max():
    # Bone-dry, hot, windy → maximum sweat-evap cooling
    bal = compute_skin_balance(t_air_c=100.0, rh=0.0, v_local=3.0, w_skin=1.0)
    # Should be clamped at -400 W/m² for any single part
    for part in bal.parts.values():
        assert part.q_evap >= -400.0 - 1e-6


def test_balance_zero_dt_overall_u_is_nan():
    bal = compute_skin_balance(t_air_c=36.0, rh=0.3, v_local=0.5)
    # T_air = T_skin → U_overall undefined
    assert math.isnan(bal.whole_body.u_overall)


def test_balance_t_mrt_independent_of_t_air():
    # User-specified T_mrt above T_air (e.g., heater nearby)
    bal_default = compute_skin_balance(t_air_c=80.0, rh=0.3, v_local=0.5)
    bal_hot_walls = compute_skin_balance(t_air_c=80.0, rh=0.3, v_local=0.5, t_mrt_c=120.0)
    assert bal_hot_walls.whole_body.q_rad > bal_default.whole_body.q_rad
    assert bal_hot_walls.whole_body.q_total > bal_default.whole_body.q_total


def test_balance_input_validation():
    with pytest.raises(ValueError):
        compute_skin_balance(t_air_c=80, rh=1.5, v_local=0.5)
    with pytest.raises(ValueError):
        compute_skin_balance(t_air_c=80, rh=0.3, v_local=-0.1)
    with pytest.raises(ValueError):
        compute_skin_balance(t_air_c=80, rh=0.3, v_local=0.5, w_skin=2.0)


def test_u_sensible_equals_h_conv_plus_h_rad():
    bal = compute_skin_balance(t_air_c=80.0, rh=0.3, v_local=1.0, t_mrt_c=85.0)
    for part in bal.parts.values():
        assert part.u_sensible == pytest.approx(part.h_conv + part.h_rad, rel=1e-9)


def test_t_operative_between_t_air_and_t_mrt():
    bal = compute_skin_balance(t_air_c=80.0, rh=0.3, v_local=0.5, t_mrt_c=110.0)
    for part in bal.parts.values():
        assert 80.0 - 1e-6 <= part.t_operative_c <= 110.0 + 1e-6


# ---------------------------------------------------------------------------
# Humidity / water addition helpers
# ---------------------------------------------------------------------------


def test_rh_humidity_ratio_round_trip():
    # Round-trip RH → w → RH
    for t in [25.0, 60.0, 80.0, 100.0]:
        for rh in [0.05, 0.30, 0.60, 0.90]:
            w = humidity_ratio_from_rh(rh, t)
            rh_back = rh_from_humidity_ratio(w, t)
            assert rh_back == pytest.approx(rh, abs=0.01)


def test_humidity_ratio_zero_for_zero_rh():
    assert humidity_ratio_from_rh(0.0, 80.0) == 0.0
    assert rh_from_humidity_ratio(0.0, 80.0) == 0.0


def test_water_addition_one_ladle_in_typical_sauna():
    # 1 ladle = 100 g into a 10 m³ sauna at 80°C
    # m_dry ≈ ρ * V = (101.325e3 / (287.05 * 353.15)) * 10 ≈ 9.99 kg
    # → w_added ≈ 0.1 / 9.99 ≈ 0.01 kg/kg = 10 g/kg
    w = humidity_ratio_from_water_addition(
        water_kg=0.1, sauna_volume_m3=10.0, t_air_c=80.0,
    )
    assert 0.009 < w < 0.011


def test_water_addition_resulting_rh_at_80c():
    # 100 g water in 10 m³ at 80°C → about 3% RH (rough check)
    w = humidity_ratio_from_water_addition(0.1, 10.0, 80.0)
    rh = rh_from_humidity_ratio(w, 80.0)
    assert 0.02 < rh < 0.05


def test_water_addition_two_ladles_double_humidity():
    w1 = humidity_ratio_from_water_addition(0.1, 10.0, 80.0)
    w2 = humidity_ratio_from_water_addition(0.2, 10.0, 80.0)
    assert w2 == pytest.approx(2.0 * w1, rel=1e-9)


def test_water_addition_with_initial_humidity():
    w0 = humidity_ratio_from_rh(0.10, 80.0)  # baseline 10% RH
    w_after = humidity_ratio_from_water_addition(
        0.1, 10.0, 80.0, initial_humidity_ratio=w0,
    )
    assert w_after > w0


def test_water_addition_zero_volume_safe():
    assert humidity_ratio_from_water_addition(0.1, 0.0, 80.0) == 0.0


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------


def test_sweep_grid_shape():
    rows = sweep_grid(
        t_air_c_values=[60.0, 80.0, 100.0],
        rh_values=[0.0, 0.5],
        v_values=[0.1, 1.0, 2.0],
    )
    # 3 * 2 * 3 = 18 rows
    assert len(rows) == 18
    # Each row has whole-body and per-part keys
    sample = rows[0]
    assert "wb_q_total" in sample
    assert "wb_u_overall" in sample
    assert "face_h_conv" in sample
    assert "back_h_conv" in sample


def test_sweep_grid_q_total_increases_with_t_air():
    rows = sweep_grid(
        t_air_c_values=[60.0, 80.0, 100.0],
        rh_values=[0.3],
        v_values=[1.0],
    )
    # Three rows in T_air order
    q_totals = [r["wb_q_total"] for r in rows]
    assert q_totals[0] < q_totals[1] < q_totals[2]
