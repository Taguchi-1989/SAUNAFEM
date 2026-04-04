"""Tests for two-zone plume model solver."""

from __future__ import annotations

from pathlib import Path

import yaml

from harness.simple_solver import _plume_entrainment, solve_two_zone


def _write_case_yaml(tmp_path: Path, **overrides) -> Path:
    """Helper to create a temporary case YAML."""
    data = {
        "case": {"name": "test", "description": "test", "type": "steady"},
        "geometry": {
            "dimensions": {"x": 3.0, "y": 2.5, "z": 2.5},
            "mesh_level": "M0",
        },
        "boundary_conditions": {
            "walls": {"temperature": 293.15, "type": "mixed"},
            "heater": {
                "power_kw": 9.0,
                "position": {"x": 0.0, "y": 0.1, "z": 1.25},
                "width": 0.6,
                "height": 0.5,
            },
        },
        "solver": {
            "name": "buoyantPimpleFoam",
            "end_time": 300,
            "write_interval": 10,
            "delta_t": 0.05,
            "averaging_start": 150,
        },
        "probes": [
            {"name": "upper_bench", "position": {"x": 1.5, "y": 2.0, "z": 1.25}},
            {"name": "lower_bench", "position": {"x": 1.5, "y": 0.8, "z": 1.25}},
            {"name": "floor_level", "position": {"x": 1.5, "y": 0.1, "z": 1.25}},
        ],
    }
    # Apply overrides
    for key, val in overrides.items():
        if key == "power_kw":
            data["boundary_conditions"]["heater"]["power_kw"] = val
        elif key == "t_wall":
            data["boundary_conditions"]["walls"]["temperature"] = val
    path = tmp_path / "case.yaml"
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(data, f)
    return path


class TestPlumeEntrainment:
    def test_zero_height_returns_high_temp(self) -> None:
        m, t = _plume_entrainment(9000.0, 0.0, 293.15)
        assert m == 0.0
        assert t > 350.0

    def test_mass_flow_increases_with_height(self) -> None:
        m1, _ = _plume_entrainment(9000.0, 1.0, 293.15)
        m2, _ = _plume_entrainment(9000.0, 2.0, 293.15)
        assert m2 > m1

    def test_plume_temp_decreases_with_height(self) -> None:
        _, t1 = _plume_entrainment(9000.0, 0.5, 293.15)
        _, t2 = _plume_entrainment(9000.0, 2.0, 293.15)
        assert t1 > t2  # more entrainment dilutes plume

    def test_higher_power_increases_temp(self) -> None:
        _, t1 = _plume_entrainment(5000.0, 1.0, 293.15)
        _, t2 = _plume_entrainment(15000.0, 1.0, 293.15)
        assert t2 > t1

    def test_zero_power(self) -> None:
        m, t = _plume_entrainment(0.0, 1.0, 293.15)
        assert m == 0.0


class TestSolveTwoZone:
    def test_converges(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert result.converged is True

    def test_thermal_stratification(self, tmp_path: Path) -> None:
        """Upper bench must be hotter than lower bench."""
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert result.probe_values["upper_bench"] > result.probe_values["lower_bench"]

    def test_interface_within_room(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert 0 < result.interface_height < 2.5

    def test_upper_layer_hotter_than_lower(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert result.upper_layer_temp > result.lower_layer_temp

    def test_plume_mass_flow_positive(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert result.plume_mass_flow > 0

    def test_higher_power_higher_temps(self, tmp_path: Path) -> None:
        (tmp_path / "lo").mkdir()
        (tmp_path / "hi").mkdir()
        path_lo = _write_case_yaml(tmp_path / "lo", power_kw=5.0)
        path_hi = _write_case_yaml(tmp_path / "hi", power_kw=15.0)
        r_lo = solve_two_zone(path_lo, max_iter=10000)
        r_hi = solve_two_zone(path_hi, max_iter=10000)
        assert r_hi.upper_layer_temp > r_lo.upper_layer_temp

    def test_profile_length_matches_n_profile(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, n_profile=40, max_iter=5000)
        assert len(result.y_positions) == 40
        assert len(result.temperatures) == 40

    def test_all_probes_present(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=5000)
        assert "upper_bench" in result.probe_values
        assert "lower_bench" in result.probe_values
        assert "floor_level" in result.probe_values

    def test_residual_history_populated(self, tmp_path: Path) -> None:
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=5000)
        assert len(result.residual_history) > 0
        # Residuals should generally decrease
        assert result.residual_history[-1] < result.residual_history[0]
