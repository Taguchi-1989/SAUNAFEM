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


def _write_aufguss_yaml(tmp_path: Path, beta_aug: float = 0.5,
                         start_time: float = 0.0,
                         duration: float = 1000.0) -> Path:
    """Helper to create a case YAML with aufguss parameters."""
    data = {
        "case": {"name": "aufguss_test", "description": "test", "type": "steady"},
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
        "aufguss": {
            "beta_aug": beta_aug,
            "start_time": start_time,
            "duration": duration,
        },
        "probes": [
            {"name": "upper_bench", "position": {"x": 1.5, "y": 2.0, "z": 1.25}},
            {"name": "lower_bench", "position": {"x": 1.5, "y": 0.8, "z": 1.25}},
            {"name": "floor_level", "position": {"x": 1.5, "y": 0.1, "z": 1.25}},
        ],
    }
    path = tmp_path / "case.yaml"
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(data, f)
    return path


class TestBetaAug:
    def test_no_aufguss_default(self, tmp_path: Path) -> None:
        """Standard case without aufguss key should have beta_aug_applied == 0.0."""
        path = _write_case_yaml(tmp_path)
        result = solve_two_zone(path, max_iter=10000)
        assert result.beta_aug_applied == 0.0

    def test_aufguss_reduces_stratification(self, tmp_path: Path) -> None:
        """With aufguss, upper-lower temp difference should be smaller than dry case."""
        (tmp_path / "dry").mkdir()
        (tmp_path / "aug").mkdir()
        path_dry = _write_case_yaml(tmp_path / "dry")
        path_aug = _write_aufguss_yaml(tmp_path / "aug", beta_aug=0.5)

        r_dry = solve_two_zone(path_dry, max_iter=10000)
        r_aug = solve_two_zone(path_aug, max_iter=10000)

        strat_dry = r_dry.upper_layer_temp - r_dry.lower_layer_temp
        strat_aug = r_aug.upper_layer_temp - r_aug.lower_layer_temp

        assert strat_aug < strat_dry, (
            f"Aufguss should reduce stratification: dry={strat_dry:.2f}, aug={strat_aug:.2f}"
        )

    def test_aufguss_energy_conservation(self, tmp_path: Path) -> None:
        """With aufguss, total energy should be approximately conserved.

        The mixing transfers heat from upper to lower but does not create energy.
        Compare total thermal energy (m*cp*T) between aufguss and dry cases:
        they should be similar since aufguss only redistributes, not adds, heat.
        """
        (tmp_path / "dry").mkdir()
        (tmp_path / "aug").mkdir()
        path_dry = _write_case_yaml(tmp_path / "dry")
        path_aug = _write_aufguss_yaml(tmp_path / "aug", beta_aug=0.5)

        r_dry = solve_two_zone(path_dry, max_iter=10000)
        r_aug = solve_two_zone(path_aug, max_iter=10000)

        # Approximate total energy as average temperature across the profile
        avg_t_dry = float(r_dry.temperatures.mean())
        avg_t_aug = float(r_aug.temperatures.mean())

        # They should be within a few degrees — mixing redistributes, not creates
        assert abs(avg_t_aug - avg_t_dry) < 10.0, (
            f"Aufguss should conserve energy: dry_avg={avg_t_dry:.2f}, aug_avg={avg_t_aug:.2f}"
        )
