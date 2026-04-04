"""Tests for case_builder module."""

from __future__ import annotations

from pathlib import Path

from harness.case_builder import (
    build_case,
    compute_heater_params,
    compute_mesh_params,
)


class TestComputeMeshParams:
    def test_m0_level(self) -> None:
        geo = {"dimensions": {"x": 3.0, "y": 2.5, "z": 2.5}, "mesh_level": "M0"}
        result = compute_mesh_params(geo)
        assert result["dim_x"] == 3.0
        assert result["nx"] == 24  # 3.0 * 8
        assert result["ny"] == 20  # 2.5 * 8
        assert result["nz"] == 20

    def test_m1_level(self) -> None:
        geo = {"dimensions": {"x": 3.0, "y": 2.5, "z": 2.5}, "mesh_level": "M1"}
        result = compute_mesh_params(geo)
        assert result["nx"] == 48  # 3.0 * 16

    def test_default_level(self) -> None:
        geo = {"dimensions": {"x": 1.0, "y": 1.0, "z": 1.0}}
        result = compute_mesh_params(geo)
        # Default is M0 density=8
        assert result["nx"] == 8

    def test_minimum_cells(self) -> None:
        geo = {"dimensions": {"x": 0.1, "y": 0.1, "z": 0.1}, "mesh_level": "M0"}
        result = compute_mesh_params(geo)
        # Should be at least 4 cells per axis
        assert result["nx"] >= 4
        assert result["ny"] >= 4


class TestComputeHeaterParams:
    def test_heat_flux_calculation(self) -> None:
        bc = {"heater": {"power_kw": 9.0, "width": 0.6, "height": 0.5},
              "walls": {"temperature": 293.15}}
        geo = {"dimensions": {"x": 3.0, "y": 2.5, "z": 2.5}}
        result = compute_heater_params(bc, geo)
        # 9000 W / (0.6 * 0.5) = 30000 W/m2
        assert result["heat_flux"] == 30000.0
        assert result["T_walls"] == 293.15
        assert result["T_initial"] == 293.15

    def test_default_wall_temperature(self) -> None:
        bc = {"heater": {"power_kw": 5.0}}
        geo = {"dimensions": {"x": 2.0, "y": 2.0, "z": 2.0}}
        result = compute_heater_params(bc, geo)
        assert result["T_walls"] == 293.15


class TestBuildCase:
    def test_creates_case_directory(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        result = build_case(sample_case_path, output_dir=out)
        assert result == out
        assert out.exists()

    def test_creates_openfoam_structure(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        build_case(sample_case_path, output_dir=out)
        # Check key directories exist
        assert (out / "0").is_dir()
        assert (out / "constant").is_dir()
        assert (out / "system").is_dir()

    def test_creates_key_files(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        build_case(sample_case_path, output_dir=out)
        assert (out / "system" / "blockMeshDict").is_file()
        assert (out / "system" / "controlDict").is_file()
        assert (out / "system" / "fvSchemes").is_file()
        assert (out / "system" / "fvSolution").is_file()
        assert (out / "0" / "T").is_file()
        assert (out / "0" / "U").is_file()
        assert (out / "0" / "p_rgh").is_file()
        assert (out / "constant" / "g").is_file()

    def test_template_rendering_content(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        build_case(sample_case_path, output_dir=out)

        # Check blockMeshDict has correct dimensions
        mesh_dict = (out / "system" / "blockMeshDict").read_text(encoding="utf-8")
        assert "3.0" in mesh_dict  # dim_x
        assert "2.5" in mesh_dict  # dim_y
        assert "heater_wall_surround" in mesh_dict
        assert "0.1" in mesh_dict
        assert "0.6" in mesh_dict

        # Check controlDict has correct solver
        control = (out / "system" / "controlDict").read_text(encoding="utf-8")
        assert "buoyantPimpleFoam" in control

        # Check T file has heat flux
        t_file = (out / "0" / "T").read_text(encoding="utf-8")
        assert "externalWallHeatFluxTemperature" in t_file
        assert "30000.0" in t_file  # 9000/(0.6*0.5)
        assert "heater_wall_surround" in t_file

    def test_probes_in_controldict(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        build_case(sample_case_path, output_dir=out)
        control = (out / "system" / "controlDict").read_text(encoding="utf-8")
        assert "upper_bench" in control
        assert "lower_bench" in control
        assert "1.5 2.0 1.25" in control

    def test_invalid_yaml_raises(self, tmp_path: Path) -> None:
        bad_yaml = tmp_path / "bad.yaml"
        bad_yaml.write_text("case:\n  name: test\n", encoding="utf-8")
        import pytest

        with pytest.raises(ValueError, match="Schema validation failed"):
            build_case(bad_yaml, output_dir=tmp_path / "out")

    def test_oversized_heater_raises(self, sample_case_path: Path, tmp_path: Path) -> None:
        import pytest

        bad_yaml = tmp_path / "bad_heater.yaml"
        text = sample_case_path.read_text(encoding="utf-8").replace("width: 0.6", "width: 3.5")
        bad_yaml.write_text(text, encoding="utf-8")

        with pytest.raises(ValueError, match="Heater width"):
            build_case(bad_yaml, output_dir=tmp_path / "out")

    def test_out_of_bounds_heater_raises(self, sample_case_path: Path, tmp_path: Path) -> None:
        import pytest

        bad_yaml = tmp_path / "bad_heater_position.yaml"
        text = sample_case_path.read_text(encoding="utf-8").replace("z: 1.25", "z: 2.4", 1)
        bad_yaml.write_text(text, encoding="utf-8")

        with pytest.raises(ValueError, match="Heater z-range"):
            build_case(bad_yaml, output_dir=tmp_path / "out")

    def test_non_wall_heater_x_raises(self, sample_case_path: Path, tmp_path: Path) -> None:
        import pytest

        bad_yaml = tmp_path / "bad_heater_x.yaml"
        text = sample_case_path.read_text(encoding="utf-8").replace("x: 0.0", "x: 0.5", 1)
        bad_yaml.write_text(text, encoding="utf-8")

        with pytest.raises(ValueError, match="boundary_conditions.heater.position.x must be 0"):
            build_case(bad_yaml, output_dir=tmp_path / "out")

    def test_overwrites_existing(self, sample_case_path: Path, tmp_path: Path) -> None:
        out = tmp_path / "test_case"
        build_case(sample_case_path, output_dir=out)
        # Build again — should not fail
        build_case(sample_case_path, output_dir=out)
        assert (out / "system" / "blockMeshDict").is_file()
