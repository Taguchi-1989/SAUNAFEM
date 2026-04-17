"""Tests for YAML schema validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from harness.schema import load_and_validate, load_yaml, validate_case


class TestLoadYaml:
    def test_loads_valid_yaml(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        assert isinstance(data, dict)
        assert "case" in data

    def test_case_name(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        assert data["case"]["name"] == "dry_sauna_steady"


class TestValidateCase:
    def test_valid_case_no_errors(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        errors = validate_case(data)
        assert errors == []

    def test_missing_required_field(self) -> None:
        data = {"case": {"name": "test", "description": "test"}}
        errors = validate_case(data)
        assert len(errors) > 0
        assert any("required" in e.lower() for e in errors)

    def test_invalid_solver_name(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        data["solver"]["name"] = "invalidSolver"
        errors = validate_case(data)
        assert len(errors) > 0

    def test_negative_dimension_rejected(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        data["geometry"]["dimensions"]["x"] = -1.0
        errors = validate_case(data)
        assert len(errors) > 0

    def test_empty_dict_rejected(self) -> None:
        errors = validate_case({})
        assert len(errors) > 0


class TestHeaterModelSchema:
    def test_surface_flux_valid(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        data["boundary_conditions"]["heater"]["model"] = "surface_flux"
        errors = validate_case(data)
        assert errors == []

    def test_volume_source_valid(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        data["boundary_conditions"]["heater"]["model"] = "volume_source"
        data["boundary_conditions"]["heater"]["depth"] = 0.3
        errors = validate_case(data)
        assert errors == []

    def test_invalid_model_rejected(self, sample_case_path: Path) -> None:
        data = load_yaml(sample_case_path)
        data["boundary_conditions"]["heater"]["model"] = "invalid"
        errors = validate_case(data)
        assert len(errors) > 0

    def test_no_model_defaults_ok(self, sample_case_path: Path) -> None:
        """Existing YAML without heater.model should still pass."""
        data = load_yaml(sample_case_path)
        assert "model" not in data["boundary_conditions"]["heater"]
        errors = validate_case(data)
        assert errors == []


class TestLoadAndValidate:
    def test_valid_file(self, sample_case_path: Path) -> None:
        errors = load_and_validate(sample_case_path)
        assert errors == []

    def test_nonexistent_file_raises(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_and_validate("/nonexistent/path.yaml")
