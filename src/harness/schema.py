"""YAML loading and JSON Schema validation."""

from __future__ import annotations

import json
from pathlib import Path

import jsonschema
import yaml

_SCHEMA_DIR = Path(__file__).resolve().parent.parent.parent / "configs" / "schemas"
_CASE_SCHEMA_PATH = _SCHEMA_DIR / "case_schema.json"


def _load_schema() -> dict:
    """Load the case JSON Schema from disk."""
    with open(_CASE_SCHEMA_PATH, encoding="utf-8") as f:
        return json.load(f)


def load_yaml(path: str | Path) -> dict:
    """Load a YAML file and return its contents as a dict."""
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def validate_case(data: dict) -> list[str]:
    """Validate case data against the JSON Schema.

    Returns a list of error messages (empty if valid).
    """
    schema = _load_schema()
    validator = jsonschema.Draft202012Validator(schema)
    return [err.message for err in validator.iter_errors(data)]


def load_and_validate(path: str | Path) -> list[str]:
    """Load a YAML case file and validate it.

    Returns a list of error messages (empty if valid).
    """
    data = load_yaml(path)
    return validate_case(data)
