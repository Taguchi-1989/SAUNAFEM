"""Mesh generation execution."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from utils.wsl import wsl_exec


@dataclass
class MeshResult:
    """Result of mesh generation."""

    cell_count: int
    success: bool
    log: str
    quality: dict[str, float] = field(default_factory=dict)


def parse_check_mesh(output: str) -> tuple[int, dict[str, float]]:
    """Parse checkMesh output for cell count and quality metrics.

    Returns (cell_count, quality_dict).
    """
    cell_count = 0
    quality: dict[str, float] = {}

    # OpenFOAM variants report either "cells:" or "Number of cells ="
    m = re.search(r"cells:\s+(\d+)", output)
    if not m:
        m = re.search(r"Number of cells\s*=\s*(\d+)", output)
    if m:
        cell_count = int(m.group(1))

    # Extract max aspect ratio
    m = re.search(r"Max aspect ratio\s*=\s*([\d.eE+-]+)", output)
    if m:
        quality["max_aspect_ratio"] = float(m.group(1))

    # Extract max non-orthogonality
    m = re.search(r"Mesh non-orthogonality.*Max:\s*([\d.eE+-]+)", output)
    if m:
        quality["max_non_orthogonality"] = float(m.group(1))

    # Extract max skewness
    m = re.search(r"Max skewness\s*=\s*([\d.eE+-]+)", output)
    if m:
        quality["max_skewness"] = float(m.group(1))

    return cell_count, quality


def run_mesh(case_dir: Path, check: bool = True) -> MeshResult:
    """Run blockMesh (and optionally checkMesh) on the case directory.

    Args:
        case_dir: Path to the OpenFOAM case directory.
        check: Whether to run checkMesh after blockMesh.

    Returns:
        MeshResult with cell count and quality information.
    """
    # Run blockMesh
    result = wsl_exec("blockMesh", cwd=case_dir)
    log = result.stdout

    cell_count = 0
    quality: dict[str, float] = {}

    if check:
        check_result = wsl_exec("checkMesh", cwd=case_dir)
        log += "\n" + check_result.stdout
        cell_count, quality = parse_check_mesh(check_result.stdout)

    return MeshResult(
        cell_count=cell_count,
        success=result.returncode == 0 and (not check or check_result.returncode == 0),
        log=log,
        quality=quality,
    )
