"""Case definition -> OpenFOAM directory structure builder."""

from __future__ import annotations

import shutil
from pathlib import Path

import jinja2

from harness.schema import load_and_validate, load_yaml

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_TEMPLATE_DIR = _PROJECT_ROOT / "foam_templates" / "base_case"
_RESULTS_DIR = _PROJECT_ROOT / "results"

# Cells per meter for each mesh level
_MESH_DENSITY: dict[str, int] = {
    "M0": 8,   # ~8/m -> 24x20x20 = ~9600 cells (coarse 2D-ish)
    "M1": 16,  # ~16/m -> ~48x40x40 = ~76800
    "M2": 28,  # ~28/m -> ~84x70x70 = ~411600
    "M3": 40,  # ~40/m -> ~120x100x100 = ~1200000
}


def _validate_positive_patch_size(size: float, limit: float, axis_name: str) -> float:
    """Validate that a heater patch size fits inside the domain."""
    if size <= 0:
        raise ValueError(f"Heater {axis_name} must be positive, got {size}.")
    if size > limit:
        raise ValueError(
            f"Heater {axis_name}={size} exceeds domain {axis_name}={limit}."
        )
    return size


def _unique_points(points: list[float]) -> list[float]:
    """Drop duplicate split points while preserving order."""
    result: list[float] = []
    for point in points:
        if not result or abs(point - result[-1]) > 1e-9:
            result.append(point)
    return result


def _allocate_segment_cells(points: list[float], total_cells: int) -> list[int]:
    """Allocate at least one cell to each segment while preserving total count."""
    spans = [points[i + 1] - points[i] for i in range(len(points) - 1)]
    total_span = sum(spans)
    if total_span <= 0:
        raise ValueError("Mesh split points must define a positive span.")

    cells = [
        max(1, round(total_cells * span / total_span))
        for span in spans
    ]

    diff = total_cells - sum(cells)
    while diff != 0:
        if diff > 0:
            idx = max(range(len(cells)), key=lambda i: spans[i])
            cells[idx] += 1
            diff -= 1
            continue

        candidates = [i for i, cell_count in enumerate(cells) if cell_count > 1]
        if not candidates:
            break
        idx = max(candidates, key=lambda i: spans[i])
        cells[idx] -= 1
        diff += 1

    return cells


def _build_block_mesh_context(geometry: dict, boundary_conditions: dict) -> dict:
    """Build a segmented mesh so the heater occupies its own wall patch."""
    mesh = compute_mesh_params(geometry)
    dims = geometry["dimensions"]
    heater = boundary_conditions.get("heater", {})

    width = _validate_positive_patch_size(
        heater.get("width", 0.5), dims["z"], "width"
    )
    height = _validate_positive_patch_size(
        heater.get("height", 0.5), dims["y"], "height"
    )

    position = heater.get("position", {})
    x_pos = position.get("x", 0.0)
    if abs(x_pos) > 1e-9:
        raise ValueError(
            "Phase 1 heater must be mounted on the x=0 wall, so heater.position.x must be 0."
        )
    y0 = position.get("y", 0.0)
    if y0 < 0 or y0 + height > dims["y"]:
        raise ValueError(
            f"Heater y-range [{y0}, {y0 + height}] must stay within [0, {dims['y']}]."
        )
    y1 = y0 + height
    z_center = position.get("z", dims["z"] / 2)
    z0 = z_center - width / 2
    z1 = z_center + width / 2
    if z0 < 0 or z1 > dims["z"]:
        raise ValueError(
            f"Heater z-range [{z0}, {z1}] must stay within [0, {dims['z']}]."
        )
    z1 = z0 + width

    y_points = _unique_points([0.0, y0, y1, dims["y"]])
    z_points = _unique_points([0.0, z0, z1, dims["z"]])
    y_cells = _allocate_segment_cells(y_points, mesh["ny"])
    z_cells = _allocate_segment_cells(z_points, mesh["nz"])

    nz_points = len(z_points)
    plane_size = len(y_points) * nz_points

    def vertex_index(x_side: int, y_idx: int, z_idx: int) -> int:
        return x_side * plane_size + y_idx * nz_points + z_idx

    vertices: list[tuple[float, float, float]] = []
    for x in (0.0, dims["x"]):
        for y in y_points:
            for z in z_points:
                vertices.append((x, y, z))

    blocks: list[dict[str, object]] = []
    floor_faces: list[tuple[int, int, int, int]] = []
    ceiling_faces: list[tuple[int, int, int, int]] = []
    front_faces: list[tuple[int, int, int, int]] = []
    back_faces: list[tuple[int, int, int, int]] = []
    heater_faces: list[tuple[int, int, int, int]] = []
    heater_surround_faces: list[tuple[int, int, int, int]] = []
    opposite_faces: list[tuple[int, int, int, int]] = []

    heater_y_range = (y0, y1)
    heater_z_range = (z0, z1)

    for y_idx in range(len(y_points) - 1):
        for z_idx in range(len(z_points) - 1):
            v0 = vertex_index(0, y_idx, z_idx)
            v1 = vertex_index(1, y_idx, z_idx)
            v2 = vertex_index(1, y_idx + 1, z_idx)
            v3 = vertex_index(0, y_idx + 1, z_idx)
            v4 = vertex_index(0, y_idx, z_idx + 1)
            v5 = vertex_index(1, y_idx, z_idx + 1)
            v6 = vertex_index(1, y_idx + 1, z_idx + 1)
            v7 = vertex_index(0, y_idx + 1, z_idx + 1)

            blocks.append({
                "vertices": (v0, v1, v2, v3, v4, v5, v6, v7),
                "cells": (mesh["nx"], y_cells[y_idx], z_cells[z_idx]),
            })

            x0_face = (v0, v4, v7, v3)
            opposite_faces.append((v1, v2, v6, v5))

            segment_y = (y_points[y_idx], y_points[y_idx + 1])
            segment_z = (z_points[z_idx], z_points[z_idx + 1])
            is_heater_face = segment_y == heater_y_range and segment_z == heater_z_range
            if is_heater_face:
                heater_faces.append(x0_face)
            else:
                heater_surround_faces.append(x0_face)

            if y_idx == 0:
                floor_faces.append((v0, v1, v5, v4))
            if y_idx == len(y_points) - 2:
                ceiling_faces.append((v3, v7, v6, v2))
            if z_idx == 0:
                front_faces.append((v0, v3, v2, v1))
            if z_idx == len(z_points) - 2:
                back_faces.append((v4, v5, v6, v7))

    return {
        **mesh,
        "heater_area": round(height * width, 6),
        "heater_y0": round(y0, 6),
        "heater_y1": round(y1, 6),
        "heater_z0": round(z0, 6),
        "heater_z1": round(z1, 6),
        "vertices": vertices,
        "blocks": blocks,
        "floor_faces": floor_faces,
        "ceiling_faces": ceiling_faces,
        "heater_faces": heater_faces,
        "heater_surround_faces": heater_surround_faces,
        "opposite_faces": opposite_faces,
        "front_faces": front_faces,
        "back_faces": back_faces,
    }


def compute_mesh_params(geometry: dict) -> dict:
    """Compute blockMesh parameters from geometry config.

    Returns dict with dim_x/y/z and nx/ny/nz cell counts.
    """
    dims = geometry["dimensions"]
    level = geometry.get("mesh_level", "M0")
    density = _MESH_DENSITY.get(level, _MESH_DENSITY["M0"])

    return {
        "dim_x": dims["x"],
        "dim_y": dims["y"],
        "dim_z": dims["z"],
        "nx": max(4, round(dims["x"] * density)),
        "ny": max(4, round(dims["y"] * density)),
        "nz": max(4, round(dims["z"] * density)),
    }


def compute_heater_params(boundary_conditions: dict, geometry: dict) -> dict:
    """Compute heater heat flux and related parameters.

    Returns dict with heat_flux (W/m2) and heater geometry.
    """
    heater = boundary_conditions.get("heater", {})
    power_kw = heater.get("power_kw", 9.0)
    width = _validate_positive_patch_size(
        heater.get("width", 0.5), geometry["dimensions"]["z"], "width"
    )
    height = _validate_positive_patch_size(
        heater.get("height", 0.5), geometry["dimensions"]["y"], "height"
    )
    heater_area = width * height
    heat_flux = (power_kw * 1000.0) / heater_area

    walls = boundary_conditions.get("walls", {})
    t_walls = walls.get("temperature", 293.15)

    return {
        "heat_flux": round(heat_flux, 2),
        "heater_width": width,
        "heater_height": height,
        "T_walls": t_walls,
        "T_initial": t_walls,
    }


def _build_probe_context(probes: list[dict]) -> list[dict]:
    """Convert YAML probe definitions to template context."""
    result = []
    for p in probes:
        pos = p["position"]
        result.append({
            "name": p["name"],
            "x": pos["x"],
            "y": pos["y"],
            "z": pos["z"],
        })
    return result


def render_templates(
    template_dir: Path,
    output_dir: Path,
    context: dict,
    skip_templates: list[str] | None = None,
) -> None:
    """Render all .j2 templates into the output directory.

    Walks template_dir, renders each .j2 file with Jinja2,
    and writes the result to the corresponding path in output_dir
    (stripping the .j2 extension). Non-.j2 files are copied as-is.

    Args:
        template_dir: Path to the template directory.
        output_dir: Path to the output directory.
        context: Template variables.
        skip_templates: List of relative template paths to skip (e.g. ["0/Y.H2O.j2"]).
    """
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_dir)),
        keep_trailing_newline=True,
        undefined=jinja2.StrictUndefined,
    )
    _skip = {s.replace("\\", "/") for s in (skip_templates or [])}

    for template_path in template_dir.rglob("*"):
        if template_path.is_dir():
            continue

        rel = template_path.relative_to(template_dir)
        rel_posix = str(rel).replace("\\", "/")

        if rel_posix in _skip:
            continue

        if template_path.suffix == ".j2":
            # Render Jinja2 template
            template = env.get_template(rel_posix)
            rendered = template.render(**context)
            out_path = output_dir / rel.with_suffix("")  # strip .j2
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(rendered, encoding="utf-8")
        else:
            # Copy non-template files as-is
            out_path = output_dir / rel
            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(template_path, out_path)


def build_case(case_yaml: Path, output_dir: Path | None = None) -> Path:
    """Build an OpenFOAM case directory from a YAML definition.

    Args:
        case_yaml: Path to the YAML case definition file.
        output_dir: Output directory. Defaults to results/{case_name}/.

    Returns:
        Path to the created case directory.

    Raises:
        ValueError: If the YAML fails schema validation.
    """
    errors = load_and_validate(case_yaml)
    if errors:
        raise ValueError(f"Schema validation failed: {'; '.join(errors)}")

    data = load_yaml(case_yaml)
    case_name = data["case"]["name"]

    if output_dir is None:
        output_dir = _RESULTS_DIR / case_name

    # Clean and recreate output directory
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    # Build template context
    mesh = _build_block_mesh_context(data["geometry"], data["boundary_conditions"])
    heater = compute_heater_params(data["boundary_conditions"], data["geometry"])
    solver = data["solver"]
    probes = _build_probe_context(data.get("probes", []))

    solver_name = solver["name"]
    end_time = solver.get("end_time", 1000)
    write_interval = solver.get("write_interval", 100)

    # Transient solver parameters
    if solver_name == "buoyantPimpleFoam":
        delta_t = solver.get("delta_t", 0.1)
        averaging_start = solver.get("averaging_start", end_time * 0.5)
    else:
        delta_t = 1
        averaging_start = 0

    # Determine mixture type from loyly config
    loyly = data.get("loyly")
    mixture_type = "multiComponent" if loyly else "pure"

    context = {
        **mesh,
        **heater,
        "solver_name": solver_name,
        "end_time": end_time,
        "write_interval": write_interval,
        "delta_t": delta_t,
        "averaging_start": averaging_start,
        "probes": probes,
        "mixture_type": mixture_type,
        "Y_H2O_initial": 0.01,
        "Y_H2O_heater": 0.01,
    }

    # Skip vapor field template for pure mixture cases
    skip: list[str] = []
    if mixture_type != "multiComponent":
        skip.append("0/H2O.j2")

    render_templates(_TEMPLATE_DIR, output_dir, context, skip_templates=skip)
    return output_dir
