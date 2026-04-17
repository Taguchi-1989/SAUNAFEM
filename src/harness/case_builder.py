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


def _build_block_mesh_context(
    geometry: dict,
    boundary_conditions: dict,
    ventilation: dict | None = None,
) -> dict:
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
    supply_vent_faces: list[tuple[int, int, int, int]] = []
    exhaust_vent_faces: list[tuple[int, int, int, int]] = []

    heater_y_range = (y0, y1)
    heater_z_range = (z0, z1)

    # Ventilation vent placement: supply on x=0 wall (lowest y, first z segment),
    # exhaust on opposite wall (highest y, first z segment).
    vent_enabled = ventilation is not None
    # Pick the first z-segment center for vent placement
    supply_vent_seg: tuple[int, int] | None = None
    exhaust_vent_seg: tuple[int, int] | None = None
    if vent_enabled:
        # Supply vent: lowest y-segment, first z-segment on x=0 wall
        supply_vent_seg = (0, 0)
        # Exhaust vent: highest y-segment, first z-segment on opposite wall
        exhaust_vent_seg = (len(y_points) - 2, 0)

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
            opp_face = (v1, v2, v6, v5)

            # Classify opposite wall face (exhaust vent or normal)
            if exhaust_vent_seg == (y_idx, z_idx):
                exhaust_vent_faces.append(opp_face)
            else:
                opposite_faces.append(opp_face)

            # Classify x=0 wall face (heater, supply vent, or surround)
            segment_y = (y_points[y_idx], y_points[y_idx + 1])
            segment_z = (z_points[z_idx], z_points[z_idx + 1])
            is_heater_face = segment_y == heater_y_range and segment_z == heater_z_range

            if is_heater_face:
                heater_faces.append(x0_face)
            elif supply_vent_seg == (y_idx, z_idx):
                supply_vent_faces.append(x0_face)
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

    result = {
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
    if vent_enabled:
        result["supply_vent_faces"] = supply_vent_faces
        result["exhaust_vent_faces"] = exhaust_vent_faces
    return result


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

    # Heater surface temperature: typical sauna stove surface 200-400°C
    # Use power-based estimate capped to realistic range
    # T_heater = T_wall + Q/(h_total * A), h_total ≈ 200 W/(m²K) (high rad+conv)
    h_eff_heater = 200.0
    t_heater_raw = t_walls + heat_flux / h_eff_heater
    t_heater = min(t_heater_raw, t_walls + 300)  # cap at wall+300K (~600K max)

    heater_model = heater.get("model", "surface_flux")
    heater_depth = heater.get("depth", 0.3)
    heater_power_w = power_kw * 1000.0
    heater_volume = heater_depth * height * width
    heater_power_density = heater_power_w / max(heater_volume, 1e-9)

    return {
        "heat_flux": round(heat_flux, 2),
        "heater_width": width,
        "heater_height": height,
        "T_walls": t_walls,
        "T_initial": t_walls,
        "T_heater": round(t_heater, 2),
        "heater_model": heater_model,
        "heater_depth": heater_depth,
        "heater_power_W": round(heater_power_w, 2),
        "heater_power_density": round(heater_power_density, 2),
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


def _generate_initial_fields(
    case_yaml: Path, output_dir: Path, geometry: dict, mesh: dict,
    boundary_conditions: dict,
) -> None:
    """Run simple solver and write nonuniform T/p_rgh initial fields.

    Uses the 2-zone plume model steady solution as initial condition
    for OpenFOAM, drastically reducing convergence time.
    """
    import numpy as np

    from harness.simple_solver import solve_two_zone

    try:
        result = solve_two_zone(case_yaml, max_iter=10000)
    except Exception:
        return  # If simple solver fails, keep uniform fields

    if not result.converged:
        # Still use the result — it's a better guess than uniform 20°C
        pass

    dims = geometry["dimensions"]
    room_h = dims["y"]

    # Get mesh y-structure from blockMeshDict context
    heater = boundary_conditions.get("heater", {})
    h_width = heater.get("width", 0.5)
    h_height = heater.get("height", 0.5)
    h_pos = heater.get("position", {})
    y0 = h_pos.get("y", 0.0)
    y1 = y0 + h_height
    z_center = h_pos.get("z", dims["z"] / 2)
    z0 = z_center - h_width / 2
    z1 = z0 + h_width

    y_points = _unique_points([0.0, y0, y1, dims["y"]])
    z_points = _unique_points([0.0, z0, z1, dims["z"]])

    density = _MESH_DENSITY.get(geometry.get("mesh_level", "M0"), 8)
    nx = max(1, round(dims["x"] * density))
    ny = max(1, round(dims["y"] * density))
    nz = max(1, round(dims["z"] * density))

    y_cells = _allocate_segment_cells(y_points, ny)
    z_cells = _allocate_segment_cells(z_points, nz)

    # Build cell centre y-coordinates for each block (ordered as blockMesh emits them)
    cell_y_list: list[float] = []
    for yi in range(len(y_points) - 1):
        for zi in range(len(z_points) - 1):
            n_y = y_cells[yi]
            n_z = z_cells[zi]
            y_lo, y_hi = y_points[yi], y_points[yi + 1]
            for j in range(n_y):
                y_c = y_lo + (j + 0.5) * (y_hi - y_lo) / n_y
                for _ in range(nx * n_z):
                    cell_y_list.append(y_c)

    cell_y = np.array(cell_y_list)
    n_cells = len(cell_y)

    # Map simple solver profile to each cell
    t_upper = result.upper_layer_temp
    t_lower = result.lower_layer_temp
    z_int = result.interface_height
    delta = 0.15 * room_h  # transition layer thickness

    # Sigmoid interpolation (same as simple_solver._build_profile_and_probes)
    T_field = t_lower + (t_upper - t_lower) / (
        1.0 + np.exp(-3.0 * (cell_y - z_int) / (delta / 2.0))
    )

    # Corresponding p_rgh: hydrostatic adjustment
    # p_rgh = p - rho * g * (y - hRef), with hRef=0
    # rho = rho_0 * T_ref / T, rho_0=1.1, T_ref=300K (approximate)
    R_air = 8314.0 / 28.96  # J/(kg·K)
    p_atm = 101325.0
    rho_field = p_atm / (R_air * T_field)
    g_mag = 9.81
    p_rgh_field = p_atm - rho_field * g_mag * cell_y

    # Write T file with nonuniform internalField
    t_file = output_dir / "0" / "T"
    t_text = t_file.read_text(encoding="utf-8")
    # Replace "internalField   uniform XXX;" with nonuniform list
    t_values = "\n".join(f"{v:.6f}" for v in T_field)
    t_text = t_text.replace(
        f"internalField   uniform {result.lower_layer_temp};",  # might not match
        f"internalField   nonuniform List<scalar>\n{n_cells}\n(\n{t_values}\n)\n;",
    )
    # More robust: replace any "internalField   uniform ...;"
    import re
    t_text = re.sub(
        r'internalField\s+uniform\s+[\d.]+\s*;',
        f'internalField   nonuniform List<scalar>\n{n_cells}\n(\n{t_values}\n)\n;',
        t_text,
        count=1,
    )
    t_file.write_text(t_text, encoding="utf-8")

    # Write p_rgh with nonuniform internalField
    prgh_file = output_dir / "0" / "p_rgh"
    prgh_text = prgh_file.read_text(encoding="utf-8")
    prgh_values = "\n".join(f"{v:.6f}" for v in p_rgh_field)
    prgh_text = re.sub(
        r'internalField\s+uniform\s+[\d.]+\s*;',
        f'internalField   nonuniform List<scalar>\n{n_cells}\n(\n{prgh_values}\n)\n;',
        prgh_text,
        count=1,
    )
    prgh_file.write_text(prgh_text, encoding="utf-8")


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
    vent_cfg = data.get("ventilation")
    vent_enabled = vent_cfg is not None and vent_cfg.get("model", "none") != "none"
    mesh = _build_block_mesh_context(
        data["geometry"],
        data["boundary_conditions"],
        ventilation=vent_cfg if vent_enabled else None,
    )
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

    # Determine mixture type: explicit solver.mixture overrides auto-detection
    loyly = data.get("loyly")
    mixture_type = solver.get("mixture", "multiComponent" if loyly else "pure")

    # Aufguss jet momentum source parameters
    aufguss = data.get("aufguss")
    aufguss_enabled = aufguss is not None

    # Buoyancy production term for k-equation (enabled by default)
    turbulence = data.get("turbulence", {})
    buoyancy_production = turbulence.get("buoyancy_production", True)

    # Species transport (needed for multiComponentMixture)
    species_transport = mixture_type == "multiComponent"

    # Radiation model: "none", "viewFactor", or "fvDOM"
    radiation_cfg = data.get("radiation", {})
    radiation_model = radiation_cfg.get("model", "none") if radiation_cfg else "none"

    # Run simple solver to get wall inner temperature for BC initialization
    t_wall_inner = heater["T_walls"]  # fallback to external wall temp
    try:
        from harness.simple_solver import solve_two_zone
        _simple = solve_two_zone(case_yaml, max_iter=10000)
        t_wall_inner = round(_simple.wall_inner_temp, 2)
    except Exception:
        pass

    # Wall heat transfer: lambda_w / delta_w gives effective U-value to outside
    walls = data["boundary_conditions"].get("walls", {})
    wall_thickness = walls.get("thickness", 0.08)
    wall_conductivity = walls.get("conductivity", 0.12)
    wall_htc_override = walls.get("wall_htc_override")
    if wall_htc_override is not None:
        wall_htc = wall_htc_override
    else:
        wall_htc = wall_conductivity / max(wall_thickness, 0.01)  # W/(m²K)

    context = {
        **mesh,
        **heater,
        "T_wall_inner": t_wall_inner,
        "wall_htc": round(wall_htc, 4),
        "solver_name": solver_name,
        "end_time": end_time,
        "write_interval": write_interval,
        "delta_t": delta_t,
        "delta_t_initial": delta_t,
        "max_delta_t": delta_t * 10,
        "averaging_start": averaging_start,
        "probes": probes,
        "mixture_type": mixture_type,
        "Y_H2O_initial": 0.01,
        "Y_H2O_heater": 0.01,
        "aufguss_enabled": aufguss_enabled,
        "aufguss_jet_velocity": aufguss.get("jet_velocity", 2.0) if aufguss_enabled else 0.0,
        "aufguss_duration": aufguss.get("duration", 1.0) if aufguss_enabled else 1.0,
        "buoyancy_production": buoyancy_production,
        "species_transport": species_transport,
        "radiation_model": radiation_model,
        "turbulence_model": data.get("turbulence", {}).get("model", "kOmegaSST"),
        "convection_scheme": data.get("turbulence", {}).get("convection_scheme", "linearUpwind"),
        "ventilation": vent_enabled,
        "T_ambient": vent_cfg.get("T_ambient", 293.15) if vent_enabled else 293.15,
    }

    # Skip vapor field template for pure mixture cases
    skip: list[str] = []
    if mixture_type != "multiComponent":
        skip.append("0/H2O.j2")

    # Skip IDefault template when not using fvDOM radiation
    if radiation_model != "fvDOM":
        skip.append("0/IDefault.j2")

    render_templates(_TEMPLATE_DIR, output_dir, context, skip_templates=skip)

    # Initialize T and p_rgh from simple solver steady solution
    uniform_init = solver.get("uniform_init", False)
    if not uniform_init:
        _generate_initial_fields(
            case_yaml, output_dir, data["geometry"], mesh, data["boundary_conditions"],
        )

    return output_dir
