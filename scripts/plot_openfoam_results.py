"""Plot OpenFOAM buoyantPimpleFoam results vs simple solver.

Usage: python scripts/plot_openfoam_results.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_openfoam_field(case_dir: Path, time: str, field: str) -> np.ndarray:
    """Parse an OpenFOAM internal field as numpy array."""
    fpath = case_dir / time / field
    text = fpath.read_text(encoding="utf-8")

    # Find the internalField block
    start = text.find("internalField")
    if start < 0:
        raise ValueError(f"No internalField in {fpath}")

    # Check if uniform
    rest = text[start:]
    first_line = rest.split("\n")[0]
    if "uniform" in first_line and "nonuniform" not in first_line:
        val_str = first_line.split("uniform")[1].split(";")[0].strip()
        return np.array([float(val_str)])

    # Non-uniform list: find count then ( ... ) block
    # Format: "nonuniform List<scalar> \n COUNT \n ( val1 val2 ... )"
    import re
    count_match = re.search(r'(\d+)\s*\(', rest)
    if not count_match:
        raise ValueError("Cannot parse field data")
    paren_start = rest.find("(", count_match.start())
    paren_end = rest.find(")", paren_start)
    data_str = rest[paren_start + 1:paren_end].strip()
    values = [float(x) for x in data_str.split()]
    return np.array(values)


def parse_cell_centres(case_dir: Path) -> np.ndarray | None:
    """Try to read cell centre y-coordinates from mesh."""
    # We'll reconstruct from blockMeshDict
    bmd = case_dir / "system" / "blockMeshDict"
    text = bmd.read_text(encoding="utf-8")

    # Extract dimensions from blocks section
    # Simple: assume uniform grid, extract from vertices
    # For now, compute y-positions from cell count
    return None


def extract_vertical_profile(case_dir: Path, time: str, field: str = "T",
                             n_bins: int = 40, room_h: float = 2.5) -> tuple[np.ndarray, np.ndarray]:
    """Extract vertical profile by running writeCellCentres or approximating.

    For multi-block meshes, read cell centres from the C file if available,
    otherwise approximate y from cell index ordering.
    """
    import subprocess, re

    values = parse_openfoam_field(case_dir, time, field)
    n_cells = len(values)

    # Try to read cell centres from constant/polyMesh or generate them
    cc_file = case_dir / "0" / "C"
    if not cc_file.exists():
        cc_file = case_dir / time / "C"

    if cc_file.exists():
        cc_text = cc_file.read_text(encoding="utf-8")
        # Parse vector field: (x y z) per cell
        vectors = re.findall(r'\(([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\)', cc_text)
        y_coords = np.array([float(v[1]) for v in vectors[:n_cells]])
    else:
        # Approximate: for structured mesh, assume uniform y distribution
        # Use n_bins approach with uniform spacing
        y_coords = np.linspace(room_h / (2 * n_cells), room_h, n_cells)

    # Bin average by y-coordinate
    bin_edges = np.linspace(0, room_h, n_bins + 1)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    profile = np.zeros(n_bins)
    counts = np.zeros(n_bins)

    for i, (y, v) in enumerate(zip(y_coords, values)):
        idx = min(int(y / room_h * n_bins), n_bins - 1)
        profile[idx] += v
        counts[idx] += 1

    mask = counts > 0
    profile[mask] /= counts[mask]
    # Fill empty bins with interpolation
    if not mask.all():
        profile = np.interp(bin_centres, bin_centres[mask], profile[mask])

    return bin_centres, profile


def main():
    case_dir = Path("results/openfoam_dry")

    # Find latest time directory
    time_dirs = []
    for d in case_dir.iterdir():
        if d.is_dir():
            try:
                float(d.name)
                if d.name != "0":
                    time_dirs.append(d.name)
            except ValueError:
                pass
    time_dirs.sort(key=float)

    if not time_dirs:
        print("No result time directories found!")
        return

    latest = time_dirs[-1]
    print(f"Latest time: {latest}")

    # Extract T profile
    y_of, T_of = extract_vertical_profile(case_dir, latest, "T")

    if len(y_of) == 0:
        print("Could not extract profile")
        return

    # Run simple solver for comparison
    import tempfile, yaml
    from harness.simple_solver import solve_two_zone

    tmp = Path(tempfile.mkdtemp())
    cfg = {
        "geometry": {"dimensions": {"x": 3.0, "y": 2.5, "z": 2.5}, "mesh_level": "M0"},
        "boundary_conditions": {
            "walls": {"temperature": 293.15, "model": "lumped", "thickness": 0.08, "conductivity": 0.12, "rho_cp": 500000},
            "heater": {"power_kw": 18.0, "position": {"x": 0, "y": 0.1, "z": 1.25}, "width": 0.6, "height": 0.5},
        },
        "case": {"name": "compare", "type": "steady"},
        "probes": [],
    }
    (tmp / "c.yaml").write_text(yaml.dump(cfg), encoding="utf-8")
    simple = solve_two_zone(tmp / "c.yaml", max_iter=10000)

    # Plot comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("OpenFOAM vs Simple Solver Comparison", fontsize=14, fontweight="bold")

    # Temperature profile
    ax = axes[0]
    ax.plot(T_of - 273.15, y_of, "b-o", markersize=3, linewidth=2, label=f"OpenFOAM (t={latest}s)")
    ax.plot(simple.temperatures - 273.15, simple.y_positions, "r--", linewidth=2, label="Simple (steady)")
    ax.axhline(simple.interface_height, color="green", linestyle=":", alpha=0.5, label=f"Interface z={simple.interface_height:.2f}m")
    ax.set_xlabel("Temperature [C]")
    ax.set_ylabel("Height [m]")
    ax.set_title("Vertical Temperature Profile")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Summary stats
    ax = axes[1]
    ax.axis("off")
    of_upper = T_of[y_of > 1.5].mean() - 273.15 if len(y_of[y_of > 1.5]) > 0 else 0
    of_lower = T_of[y_of < 1.0].mean() - 273.15 if len(y_of[y_of < 1.0]) > 0 else 0

    text = (
        f"OpenFOAM (t={latest}s, M0 mesh):\n"
        f"  T_upper (y>1.5m avg) = {of_upper:.1f} C\n"
        f"  T_lower (y<1.0m avg) = {of_lower:.1f} C\n"
        f"  Delta T             = {of_upper - of_lower:.1f} K\n"
        f"\n"
        f"Simple Solver (steady):\n"
        f"  T_upper             = {simple.upper_layer_temp - 273.15:.1f} C\n"
        f"  T_lower             = {simple.lower_layer_temp - 273.15:.1f} C\n"
        f"  Delta T             = {simple.upper_layer_temp - simple.lower_layer_temp:.1f} K\n"
        f"  z_interface         = {simple.interface_height:.2f} m\n"
        f"\n"
        f"Note: OpenFOAM ran only {latest}s (not steady).\n"
        f"Full run needs ~300s for thermal equilibrium."
    )
    ax.text(0.05, 0.95, text, transform=ax.transAxes, verticalalignment="top",
            fontsize=11, family="monospace",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax.set_title("Comparison Summary")

    plt.tight_layout()
    out_path = Path("results/openfoam_vs_simple.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
