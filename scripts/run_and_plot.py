"""Run 3 simulation scenarios and generate comparison plots.

Usage: python scripts/run_and_plot.py
Output: results/simulation_results.png
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from harness.simple_solver import solve_transient, solve_two_zone

# --- Common room config ---
BASE_CFG = {
    "geometry": {
        "dimensions": {"x": 3.0, "y": 2.5, "z": 2.5},
        "mesh_level": "M0",
    },
    "boundary_conditions": {
        "walls": {
            "temperature": 293.15,
            "model": "lumped",
            "thickness": 0.08,       # 80mm insulated wood panel
            "conductivity": 0.12,
            "rho_cp": 500000,
        },
        "heater": {
            "power_kw": 18.0,        # ~1 kW/m3 for 18.75 m3 room
            "position": {"x": 0.0, "y": 0.1, "z": 1.25},
            "width": 0.6,
            "height": 0.5,
        },
    },
    "probes": [
        {"name": "upper_bench", "position": {"x": 1.5, "y": 2.0, "z": 1.25}, "fields": ["T"]},
        {"name": "lower_bench", "position": {"x": 1.5, "y": 0.8, "z": 1.25}, "fields": ["T"]},
    ],
}


def write_yaml(cfg: dict, directory: Path) -> Path:
    p = directory / "case.yaml"
    p.write_text(yaml.dump(cfg, allow_unicode=True), encoding="utf-8")
    return p


def run_scenarios():
    tmp = Path(tempfile.mkdtemp())

    # --- Scenario 1: Dry steady-state ---
    print("Running: Dry steady-state...")
    dry_dir = tmp / "dry"
    dry_dir.mkdir()
    dry_cfg = {**BASE_CFG, "case": {"name": "dry_steady", "type": "steady"}}
    dry_path = write_yaml(dry_cfg, dry_dir)
    dry_result = solve_two_zone(dry_path, max_iter=10000)

    # --- Scenario 2: Löyly transient ---
    print("Running: Loyly transient (500 mL)...")
    loyly_dir = tmp / "loyly"
    loyly_dir.mkdir()
    loyly_cfg = {
        **BASE_CFG,
        "case": {"name": "loyly_transient", "type": "transient"},
        "loyly": {"water_ml": 500, "time": 30.0, "tau_evap": 5.0},
    }
    loyly_path = write_yaml(loyly_cfg, loyly_dir)
    loyly_result = solve_transient(loyly_path, end_time=300.0, physical_dt=0.5, record_interval=1.0)

    # --- Scenario 3: Aufguss transient ---
    print("Running: Aufguss transient (beta=0.5)...")
    aufguss_dir = tmp / "aufguss"
    aufguss_dir.mkdir()
    aufguss_cfg = {
        **BASE_CFG,
        "case": {"name": "aufguss_transient", "type": "transient"},
        "aufguss": {"beta_aug": 0.5, "start_time": 60.0, "duration": 30.0},
    }
    aufguss_path = write_yaml(aufguss_cfg, aufguss_dir)
    aufguss_result = solve_transient(aufguss_path, end_time=200.0, physical_dt=0.5, record_interval=1.0)

    return dry_result, loyly_result, aufguss_result


def plot_results(dry, loyly, aufguss):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("SaunaFlow Simulation Results", fontsize=16, fontweight="bold")

    # Color scheme
    c_upper = "#d62728"
    c_lower = "#1f77b4"
    c_interface = "#2ca02c"
    c_humidity = "#9467bd"

    # --- (0,0) Dry steady-state temperature profile ---
    ax = axes[0, 0]
    ax.plot(dry.temperatures - 273.15, dry.y_positions, color=c_upper, linewidth=2)
    ax.axhline(dry.interface_height, color=c_interface, linestyle="--", alpha=0.7, label=f"Interface z={dry.interface_height:.2f}m")
    ax.set_xlabel("Temperature [°C]")
    ax.set_ylabel("Height [m]")
    ax.set_title("Dry Steady-State Profile")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (0,1) Löyly time series: T_upper, T_lower ---
    ax = axes[0, 1]
    ax.plot(loyly.time, loyly.t_upper_series - 273.15, color=c_upper, linewidth=1.5, label="T_upper")
    ax.plot(loyly.time, loyly.t_lower_series - 273.15, color=c_lower, linewidth=1.5, label="T_lower")
    ax.axvline(30, color="gray", linestyle=":", alpha=0.5, label="Löyly @ 30s")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.set_title("Löyly Transient (500 mL)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (0,2) Löyly humidity ---
    ax = axes[0, 2]
    ax.plot(loyly.time, loyly.humidity_series * 1000, color=c_humidity, linewidth=1.5)
    ax.axvline(30, color="gray", linestyle=":", alpha=0.5, label="Löyly @ 30s")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Humidity ratio [g/kg]")
    ax.set_title("Löyly Humidity Response")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (1,0) Aufguss time series ---
    ax = axes[1, 0]
    ax.plot(aufguss.time, aufguss.t_upper_series - 273.15, color=c_upper, linewidth=1.5, label="T_upper")
    ax.plot(aufguss.time, aufguss.t_lower_series - 273.15, color=c_lower, linewidth=1.5, label="T_lower")
    ax.axvspan(60, 90, color="orange", alpha=0.15, label="Aufguss window")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Temperature [°C]")
    ax.set_title("Aufguss Transient (β=0.5)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (1,1) Interface height comparison ---
    ax = axes[1, 1]
    ax.plot(loyly.time, loyly.z_int_series, color=c_interface, linewidth=1.5, label="Löyly")
    ax.plot(aufguss.time, aufguss.z_int_series, color="orange", linewidth=1.5, label="Aufguss")
    ax.axhline(dry.interface_height, color="gray", linestyle="--", alpha=0.5, label=f"Dry steady ({dry.interface_height:.2f}m)")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Interface height [m]")
    ax.set_title("Interface Height Evolution")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- (1,2) Perceived temperature (Löyly) ---
    ax = axes[1, 2]
    ax.plot(loyly.time, loyly.perceived_upper_series, color=c_upper, linewidth=1.5, label="Perceived T (upper)")
    ax.axvline(30, color="gray", linestyle=":", alpha=0.5, label="Löyly @ 30s")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Perceived Temperature [°C]")
    ax.set_title("Skin Heat Balance (K-06)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Summary text ---
    fig.text(0.02, 0.02,
        f"Dry: T_upper={dry.upper_layer_temp-273.15:.1f}°C  T_lower={dry.lower_layer_temp-273.15:.1f}°C  "
        f"ΔT={dry.upper_layer_temp-dry.lower_layer_temp:.1f}K  z_int={dry.interface_height:.2f}m  "
        f"Perceived={dry.perceived_temp_upper:.1f}°C",
        fontsize=9, family="monospace", color="gray")

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])

    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "simulation_results.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nSaved: {out_path}")
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    dry, loyly, aufguss = run_scenarios()
    out = plot_results(dry, loyly, aufguss)
    print("Done.")
