"""SaunaFlow Visualization UI.

Run with: streamlit run ui/app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from harness.kpi import evaluate_phase1_kpis
from harness.simple_solver import solve_two_zone

st.set_page_config(page_title="SaunaFlow", page_icon="♨", layout="wide")

st.title("SaunaFlow - Sauna Thermal Simulation")
st.markdown("*Two-zone plume model (Morton-Taylor-Turner entrainment)*")

# --- Sidebar: Parameters ---
st.sidebar.header("Case Parameters")

case_yaml = Path(__file__).resolve().parent.parent / "configs" / "cases" / "dry_sauna_steady.yaml"

st.sidebar.subheader("Room Dimensions")
room_x = st.sidebar.slider("Width (x) [m]", 1.0, 6.0, 3.0, 0.1)
room_y = st.sidebar.slider("Height (y) [m]", 1.5, 4.0, 2.5, 0.1)
room_z = st.sidebar.slider("Depth (z) [m]", 1.0, 6.0, 2.5, 0.1)

st.sidebar.subheader("Heater")
heater_power = st.sidebar.slider("Power [kW]", 1.0, 20.0, 9.0, 0.5)
heater_y = st.sidebar.slider("Heater Y position [m]", 0.0, room_y - 0.5, 0.1, 0.1)
heater_height = st.sidebar.slider("Heater height [m]", 0.2, 1.0, 0.5, 0.1)
heater_width = st.sidebar.slider("Heater width [m]", 0.2, 1.5, 0.6, 0.1)

st.sidebar.subheader("Environment")
t_wall_c = st.sidebar.slider("Wall temperature [C]", 10.0, 30.0, 20.0, 1.0)
t_wall_k = t_wall_c + 273.15

st.sidebar.subheader("Probes")
upper_y = st.sidebar.slider("Upper bench Y [m]", 0.5, room_y, 2.0, 0.1)
lower_y = st.sidebar.slider("Lower bench Y [m]", 0.2, room_y, 0.8, 0.1)

# --- Build dynamic YAML-like data and solve ---
# We'll create a temporary YAML with the sidebar values
import tempfile

import yaml

case_data = {
    "case": {
        "name": "ui_interactive",
        "description": "Interactive UI case",
        "type": "steady",
    },
    "geometry": {
        "dimensions": {"x": room_x, "y": room_y, "z": room_z},
        "mesh_level": "M0",
    },
    "boundary_conditions": {
        "walls": {"temperature": t_wall_k, "type": "mixed"},
        "heater": {
            "power_kw": heater_power,
            "position": {"x": 0.0, "y": heater_y, "z": room_z / 2},
            "width": heater_width,
            "height": heater_height,
        },
    },
    "solver": {
        "name": "buoyantSimpleFoam",
        "end_time": 1000,
        "write_interval": 100,
    },
    "probes": [
        {"name": "upper_bench", "position": {"x": room_x / 2, "y": upper_y, "z": room_z / 2}, "fields": ["T"]},
        {"name": "lower_bench", "position": {"x": room_x / 2, "y": lower_y, "z": room_z / 2}, "fields": ["T"]},
        {"name": "floor_level", "position": {"x": room_x / 2, "y": 0.1, "z": room_z / 2}, "fields": ["T"]},
    ],
}

# Write temp YAML and solve
with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False, encoding="utf-8") as f:
    yaml.dump(case_data, f, default_flow_style=False)
    tmp_yaml = Path(f.name)

try:
    result = solve_two_zone(tmp_yaml, n_profile=80, max_iter=10000)
finally:
    tmp_yaml.unlink(missing_ok=True)

# --- Results ---
col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("Vertical Temperature Profile")

    fig, ax = plt.subplots(figsize=(8, 6))

    # Temperature profile
    t_celsius = result.temperatures - 273.15
    ax.plot(t_celsius, result.y_positions, "b-", linewidth=2.5, label="Temperature")

    # Probe markers
    colors = {"upper_bench": "#e74c3c", "lower_bench": "#f39c12", "floor_level": "#3498db"}
    labels = {"upper_bench": "Upper Bench", "lower_bench": "Lower Bench", "floor_level": "Floor Level"}
    for name, temp_k in result.probe_values.items():
        probe_y = next(p["position"]["y"] for p in case_data["probes"] if p["name"] == name)
        ax.plot(temp_k - 273.15, probe_y, "o", color=colors.get(name, "gray"),
                markersize=12, label=f"{labels.get(name, name)}: {temp_k - 273.15:.1f}C", zorder=5)

    # Heater zone
    ax.axhspan(heater_y, heater_y + heater_height, alpha=0.15, color="red", label="Heater zone")

    # Interface line
    ax.axhline(y=result.interface_height, color="green", linestyle="--",
               linewidth=1.5, alpha=0.7, label=f"Interface: {result.interface_height:.2f} m")

    ax.set_xlabel("Temperature [C]", fontsize=13)
    ax.set_ylabel("Height [m]", fontsize=13)
    ax.set_title("Sauna Temperature Stratification", fontsize=15)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(t_wall_c - 5, max(t_celsius) + 10)
    ax.set_ylim(0, room_y)
    fig.tight_layout()
    st.pyplot(fig)
    plt.close()

with col2:
    st.subheader("KPI Results")

    kpis = evaluate_phase1_kpis(result.probe_values)

    for kpi in kpis:
        status_icon = "PASS" if kpi.pass_fail == "pass" else ("FAIL" if kpi.pass_fail == "fail" else "-")
        if kpi.pass_fail == "pass":
            st.success(f"**{kpi.kpi_id}**: {kpi.name}")
        elif kpi.pass_fail == "fail":
            st.error(f"**{kpi.kpi_id}**: {kpi.name}")
        else:
            st.info(f"**{kpi.kpi_id}**: {kpi.name}")
        st.metric(label=f"{kpi.kpi_id} Value", value=f"{kpi.value:.2f} {kpi.unit}")

    st.divider()
    st.subheader("Probe Temperatures")
    for name, temp_k in result.probe_values.items():
        label = labels.get(name, name)
        st.metric(label=label, value=f"{temp_k - 273.15:.1f} C")

    st.divider()
    st.subheader("Two-Zone Model")
    st.write(f"**Interface height:** {result.interface_height:.2f} m")
    st.write(f"**Upper layer:** {result.upper_layer_temp - 273.15:.1f} C")
    st.write(f"**Lower layer:** {result.lower_layer_temp - 273.15:.1f} C")
    st.write(f"**Plume flow:** {result.plume_mass_flow:.3f} kg/s")

    st.divider()
    st.subheader("Solver Info")
    st.write(f"**Converged:** {'Yes' if result.converged else 'No'}")
    st.write(f"**Iterations:** {result.iterations}")
    st.write(f"**Final Residual:** {result.residual_history[-1]:.2e}")

# --- Convergence Plot ---
st.subheader("Convergence History")
fig2, ax2 = plt.subplots(figsize=(10, 3))
ax2.semilogy(result.residual_history, "b-", linewidth=1)
ax2.axhline(y=1e-6, color="r", linestyle="--", alpha=0.5, label="Tolerance")
ax2.set_xlabel("Iteration")
ax2.set_ylabel("Max Residual")
ax2.set_title("Solver Convergence")
ax2.legend()
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
st.pyplot(fig2)
plt.close()

# --- Room Cross-Section Visualization ---
st.subheader("Room Cross-Section (Y-Z plane)")
fig3, ax3 = plt.subplots(figsize=(8, 6))

# Create a 2D temperature field (approximate: same T across width)
nz_vis = 30
z_vis = np.linspace(0, room_z, nz_vis)
Y, Z = np.meshgrid(result.y_positions, z_vis, indexing="ij")
T_2d = np.tile(result.temperatures[:, np.newaxis], (1, nz_vis))

# Add some lateral variation (cooler near walls)
for j in range(nz_vis):
    wall_factor = 1.0 - 0.3 * (1.0 - np.sin(np.pi * z_vis[j] / room_z))
    T_2d[:, j] = t_wall_k + (T_2d[:, j] - t_wall_k) * wall_factor

contour = ax3.contourf(Z, Y, T_2d - 273.15, levels=20, cmap="RdYlBu_r")
plt.colorbar(contour, ax=ax3, label="Temperature [C]")

# Draw heater
heater_z0 = room_z / 2 - heater_width / 2
ax3.add_patch(plt.Rectangle((heater_z0, heater_y), heater_width, heater_height,
              linewidth=2, edgecolor="red", facecolor="red", alpha=0.5))
ax3.text(heater_z0 + heater_width / 2, heater_y + heater_height / 2, "Heater",
         ha="center", va="center", fontsize=9, color="white", fontweight="bold")

# Draw probes
for name, temp_k in result.probe_values.items():
    probe_data = next(p for p in case_data["probes"] if p["name"] == name)
    pz = probe_data["position"]["z"]
    py = probe_data["position"]["y"]
    ax3.plot(pz, py, "ko", markersize=8, zorder=5)
    ax3.annotate(f"{labels.get(name, name)}\n{temp_k - 273.15:.1f}C",
                 (pz, py), textcoords="offset points", xytext=(10, 5),
                 fontsize=9, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

ax3.set_xlabel("Depth (z) [m]", fontsize=12)
ax3.set_ylabel("Height (y) [m]", fontsize=12)
ax3.set_title("Temperature Distribution (cross-section at x = midplane)", fontsize=13)
ax3.set_xlim(0, room_z)
ax3.set_ylim(0, room_y)
fig3.tight_layout()
st.pyplot(fig3)
plt.close()
