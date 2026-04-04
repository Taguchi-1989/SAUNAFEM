"""SaunaFlow Visualization UI.

Run with: streamlit run ui/app.py
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import yaml

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

st.sidebar.subheader("Benches")
bench_depth = st.sidebar.slider("Bench depth [m]", 0.3, 1.0, 0.6, 0.05)
upper_bench_y = st.sidebar.slider("Upper bench height [m]", 0.5, room_y - 0.3, 1.2, 0.05)
lower_bench_y = st.sidebar.slider("Lower bench height [m]", 0.2, room_y - 0.3, 0.45, 0.05)
bench_wall_side = st.sidebar.selectbox("Bench wall", ["opposite (right)", "heater (left)"], index=0)

st.sidebar.subheader("Person")
person_height_cm = st.sidebar.slider("Person height [cm]", 140, 200, 170, 1)
person_seat = st.sidebar.selectbox("Sitting on", ["Upper bench", "Lower bench"])

st.sidebar.subheader("Environment")
t_wall_c = st.sidebar.slider("Wall temperature [C]", 10.0, 30.0, 20.0, 1.0)
t_wall_k = t_wall_c + 273.15

st.sidebar.subheader("Löyly (Steam)")
loyly_enabled = st.sidebar.checkbox("Enable löyly", value=False)
if loyly_enabled:
    water_ml = st.sidebar.slider("Water [mL]", 10, 500, 100, 10)
    tau_evap = st.sidebar.slider("Evaporation tau [s]", 1.0, 20.0, 5.0, 0.5)
else:
    water_ml = 0
    tau_evap = 5.0

st.sidebar.subheader("Aufguss")
aufguss_enabled = st.sidebar.checkbox("Enable aufguss", value=False)
if aufguss_enabled:
    beta_aug = st.sidebar.slider("Mixing coefficient β", 0.1, 2.0, 0.5, 0.1)
else:
    beta_aug = 0.0

# Probe positions follow bench positions
upper_probe_y = upper_bench_y + 0.6
lower_probe_y = lower_bench_y + 0.6

# Shared display constants
PROBE_LABELS = {"upper_bench": "Upper Bench", "lower_bench": "Lower Bench", "floor_level": "Floor Level"}
PROBE_COLORS = {"upper_bench": "#e74c3c", "lower_bench": "#f39c12", "floor_level": "#3498db"}
SEATED_HEIGHT_RATIO = 0.52  # seated eye-height as fraction of standing height

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
        "walls": {
            "temperature": t_wall_k,
            "type": "mixed",
            "model": "lumped",
            "thickness": 0.015,
            "conductivity": 0.12,
            "rho_cp": 500000,
        },
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
        {"name": "upper_bench", "position": {"x": room_x / 2, "y": upper_probe_y, "z": room_z / 2}, "fields": ["T"]},
        {"name": "lower_bench", "position": {"x": room_x / 2, "y": lower_probe_y, "z": room_z / 2}, "fields": ["T"]},
        {"name": "floor_level", "position": {"x": room_x / 2, "y": 0.1, "z": room_z / 2}, "fields": ["T"]},
    ],
}

if loyly_enabled:
    case_data["loyly"] = {"water_ml": water_ml, "time": 0.0, "tau_evap": tau_evap}
if aufguss_enabled:
    case_data["aufguss"] = {"beta_aug": beta_aug, "start_time": 0.0, "duration": 1000.0}

# Write temp YAML and solve
with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False, encoding="utf-8") as f:
    yaml.dump(case_data, f, default_flow_style=False)
    tmp_yaml = Path(f.name)

try:
    result = solve_two_zone(tmp_yaml, n_profile=80, max_iter=10000)
finally:
    tmp_yaml.unlink(missing_ok=True)


# ==========================================================
# SVG Room Cross-Section Diagram
# ==========================================================
def build_room_svg(
    room_w: float,
    room_h: float,
    heater_x0: float,
    heater_y0: float,
    heater_w: float,
    heater_h: float,
    heater_power_kw: float,
    upper_bench_h: float,
    lower_bench_h: float,
    bench_d: float,
    bench_side: str,
    person_h_cm: int,
    person_on: str,
    interface_h: float,
    t_upper_c: float,
    t_lower_c: float,
    t_wall_c: float,
    t_wall_inner_c: float = 0.0,
) -> str:
    """Generate an SVG cross-section of the sauna room.

    All physical dimensions are in metres; the SVG is scaled to fit a
    fixed viewport while preserving aspect ratio.
    """
    # --- SVG coordinate system ---
    # Physical origin: bottom-left of room interior.
    # SVG canvas: 800 x 600 with padding for labels.
    pad_l, pad_r, pad_t, pad_b = 80, 80, 50, 60
    canvas_w, canvas_h = 800, 600
    draw_w = canvas_w - pad_l - pad_r
    draw_h = canvas_h - pad_t - pad_b

    # Scale factor (m -> px), preserving aspect ratio
    sx = draw_w / room_w
    sy = draw_h / room_h
    scale = min(sx, sy)
    rw = room_w * scale
    rh = room_h * scale

    # Offset to centre the room drawing in the canvas
    ox = pad_l + (draw_w - rw) / 2
    oy = pad_t + (draw_h - rh) / 2

    def px(m_x: float) -> float:
        return ox + m_x * scale

    def py(m_y: float) -> float:
        return oy + rh - m_y * scale  # flip Y (SVG Y is top-down)

    # --- Colour helpers ---
    def temp_color(t: float) -> str:
        """Map temperature (C) to a colour: blue(20) -> red(100)."""
        frac = max(0.0, min(1.0, (t - 20) / 80))
        r = int(40 + 215 * frac)
        g = int(80 + 100 * (1 - abs(frac - 0.5) * 2))
        b = int(220 - 200 * frac)
        return f"#{r:02x}{g:02x}{b:02x}"

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {canvas_w} {canvas_h}" '
        f'width="100%" height="auto" '
        f'style="max-width:{canvas_w}px;font-family:sans-serif;">'
    )

    # --- All SVG defs in one block ---
    upper_col = temp_color(t_upper_c)
    lower_col = temp_color(t_lower_c)
    interface_frac = interface_h / room_h
    parts.append(
        "<defs>"
        f'<linearGradient id="thermalGrad" x1="0" y1="0" x2="0" y2="1">'
        f'<stop offset="0%" stop-color="{upper_col}" stop-opacity="0.45"/>'
        f'<stop offset="{(1 - interface_frac) * 100:.0f}%" stop-color="{upper_col}" stop-opacity="0.35"/>'
        f'<stop offset="{(1 - interface_frac) * 100 + 5:.0f}%" stop-color="{lower_col}" stop-opacity="0.25"/>'
        f'<stop offset="100%" stop-color="{lower_col}" stop-opacity="0.18"/>'
        f'</linearGradient>'
        '<marker id="arrowRed" markerWidth="6" markerHeight="6" '
        'refX="3" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="#e63946"/></marker>'
        '<marker id="dimArrowS" markerWidth="6" markerHeight="6" refX="1" refY="3" orient="auto">'
        '<path d="M6,0 L0,3 L6,6" fill="none" stroke="#666" stroke-width="1"/></marker>'
        '<marker id="dimArrowE" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto">'
        '<path d="M0,0 L6,3 L0,6" fill="none" stroke="#666" stroke-width="1"/></marker>'
        "</defs>"
    )

    # --- Room outline (filled with thermal gradient) ---
    rx0, ry0 = px(0), py(room_h)
    parts.append(
        f'<rect x="{rx0}" y="{ry0}" width="{rw}" height="{rh}" '
        f'fill="url(#thermalGrad)" stroke="#333" stroke-width="3" rx="2"/>'
    )

    # --- Wall labels ---
    fs = 11
    parts.append(f'<text x="{px(room_w / 2)}" y="{py(room_h) - 8}" text-anchor="middle" font-size="{fs}" fill="#555">ceiling</text>')
    parts.append(f'<text x="{px(room_w / 2)}" y="{py(0) + 16}" text-anchor="middle" font-size="{fs}" fill="#555">floor</text>')

    # --- Heater (on the left wall) ---
    hx = px(0)
    hy = py(heater_y0 + heater_h)
    hw_px = heater_w * scale * 0.6  # visual width (exaggerated a bit)
    hh_px = heater_h * scale
    parts.append(
        f'<rect x="{hx}" y="{hy}" width="{hw_px}" height="{hh_px}" '
        f'fill="#e63946" stroke="#b71c1c" stroke-width="2" rx="4"/>'
    )
    parts.append(
        f'<text x="{hx + hw_px / 2}" y="{hy + hh_px / 2 + 4}" '
        f'text-anchor="middle" font-size="11" fill="#fff" font-weight="bold">Stove</text>'
    )
    parts.append(
        f'<text x="{hx + hw_px / 2}" y="{hy + hh_px / 2 + 16}" '
        f'text-anchor="middle" font-size="9" fill="#fdd">{heater_power_kw:.0f} kW</text>'
    )

    # Heat arrows (rising from stove)
    arrow_x = hx + hw_px + 8
    for i in range(4):
        ay_base = hy - 10 - i * 24
        opacity = 0.7 - i * 0.12
        parts.append(
            f'<path d="M{arrow_x},{ay_base} C{arrow_x + 6},{ay_base - 10} '
            f'{arrow_x - 6},{ay_base - 18} {arrow_x},{ay_base - 26}" '
            f'stroke="#e63946" stroke-width="2" fill="none" opacity="{opacity:.2f}" '
            f'marker-end="url(#arrowRed)"/>'
        )
    # --- Benches ---
    bench_on_right = "opposite" in bench_side
    if bench_on_right:
        bx_start = px(room_w) - bench_d * scale
        bx_label = px(room_w) - bench_d * scale / 2
    else:
        bx_start = px(0)
        bx_label = px(0) + bench_d * scale / 2

    bench_thick = 6
    bd_px = bench_d * scale

    def _draw_bench(
        label: str, bench_h_m: float, leg_target_m: float, fill: str,
    ) -> float:
        """Draw a bench with legs and label. Returns bench Y in SVG coords."""
        by = py(bench_h_m)
        parts.append(
            f'<rect x="{bx_start}" y="{by}" width="{bd_px}" height="{bench_thick}" '
            f'fill="{fill}" stroke="#5d4037" stroke-width="1.5" rx="2"/>'
        )
        for leg_off in [0.08, 0.92]:
            lx = bx_start + bd_px * leg_off
            parts.append(
                f'<line x1="{lx}" y1="{by + bench_thick}" x2="{lx}" y2="{py(leg_target_m)}" '
                f'stroke="#5d4037" stroke-width="2"/>'
            )
        parts.append(
            f'<text x="{bx_label}" y="{by - 6}" text-anchor="middle" '
            f'font-size="10" fill="#5d4037" font-weight="bold">{label}</text>'
        )
        return by

    ub_y = _draw_bench("Upper bench", upper_bench_h, lower_bench_h, "#8d6e63")
    lb_y = _draw_bench("Lower bench", lower_bench_h, 0.0, "#a1887f")

    # --- Person (stick figure, sitting) ---
    person_h_m = person_h_cm / 100.0
    sitting_h = person_h_m * 0.52  # seated height ~ 52% of standing
    torso_h = sitting_h * 0.65
    head_r_m = 0.10  # head radius in m

    if person_on == "Upper bench":
        seat_y_m = upper_bench_h
    else:
        seat_y_m = lower_bench_h

    person_cx = bx_start + bd_px * 0.5

    # Torso
    torso_bot = py(seat_y_m)
    torso_top = py(seat_y_m + torso_h)
    parts.append(
        f'<line x1="{person_cx}" y1="{torso_bot}" x2="{person_cx}" y2="{torso_top}" '
        f'stroke="#1565c0" stroke-width="3" stroke-linecap="round"/>'
    )

    # Head
    head_cy = py(seat_y_m + torso_h + head_r_m)
    head_r_px = head_r_m * scale
    parts.append(
        f'<circle cx="{person_cx}" cy="{head_cy}" r="{head_r_px}" '
        f'fill="#bbdefb" stroke="#1565c0" stroke-width="2"/>'
    )

    # Arms (slightly out)
    arm_y = py(seat_y_m + torso_h * 0.65)
    arm_len = 0.22 * scale
    parts.append(
        f'<line x1="{person_cx - arm_len}" y1="{arm_y + 6}" '
        f'x2="{person_cx + arm_len}" y2="{arm_y + 6}" '
        f'stroke="#1565c0" stroke-width="2.5" stroke-linecap="round"/>'
    )

    # Legs (bent, sitting)
    knee_y = py(seat_y_m)
    foot_y = py(seat_y_m - 0.42)  # lower legs hang or rest
    leg_spread = 0.08 * scale
    for sign in [-1, 1]:
        kx = person_cx + sign * leg_spread
        fx = person_cx + sign * leg_spread * 2.5
        parts.append(
            f'<line x1="{kx}" y1="{knee_y}" x2="{fx}" y2="{foot_y}" '
            f'stroke="#1565c0" stroke-width="2.5" stroke-linecap="round"/>'
        )

    # Person height label
    person_label_x = person_cx + 20
    parts.append(
        f'<text x="{person_label_x}" y="{head_cy - head_r_px - 4}" '
        f'font-size="10" fill="#1565c0">{person_h_cm}cm</text>'
    )

    # --- Interface height (dashed line) ---
    iy = py(interface_h)
    parts.append(
        f'<line x1="{px(0)}" y1="{iy}" x2="{px(room_w)}" y2="{iy}" '
        f'stroke="#2e7d32" stroke-width="1.5" stroke-dasharray="8,4"/>'
    )
    parts.append(
        f'<text x="{px(room_w) + 4}" y="{iy + 4}" font-size="10" fill="#2e7d32">'
        f'Interface {interface_h:.2f}m</text>'
    )

    # --- Temperature labels in zones ---
    upper_mid_y = py((interface_h + room_h) / 2)
    lower_mid_y = py(interface_h / 2)
    parts.append(
        f'<text x="{px(room_w / 2)}" y="{upper_mid_y}" text-anchor="middle" '
        f'font-size="14" fill="{temp_color(t_upper_c)}" font-weight="bold" opacity="0.85">'
        f'Upper: {t_upper_c:.1f} &#176;C</text>'
    )
    parts.append(
        f'<text x="{px(room_w / 2)}" y="{lower_mid_y}" text-anchor="middle" '
        f'font-size="14" fill="{temp_color(t_lower_c)}" font-weight="bold" opacity="0.85">'
        f'Lower: {t_lower_c:.1f} &#176;C</text>'
    )

    # --- Dimension arrows ---
    dim_col = "#666"

    # Room height (right side)
    dx_r = px(room_w) + 28
    parts.append(
        f'<line x1="{dx_r}" y1="{py(0)}" x2="{dx_r}" y2="{py(room_h)}" '
        f'stroke="{dim_col}" stroke-width="1" marker-start="url(#dimArrowS)" marker-end="url(#dimArrowE)"/>'
    )
    parts.append(
        f'<text x="{dx_r + 4}" y="{py(room_h / 2) + 4}" font-size="11" fill="{dim_col}" '
        f'writing-mode="tb">{room_y:.1f}m</text>'
    )

    # Room width (bottom)
    dy_b = py(0) + 28
    parts.append(
        f'<line x1="{px(0)}" y1="{dy_b}" x2="{px(room_w)}" y2="{dy_b}" '
        f'stroke="{dim_col}" stroke-width="1" marker-start="url(#dimArrowS)" marker-end="url(#dimArrowE)"/>'
    )
    parts.append(
        f'<text x="{px(room_w / 2)}" y="{dy_b + 16}" text-anchor="middle" '
        f'font-size="11" fill="{dim_col}">{room_x:.1f}m (width) / {room_z:.1f}m (depth)</text>'
    )

    # Heater height dimension (left side)
    dx_l = px(0) - 16
    parts.append(
        f'<line x1="{dx_l}" y1="{py(heater_y0)}" x2="{dx_l}" y2="{py(heater_y0 + heater_h)}" '
        f'stroke="#e63946" stroke-width="1" marker-start="url(#dimArrowS)" marker-end="url(#dimArrowE)"/>'
    )
    parts.append(
        f'<text x="{dx_l - 4}" y="{py(heater_y0 + heater_h / 2) + 4}" '
        f'text-anchor="end" font-size="9" fill="#e63946">{heater_h:.1f}m</text>'
    )

    # Bench-to-stove distance
    if bench_on_right:
        dist_stove = room_w - bench_d
    else:
        dist_stove = 0.0  # bench on same wall as stove
    if dist_stove > 0.1:
        dist_y_px = py(upper_bench_h) + 20
        parts.append(
            f'<line x1="{px(0) + hw_px + 2}" y1="{dist_y_px}" '
            f'x2="{bx_start - 2}" y2="{dist_y_px}" '
            f'stroke="#ff9800" stroke-width="1" stroke-dasharray="4,3" '
            f'marker-start="url(#dimArrowS)" marker-end="url(#dimArrowE)"/>'
        )
        parts.append(
            f'<text x="{(px(0) + hw_px + bx_start) / 2}" y="{dist_y_px - 4}" '
            f'text-anchor="middle" font-size="10" fill="#ff9800">{dist_stove:.1f}m</text>'
        )

    # Upper bench height dimension
    ubdx = bx_start - 12
    parts.append(
        f'<line x1="{ubdx}" y1="{py(0)}" x2="{ubdx}" y2="{py(upper_bench_h)}" '
        f'stroke="#5d4037" stroke-width="1" stroke-dasharray="3,3"/>'
    )
    parts.append(
        f'<text x="{ubdx - 3}" y="{py(upper_bench_h / 2) + 3}" '
        f'text-anchor="end" font-size="9" fill="#5d4037">{upper_bench_h:.2f}m</text>'
    )

    # --- Wall temperature label ---
    parts.append(
        f'<text x="{px(0) - 4}" y="{py(room_h / 2)}" text-anchor="end" '
        f'font-size="9" fill="#888" transform="rotate(-90,{px(0) - 4},{py(room_h / 2)})">'
        f'T_wall: {t_wall_c:.0f}&#8594;{t_wall_inner_c:.0f}&#176;C</text>'
    )

    parts.append("</svg>")
    return "\n".join(parts)


# ==========================================================
# Main layout
# ==========================================================
tab_diagram, tab_profile, tab_contour, tab_comfort = st.tabs(
    ["Room Diagram", "Temperature Profile", "Cross-Section Contour", "Thermal Comfort"]
)

# --- Tab 1: Room Diagram (SVG) ---
with tab_diagram:
    st.subheader("Sauna Room Cross-Section")

    svg = build_room_svg(
        room_w=room_x,
        room_h=room_y,
        heater_x0=0.0,
        heater_y0=heater_y,
        heater_w=heater_width,
        heater_h=heater_height,
        heater_power_kw=heater_power,
        upper_bench_h=upper_bench_y,
        lower_bench_h=lower_bench_y,
        bench_d=bench_depth,
        bench_side=bench_wall_side,
        person_h_cm=person_height_cm,
        person_on=person_seat,
        interface_h=result.interface_height,
        t_upper_c=result.upper_layer_temp - 273.15,
        t_lower_c=result.lower_layer_temp - 273.15,
        t_wall_c=t_wall_c,
        t_wall_inner_c=result.wall_inner_temp - 273.15,
    )

    st.markdown(svg, unsafe_allow_html=True)

    # Key metrics below the diagram
    col_a, col_b, col_c, col_d, col_e = st.columns(5)
    with col_a:
        st.metric("Upper Layer", f"{result.upper_layer_temp - 273.15:.1f} °C")
    with col_b:
        st.metric("Lower Layer", f"{result.lower_layer_temp - 273.15:.1f} °C")
    with col_c:
        st.metric("Wall Surface", f"{result.wall_inner_temp - 273.15:.1f} °C")
    with col_d:
        st.metric("Interface", f"{result.interface_height:.2f} m")
    with col_e:
        st.metric("Humidity", f"{result.humidity_ratio * 1000:.1f} g/kg")

    seat_base = upper_bench_y if person_seat == "Upper bench" else lower_bench_y
    head_y = seat_base + person_height_cm / 100.0 * SEATED_HEIGHT_RATIO

    # Interpolate temperature at head height
    idx = int(np.clip(head_y / (room_y / 80), 0, 79))
    head_temp = result.temperatures[idx] - 273.15

    st.info(
        f"**Head position:** {head_y:.2f} m from floor "
        f"({'above' if head_y > result.interface_height else 'below'} interface) "
        f"  |  **Temperature at head:** {head_temp:.1f} C"
    )

# --- Tab 2: Temperature Profile ---
with tab_profile:
    st.subheader("Vertical Temperature Profile")

    fig, ax = plt.subplots(figsize=(8, 6))

    t_celsius = result.temperatures - 273.15
    ax.plot(t_celsius, result.y_positions, "b-", linewidth=2.5, label="Temperature")

    for name, temp_k in result.probe_values.items():
        probe_y = next(p["position"]["y"] for p in case_data["probes"] if p["name"] == name)
        ax.plot(temp_k - 273.15, probe_y, "o", color=PROBE_COLORS.get(name, "gray"),
                markersize=12, label=f"{PROBE_LABELS.get(name, name)}: {temp_k - 273.15:.1f}C", zorder=5)

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

# --- Tab 3: Cross-Section Contour ---
with tab_contour:
    st.subheader("Room Cross-Section (Y-Z plane)")
    fig3, ax3 = plt.subplots(figsize=(8, 6))

    nz_vis = 30
    z_vis = np.linspace(0, room_z, nz_vis)
    Y, Z = np.meshgrid(result.y_positions, z_vis, indexing="ij")
    T_2d = np.tile(result.temperatures[:, np.newaxis], (1, nz_vis))

    wall_factors = 1.0 - 0.3 * (1.0 - np.sin(np.pi * z_vis / room_z))
    T_2d = t_wall_k + (T_2d - t_wall_k) * wall_factors[np.newaxis, :]

    contour = ax3.contourf(Z, Y, T_2d - 273.15, levels=20, cmap="RdYlBu_r")
    plt.colorbar(contour, ax=ax3, label="Temperature [C]")

    heater_z0 = room_z / 2 - heater_width / 2
    ax3.add_patch(plt.Rectangle((heater_z0, heater_y), heater_width, heater_height,
                  linewidth=2, edgecolor="red", facecolor="red", alpha=0.5))
    ax3.text(heater_z0 + heater_width / 2, heater_y + heater_height / 2, "Heater",
             ha="center", va="center", fontsize=9, color="white", fontweight="bold")

    for name, temp_k in result.probe_values.items():
        probe_data = next(p for p in case_data["probes"] if p["name"] == name)
        pz = probe_data["position"]["z"]
        p_y = probe_data["position"]["y"]
        ax3.plot(pz, p_y, "ko", markersize=8, zorder=5)
        ax3.annotate(f"{PROBE_LABELS.get(name, name)}\n{temp_k - 273.15:.1f}C",
                     (pz, p_y), textcoords="offset points", xytext=(10, 5),
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

# --- Tab 4: Thermal Comfort ---
with tab_comfort:
    st.subheader("Thermal Comfort Analysis")

    col_l, col_r = st.columns(2)

    with col_l:
        st.markdown("### Temperature Summary")
        st.metric("Upper Layer (dry-bulb)", f"{result.upper_layer_temp - 273.15:.1f} °C")
        st.metric("Lower Layer (dry-bulb)", f"{result.lower_layer_temp - 273.15:.1f} °C")
        st.metric("Wall Surface", f"{result.wall_inner_temp - 273.15:.1f} °C")

    with col_r:
        st.markdown("### Humidity & Perceived Temperature")
        st.metric("Absolute Humidity", f"{result.humidity_ratio * 1000:.1f} g/kg")
        st.metric("Relative Humidity", f"{result.relative_humidity:.1%}")
        st.metric("Perceived Temp (Upper)", f"{result.perceived_temp_upper:.1f} °C")
        st.metric("Perceived Temp (Lower)", f"{result.perceived_temp_lower:.1f} °C")

    # Thermal stress gauge
    perceived = result.perceived_temp_upper
    if perceived < 60:
        level, color = "Comfortable", "green"
    elif perceived < 80:
        level, color = "Moderate", "orange"
    elif perceived < 100:
        level, color = "Intense", "red"
    else:
        level, color = "Extreme", "darkred"

    st.markdown(f"### Thermal Stress Level: :{color}[**{level}**] ({perceived:.0f}°C perceived)")

# --- KPI & Solver Info (below tabs) ---
st.divider()
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("KPI Results")
    kpis = evaluate_phase1_kpis(result.probe_values)
    for kpi in kpis:
        if kpi.pass_fail == "pass":
            st.success(f"**{kpi.kpi_id}**: {kpi.name}")
        elif kpi.pass_fail == "fail":
            st.error(f"**{kpi.kpi_id}**: {kpi.name}")
        else:
            st.info(f"**{kpi.kpi_id}**: {kpi.name}")
        st.metric(label=f"{kpi.kpi_id} Value", value=f"{kpi.value:.2f} {kpi.unit}")

with col2:
    st.subheader("Solver Convergence")
    fig2, ax2 = plt.subplots(figsize=(6, 3))
    ax2.semilogy(result.residual_history, "b-", linewidth=1)
    ax2.axhline(y=1e-4, color="r", linestyle="--", alpha=0.5, label="Tolerance")
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Max Residual")
    ax2.set_title(f"{'Converged' if result.converged else 'NOT converged'} ({result.iterations} iter)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    st.pyplot(fig2)
    plt.close()
