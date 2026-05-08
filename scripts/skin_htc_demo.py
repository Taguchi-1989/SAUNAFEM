"""Sweep skin overall heat transfer coefficient over (T_air, RH, V).

Standalone driver for `harness.skin_htc`. Produces a CSV table and a
4-panel plot characterising "perceived heat" as a function of dry-bulb
temperature, relative humidity, and local air velocity.

Usage (from repo root):
    PYTHONPATH=src python scripts/skin_htc_demo.py
    PYTHONPATH=src python scripts/skin_htc_demo.py --t-mrt 110 --output-dir results/skin_htc

Outputs (under --output-dir, default results/skin_htc/):
    sweep.csv               flat table of every (T, RH, V) combination
    u_overall_vs_T_V.png    heat map: U_overall(T_air, V) at fixed RH
    q_total_vs_T_V.png      heat map: q_total(T_air, V) at fixed RH
    components_vs_V.png     line plot: q_conv, q_rad, q_evap vs V
    u_vs_RH.png             line plot: U_overall vs RH at fixed (T, V)
    sweat_q_total_vs_RH.png line plot: q_total vs RH for several w_skin
    sweat_q_evap_vs_RH.png  line plot: q_evap vs RH for several w_skin
    sweat_dry_vs_wet_bars.png  stacked bars: components for dry vs wet skin
    per_part_bars.png       per-body-part q_total in canonical scenarios
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from harness.skin_htc import BODY_PARTS, compute_skin_balance, sweep_grid

# Default sweep ranges spanning typical Finnish/Aufguss sauna conditions
DEFAULT_T_AIR = [60.0, 70.0, 80.0, 90.0, 100.0, 110.0]    # °C
DEFAULT_RH = [0.05, 0.10, 0.20, 0.30, 0.50, 0.80]          # -
DEFAULT_V = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]            # m/s


def write_csv(rows: list[dict], path: Path) -> None:
    """Write sweep rows as CSV."""
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {path} ({len(rows)} rows)")


def plot_u_overall_heatmap(
    t_air: list[float],
    v_values: list[float],
    rh_fixed: float,
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """Heat map of whole-body U_overall as a function of (T_air, V) at fixed RH."""
    grid = np.zeros((len(v_values), len(t_air)))
    for i, v in enumerate(v_values):
        for j, t in enumerate(t_air):
            bal = compute_skin_balance(t_air_c=t, rh=rh_fixed, v_local=v, t_mrt_c=t_mrt_c)
            grid[i, j] = bal.whole_body.u_overall

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(
        grid,
        origin="lower",
        aspect="auto",
        extent=[min(t_air), max(t_air), min(v_values), max(v_values)],
        cmap="inferno",
    )
    ax.set_xlabel("T_air [°C]")
    ax.set_ylabel("V_local [m/s]")
    ax.set_title(f"Whole-body U_overall  (RH={rh_fixed:.0%})")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("U_overall [W/(m²·K)]")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_q_total_heatmap(
    t_air: list[float],
    v_values: list[float],
    rh_fixed: float,
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """Heat map of whole-body q_total as a function of (T_air, V) at fixed RH."""
    grid = np.zeros((len(v_values), len(t_air)))
    for i, v in enumerate(v_values):
        for j, t in enumerate(t_air):
            bal = compute_skin_balance(t_air_c=t, rh=rh_fixed, v_local=v, t_mrt_c=t_mrt_c)
            grid[i, j] = bal.whole_body.q_total

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(
        grid,
        origin="lower",
        aspect="auto",
        extent=[min(t_air), max(t_air), min(v_values), max(v_values)],
        cmap="hot",
    )
    ax.set_xlabel("T_air [°C]")
    ax.set_ylabel("V_local [m/s]")
    ax.set_title(f"Whole-body q_total  (RH={rh_fixed:.0%})")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("q_total [W/m²]")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_components_vs_velocity(
    t_air_fixed: float,
    rh_fixed: float,
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """Line plot of q_conv, q_rad, q_evap, q_total vs V at fixed (T, RH)."""
    v_values = np.linspace(0.05, 3.0, 60)
    q_conv, q_rad, q_evap, q_total = [], [], [], []
    for v in v_values:
        bal = compute_skin_balance(
            t_air_c=t_air_fixed, rh=rh_fixed, v_local=float(v), t_mrt_c=t_mrt_c,
        )
        wb = bal.whole_body
        q_conv.append(wb.q_conv)
        q_rad.append(wb.q_rad)
        q_evap.append(wb.q_evap)
        q_total.append(wb.q_total)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(v_values, q_conv, label="q_conv (convection)")
    ax.plot(v_values, q_rad, label="q_rad (radiation)", linestyle="--")
    ax.plot(v_values, q_evap, label="q_evap (evap/cond)", linestyle=":")
    ax.plot(v_values, q_total, label="q_total", color="black", linewidth=2)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("V_local [m/s]")
    ax.set_ylabel("Heat flux [W/m²]")
    ax.set_title(
        f"Skin heat flux components  (T_air={t_air_fixed:.0f}°C, RH={rh_fixed:.0%})"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_u_vs_rh(
    t_air_fixed: float,
    v_values: list[float],
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """U_overall vs RH at several velocities, fixed T_air."""
    rh_values = np.linspace(0.0, 0.95, 50)
    fig, ax = plt.subplots(figsize=(8, 5))
    for v in v_values:
        u_arr = []
        for rh in rh_values:
            bal = compute_skin_balance(
                t_air_c=t_air_fixed, rh=float(rh), v_local=v, t_mrt_c=t_mrt_c,
            )
            u_arr.append(bal.whole_body.u_overall)
        ax.plot(rh_values, u_arr, label=f"V = {v:.1f} m/s")
    ax.set_xlabel("Relative humidity [-]")
    ax.set_ylabel("U_overall [W/(m²·K)]")
    ax.set_title(f"Whole-body U_overall vs RH  (T_air={t_air_fixed:.0f}°C)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def _crossover_rh(t_air_c: float, t_skin_c: float = 36.0) -> float:
    """RH at which P_vapor = P_sat,skin (transition from sweat-cool → condensation)."""
    from harness.skin_htc import saturation_pressure_kpa
    p_sat_air = saturation_pressure_kpa(t_air_c)
    p_sat_skin = saturation_pressure_kpa(t_skin_c)
    rh = p_sat_skin / p_sat_air if p_sat_air > 0 else 1.0
    return float(min(rh, 1.0))


def plot_sweat_q_total_vs_rh(
    t_air_fixed: float,
    v_fixed: float,
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """q_total vs RH for several skin wettedness values.

    Two-panel layout:
      Top:    full RH range (0-0.95) — shows condensation-dominated regime
      Bottom: zoomed to RH < crossover — where sweating actually matters
    """
    w_skin_values = [0.0, 0.2, 0.4, 0.7, 1.0]
    rh_full = np.linspace(0.0, 0.95, 80)
    rh_cross = _crossover_rh(t_air_fixed)
    rh_zoom = np.linspace(0.0, max(rh_cross * 1.4, 0.15), 80)

    def q_total_curve(w: float, rh_array: np.ndarray) -> list[float]:
        out = []
        for rh in rh_array:
            bal = compute_skin_balance(
                t_air_c=t_air_fixed, rh=float(rh), v_local=v_fixed,
                t_mrt_c=t_mrt_c, w_skin=w,
            )
            out.append(bal.whole_body.q_total)
        return out

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(9, 8))
    for w in w_skin_values:
        label = f"w_skin = {w:.1f}"
        if w == 0.0:
            label += " (no sweat)"
        elif w == 1.0:
            label += " (max sweat)"
        ax_top.plot(rh_full, q_total_curve(w, rh_full), label=label, linewidth=1.8)
        ax_bot.plot(rh_zoom, q_total_curve(w, rh_zoom), label=label, linewidth=1.8)

    for ax in (ax_top, ax_bot):
        ax.axvline(rh_cross, color="gray", linewidth=0.7, linestyle="--",
                   label=f"crossover RH = {rh_cross:.3f}")
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.set_xlabel("Relative humidity [-]")
        ax.set_ylabel("q_total [W/m²]")
        ax.grid(True, alpha=0.3)
    ax_top.set_title(
        f"Sweating effect on q_total  "
        f"(T_air={t_air_fixed:.0f}°C, V={v_fixed:.1f} m/s)\n"
        f"top: full RH range — bottom: zoom to sweating-relevant region"
    )
    ax_top.legend(fontsize=8, loc="upper left")
    ax_bot.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_sweat_q_evap_vs_rh(
    t_air_fixed: float,
    v_fixed: float,
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """q_evap vs RH for several w_skin.

    Two-panel zoom: full range on top, crossover region on bottom.
    Shows how sweating only affects the negative (cooling) branch — the
    condensation branch is wettedness-independent.
    """
    w_skin_values = [0.0, 0.2, 0.4, 0.7, 1.0]
    rh_full = np.linspace(0.0, 0.95, 80)
    rh_cross = _crossover_rh(t_air_fixed)
    rh_zoom = np.linspace(0.0, max(rh_cross * 1.4, 0.15), 80)

    def q_evap_curve(w: float, rh_array: np.ndarray) -> list[float]:
        out = []
        for rh in rh_array:
            bal = compute_skin_balance(
                t_air_c=t_air_fixed, rh=float(rh), v_local=v_fixed,
                t_mrt_c=t_mrt_c, w_skin=w,
            )
            out.append(bal.whole_body.q_evap)
        return out

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(9, 8))
    for w in w_skin_values:
        ax_top.plot(rh_full, q_evap_curve(w, rh_full), label=f"w_skin = {w:.1f}",
                    linewidth=1.8)
        ax_bot.plot(rh_zoom, q_evap_curve(w, rh_zoom), label=f"w_skin = {w:.1f}",
                    linewidth=1.8)

    for ax in (ax_top, ax_bot):
        ax.axvline(rh_cross, color="gray", linewidth=0.7, linestyle="--",
                   label=f"crossover RH = {rh_cross:.3f}")
        ax.axhline(0, color="gray", linewidth=0.4, linestyle="--")
        ax.set_xlabel("Relative humidity [-]")
        ax.set_ylabel("q_evap [W/m²]")
        ax.grid(True, alpha=0.3)

    ax_top.set_title(
        f"Evaporative/condensation flux vs RH and sweating  "
        f"(T_air={t_air_fixed:.0f}°C, V={v_fixed:.1f} m/s)\n"
        f"negative = sweat cooling, positive = vapor condensing on skin"
    )
    ax_top.legend(fontsize=8, loc="upper left")
    ax_bot.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_dry_vs_wet_bars(
    t_mrt_c: float | None,
    out_path: Path,
) -> None:
    """Side-by-side stacked bars of (q_conv, q_rad, q_evap) for dry vs sweating skin
    across canonical sauna scenarios.
    """
    scenarios = [
        ("Mild dry\ncalm",        70.0, 0.10, 0.10),
        ("Hot dry\ncalm",         90.0, 0.05, 0.10),
        ("Hot dry\nAufguss",      90.0, 0.05, 2.00),
        ("Löyly humid\ncalm",     85.0, 0.50, 0.20),
        ("Löyly humid\nAufguss",  85.0, 0.50, 2.00),
    ]
    n = len(scenarios)
    width = 0.35
    x = np.arange(n)

    def collect(w_skin: float) -> tuple[list[float], list[float], list[float]]:
        qc, qr, qe = [], [], []
        for _, t_air, rh, v in scenarios:
            bal = compute_skin_balance(
                t_air_c=t_air, rh=rh, v_local=v, t_mrt_c=t_mrt_c, w_skin=w_skin,
            )
            qc.append(bal.whole_body.q_conv)
            qr.append(bal.whole_body.q_rad)
            qe.append(bal.whole_body.q_evap)
        return qc, qr, qe

    qc_dry, qr_dry, qe_dry = collect(0.0)
    qc_wet, qr_wet, qe_wet = collect(1.0)

    fig, ax = plt.subplots(figsize=(10, 5))
    # Stack: positive components stacked up from 0; negative stacked down
    def stack(ax, x, qc, qr, qe, width, label_prefix, hatch):
        bottom_pos = np.zeros(len(qc))
        bottom_neg = np.zeros(len(qc))
        for vals, color, lbl in [
            (qc, "tab:red", "q_conv"),
            (qr, "tab:orange", "q_rad"),
            (qe, "tab:blue", "q_evap"),
        ]:
            arr = np.array(vals)
            pos = np.where(arr > 0, arr, 0)
            neg = np.where(arr < 0, arr, 0)
            ax.bar(
                x, pos, width, bottom=bottom_pos, color=color, alpha=0.85,
                hatch=hatch, edgecolor="black", linewidth=0.4,
                label=f"{label_prefix} {lbl}" if hatch == "" else None,
            )
            ax.bar(
                x, neg, width, bottom=bottom_neg, color=color, alpha=0.85,
                hatch=hatch, edgecolor="black", linewidth=0.4,
            )
            bottom_pos += pos
            bottom_neg += neg

    stack(ax, x - width / 2, qc_dry, qr_dry, qe_dry, width, "dry", "")
    stack(ax, x + width / 2, qc_wet, qr_wet, qe_wet, width, "wet", "//")

    # Total markers
    totals_dry = [a + b + c for a, b, c in zip(qc_dry, qr_dry, qe_dry, strict=True)]
    totals_wet = [a + b + c for a, b, c in zip(qc_wet, qr_wet, qe_wet, strict=True)]
    ax.scatter(x - width / 2, totals_dry, marker="D", color="black", s=40,
               zorder=5, label="q_total (dry, w=0)")
    ax.scatter(x + width / 2, totals_wet, marker="o", color="black", s=40,
               zorder=5, label="q_total (wet, w=1)")

    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([s[0] for s in scenarios], fontsize=9)
    ax.set_ylabel("Heat flux [W/m²]")
    ax.set_title(
        "Skin heat flux components: dry skin (left, w=0) vs sweating skin (right, w=1, hatched)"
    )
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_per_part_bars(
    t_mrt_c: float | None,
    w_skin: float,
    out_path: Path,
) -> None:
    """Bar chart of q_total per body part across scenarios."""
    scenarios = [
        ("Hot dry\ncalm",         90.0, 0.05, 0.10),
        ("Hot dry\nAufguss",      90.0, 0.05, 2.00),
        ("Löyly humid\ncalm",     85.0, 0.50, 0.20),
        ("Löyly humid\nAufguss",  85.0, 0.50, 2.00),
    ]
    parts = ("face", "chest", "back", "arm", "thigh", "calf")
    width = 0.13
    x = np.arange(len(scenarios))

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, part in enumerate(parts):
        vals = []
        for _, t_air, rh, v in scenarios:
            bal = compute_skin_balance(
                t_air_c=t_air, rh=rh, v_local=v, t_mrt_c=t_mrt_c, w_skin=w_skin,
            )
            vals.append(bal.parts[part].q_total)
        ax.bar(x + (i - 2.5) * width, vals, width, label=part)
    ax.set_xticks(x)
    ax.set_xticklabels([s[0] for s in scenarios], fontsize=9)
    ax.set_ylabel("q_total [W/m²]")
    ax.set_title(f"Per-body-part heat flux (w_skin = {w_skin:.1f})")
    ax.legend(ncol=6, fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def print_summary_table(t_mrt_c: float | None) -> None:
    """Print a small comparison table for a few canonical sauna scenarios."""
    scenarios = [
        ("Mild dry, calm",      70.0, 0.10, 0.10),
        ("Hot dry, calm",       90.0, 0.05, 0.10),
        ("Hot dry, Aufguss",    90.0, 0.05, 2.00),
        ("Löyly humid, calm",   85.0, 0.50, 0.20),
        ("Löyly humid, Aufguss",85.0, 0.50, 2.00),
        ("Extreme",            105.0, 0.30, 2.50),
    ]
    print()
    print("Scenario summary (whole body):")
    header = (
        f"  {'name':<22}{'T':>6}{'RH':>6}{'V':>6}"
        f"{'h_c':>8}{'h_r':>8}{'q_c':>9}{'q_r':>9}{'q_e':>9}{'q_tot':>9}{'U_ov':>9}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    for name, t_air, rh, v in scenarios:
        bal = compute_skin_balance(t_air_c=t_air, rh=rh, v_local=v, t_mrt_c=t_mrt_c)
        wb = bal.whole_body
        u_str = f"{wb.u_overall:>9.2f}" if not np.isnan(wb.u_overall) else "      nan"
        print(
            f"  {name:<22}{t_air:>6.0f}{rh:>6.2f}{v:>6.2f}"
            f"{wb.h_conv:>8.2f}{wb.h_rad:>8.2f}"
            f"{wb.q_conv:>9.0f}{wb.q_rad:>9.0f}{wb.q_evap:>9.0f}"
            f"{wb.q_total:>9.0f}{u_str}"
        )
    print()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--t-mrt", type=float, default=None,
                   help="Mean radiant temperature [°C]. Default: T_air per row.")
    p.add_argument("--t-skin", type=float, default=36.0, help="Skin temperature [°C].")
    p.add_argument("--w-skin", type=float, default=0.4, help="Skin wettedness [0-1].")
    p.add_argument("--output-dir", type=Path, default=Path("results/skin_htc"))
    p.add_argument("--rh-for-heatmap", type=float, default=0.30,
                   help="RH at which to render U/q heat maps.")
    p.add_argument("--t-for-components", type=float, default=90.0,
                   help="T_air [°C] for components-vs-V plot.")
    p.add_argument("--rh-for-components", type=float, default=0.30,
                   help="RH for components-vs-V plot.")
    p.add_argument("--t-for-rh-curve", type=float, default=85.0,
                   help="T_air [°C] for U-vs-RH plot.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Writing outputs to {out_dir}/")

    rows = sweep_grid(
        t_air_c_values=DEFAULT_T_AIR,
        rh_values=DEFAULT_RH,
        v_values=DEFAULT_V,
        t_mrt_c=args.t_mrt,
        t_skin_c=args.t_skin,
        w_skin=args.w_skin,
        parts=BODY_PARTS,
    )
    write_csv(rows, out_dir / "sweep.csv")

    plot_u_overall_heatmap(
        DEFAULT_T_AIR, DEFAULT_V, args.rh_for_heatmap, args.t_mrt,
        out_dir / "u_overall_vs_T_V.png",
    )
    plot_q_total_heatmap(
        DEFAULT_T_AIR, DEFAULT_V, args.rh_for_heatmap, args.t_mrt,
        out_dir / "q_total_vs_T_V.png",
    )
    plot_components_vs_velocity(
        args.t_for_components, args.rh_for_components, args.t_mrt,
        out_dir / "components_vs_V.png",
    )
    plot_u_vs_rh(
        args.t_for_rh_curve, [0.1, 0.5, 1.0, 2.0], args.t_mrt,
        out_dir / "u_vs_RH.png",
    )
    # Sweating-effect plots
    plot_sweat_q_total_vs_rh(
        args.t_for_rh_curve, 0.5, args.t_mrt,
        out_dir / "sweat_q_total_vs_RH.png",
    )
    plot_sweat_q_evap_vs_rh(
        args.t_for_rh_curve, 0.5, args.t_mrt,
        out_dir / "sweat_q_evap_vs_RH.png",
    )
    plot_dry_vs_wet_bars(args.t_mrt, out_dir / "sweat_dry_vs_wet_bars.png")
    plot_per_part_bars(args.t_mrt, args.w_skin, out_dir / "per_part_bars.png")

    print_summary_table(args.t_mrt)


if __name__ == "__main__":
    main()
