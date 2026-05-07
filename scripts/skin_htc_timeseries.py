"""Time-series sauna events: q(t), T_eq(t), and cumulative Q.

Real-world sauna intensity isn't steady — it's a sequence of pulses
(löyly throws, Aufguss towel waves) on top of a varying background.
This script defines named time-evolving (T_air, RH, V) profiles and
integrates the skin heat balance through them.

Key quantity: Q [J/m²] = ∫ q(t) dt — total energy delivered to skin
over the event. Lets us compare "how much heat does 1 ladle deliver?"
on the same axis as steady-state q.

Reference values:
- Latent heat per ladle (100 mL water):
    L * m = 2.26 MJ/kg * 0.1 kg = 226 kJ total
    on 1.8 m² body if fully condensed: 125 kJ/m² upper bound
- Steady "hot dry calm" sauna: q ~ 700 W/m²; 60 s of exposure = 42 kJ/m²
  → 1 ladle ≈ 3 minutes of normal exposure delivered in ~30 s

Usage:
    PYTHONPATH=src python scripts/skin_htc_timeseries.py
    PYTHONPATH=src python scripts/skin_htc_timeseries.py --output-dir results/skin_htc_timeseries

Outputs:
    timeseries.csv         per-event q(t), T_eq(t), Q_cum(t)
    events_q_t.png         q(t) for each named event, panelled
    events_t_eq.png        T_equivalent(t) for each event
    events_Q_cumulative.png cumulative Q(t) — total energy delivered
    energy_per_ladle.png   bar chart: Q absorbed per ladle vs ladle latent heat
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from harness.skin_htc import (
    BODY_SURFACE_AREA_M2,
    LATENT_HEAT_VAPORIZATION,
    compute_skin_balance,
    humidity_ratio_from_rh,
    humidity_ratio_from_water_addition,
    rh_from_humidity_ratio,
)


@dataclass
class SaunaEvent:
    """A time-evolving sauna scenario."""
    name: str
    duration_s: float
    base_t_air_c: float
    base_rh: float
    base_v: float
    base_t_mrt_c: float
    volume_m3: float
    pour_times_s: list[float] = field(default_factory=list)
    pour_water_kg: float = 0.1                # per pour (1 ladle = 0.1 kg)
    pour_decay_tau_s: float = 60.0            # humidity exp decay after pour
    aufguss_windows: list[tuple[float, float, float]] = field(default_factory=list)
    # each window: (t_start, t_end, V_during_aufguss [m/s])
    t_air_ramp: tuple[float, float] | None = None  # (t_start, t_end) for linear ramp
    t_air_target_c: float = 0.0


# ---------------------------------------------------------------------------
# Named events
# ---------------------------------------------------------------------------

EVENTS: list[SaunaEvent] = [
    # (1) Single löyly: hot dry baseline, one pour at t=15s, observe decay
    SaunaEvent(
        name="Single löyly (1 ladle)",
        duration_s=180.0,
        base_t_air_c=85.0, base_rh=0.07, base_v=0.10, base_t_mrt_c=90.0,
        volume_m3=10.0,
        pour_times_s=[15.0],
        pour_water_kg=0.1,
        pour_decay_tau_s=70.0,
    ),

    # (2) Aufguss session (3 pours + escalating fan/towel velocity)
    SaunaEvent(
        name="Aufguss session",
        duration_s=300.0,
        base_t_air_c=88.0, base_rh=0.08, base_v=0.10, base_t_mrt_c=92.0,
        volume_m3=15.0,
        pour_times_s=[20.0, 90.0, 170.0],
        pour_water_kg=0.1,
        pour_decay_tau_s=80.0,
        aufguss_windows=[
            ( 30.0,  60.0, 1.5),   # gentle wave after 1st pour
            (100.0, 140.0, 2.5),   # standard Aufguss after 2nd pour
            (180.0, 230.0, 3.5),   # intense finale after 3rd pour
        ],
    ),

    # (3) Shikiji herbal sauna — sustained high humidity, low T, calm.
    # The "intensity" comes from continuous condensation, not events.
    SaunaEvent(
        name="Shikiji herbal",
        duration_s=300.0,
        base_t_air_c=60.0, base_rh=0.95, base_v=0.10, base_t_mrt_c=62.0,
        volume_m3=6.0,
        pour_times_s=[],
    ),

    # (4) Kagaya-style progressive: temperature ramps up and pours are
    # spaced to build cumulative humidity → escalating intensity.
    SaunaEvent(
        name="Kagaya progressive",
        duration_s=420.0,
        base_t_air_c=80.0, base_rh=0.10, base_v=0.10, base_t_mrt_c=88.0,
        volume_m3=12.0,
        pour_times_s=[60.0, 150.0, 240.0, 330.0],
        pour_water_kg=0.1,
        pour_decay_tau_s=120.0,  # slower decay → cumulative build-up
        t_air_ramp=(0.0, 420.0),
        t_air_target_c=98.0,
    ),
]


def build_profile(event: SaunaEvent, dt: float = 1.0) -> dict:
    """Return arrays of (t, T_air, RH, V, T_mrt) for one event."""
    n = int(event.duration_s / dt) + 1
    t = np.linspace(0.0, event.duration_s, n)

    # T_air: optional linear ramp from base to target over a window
    t_air = np.full(n, event.base_t_air_c)
    if event.t_air_ramp is not None:
        ts, te = event.t_air_ramp
        ramp_frac = np.clip((t - ts) / max(te - ts, 1e-6), 0.0, 1.0)
        t_air = event.base_t_air_c + ramp_frac * (
            event.t_air_target_c - event.base_t_air_c
        )

    # T_mrt tracks T_air with the original offset
    mrt_offset = event.base_t_mrt_c - event.base_t_air_c
    t_mrt = t_air + mrt_offset

    # Humidity: start with baseline, add exp-decay pulses for each pour
    rh = np.full(n, event.base_rh)
    # Convert to humidity ratio for proper accounting, then back to RH
    w_arr = np.zeros(n)
    for i in range(n):
        w_arr[i] = humidity_ratio_from_rh(event.base_rh, t_air[i])
    for pt in event.pour_times_s:
        # peak humidity ratio added by this pour, decaying over time
        # use t_air at pour time for the dry-air mass calc
        idx_pour = min(int(pt / dt), n - 1)
        delta_w_peak = humidity_ratio_from_water_addition(
            event.pour_water_kg, event.volume_m3, t_air[idx_pour],
        ) - humidity_ratio_from_rh(0.0, t_air[idx_pour])
        # exp ramp-up over 3 s, then exp decay
        rise_tau = 3.0
        for i in range(n):
            if t[i] >= pt:
                age = t[i] - pt
                contribution = delta_w_peak * (1.0 - np.exp(-age / rise_tau)) \
                    * np.exp(-age / event.pour_decay_tau_s)
                w_arr[i] += contribution
    # Convert back to RH (clipped to 1.0)
    for i in range(n):
        rh[i] = rh_from_humidity_ratio(w_arr[i], t_air[i])

    # Velocity: baseline + Aufguss windows
    v = np.full(n, event.base_v)
    for ts, te, v_aug in event.aufguss_windows:
        # smooth ramp at edges (5 s rise/fall)
        ramp = 5.0
        for i in range(n):
            if ts <= t[i] <= te:
                v[i] = max(v[i], v_aug)
            elif ts - ramp < t[i] < ts:
                frac = (t[i] - (ts - ramp)) / ramp
                v[i] = max(v[i], event.base_v + frac * (v_aug - event.base_v))
            elif te < t[i] < te + ramp:
                frac = 1.0 - (t[i] - te) / ramp
                v[i] = max(v[i], event.base_v + frac * (v_aug - event.base_v))

    return {
        "t": t, "t_air_c": t_air, "rh": rh, "v": v, "t_mrt_c": t_mrt,
        "w": w_arr,
    }


def integrate_event(
    event: SaunaEvent,
    *,
    w_skin: float = 0.4,
    t_skin_c: float = 36.0,
    dt: float = 1.0,
) -> dict:
    """Run skin balance through the time profile, return all series + Q_cum."""
    p = build_profile(event, dt=dt)
    n = len(p["t"])

    q_arr = np.zeros(n)
    qc_arr = np.zeros(n)
    qr_arr = np.zeros(n)
    qe_arr = np.zeros(n)
    teq_arr = np.zeros(n)
    h_conv_arr = np.zeros(n)
    face_q_arr = np.zeros(n)

    for i in range(n):
        bal = compute_skin_balance(
            t_air_c=float(p["t_air_c"][i]),
            rh=float(p["rh"][i]),
            v_local=float(p["v"][i]),
            t_mrt_c=float(p["t_mrt_c"][i]),
            w_skin=w_skin,
            t_skin_c=t_skin_c,
        )
        wb = bal.whole_body
        q_arr[i] = wb.q_total
        qc_arr[i] = wb.q_conv
        qr_arr[i] = wb.q_rad
        qe_arr[i] = wb.q_evap
        teq_arr[i] = wb.t_equivalent_c
        h_conv_arr[i] = wb.h_conv
        face_q_arr[i] = bal.parts["face"].q_total

    # Trapezoid cumulative integration: Q_cum [J/m²]
    Q_cum = np.zeros(n)
    for i in range(1, n):
        dt_i = p["t"][i] - p["t"][i - 1]
        Q_cum[i] = Q_cum[i - 1] + 0.5 * (q_arr[i] + q_arr[i - 1]) * dt_i

    return {
        **p,
        "q": q_arr, "q_conv": qc_arr, "q_rad": qr_arr, "q_evap": qe_arr,
        "t_eq": teq_arr, "h_conv": h_conv_arr, "face_q": face_q_arr,
        "Q_cum": Q_cum,
    }


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def plot_events_panel(
    series: list[tuple[SaunaEvent, dict]],
    key: str,
    ylabel: str,
    title: str,
    out_path: Path,
    extra_lines: dict | None = None,
) -> None:
    """One panel per event, plotting `series[i][1][key]` vs time."""
    n = len(series)
    cols = 2
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(13, 3.0 * rows), sharex=False)
    axes = np.array(axes).flatten()
    for idx, (event, s) in enumerate(series):
        ax = axes[idx]
        ax.plot(s["t"], s[key], color="tab:red", linewidth=1.7, label=key)
        # Mark pours
        for pt in event.pour_times_s:
            ax.axvline(pt, color="blue", linestyle=":", alpha=0.7)
        # Mark Aufguss windows
        for ts, te, _ in event.aufguss_windows:
            ax.axvspan(ts, te, color="orange", alpha=0.15)
        if extra_lines:
            for label, value in extra_lines.items():
                ax.axhline(value, color="gray", linestyle="--", alpha=0.5, label=label)
        ax.set_title(f"{event.name}  (T_base={event.base_t_air_c:.0f}°C)", fontsize=10)
        ax.set_xlabel("time [s]")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=8, loc="upper right")
    for k in range(n, len(axes)):
        axes[k].axis("off")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_q_components_for_aufguss(s: dict, event: SaunaEvent, out_path: Path) -> None:
    """Component breakdown of q(t) during the Aufguss session."""
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(s["t"], s["q_conv"], label="q_conv", linewidth=1.4)
    ax.plot(s["t"], s["q_rad"],  label="q_rad",  linewidth=1.4, linestyle="--")
    ax.plot(s["t"], s["q_evap"], label="q_evap (cond>0 / sweat<0)",
            linewidth=1.4, linestyle=":")
    ax.plot(s["t"], s["q"],      label="q_total", linewidth=2.3, color="black")
    for pt in event.pour_times_s:
        ax.axvline(pt, color="blue", linestyle=":", alpha=0.7,
                   label="pour" if pt == event.pour_times_s[0] else None)
    for ts, te, v in event.aufguss_windows:
        is_first = (ts, te) == event.aufguss_windows[0][:2]
        lbl = f"Aufguss V={v:.1f}" if is_first else None
        ax.axvspan(ts, te, color="orange", alpha=0.15, label=lbl)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("Heat flux [W/m²]")
    ax.set_title(f"{event.name} — component breakdown")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_energy_per_ladle(series: list[tuple[SaunaEvent, dict]], out_path: Path) -> None:
    """Bar chart: total Q delivered vs latent heat budget per ladle."""
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 5))
    names = [e.name for e, _ in series]

    # Left: total Q absorbed [kJ/m²] and equivalent total energy on body [kJ]
    Q_per_m2 = np.array([s["Q_cum"][-1] / 1000.0 for _, s in series])  # kJ/m²
    Q_total_body = Q_per_m2 * BODY_SURFACE_AREA_M2  # kJ

    x = np.arange(len(series))
    width = 0.35
    ax_left.bar(x - width / 2, Q_per_m2, width, color="tab:red",
                label="kJ/m² absorbed")
    ax_left.bar(x + width / 2, Q_total_body, width, color="tab:orange",
                label=f"kJ total (× {BODY_SURFACE_AREA_M2} m² body)")
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(names, fontsize=8, rotation=20, ha="right")
    ax_left.set_ylabel("Energy delivered to skin")
    ax_left.set_title("Cumulative Q over the whole event")
    ax_left.legend(fontsize=8)
    ax_left.grid(True, axis="y", alpha=0.3)

    # Right: per-ladle breakdown — Q absorbed per pour vs latent heat per pour
    latent_per_ladle_kj = LATENT_HEAT_VAPORIZATION * 0.1 / 1000.0  # 226 kJ
    latent_per_ladle_per_m2 = (LATENT_HEAT_VAPORIZATION * 0.1) / BODY_SURFACE_AREA_M2 / 1000.0

    # Per-pour comparison: only meaningful for events that include pours
    pour_mask = np.array([len(e.pour_times_s) > 0 for e, _ in series])
    n_pours = np.array([max(len(e.pour_times_s), 1) for e, _ in series])
    Q_per_pour_per_m2 = np.where(pour_mask, Q_per_m2 / n_pours, np.nan)

    bar_x = x[pour_mask]
    bar_vals = Q_per_pour_per_m2[pour_mask]
    bar_labels = [n for n, m in zip(names, pour_mask, strict=True) if m]
    ax_right.bar(bar_x, bar_vals, color="tab:red", alpha=0.85,
                 label="actual Q absorbed per pour [kJ/m²]")
    latent_label = (
        f"latent heat in 1 ladle / 1.8 m² ({latent_per_ladle_per_m2:.0f} kJ/m²)"
    )
    ax_right.axhline(latent_per_ladle_per_m2, color="black", linestyle="--",
                     linewidth=1.4, label=latent_label)
    # Annotate skipped (no-pour) events
    skipped = [n for n, m in zip(names, pour_mask, strict=True) if not m]
    if skipped:
        ax_right.text(0.98, 0.02, "no pours: " + ", ".join(skipped),
                      transform=ax_right.transAxes, fontsize=8, ha="right",
                      style="italic", color="gray")
    ax_right.set_xticks(bar_x)
    ax_right.set_xticklabels(bar_labels, fontsize=8, rotation=20, ha="right")
    ax_right.set_ylabel("kJ/m² per pour")
    ax_right.set_title("How much of the löyly latent heat actually reaches skin?")
    ax_right.legend(fontsize=8, loc="upper right")
    ax_right.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        f"Q (time-integrated heat absorbed by skin) — 1 ladle latent heat = "
        f"{latent_per_ladle_kj:.0f} kJ"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def write_csv(series: list[tuple[SaunaEvent, dict]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for event, s in series:
        for i in range(len(s["t"])):
            rows.append({
                "event": event.name,
                "t_s": float(s["t"][i]),
                "T_air_C": float(s["t_air_c"][i]),
                "T_mrt_C": float(s["t_mrt_c"][i]),
                "RH": float(s["rh"][i]),
                "V_m_s": float(s["v"][i]),
                "humidity_ratio_g_kg": float(s["w"][i]) * 1000.0,
                "q_conv": float(s["q_conv"][i]),
                "q_rad": float(s["q_rad"][i]),
                "q_evap": float(s["q_evap"][i]),
                "q_total": float(s["q"][i]),
                "T_equivalent_C": float(s["t_eq"][i]),
                "Q_cumulative_J_per_m2": float(s["Q_cum"][i]),
                "face_q_total": float(s["face_q"][i]),
            })
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {out_path} ({len(rows)} rows)")


def print_summary(series: list[tuple[SaunaEvent, dict]]) -> None:
    print(f"\n{'Event':<26}{'dur':>6}{'q_max':>9}{'T_eq_pk':>9}"
          f"{'Q_total':>11}{'Q/ladle':>11}{'Q/body':>11}")
    print("  " + "-" * 86)
    for event, s in series:
        Q_total_kj_m2 = s["Q_cum"][-1] / 1000.0
        Q_total_body_kj = Q_total_kj_m2 * BODY_SURFACE_AREA_M2
        n_pours = max(len(event.pour_times_s), 1)
        Q_per_pour = Q_total_kj_m2 / n_pours
        print(f"  {event.name:<24}{event.duration_s:>6.0f}"
              f"{s['q'].max():>9.0f}{s['t_eq'].max():>9.1f}"
              f"{Q_total_kj_m2:>9.0f} kJ/m²"
              f"{Q_per_pour:>9.1f} kJ/m²"
              f"{Q_total_body_kj:>9.0f} kJ")
    print("\n  Reference: latent heat in 1 ladle (100 mL) = 226 kJ total,")
    print("             ÷ 1.8 m² body = 125 kJ/m² upper bound.\n")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/skin_htc_timeseries"))
    p.add_argument("--w-skin", type=float, default=0.4)
    p.add_argument("--t-skin", type=float, default=36.0)
    p.add_argument("--dt", type=float, default=1.0, help="Time step [s]")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing outputs to {out_dir}/")

    series: list[tuple[SaunaEvent, dict]] = []
    for event in EVENTS:
        s = integrate_event(event, w_skin=args.w_skin, t_skin_c=args.t_skin,
                            dt=args.dt)
        series.append((event, s))

    write_csv(series, out_dir / "timeseries.csv")
    plot_events_panel(
        series, "q", "q [W/m²]",
        "Whole-body skin heat flux q(t)  (blue|löyly  orange|Aufguss)",
        out_dir / "events_q_t.png",
    )
    plot_events_panel(
        series, "t_eq", "T_equivalent [°C]",
        "Perceived (equivalent) skin temperature T_eq(t)",
        out_dir / "events_t_eq.png",
    )
    plot_events_panel(
        series, "Q_cum", "Q_cum [J/m²]",
        "Cumulative absorbed energy Q(t) = ∫ q dt",
        out_dir / "events_Q_cumulative.png",
    )
    # Find the Aufguss event for component breakdown
    for event, s in series:
        if "Aufguss" in event.name:
            plot_q_components_for_aufguss(
                s, event, out_dir / "aufguss_components.png",
            )
            break
    plot_energy_per_ladle(series, out_dir / "energy_per_ladle.png")

    print_summary(series)


if __name__ == "__main__":
    main()
