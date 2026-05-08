"""Sauna scenario explorer: compare types, löyly amounts, and Aufguss effects.

Builds named real-world sauna scenarios from measurable inputs
(T_air, room volume, water poured, fan/towel velocity) using the
standalone skin heat-balance calculator. Independent of the CFD/zone
solver.

Usage:
    PYTHONPATH=src python scripts/skin_htc_scenarios.py
    PYTHONPATH=src python scripts/skin_htc_scenarios.py --output-dir results/skin_htc_scenarios

Outputs (under --output-dir):
    scenarios.csv                      flat table of every scenario × treatment
    sauna_types_baseline.png           5 sauna types side-by-side (baseline)
    loyly_progression.png              q_total vs # ladles, with/without Aufguss
    treatment_matrix.png               heatmap [scenario × treatment] of q_total
    face_focus.png                     face vs whole-body across treatments
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from harness.skin_htc import (
    compute_skin_balance,
    humidity_ratio_from_rh,
    humidity_ratio_from_water_addition,
    rh_from_humidity_ratio,
)

# 1 standard ladle ≈ 100 mL ≈ 0.1 kg water
LADLE_KG = 0.1


@dataclass(frozen=True)
class SaunaType:
    """A real-world sauna profile defined by measurable parameters."""
    name: str
    t_air_c: float           # measured dry-bulb temperature [°C]
    baseline_rh: float       # ambient humidity before löyly [-]
    t_mrt_c: float           # mean radiant temperature [°C]; ≈ wall temp
    v_calm: float            # background air velocity, no Aufguss [m/s]
    v_aufguss: float         # face-level velocity during Aufguss [m/s]
    volume_m3: float         # room volume for löyly humidity calc [m³]


# Canonical sauna types (rough but realistic).  Two Japanese well-known
# venues are included: Shikiji's herbal sauna (低温高湿度・薬草スチーム) and
# a generic "Kagaya-style" progressive löyly room used for the timeseries
# study — its baseline is conservative, the magic happens over time.
SAUNA_TYPES: list[SaunaType] = [
    SaunaType("Finnish dry",     90.0, 0.05,  95.0, 0.10, 2.5, 10.0),
    SaunaType("Russian banya",   65.0, 0.30,  68.0, 0.10, 1.5,  8.0),
    SaunaType("German Aufguss",  85.0, 0.08,  90.0, 0.10, 3.0, 15.0),
    SaunaType("Japanese rock",   65.0, 0.15,  75.0, 0.10, 1.0, 10.0),
    SaunaType("Turkish hammam",  45.0, 0.95,  46.0, 0.10, 0.5, 12.0),
    # Shikiji 薬草サウナ (Shizuoka): ~60°C steam-saturated herbal mist, calm
    SaunaType("Shikiji herbal",  60.0, 0.95,  62.0, 0.10, 0.3,  6.0),
    # Kagaya-style progressive: hot dry baseline, intense löyly events
    SaunaType("Kagaya prog.",    88.0, 0.10,  98.0, 0.10, 2.0, 12.0),
]


@dataclass(frozen=True)
class Treatment:
    label: str
    ladles: int
    aufguss: bool


TREATMENTS: list[Treatment] = [
    Treatment("baseline",         0, False),
    Treatment("+1 ladle",         1, False),
    Treatment("+2 ladles",        2, False),
    Treatment("Aufguss only",     0, True),
    Treatment("+1 ladle\n+Aufguss", 1, True),
    Treatment("+2 ladles\n+Aufguss", 2, True),
]


def evaluate(
    sauna: SaunaType,
    treatment: Treatment,
    *,
    w_skin: float = 0.4,
    t_skin_c: float = 36.0,
):
    """Run one (sauna, treatment) and return the SkinHeatBalance plus inputs."""
    # Baseline humidity from measured RH
    w0 = humidity_ratio_from_rh(sauna.baseline_rh, sauna.t_air_c)
    # Add löyly water (well-mixed)
    water_kg = treatment.ladles * LADLE_KG
    w_total = humidity_ratio_from_water_addition(
        water_kg, sauna.volume_m3, sauna.t_air_c, initial_humidity_ratio=w0,
    )
    rh_total = rh_from_humidity_ratio(w_total, sauna.t_air_c)
    v = sauna.v_aufguss if treatment.aufguss else sauna.v_calm
    bal = compute_skin_balance(
        t_air_c=sauna.t_air_c,
        rh=rh_total,
        v_local=v,
        t_mrt_c=sauna.t_mrt_c,
        t_skin_c=t_skin_c,
        w_skin=w_skin,
    )
    return bal, dict(rh_total=rh_total, v=v, water_kg=water_kg, w_humidity=w_total)


def write_csv(out_path: Path) -> None:
    """Run the full scenario × treatment matrix and dump as CSV."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for sauna in SAUNA_TYPES:
        for tr in TREATMENTS:
            bal, ctx = evaluate(sauna, tr)
            wb = bal.whole_body
            face = bal.parts["face"]
            rows.append({
                "sauna": sauna.name,
                "treatment": tr.label.replace("\n", " "),
                "ladles": tr.ladles,
                "aufguss": tr.aufguss,
                "T_air_C": sauna.t_air_c,
                "T_mrt_C": sauna.t_mrt_c,
                "V_used_m_s": ctx["v"],
                "RH_after": ctx["rh_total"],
                "humidity_ratio_g_kg": ctx["w_humidity"] * 1000.0,
                "wb_h_conv": wb.h_conv,
                "wb_h_rad": wb.h_rad,
                "wb_q_conv": wb.q_conv,
                "wb_q_rad": wb.q_rad,
                "wb_q_evap": wb.q_evap,
                "wb_q_total": wb.q_total,
                "wb_u_overall": wb.u_overall,
                "wb_t_equivalent_C": wb.t_equivalent_c,
                "face_q_total": face.q_total,
                "face_h_conv": face.h_conv,
            })
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote {out_path} ({len(rows)} rows)")


def plot_sauna_types_baseline(out_path: Path) -> None:
    """Stacked-bar comparison of components across sauna types at baseline."""
    qc, qr, qe, qt = [], [], [], []
    rh_after, v_used = [], []
    for sauna in SAUNA_TYPES:
        bal, ctx = evaluate(sauna, TREATMENTS[0])
        wb = bal.whole_body
        qc.append(wb.q_conv)
        qr.append(wb.q_rad)
        qe.append(wb.q_evap)
        qt.append(wb.q_total)
        rh_after.append(ctx["rh_total"])
        v_used.append(ctx["v"])

    x = np.arange(len(SAUNA_TYPES))
    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.7
    # Stack positives
    pos_bottom = np.zeros(len(x))
    neg_bottom = np.zeros(len(x))
    for vals, color, label in [
        (qc, "tab:red", "q_conv"),
        (qr, "tab:orange", "q_rad"),
        (qe, "tab:blue", "q_evap"),
    ]:
        arr = np.array(vals)
        pos = np.where(arr > 0, arr, 0)
        neg = np.where(arr < 0, arr, 0)
        ax.bar(x, pos, width, bottom=pos_bottom, color=color, alpha=0.85,
               edgecolor="black", linewidth=0.4, label=label)
        ax.bar(x, neg, width, bottom=neg_bottom, color=color, alpha=0.85,
               edgecolor="black", linewidth=0.4)
        pos_bottom += pos
        neg_bottom += neg

    ax.scatter(x, qt, color="black", s=50, marker="D", zorder=5,
               label="q_total (sum)")
    ax.axhline(0, color="black", linewidth=0.5)
    labels = [
        f"{s.name}\nT={s.t_air_c:.0f}°C\nRH={rh:.1%}, V={v:.1f}"
        for s, rh, v in zip(SAUNA_TYPES, rh_after, v_used, strict=True)
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Heat flux [W/m²]")
    ax.set_title("Sauna type comparison at baseline (no löyly, no Aufguss)")
    ax.legend(loc="upper left", ncol=4, fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_loyly_progression(out_path: Path) -> None:
    """For each sauna type, q_total vs ladle count, with and without Aufguss."""
    n_ladles = list(range(0, 7))  # 0, 1, ..., 6 ladles

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=False)
    axes = axes.flatten()
    for idx, sauna in enumerate(SAUNA_TYPES):
        ax = axes[idx]
        for aufguss in (False, True):
            v = sauna.v_aufguss if aufguss else sauna.v_calm
            qts, rhs = [], []
            for n in n_ladles:
                tr = Treatment(f"+{n} ladles", n, aufguss)
                bal, ctx = evaluate(sauna, tr)
                qts.append(bal.whole_body.q_total)
                rhs.append(ctx["rh_total"])
            label = f"{'Aufguss' if aufguss else 'calm'} (V={v:.1f} m/s)"
            ax.plot(n_ladles, qts, "-o", label=label, linewidth=1.8)
        ax.set_title(f"{sauna.name}\nT={sauna.t_air_c:.0f}°C, V={sauna.volume_m3:.0f} m³")
        ax.set_xlabel("# ladles (100 mL each)")
        ax.set_ylabel("q_total [W/m²]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    # hide unused subplot
    for k in range(len(SAUNA_TYPES), len(axes)):
        axes[k].axis("off")
    fig.suptitle("Löyly dose-response — whole-body heat flux vs # ladles", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_treatment_matrix(out_path: Path) -> None:
    """Heatmap of q_total over (sauna × treatment) and a second of T_equivalent."""
    qt_grid = np.zeros((len(SAUNA_TYPES), len(TREATMENTS)))
    teq_grid = np.zeros_like(qt_grid)
    for i, sauna in enumerate(SAUNA_TYPES):
        for j, tr in enumerate(TREATMENTS):
            bal, _ = evaluate(sauna, tr)
            qt_grid[i, j] = bal.whole_body.q_total
            teq_grid[i, j] = bal.whole_body.t_equivalent_c

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, grid, title, cmap, unit in [
        (axes[0], qt_grid, "q_total [W/m²]", "hot", "W/m²"),
        (axes[1], teq_grid, "T_equivalent [°C]", "magma", "°C"),
    ]:
        im = ax.imshow(grid, aspect="auto", cmap=cmap)
        ax.set_xticks(np.arange(len(TREATMENTS)))
        ax.set_xticklabels([t.label for t in TREATMENTS], fontsize=9)
        ax.set_yticks(np.arange(len(SAUNA_TYPES)))
        ax.set_yticklabels([s.name for s in SAUNA_TYPES], fontsize=9)
        # Annotate cells
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                ax.text(j, i, f"{grid[i, j]:.0f}",
                        ha="center", va="center", fontsize=8,
                        color="white" if grid[i, j] > grid.mean() else "black")
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(unit)
        ax.set_title(title)
    fig.suptitle("Treatment matrix: each cell is one (sauna type × treatment) result",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_face_focus(out_path: Path) -> None:
    """Face vs whole-body q_total under each treatment, for each sauna."""
    n_sauna = len(SAUNA_TYPES)
    n_tr = len(TREATMENTS)
    width = 0.35
    x = np.arange(n_tr)

    fig, axes = plt.subplots(n_sauna, 1, figsize=(11, 2.5 * n_sauna), sharex=True)
    for idx, sauna in enumerate(SAUNA_TYPES):
        ax = axes[idx]
        face_q, wb_q = [], []
        for tr in TREATMENTS:
            bal, _ = evaluate(sauna, tr)
            face_q.append(bal.parts["face"].q_total)
            wb_q.append(bal.whole_body.q_total)
        ax.bar(x - width / 2, wb_q, width, label="whole body", color="tab:blue")
        ax.bar(x + width / 2, face_q, width, label="face", color="tab:red")
        ax.set_title(f"{sauna.name}  (T={sauna.t_air_c:.0f}°C, V_aug={sauna.v_aufguss:.1f})",
                     fontsize=10)
        ax.set_ylabel("q_total [W/m²]")
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([t.label for t in TREATMENTS], fontsize=9)
    fig.suptitle("Face vs whole-body — face takes the brunt of Aufguss + löyly",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"  wrote {out_path}")


def print_summary(t_skin_c: float, w_skin: float) -> None:
    """Compact text table of all (sauna × treatment) cells."""
    print(f"\n(t_skin = {t_skin_c:.1f}°C, w_skin = {w_skin:.2f})")
    print(f"\n{'Sauna':<18}{'Treatment':<22}{'T':>5}{'V':>6}{'RH':>7}"
          f"{'q_c':>9}{'q_r':>9}{'q_e':>9}{'q_tot':>9}{'T_eq':>8}")
    print("  " + "-" * 100)
    for sauna in SAUNA_TYPES:
        for tr in TREATMENTS:
            bal, ctx = evaluate(sauna, tr, w_skin=w_skin, t_skin_c=t_skin_c)
            wb = bal.whole_body
            label = tr.label.replace("\n", " ")
            print(f"{sauna.name:<18}{label:<22}{sauna.t_air_c:>5.0f}{ctx['v']:>6.2f}"
                  f"{ctx['rh_total']:>7.2%}"
                  f"{wb.q_conv:>9.0f}{wb.q_rad:>9.0f}{wb.q_evap:>9.0f}"
                  f"{wb.q_total:>9.0f}{wb.t_equivalent_c:>8.1f}")
        print()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output-dir", type=Path,
                   default=Path("results/skin_htc_scenarios"))
    p.add_argument("--w-skin", type=float, default=0.4,
                   help="Skin wettedness (sweating intensity), 0=dry, 1=fully wet")
    p.add_argument("--t-skin", type=float, default=36.0,
                   help="Mean skin temperature [°C]")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Writing outputs to {out_dir}/")
    write_csv(out_dir / "scenarios.csv")
    plot_sauna_types_baseline(out_dir / "sauna_types_baseline.png")
    plot_loyly_progression(out_dir / "loyly_progression.png")
    plot_treatment_matrix(out_dir / "treatment_matrix.png")
    plot_face_focus(out_dir / "face_focus.png")
    print_summary(args.t_skin, args.w_skin)


if __name__ == "__main__":
    main()
