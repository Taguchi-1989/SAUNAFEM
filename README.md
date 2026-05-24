# SaunaFlow

**CFD-driven thermal simulation of sauna environments — from a zero-dimensional browser calculator to full 3D OpenFOAM analysis.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.11+-3776ab.svg)](https://python.org)
[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-v2312-c00?logo=openfoam)](https://openfoam.org)
[![GitHub Pages](https://img.shields.io/badge/Live_Demo-GitHub_Pages-222?logo=github)](https://taguchi-1989.github.io/SAUNAFEM/)

---

> *Why does the same sauna feel completely different depending on who's running it?*
> The answer is in the fluid dynamics.

---

## Live Demo — No install required

**[→ Löyly Calculator (interactive)](https://taguchi-1989.github.io/SAUNAFEM/)**

A zero-dimensional heat balance simulator running entirely in your browser:
- Adjust heater power, room humidity, water volume, and ventilation
- Instant Humidex body-feel index, latent/convective/radiant heat split
- Dependency graphs and sensitivity popovers with governing equations

---

## What is SaunaFlow?

A sauna is a deceptively complex thermal system:

- An 18 kW heater with non-uniform surface temperatures
- Buoyancy-driven stratification (upper bench 95 °C, lower bench 54 °C in the same room)
- Transient steam from löyly — a pulse of latent heat and vapor that travels upward in seconds
- Aerodynamic disruption from Aufguss towel waves

None of this is captured by rule-of-thumb temperature guidelines.

**SaunaFlow** is a Python harness that turns declarative YAML case definitions into fully reproducible OpenFOAM simulations, computes sauna-specific KPIs (K-01 through K-07), and is being built toward experimental sensor validation.

---

## CFD Results

13 parametric cases run with `buoyantPimpleFoam` on a Finnish-style sauna geometry:

| Case | Heater | Upper bench | Lower bench | Target range |
|------|--------|-------------|-------------|--------------|
| L-1 | 13 kW | **95 °C** | **54 °C** | 80–100 / 40–60 °C ✅ |
| K-1 | 18 kW | 102 °C | 61 °C | — |
| K-2 | 18 kW + vent | 88 °C | 49 °C | — |
| K-3 | 18 kW + löyly | peak +12 °C | peak +7 °C | — |

**Key finding:** `buoyantPimpleFoam` is essential. `simpleFoam`-based solvers do not converge for buoyancy-dominated flows even at 50,000 iterations — the pressure-velocity coupling cannot resolve the density-driven feedback loop.

---

## Physics modeled

| Phenomenon | Model |
|-----------|-------|
| Buoyancy-driven stratification | Boussinesq approximation, k-ε turbulence |
| Heater radiation | View factor (surface-to-surface) |
| Löyly steam | Transient vapor volume source with latent heat |
| Aufguss airflow | Momentum source (towel wave), local convective enhancement |
| Skin heat transfer | 0D: convection + radiation + evaporative (sweating) |

---

## Architecture

```
configs/cases/dry_sauna_steady.yaml   ← declarative input
          │
          ▼
   case_builder.py    ← expands YAML → OpenFOAM directory
          │
          ▼
   solver_runner.py   ← runs buoyantPimpleFoam (local or WSL2)
          │
          ▼
   probe_parser.py    ← reads time-series probe output
          │
          ▼
        kpi.py        ← K-01 … K-07
          │
          ▼
    reporting.py      ← Markdown + HTML report
```

Swap a YAML file, re-run — identical inputs produce identical results.

---

## KPIs

| ID | Metric |
|----|--------|
| K-01 | Steady-state temperature differential (upper / lower bench) |
| K-02 | Post-löyly peak temperature |
| K-03 | Post-löyly peak humidity |
| K-04 | Steam peak arrival time |
| K-05 | Face-level wind speed peak (Aufguss) |
| K-06 | Simplified thermal stress index |
| K-07 | Upper / lower relative temperature difference |

---

## Quickstart

**Prerequisites:** Python 3.11+, OpenFOAM v2312 (Ubuntu native or WSL2)

```bash
git clone https://github.com/Taguchi-1989/SAUNAFEM.git
cd SAUNAFEM
pip install -e .
```

```bash
# Build OpenFOAM case from YAML
PYTHONPATH=src python -m harness.cli build configs/cases/dry_sauna_steady.yaml

# Run solver (WSL2)
wsl -- /usr/bin/openfoam2312 bash scripts/run_openfoam_wsl.sh

# Generate report
PYTHONPATH=src python -m harness.cli report results/dry_sauna_steady/

# Tests
pytest tests/ -x -q
```

---

## Parametric study coverage

29 YAML-defined cases across:

- Heater power: 8–18 kW
- Löyly water: 0.1–0.5 L per throw
- Ventilation: none / natural / forced
- Wall material: spruce / cedar, 45–90 mm
- Occupancy: empty / 2-person

---

## Repository layout

```
src/harness/       Python orchestration — CLI, case builder, solver runner, KPI, reporting
configs/cases/     29 YAML case definitions
foam_templates/    OpenFOAM base case templates (generated, never hand-edited)
tools/             Standalone HTML calculators (no server required)
experiments/       Sensor data — raw CSV + processed (Phase 4)
firmware/          Sensor firmware for the experimental rig
docs/              Governing equations (PDF/TeX), parametric study reports, field notes
```

---

## Project status

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Repo, schema, CLI | ✅ |
| 1 | Dry sauna steady-state CFD | ✅ |
| 2 | Löyly transient | ✅ |
| 3 | Aufguss momentum source | 🔄 |
| 4 | Experimental sensor validation | 📋 |
| 5 | Automated batch comparison | 📋 |

---

## Contributing

Issues and PRs are welcome. This is likely the only public repository combining OpenFOAM with sauna-specific thermal modeling, so if you're working on:

- Thermal comfort in small enclosures
- Transient buoyancy-driven flows
- Building physics / HVAC with OpenFOAM
- Sauna design or research

...you're in the right place.

---

---

## 日本語ガイド

**SaunaFlow** は、サウナ熱環境を物理シミュレーションで解析するオープンソースプロジェクトです。

### なぜ作ったか

「同じサウナなのに、誰が管理するかで体感が全然違う」という疑問が出発点です。温度計の数値だけでは説明できない — ロウリュの蒸気がどう広がるか、上下の温度成層がどう形成されるか、Aufgussの風がどう熱流束を変えるか。これを流体力学で定量化することが目的です。

### できること

| | ツール |
|---|---|
| **ブラウザだけで今すぐ試す** | [ロウリュ体感シミュレータ](https://taguchi-1989.github.io/SAUNAFEM/) — スライダーを動かすと即座にHumidex・熱収支・感度グラフが更新 |
| **OpenFOAMでCFD計算** | YAMLで条件定義 → 自動でOpenFOAMケース生成 → ソルバー実行 → KPI自動計算 |
| **パラメトリックスタディ** | ヒーター出力・ロウリュ水量・換気・壁材の組み合わせを29ケース実施済み |

### 主な発見

- **上段95°C / 下段54°C** の温度成層を13kWヒーターで再現（実測目標値内）
- `buoyantPimpleFoam` が浮力主導流れに必須。`simpleFoam` は5万イテレーションでも収束しない
- ヒーター壁温の感度は低い（473→673K で気温差 <1°C）
- 換気の境界条件設定が結果に大きく影響する

### 参加・フィードバック

Issueで質問・提案を歓迎します。サウナ設計、建築熱環境、OpenFOAM活用に興味のある方はぜひ。

---

*MIT License*
