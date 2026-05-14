"""Generate a self-contained HTML report from an OpenFOAM case result directory.

Reads:
    <case>/postProcessing/probes/0/T
    <case>/postProcessing/volAverageT/0/volFieldValue.dat (or similar)
    <case>/postProcessing/wallHeatFlux/0/wallHeatFlux.dat

Writes:
    <case>/report.html (standalone, opens in any browser; uses Chart.js via CDN)

Usage:
    python scripts/generate_html_report.py results/L1
    python scripts/generate_html_report.py results/parametric_J1 --title "J1 case"
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path


def parse_probe_T(probe_dir: Path) -> tuple[list[float], list[list[float]], list[str]]:
    """Return (times, [series per probe], probe_labels)."""
    if not probe_dir.exists():
        return [], [], []
    sub = next(probe_dir.iterdir(), None)
    if sub is None:
        return [], [], []
    f = sub / "T"
    if not f.exists():
        return [], [], []
    text = f.read_text(encoding="utf-8", errors="ignore")
    labels: list[str] = []
    times: list[float] = []
    cols: list[list[float]] = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("# Probe"):
            m = re.match(r"# Probe \d+ \(([^)]+)\)", s)
            labels.append(m.group(1) if m else s)
            continue
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if not parts:
            continue
        try:
            t = float(parts[0])
        except ValueError:
            continue
        vals = [float(x) for x in parts[1:]]
        if not cols:
            cols = [[] for _ in vals]
        times.append(t)
        for i, v in enumerate(vals):
            if i < len(cols):
                cols[i].append(v - 273.15)  # K -> °C
    if not labels:
        labels = [f"probe{i}" for i in range(len(cols))]
    return times, cols, labels[: len(cols)]


def parse_vol_average(vol_dir: Path) -> tuple[list[float], list[float]]:
    if not vol_dir.exists():
        return [], []
    sub = next(vol_dir.iterdir(), None)
    if sub is None:
        return [], []
    files = list(sub.glob("*"))
    if not files:
        return [], []
    text = files[0].read_text(encoding="utf-8", errors="ignore")
    times, vals = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            times.append(float(parts[0]))
            vals.append(float(parts[1]) - 273.15)
        except ValueError:
            continue
    return times, vals


def parse_wall_heat_flux(wf_dir: Path) -> dict:
    """Return {patches: [...], times: [...], min: {p:[...]}, max:{}, integral:{}}."""
    if not wf_dir.exists():
        return {}
    sub = next(wf_dir.iterdir(), None)
    if sub is None:
        return {}
    files = list(sub.glob("*"))
    if not files:
        return {}
    text = files[0].read_text(encoding="utf-8", errors="ignore")
    series_min: dict[str, list[float]] = defaultdict(list)
    series_max: dict[str, list[float]] = defaultdict(list)
    series_int: dict[str, list[float]] = defaultdict(list)
    times_per_patch: dict[str, list[float]] = defaultdict(list)
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 5:
            continue
        try:
            t = float(parts[0])
            patch = parts[1]
            mn = float(parts[2])
            mx = float(parts[3])
            integral = float(parts[4])
        except ValueError:
            continue
        times_per_patch[patch].append(t)
        series_min[patch].append(mn)
        series_max[patch].append(mx)
        series_int[patch].append(integral)
    patches = list(times_per_patch.keys())
    if not patches:
        return {}
    times = times_per_patch[patches[0]]
    return {
        "patches": patches,
        "times": times,
        "min": dict(series_min),
        "max": dict(series_max),
        "integral": dict(series_int),
    }


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<title>SaunaFlow CFD Report — {title}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
  body {{ font-family: -apple-system, "Segoe UI", "Hiragino Sans", sans-serif; margin: 0; padding: 20px; background: #1a1a1a; color: #e8e8e8; }}
  h1 {{ font-size: 20px; margin: 0 0 6px; color: #ffb96b; }}
  h2 {{ font-size: 14px; margin: 24px 0 8px; color: #ffb96b; border-bottom: 1px solid #333; padding-bottom: 4px; }}
  .meta {{ color: #888; font-size: 12px; margin-bottom: 18px; }}
  .kpis {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(170px, 1fr)); gap: 10px; margin-bottom: 14px; }}
  .kpi {{ background: #242424; border: 1px solid #333; border-radius: 6px; padding: 10px 12px; }}
  .kpi .lbl {{ font-size: 10px; color: #999; text-transform: uppercase; letter-spacing: 0.5px; }}
  .kpi .num {{ font-size: 20px; font-weight: 700; color: #ffb96b; font-variant-numeric: tabular-nums; }}
  .kpi .unit {{ font-size: 11px; color: #aaa; margin-left: 3px; }}
  .panel {{ background: #242424; border: 1px solid #333; border-radius: 6px; padding: 14px; margin-bottom: 14px; }}
  canvas {{ background: #1e1e1e; border-radius: 4px; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
  th, td {{ padding: 6px 10px; text-align: right; border-bottom: 1px solid #333; font-variant-numeric: tabular-nums; }}
  th {{ color: #999; font-weight: 600; text-align: left; }}
  td:first-child, th:first-child {{ text-align: left; color: #ddd; }}
  .note {{ font-size: 10px; color: #888; margin-top: 8px; }}
</style>
</head>
<body>
<h1>SaunaFlow CFD Report — {title}</h1>
<div class="meta">case dir: <code>{case_dir}</code> · times: {t_start:.1f} – {t_end:.1f} s · probes: {n_probes}</div>

<div class="kpis">
{kpi_html}
</div>

<h2>プローブ温度の時間推移</h2>
<div class="panel"><div style="height:320px"><canvas id="chartProbes"></canvas></div></div>

<h2>体積平均温度の収束</h2>
<div class="panel"><div style="height:240px"><canvas id="chartVol"></canvas></div></div>

<h2>壁面熱流束（パッチ別 integral）</h2>
<div class="panel"><div style="height:300px"><canvas id="chartWall"></canvas></div></div>

<h2>最終時刻パッチサマリ</h2>
<div class="panel">
<table id="patchTable"></table>
<div class="note">符号: 正=外部から流入, 負=外部へ放熱（OpenFOAM wallHeatFlux 規約）</div>
</div>

<script>
const DATA = {data_json};

function makeChart(id, datasets, ylabel, xlabel='時間 [s]') {{
  return new Chart(document.getElementById(id), {{
    type: 'line',
    data: {{ datasets }},
    options: {{
      responsive: true, maintainAspectRatio: false, animation: false,
      plugins: {{ legend: {{ labels: {{ color: '#ccc', font:{{size:11}} }} }} }},
      scales: {{
        x: {{ type:'linear', ticks:{{color:'#999'}}, grid:{{color:'#333'}}, title:{{display:true,text:xlabel,color:'#999'}} }},
        y: {{ ticks:{{color:'#999'}}, grid:{{color:'#333'}}, title:{{display:true,text:ylabel,color:'#999'}} }}
      }}
    }}
  }});
}}

const COLORS = ['#ff8040','#40b0ff','#80ffc0','#ff80c0','#ffd060','#a0a0ff','#ff6060','#60d0d0'];

// Probes
const probeDs = DATA.probes.labels.map((lab, i) => ({{
  label: lab, borderColor: COLORS[i % COLORS.length], backgroundColor: 'transparent',
  data: DATA.probes.times.map((t,k)=>({{x:t, y:DATA.probes.series[i][k]}})),
  tension: 0.2, pointRadius: 0
}}));
makeChart('chartProbes', probeDs, '温度 [°C]');

// Vol avg
makeChart('chartVol', [{{
  label:'volume average T', borderColor:'#ffb96b', backgroundColor:'transparent',
  data: DATA.volAvg.times.map((t,k)=>({{x:t, y:DATA.volAvg.values[k]}})),
  tension: 0.2, pointRadius: 0
}}], '体積平均温度 [°C]');

// Wall heat flux per-patch integrals
if (DATA.wall && DATA.wall.patches) {{
  const wDs = DATA.wall.patches.map((p, i) => ({{
    label: p, borderColor: COLORS[i % COLORS.length], backgroundColor: 'transparent',
    data: DATA.wall.times.map((t,k)=>({{x:t, y:DATA.wall.integral[p][k]}})),
    tension: 0.2, pointRadius: 0
  }}));
  makeChart('chartWall', wDs, '熱流束 integral [W]');

  // Final-time table
  const t = document.getElementById('patchTable');
  let html = '<thead><tr><th>patch</th><th>min [W/m²]</th><th>max [W/m²]</th><th>integral [W]</th></tr></thead><tbody>';
  for (const p of DATA.wall.patches) {{
    const last = DATA.wall.times.length - 1;
    html += `<tr><td>${{p}}</td><td>${{DATA.wall.min[p][last].toFixed(1)}}</td><td>${{DATA.wall.max[p][last].toFixed(1)}}</td><td>${{DATA.wall.integral[p][last].toFixed(1)}}</td></tr>`;
  }}
  html += '</tbody>';
  t.innerHTML = html;
}}
</script>
</body>
</html>
"""


def kpi_card(label: str, value: str, unit: str = "") -> str:
    return f'<div class="kpi"><div class="lbl">{label}</div><div class="num">{value}<span class="unit">{unit}</span></div></div>'


def _subsample_idx(n: int, max_pts: int) -> list[int]:
    if n <= max_pts:
        return list(range(n))
    step = n / max_pts
    return [int(i * step) for i in range(max_pts)]


def generate(case_dir: Path, title: str | None = None, output: Path | None = None,
             max_pts: int = 1500) -> Path:
    pp = case_dir / "postProcessing"
    times, probe_series, labels = parse_probe_T(pp / "probes")
    vol_t, vol_v = parse_vol_average(pp / "volAverageT")
    wall = parse_wall_heat_flux(pp / "wallHeatFlux")

    # Subsample long time series to keep HTML size reasonable
    if times and len(times) > max_pts:
        idx = _subsample_idx(len(times), max_pts)
        times = [times[i] for i in idx]
        probe_series = [[s[i] for i in idx] for s in probe_series]
    if vol_t and len(vol_t) > max_pts:
        idx = _subsample_idx(len(vol_t), max_pts)
        vol_t = [vol_t[i] for i in idx]
        vol_v = [vol_v[i] for i in idx]
    if wall and wall.get("times") and len(wall["times"]) > max_pts:
        idx = _subsample_idx(len(wall["times"]), max_pts)
        wall["times"] = [wall["times"][i] for i in idx]
        for p in wall["patches"]:
            wall["min"][p] = [wall["min"][p][i] for i in idx]
            wall["max"][p] = [wall["max"][p][i] for i in idx]
            wall["integral"][p] = [wall["integral"][p][i] for i in idx]

    # KPIs
    kpis = []
    if probe_series and probe_series[0]:
        finals = [s[-1] for s in probe_series if s]
        kpis.append(kpi_card("プローブ最高(最終)", f"{max(finals):.1f}", "°C"))
        kpis.append(kpi_card("プローブ最低(最終)", f"{min(finals):.1f}", "°C"))
        kpis.append(kpi_card("成層差 ΔT", f"{max(finals)-min(finals):.1f}", "°C"))
    if vol_v:
        kpis.append(kpi_card("体積平均T(最終)", f"{vol_v[-1]:.1f}", "°C"))
    if wall and wall.get("patches"):
        last_idx = len(wall["times"]) - 1
        total_in = sum(wall["integral"][p][last_idx] for p in wall["patches"] if wall["integral"][p][last_idx] > 0)
        total_out = sum(wall["integral"][p][last_idx] for p in wall["patches"] if wall["integral"][p][last_idx] < 0)
        kpis.append(kpi_card("壁面流入合計", f"{total_in:.0f}", "W"))
        kpis.append(kpi_card("壁面放熱合計", f"{total_out:.0f}", "W"))
    t_all = times or vol_t or (wall.get("times") if wall else []) or [0.0]
    kpis.append(kpi_card("シミュレーション終了", f"{t_all[-1]:.1f}", "s"))

    data = {
        "probes": {"times": times, "series": probe_series, "labels": labels},
        "volAvg": {"times": vol_t, "values": vol_v},
        "wall": wall,
    }

    title_str = title or case_dir.name
    html = HTML_TEMPLATE.format(
        title=title_str,
        case_dir=str(case_dir),
        t_start=t_all[0] if t_all else 0,
        t_end=t_all[-1] if t_all else 0,
        n_probes=len(probe_series),
        kpi_html="\n".join(kpis),
        data_json=json.dumps(data),
    )

    out = output or (case_dir / "report.html")
    out.write_text(html, encoding="utf-8")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case_dir", type=Path, help="OpenFOAM case directory (results/<case>)")
    ap.add_argument("--title", default=None, help="Report title (default: case name)")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Output HTML path")
    args = ap.parse_args()
    out = generate(args.case_dir, title=args.title, output=args.output)
    print(f"Report written to: {out}")


if __name__ == "__main__":
    main()
