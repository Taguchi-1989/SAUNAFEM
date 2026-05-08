"""Build a single self-contained HTML report of the skin heat-flux work.

Re-runs (or reuses outputs of) the three demos and assembles the result
into one portable HTML file with all plots base64-embedded — no external
dependencies, opens offline in any browser.

Usage:
    PYTHONPATH=src python scripts/skin_htc_report.py
    PYTHONPATH=src python scripts/skin_htc_report.py --output results/report.html
    PYTHONPATH=src python scripts/skin_htc_report.py --skip-rebuild   # reuse existing PNGs
"""

from __future__ import annotations

import argparse
import base64
import datetime as _dt
import subprocess
import sys
from html import escape
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

SWEEP_DIR = REPO_ROOT / "results" / "skin_htc"
SCEN_DIR = REPO_ROOT / "results" / "skin_htc_scenarios"
TS_DIR = REPO_ROOT / "results" / "skin_htc_timeseries"


def _embed_png(path: Path) -> str:
    if not path.exists():
        return f"<p><em>(missing: {escape(str(path))})</em></p>"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" alt="{escape(path.name)}" />'


def _figure(path: Path, caption: str) -> str:
    return (
        '<figure>'
        f'{_embed_png(path)}'
        f'<figcaption>{caption}</figcaption>'
        '</figure>'
    )


def _rebuild_outputs() -> None:
    scripts = [
        "scripts/skin_htc_demo.py",
        "scripts/skin_htc_scenarios.py",
        "scripts/skin_htc_timeseries.py",
    ]
    env_cmd_prefix = ["env", f"PYTHONPATH={REPO_ROOT / 'src'}", sys.executable]
    for script in scripts:
        print(f"Rebuilding via {script} …")
        subprocess.run(env_cmd_prefix + [str(REPO_ROOT / script)],
                       check=True, cwd=REPO_ROOT)


CSS = """
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Hiragino Sans", "Yu Gothic",
               "Segoe UI", Roboto, sans-serif;
  line-height: 1.65;
  color: #222;
  background: #fafaf8;
  margin: 0;
}
.container { max-width: 1100px; margin: 0 auto; padding: 2rem 1.5rem 4rem; }
header { border-bottom: 3px solid #444; padding-bottom: 1rem; margin-bottom: 2rem; }
h1 { font-size: 1.9rem; margin: 0 0 0.4rem; }
header .meta { color: #666; font-size: 0.9rem; }
nav { background: #f1ede4; border-radius: 6px; padding: 0.8rem 1.2rem; margin: 1.2rem 0 2rem; }
nav ol { margin: 0; padding-left: 1.2rem; columns: 2 220px; column-gap: 1.5rem; }
nav a { color: #1d4d6b; text-decoration: none; }
nav a:hover { text-decoration: underline; }
section { margin: 2.4rem 0; }
h2 { color: #1d4d6b; border-left: 4px solid #1d4d6b; padding-left: 0.6rem;
     margin-top: 2.4rem; }
h3 { color: #294b29; margin-top: 1.6rem; }
figure { margin: 1rem 0; text-align: center; background: #fff;
         border: 1px solid #e5e0d4; border-radius: 6px; padding: 0.5rem; }
figure img { max-width: 100%; height: auto; display: block; margin: 0 auto; }
figcaption { font-size: 0.85rem; color: #555; margin-top: 0.4rem; padding: 0 0.5rem 0.3rem; }
table { border-collapse: collapse; width: 100%; margin: 1rem 0;
        font-size: 0.92rem; background: #fff; }
th, td { border: 1px solid #d8d3c4; padding: 0.4rem 0.6rem; text-align: right; }
th { background: #efe9d8; }
td.left, th.left { text-align: left; }
.formula { background: #fff; border-left: 3px solid #c8a04a;
           padding: 0.6rem 1rem; margin: 1rem 0; font-family: ui-monospace, monospace;
           overflow-x: auto; }
.note { background: #fff7e0; border-left: 3px solid #c8a04a;
        padding: 0.6rem 1rem; margin: 1rem 0; font-size: 0.92rem; }
.callout { background: #e8f1ee; border-left: 3px solid #2a8b6e;
           padding: 0.8rem 1.2rem; margin: 1.2rem 0; }
code { font-family: ui-monospace, "SF Mono", monospace; font-size: 0.9em;
       background: #f1ede4; padding: 0.05rem 0.3rem; border-radius: 3px; }
.section-summary { font-style: italic; color: #555; margin-bottom: 1rem; }
ul.kv { list-style: none; padding-left: 0; }
ul.kv li { margin: 0.15rem 0; }
ul.kv li code { display: inline-block; min-width: 8em; }
.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }
@media (max-width: 700px) { nav ol { columns: 1; } .two-col { grid-template-columns: 1fr; } }
footer { margin-top: 3rem; border-top: 1px solid #d8d3c4; padding-top: 1rem;
         color: #888; font-size: 0.85rem; }
"""


def _section_physics() -> str:
    return """
<section id="physics">
  <h2>1. 物理モデル</h2>
  <p class="section-summary">皮膚表面 1 点での熱収支を 3 成分に分解。
     入力は乾球温度・相対湿度・局所風速・平均放射温度の 4 つだけ。</p>

  <h3>1.1 熱収支</h3>
  <p>サインルール: q &gt; 0 = 皮膚に熱が入る方向。</p>
  <div class="formula">
    q_total = q_conv + q_rad + q_evap   [W/m²]<br>
    q_conv  = h_c · (T_air − T_skin)<br>
    q_rad   = h_r · (T_mrt − T_skin)<br>
    q_evap  = h_e · w_skin · (P_air − P_sat,skin)
  </div>

  <h3>1.2 各熱伝達係数</h3>
  <ul class="kv">
    <li><code>h_c(V, 部位)</code> 自然対流 (Mitchell 1974) と強制対流 (de Dear 1997)
        の Churchill ブレンド (n=3)</li>
    <li><code>h_r(T_mrt, T_skin)</code> 線形化グレイボディ放射
        h_r = ε σ (T_mrt² + T_skin²)(T_mrt + T_skin)</li>
    <li><code>h_e = 16.5 · h_c</code> Lewis 関係 [W/(m²·kPa)]</li>
  </ul>
  <p>典型値（87°C 環境・無風）:
     h_c ≈ 6.5、h_r ≈ 8.4、h_e ≈ 107。</p>

  <h3>1.3 部位別係数 (de Dear 1997)</h3>
  <p>強制対流 h_forc = a · V<sup>b</sup>:</p>
  <table>
    <tr><th class="left">部位</th><th>a</th><th>b</th></tr>
    <tr><td class="left">face</td><td>14.0</td><td>0.61</td></tr>
    <tr><td class="left">chest</td><td>10.4</td><td>0.56</td></tr>
    <tr><td class="left">back</td><td>9.1</td><td>0.61</td></tr>
    <tr><td class="left">arm</td><td>13.4</td><td>0.55</td></tr>
    <tr><td class="left">thigh</td><td>10.6</td><td>0.59</td></tr>
    <tr><td class="left">calf</td><td>11.6</td><td>0.59</td></tr>
    <tr><td class="left">average</td><td>10.4</td><td>0.60</td></tr>
  </table>

  <h3>1.4 q [W/m²] と Q [J/m²] の関係</h3>
  <div class="formula">
    Q(t) = ∫₀ᵗ q(τ) dτ   [J/m²]<br>
    1 杯 (100 mL) 蒸発潜熱: L · m = 2.26 MJ/kg × 0.1 kg = 226 kJ<br>
    体表 1.8 m² に全凝縮なら 126 kJ/m² が「ロウリュ 1 杯の物理上限」
  </div>
  <div class="callout">
    <strong>瞬時量 q だけでは「1 杯のロウリュが何を意味するか」は分からない。</strong>
    時間積分 Q を取ると、各イベントが潜熱バジェット（126 kJ/m²/杯）に対して
    どこに乗るかが見え、エネルギー収支として読める。
  </div>
</section>
"""


def _section_static() -> str:
    return f"""
<section id="static">
  <h2>2. 静的シナリオ — サウナタイプ × 処置</h2>
  <p class="section-summary">7 種のサウナタイプ × 6 処置（無処置 / 1〜2 杯 / Aufguss /
     組合せ）の定常状態行列。</p>

  {_figure(SCEN_DIR / "sauna_types_baseline.png",
           "図 2-1. 5 タイプのサウナ素のままでの皮膚熱流束分解。"
           "Turkish hammam は気温が低くても凝縮が支配的。")}

  {_figure(SCEN_DIR / "treatment_matrix.png",
           "図 2-2. サウナタイプ × 処置 のヒートマップ。"
           "左: 全身熱流束 q_total。右: 体感等価温度 T_equivalent。"
           "ロウリュ + Aufguss の組合せで両軸とも最大値。")}

  {_figure(SCEN_DIR / "loyly_progression.png",
           "図 2-3. 1 杯増えると q_total が何 W/m² 上がるか。"
           "calm（青）と Aufguss（橙）の差は風速依存。")}

  {_figure(SCEN_DIR / "face_focus.png",
           "図 2-4. 顔面 vs 全身。Aufguss + 高湿度時に顔面が突出する。")}
</section>
"""


def _section_sweating() -> str:
    return f"""
<section id="sweating">
  <h2>3. 発汗（皮膚濡れ率 w_skin）の影響</h2>
  <p class="section-summary">発汗は P_vapor &lt; P_sat,skin（低湿度域）でのみ
     冷却として効く。湿度が高くなり凝縮側に入ると発汗の有無は無関係。</p>

  {_figure(SWEEP_DIR / "sweat_q_total_vs_RH.png",
           "図 3-1. q_total vs RH を w_skin 5 段階で。下段ズームで"
           "クロスオーバー RH（凝縮開始点）と発汗冷却の差が明瞭。")}

  {_figure(SWEEP_DIR / "sweat_q_evap_vs_RH.png",
           "図 3-2. q_evap だけ抜き出し。下段で物理上限 −400 W/m² で"
           "頭打ちになる「発汗能力の限界」が見える。")}

  {_figure(SWEEP_DIR / "sweat_dry_vs_wet_bars.png",
           "図 3-3. 5 シナリオで w_skin = 0（皮膚乾燥）vs 1（全身発汗）の比較。"
           "低湿サウナでは発汗が q_total を 30〜50% 削るが、"
           "ロウリュ後は同値に収束。")}

  {_figure(SWEEP_DIR / "per_part_bars.png",
           "図 3-4. 部位別 q_total。Aufguss 時は顔面が他部位より明確に高い。")}
</section>
"""


def _section_timeseries() -> str:
    return f"""
<section id="timeseries">
  <h2>4. 時系列 — q(t), T_eq(t), Q(t)</h2>
  <p class="section-summary">5 つの時間発展イベント（純ベースライン、1 杯ロウリュ、
     3 段 Aufguss、しきじ薬草、加賀屋プログレッシブ）を整数秒刻みで
     積分。Q = ∫ q dt で「セッション全体の被熱量」を定量化。</p>

  <h3>4.1 比較ダッシュボード（主力画 1）</h3>
  {_figure(TS_DIR / "comparison_dashboard.png",
           "図 4-1. 5 イベントを軸共有で重ね描き。一画でパターン比較が可能 — "
           "Shikiji の平坦高原、Aufguss のスパイク、加賀屋の階段、"
           "Single löyly の小山、Finnish dry の横ばい。"
           "Q_cum パネルの黒点線は 1 杯の潜熱バジェット 126 kJ/m²。")}

  <h3>4.2 湿度パラドックス（主力画 2）</h3>
  {_figure(TS_DIR / "humidity_paradox.png",
           "図 4-2. 「Hot dry 90°C / RH 5%」vs「Cool humid 60°C / RH 95%」。"
           "気温が 30°C 低いのに皮膚熱流束は 2 倍、体感温度は +68°C。"
           "凝縮が支配する物理的説明をテキストパネルで明記。")}

  <h3>4.3 詳細パネル</h3>
  {_figure(TS_DIR / "events_q_t.png",
           "図 4-3. q(t) の個別パネル。Single löyly: 700 → 747 で減衰。"
           "Aufguss: 800 → 1260 → 1620 → 1898 W/m² と階段状にエスカレート。")}

  {_figure(TS_DIR / "events_t_eq.png",
           "図 4-4. T_equivalent(t)。しきじ薬草の体感 153.5°C は 60°C 環境で"
           "凝縮潜熱が支配することの帰結。")}

  {_figure(TS_DIR / "events_Q_cumulative.png",
           "図 4-5. Q_cum(t) = ∫ q dt。傾きが大きいほど被熱率が高い。"
           "しきじが最も急、加賀屋が二次関数的に加速。")}

  {_figure(TS_DIR / "aufguss_components.png",
           "図 4-6. Aufguss セッションの成分分解。風（Aufguss 窓）が来るたび "
           "q_conv が四角く立ち上がり、ロウリュの凝縮成分（緑）が"
           "それに上乗せ。q_rad はほぼ一定。")}

  <h3>4.4 1 杯あたりのエネルギー</h3>
  {_figure(TS_DIR / "energy_per_ladle.png",
           "図 4-7. 左: 全イベント合計 Q [kJ/m²] と 1.8 m² 換算の総 kJ。"
           "右: 1 杯あたりの実吸収 Q (kJ/m²) を潜熱バジェット 126 kJ/m² と比較。"
           "実測 115〜129 kJ/m²/杯 で物理上限ぎりぎりに張り付く — "
           "ロウリュの蒸発潜熱はほぼすべて皮膚（と壁）に降ってきている。")}
</section>
"""


def _section_results_table() -> str:
    rows = [
        # (event, dur, q_max, T_eq_pk, Q_total, Q_per_pour, Q_body)
        ("Finnish dry 90°C",       300, 761,  85.9, 228, "—",   411),
        ("Single löyly (1 ladle)", 180, 747,  86.0, 129, "128.6",232),
        ("Aufguss session",        300, 1898, 96.0, 345, "115.1",622),
        ("Shikiji herbal",         300, 1520, 153.5,456, "—",   821),
        ("Kagaya progressive",     420, 1575, 134.8,467, "116.9",841),
    ]
    body = "".join(
        f"<tr><td class='left'>{escape(name)}</td>"
        f"<td>{d}</td><td>{q}</td><td>{te}</td>"
        f"<td>{Q}</td><td>{qp}</td><td>{qb}</td></tr>"
        for (name, d, q, te, Q, qp, qb) in rows
    )
    return f"""
<section id="summary">
  <h2>5. 主要数値サマリ</h2>
  <table>
    <tr>
      <th class="left">イベント</th>
      <th>持続 [s]</th>
      <th>q_max [W/m²]</th>
      <th>T_eq ピーク [°C]</th>
      <th>Q_total [kJ/m²]</th>
      <th>Q / 杯 [kJ/m²]</th>
      <th>Q (1.8 m²) [kJ]</th>
    </tr>
    {body}
  </table>
  <p class="note"><strong>参照</strong>: 1 杯（100 mL）潜熱 = 226 kJ、
     体表 1.8 m² 換算で 126 kJ/m² が物理上限。</p>
</section>
"""


def _section_conclusions() -> str:
    return """
<section id="conclusions">
  <h2>6. まとめ</h2>
  <ul>
    <li><strong>体感熱流束は 3 成分の合成</strong>: 対流 (h_c · ΔT)、
        放射 (h_r · ΔT_mrt)、湿度 (h_e · w · ΔP)。
        どれか 1 つでは説明できない。</li>
    <li><strong>湿度パラドックス</strong>: 60°C / 95% RH は 90°C / 5% RH の
        <strong>2 倍の熱流束</strong>を皮膚に与える。
        体感等価温度では <strong>+68°C</strong>。
        凝縮潜熱の威力が定量化された。</li>
    <li><strong>Aufguss の体感上昇</strong>: 風で h_c が 6 → 17 W/(m²·K) と
        2.8 倍になり、対流成分が 340 → 1200 W/m² に。
        ロウリュ重ねがけで T_eq が 88 → 96°C と階段状に上がる。</li>
    <li><strong>Q は正しい単位</strong>: 1 杯の潜熱バジェット
        126 kJ/m² に対し、実吸収 Q は 115〜129 kJ/m² で張り付く。
        瞬時 q だけでは「1 杯」の意味は計れない。</li>
    <li><strong>発汗の有効域は限定的</strong>: 低湿度で −350〜−400 W/m² の
        冷却を提供するが、凝縮側に入ると無効化される。
        サウナ内部設計では、発汗能を超えない湿度・風速組合せを意識する余地。</li>
  </ul>

  <h2>7. ファイル / 再現方法</h2>
  <ul class="kv">
    <li><code>src/harness/skin_htc.py</code> 純関数モジュール</li>
    <li><code>tests/unit/test_skin_htc.py</code> 単体テスト 29 件</li>
    <li><code>scripts/skin_htc_demo.py</code> パラメータスイープ</li>
    <li><code>scripts/skin_htc_scenarios.py</code> 静的シナリオマトリクス</li>
    <li><code>scripts/skin_htc_timeseries.py</code> 時間発展イベント</li>
    <li><code>scripts/skin_htc_report.py</code> このレポート生成器</li>
  </ul>
  <p>再現:</p>
  <div class="formula">
    PYTHONPATH=src python scripts/skin_htc_report.py<br>
    # → results/skin_htc_report.html
  </div>
</section>
"""


def build_html() -> str:
    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    nav = """
<nav>
  <ol>
    <li><a href="#physics">物理モデル</a></li>
    <li><a href="#static">静的シナリオ</a></li>
    <li><a href="#sweating">発汗の影響</a></li>
    <li><a href="#timeseries">時系列イベント</a></li>
    <li><a href="#summary">主要数値サマリ</a></li>
    <li><a href="#conclusions">まとめ</a></li>
  </ol>
</nav>
"""
    body = (
        _section_physics()
        + _section_static()
        + _section_sweating()
        + _section_timeseries()
        + _section_results_table()
        + _section_conclusions()
    )
    return f"""<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>サウナ皮膚熱流束レポート</title>
<style>{CSS}</style>
</head>
<body>
<div class="container">
  <header>
    <h1>サウナ皮膚熱流束レポート</h1>
    <div class="meta">SaunaFlow PoC — skin_htc モジュール / 生成: {now}</div>
    <p>サウナ室内の人体皮膚にどれだけの熱が、どのような物理機構で
       入ってくるかを、温度・湿度・風速・放射の 4 入力だけから
       スタンドアロン計算するモジュールの結果集。CFD ソルバには依存せず、
       純関数として独立。</p>
  </header>
  {nav}
  {body}
  <footer>
    Generated by <code>scripts/skin_htc_report.py</code>.
    All plots base64-embedded — this HTML is a single self-contained file.
  </footer>
</div>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output", type=Path,
                   default=REPO_ROOT / "results" / "skin_htc_report.html")
    p.add_argument("--skip-rebuild", action="store_true",
                   help="Reuse existing PNG outputs without re-running demos.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.skip_rebuild:
        _rebuild_outputs()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    html = build_html()
    args.output.write_text(html, encoding="utf-8")
    size_kb = args.output.stat().st_size / 1024
    print(f"\nWrote {args.output} ({size_kb:.0f} kB)")


if __name__ == "__main__":
    main()
