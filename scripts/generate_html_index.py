"""Generate an index page listing all per-case HTML reports under results/.

Usage:
    python scripts/generate_html_index.py
    # writes results/index.html
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def main() -> None:
    cases = sorted(p.parent.name for p in RESULTS.glob("*/report.html"))
    rows = []
    for c in cases:
        # group by prefix
        rows.append(
            f'<li><a href="{c}/report.html">{c}</a></li>'
        )
    html = f"""<!DOCTYPE html>
<html lang="ja"><head><meta charset="UTF-8"><title>SaunaFlow CFD Reports</title>
<style>
body {{ font-family: -apple-system,"Segoe UI","Hiragino Sans",sans-serif; background:#1a1a1a; color:#e8e8e8; margin:0; padding:24px; }}
h1 {{ color:#ffb96b; font-size:20px; margin:0 0 14px; }}
.meta {{ color:#888; font-size:12px; margin-bottom:18px; }}
ul {{ list-style:none; padding:0; display:grid; grid-template-columns:repeat(auto-fill,minmax(220px,1fr)); gap:8px; }}
li a {{ display:block; background:#242424; border:1px solid #333; border-radius:6px; padding:10px 14px;
       color:#ffb96b; text-decoration:none; font-variant-numeric:tabular-nums; }}
li a:hover {{ background:#2d2d2d; border-color:#555; }}
.tools {{ margin-top:24px; padding-top:14px; border-top:1px solid #333; }}
.tools a {{ color:#80b8ff; }}
</style></head><body>
<h1>SaunaFlow CFD Reports</h1>
<div class="meta">Generated reports under <code>results/</code> ({len(cases)} cases)</div>
<ul>
{chr(10).join(rows)}
</ul>
<div class="tools">
  <h1 style="font-size:14px">Tools</h1>
  <ul>
    <li><a href="../tools/loyly_calculator.html">ロウリュ体感温度シミュレータ (0次元)</a></li>
  </ul>
</div>
</body></html>
"""
    out = RESULTS / "index.html"
    out.write_text(html, encoding="utf-8")
    print(f"Index written to: {out}")


if __name__ == "__main__":
    main()
