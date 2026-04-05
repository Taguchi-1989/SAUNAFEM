---
name: openfoam-debug
description: Debug and fix OpenFOAM buoyantPimpleFoam case setup errors
triggers:
  - "openfoam"
  - "buoyantPimpleFoam"
  - "foam error"
  - "NaN"
  - "p_rgh"
---

# OpenFOAM buoyantPimpleFoam デバッグスキル

buoyantPimpleFoam ケースのエラーを診断・修正するためのスキル。

## 既知のエラーパターンと修正

### NaN in p_rgh
1. `p_rgh` 初期値が 0 → **101325 Pa** に変更
2. `pRefValue` が 0 → **101325** に変更
3. `hRef` ファイルが無い → **作成** (value = 0)

### wallDist method not found
→ `fvSchemes` に `wallDist { method meshWave; }` を追加

### div scheme not found
- `div(phi,h)` — sensibleEnthalpy の場合に必要
- `div(((rho*nuEff)*dev2(T(grad(U)))))` — 圧縮性用
- `div(phi,Yi_h)` — multiComponent の場合

### codedSource name not found
→ トップレベルに `name` キーを追加（Coeffs 内だけでは不足）

### rho solver not found
→ `fvSolution` に `rho { solver PCG; preconditioner DIC; }` を追加

### FPE / sigFpe
1. `FOAM_SIGFPE=false` で実行
2. codedSource でゼロ除算ガード: `max(rho, 0.5)`, `k + SMALL`
3. 初期 deltaT を 0.001 以下に

### externalWallHeatFluxTemperature で発散
→ `fixedValue` に変更して安定化後に戻す

## ケースビルド→実行パイプライン

```bash
# 1. Python でケースビルド
PYTHONPATH=src python -c "
from harness.case_builder import build_case
from pathlib import Path
build_case(Path('configs/cases/dry_sauna_steady.yaml'), output_dir=Path('results/openfoam_dry'))
"

# 2. WSL で実行
wsl -d Ubuntu -- /usr/bin/openfoam2312 bash -c "
cd /mnt/d/dev/SaunaFEM/results/openfoam_dry
export FOAM_SIGFPE=false
blockMesh
checkMesh
foamDictionary system/controlDict -entry endTime -set 300
buoyantPimpleFoam
"

# 3. 結果プロット
python scripts/plot_openfoam_results.py
```

## 参照
- `docs/openfoam_troubleshooting.md` — 詳細なエラー記録
- `foam_templates/base_case/` — Jinja2 テンプレート
- `src/harness/case_builder.py` — ケースビルダー
