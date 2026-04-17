# SaunaFlow OpenFOAM計算キャンペーン技術報告書

## 1. プロジェクト概要

**SaunaFlow** は、フィンランド式ドライサウナの熱環境をCFDシミュレーション可能にするPython駆動型ハーネスです。

**計算条件**
- **計算エンジン**: OpenFOAM v2312 (buoyantSimpleFoam / buoyantPimpleFoam)
- **ハーネス**: Python 3.11+ 駆動の宣言的YAML/JSON定義
- **対象環境**: Ubuntu Linux (WSL2経由でWindows実行可)
- **目標**: 18kW電気ヒーター搭載の3.0m×2.5m×2.5mドライサウナ室における温度成層化の再現

**基本物理モデル**
- 密度駆動対流(buoyancy-driven flow)
- 壁面熱伝達(externalWallHeatFluxTemperature BC)
- ヒーター熱源(surface flux またはvolume source)
- 乱流熱輸送(k-epsilon, enthalpy-based)
- オプション: 放射伝熱、通風システム、Aufguss(時間的蒸気導入)

---

## 2. 本セッションで実装・検証した機能

### 2.1 熱収支自動集約パイプライン (Phase 1)

OpenFOAM postProcessing出力から熱エネルギー収支を自動抽出・集計するパイプラインを実装しました。

**実装内容**
- **controlDict.j2**: 2つのfunction object追加
  - `wallHeatFlux`: 各パッチの壁面熱フラックス [W] を統計
  - `volAverageT`: ドメイン全体の体積平均温度 [K] を計算

- **heat_balance_parser.py** (`src/harness/heat_balance_parser.py`)
  - OpenFOAM v2312 postProcessing形式に対応
    - Combined format: `wallHeatFlux.dat` (行形式: time, patch, min, max, integral)
    - Per-patch形式: `floor_wallHeatFlux.dat`, `ceiling_wallHeatFlux.dat` など
    - Legacy形式: `surfaceFieldValue.dat`
  - `parse_wall_heat_flux()`: パッチごとのフラックス時系列を抽出
  - `parse_vol_average_t()`: `volFieldValue.dat` から体積平均T時系列を抽出

- **HeatBalance データクラス**
  - `heater_input_W`: ヒーター面からの入熱 [W]
  - `wall_loss_W`: 壁面からの出熱 [W]
  - `vent_loss_W`: 通風口からの出熱 [W]
  - `imbalance_W`, `imbalance_pct`: エネルギー収支の不均衡
  - `vol_avg_T`: 体積平均温度 [K]
  - `patch_fluxes`: パッチ別フラックス辞書

- **reporting.py**: `heat_balance_to_markdown()`関数
  - HeatBalanceオブジェクトをMarkdownテーブルに変換
  - 各パッチのフラックス内訳を表示
  - 不均衡度をパーセンテージで表示

**テスト** (14テスト合格)
- `tests/unit/test_heat_balance_parser.py`
  - HeatBalance計算(imbalance, zero-division safety)
  - Per-patch形式と combined形式の両方をパース
  - `_wallHeatFlux` suffix自動除去(v2312互換)
  - 複数time directoriesの処理

### 2.2 ヒーターモデルA/B比較機能 (Phase 2)

異なるヒーター熱源モデルの効果を比較可能にしました。

**ヒーターモデル選択肢**

| モデル | BC種別 | 実装場所 | 特徴 |
|--------|--------|---------|------|
| `surface_flux` | externalWallHeatFluxTemperature mode=flux | T.j2 | 60,000 W/m² の flux BC を直接指定 |
| `volume_source` | scalarSemiImplicitSource (enthalpy) | fvOptions.j2 | cellZone 内に distributed volume heat source |

**スキーマ拡張** (`configs/schemas/case_schema.json`)

heater.model オプション: surface_flux | volume_source、heater.depth 追加

**テンプレート変更**

1. **T.j2** (0/T.j2) - heater_wall パッチで条件分岐
   - surface_flux → flux BC (q_uniform = heat_flux)
   - volume_source → coefficient BC (h=1.0, Ta=293.15)

2. **fvOptions.j2** (constant/fvOptions.j2)
   - heaterSource with scalarSemiImplicitSource
   - cellZone heaterZone, volumeMode specific

3. **topoSetDict.j2** (system/topoSetDict.j2) - 新規ファイル
   - volume_source mode の場合、boxToCell で heaterZone cellZone を定義

4. **case_builder.py** - コンテキスト変数生成
   - heater_power_density = heater_power_W / (heater_depth * heater_height * heater_width)

**YAML設定ファイル**

- dry_sauna_steady_surfflux.yaml (Case A): surface_flux モデル
- dry_sauna_steady_volsource.yaml (Case B): volume_source、通風なし
- dry_sauna_steady_vent.yaml (Case C): volume_source + ventilation

**テスト** (46テスト合格)
- heater_model context variable生成
- heater_power_density計算
- topoSetDict生成(cellZone定義)
- T.j2条件分岐検証

### 2.3 安定性改善

| 問題 | 根本原因 | 対策 |
|------|---------|------|
| SIMPLE solver発散 | T(温度)の緩和係数が高い | enthalpy h に変更、緩和係数を0.2に低下 |
| 圧力参照値エラー | pRefValue=0 (相対圧) | pRefValue=101325 (絶対圧) に変更 |
| Residual control対象ミス | residualControl T を監視 | residualControl h に変更 |
| surface_flux mesh不安定 | M0では局所勾配過大 | uniform_init: true で非均一初期化を無効化 |

### 2.4 通風システム対応

Case C: supply_vent + exhaust_vent
- T.j2 の ventilation ブロック
- batch.py での OpenFOAM A/B 比較機能

---

## 3. 計算結果と知見

### 3.1 Case A: surface_flux (60,000 W/m²)

| 項目 | 値 |
|------|-----|
| ソルバー | buoyantSimpleFoam |
| メッシュ | M0 (76,800 cells) |
| 結果 | **発散** (iter 4-13) |

**原因**: 60,000 W/m² flux が M0 メッシュの局所セル境界層で extreme gradient を生成
**結論**: surface_flux モデルは M1+ mesh (16 cells/m) が必須

### 3.2 Case B: volume_source (18 kW, no ventilation)

| 項目 | 値 |
|------|-----|
| ソルバー | buoyantSimpleFoam |
| メッシュ | M0 (76,800 cells) |
| 反復数 | 19,293 (iter 19k で発散) |

**iter 15,000 での結果**

| 計測点 | 値 | 備考 |
|--------|-----|------|
| upper_bench | 381 K (108°C) | 測定値 80-100°C より高い |
| lower_bench | 341 K (68°C) | 温度成層化は再現 |
| floor_level | 304 K (31°C) | - |
| Vol-avg T | 356 K (83°C) | - |

**熱収支 (iter 15,000)**

| コンポーネント | 値 [W] | % of 18kW |
|--------|------|-----------|
| Heater input | +18,000 | 100% |
| Wall loss | -6,717 | -37% |
| Imbalance | +11,283 | +63% |

**解釈**: Wall loss 37% → **steady state 未達**、残り63%は air mass を加熱中

### 3.3 Case C: volume_source + 通風

| 項目 | 値 |
|------|-----|
| 反復数 | 30,000 **安定実行** (発散なし) |

**iter 30,000 での結果**

| 計測点 | 値 | vs. Case B |
|--------|-----|-----------|
| upper_bench | 344 K (70°C) | -37°C |
| Vol-avg T | 316 K (43°C) | -40°C |

**熱収支 (iter 30,000)**

| コンポーネント | 値 [W] | % of 18kW |
|--------|------|-----------|
| Heater input | +18,000 | 100% |
| Wall loss | -7,500 | -42% |
| Vent loss | -1,850 | -10% |
| Imbalance | +8,650 | +48% |

**観察**: 通風により -38°C の大幅な温度低下、最大反復数まで安定

---

## 4. 欠けている物理現象と定量的見積

| 物理現象 | 推定大きさ [W] | % of 18kW | 空気温度への効果 | 実装状況 |
|---------|--------------|----------|----------------|---------|
| **放射伝熱** (ヒーター→壁) | 3,600-5,200 | 20-29% | 空気T低下 | Template exists, disabled |
| **壁面熱容量** (transient) | 5,000-7,000 | 28-39% | warm-up中の冷却 | CFD では未実装 |
| **漏気 + 通風** | 1,700-2,800 | 9-16% | 空気T低下 | Case C で実装 |
| **木材水分蒸発** | 3,000-5,400 | 17-30% | Transient のみ | 未実装 |

**主な観察**
- 放射伝熱 20-29% が最大の欠落 → air temperature 低下へ
- Case B 108°C で放射を無視することが実測に近い可能性
- Wall thermal mass は steady state CFD では考慮外

---

## 5. 現在の問題と改善提案

### 問題1: Steady state 到達せず
- Case B: iter 15k で wall loss 37%
- Case C: iter 30k でも heat loss 52%
- 原因: SIMPLE + h=0.2 relaxation の slow convergence
- **対策**: h=0.1, 50k+ iter 実行、またはbuoyantPimpleFoam transient

### 問題2: Case A (surface_flux) M0 で不安定
- **対策**: M1 mesh に upgrade、または flux density 低下

### 問題3: 放射伝熱モデル未実装
- 実測値との最大乖離要因
- **対策**: Case D で P1 model 追加

### 問題4: scalarCodedSource (buoyancy_production) wmake 依存
- Package OpenFOAM では pre-compiled library として提供なし
- **対策**: Built-in type への置き換え

---

## 6. 次のアクション (優先度順)

### Priority A: 数値的収束の確保

1. Case C with h=0.1: iter 50,000+ 実行
2. buoyantPimpleFoam transient試験

### Priority B: 欠落物理の追加

3. Radiation model 実装 (Case D)
4. Case A on M1 mesh

### Priority C: 検証・最適化

5. 通風ACH検証、vent area調整
6. Wall thermal mass モデル
7. scalarCodedSource 代替

---

## 7. 実装ファイル一覧

### 新規作成ファイル

| ファイル | 説明 |
|---------|------|
| `src/harness/heat_balance_parser.py` | Heat balance aggregation pipeline |
| `tests/unit/test_heat_balance_parser.py` | Heat balance parser tests (14 tests) |
| `foam_templates/base_case/system/topoSetDict.j2` | cellZone heaterZone definition |
| `configs/cases/dry_sauna_steady_surfflux.yaml` | Case A config (surface flux) |
| `configs/cases/dry_sauna_steady_volsource.yaml` | Case B config (volume source) |
| `configs/cases/dry_sauna_steady_vent.yaml` | Case C config (volume + ventilation) |

### 変更ファイル

| ファイル | 変更内容 |
|---------|---------|
| `foam_templates/base_case/system/controlDict.j2` | wallHeatFlux, volAverageT functionObjects追加 |
| `foam_templates/base_case/system/fvSolution.j2` | relaxation h=0.2, pRefValue=101325 |
| `foam_templates/base_case/0/T.j2` | heater_wall conditional (surface_flux vs volume_source) |
| `foam_templates/base_case/constant/fvOptions.j2` | heaterSource scalarSemiImplicitSource 追加 |
| `src/harness/case_builder.py` | heater_model, heater_power_density context |
| `src/harness/batch.py` | parse_openfoam_case, compare_openfoam_results 追加 |
| `src/harness/reporting.py` | heat_balance_to_markdown() 追加 |
| `configs/schemas/case_schema.json` | heater.model enum, uniform_init, wall_htc_override 追加 |
| `tests/unit/test_case_builder.py` | heater model tests (46 tests total) |
| `tests/unit/test_batch.py` | Batch runner tests (A/B comparison) |
| `tests/unit/test_reporting.py` | heat_balance_to_markdown tests (6 tests) |
| `scripts/run_openfoam_wsl.sh` | Generic OpenFOAM runner with heat balance output |
| `scripts/run_case_bc.sh` | Case B vs C comparison runner |

### テスト統計

```
tests/unit/test_heat_balance_parser.py     14 tests ✓
tests/unit/test_reporting.py               6 tests  ✓
tests/unit/test_case_builder.py            46 tests ✓
Total unit tests: 211+ passing
```

---

## 8. デプロイメント・実行手順

### ケース構築

```bash
# Case B
PYTHONPATH=src python -c "
from harness.case_builder import build_case
from pathlib import Path
build_case(
    Path('configs/cases/dry_sauna_steady_volsource.yaml'),
    output_dir=Path('results/openfoam_case_b')
)"
```

### WSL2 上での計算実行

```bash
wsl -d Ubuntu -- bash /mnt/d/dev/SaunaFEM/scripts/run_openfoam_wsl.sh \
  /mnt/d/dev/SaunaFEM/results/openfoam_case_b \
  --toposet
```

### テスト実行

```bash
pytest tests/unit/ -q
pytest tests/unit/test_heat_balance_parser.py -v
```

---

## 9. 技術的な注記と制限

### OpenFOAM v2312 互換性

**対応フォーマット**
- `wallHeatFlux.dat`: combined format
- Per-patch files: `floor_wallHeatFlux.dat` など
- `volFieldValue.dat`: time-value 形式

**既知の制限**
- scalarCodedSource は pre-compiled library 依存
- Package OpenFOAM では disabled by default

### メッシュレベル

| Level | Cells | Resolution | Use case |
|-------|-------|-------------|----------|
| M0 | 76,800 | 8 cells/m | Quick iteration |
| M1 | 500,000-600,000 | 16 cells/m | Production |
| M2 | 1,500,000 | 32 cells/m | High-fidelity |

### Relaxation factors

**Recommended for stable SIMPLE convergence**
- h (enthalpy): 0.1 - 0.2 (default 0.2)
- p (pressure): 0.3
- U (velocity): 0.6

---

## 10. 参照・アーティファクト

**主要ドキュメント**
- `CLAUDE.md` (Project Instructions)
- `docs/openfoam_troubleshooting.md` (Error log, checklist)

**スクリプト**
- `scripts/run_openfoam_wsl.sh`: WSL2 generalized OpenFOAM runner
- `scripts/run_case_bc.sh`: Case B vs C comparison automation

---

## 11. 結語

**成果**
- Heat balance auto-aggregation pipeline の完全実装・テスト (14 tests)
- Heater model A/B comparison framework の構築 (3 YAML cases, 46 tests)
- OpenFOAM v2312 postProcessing 統合
- 安定性改善 (enthalpy, relaxation, pRefValue)

**現在のボトルネック**
- SIMPLE solver の slow convergence (steady state 未達)
- 放射伝熱モデル不在 (20-29% missing)
- Surface flux model の mesh依存性 (M0では不安定)

**次ステップ**
1. Case C extended run (h=0.1, 50k iter) で steady state 接近度確認
2. Radiation model 追加 (Case D)
3. buoyantPimpleFoam transient alternative 検証
4. M1 mesh upgrade で surface_flux viability test

実装は全て version control下にあり、reproducible・extendable な構造です。

---

**作成日**: 2026-04-17  
**ハーネスバージョン**: SaunaFlow Phase 2 (Heat Balance + Heater A/B)  
**計算ドメイン**: 3.0m × 2.5m × 2.5m ドライサウナ, 18kW electric heater
