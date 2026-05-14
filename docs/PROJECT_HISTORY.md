# SaunaFlow プロジェクト経緯メモ

このドキュメントは「これまで何をやってきたか」をざっと振り返るためのインデックス。
詳細は各章末尾のリンク先 doc にある。日付は最終更新コミットの目安。

---

## タイムライン

### Phase 0 — 基盤（2026-04-03 〜 04-04）
- リポジトリ初期化、ディレクトリ構成、YAML スキーマ、CLI 雛形、pytest 基盤。
- 「ハーネスが主、OpenFOAM は計算エンジン」という設計方針を確定。
- 詳細: [CLAUDE.md](../CLAUDE.md), [要件定義書.md](../要件定義書.md)

### Phase 1 — ドライサウナ定常（2026-04-04 〜 04-05）
- `simple_solver.py`（2-Zone プルーム模型）と OpenFOAM `buoyantSimpleFoam` の両系統を実装。
- 物理モデル刷新（Zukoski プルーム + SST k-omega + 2-zone）。
- 詳細: [technical_reference.md](technical_reference.md), [governing_equations.md](governing_equations.md)

### Phase 2 タスク A/B — ロウリュ・放射（2026-04-05 〜 04-09）
- 蒸気体積膨張（löyly）と幾何 viewFactor 放射モデルを導入。
- ヒーター BC を何度も入れ替え（fixedValue → externalWallHeatFlux → scalarCodedSource → fixedT サロゲート）。
- fvDOM / viewFactor / boundaryRadiationProperties を整備。
- 詳細: [openfoam_troubleshooting.md](openfoam_troubleshooting.md), [openfoam_current_issues.md](openfoam_current_issues.md)

### 過渡・熱収支セッション（2026-04-10 〜 04-15）
- `buoyantPimpleFoam` 過渡対応、時系列レポート。
- 熱収支の自動集計、ヒーター A/B 比較、換気ケース、kMin=0.01 で乱流崩壊対策。
- 詳細: コミット `6ca57dc`〜`2c9d694`、`MEMORY.md` の heat_balance_session メモ。

### パラメトリックスタディ（2026-04-15 〜 04-22）
- Case G/H/I（壁厚、ヒーター温度 473/573/673 K、放射 ON/OFF）。
- Case J/K（換気: pressIO → fixedValue、両側/片側、kW スケール）。
- Case L-1（PIMPLE 過渡 13 kW）で **上層 95℃ / 下層 54℃** と実測目標域に到達。
- 詳細: [parametric_study_wall_heater.md](parametric_study_wall_heater.md), [parametric_study_volume_source_vent.md](parametric_study_volume_source_vent.md), [parametric_study_transient_L1.md](parametric_study_transient_L1.md), [computational_insights_summary.md](computational_insights_summary.md)

### 実測フェーズ準備（2026-05-05 〜 05-13）
- L-2 (18kW) ケース設定、計測フィロソフィー策定、センサ・ファーム実装（`firmware/`）。
- ハードウェア調達計画と GPT 向け調達リサーチプロンプトを整備。
- 詳細: [measurement_philosophy.md](measurement_philosophy.md), [measurement_implementation_plan.md](measurement_implementation_plan.md), [procurement_plan.md](procurement_plan.md)

---

## 既知の未解決事項

- viewFactor 放射と SIMPLE+換気の組み合わせで温度がやや高めに張り付くケースあり（Case I 系）。
- CFD ↔ 実測の突き合わせは未着手（Phase 4）。センサー設置後に再開予定。

---

## ドキュメント索引

| 種類 | ファイル | 内容 |
|------|---------|------|
| 計算 | [computational_insights_summary.md](computational_insights_summary.md) | 全 Case (A〜L) の温度結果まとめ |
| 計算 | [openfoam_computation_summary.md](openfoam_computation_summary.md) | OpenFOAM 計算ログの集約 |
| 物理 | [governing_equations.md](governing_equations.md) / `.tex` / `.pdf` | 支配方程式（簡易版 / OpenFOAM 両系統） |
| 物理 | [technical_reference.md](technical_reference.md) | 物理モデル / 制約 / 室内ダイアグラム |
| 運用 | [openfoam_troubleshooting.md](openfoam_troubleshooting.md) | OpenFOAM エラー記録と修正チェックリスト |
| 運用 | [openfoam_current_issues.md](openfoam_current_issues.md) | 未解決事項（古い・要更新） |
| 状況 | [project_status_report.md](project_status_report.md) | 2026-04-05 時点の総合状況（やや古い） |
| パラメトリック | parametric_study_*.md | 各パラスタディの個別レポート |
| 計測 | measurement_*.md, procurement_*.md | センサー設置と調達計画 |
| 振り返り | [review_response_2026_04_06.md](review_response_2026_04_06.md) | 外部レビューへの応答メモ |

## 重要スクリプト

- `scripts/run_openfoam_wsl.sh` — WSL2 経由 OpenFOAM 実行のエントリ
- `scripts/run_and_plot.py` — 簡易版 3 シナリオ比較プロット
- `scripts/plot_openfoam_results.py` — OpenFOAM vs 簡易版比較
- `scripts/run_parametric_*.sh` — 各パラスタディのバッチ
- `scripts/generate_html_report.py`, `generate_html_index.py` — 結果 HTML 化
- `tools/loyly_calculator.html` — ロウリュ蒸気量の電卓 UI
