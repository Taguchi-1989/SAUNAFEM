# 皮膚への局所熱流束モデル設計案

## 1. 背景と課題

### 1.1 現状

`simple_solver._perceived_temperature` で皮膚表面の熱収支を計算しているが、
対流熱伝達率 (HTC) が固定値:

```python
# src/harness/simple_solver.py:248
H_CONV = 8.0     # convective HTC [W/(m2*K)]
...
q_conv = H_CONV * (t_c - T_SKIN)
```

これは座位・無風・自然対流環境を想定した値であり、以下の局所現象が
表現できない:

- アウフグース時の顔面風速ピーク（K-05 の対象現象）
- 人体プルーム（自身の体温で立ち上る上昇流）
- ヒーター対面 vs 背面で受ける熱流束の差
- 上段/下段ベンチの対流熱伝達率の違い

### 1.2 求めたい量

サウナ室内の人体上の点 $\mathbf{x}$ における対流熱流束:

$$
q_{\text{conv}}(\mathbf{x}) = u(\mathbf{x}) \cdot A(\mathbf{x}) \cdot \Delta T(\mathbf{x})
$$

- $q$: 対流熱流束 [W]（あるいは熱流束密度 [W/m²]）
- $u$: 局所対流熱伝達率 [W/(m²·K)] — **本提案の主題**
- $A$: 微小表面積 [m²]
- $\Delta T = T_{\text{air,local}} - T_{\text{skin}}$ [K]

---

## 2. $u$ の物理

サウナでは 3 種類の対流が共存する。$u$ はそれぞれで桁が違う。

| 領域 | 駆動力 | スケール | $u$ の典型値 |
|------|--------|----------|--------------|
| 自然対流（静座時） | 体表-空気の温度差で生じる人体プルーム | $V \approx 0.1$ m/s | 3–6 W/(m²·K) |
| 強制対流（アウフグース） | タオル扇動による誘起流 | $V \approx 1$–3 m/s | 15–40 W/(m²·K) |
| 混合対流（中間） | 両者が同オーダー | $V \approx 0.3$ m/s | 8–15 W/(m²·K) |

### 2.1 自然対流側 (Mitchell 1974, ASHRAE 55)

人体の鉛直部位（胴・腕）では:

$$
u_{\text{nat}} = 2.38 \, |T_{\text{skin}} - T_{\text{air}}|^{0.25} \quad [\text{W/(m²·K)}]
$$

水平上向き面（肩・頭頂部）はやや増し、下向き面（足裏・座面接触部）は減る。

### 2.2 強制対流側 (de Dear et al. 1997, Watanabe et al. 1993)

裸体の人体部位に対する実測フィット:

| 部位 | 相関式 $u_{\text{forc}}(V)$ [W/(m²·K)] | 適用範囲 V [m/s] |
|------|---------------------------------------|------------------|
| 顔面 | $14.0 \, V^{0.61}$ | 0.1–4 |
| 胸部 | $10.4 \, V^{0.56}$ | 0.1–4 |
| 背中 | $9.1  \, V^{0.61}$ | 0.1–4 |
| 上腕 | $13.4 \, V^{0.55}$ | 0.1–4 |
| 全身平均 | $10.4 \, V^{0.6}$ | 0.1–4 |

これらは標準的な ISO 9920/ASHRAE 55 系の人体熱伝達係数ライブラリ。

### 2.3 自然+強制の合成 (混合対流)

ASHRAE が採用する Churchill 型ブレンド:

$$
u = \left( u_{\text{nat}}^{n} + u_{\text{forc}}^{n} \right)^{1/n}, \quad n \approx 3
$$

$u_{\text{nat}}$ が支配的な低風速域から、$u_{\text{forc}}$ が支配的な高風速域へ
滑らかに遷移する。実装上は単純な `max(u_nat, u_forc)` でも 10% 以内の誤差。

---

## 3. $u$ にどの $V$ を使うか — 4 案

これがこの提案の核心。**$u$ そのものより「どの速度を入れるか」が難しい。**

### 案 A: 室内代表速度（最も簡単・低精度）

`SimpleSolverResult` に「アウフグース時の代表面風速 $V_{\text{face}}$」を持たせ、
全身一律でその $V$ を使う。

- $V_{\text{face}}$ は K-05 (`compute_k05`) で既に計算済 (`beta_aug / (rho * a_face)`)
- 自然対流時は $V = 0.1$ m/s 固定
- 実装コスト: 数行
- 限界: 部位ごとの差・ヒーター対面/背面の差が出ない

### 案 B: 部位セグメント分解（中精度）

人体を $N$ セグメント（顔・胸・背・腕・脚 ×2、計 8〜16 部位）に分け、
各部位ごとに:

- 露出面積 $A_i$
- 表面温度 $T_{\text{skin},i}$（ISO 9920 の部位別値）
- 局所気流速度 $V_i$
- 局所空気温度 $T_{\text{air},i}$（上層/下層に応じて $T_{\text{upper}}$ or $T_{\text{lower}}$）

部位別 $V_i$ の与え方:
- アウフグース時、顔面のみ $V_{\text{aug}}$、他は $0.5 V_{\text{aug}}$ 程度
- 静座時、頭頂・肩は人体プルーム速度 $V_p \approx 0.2$ m/s、下半身は $0.1$ m/s

総熱流束:
$$
Q_{\text{conv,total}} = \sum_i u_i \, A_i \, (T_{\text{air},i} - T_{\text{skin},i})
$$

- 実装コスト: 設定 YAML に `body_segments` セクション追加 + 集計コード
- 限界: $V_i$ は経験則ベース（CFD なし）

### 案 C: CFD プローブ駆動（高精度・OpenFOAM 必要）

OpenFOAM ケースに「人体プローブ点」を配置し、解析後に各点の $|U|$ と $T$ を読み取り、
案 B と同じ集計を後処理で行う。

- プローブ位置: 顔（座面 +1.0 m）、胸（+0.7 m）、背（+0.7 m, 壁側）、膝（+0.3 m）など
- 既存 `probe_parser.py` を流用可能
- 実装コスト: テンプレートに人体プローブ追加 + `kpi.py` に集計関数
- 限界: 人体ジオメトリ自体は無いので、人体プルームは捉えられない（ヒーター由来流のみ）

### 案 D: 共役熱伝達 / マニキン CFD（将来検討）

人体形状を OpenFOAM メッシュに含め、皮膚を熱境界として共役解析する。
- コスト極大、PoC スコープ外。`要件定義書.md` の Phase 6 以降で検討。

---

## 4. 推奨ステップ

| Step | 内容 | 目的 | 実装コスト |
|------|------|------|------------|
| **S1** | 案 A を `_perceived_temperature` に適用。`H_CONV` を $u_{\text{nat}}$ + $u_{\text{forc}}(V_{\text{face}})$ の合成に置換 | アウフグース時の体感温度を風速応答に | 0.5d |
| **S2** | 部位別係数（案 B）を導入。設定は YAML の `body_segments`。デフォルト値は ISO 9920 から | K-06 を「平均体感」から「部位別熱流束分布」に拡張 | 1–2d |
| **S3** | OpenFOAM ケースに人体プローブを追加（案 C）。`probe_parser` で $V$, $T$ を読んで $u_i$ を後処理計算 | CFD と整合した皮膚熱流束 | 2–3d |

S1 だけで K-05 と K-06 が物理的に連動する。S2/S3 は精度向上として段階的に。

---

## 5. 推奨実装案 (S1 詳細)

### 5.1 関数シグネチャ変更

```python
# 現状
def _perceived_temperature(t_c: float, rh: float, q_rad_body: float = 0.0) -> float:
    H_CONV = 8.0
    q_conv = H_CONV * (t_c - T_SKIN)
    ...

# 変更後
def _convective_htc(
    t_air_c: float,
    t_skin_c: float = 36.0,
    v_local: float = 0.1,
    body_part: str = "average",  # "face" | "chest" | "back" | "average"
) -> float:
    """Local convective heat transfer coefficient [W/(m²·K)]."""
    # Natural convection (Mitchell 1974)
    u_nat = 2.38 * abs(t_skin_c - t_air_c) ** 0.25
    # Forced convection (de Dear 1997)
    coeffs = {
        "face":    (14.0, 0.61),
        "chest":   (10.4, 0.56),
        "back":    ( 9.1, 0.61),
        "arm":     (13.4, 0.55),
        "average": (10.4, 0.60),
    }
    a, b = coeffs.get(body_part, coeffs["average"])
    u_forc = a * max(v_local, 0.05) ** b
    # Churchill blend (n=3)
    return (u_nat ** 3 + u_forc ** 3) ** (1.0 / 3.0)


def _perceived_temperature(
    t_c: float,
    rh: float,
    q_rad_body: float = 0.0,
    v_local: float = 0.1,           # NEW
    body_part: str = "average",     # NEW
) -> float:
    h_conv = _convective_htc(t_c, T_SKIN, v_local, body_part)
    q_conv = h_conv * (t_c - T_SKIN)
    # Lewis 関係の蒸発項も h_conv 連動に変更:
    q_evap = 16.5 * h_conv * (...)
    ...
```

### 5.2 呼び出し側 ($V_{\text{face}}$ の供給)

`solve_two_zone` / `solve_transient` の `_humid_air_properties` 呼び出し時、
アウフグース active 期間内なら $V_{\text{face}}$ を渡す:

```python
v_face = (
    beta_aug / (rho_upper * 0.05)  # K-05 と同じ係数
    if (beta_aug > 0 and aufguss_active)
    else 0.1                         # 自然対流
)
props = _humid_air_properties(t_upper, humidity_ratio, q_rad_body_val, v_face=v_face)
```

### 5.3 新 KPI / 出力フィールド

`SimpleSolverResult` に以下を追加:

```python
h_conv_face: float = 0.0       # 顔面対流 HTC [W/(m²·K)]
q_conv_face: float = 0.0       # 顔面対流熱流束 [W/m²]
q_total_skin: float = 0.0      # 皮膚総熱流束 (conv + rad + evap) [W/m²]
```

これらは既存の K-06 体感温度を分解した中身なので破壊的変更にはならない。

### 5.4 テスト

`tests/unit/test_skin_heat_flux.py`（新規）:

- 自然対流域 ($V = 0.1$): $u \approx 5$–6 W/(m²·K) 範囲
- 強制対流域 ($V = 2.0$): $u \approx 15$–20 W/(m²·K) 範囲
- $\Delta T = 0$ で $q_{\text{conv}} = 0$（境界条件）
- $V$ 増加で $u$ 単調増加
- 顔面 vs 背中で $u_{\text{face}} > u_{\text{back}}$ at same $V$

### 5.5 検証ケース (新規 YAML)

`configs/cases/skin_heat_flux_aufguss.yaml`:

- ベースは `aufguss_test.yaml`
- 出力: 静座時 vs アウフグース時の $q_{\text{conv,face}}$ をプロット
- 期待: $q_{\text{conv,face}}$ が 100–200 W/m² から 600–1000 W/m² に跳ね上がる

---

## 6. リスク / 留意点

| 項目 | 内容 | 対処 |
|------|------|------|
| 速度プロキシの誤差 | $V_{\text{face}} = \beta_{\text{aug}}/(\rho \cdot 0.05)$ は粗い近似 | S3 で OpenFOAM プローブ駆動に置換 |
| 蒸発項の連動 | $q_{\text{evap}}$ も Lewis 関係で $h$ 比例。$h$ が 2 倍になれば蒸発冷却も 2 倍。発汗上限 400 W/m² で頭打ちのため過大評価は限定的 | 既存の `Q_EVAP_MAX` クランプを維持 |
| 体感温度の解釈 | $u$ が動的になると `T_eq` が風速で振れる。報告値の比較条件を統一する | レポートに $V_{\text{used}}$ と $h_{\text{used}}$ を明記 |
| 人体プルームの欠落 | 案 A〜C は人体自身のプルームを陽に解かない | $V_{\text{nat}} = 0.1$ m/s をフロアとして与え、$u_{\text{nat}}$ の温度差項で実質補償 |

---

## 7. 直近のアクション

1. 本提案をレビュー → S1 のスコープ確定
2. `tests/unit/test_skin_heat_flux.py` を先に書く（テストファースト）
3. `_convective_htc` を実装、`_perceived_temperature` を移行
4. `aufguss_test.yaml` で既存 K-06 の値が妥当範囲に収まるか確認
5. `docs/governing_equations.md` の「皮膚熱収支モデル」節に式追記
