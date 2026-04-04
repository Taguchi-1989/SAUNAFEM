# SaunaFlow 支配方程式ドキュメント

本ドキュメントでは、SaunaFlow で使用する **2 つのソルバー** それぞれの支配方程式、
離散化手法、ゾーン構成、数値解法を整理する。

---

## 目次

1. [ソルバー概要](#1-ソルバー概要)
2. [正確バージョン — OpenFOAM (buoyantPimpleFoam)](#2-正確バージョン--openfoam-buoyantpimplefoam)
3. [簡易バージョン — 2-Zone プルームモデル (simple_solver.py)](#3-簡易バージョン--2-zone-プルームモデル)
4. [ゾーン構成と熱収支](#4-ゾーン構成と熱収支)
5. [離散化スキームの対応表](#5-離散化スキームの対応表)
6. [数値解法と収束制御](#6-数値解法と収束制御)
7. [物性値](#7-物性値)
8. [簡易版と正確版の対応関係](#8-簡易版と正確版の対応関係)

---

## 1. ソルバー概要

| 項目 | 正確版 (OpenFOAM) | 簡易版 (simple_solver.py) |
|------|-------------------|--------------------------|
| 次元 | 3D | 0D (2ゾーン集中定数) |
| 方程式 | N-S + エネルギー + 乱流 | ゾーン質量保存 + エネルギー保存 + プルーム相関式 |
| 乱流モデル | SST k-omega | なし (プルーム相関式に内包) |
| 圧力-速度連成 | PIMPLE 法 (非定常) | なし (速度場を解かない) |
| 時間積分 | backward (2次精度陰解法) | 前進 Euler (擬似時間進行) |
| 用途 | 本計算・検証 | 可視化 UI・パラメータスタディ |

---

## 2. 正確バージョン — OpenFOAM (buoyantPimpleFoam)

### 2.1 連続の式 (質量保存)

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{U}) = 0
$$

- $\rho$: 密度 [kg/m³]（完全ガス状態方程式で計算）
- $\mathbf{U}$: 速度ベクトル [m/s]
- 非定常計算 (buoyantPimpleFoam) のため時間微分項を保持

### 2.2 運動量保存 (Navier-Stokes)

$$
\frac{\partial (\rho \mathbf{U})}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} \otimes \mathbf{U})
= -\nabla p_{rgh}
+ \nabla \cdot \left[ \mu_{\text{eff}} \left( \nabla \mathbf{U} + (\nabla \mathbf{U})^T - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I} \right) \right]
+ \rho \mathbf{g}
$$

ここで:

- $p_{rgh} = p - \rho \mathbf{g} \cdot \mathbf{r}$: 静水圧を除いた修正圧力 [Pa]
- $\mu_{\text{eff}} = \mu + \mu_t$: 実効粘性 [Pa·s]
- $\mathbf{g} = (0, -9.81, 0)$ m/s²: 重力加速度
- **浮力駆動流**: 密度差 → 圧力勾配 → プルーム上昇

離散化 (`fvSchemes`):

```
div(phi,U):  bounded Gauss linearUpwind grad(U)   ← 2次風上 (数値拡散低減)
grad(U):     cellLimited Gauss linear 1            ← セル制限付き線形
ddt:         backward                              ← 2次精度陰的時間積分
```

### 2.3 エネルギー保存

$$
\frac{\partial (\rho h)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} h) = \nabla \cdot (\alpha_{\text{eff}} \nabla T) + q_{\text{source}}
$$

- $h = c_p T$: 比エンタルピー (sensibleEnthalpy) [J/kg]
- $\alpha_{\text{eff}} = \frac{\mu}{Pr} + \frac{\mu_t}{Pr_t}$: 実効熱拡散率 [W/(m·K)]
- $q_{\text{source}}$: ヒーター熱流束 [W/m²]（境界条件として付与）
- $c_p = 1005$ J/(kg·K), $Pr = 0.7$

離散化:

```
div(phi,T):  bounded Gauss linearUpwind default    ← 2次風上
laplacian:   Gauss linear corrected                ← 2次中心差分 + 非直交補正
```

### 2.4 乱流モデル (SST k-omega)

Menter の SST (Shear Stress Transport) k-omega モデルを採用。
標準 k-epsilon に比べ、壁面近傍の低レイノルズ数領域と
自然対流における熱伝達率の予測精度が優れる。

**乱流運動エネルギー k:**

$$
\frac{\partial (\rho k)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} k)
= \nabla \cdot \left[ (\mu + \sigma_k \mu_t) \nabla k \right]
+ P_k - \beta^* \rho k \omega
$$

**比散逸率 omega:**

$$
\frac{\partial (\rho \omega)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} \omega)
= \nabla \cdot \left[ (\mu + \sigma_\omega \mu_t) \nabla \omega \right]
+ \frac{\gamma}{\nu_t} P_k - \beta \rho \omega^2
+ 2(1 - F_1) \frac{\rho \sigma_{\omega 2}}{\omega} \nabla k \cdot \nabla \omega
$$

**乱流粘性:**

$$
\mu_t = \frac{\rho a_1 k}{\max(a_1 \omega, \; S F_2)}
$$

ここで $F_1$, $F_2$ はブレンディング関数、$S$ はひずみ速度テンソルの大きさ。

**SST k-omega の利点 (サウナ解析において):**

- 壁面近傍: k-omega モデルが自動的に適用され、壁面関数への依存が軽減
- 主流域: k-epsilon 相当の挙動にブレンドされ安定
- 自然対流の壁面熱伝達率を k-epsilon より正確に予測

**壁面境界条件:**

| 変数 | 壁面 BC |
|------|---------|
| k | kqRWallFunction |
| omega | omegaWallFunction |
| nut | nutkWallFunction |
| alphat | compressible::alphatWallFunction (Prt = 0.85) |

離散化:

```
div(phi,k):      bounded Gauss upwind    ← 1次風上 (安定性重視)
div(phi,omega):  bounded Gauss upwind    ← 1次風上
```

### 2.5 状態方程式

$$
\rho = \frac{p M_w}{R T}
$$

- $M_w = 28.96$ g/mol (空気の分子量)
- $R = 8314$ J/(kmol·K)
- 完全ガス (perfectGas) 仮定
- 熱物性モデル: `heRhoThermo` + `pureMixture` + `hConst` + `const` transport

**浮力の発生メカニズム:**
ヒーター近傍で温度上昇 → 密度低下 → 浮力 → 上昇プルーム形成

### 2.6 時間平均場

非定常計算のため、統計量として時間平均場を出力する。
`fieldAverage` 関数オブジェクトにより、指定開始時刻以降の
$\bar{T}$, $\bar{U}$, $\overline{T'^2}$, $\overline{U'^2}$ を計算。

```
averaging_start: end_time * 0.5  (デフォルト: 後半50%を平均化)
```

---

## 3. 簡易バージョン — 2-Zone プルームモデル

### 3.1 モデル概要

建築環境工学で標準的な **2層ゾーンモデル** を採用。
Morton-Taylor-Turner (MTT) のエントレインメント理論と
Zukoski のプルーム相関式に基づく。

```
 ┌─────────────────────────────────┐  ← ceiling
 │                                 │
 │   上層 (Upper hot layer)         │  温度: T_upper
 │   密度: rho_upper = rho_0 * T_wall / T_upper
 │                                 │
 ├ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┤  ← 界面高さ z_int
 │                                 │
 │   下層 (Lower cool layer)        │  温度: T_lower
 │   密度: rho_lower = rho_0 * T_wall / T_lower
 │          ┌───┐                  │
 │          │ H │ ← ヒーター         │
 │          │   │   (プルーム発生源)   │
 └──────────┴───┴──────────────────┘  ← floor
              ↑
           プルーム上昇
         (エントレインメント)
```

**1D 拡散モデルとの根本的な違い:**

- 対流を「大きな熱拡散」で近似する代わりに、プルームの質量・エネルギー輸送を直接モデル化
- 質量保存則を遵守（プルーム上昇量 = 壁面沿い下降量）
- 物理的根拠のある相関式（Zukoski のプルーム相関）を使用

### 3.2 プルームモデル (Zukoski 相関式)

ヒーター中心から高さ $z$ の位置におけるプルーム質量流量:

$$
\dot{m}_p(z) = 0.071 \, Q_c^{1/3} \, z^{5/3} + 0.0018 \, Q_c
$$

ここで:

- $Q_c$: 対流熱出力 [kW] ($Q_c = f_{\text{conv}} \times Q_{\text{heater}}$, $f_{\text{conv}} = 0.7$)
- $z$: ヒーター中心からの鉛直距離 [m]
- $\dot{m}_p$: プルーム質量流量 [kg/s]

プルーム温度（エネルギー保存から）:

$$
T_{\text{plume}} = T_{\text{lower}} + \frac{Q_c}{\dot{m}_p \, c_p}
$$

**物理的意味:** プルームが上昇するにつれ周囲空気をエントレイン（巻き込み）し、
質量流量は増加するが温度は低下する。これは MTT 理論の本質的な特徴。

### 3.3 支配方程式

**状態変数:** $T_{\text{upper}}$ (上層温度), $T_{\text{lower}}$ (下層温度), $z_{\text{int}}$ (界面高さ)

#### 上層エネルギー保存

$$
\rho_{\text{upper}} \, c_p \, V_{\text{upper}} \, \frac{dT_{\text{upper}}}{dt}
= \underbrace{\dot{m}_p \, c_p \, (T_{\text{plume}} - T_{\text{upper}})}_{\text{プルームからの熱入力}}
- \underbrace{h_{\text{wall}} \, A_{\text{wall,upper}} \, (T_{\text{upper}} - T_{\text{wall}})}_{\text{壁面熱損失}}
$$

ここで:

- $V_{\text{upper}} = A_{\text{floor}} \times (H - z_{\text{int}})$: 上層体積 [m³]
- $A_{\text{wall,upper}} = P \times (H - z_{\text{int}}) + A_{\text{floor}}$: 上層が接する壁面積 (側壁 + 天井) [m²]
- $P = 2(W + D)$: 水平断面の周長 [m]

#### 下層エネルギー保存

$$
\rho_{\text{lower}} \, c_p \, V_{\text{lower}} \, \frac{dT_{\text{lower}}}{dt}
= \underbrace{Q_{\text{rad}} \times 0.3}_{\text{輻射受熱}}
+ \underbrace{k_{\text{int}} \, A_{\text{floor}} \, \frac{T_{\text{upper}} - T_{\text{lower}}}{0.1 H}}_{\text{界面伝導}}
- \underbrace{h_{\text{wall}} \, A_{\text{wall,lower}} \, (T_{\text{lower}} - T_{\text{wall}})}_{\text{壁面熱損失}}
$$

ここで:

- $Q_{\text{rad}} = (1 - f_{\text{conv}}) \times Q_{\text{heater}}$: ヒーターからの輻射熱 [W]
- $k_{\text{int}} = 0.5$ W/(m·K): 界面の実効熱伝導率

#### 界面質量保存

$$
A_{\text{floor}} \, \frac{dz_{\text{int}}}{dt}
= -\frac{\dot{m}_p}{\rho_{\text{upper}}} + \dot{V}_{\text{return}}
$$

ここで:

- $\dot{m}_p / \rho_{\text{upper}}$: プルームが上層に供給する体積流量（界面を押し下げる）
- $\dot{V}_{\text{return}}$: 壁面沿い下降流による体積流量（界面を押し上げる）

$$
\dot{V}_{\text{return}} = \frac{h_{\text{wall}} \, A_{\text{wall,upper}} \, (T_{\text{upper}} - T_{\text{wall}})}{\rho_{\text{upper}} \, c_p \, \max(T_{\text{upper}} - T_{\text{lower}}, \, 1)}
$$

定常状態ではこの 2 つがバランスし、界面高さが安定する。

#### 密度の温度依存性

$$
\rho = \rho_0 \frac{T_{\text{wall}}}{T}
$$

- $\rho_0 = 1.1$ kg/m³ (基準密度, ~300K)

### 3.4 温度プロファイルの構築

2-Zone の集中定数 ($T_{\text{upper}}$, $T_{\text{lower}}$) から
鉛直方向の連続プロファイルをシグモイド関数で補間:

$$
T(y) = T_{\text{lower}} + \frac{T_{\text{upper}} - T_{\text{lower}}}{1 + \exp\left(-3 \cdot \frac{y - z_{\text{int}}}{\delta / 2}\right)}
$$

ここで $\delta = 0.15 H$ は遷移層の厚さ。

### 3.5 境界条件と制約

| 境界/制約 | 条件 |
|-----------|------|
| 界面高さ | $0.05 H \le z_{\text{int}} \le 0.95 H$ |
| 上層温度 | $T_{\text{wall}} \le T_{\text{upper}} \le T_{\text{wall}} + 200$ K |
| 下層温度 | $T_{\text{wall}} - 1 \le T_{\text{lower}} \le T_{\text{upper}}$ |

### 3.6 参考文献

- Morton, Taylor & Turner (1956), "Turbulent Gravitational Convection from Maintained and Instantaneous Sources", Proc. R. Soc. A 234:1-23
- Zukoski (1978), "Development of a Stratified Ceiling Layer in the Early Stages of a Closed-Room Fire", NBS-GCR-78-150
- Cooper (1982), "A Mathematical Model for Estimating Available Safe Egress Time in Fires", NBSIR 82-2612

---

## 4. ゾーン構成と熱収支

### 4.1 正確版 (OpenFOAM) のゾーン構成

```
              fixedValue T_wall
         ┌──────────────────────────┐
         │        ceiling            │
         │                          │
         │     ┌──────┐             │
  heater │     │ 内部 │             │ fixedValue
  _wall  │     │ 領域 │             │ T_wall
 (熱流束)│     │      │             │ (opposite_wall)
         │     │      │             │
         │     └──────┘             │
         │  heater_wall_surround    │
         │   (fixedValue T_wall)    │
         └──────────────────────────┘
                  floor
              fixedValue T_wall

  front/back: fixedValue T_wall
```

**各パッチの境界条件:**

| パッチ名 | 温度 BC | 速度 BC | 圧力 BC |
|----------|---------|---------|---------|
| heater_wall | externalWallHeatFluxTemperature (flux) | noSlip | fixedFluxPressure |
| heater_wall_surround | fixedValue $T_{\text{wall}}$ | noSlip | fixedFluxPressure |
| floor | fixedValue $T_{\text{wall}}$ | noSlip | fixedFluxPressure |
| ceiling | fixedValue $T_{\text{wall}}$ | noSlip | fixedFluxPressure |
| opposite_wall | fixedValue $T_{\text{wall}}$ | noSlip | fixedFluxPressure |
| front, back | fixedValue $T_{\text{wall}}$ | noSlip | fixedFluxPressure |

**ヒーター熱流束の計算:**

$$
q_{\text{heater}} = \frac{Q_{\text{kW}} \times 1000}{W_{\text{heater}} \times H_{\text{heater}}} \quad [\text{W/m}^2]
$$

### 4.2 セル単位の熱収支 (有限体積法)

OpenFOAM は有限体積法 (FVM) で各セルのエネルギー収支を解く:

$$
\underbrace{\frac{\partial}{\partial t} \int_{V_P} \rho h \, dV}_{\text{蓄熱}}
+ \underbrace{\sum_f (\rho \mathbf{U} h)_f \cdot \mathbf{S}_f}_{\text{対流フラックス}}
= \underbrace{\sum_f (\alpha_{\text{eff}} \nabla T)_f \cdot \mathbf{S}_f}_{\text{拡散フラックス}}
+ \underbrace{q \cdot V_P}_{\text{体積熱源}}
$$

- $f$: セル面, $\mathbf{S}_f$: 面積ベクトル
- $V_P$: セル体積
- **蓄熱**: 非定常項（backward スキームで 2 次精度離散化）
- **対流フラックス**: 流れによる熱輸送（linearUpwind で離散化）
- **拡散フラックス**: 伝導+乱流拡散（Gauss linear corrected で離散化）

### 4.3 簡易版の熱収支模式図

```
     ┌─────────────────────────────────┐
     │         上層 (T_upper)           │
     │                                 │
     │  入力: plume enthalpy           │
     │  出力: 壁面損失 (天井 + 上部側壁) │
     ├ ─ ─ ─ ─ ─ ─ z_int ─ ─ ─ ─ ─ ─ ┤  ← 界面 (質量保存で決定)
     │         下層 (T_lower)           │
     │                                 │
     │  入力: 輻射 + 界面伝導           │
     │  出力: 壁面損失 (床 + 下部側壁)   │
     │  出力: プルームへの空気供給       │
     │          ↑ プルーム              │
     │         [H] ヒーター              │
     └─────────────────────────────────┘
```

---

## 5. 離散化スキームの対応表

| 微分項 | 物理的意味 | OpenFOAM スキーム | 簡易版での扱い |
|--------|-----------|-------------------|---------------|
| $\partial/\partial t$ | 時間変化 | backward (2次精度陰解法) | 前進 Euler ($\Delta t = 0.5$) |
| $\nabla \cdot (\rho \mathbf{U} \phi)$ | 対流 | linearUpwind (2次風上) | プルーム相関式で代替 |
| $\nabla \cdot (k \nabla T)$ | 拡散 | Gauss linear corrected | 界面伝導 ($k_{\text{int}}$) |
| $\nabla p$ | 圧力勾配 | Gauss linear | 解かない |
| $q_{\text{source}}$ | 熱源 | 境界条件 (面熱流束) | プルームエンタルピー |
| 壁面熱伝達 | 壁損失 | 壁面関数 + fixedValue | Newton 冷却則 ($h_{\text{wall}}$) |

---

## 6. 数値解法と収束制御

### 6.1 正確版 — PIMPLE 法

PIMPLE (PISO + SIMPLE の融合) アルゴリズム。
密閉空間の強い自然対流は本質的に非定常であるため、
定常ソルバー (SIMPLE) ではなく非定常ソルバーを採用し、
時間平均で統計量を取得する。

```
各時間ステップ:
  外側ループ (nOuterCorrectors = 2):
    1. 運動量方程式を解く → U* (仮の速度場)
    2. エネルギー方程式を解く → T
    3. 乱流方程式を解く → k, omega
    4. 密度を更新 → rho = f(p, T)
    内側ループ (nCorrectors = 1):
      5. 圧力補正方程式を解く → p'
      6. 速度を補正 → U = U* + U'
  時間ステップ調整 (maxCo = 0.5)
```

**線形ソルバー:**

| 変数 | ソルバー | 前処理 | 許容残差 |
|------|---------|--------|---------|
| $p_{rgh}$ | GAMG (代数マルチグリッド) | Gauss-Seidel | 1e-7 (rel 0.01) |
| $p_{rgh}$ Final | GAMG | Gauss-Seidel | 1e-7 (rel 0) |
| $U, T, k, \omega$ | PBiCGStab | DILU | 1e-7 (rel 0.1) |
| $U, T, k, \omega$ Final | PBiCGStab | DILU | 1e-7 (rel 0) |

**PIMPLE 残差制御:**

| 変数 | 閾値 |
|------|------|
| $p_{rgh}$ | 1e-4 |
| $U$ | 1e-4 |
| $T$ | 1e-5 (より厳密) |

**緩和係数 (Under-relaxation):**

| 変数 | 緩和係数 | 備考 |
|------|---------|------|
| $p_{rgh}$ | 0.3 | 保守的 (圧力は不安定になりやすい) |
| $U$ | 0.7 | 中程度 |
| $T$ | 0.5 | やや保守的 (浮力との結合が強い) |
| $k, \omega$ | 0.7 | 中程度 |

**時間ステップ制御:**

| パラメータ | 値 | 備考 |
|-----------|-----|------|
| deltaT | 0.05 s | 初期時間ステップ |
| maxCo | 0.5 | Courant 数上限 |
| adjustTimeStep | yes | 適応的時間ステップ調整 |
| endTime | 300 s | 計算終了時刻 |

### 6.2 簡易版 — 擬似時間進行

| 項目 | 値 |
|------|-----|
| 解法 | 前進 Euler（擬似時間ステップ $\Delta t = 0.5$） |
| 状態変数 | $T_{\text{upper}}$, $T_{\text{lower}}$, $z_{\text{int}}$ |
| 収束判定 | $\max(|T_{\text{upper}}^{n+1} - T_{\text{upper}}^n|, \; 10 |z_{\text{int}}^{n+1} - z_{\text{int}}^n|) < 10^{-4}$ |
| 最大反復数 | 10,000 |
| 物理クリッピング | $T_{\text{upper}} \in [T_{\text{wall}},\; T_{\text{wall}}+200]$ K |
| | $z_{\text{int}} \in [0.05H,\; 0.95H]$ |

---

## 7. 物性値

### 空気 (全ソルバー共通)

| 物性 | 記号 | 値 | 単位 |
|------|------|-----|------|
| 分子量 | $M_w$ | 28.96 | g/mol |
| 定圧比熱 | $c_p$ | 1005 | J/(kg·K) |
| 動粘性係数 | $\mu$ | $1.8 \times 10^{-5}$ | Pa·s |
| プラントル数 | $Pr$ | 0.7 | - |
| 熱伝導率 | $\lambda = \mu c_p / Pr$ | 0.0258 | W/(m·K) |

### 簡易版追加パラメータ

| パラメータ | 記号 | 値 | 根拠 |
|-----------|------|-----|------|
| 基準密度 | $\rho_0$ | 1.1 kg/m³ | ~300K での空気密度 |
| 壁面熱伝達率 | $h_{\text{wall}}$ | 8.0 W/(m²·K) | 自然対流の標準値 |
| 対流熱割合 | $f_{\text{conv}}$ | 0.7 | ヒーター出力の70%が対流 |
| 界面伝導率 | $k_{\text{int}}$ | 0.5 W/(m·K) | 界面を介した弱い伝導 |

---

## 8. 簡易版と正確版の対応関係

### 物理モデルの対応

| 物理現象 | 正確版 (OpenFOAM) | 簡易版 (2-Zone) | 簡易化による影響 |
|---------|-------------------|-----------------|----------------|
| **対流** | N-S 方程式で速度場を解く | MTT プルーム相関式で質量流量を算出 | 詳細な速度場は不明だが、プルームの総輸送量は再現 |
| **乱流** | SST k-omega で乱流粘性を計算 | プルーム相関式に乱流効果が内包 | 局所的な乱流強度の変化を再現不可 |
| **浮力** | 密度差 → 圧力勾配 → 速度場 | 密度の温度依存 + プルーム浮力相関 | 浮力-運動量の直接結合はないが、エネルギー輸送は再現 |
| **質量保存** | 連続の式で厳密に保存 | プルーム上昇 = 壁面沿い下降で保存 | ゾーン間の質量収支は保存される |
| **壁面境界層** | omegaWallFunction + メッシュ解像 | Newton 冷却 ($h_{\text{wall}}$ 固定) | 壁面近傍の詳細を無視 |
| **水平方向分布** | 3D で空間分解 | 各層内で均一を仮定 | 水平方向の温度むらを再現不可 |
| **密度変化** | 完全ガス $\rho = pM/RT$ | $\rho = \rho_0 T_{\text{wall}} / T$ | 簡易的だが温度依存あり |

### 簡易版の限界

1. **速度場の不在**: プルーム相関式は総質量流量を与えるが、空間的な速度分布は得られない。
   アウフグース（Phase 3）のジェット流れは原理的にモデル化不可能。
2. **水平方向の均一仮定**: 各層内の水平温度分布が均一のため、
   ヒーター直上と壁面近傍の温度差は再現できない。
3. **蒸気輸送の困難**: ロウリュ（Phase 2）の蒸気は対流に乗って輸送されるが、
   速度場がないため拡散モデルでの近似が必要になり精度に限界がある。

### 簡易版で正しく再現できること

- **温度成層**: 上層高温・下層低温の 2 層構造（物理的に正しいメカニズムで再現）
- **界面高さ**: プルームの強さと壁面損失のバランスで決定される界面高さ
- **KPI の傾向**: K-01 (上下温度差) > 0 の判定は信頼できる
- **パラメータ感度**: ヒーター出力・壁温・部屋サイズを変えたときの応答傾向
- **エントレインメント効果**: プルーム高さが大きいほど巻き込みが増え温度が低下する現象

---

## 付録: フェーズ別の方程式拡張予定

| フェーズ | OpenFOAM | 簡易版 | 備考 |
|---------|----------|--------|------|
| Phase 1 (現在) | 上記の非定常方程式 (buoyantPimpleFoam + SST k-omega) | 2-Zone プルームモデル | 温度成層の再現 |
| Phase 2 (lolylu) | 蒸気輸送方程式 $\frac{\partial (\rho Y)}{\partial t} + \nabla \cdot (\rho \mathbf{U} Y) = \nabla \cdot (D_{\text{eff}} \nabla Y)$ | 拡張検討中 (2-Zone への蒸気層追加は限界あり) | 蒸気は対流支配のため簡易版での再現は困難 |
| Phase 3 (Aufguss) | 運動量ソース項 (ジェット) | 簡易版では対応不可 | 速度場が必要 |
| Phase 4-5 | 方程式追加なし | - | 実験データとの比較・自動化 |
