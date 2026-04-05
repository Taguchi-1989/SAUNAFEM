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
| 方程式 | N-S + エネルギー + 乱流 + 蒸気輸送 | ゾーン質量/エネルギー保存 + プルーム + 壁面 + 湿度 |
| 状態変数 | $\mathbf{U}, p, T, k, \omega, Y$ | $T_{\text{upper}}, T_{\text{lower}}, z_{\text{int}}, T_{\text{wall,inner}}, w$ |
| 乱流モデル | SST k-omega | なし (プルーム相関式に内包) |
| 壁面モデル | 壁面関数 + fixedValue | 集中定数 (lumped) または固定温度 |
| 輻射モデル | なし (Phase 2 で fvDOM 導入予定) | View Factor (幾何近似) |
| 湿度連成 | multiComponentMixture (H₂O) | 混合物性 ($c_{p,\text{mix}}, \lambda_{\text{mix}}, h_{\text{wall,eff}}$) |
| 圧力-速度連成 | PIMPLE 法 (非定常) | なし (速度場を解かない) |
| 時間積分 | backward (2次精度陰解法) | 定常: 擬似時間 / 過渡: 前進 Euler (実時間) |
| 用途 | 本計算・検証 | 可視化 UI・パラメータスタディ・KPI 算出 |

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

### 3.3 支配方程式（現行実装）

**状態変数:** $T_{\text{upper}}$, $T_{\text{lower}}$, $z_{\text{int}}$, $T_{\text{wall,inner}}$ (lumped壁面), $w$ (湿度比)

#### 湿度連成物性

全ての熱物性は湿度比 $w$ [kg vapor / kg dry air] に依存する:

$$
y_v = \frac{w}{1+w}, \quad
c_{p,\text{mix}} = (1-y_v) c_{p,\text{air}} + y_v c_{p,\text{vapor}}, \quad
\lambda_{\text{mix}} = (1-y_v) \lambda_{\text{air}} + y_v \lambda_{\text{vapor}}
$$

$$
h_{\text{wall,eff}} = h_{\text{wall,base}} \left(\frac{c_{p,\text{mix}}}{c_{p,\text{air}}}\right)^{0.25} \left(\frac{\lambda_{\text{mix}}}{\lambda_{\text{air}}}\right)^{0.75}
$$

以下では $c_p = c_{p,\text{mix}}$, $h_{\text{wall}} = h_{\text{wall,eff}}$ と略記する。

#### 上層エネルギー保存

$$
\rho_{\text{upper}} \, c_p \, V_{\text{upper}} \, \frac{dT_{\text{upper}}}{dt}
= \underbrace{\dot{m}_p \, c_p \, (T_{\text{plume}} - T_{\text{upper}})}_{\text{プルームからの熱入力}}
- \underbrace{h_{\text{wall}} \, A_{\text{wall,upper}} \, (T_{\text{upper}} - T_{\text{wall,inner}})}_{\text{壁面への熱損失}}
$$

ここで:

- $V_{\text{upper}} = A_{\text{floor}} \times (H - z_{\text{int}})$: 上層体積 [m³]
- $A_{\text{wall,upper}} = P \times (H - z_{\text{int}}) + A_{\text{floor}}$: 上層が接する壁面積 (側壁 + 天井) [m²]
- $P = 2(W + D)$: 水平断面の周長 [m]
- $T_{\text{wall,inner}}$: 壁面内面温度（lumped壁面モデルで計算。`model: fixed` の場合は $T_{\text{wall,outer}}$ 固定）

**注:** ロウリュ（蒸気）は上層の湿度 $w$ を変化させるが、
蒸発潜熱はヒーター石から奪われるため、上層空気に正味のエネルギー追加はない。

#### 下層エネルギー保存

$$
\rho_{\text{lower}} \, c_p \, V_{\text{lower}} \, \frac{dT_{\text{lower}}}{dt}
= \underbrace{Q_{\text{rad,lower}}}_{\text{輻射受熱}}
+ \underbrace{k_{\text{int}} \, A_{\text{floor}} \, \frac{T_{\text{upper}} - T_{\text{lower}}}{0.1 H}}_{\text{界面伝導}}
- \underbrace{h_{\text{wall}} \, A_{\text{wall,lower}} \, (T_{\text{lower}} - T_{\text{wall,inner}})}_{\text{壁面への熱損失}}
$$

**輻射経路（壁面モデルに依存）:**

- **`model: lumped`（推奨）:** ヒーター輻射は全て壁面に吸収される。
  壁面が温まり、壁面から空気への対流伝熱で間接的に加熱。
  → $Q_{\text{rad,lower}} = 0$（壁面モデルが輻射を処理する）
- **`model: fixed`:** 壁温固定のため壁面の蓄熱なし。
  ヒーター輻射の一部が直接下層空気に入力される。
  → $Q_{\text{rad,lower}} = Q_{\text{rad}} \times f_{\text{rad,lower}}$

ここで:

- $Q_{\text{rad}} = (1 - f_{\text{conv}}) \times Q_{\text{heater}}$: ヒーターからの輻射熱 [W]
- $f_{\text{rad,lower}} = F_{\text{heater→floor}} + F_{\text{heater→lower walls}}$: 幾何学的 View Factor
- $k_{\text{int}} = 0.5$ W/(m·K): 界面の実効熱伝導率

#### 壁面集中定数モデル（`model: lumped`）

壁面内面温度の時間発展:

$$
(\rho c_p)_w \, \delta_w \, A_{\text{wall,total}} \, \frac{dT_{\text{wall,inner}}}{dt}
= \underbrace{Q_{\text{conv→wall}}}_{\text{空気→壁面}}
+ \underbrace{Q_{\text{rad}}}_{\text{ヒーター輻射全量}}
- \underbrace{\frac{\lambda_w}{\delta_w} A_{\text{wall,total}} (T_{\text{wall,inner}} - T_{\text{wall,outer}})}_{\text{壁面→外部}}
$$

ここで:

- $Q_{\text{conv→wall}} = h_{\text{wall}} A_{\text{wall,upper}} (T_{\text{upper}} - T_{\text{wall,inner}}) + h_{\text{wall}} A_{\text{wall,lower}} (T_{\text{lower}} - T_{\text{wall,inner}})$
- $Q_{\text{rad}} = (1 - f_{\text{conv}}) Q_{\text{heater}}$: ヒーター輻射全量が壁面に吸収
- $\delta_w = 0.015$ m, $\lambda_w = 0.12$ W/(m·K), $(\rho c_p)_w = 5 \times 10^5$ J/(m³·K)

#### 界面質量保存

$$
A_{\text{floor}} \, \frac{dz_{\text{int}}}{dt}
= -\frac{\dot{m}_p}{\rho_{\text{upper}}} + \dot{V}_{\text{return}} - \dot{V}_{\text{steam}}
$$

ここで:

- $\dot{m}_p / \rho_{\text{upper}}$: プルームが上層に供給する体積流量（界面を押し下げる）
- $\dot{V}_{\text{return}}$: 壁面沿い下降流による体積流量（界面を押し上げる）
- $\dot{V}_{\text{steam}}$: 蒸気体積膨張（ロウリュ時のみ、過渡ソルバーで使用）

$$
\dot{V}_{\text{return}} = \frac{h_{\text{wall}} \, A_{\text{wall,upper}} \, (T_{\text{upper}} - T_{\text{wall,inner}})}{\rho_{\text{upper}} \, c_p \, \max(T_{\text{upper}} - T_{\text{lower}}, \, 1)}
$$

定常状態ではこれらがバランスし、界面高さが安定する。

#### アウフグース強制混合（`aufguss` 設定時のみ）

定常ソルバーでは常時有効、過渡ソルバーでは指定時間窓で有効:

$$
\frac{dT_{\text{upper}}}{dt} \mathrel{-}= \frac{\beta_{\text{aug}} \, c_p \, (T_{\text{upper}} - T_{\text{lower}})}{m_{\text{upper}} \, c_p}
$$

$$
\frac{dT_{\text{lower}}}{dt} \mathrel{+}= \frac{\beta_{\text{aug}} \, c_p \, (T_{\text{upper}} - T_{\text{lower}})}{m_{\text{lower}} \, c_p}
$$

$\beta_{\text{aug}}$ [kg/s] は OpenFOAM の事前計算から抽出する ROM パラメータ。

#### 蒸気投入（ロウリュ）

**定常ソルバー:** 全水量が蒸発した後の平衡状態を計算する。
蒸発潜熱は石/ヒーター側から奪われるため、空気への正味エネルギー追加はない。
湿度を一括設定: $w = m_{\text{water}} / m_{\text{upper}}$

**過渡ソルバー:** Spalding 型の指数減衰蒸発モデルを区間積分で評価:

$$
\Delta m_{\text{steam}} = m_{\text{water}} \left( e^{-t_0/\tau} - e^{-t_1/\tau} \right)
$$

蒸発した蒸気は湿度 $w$ を増加させ、体積膨張 $\dot{V}_{\text{steam}}$ として界面を押し下げる。

#### 密度の温度依存性

$$
\rho = \rho_0 \frac{T_{\text{wall,outer}}}{T}
$$

- $\rho_0 = 1.1$ kg/m³ (基準密度, ~300K)

#### 収束判定（定常ソルバー）

$$
\max\left( \lvert \Delta T_{\text{upper}} \rvert, \; \lvert \Delta T_{\text{lower}} \rvert, \; 10 \lvert \Delta z_{\text{int}} \rvert, \; \lvert \Delta T_{\text{wall,inner}} \rvert \right) < 10^{-4}
$$

4つの状態変数すべてが収束するまで反復する（最低100反復）。

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
| 収束判定 | $\max(\lvert T_{\text{upper}}^{n+1} - T_{\text{upper}}^n \rvert, \; 10 \lvert z_{\text{int}}^{n+1} - z_{\text{int}}^n \rvert) < 10^{-4}$ |
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
| 壁面熱伝達率（基準） | $h_{\text{wall,base}}$ | 8.0 W/(m²·K) | 自然対流の標準値（湿度で補正） |
| 対流熱割合 | $f_{\text{conv}}$ | 0.7 | ヒーター出力の70%が対流 |
| 界面伝導率 | $k_{\text{int}}$ | 0.5 W/(m·K) | 界面を介した弱い伝導 |
| 水蒸気比熱 | $c_{p,\text{vapor}}$ | 1860 J/(kg·K) | NIST (373K) |
| 水蒸気熱伝導率 | $\lambda_{\text{vapor}}$ | 0.025 W/(m·K) | Incropera (373K) |
| 蒸発潜熱 | $L_v$ | $2.26 \times 10^6$ J/kg | 100°C |

### 壁面モデルパラメータ（`model: lumped` 使用時）

| パラメータ | 記号 | 値 | 根拠 |
| --------- | ---- | --- | ---- |
| 木材パネル厚さ | $\delta_w$ | 0.015 m | サウナ内装板の標準厚 |
| 木材熱伝導率 | $\lambda_w$ | 0.12 W/(m·K) | 針葉樹（杉/ヒノキ） |
| 木材体積熱容量 | $(\rho c_p)_w$ | $5 \times 10^5$ J/(m³·K) | $\rho \approx 450$ kg/m³, $c_p \approx 1100$ J/(kg·K) |

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

## 9. 既知の物理的制約とフェーズ拡張計画

本節では、Phase 1 の方程式系に内在する物理的制約を明示し、
Phase 2–3 で必要となる方程式拡張の具体的方針を示す。

### 9.1 蒸気化による体積膨張の無視 (Phase 2 で致命的)

**問題:** 現在の界面質量保存式は、乾燥空��のプルームのみを扱う。
ロウリュ時に水 100 mL が気化すると、理想気体近似で約 170 L の蒸気が瞬時に生成される。
一般的なサウナ室（体積 $V \sim 20$ m³）においてこの体積は無視できず、
界面高さ $z_{\text{int}}$ を急激に押し下げる「ピストン効果」が生じる。
現在の方程式にはこのソース項が存在しない。

**拡張方針:** 界面質量保存式に気化体積生成速度を追加する:

$$
A_{\text{floor}} \, \frac{dz_{\text{int}}}{dt}
= -\frac{\dot{m}_p}{\rho_{\text{upper}}} + \dot{V}_{\text{return}} - \dot{V}_{\text{steam}}
$$

ここで:

$$
\dot{V}_{\text{steam}} = \frac{\dot{m}_{\text{water}} \, R \, T_{\text{upper}}}{p \, M_{w,\text{H}_2\text{O}}}
$$

- $\dot{m}_{\text{water}}$: 水の蒸発速度 [kg/s]（ロウリュ投入量と蒸発モデルから決定）
- $M_{w,\text{H}_2\text{O}} = 18.015$ g/mol
- 蒸発速度のモデルは Spalding の質量移動理論を使用予定

上層エネルギー保存式にも蒸発潜熱の寄与を追加する:

$$
Q_{\text{steam}} = \dot{m}_{\text{water}} \, L_v
$$

ここで $L_v = 2.26 \times 10^6$ J/kg（水の蒸発潜熱, 100°C）。
この項は上層の温度を一時的に上昇させつつ、界面を押し下げ、
湿度ピーク（K-03）と到達時間（K-04）の両方に影響する。

### 9.2 アウフグース: ROM (低次元化モデル) によるハイブリッド・アプローチ (Phase 3)

**問題:** 「簡易版では対応不可」では、シミュレータ���としての価値が不十分。
0D モデルで速度場（タオルによるジェット流れ）を解くことは不可能だが、
「上層の熱を強制的に下層へ引きずり下ろす」効果をパラメータ化することは可能。

#### 拡張方針: ハイブリッド ROM (Reduced Order Model)

1. **事前計算:** OpenFOAM で複数パターンのアウフグース（タオル振り周波数、振幅）を計算
2. **パラメータ抽出:** 各条件における **層間強制混合係数 $\beta_{\text{aug}}$** [kg/s] を抽出。
   $\beta_{\text{aug}}$ はタオルが生成する強制循環による層間質量交換率を表す
3. **簡易版への統合:** アウフグース実行時のみ有効な強制混合項を追加:

上層エネルギー保存:

$$
\rho_{\text{upper}} c_p V_{\text{upper}} \frac{dT_{\text{upper}}}{dt}
= \cdots - \beta_{\text{aug}} \, c_p \, (T_{\text{upper}} - T_{\text{lower}})
$$

下層エネルギー保存:

$$
\rho_{\text{lower}} c_p V_{\text{lower}} \frac{dT_{\text{lower}}}{dt}
= \cdots + \beta_{\text{aug}} \, c_p \, (T_{\text{upper}} - T_{\text{lower}})
$$

$\beta_{\text{aug}}$ は OpenFOAM の結果から次のように抽出:

$$
\beta_{\text{aug}} = \frac{\dot{Q}_{\text{forced}}}{ c_p (T_{\text{upper}} - T_{\text{lower}})}
$$

ここで $\dot{Q}_{\text{forced}}$ はアウフグースにより追加的に輸送された熱量（CFD 結果から計測）。
これにより、計算負荷を上げずに「熱波が来る」現象を K-05（顔面風速ピーク）と
K-06（簡易熱ストレス指標）に反映できる。

### 9.3 正確版の多成分化: 組成浮力の導入 (Phase 2 で必須)

**問題:** 現在の `pureMixture` + `perfectGas`（$M_w = 28.96$ 固定）は、
水蒸気（$M_w = 18.015$）と乾燥空気の混合による密度変化を表現できない。
ロウリュ時のプルーム急上昇は、温度差による **熱的浮力** に加え、
蒸気濃度差による **組成浮力（solutal buoyancy）** の寄与が大きい。
水蒸気は空気より軽い（$18/29 \approx 0.62$）ため、蒸気リッチなプルームは
温度差だけで予測されるより速く上昇する。

**拡張方針:** Phase 2 で熱物理モデルを変更:

```c
thermoType
{
    type            heRhoThermo;
    mixture         multiComponentMixture;   // ← pureMixture から変更
    transport       sutherland;              // ← 温度依存粘性
    thermo          janaf;                   // ← 温度依存 Cp
    equationOfState perfectGas;
    specie          specie;
    energy          sensibleEnthalpy;
}
```

蒸気質量分率 $Y$ の輸送方程式を追加:

$$
\frac{\partial (\rho Y)}{\partial t} + \nabla \cdot (\rho \mathbf{U} Y)
= \nabla \cdot (D_{\text{eff}} \nabla Y) + S_Y
$$

混合密度:

$$
\rho = \frac{p}{R T} \left( \frac{Y}{M_{w,\text{H}_2\text{O}}} + \frac{1-Y}{M_{w,\text{air}}} \right)^{-1}
$$

ここで $D_{\text{eff}} = D_{\text{mol}} + D_t$（分子拡散 + 乱流拡散）。
$D_{\text{mol}} \approx 2.5 \times 10^{-5}$ m²/s（空気中の水蒸気拡散係数, 373K）。

これにより密度が温度 $T$ と蒸気質量分率 $Y$ の両方に依存し、
組成浮力が自動的に運動量方程式に反映される。

### 9.4 輻射モデルの改善

**問題:** 簡易版の下層エネルギー保存における $Q_{\text{rad}} \times 0.3$（固定係数）は
以下の感度を失っている:

- ストーブの配置（部屋の中央 vs 隅）
- 部屋の形状・アスペクト比
- ベンチや障害物による遮蔽

サウナにおける人体の受熱量の半分以上は輻射によるものであり（特に高温域）、
固定係数では熱ストレス指標（K-06）の空間的分布を正しく予測できない。

**拡張方針 (簡易版):**

ヒーター面を放射源、各受熱面（上部壁、下部壁、床、ベンチ）を受熱体とした
**形態係数（View Factor）** を幾何学的に概算する:

$$
F_{ij} = \frac{1}{A_i} \int_{A_i} \int_{A_j} \frac{\cos\theta_i \cos\theta_j}{\pi r^2} dA_j dA_i
$$

簡易な直方体近似では、対向する平行平板の形態係数の解析解を利用できる
（Hottel & Sarofim, 1967; Howell et al., "Thermal Radiation Heat Transfer", 6th ed.）。

下層の輻射受熱を形態係数で置き換える:

$$
Q_{\text{rad,lower}} = Q_{\text{rad}} \times (F_{\text{heater} \to \text{floor}} + F_{\text{heater} \to \text{lower walls}})
$$

$Q_{\text{rad}} \times 0.3$ の代わりにこれを使用すれば、部屋の形状やヒーター配置の変更に対して
自動的に輻射配分が変化する。

**拡張方針 (OpenFOAM):**

Phase 2 以降で輻射モデルを導入する:

```c
radiationModel  fvDOM;    // Finite Volume Discrete Ordinates Method
// または
radiationModel  viewFactor;  // 壁面間のみの輻射交換（計算コスト低）

absorptionEmissionModel  greyMeanAbsorptionEmission;
scatterModel    none;
```

特に高温域（$T > 100$ °C）では壁面輻射による二次加熱（ストーブ → 壁 → 空気層 → 反対側の壁）
が自然対流の境界層形成に影響するため、Phase 2 では必須。
`viewFactor` モデルは壁面間のみの輻射交換を計算するため、計算コストが低く
PoC フェーズに適している。

---

## 10. 既知の物理的制約（未解決）

本プロジェクトの現行モデルに存在する物理的制約のうち、
認識済みだが未修正のものを以下に列挙する。
これらは今後のフェーズで段階的に対処する予定。

### 10.1 換気モデルの欠落

**現状:** サウナ室を完全密閉空間として扱っている。

**物理的問題:** 実際のサウナには給排気口が存在し、ロウリュ時の体積膨張は主に
排気口から外部へ放出される。密閉仮定では界面低下を過大評価する。
また、定常時の新鮮空気導入は下層温度と湿度に決定的な影響を与える。

**今後の方針:** オリフィス流量式
$\dot{m}_{\text{vent}} = C_d A \sqrt{2 \rho \Delta p}$
に基づく換気モデルを質量・エネルギー保存式に組み込む。
サウナの給排気口サイズ（通常 $A \sim 0.01\text{--}0.04$ m²）と
位置（吸気口: 低位、排気口: 高位）をYAMLパラメータ化する。

### 10.2 Aufguss ROM の質量保存との不整合

**現状:** $\beta_{\text{aug}}$ によるエネルギー交換のみをモデル化。

**物理的問題:** 対流による熱輸送は必ず質量の移動を伴う。
上層から下層へ空気が移動すれば、下層から上層へも同質量が押し出される。
エネルギー保存式だけに混合項を入れて質量保存式（界面方程式）に
反映しないと、ゾーンの質量とエンタルピーの整合性が破綻する。

**今後の方針:** 界面質量保存式にも $\dot{m}_{\text{mix}}$ に基づく
双方向質量交換項を追加する。

### 10.3 Zukoski プルームモデルの仮想原点未補正

**現状:** ヒーター中心位置を直接 $z=0$ として Zukoski 相関式に入力。

**物理的問題:** Zukoski の式は「点熱源」を仮定している。
有限サイズのサウナストーブでは、物理的なヒーター位置と
プルーム計算上の仮想原点 $z_0$ にずれが生じる。
火災（数百～千℃）とサウナ（数十～数百℃）では密度差のスケールも異なり、
Boussinesq 近似の妥当性限界を超える可能性がある。

**今後の方針:** Heskestad (1984) の仮想原点補正
$z_0 = -1.02 D + 0.083 Q^{2/5}$ を導入し、
$z$ → $z - z_0$ に置き換える。

### 10.4 体感温度モデル (皮膚熱収支モデルに置換済)

**旧実装:** K-06 に Steadman の heat index 近似式を使用していたが、
気温 20-50℃ 程度の屋外条件のみ有効であり、サウナ温度域 (60-120℃) では適用外だった。

**現在の実装:** 皮膚表面の熱収支を直接評価する等価温度モデルに置換済。

- 対流熱流束: $q_{\text{conv}} = h_{\text{conv}} (T_{\text{air}} - T_{\text{skin}})$
- 蒸発/凝縮項: Lewis 関係により湿度の影響を評価
- ヒーター直接放射: $q_{\text{rad,body}} = \varepsilon \sigma F_{\text{body}} (T_{\text{heater}}^4 - T_{\text{skin}}^4)$
- 等価温度: $T_{\text{eq}} = T_{\text{skin}} + q_{\text{total}} / h_{\text{ref}}$

60-120℃ の範囲で物理的に妥当な値を返す。

### 10.5 SST k-omega の浮力生成項

**現状:** 標準の SST k-omega モデルを使用（`kOmegaSST`）。

**物理的問題:** 標準の SST k-omega は浮力による乱流生成/消散の項
（buoyancy production term $G_b = -\frac{\mu_t}{Pr_t \rho} \mathbf{g} \cdot \nabla\rho$）
を含んでいない。安定成層（上層が高温）では浮力が乱流を抑制するが、
この効果が無視されると温度成層の鋭さを過小予測する。

**今後の方針:** OpenFOAM の乱流モデルに浮力項を追加するか、
`buoyantKOmegaSST` 相当のモデルを採用する。

### 10.6 輻射経路の設計判断

**現状:** `wall_cfg == "lumped"` の場合、全輻射は壁面に吸収され、
壁面から空気への対流伝熱で間接的に空気を加熱する。
`wall_cfg == "fixed"` の場合は下層空気に直接入力される。

**設計判断:** 実際のサウナでは輻射は壁面・床・天井に吸収され、
壁面温度を上げて二次的に空気を加熱する。
現在の lumped モデルはこの物理を近似的に表現している。
人体が直接受ける輻射（ヒーターからの直達輻射）は
皮膚熱収支モデルの $q_{\text{rad,body}}$ 項として実装済であり、
ヒーター--人体間の形態係数 $F_{\text{body}}$ を通じて
K-06（熱ストレス）の計算に反映される。

---

## 付録 A: フェーズ別の実装状況

| フェーズ | OpenFOAM | 簡易版 | 状態 |
|---------|----------|--------|------|
| Phase 1 | buoyantPimpleFoam + SST k-omega, pureMixture | 2-Zone プルームモデル (定常) | **実装済** |
| Phase 2 | multiComponentMixture + 蒸気輸送 ($Y$) + H2O テンプレート | 蒸気体積膨張 $\dot{V}_{\text{steam}}$ + 蒸発潜熱 $Q_{\text{steam}}$ + 形態係数輻射 + 壁面昇温モデル + 湿度連成物性 + 非定常ソルバー (`solve_transient`) | **実装済** |
| Phase 3 | 運動量ソース項 (ジェット) — 未実装 | ROM: 強制混合係数 $\beta_{\text{aug}}$ + 抽出スクリプト (`extract_beta_aug.py`) | **枠組み実装済** |
| Phase 4-5 | — | CSV 取込 (`validation.py`) + プローブ比較 (`compare_probes`) + レポート生成 (`reporting.py`) | **枠組み実装済** |

## 付録 B: KPI 一覧と実装状況

| KPI ID | 名称 | 入力 | 実装 |
| ------ | ---- | ---- | ---- |
| K-01 | 定常温度差 (上段-下段) | プローブ定常値 | **済** (`kpi.py`) |
| K-02 | Löyly 後ピーク温度上昇 | 時系列 $T_{\text{upper}}(t)$ | **済** |
| K-03 | Löyly 後ピーク絶対湿度 | 時系列 $w(t)$ | **済** |
| K-04 | ピーク到達時間 | 時系列 $T_{\text{upper}}(t)$ + イベント時刻 | **済** |
| K-05 | 顔面風速ピーク (proxy) | $\beta_{\text{aug}}$ から推定 | **済** (ROM proxy) |
| K-06 | 簡易熱ストレス指標 | 体感温度 (皮膚熱収支モデル+放射) | **済** |
| K-07 | 上下相対温度差 | プローブ定常値 | **済** (`kpi.py`) |

## 付録 C: 簡易版の追加モデル (Phase 2 実装分)

### 壁面集中定数モデル (lumped wall)

壁面内面温度 $T_{\text{wall,inner}}$ を状態変数として追加:

$$
(\rho c_p)_w \, \delta_w \, A_{\text{wall}} \, \frac{dT_{\text{wall,inner}}}{dt}
= \underbrace{Q_{\text{conv→wall}}}_{\text{空気→壁面対流}}
+ \underbrace{Q_{\text{rad→wall}}}_{\text{ヒーター輻射}}
- \underbrace{\frac{\lambda_w}{\delta_w} A_{\text{wall}} (T_{\text{wall,inner}} - T_{\text{wall,outer}})}_{\text{壁面→外部熱伝導}}
$$

パラメータ:

- $\delta_w = 0.015$ m (木製パネル厚さ)
- $\lambda_w = 0.12$ W/(m·K) (木の熱伝導率)
- $(\rho c_p)_w = 5 \times 10^5$ J/(m³·K) (木の体積熱容量)

`model: fixed` で従来の固定壁温、`model: lumped` で昇温モデルを選択。

### 湿度連成物性モデル

混合気体の物性を湿度比 $w$ [kg/kg] に応じて補正:

$$
c_{p,\text{mix}} = (1 - y_v) \, c_{p,\text{air}} + y_v \, c_{p,\text{vapor}}
$$

$$
h_{\text{wall,eff}} = h_{\text{wall,base}} \left(\frac{c_{p,\text{mix}}}{c_{p,\text{air}}}\right)^{0.25} \left(\frac{\lambda_{\text{mix}}}{\lambda_{\text{air}}}\right)^{0.75}
$$

ここで $y_v = w / (1+w)$ は蒸気質量分率。

### 体感温度 (皮膚熱収支モデル)

従来の Steadman (1979) 近似は屋外条件 (20-50℃) のみ有効であり、
サウナ温度域 (60-120℃) では適用外となる。
本モデルでは皮膚表面の熱収支を直接評価し、等価温度 $T_{\text{eq}}$ を算出する。

**対流熱流束:**

$$
q_{\text{conv}} = h_{\text{conv}} (T_{\text{air}} - T_{\text{skin}})
$$

ここで $T_{\text{skin}} = 36$℃、$h_{\text{conv}} = 8 \; \text{W/(m²K)}$。

**蒸発 / 凝縮項 (Lewis 関係):**

$p_{\text{vapor}} > p_{\text{sat,skin}}$ (凝縮) のとき:

$$
q_{\text{evap}} = \frac{16.5 \, h_{\text{conv}} (p_{\text{vapor}} - p_{\text{sat,skin}})}{1000}
$$

$p_{\text{vapor}} \leq p_{\text{sat,skin}}$ (蒸発冷却) のとき:

$$
q_{\text{evap}} = -\frac{16.5 \, h_{\text{conv}} (p_{\text{sat,skin}} - p_{\text{vapor}})}{1000}
$$

**ヒーターからの直接放射:**

$$
q_{\text{rad,body}} = \varepsilon_{\text{body}} \sigma F_{\text{body}} (T_{\text{heater}}^4 - T_{\text{skin}}^4)
$$

ここで $\varepsilon_{\text{body}} = 0.97$、$F_{\text{body}}$ はヒーター--人体間の形態係数
(小面積近似 $F \approx A_{\text{body}} / (2\pi d^2)$ で算出)。

**等価温度:**

$$
T_{\text{eq}} = T_{\text{skin}} + \frac{q_{\text{conv}} + q_{\text{rad,body}} + q_{\text{evap}}}{h_{\text{ref}}}
$$

ここで $h_{\text{ref}} = 10 \; \text{W/(m²K)}$。この等価温度は K-06 (簡易熱ストレス指標) に直接使用される。

### ヒーター容量の目安

業界標準 (Harvia, HELO): 約 **1 kW / m³** のサウナ室体積。

| 部屋サイズ | 体積 | 推奨ヒーター |
| --------- | ---- | ----------- |
| 2.0×1.5×2.2 m (1-2人) | 6.6 m³ | 6-8 kW |
| 3.0×2.5×2.5 m (4-6人) | 18.8 m³ | 18-20 kW |
| 4.0×3.0×2.5 m (商業用) | 30 m³ | 30 kW |
