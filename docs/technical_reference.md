# SaunaFlow 技術リファレンス — 定数・相関式・支配方程式の根拠と離散化

本ドキュメントでは、SaunaFlow で使用するすべての物理定数、相関式、支配方程式について
**出典・分野・理論的根拠・本計算での変形**を体系的に整理する。

---

## 目次

1. [物理定数の出典と根拠](#1-物理定数の出典と根拠)
2. [相関式・経験式の理論的背景](#2-相関式経験式の理論的背景)
3. [簡易版 (2-Zone) の方程式導出](#3-簡易版-2-zone-の方程式導出)
4. [正確版 (OpenFOAM) の支配方程式](#4-正確版-openfoam-の支配方程式)
5. [有限体積法 (FVM) による離散化](#5-有限体積法-fvm-による離散化)
6. [離散化スキームの詳細](#6-離散化スキームの詳細)
7. [圧力-速度連成アルゴリズム](#7-圧力-速度連成アルゴリズム)
8. [乱流モデル SST k-omega の詳細](#8-乱流モデル-sst-k-omega-の詳細)
9. [壁面関数の理論](#9-壁面関数の理論)
10. [線形ソルバーと緩和係数の根拠](#10-線形ソルバーと緩和係数の根拠)

---

## 1. 物理定数の出典と根拠

### 1.1 空気の熱力学物性

| 定数 | 記号 | 値 | 単位 | 出典 | 分野 | 根拠 |
|------|------|-----|------|------|------|------|
| 分子量 | $M_w$ | 28.96 | g/mol | NIST Chemistry WebBook | 化学物性 | 乾燥空気の平均分子量。N₂ (28.01, 78%) + O₂ (32.00, 21%) + Ar (39.95, 0.93%) の組成加重平均 |
| 定圧比熱 | $c_p$ | 1005 | J/(kg·K) | Incropera et al., "Fundamentals of Heat and Mass Transfer", 7th ed., Table A.4 | 熱力学 | 300K における乾燥空気の定圧比熱。サウナ温度域 (293-373K) での変動は ±2% 以内のため定数扱いが妥当 |
| 動粘性係数 | $\mu$ | 1.8×10⁻⁵ | Pa·s | Sutherland の式から 300K の値。Incropera, Table A.4 | 流体力学 | Sutherland の式: $\mu = \mu_0 (T/T_0)^{3/2} (T_0 + S)/(T + S)$ ($\mu_0 = 1.716 \times 10^{-5}$, $T_0 = 273.15$ K, $S = 110.4$ K)。本計算では温度依存性を無視し const transport で代用 |
| プラントル数 | $Pr$ | 0.7 | — | 空気の実験値。Schlichting, "Boundary Layer Theory", 8th ed. | 伝熱工学 | $Pr = \nu / \alpha = \mu c_p / \lambda$。気体は分子運動論から $Pr \approx 4\gamma/(9\gamma-5)$（$\gamma = 1.4$ → $Pr \approx 0.74$）。実験値は 0.707 (300K)。温度依存小 |
| 熱伝導率 | $\lambda$ | 0.026 | W/(m·K) | $\lambda = \mu c_p / Pr$ から導出 | 伝熱工学 | 厳密には $0.0258$ だが、簡易版コードでは $0.026$ に丸めて使用。OpenFOAM では $\mu$ と $Pr$ から内部計算 |
| 重力加速度 | $g$ | 9.81 | m/s² | WGS 84 標準重力 | 力学 | 日本の緯度 (35°N) での値は 9.798 m/s² だが、工学的に 9.81 を使用 |

**本計算での使用:**
- OpenFOAM: `thermophysicalProperties` に $M_w$, $c_p$, $\mu$, $Pr$ を設定 → 内部で $\lambda$, $\alpha$ を自動算出
- 簡易版: $\rho_0$, $c_p$, $g$ を直接使用

### 1.2 簡易版固有のパラメータ

| 定数 | 記号 | 値 | 単位 | 出典 | 分野 | 根拠 |
|------|------|-----|------|------|------|------|
| 基準密度 | $\rho_0$ | 1.1 | kg/m³ | 理想気体の状態方程式 | 気体物性 | $\rho = p M_w / (R T) = 101325 \times 0.02896 / (8.314 \times 300) = 1.177$。$\rho_0 = 1.1$ は ~300K での丸め値。厳密ではないが簡易モデルの精度要求に対して十分 |
| 壁面熱伝達率 | $h_{\text{wall}}$ | 8.0 | W/(m²·K) | Churchill & Chu (1975), "Correlating equations for laminar and turbulent free convection from a vertical plate", Int. J. Heat Mass Transfer 18:1323 | 伝熱工学（自然対流） | 鉛直平板の自然対流相関式: $\overline{Nu} = 0.68 + 0.670 Ra^{1/4} / [1+(0.492/Pr)^{9/16}]^{4/9}$。サウナ条件 ($\Delta T \approx 50$ K, $L = 2.5$ m) → $Ra \approx 10^{10}$ → $\overline{Nu} \approx 200$ → $h \approx Nu \cdot \lambda / L \approx 2$ W/(m²·K) だが、乱流遷移後で $h = 5\text{--}10$ 程度。$h = 8$ は高めの見積もり（安全側） |
| 対流熱割合 | $f_{\text{conv}}$ | 0.7 | — | Drysdale, "An Introduction to Fire Dynamics", 3rd ed., Ch.10 | 火災工学・燃焼工学 | ヒーター/火炎の全発熱量のうち対流で上昇気流に乗る割合。電気ヒーター（サウナストーブ）は輻射率が比較的低く、対流優勢。文献では 0.6–0.8 の範囲。$f_{\text{conv}} = 0.7$ は Zukoski (1978) が防火工学で採用した標準値 |
| 界面伝導率 | $k_{\text{int}}$ | 0.5 | W/(m·K) | モデルパラメータ（理論値ではない） | ゾーンモデル工学 | 2層界面を横切る熱輸送を実効熱伝導で近似。空気の真の熱伝導率 (0.026) の約20倍 → 界面での乱流混合・エントレインメント効果を吸収する経験的パラメータ。Cooper (1982) の2層モデルでも同様の扱い |
| 乱流プラントル数 | $Pr_t$ | 0.85 | — | Reynolds analogy の修正。Kays & Crawford, "Convective Heat and Mass Transfer", 4th ed. | 乱流伝熱 | 乱流における運動量輸送と熱輸送の比。$Pr_t = \nu_t / \alpha_t$。理論的には Prt ≈ 1（Reynolds analogy）だが、実験的に $Pr_t = 0.85\text{--}0.9$ が最適。OpenFOAM の `alphatWallFunction` で使用 |

### 1.3 壁面・湿度モデルパラメータ (Phase 2 追加分)

| 定数 | 記号 | 値 | 単位 | 出典 | 分野 | 根拠 |
|------|------|-----|------|------|------|------|
| 木材パネル厚さ | $\delta_w$ | 0.015 | m | 建築仕様 | 建築工学 | サウナ内装の杉/ヒノキ板厚。15-20mm が標準 |
| 木材熱伝導率 | $\lambda_w$ | 0.12 | W/(m·K) | Incropera, Table A.3 | 伝熱工学 | 針葉樹の繊維方向。樹種により 0.10-0.16 |
| 木材体積熱容量 | $(\rho c_p)_w$ | 5×10⁵ | J/(m³·K) | 木材工学ハンドブック | 材料工学 | $\rho \approx 450$ kg/m³, $c_p \approx 1100$ J/(kg·K) → $\rho c_p \approx 5 \times 10^5$ |
| 水蒸気比熱 | $c_{p,\text{vapor}}$ | 1860 | J/(kg·K) | NIST Chemistry WebBook | 熱力学 | 373K における水蒸気定圧比熱 |
| 水蒸気熱伝導率 | $\lambda_{\text{vapor}}$ | 0.025 | W/(m·K) | Incropera, Table A.6 | 伝熱工学 | 373K における水蒸気。空気(0.026)と同程度 |
| 蒸発潜熱 | $L_v$ | 2.26×10⁶ | J/kg | CRC Handbook | 熱力学 | 100°C における水の蒸発潜熱 |
| ヒーター容量目安 | — | ~1 | kW/m³ | Harvia サウナカリキュレーター | 設備工学 | 業界標準 (Harvia, HELO, HUUM)。3.0×2.5×2.5m で 18kW |

### 1.4 体感温度・輻射・換気モデルパラメータ

| 定数 | 記号 | 値 | 単位 | 出典 | 分野 | 根拠 |
|------|------|-----|------|------|------|------|
| 平均皮膚温度 | $T_{\text{skin}}$ | 36.0 | °C | ISO 7933 (予測熱ひずみ) | 生理学 | 温熱環境下の平均皮膚表面温度。個人差あり（34-37°C） |
| 皮膚濡れ率 | $w_{\text{skin}}$ | 0.4 | — | Gagge et al. (1971) | 温熱生理学 | 発汗による皮膚表面の湿潤度。0=完全乾燥、1=完全湿潤。サウナでは0.3-0.6程度 |
| 最大蒸発冷却 | $Q_{\text{evap,max}}$ | 400 | W/m² | Kerslake (1972), "The Stress of Hot Environments" | 温熱生理学 | 生理的限界。汗腺の最大発汗能力と皮膚面積から |
| Lewis 係数 | — | 16.5 | W/(m²·kPa) per h | Lewis (1922) の関係 | 伝熱工学 | 空気-水系の Lewis 数 ≈ 1 から $h_e \approx 16.5 h_c$ |
| ヒーター射出率 | $\varepsilon_{\text{heater}}$ | 0.90 | — | 石/金属表面の文献値 | 輻射伝熱 | サウナストーンの表面放射率。粗面石は0.85-0.95 |
| 人体射出率 | $\varepsilon_{\text{body}}$ | 0.97 | — | Hardy & Muschenheim (1934) | 生物物理学 | 人体皮膚の遠赤外域放射率。ほぼ黒体 |
| Stefan-Boltzmann定数 | $\sigma$ | 5.67×10⁻⁸ | W/(m²·K⁴) | CODATA 2018 | 物理学 | 黒体輻射の基本定数 |
| オリフィス流量係数 | $C_d$ | 0.6 | — | Idelchik, "Handbook of Hydraulic Resistance" | 流体力学 | 鋭縁オリフィスの標準値。丸孔で0.6-0.65 |

### 1.5 状態方程式

**理論根拠: 理想気体の状態方程式**

$$
pV = nRT \quad \Rightarrow \quad \rho = \frac{p M_w}{R T}
$$

- **出典:** Clausius (1857) 以降の気体運動論。熱力学の基本法則
- **分野:** 熱力学
- **適用条件:** 圧力が臨界圧力より十分低く（空気: $p_c = 3.77$ MPa）、温度が臨界温度より十分高い（空気: $T_c = 132.5$ K）場合に成立。サウナ条件（1 atm, 293–373K）では誤差 < 0.1%

**OpenFOAM での実装:**

```
equationOfState  perfectGas;   // rho = p * Mw / (R * T)
```

密度が温度と圧力の関数 → 温度上昇 → 密度低下 → 浮力発生。これが自然対流の駆動力。

**簡易版での変形:**

$$
\rho = \rho_0 \frac{T_{\text{ref}}}{T}
$$

圧力一定（密閉空間だが圧力変動は $O(10)$ Pa、大気圧 $10^5$ Pa に対して微小）の仮定の下、
$\rho \propto 1/T$ → $\rho / \rho_0 = T_0 / T$。ここで $T_{\text{ref}} = T_{\text{wall}}$, $\rho_0 = 1.1$ (対応温度 ~300K)。

---

## 2. 相関式・経験式の理論的背景

### 2.1 Zukoski プルーム相関式

**原式:**

$$
\dot{m}_p(z) = 0.071 \, Q_c^{1/3} \, z^{5/3} + 0.0018 \, Q_c
$$

($Q_c$ [kW], $z$ [m], $\dot{m}_p$ [kg/s])

**出典:**
- Zukoski, E.E. (1978), "Development of a Stratified Ceiling Layer in the Early Stages of a Closed-Room Fire", NBS-GCR-78-150
- Zukoski, E.E., Kubota, T., & Cetegen, B. (1981), "Entrainment in fire plumes", Fire Safety Journal 3(3):107-121

**分野:** 火災工学（Fire Plume Dynamics）、建築環境工学

**理論的導出:**

この式は Morton, Taylor & Turner (1956) のエントレインメント仮説から導かれる。

**Step 1: MTT エントレインメント仮説**

プルーム縁での巻き込み速度 $v_e$ が上昇速度 $w$ に比例する:

$$
v_e = \alpha_e \, w
$$

ここで $\alpha_e \approx 0.1\text{--}0.15$ がエントレインメント係数。

- **出典:** Morton, B.R., Taylor, G.I. & Turner, J.S. (1956), "Turbulent Gravitational Convection from Maintained and Instantaneous Sources", Proc. R. Soc. A 234:1-23
- **根拠:** Taylor のself-similarity仮説 — 乱流プルーム断面の速度・温度プロファイルが高さによらず相似形を保つ。実験で広く検証済み

**Step 2: プルームの自己相似解**

MTT モデルの支配方程式（質量・運動量・浮力フラックス保存）:

$$
\frac{d}{dz}(\rho_\infty b^2 w) = 2 \alpha_e \rho_\infty b w \quad \text{(質量)}
$$

$$
\frac{d}{dz}(\rho_\infty b^2 w^2) = g \frac{\Delta\rho}{\rho_\infty} b^2 \quad \text{(運動量)}
$$

$$
\frac{d}{dz}(\rho_\infty b^2 w \Delta T) = 0 \quad \text{(浮力フラックス保存)}
$$

ここで $b$ はプルーム半径。Boussinesq 近似下で $\Delta\rho/\rho \approx -\Delta T/T$。

**Step 3: べき乗則解**

上記 ODE 系のべき乗則解 $b \propto z$, $w \propto z^{-1/3}$ から:

$$
\dot{m}_p = C \, Q_c^{1/3} \, z^{5/3}
$$

ここで $C$ は $\alpha_e$, $g$, $T_\infty$, $c_p$ の関数。

**Step 4: Zukoski の実験的修正**

Zukoski は実寸大の火災プルーム実験から:
- 第1項 $0.071 Q_c^{1/3} z^{5/3}$: MTT 理論のべき乗則（遠方場）
- 第2項 $0.0018 Q_c$: 火源近傍の有限サイズ効果の補正（$z$ が小さい時に効く）

$C = 0.071$ は Zukoski の実験データへのフィッティングから決定。
Heskestad (1984) は同様の理論から $C = 0.071$ を独立に導出しており、この値の信頼性は高い。

**本計算での使用** (`simple_solver.py:78-81`):

```python
q_kw = q_conv_w / 1000.0   # W → kW に変換
m_dot = 0.071 * q_kw ** (1/3) * z ** (5/3) + 0.0018 * q_kw
```

入力は対流熱出力 $Q_c = f_{\text{conv}} \times Q_{\text{heater}}$ (kW単位) と
ヒーター中心からの高さ $z$ (m)。出力はプルーム質量流量 (kg/s)。

**サウナへの適用妥当性:**
- Zukoski の式は point source plume（点源プルーム）が前提。サウナヒーター (0.6m×0.5m) の
  寸法は室内高さ (2.5m) に対して十分小さいため、界面高さ (~1.5m) での適用は妥当
- 密閉空間効果（ceiling jet, 壁面反射）は本式に含まれないが、2-Zone モデルの
  壁面熱損失項で間接的にモデル化

### 2.2 プルーム温度（エネルギー保存）

$$
T_{\text{plume}} = T_{\text{lower}} + \frac{Q_c}{\dot{m}_p \, c_p}
$$

**理論根拠:** プルーム断面積分したエネルギー保存則（第一法則）

$$
\dot{m}_p c_p T_{\text{plume}} = \dot{m}_p c_p T_{\text{lower}} + Q_c
$$

→ プルームが下層空気を $\dot{m}_p$ [kg/s] の速度でエントレインし、
$Q_c$ [W] の対流熱を受け取る。巻き込みが多い（$z$ が大きい）ほど
$\dot{m}_p$ が増え $T_{\text{plume}}$ は低下する → MTT 理論の本質的特徴。

- **出典:** 熱力学第一法則の定常開放系への適用。成書では Bird, Stewart & Lightfoot, "Transport Phenomena", 2nd ed., Ch.15
- **分野:** 熱流体力学

### 2.3 Newton の冷却則（壁面熱損失）

$$
Q_{\text{wall}} = h_{\text{wall}} \cdot A_{\text{wall}} \cdot (T_{\text{fluid}} - T_{\text{wall}})
$$

- **出典:** Newton (1701), "Scala graduum caloris" (近似法則として。厳密な導出は対流伝熱理論)
- **分野:** 伝熱工学
- **理論:** 壁面近傍の温度境界層を通しての熱輸送を、熱伝達率 $h$ で集約。
  $h$ は自然対流の場合 Nusselt 数相関式 (Churchill & Chu, 1975) から得られるが、
  簡易版では $h = 8.0$ W/(m²·K) に固定

**Churchill-Chu 相関式（$h_{\text{wall}} = 8.0$ の根拠）:**

鉛直平板の自然対流（全域相関式）:

$$
\overline{Nu}_L = \left\{ 0.825 + \frac{0.387 Ra_L^{1/6}}{[1 + (0.492/Pr)^{9/16}]^{8/27}} \right\}^2
$$

サウナ条件での見積もり:
- $L = 2.5$ m（壁面高さ）
- $\Delta T = 50$ K（上層温度 - 壁温）
- $\beta = 1/T \approx 1/320 = 3.1 \times 10^{-3}$ K⁻¹
- $\nu = 1.6 \times 10^{-5}$ m²/s
- $Ra = g \beta \Delta T L^3 / (\nu \alpha) = 9.81 \times 3.1 \times 10^{-3} \times 50 \times 2.5^3 / (1.6 \times 10^{-5} \times 2.3 \times 10^{-5}) \approx 6.5 \times 10^{10}$
- $\overline{Nu} \approx 370$
- $h = Nu \cdot \lambda / L = 370 \times 0.026 / 2.5 \approx 3.8$ W/(m²·K)

ただし実際のサウナは:
- 壁面が木材（表面粗さによる乱流促進）
- 四方の壁が同時に冷却（角部の流れ干渉）
- 天井面の自然対流（加熱面が下向き ≠ 鉛直平板）

→ 有効 $h$ は純粋な鉛直平板よりやや大きく、$h = 5\text{--}10$ W/(m²·K) が妥当範囲。
$h = 8.0$ はこの範囲の中間値を採用したもの。

### 2.4 シグモイド補間（温度プロファイル構築）

$$
T(y) = T_{\text{lower}} + \frac{T_{\text{upper}} - T_{\text{lower}}}{1 + \exp\left(-3 \cdot \frac{y - z_{\text{int}}}{\delta / 2}\right)}
$$

- **出典:** 数学的ツール（ロジスティック関数）。物理的根拠は Cooper (1982) の温度遷移層モデル
- **分野:** 数値解析・ゾーンモデル工学
- **根拠:** 2層ゾーンモデルの集中定数 ($T_{\text{upper}}$, $T_{\text{lower}}$) を
  鉛直方向に連続分布化するための後処理。物理的には界面近傍で急峻な温度変化が
  有限の遷移層厚さ $\delta = 0.15H$ で生じることを表現。
  係数 "3" はシグモイドの傾き制御パラメータ（$\pm\delta/2$ の範囲で 5%~95% の遷移）

### 2.5 体感温度モデル（皮膚熱収支）

Steadman (1979) の heat index はサウナ温度域 (>80°C) で破綻するため、
皮膚表面の熱収支に基づく等価温度モデルに置き換えた。

**皮膚表面の熱収支:**

$$
q_{\text{total}} = \underbrace{h_c (T_{\text{air}} - T_{\text{skin}})}_{\text{対流}}
+ \underbrace{q_{\text{rad,body}}}_{\text{直達輻射}}
+ \underbrace{q_{\text{evap}}}_{\text{蒸発/結露}}
$$

**蒸発・結露項 (Lewis 関係):**

$$
q_{\text{evap}} = \begin{cases}
16.5 \, h_c \, (p_v - p_{\text{sat,skin}}) & \text{if } p_v > p_{\text{sat,skin}} \text{ (結露→加熱)} \\
-\min\left(w_{\text{skin}} \cdot 16.5 \, h_c \, (p_{\text{sat,skin}} - p_v), \; Q_{\text{evap,max}}\right) & \text{(蒸発→冷却)}
\end{cases}
$$

- $p_v = RH \cdot p_{\text{sat}}(T_{\text{air}})$: 空気中の水蒸気分圧 [kPa]
- $p_{\text{sat,skin}} = p_{\text{sat}}(T_{\text{skin}})$: 皮膚温度での飽和蒸気圧 [kPa]
- $p_{\text{sat}}$ は Magnus 式: $p = 0.61078 \exp(17.27 T / (T + 237.3))$ [kPa]

**等価温度:**

$$
T_{\text{perceived}} = T_{\text{skin}} + \frac{q_{\text{total}}}{h_{\text{ref}}}
$$

$h_{\text{ref}} = 10$ W/(m²·K) は正規化用の基準伝熱係数。

- **出典:** ISO 7933 (2004) "Ergonomics of the thermal environment" の皮膚熱収支モデルを簡略化
- **分野:** 温熱生理学・人間工学
- **Steadman からの改善点:** 結露加熱を明示的にモデル化（サウナ蒸気による皮膚結露で体感温度が急上昇する現象を表現可能）

### 2.6 ヒーターから人体への直達輻射

**Stefan-Boltzmann の法則:**

$$
q_{\text{rad,body}} = \varepsilon_{\text{body}} \, \sigma \, F_{\text{heater→body}} \, (T_{\text{heater}}^4 - T_{\text{skin}}^4)
$$

ヒーター表面温度は放射出力から逆算:

$$
T_{\text{heater}} = \left(\frac{Q_{\text{rad}} / A_{\text{heater}}}{\varepsilon_{\text{heater}} \, \sigma} + T_{\text{wall,inner}}^4 \right)^{1/4}
$$

- **出典:** Siegel & Howell, "Thermal Radiation Heat Transfer", 5th ed.
- **分野:** 輻射伝熱
- **本計算での意義:** ヒーターから人体への直達輻射は、壁面経由の間接輻射とは独立に体感温度に寄与する。これは「同じ気温でもヒーターの近くは熱い」という実体験と一致

### 2.7 換気モデル（スタック効果）

密閉仮定を撤廃し、サウナ室の給排気口を通じた自然換気をモデル化。

**スタック効果による圧力差:**

$$
\Delta p = \rho_{\text{amb}} \, g \, (z_{\text{exhaust}} - z_{\text{supply}}) \, \frac{T_{\text{col}} - T_{\text{amb}}}{T_{\text{col}}}
$$

ここで $T_{\text{col}}$ は給排気口間の空気柱の加重平均温度。

**オリフィス流量:**

$$
\dot{m}_{\text{vent}} = A_{\text{eff}} \, \sqrt{2 \rho \lvert \Delta p \rvert} \cdot \text{sign}(\Delta p)
$$

- $A_{\text{eff}} = \min(C_{d,s} A_{\text{supply}}, C_{d,e} A_{\text{exhaust}})$: $C_d$ を含む実効面積
- $\rho$: 上流側の密度（流入時は外気、流出時は室内空気）

- **出典:** ASHRAE Handbook — Fundamentals (2017), Ch.16 "Ventilation and Infiltration"
- **分野:** 建築環境工学
- **サウナへの適用:** フィンランドサウナは床近くの吸気口と天井近くの排気口を持つ設計が標準。ロウリュ時の体積膨張は主にこの排気口から放出される

---

## 3. 簡易版 (2-Zone) の方程式導出

### 3.1 出発点: 熱力学第一法則の制御体積形

閉じた制御体積（ゾーン）に対する第一法則:

$$
\frac{dE}{dt} = \dot{Q}_{\text{in}} - \dot{Q}_{\text{out}} + \dot{m}_{\text{in}} h_{\text{in}} - \dot{m}_{\text{out}} h_{\text{out}}
$$

- **出典:** Cengel & Boles, "Thermodynamics: An Engineering Approach", 8th ed., Ch.5
- **分野:** 工学熱力学

### 3.2 上層ゾーンへの適用

制御体積: 界面から天井まで。入口: プルーム（底部から上昇）。出口: 壁面沿い下降流。

$$
\rho_{\text{upper}} c_p V_{\text{upper}} \frac{dT_{\text{upper}}}{dt}
= \underbrace{\dot{m}_p c_p (T_{\text{plume}} - T_{\text{upper}})}_{\text{プルームエンタルピー流入}}
- \underbrace{h_{\text{wall}} A_{\text{wall,upper}} (T_{\text{upper}} - T_{\text{wall}})}_{\text{壁面損失}}
$$

**変形の過程:**
1. $dE/dt = \rho c_p V \, dT/dt$ （内部エネルギー変化 ≈ エンタルピー変化、$p$ 一定）
2. プルーム入力 = $\dot{m}_p c_p T_{\text{plume}} - \dot{m}_p c_p T_{\text{upper}}$
   （質量 $\dot{m}_p$ が $T_{\text{plume}}$ で流入し、上層の $T_{\text{upper}}$ と混合）
3. 壁面損失 = Newton の冷却則を上層が接する全壁面に適用
4. $A_{\text{wall,upper}} = P(H - z_{\text{int}}) + A_{\text{floor}}$
   — 周長 $P = 2(W+D)$ × 上層高さ（側壁） + 底面積（天井）

### 3.3 下層ゾーンへの適用

$$
\rho_{\text{lower}} c_p V_{\text{lower}} \frac{dT_{\text{lower}}}{dt}
= \underbrace{Q_{\text{rad,lower}}}_{\text{輻射受熱}}
+ \underbrace{k_{\text{int}} A_{\text{floor}} \frac{T_{\text{upper}} - T_{\text{lower}}}{0.1H}}_{\text{界面伝導}}
- \underbrace{h_{\text{wall,eff}} A_{\text{wall,lower}} (T_{\text{lower}} - T_{\text{wall,inner}})}_{\text{壁面損失}}
$$

**輻射経路の設計（壁面モデル依存）:**

- **`model: lumped`（推奨）:** ヒーター輻射は壁面集中定数モデルに全量吸収。
  壁面が温まり、対流で空気へ返す。→ $Q_{\text{rad,lower}} = 0$（二重計上を防止）
- **`model: fixed`:** 壁面蓄熱なし。ヒーター輻射の一部が直接下層空気へ。
  → $Q_{\text{rad,lower}} = Q_{\text{rad}} \times f_{\text{rad,lower}}$

$f_{\text{rad,lower}}$ は幾何学的 View Factor で計算（セクション 2.3 の形態係数近似）。
固定係数 0.3 は撤廃。

**各項の根拠:**

- **界面伝導**: Fourier の法則をゾーン間に適用。
  実効距離 $0.1H$、$k_{\text{int}} = 0.5$ W/(m·K) は界面乱流混合の経験的パラメータ
- **壁面損失**: $h_{\text{wall,eff}}$ は湿度連成（セクション 2.5 参照）。
  $T_{\text{wall}}$ は壁面内面温度（lumped モデル時）

### 3.4 壁面集中定数モデル

$$
(\rho c_p)_w \delta_w A_{\text{wall}} \frac{dT_{\text{wall,inner}}}{dt}
= Q_{\text{conv→wall}} + Q_{\text{rad}} - \frac{\lambda_w}{\delta_w} A_{\text{wall}} (T_{\text{wall,inner}} - T_{\text{wall,outer}})
$$

**導出:** 壁面パネルを厚さ $\delta_w$ の均一温度板として lumped capacitance 近似。
Biot 数 $Bi = h \delta_w / \lambda_w = 8 \times 0.015 / 0.12 = 1.0$ で、
厳密には lumped 近似の限界 ($Bi < 0.1$) を超えるが、木材内部の温度分布は
本モデルの主要な不確実性に比べて小さいため許容する。

- **出典:** Incropera, Ch.5 "Transient Conduction" — Lumped Capacitance Method
- **分野:** 伝熱工学

### 3.5 換気の質量・エネルギー保存への寄与

換気モデル（セクション 2.7）が有効な場合、下層エネルギー保存に以下を追加:

$$
\dot{Q}_{\text{vent}} = \dot{m}_{\text{vent}} \, c_p \, (T_{\text{ambient}} - T_{\text{lower}})
$$

界面質量保存にも換気による体積変化を追加:

$$
\dot{V}_{\text{vent}} = \dot{m}_{\text{vent}} / \rho_{\text{lower}}
$$

### 3.6 界面の質量保存

$$
A_{\text{floor}} \frac{dz_{\text{int}}}{dt} = -\frac{\dot{m}_p}{\rho_{\text{upper}}} + \dot{V}_{\text{return}} - \dot{V}_{\text{steam}} + \dot{V}_{\text{mix}} + \dot{V}_{\text{vent}}
$$

**導出:**

1. 体積保存 $V_{\text{upper}} + V_{\text{lower}} = V_{\text{total}}$ (一定、換気なし時)
2. $dV_{\text{upper}}/dt = -A_{\text{floor}} \, dz_{\text{int}}/dt$
3. 上層への体積流入 = プルーム $\dot{m}_p / \rho_{\text{upper}}$
4. 上層からの体積流出 = 壁面下降流 $\dot{V}_{\text{return}}$
5. 蒸気膨張 = $\dot{V}_{\text{steam}}$（界面を押し下げる）
6. Aufguss 混合 = $\dot{V}_{\text{mix}} = \beta_{\text{aug}} (1/\rho_{\text{lower}} - 1/\rho_{\text{upper}})$（密度差による正味体積効果）
7. 換気 = $\dot{V}_{\text{vent}} = \dot{m}_{\text{vent}} / \rho_{\text{lower}}$（新鮮空気流入で界面を押し上げる）

**リターンフロー $\dot{V}_{\text{return}}$ の根拠:**

$$
\dot{V}_{\text{return}} = \frac{h_{\text{wall}} A_{\text{wall,upper}} (T_{\text{upper}} - T_{\text{wall}})}{\rho_{\text{upper}} c_p \max(T_{\text{upper}} - T_{\text{lower}}, 1)}
$$

これは壁面冷却による密度増加 → 下降流の簡易モデル:
1. 壁面での冷却率: $\dot{Q} = h_{\text{wall}} A_{\text{wall}} \Delta T$
2. 冷却された空気が $\Delta T_{\text{layers}} = T_{\text{upper}} - T_{\text{lower}}$ だけ温度低下すると下層に沈降
3. 沈降質量流量: $\dot{m}_{\text{return}} = \dot{Q} / (c_p \Delta T_{\text{layers}})$
4. 体積流量: $\dot{V}_{\text{return}} = \dot{m}_{\text{return}} / \rho_{\text{upper}}$

- **出典:** Cooper, L.Y. (1982), "A Mathematical Model for Estimating Available Safe Egress Time in Fires", NBSIR 82-2612 — 2層モデルの壁面プルーム/下降流の扱い
- **分野:** 防火工学

### 3.5 時間積分（前進 Euler 法）

$$
T^{n+1} = T^n + \Delta t \cdot f(T^n)
$$

- **出典:** Euler (1768)。数値解析の最も基本的な ODE 積分法
- **分野:** 数値解析
- **本計算:** $\Delta t = 0.5$（擬似時間ステップ）。定常解のみが必要なため、
  物理時間の精度は不要。安定性条件 ($\Delta t < \tau_{\text{min}}$) を満たす程度に小さければよい
- **精度:** 1次精度 $O(\Delta t)$。定常解には影響しない（過渡応答の精度のみに影響）

---

## 4. 正確版 (OpenFOAM) の支配方程式

buoyantPimpleFoam が解く方程式系を、**連続形→積分形→離散形**の順に記述する。

### 4.1 連続の式（質量保存）

**微分形（連続形）:**

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{U}) = 0
$$

- **出典:** Euler (1757) の流体運動方程式、Reynolds の輸送定理 (1903) から導出
- **分野:** 流体力学（基本保存則）
- **物理:** 任意の検査体積への質量流入 = 質量蓄積。圧縮性流体では密度変化を含む

**積分形（制御体積）:**

$$
\frac{\partial}{\partial t} \int_V \rho \, dV + \oint_S \rho \mathbf{U} \cdot d\mathbf{S} = 0
$$

Gauss の発散定理 $\int_V \nabla \cdot \mathbf{F} \, dV = \oint_S \mathbf{F} \cdot d\mathbf{S}$ を用いた変形。

**OpenFOAM での離散形:**

$$
\frac{(\rho_P V_P)^{n+1} - (\rho_P V_P)^{n}}{\Delta t} + \sum_f (\rho \mathbf{U})_f \cdot \mathbf{S}_f = 0
$$

- $P$: セル中心、$f$: セル面
- $\mathbf{S}_f$: 面の外向き法線ベクトル（面積 × 方向）
- $(\rho \mathbf{U})_f \cdot \mathbf{S}_f$: 面を通る質量フラックス $\phi_f = \dot{m}_f$
- buoyantPimpleFoam では密度は状態方程式 $\rho = pM_w/(RT)$ から更新

### 4.2 運動量保存（Navier-Stokes 方程式）

**微分形:**

$$
\frac{\partial (\rho \mathbf{U})}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} \otimes \mathbf{U})
= -\nabla p_{\text{rgh}}
+ \nabla \cdot \boldsymbol{\tau}_{\text{eff}}
+ \rho \mathbf{g}
$$

ここで:

$$
\boldsymbol{\tau}_{\text{eff}} = (\mu + \mu_t) \left[ \nabla \mathbf{U} + (\nabla \mathbf{U})^T - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I} \right]
$$

- **出典:** Navier (1822) & Stokes (1845)。Newton の粘性法則 + Cauchy の運動方程式
- **分野:** 流体力学（基礎方程式）

**各項の物理的意味:**

| 項 | 物理 | 数学 |
|----|------|------|
| $\partial(\rho\mathbf{U})/\partial t$ | 運動量の時間変化率 | 非定常項 |
| $\nabla \cdot (\rho \mathbf{U} \otimes \mathbf{U})$ | 対流による運動量輸送 | 非線形移流項 |
| $-\nabla p_{\text{rgh}}$ | 圧力勾配力 | 体積力（静水圧除去済み） |
| $\nabla \cdot \boldsymbol{\tau}_{\text{eff}}$ | 粘性応力（分子+乱流） | 拡散項 |
| $\rho \mathbf{g}$ | 重力（浮力の源泉） | 体積力ソース |

**修正圧力 $p_{\text{rgh}}$ の導入:**

$$
p_{\text{rgh}} = p - \rho \mathbf{g} \cdot \mathbf{r}
$$

- **根拠:** 静水圧成分を除去することで、圧力の数値精度を改善。
  浮力駆動流では $p$ の変動が静水圧に比べて微小 ($O(10)$ Pa vs $O(10^4)$ Pa)
  のため、直接 $p$ を解くと桁落ちが生じる
- **出典:** Patankar (1980), "Numerical Heat Transfer and Fluid Flow", Ch.6
- 全壁面で `fixedFluxPressure` を適用 — $\nabla p_{\text{rgh}} \cdot \mathbf{n} = -\rho \mathbf{g} \cdot \mathbf{n}$

**積分形:**

$$
\frac{\partial}{\partial t} \int_V \rho \mathbf{U} \, dV
+ \oint_S (\rho \mathbf{U} \otimes \mathbf{U}) \cdot d\mathbf{S}
= -\oint_S p_{\text{rgh}} \, d\mathbf{S}
+ \oint_S \boldsymbol{\tau}_{\text{eff}} \cdot d\mathbf{S}
+ \int_V \rho \mathbf{g} \, dV
$$

**OpenFOAM 離散形:**

$$
\frac{(\rho \mathbf{U})_P^{n+1} V_P - (\rho \mathbf{U})_P^n V_P}{\Delta t}
+ \sum_f \phi_f \mathbf{U}_f
= -\sum_f (p_{\text{rgh}})_f \mathbf{S}_f
+ \sum_f (\mu_{\text{eff}})_f (\nabla \mathbf{U})_f \cdot \mathbf{S}_f
+ \rho_P \mathbf{g} V_P
$$

- $\phi_f = (\rho \mathbf{U})_f \cdot \mathbf{S}_f$: 質量フラックス
- $\mathbf{U}_f$: 面補間された速度（linearUpwind スキームで評価）

### 4.3 エネルギー保存

**微分形（sensibleEnthalpy 形式）:**

$$
\frac{\partial (\rho h)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} h)
= \nabla \cdot (\alpha_{\text{eff}} \nabla T) + S_h
$$

ここで:
- $h = c_p T$ （sensibleEnthalpy。`hConst` なので $c_p$ 一定）
- $\alpha_{\text{eff}} = \frac{\mu}{Pr} + \frac{\mu_t}{Pr_t} = \frac{\lambda}{c_p} + \frac{\mu_t}{Pr_t}$
- $S_h$: 熱源項（本計算ではヒーター壁面の熱流束として境界条件で付与）

- **出典:** エネルギー保存 — 熱力学第一法則の流体版。Bird, Stewart & Lightfoot, "Transport Phenomena", Ch.11
- **分野:** 伝熱流体力学

**$h = c_p T$ とする理由:**
OpenFOAM の `sensibleEnthalpy` は $h = \int c_p dT$ を使用。`hConst` ($c_p$ 一定) なので $h = c_p T$ + const。
エネルギー方程式を $T$ について解くと:

$$
\frac{\partial (\rho c_p T)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} c_p T)
= \nabla \cdot (\alpha_{\text{eff}} c_p \nabla T) + S_h
$$

$c_p$ 一定で両辺を $c_p$ で割れば温度方程式になるが、OpenFOAM ではエンタルピー形式のまま解く。

**ヒーター境界条件 — `externalWallHeatFluxTemperature`:**

$$
q_{\text{heater}} = \frac{Q_{\text{heater}}}{A_{\text{heater}}} = \frac{P_{\text{kW}} \times 1000}{W_{\text{heater}} \times H_{\text{heater}}} \quad [\text{W/m}^2]
$$

例: $P = 9$ kW, $W = 0.6$ m, $H = 0.5$ m → $q = 9000 / 0.3 = 30{,}000$ W/m²

- **分野:** 伝熱工学（境界条件）
- 壁面セルでは $q = -\lambda (\partial T / \partial n)_{\text{wall}}$ が指定される

**離散形:**

$$
\frac{(\rho h)_P^{n+1} V_P - (\rho h)_P^n V_P}{\Delta t}
+ \sum_f \phi_f h_f
= \sum_f (\alpha_{\text{eff}})_f (\nabla T)_f \cdot \mathbf{S}_f
+ S_h V_P
$$

### 4.4 方程式系のまとめ（OpenFOAM が解く 7 変数）

| # | 変数 | 方程式 | 未知数 |
|---|------|--------|--------|
| 1 | $\rho$ | 状態方程式 $\rho = pM_w/(RT)$ | 密度（代数的に決定） |
| 2-4 | $U_x, U_y, U_z$ | 運動量保存 (N-S) × 3成分 | 速度場 |
| 5 | $p_{\text{rgh}}$ | 圧力補正方程式（連続の式から導出） | 修正圧力 |
| 6 | $T$ | エネルギー保存 | 温度 |
| 7 | $k$ | 乱流運動エネルギー方程式 | 乱流 k |
| 8 | $\omega$ | 比散逸率方程式 | 乱流 omega |

→ 8 変数を各時間ステップで反復的に解く（PIMPLE ループ）

---

## 5. 有限体積法 (FVM) による離散化

### 5.1 FVM の基本原理

一般の保存則:

$$
\frac{\partial (\rho \phi)}{\partial t} + \nabla \cdot (\rho \mathbf{U} \phi) = \nabla \cdot (\Gamma \nabla \phi) + S_\phi
$$

を制御体積 $V_P$ で積分する。

- **出典:** Patankar, S.V. (1980), "Numerical Heat Transfer and Fluid Flow", Hemisphere Publishing
- **分野:** 数値流体力学 (CFD)
- **利点:** 離散レベルで保存則が厳密に成立（差分法では必ずしも成立しない）

### 5.2 時間微分項の離散化

**backward スキーム（2次精度陰的3点後退差分）:**

$$
\frac{\partial (\rho \phi)}{\partial t} \approx \frac{3(\rho\phi)_P^{n+1} - 4(\rho\phi)_P^n + (\rho\phi)_P^{n-1}}{2\Delta t} V_P
$$

- **出典:** Ferziger & Perić, "Computational Methods for Fluid Dynamics", 3rd ed., Ch.6
- **分野:** 数値解析
- **導出:** Taylor 展開 $\phi^n = \phi^{n+1} - \Delta t \phi' + \frac{\Delta t^2}{2} \phi'' - \cdots$ を
  3点で組み合わせて $\phi'$ の $O(\Delta t^2)$ 近似を構成
- **特性:** 陰的 → 無条件安定（CFL 制約なし）。2次精度 → Euler (1次) より時間精度が高い
- **本計算:** buoyantPimpleFoam（非定常）で使用。自然対流の振動解を捉えるため2次精度が必要

**比較 — 前進 Euler（簡易版で使用）:**

$$
\frac{\phi^{n+1} - \phi^n}{\Delta t} = f(\phi^n)
$$

1次精度、陽的 → CFL 制約あり。簡易版では定常解のみ求めるため十分。

### 5.3 対流項の離散化

一般形:

$$
\oint_S (\rho \mathbf{U} \phi) \cdot d\mathbf{S} \approx \sum_f \phi_f \, \underbrace{(\rho \mathbf{U})_f \cdot \mathbf{S}_f}_{\phi_f = \text{面での} \phi}
$$

問題は **$\phi_f$ の評価方法**（セル中心値 $\phi_P$, $\phi_N$ からの補間）。

#### 5.3.1 1次風上 (upwind) — k, omega に使用

$$
\phi_f = \begin{cases} \phi_P & \text{if } \dot{m}_f > 0 \\ \phi_N & \text{if } \dot{m}_f < 0 \end{cases}
$$

- **出典:** Courant, Isaacson & Rees (1952)
- **特性:** 1次精度、数値拡散大、無条件有界（非物理的振動なし）
- **本計算:** 乱流量 $k$, $\omega$ は正値であることが必須 → 安定性重視で upwind
- 対応 fvSchemes: `div(phi,k) bounded Gauss upwind;`

#### 5.3.2 2次風上 (linearUpwind) — U, T に使用

$$
\phi_f = \phi_P + (\nabla \phi)_P \cdot \mathbf{d}_{Pf}
$$

ここで $\mathbf{d}_{Pf}$ はセル中心 $P$ から面 $f$ へのベクトル。

- **出典:** Warming & Beam (1976), "Upwind second-order difference schemes", AIAA Journal
- **分野:** 数値流体力学
- **導出:** Taylor 展開で $\phi_f = \phi_P + \nabla\phi \cdot \delta\mathbf{x} + O(\delta x^2)$。
  勾配 $\nabla\phi$ をセル中心で評価し、1次風上に勾配補正を加える
- **特性:** 2次精度、数値拡散低減、有界性は保証されない（`bounded` で制限）
- 対応 fvSchemes: `div(phi,U) bounded Gauss linearUpwind grad(U);`

**`bounded` の意味:**
面フラックスの発散 $\sum_f \phi_f \dot{m}_f$ から連続の式の離散形を引くことで、
定常状態で $\sum_f \dot{m}_f \neq 0$ の場合の不安定性を除去する限定操作。

**`grad(U)` の意味:**
linearUpwind で使用する勾配場を指定。`cellLimited Gauss linear 1` が使われる。

#### 5.3.3 セル制限付き勾配 (cellLimited)

$$
(\nabla \phi)_{\text{limited}} = \min(1, \psi) \cdot (\nabla \phi)_{\text{unlimited}}
$$

$\psi$ はリミッター（セル中心値が隣接セル値の範囲を超えないよう制限）。

- **出典:** Barth & Jespersen (1989), "The design and application of upwind schemes on unstructured meshes"
- **分野:** 数値流体力学
- 対応 fvSchemes: `grad(U) cellLimited Gauss linear 1;`（リミッター係数 = 1 = 完全制限）

### 5.4 拡散項の離散化

$$
\oint_S (\Gamma \nabla \phi) \cdot d\mathbf{S} \approx \sum_f \Gamma_f \frac{\phi_N - \phi_P}{|\mathbf{d}_{PN}|} |\mathbf{S}_f| + \text{non-orthogonal correction}
$$

- **出典:** Jasak, H. (1996), "Error Analysis and Estimation for the Finite Volume Method with Applications to Fluid Flows", PhD thesis, Imperial College London
- **分野:** 数値流体力学

**非直交補正 (corrected):**

直交メッシュでは $\mathbf{S}_f \parallel \mathbf{d}_{PN}$ だが、非直交メッシュでは
$\mathbf{S}_f$ と $\mathbf{d}_{PN}$ にずれが生じる。`corrected` は:

$$
(\nabla \phi)_f \cdot \mathbf{S}_f = \frac{\phi_N - \phi_P}{|\mathbf{d}_{PN}|} |\mathbf{S}_f| \cos\theta + (\nabla \phi)_f \cdot \mathbf{k}_f
$$

ここで $\mathbf{k}_f = \mathbf{S}_f - \frac{|\mathbf{S}_f|}{|\mathbf{d}_{PN}|} \mathbf{d}_{PN}$ が非直交補正ベクトル。

- 対応 fvSchemes: `laplacian default Gauss linear corrected;`

### 5.5 面法線勾配 (snGrad)

$$
\frac{\partial \phi}{\partial n}\bigg|_f \approx \frac{\phi_N - \phi_P}{|\mathbf{d}_{PN}|} + \text{correction}
$$

`corrected` は laplacian と同様の非直交補正を適用。

### 5.6 離散方程式の最終形

各セル $P$ について、離散化後の方程式は線形代数系:

$$
a_P \phi_P = \sum_N a_N \phi_N + b_P
$$

ここで:
- $a_P$: 中心セル係数（時間項 + 対流 + 拡散の寄与）
- $a_N$: 隣接セル係数（対流 + 拡散の寄与）
- $b_P$: ソース項（重力、熱源など）

全セルで組み立てると疎行列方程式 $\mathbf{A}\boldsymbol{\phi} = \mathbf{b}$ になる。

---

## 6. 離散化スキームの詳細（fvSchemes 対応表）

| fvSchemes エントリ | 数学的操作 | スキーム | 精度 | 安定性 | 根拠 |
|---|---|---|---|---|---|
| `ddt: backward` | $\partial/\partial t$ | 3点後退差分 | 2次 | 無条件安定 | 非定常自然対流の振動を捕捉 |
| `grad: Gauss linear` | $\nabla \phi$ | Green-Gauss + 線形補間 | 2次 | — | 標準的な勾配評価法 |
| `grad(U): cellLimited Gauss linear 1` | $\nabla \mathbf{U}$ | 制限付き Green-Gauss | 2次（制限時低下） | 有界 | linearUpwind の勾配。非物理的振動防止 |
| `div(phi,U): bounded Gauss linearUpwind grad(U)` | $\nabla \cdot (\rho\mathbf{U}\mathbf{U})$ | 2次風上 | 2次 | 有界 | 速度場の数値拡散低減。自然対流では対流項が支配的 |
| `div(phi,T): bounded Gauss linearUpwind default` | $\nabla \cdot (\rho\mathbf{U}T)$ | 2次風上 | 2次 | 有界 | 温度場の数値拡散低減。鋭い温度勾配の解像に必要 |
| `div(phi,k): bounded Gauss upwind` | $\nabla \cdot (\rho\mathbf{U}k)$ | 1次風上 | 1次 | 安定 | $k > 0$ の物理的制約を維持。乱流量の安定性優先 |
| `div(phi,omega): bounded Gauss upwind` | $\nabla \cdot (\rho\mathbf{U}\omega)$ | 1次風上 | 1次 | 安定 | $\omega > 0$ の物理的制約を維持 |
| `laplacian: Gauss linear corrected` | $\nabla \cdot (\Gamma \nabla \phi)$ | Green-Gauss + 非直交補正 | 2次 | — | 非直交メッシュでも2次精度を維持 |
| `snGrad: corrected` | $\partial\phi/\partial n |_f$ | 非直交補正付き | 2次 | — | laplacian と整合 |

---

## 7. 圧力-速度連成アルゴリズム

### 7.1 PIMPLE アルゴリズム

PISO (Pressure Implicit with Splitting of Operators) と SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) を融合した非定常圧力-速度連成法。

- **SIMPLE の出典:** Patankar & Spalding (1972), "A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows", Int. J. Heat Mass Transfer 15:1787
- **PISO の出典:** Issa, R.I. (1986), "Solution of the implicitly discretised fluid flow equations by operator-splitting", J. Comput. Phys. 62:40-65
- **PIMPLE の出典:** OpenFOAM 実装。Jasak (1996) の doctoral thesis に基盤
- **分野:** 数値流体力学（圧力-速度連成）

**アルゴリズムの流れ（各時間ステップ）:**

```
for each time step n → n+1:
    for outer = 1 to nOuterCorrectors (=2):     ← SIMPLE 的外側反復
        ① 運動量予測: H(U) で U* を計算
           a_P U*_P = H(U) - ∇p^n_rgh
           ここで H(U) = -Σ_N a_N U_N + source

        ② エネルギー方程式を解く → T^{n+1}

        ③ 乱流方程式を解く → k, omega

        ④ 密度更新: rho = p M_w / (R T)

        for inner = 1 to nCorrectors (=1):       ← PISO 的内側補正
            ⑤ 圧力方程式を組み立て・求解:
               ∇·(1/a_P ∇p'_rgh) = ∇·(H(U)/a_P)
               p^{n+1}_rgh = p^n_rgh + p'_rgh

            ⑥ 速度補正:
               U^{n+1} = H(U)/a_P - (1/a_P)∇p^{n+1}_rgh

            ⑦ 質量フラックス補正:
               phi^{n+1} = phi* - (1/a_P)_f |S_f| (∇p'_rgh)_f
        end inner
    end outer

    ⑧ 時間ステップ調整: Δt = min(Δt_max, maxCo * Δx_min / |U|_max)
end
```

**なぜ PIMPLE か（SIMPLE でない理由）:**

密閉空間の自然対流 ($Ra > 10^9$) は本質的に非定常（振動解）。
定常ソルバー (buoyantSimpleFoam) では収束せず残差が振動する場合がある。
buoyantPimpleFoam で非定常計算し、時間平均で統計量を取得する方が物理的に正しい。

### 7.2 圧力方程式の導出

運動量方程式の離散形:

$$
a_P \mathbf{U}_P = \mathbf{H}(\mathbf{U}) - \nabla p_{\text{rgh}}
$$

$\mathbf{U}_P$ について解くと:

$$
\mathbf{U}_P = \frac{\mathbf{H}(\mathbf{U})}{a_P} - \frac{1}{a_P} \nabla p_{\text{rgh}}
$$

これを連続の式 $\nabla \cdot (\rho \mathbf{U}) = 0$ に代入:

$$
\nabla \cdot \left( \rho \frac{\mathbf{H}(\mathbf{U})}{a_P} \right) = \nabla \cdot \left( \frac{\rho}{a_P} \nabla p_{\text{rgh}} \right)
$$

→ これが**圧力のポアソン方程式**（楕円型 PDE）。GAMG で解く。

---

## 8. 乱流モデル SST k-omega の詳細

### 8.1 モデルの出典と位置づけ

- **出典:** Menter, F.R. (1994), "Two-equation eddy-viscosity turbulence models for engineering applications", AIAA Journal 32(8):1598-1605
- **分野:** 乱流工学
- **位置づけ:** RANS (Reynolds-Averaged Navier-Stokes) の2方程式モデル。
  壁面近傍で k-omega モデル、主流域で k-epsilon モデルに自動的にブレンドする

### 8.2 なぜ SST k-omega を選択したか

| モデル | 壁面近傍 | 主流域 | 自然対流 | 計算コスト |
|--------|---------|--------|---------|-----------|
| 標準 k-epsilon | 壁面関数必須（y+ > 30） | 良好 | 熱伝達率を過大予測 | 低 |
| 標準 k-omega | 直接解法可能（y+ ~ 1） | 入口条件に敏感 | 良好 | 中 |
| SST k-omega | 直接解法可能 | k-epsilon 相当 | 良好 | 中 |
| LES/DES | 壁面モデリング | 大規模渦を直接解く | 最良 | 極めて高 |

- **根拠:** Heschl et al. (2005) "CFD simulation of the thermal comfort in an office room heated by a radiator" — 室内自然対流で SST k-omega が k-epsilon より壁面熱伝達を正確に予測
- サウナは壁面からの自然対流が支配的 → 壁面熱伝達の精度が全体精度を左右

### 8.3 輸送方程式

**乱流運動エネルギー $k$:**

$$
\frac{\partial (\rho k)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} k)
= \nabla \cdot \left[ (\mu + \sigma_k \mu_t) \nabla k \right]
+ P_k - \beta^* \rho k \omega
$$

| 項 | 物理 | 数学的役割 |
|----|------|-----------|
| $\partial(\rho k)/\partial t$ | 乱流エネルギーの時間変化 | 非定常蓄積 |
| $\nabla \cdot (\rho \mathbf{U} k)$ | 平均流れによる乱流エネルギー輸送 | 対流 |
| $(\mu + \sigma_k \mu_t) \nabla k$ | 分子+乱流拡散 | 拡散 |
| $P_k = \mu_t S^2$ | 平均速度勾配からの生成 | ソース（+） |
| $\beta^* \rho k \omega$ | 乱流の散逸 | シンク（−） |

- $S = \sqrt{2 S_{ij} S_{ij}}$: ひずみ速度テンソルの大きさ
- $P_k$ は `min(P_k, 10 * beta* * rho * k * omega)` で制限（stagnation point anomaly 防止）

**比散逸率 $\omega$:**

$$
\frac{\partial (\rho \omega)}{\partial t}
+ \nabla \cdot (\rho \mathbf{U} \omega)
= \nabla \cdot \left[ (\mu + \sigma_\omega \mu_t) \nabla \omega \right]
+ \frac{\gamma}{\nu_t} P_k - \beta \rho \omega^2
+ 2(1 - F_1) \frac{\rho \sigma_{\omega 2}}{\omega} \nabla k \cdot \nabla \omega
$$

最後の項 (cross-diffusion) が k-epsilon と k-omega のブレンドを実現する。

### 8.4 モデル定数

| 定数 | k-omega 側 ($\phi_1$) | k-epsilon 側 ($\phi_2$) | 出典 |
|------|----------------------|------------------------|------|
| $\sigma_k$ | 0.85 | 1.0 | Menter (1994) |
| $\sigma_\omega$ | 0.5 | 0.856 | Menter (1994) |
| $\beta$ | 0.075 | 0.0828 | Wilcox (1988) / Launder-Sharma (1974) から変換 |
| $\gamma$ | $\beta/\beta^* - \sigma_\omega \kappa^2/\sqrt{\beta^*}$ | 同左 | 対数則との整合性条件から導出 |
| $\beta^*$ | 0.09 | 0.09 | 実験値 ($\beta^* = C_\mu$) |
| $a_1$ | 0.31 | — | Bradshaw の仮説 $-\overline{u'v'} = a_1 k$ |
| $\kappa$ | 0.41 | — | von Kármán 定数（壁面乱流の対数則） |

**ブレンド:** $\phi = F_1 \phi_1 + (1-F_1) \phi_2$

### 8.5 ブレンディング関数

$$
F_1 = \tanh(\arg_1^4)
$$

$$
\arg_1 = \min\left[ \max\left( \frac{\sqrt{k}}{\beta^* \omega d}, \frac{500\nu}{d^2 \omega} \right), \frac{4\rho\sigma_{\omega 2} k}{CD_{k\omega} d^2} \right]
$$

- $d$: 最近壁面からの距離
- $CD_{k\omega} = \max(2\rho\sigma_{\omega 2} \nabla k \cdot \nabla \omega / \omega, \; 10^{-10})$

壁面近傍 ($d \to 0$) で $F_1 \to 1$（k-omega）、遠方で $F_1 \to 0$（k-epsilon）。

**乱流粘性:**

$$
\mu_t = \frac{\rho a_1 k}{\max(a_1 \omega, \; S F_2)}
$$

$F_2$ は第2ブレンディング関数。分母の $S F_2$ が Bradshaw の仮説
（逆圧力勾配領域での乱流応力制限）を実現 → SST の "Shear Stress Transport" の名前の由来。

---

## 9. 壁面関数の理論

### 9.1 壁面境界層構造

壁面近傍の速度分布（無次元化）:

$$
u^+ = \frac{U}{u_\tau}, \quad y^+ = \frac{y u_\tau}{\nu}, \quad u_\tau = \sqrt{\tau_w / \rho}
$$

| 領域 | $y^+$ 範囲 | 速度則 |
|------|-----------|--------|
| 粘性底層 | $y^+ < 5$ | $u^+ = y^+$（線形） |
| バッファ層 | $5 < y^+ < 30$ | 遷移（理論式なし） |
| 対数則層 | $y^+ > 30$ | $u^+ = \frac{1}{\kappa} \ln(y^+) + B$ ($\kappa = 0.41$, $B = 5.2$) |

- **対数則の出典:** von Kármán (1930), "Mechanische Aehnlichkeit und Turbulenz", Proc. 3rd Int. Congr. Applied Mechanics
- **分野:** 乱流境界層理論

### 9.2 本計算で使用する壁面関数

| 変数 | 壁面関数 | 役割 | 理論根拠 |
|------|---------|------|---------|
| $k$ | `kqRWallFunction` | 壁面での $k$ に zero-gradient を適用 | $k$ は壁面近傍で有限値を持つ（壁面で $\partial k / \partial n \approx 0$） |
| $\omega$ | `omegaWallFunction` | 壁面セルで $\omega = 6\nu/(\beta_1 y^2)$ を強制 | Wilcox (1988) の壁面漸近解。粘性底層で $\omega \propto 1/y^2$ |
| $\nu_t$ | `nutkWallFunction` | $y^+$ に応じて $\nu_t$ を壁面関数値で設定 | 対数則層: $\nu_t = \kappa u_\tau y - \nu$、粘性底層: $\nu_t = 0$ |
| $\alpha_t$ | `compressible::alphatWallFunction` | $\alpha_t = \nu_t / Pr_t$ を壁面で設定 | Reynolds analogy の修正。$Pr_t = 0.85$ で温度の乱流輸送をモデル化 |

### 9.3 $y^+$ の影響と本計算の方針

SST k-omega は壁面関数なし ($y^+ \sim 1$) でも壁面関数あり ($y^+ \sim 30\text{--}50$) でも動作可能（automatic wall treatment）。

本計算の M0 メッシュ (2-10万セル) では壁面セルの $y^+$ が 5-50 程度になり、
壁面関数が自動的にブレンドする。M2 以上の高精細メッシュでは $y^+ \sim 1$ を目標に
壁面近傍を細分化し、壁面関数への依存を低減する。

---

## 10. 線形ソルバーと緩和係数の根拠

### 10.1 線形ソルバー選択

| 変数 | ソルバー | 根拠 |
|------|---------|------|
| $p_{\text{rgh}}$ | GAMG (Geometric-Algebraic Multi-Grid) | 圧力方程式は楕円型 → マルチグリッドが最適。出典: Wesseling (1992), "An Introduction to Multigrid Methods" |
| $U, T, k, \omega$ | PBiCGStab (Preconditioned Bi-Conjugate Gradient Stabilized) | 非対称行列（対流項が非対称性を導入）に適した Krylov 部分空間法。出典: van der Vorst (1992), "Bi-CGSTAB", SIAM J. Sci. Stat. Comput. 13:631 |

**前処理:**
- GAMG: Gauss-Seidel スムーザ — マルチグリッドの標準的スムーザ
- PBiCGStab: DILU (Diagonal Incomplete LU) — 非対称行列用の不完全分解前処理

### 10.2 緩和係数

$$
\phi^{k+1} = \phi^k + \alpha (\phi^{k+1}_{\text{solved}} - \phi^k)
$$

$\alpha < 1$ で Under-relaxation（安定化）、$\alpha = 1$ で無緩和。

| 変数 | $\alpha$ | 根拠 |
|------|----------|------|
| $p_{\text{rgh}}$ | 0.3 | 圧力は運動量と強く結合し不安定になりやすい。SIMPLE 法では $\alpha_p = 0.2\text{--}0.3$ が標準。出典: Patankar (1980) |
| $U$ | 0.7 | 速度は圧力ほど敏感でない。$\alpha_U = 0.5\text{--}0.8$ が標準 |
| $T$ | 0.5 | 温度は密度（→浮力→速度→圧力）と連鎖結合する。自然対流では $\alpha_T = 0.3\text{--}0.7$ が推奨 |
| $k, \omega$ | 0.7 | 乱流方程式は比較的安定。$0.5\text{--}0.8$ が標準 |

**注:** PIMPLE（非定常）では外側反復の最終パスで `Final` フィールドの `relTol 0` が適用され、
各時間ステップの最終解は相対許容残差なしで求解される。

### 10.3 CFL 数と適応的時間ステップ

$$
Co = \frac{|\mathbf{U}| \Delta t}{\Delta x}
$$

- $Co_{\max} = 0.5$: backward スキーム（陰的）は理論上無条件安定だが、
  精度確保のため $Co < 1$ が推奨。自然対流の振動的解では $Co = 0.5$ が安全
- `adjustTimeStep yes`: 速度場に応じて $\Delta t$ を動的調整
- $\Delta t_{\text{init}} = 0.05$ s: 初期の安全な時間刻み
- $\Delta t_{\text{max}} = 0.05$ s: 上限（YAML の `delta_t`）

---

## 付録 A: 方程式番号と実装ファイルの対応表

| 方程式 | 理論出典 | OpenFOAM 実装 | 簡易版コード |
|--------|---------|---------------|-------------|
| 連続の式 | Euler (1757) | buoyantPimpleFoam 内部 | — (0D, 界面質量保存で代替) |
| N-S | Navier-Stokes (1822/1845) | buoyantPimpleFoam 内部 | — (速度を解かない) |
| エネルギー保存 | 第一法則 | buoyantPimpleFoam 内部 | `simple_solver.py` `solve_two_zone()` L707-776 |
| 状態方程式 | 理想気体 | `thermophysicalProperties` | `simple_solver.py` L683-684 |
| SST k-omega | Menter (1994) | `turbulenceProperties` | — (相関式に内包) |
| Zukoski プルーム | Zukoski (1978) | — | `simple_solver.py` `_plume_entrainment()` L164-209 |
| Newton 冷却 | Newton (1701) | — (壁面関数+BC) | `simple_solver.py` L714, L754 |
| 界面質量保存 | Cooper (1982) | — | `simple_solver.py` L778-789 |
| 壁面 lumped | Incropera Ch.5 | — | `simple_solver.py` L796-806 |
| 換気 (stack effect) | ASHRAE 2017 | — | `simple_solver.py` `_ventilation_flow()` L282-355 |
| 体感温度 (皮膚熱収支) | ISO 7933 | — | `simple_solver.py` `_perceived_temperature()` L227-279 |
| 人体輻射 | Stefan-Boltzmann | — | `simple_solver.py` `_q_rad_body()` L415-456 |
| View Factor | Hottel & Sarofim | — | `simple_solver.py` `_compute_view_factors()` L76-161 |

## 付録 B: 参考文献一覧

1. **Morton, Taylor & Turner (1956)** — "Turbulent Gravitational Convection from Maintained and Instantaneous Sources", Proc. R. Soc. A 234:1-23. *エントレインメント仮説の原論文*
2. **Zukoski (1978)** — "Development of a Stratified Ceiling Layer in the Early Stages of a Closed-Room Fire", NBS-GCR-78-150. *プルーム質量流量相関式*
3. **Cooper (1982)** — "A Mathematical Model for Estimating Available Safe Egress Time in Fires", NBSIR 82-2612. *2層ゾーンモデルの体系化*
4. **Patankar & Spalding (1972)** — SIMPLE アルゴリズムの原論文
5. **Patankar (1980)** — "Numerical Heat Transfer and Fluid Flow", Hemisphere Publishing. *FVM・SIMPLE の教科書*
6. **Issa (1986)** — PISO アルゴリズムの原論文
7. **Menter (1994)** — SST k-omega モデルの原論文
8. **Jasak (1996)** — PhD thesis, Imperial College. *OpenFOAM の FVM 実装の理論的基盤*
9. **Ferziger & Perić (2002)** — "Computational Methods for Fluid Dynamics", 3rd ed. *CFD 離散化の教科書*
10. **Churchill & Chu (1975)** — 鉛直平板自然対流の相関式
11. **Incropera et al.** — "Fundamentals of Heat and Mass Transfer", 7th ed. *空気物性値の出典*
12. **Bird, Stewart & Lightfoot** — "Transport Phenomena", 2nd ed. *輸送方程式の導出*
13. **Drysdale** — "An Introduction to Fire Dynamics", 3rd ed. *対流熱割合 $f_{\text{conv}}$*
14. **Wilcox (1988)** — k-omega モデルの原論文. *壁面境界条件の漸近解*
15. **von Kármán (1930)** — 対数則の原論文
16. **Warming & Beam (1976)** — linearUpwind スキームの理論
17. **Barth & Jespersen (1989)** — セル制限付き勾配の理論
18. **van der Vorst (1992)** — BiCGStab 法の原論文
19. **Wesseling (1992)** — マルチグリッド法の教科書
20. **Kays & Crawford** — "Convective Heat and Mass Transfer", 4th ed. *乱流プラントル数*
