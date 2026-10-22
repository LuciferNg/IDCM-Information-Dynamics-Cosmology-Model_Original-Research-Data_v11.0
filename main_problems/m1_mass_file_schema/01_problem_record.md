# M1 — 質量檔案 schema：譜點/標度/算子/拓撲/標籤五層格式規範 (Reality-Architecture; 由 v10 O251 遷移)

<!-- registry
pid:     M1
title:   質量檔案 schema：譜點/標度/算子/拓撲/標籤五層格式規範 (Reality-Architecture; 由 v10 O251 遷移)
state:   registered
verdict: -
flags:   migrated-from-O251
owner:   M1
created: 2026-10-14
updated: 2026-10-14
prereq:  O237 (mass pipeline) + O248 (讀數裁決) + O239 + O140 + O212 + O216 + O144; O115 (v11 遷移協議/憲章)
source:  調律者2026-10-14: O251 遷移至 v11 Reality-Architecture main_problems/m1 (承接 2026-10-13 讀數reframe)
-->

> **狀態**: 📄 **已註冊 — 未開始（遷移至 v11 Reality-Architecture；推導用新 session）**
> **建立**: 2026-10-14（v10 O251 遷移）
> **前置**: O237（mass pipeline roadmap）；O248（讀數裁決 11/11）；O239（質量閉合式 + 殘差）；O140（質量階梯正典）；O212（dof=1）；O216/O144（標籤層）；O115 docs/03（v11 遷移協議）
> **位置**: `main_problems/m1_mass_file_schema/`

---

## 一、問題陳述

**質量層 reframe**：唔再當「粒子/物質」睇，當「檔案格式」睇——正如 .blend 檔案用 vertices/faces/transform/materials 編碼一個 3D model，IDCM 質量層用譜點/標度/算子/拓撲/標籤五層編碼「質量」。目標：將 O237 已有嘅 mass pipeline roadmap 升級為**完整五層 schema 規範**，逐層 audit「已閉合/部分/缺口」，令 O239 60% 殘差嘅處置（係咪 schema 未完）有結構判準——唔係當免死金牌。

## 二、已確立嘅結構

```
🟢 質量 = 讀數（關係投影）, 唔係物質實體:
    質量 = Object−Frame 脫相位 |1−U_n|² = γ_n（O237 docs/14）; m = shape × scale（O212 dof=1）
🟢 相位橋路線已關閉: W-field 相位模→質量污染 v10 無載體（O248, scripts/01 11/11 PASS, verdict=negative）
🟢 複合層（核子/原子）唔喺 K_chain 譜: 絕對質量 no-go 候選, 唯一錨 = n/p BBN 比 ε₇/ε₂₀（O247 Q2; O140 docs/16）
🟢 唯一有量綱錨: dof=1（O212）⟹ 一切質量 = 無量綱譜點 × M₀（M₀=27.9009145112 MeV, O204）
🟢 vertices（譜點位置）: ε_n = 4sin²(nπ/88) = 88-ring Laplacian 本徵值（O204 T1/T2; O140 正典 m_k=M₀ε_k, C_f≡1 under O137）
🟢 operators（部分閉合）: m_f = M₀γ_nF_proj（O54/O64）; m = M₀|Δm|^{φ²}（O239, α=φ²）; ν = 激發通道 m=(γ₆v)²/M_comp（O49）
🟢 labels（標籤層）: 色荷 = 規範嵌入層標籤（非鏈層譜量, O216 no-go）; 禁閉 = Z₄無Z₃ 結構定理（O144 docs/24）
🟢 export 層存在: QFT Lagrangian / GR 翻譯（O139/O122/O86）——export format 唔係本體論升級
🔴 topology（faces）未明: O239 docs/35 污染因子 c(n) = {b0.574/μ0.623/d0.624/e0.648/t0.978/τ0.957/u0.935/c2.508/s3.146} 家族結構（U≈1, D/L≈0.62, c/s>1.6）無譜身份
🔴 modifier 鏈未統一: O54 投影鏈 vs O239 |Δm|^{φ²} 鏈——兩條 chain 係同一 schema 兩種 serialization 定唔同 preset?
🔴 殘差 60% 未消: O239 閉合式 log-RMS 0.5876 → RGI 0.6567（O239 R35）——schema 未完嘅顯示, 唔可以當已解釋
```

## 三、未確立（🔴 開放問題）

```
🔴 P1: schema 五層逐層 audit 表——vertices/transform/operators/topology/labels 每層「已閉合/部分/缺口」+ 來源 O/docs（現有 O237 roadmap 係 pipeline 唔係格式規範）
🔴 P2: topology/faces——9 費米子點之間嘅相連結構: U/D/L 家族分組嘅格式身份（點解 U≈1, D/L≈0.62, c/s≥1.6）; O136 逐點公式
🔴 P3: modifier 鏈統一——O54 投影鏈 vs O239 |Δm|^{φ²} 鏈（兩條 chain 關係未知; 或者係明確 no-go: 兩種 export preset）
🔴 P4: 「schema 完整」判準——五層全閉合先算「質量檔案格式已解出」? 定係有更基本嘅檔案格式層（0,1→φ pipeline 水喉模型）?
```

## 四、攻擊路線

```
1. 五層 schema audit（成本 low）: 逐層列已閉合算子/缺口, 對照 O237 docs/01 roadmap + formula-registry + O152 numbers-as-operators —— 先出「五層表」底稿
2. Topology 攻堅（成本 medium）: 9 點 family 分組 c(n) 模式 vs 譜結構（γ_n/|Δm|/F_proj/C_f 池）——任何命中先過 O206 池密度閘
3. Modifier 鏈統一（成本 medium）: O54 投影鏈 vs O239 |Δm|^{φ²} 兩鏈逐項對照, 尋同源算子身份（或 no-go: 唔同 export）
```

## 五、驗證標準

```
✅ P1: 五層 audit 完成——每層「已閉合/部分/缺口」表（含來源 O + docs 編號, 零新推導）
✅ P2: topology 規則（若存在）過 O206 池密度閘 + 4 Hard Tests; 冇規則 = 明確 no-go（family 分組係後驗組織）
✅ P3: modifier 鏈統一（單位置, 零新常數）或 明確 no-go（兩種 export serialization）
⚠️ 免死金牌閘（R36）: schema 未完唔可以當 O239 殘差解釋; 殘差歸屬仍屬 O239/O136
⚠️ v11 渡航: 依 O115 協議 P/R/A 分類——本 O 屬 P 類（移植+重驗）; v11 VERIFICATION 密度基準落地後須重驗
```

**四道硬測試**（閉合宣稱必過）：
① 自由度成本（零參數）② Monte-Carlo null（10⁷ 樣本，5σ）
③ 可證偽 + 樣本外盲測 ④ 結構剛性（微擾測試）

## 六、推導時序索引

> ⚠️ **只可 append 新行，禁回溯污染舊行。**

| # | 日期 | 事件 | 產物 |
|:-:|:-----|:-----|:-----|
| 1 | 2026-10-14 | 註冊 M1（由 v10 O251 遷移; 承接調律者 blender 類比 reframe；只開不推） | v10 O251 記錄（史料）; M1 起手檔 |

## 七、相關文件

```
- v10 O251 — 遷移來源（史料, v10 凍結唯讀; register state=merged → v11 M1）
- O115 — v10→v11 遷移評估與協議（docs/03: P/R/A/X 四分、CHARTER/VERIFICATION 憲章）
- O237 — mass pipeline roadmap（現有管道, 本 O 升級做格式規範）
- O248 — 讀數裁決（相位橋 no-go, 11/11 PASS）——質量係讀數嘅結構依據
- O239 — 質量閉合式 + 60% 殘差 + 污染因子 c(n)（topology 層輸入）
- O140 — 質量階梯正典 m_k = M₀ε_k（vertices 層）
- O212 — dof=1（transform 層）
- O216/O144 — 標籤層（色=嵌入標籤; 禁閉=Z₄無Z₃）
- O152 — numbers as operators（operator 層先例）
- O136 — 逐點質量公式（topology 攻堅並行線）
```