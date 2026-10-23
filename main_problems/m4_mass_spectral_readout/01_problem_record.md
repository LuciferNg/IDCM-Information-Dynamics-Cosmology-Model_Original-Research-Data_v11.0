# M4 — 質量 = 譜讀數差值：全費米子譜差閉合（只開不推）

<!-- registry
pid:     M4
title:   質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2, R7 定案承接；只開不推）
state:   partial
verdict: -
flags:   running; r1-gate-pass; r2-bridge-exact; r3-lattice-rejudge; r4-negative
owner:   M4
created: 2026-10-21
updated: 2026-10-24
prereq:  O257 (R7 docs/07 定案) + O66 (M_comp 降級) + O89 (γ₁ 冪階梯, O85 污染通報後更新版) + O57 (m_ν 盲測) + P6 (§10 差值閉合) + O73 (c=Goldstone) + O212 (dof=1)
source:  調律者 2026-10-21: 開v11 新M, 將/tmp全部py移入scripts；新session重建推導鏈及報告文檔
-->

> **狀態**: 📄 **已註冊 — 未開始（只開不推；推導鏈＋報告文檔由新 session 完成）**
> **2026-10-24 追加**: 🟡 **partial — R1 閘通過、R2 橋精確、R3 lattice 重判、R4 A 攻堅（negative）；未 close（見 §八／docs/06）**
> **建立**: 2026-10-21（v11 Reality-Architecture）
> **位置**: `main_problems/m4_mass_spectral_readout/`
> **來源**: O257 R7 定案（`open_problems_fresh/o257_mcomp_log10_mechanism/docs/07_round7_difference_closure.md`）＋ /tmp 8 支測試 script 全量遷移

---

## 〇、歸屬審計（fence）

| 攻擊項 | 已有主 | 裁定 |
|:--|:--|:--|
| 中微子 Δm² 差值閉合 | O257 R7 🟢 已定案 | 本 M 只承接（唔重推） |
| M_comp 絕對標度 | O66 已降級（量級 ±3%）＋O257 negative | 🔴 唔重推 |
| GUT 標度 γ₁ 冪階梯 | O89（O85 污染通報後更新版） | 引用需核證；本 M 做 cross-check |
| 夸克/輕子/hadron 譜差掃描 | /tmp 8 scripts（本 M 唯一新產物搬入） | 本 M 驗證＋形式化 |

## 一、問題陳述

**核心問題（調律者 framing）**：質量 = 駐波態之間交互嘅**同步邊界譜讀數差值**。
R7 喺中微子 sector 確立「差值閉合力學」：實驗差譜（Δm²）+ 譜比 R 直接出質量譜，M_comp 係多餘中間人；m₃ = γ₆²·A（讀數×譜強度）。

**本 M 任務**：
1. 重建推導鏈（8 支遷入 scripts 由新 session 逐支驗證）
2. 全費米子譜差結構：quark/lepton/hadron 掃描已做（全部負向/撞池），需形式化記錄
3. GUT 標度錨（γ₁ 冪階梯）可否加入差值閉合——O89 已有 M_GUT/PS/R/EW = M_Pl×γ₁ᵏ；γ₁↔γ₆ 結構橋未解
4. 推導鏈重建後出報告文檔

## 二、已確立嘅結構（遷入 scripts 結論快照）

### 2.1 中微子（O257 R7 定案，🟢）
- Δm²(7.53e-5, 2.434e-3) + R(≈γ30/γ17) ⇒ m₁=0.020892/m₂=0.022623/m₃=0.053577/Σm=0.09709 eV（dev 0.02-0.08%）
- **m₃ = γ₆²·A**（γ₆²=0.0326662, dev 0.002%）；A = v²/M_comp×1e9 = 1.64044 eV/γ²（A 身份 = 唯一未閉合標度）
- 「質量比=γ譜比」直路：🔴 死（MC null A∪B 池 0.1% 級 63.4%）

### 2.2 夸克（scripts/01-03，🔴 全負向）
- m_c/m_s vs γ20/γ5 dev 0.042%——MC null **97.9%** 撞池（γ比池密度=隨機級）→ 後驗，唔升格
- 讀數×譜強度 m_q=γ²·A：全部 >6%；m_q=γ·v：m_c/v≈γ1(1.29%)/m_t/v≈γ12(1.67%)，事後揀 mode 唔過閘
- u/d/s 不確定度 10-20%：0.04% 級「閉合」係 fit 入測量誤差，零預言力

### 2.3 帶電輕子（scripts/04，🔴 全負向）
- μ/e τ/μ 直接 γ比 dev 1.04-1.15%、τ/e 77%；MC null 16.8%（輕子池細過夸克但仍 >5%）
- memory μ/e=(γ₇/γ₁)(F_μ/F_e)=206.705(0.031%)：F_μ/F_e≈4.3069 係**後驗 permissive 因子**（複核坐實），唔係 EXACT 結構

### 2.4 hadron 質量差（scripts/05，🔴 負向＋🟢 物理 insight）
- π⁺−π⁰=+4.5936 MeV（em probe) vs K⁺−K⁰=−3.937 MeV（strong probe）**反相** = 「雙 probe 讀數差」實驗事實（framing 支持）
- Δm 比 γ池：全 >0.01%；n−p/m_e=γ29/γ16 dev 0.0195%（γ29=L₇/γ16=B² 身份靚）但 MC null 8.7%（0.01% 級）→ density 殺
- d−u 不確定度 ±10%：anchor 唔夠格

### 2.5 hop 階梯 + 譜跨度（scripts/06-08，🟡 已判）
- s/d=20.0 EXACT（hop 池零偏差）——**但 m_d ±10%** → 誤差帶 [17.4,22.6]，PDG 引用值巧合
- quark/neutrino log-span 比 = **11.994978 ≈ 12 = hop_max = 88−L₉**（dev 0.042%，MC null 0.166% 低）——但 m_u ±10% → R∈[11.89,12.11]，誤差帶 kill
- 微子 log(m₂/m₁)=0.0347 vs hop 11→12 微間距 0.0378：dev 8.3% 唔過

### 2.6 GUT 標度錨（O89 跨查＋實測，🟡 結構實例成立）
- **γ₁ 冪階梯**：M_GUT/PS/R/EW = M_Pl×γ₁ᵏ×{1,1,1,N_d}（6.22e16/3.17e14/1.62e12/244.24 GeV）
- M_EW = M_Pl·γ₁⁸·N_d = 244.24 vs v=246.22：dev **0.80%**（O89 已錄）
- **「譜讀數×譜強度」統一形態第二實例**：γ₁（gauge 邊界）同 γ₆（fermion 邊界）各自內部閉合
- ❌ γ₁ᵏ×F 對 9 費米子零命中 <2%；M_comp log₁₀ k=2.407 非整數——γ₁↔γ₆ 冇直接耦合
- ⚠️ O89 有 O85 污染通報（舊 1e16/1e14 硬編碼）；現 scripts 已更新為 M_Pl×γ₁ᵏ，引用需 RAG 核證

## 三、未確立（🔴 開放問題）

1. **A 嘅譜身份**（A=v²/M_comp×1e9=1.64044 eV/γ²）——唯一未閉合標度，鏈返 M_comp
2. **γ₁↔γ₆ 結構橋**（「同步邊界之間嘅邊界」）——有冇 γ₆/γ₁=35.47≈? 或混合形式
3. **全費米子譜差結構** 嘅正確形式（單點質量比已全滅；差譜 anchor 先行原則）
4. **quark span≈12=hop_max** 係真結構定巧合（需 lattice 高精度 m_u/m_d/m_s 重測——FLAG 值代入重做 MC null）
5. **新實驗 anchor 落地後重測**（KATRIN/LEGEND m_ββ、DUNE、lattice QCD 1% 級）

## 四、攻擊路線（新 session，僅列）

- R1: 逐支重跑 8 支遷入 scripts 重建推導鏈（失敗即停）
- R2: γ₁↔γ₆ 橋掃描（γ₆/γ₁、γ₁^kγ₆^l·A 混合形式，density 閘）
- R3: lattice 值（FLAG 2024 m_u/m_d/m_s）代入重做 span/比掃描
- R4: A 譜身份攻堅（v²/M_comp→純譜量重寫）
- R5: 報告文檔（docs/01-03：定案承接／掃描全滅／GUT 錨分層）

## 五、驗證標準

- 每支 script 重跑 PASS 數一致（重建閘）
- 任何新閉合過四硬測（零參數／盲測／樣本外／MC null density 閘）
- u/d/s 誤差帶內命中一律標後驗，唔升格
- 引 O89 γ₁ 階梯前先 RAG 核證 O85 通報後狀態

## 六、狀態

- 📄 **已註冊 — 只開不推**（推導鏈＋報告文檔 = 新 session 第一步）
- **2026-10-24**：🟡 partial — R1 閘 8/8、R2 結構橋（Chebyshev 精確）、R3 lattice 重判（span≈12 升 🟡）、R4 A（🔴 negative）；詳 §八＋`docs/01–06`

## 七、相關文件

- `scripts/01_quark_scan.py` … `scripts/08_span_12_judge.py`（8 支遷入）
- `scripts/09_density_gate.py` … `scripts/12_A_spectral_identity.py`（本輪新增）
- O257 docs/07（R7 定案，上游）
- O89 docs/01_gamma_cascade_closure.md（γ₁ 階梯，O85 通報後更新版）
- O57 docs/16（2.96% 先例，誠實負面）
- O66 docs/02（M_comp 1/76 循環，已降級）

---

## 八、R1–R4 執行結果（2026-10-24 追加；原文 §〇–§七 一字未改）

**推導時序索引**：R1 閘（04:56）→ F1 修 01 [C]（04:58）→ F2 05 MC 判別（05:00）→ R3 lattice（05:03，含 FLAG PDF 實檔）→ R2 橋（05:06）→ R4 A（05:07）→ 文檔 docs/01–06（05:15）

| 輪 | 結果 | 檔 |
|:--|:--|:--|
| R1 閘 | 8/8 跑通；7/8 逐位吻合；01 首跑 crash（TypeError，/tmp 原版 bug）已修；05 MC 43.3% ≠ 記錄 8.7%（不可復現） | docs/01 |
| R1-corr | 統一 density protocol（thr = 命中自身 dev；域 = 受測比值域）：**全部 sector PN ≥ 31%** ⟹ 夸克/輕子/hadron 命中全滅 | docs/02 |
| R2 | γ_n = 2−2T_n(c) 精確；遞推 `γ_{n+1}=2cγ_n−γ_{n−1}+2γ₁`；γ₆/γ₁ = deg-20 本原元（ℚ(ζ₈₈)⁺）⟹ **「橋」唔缺常數** | docs/03 |
| R3 | FLAG 2024 lattice：c/s 🔴（PN 78.75%）、s/d=20 🔴（19.8866±0.2452，PN 99.69%） | docs/04 |
| R3 ★ | **span_Q/span_ν = 12 由 🔴 升 🟡**：R=12.0055±0.0399（0.14σ）；預言 **m_u = 2.1498 MeV**（FLAG 2.14±0.08，+0.12σ） | docs/04 |
| R3 🟡 | 新命中 **m_d/m_u = L₅/F₅ = 2.19929±0.08605**（dev 0.032%, 0.01σ；密度 P1 1.85%；後驗） | docs/04 |
| R4 | A = 1.64337 eV/γ²（v_SM）；族 13104 組合 **0 命中 <0.1%**，密度 5.21% ⟹ 🔴 未閉合 | docs/05 |
| R2 🔴 | O89 M_GUT/PS/R 三階梯 **不可引用**（O85 由零重建證階梯唔存在、目標硬編碼） | docs/03 |

**狀態**：M4 = partial（未 close）。剩 🟡 兩項：span≈12（等 m_u <0.5% 精度 ＋ sample-out blind）、m_d/m_u=11/5（等第二實例）。

---

## 九、§2.4 MC 數字更正（append-only 附註，原文不動）

- §2.4 寫「n−p/m_e=γ29/γ16 dev 0.0195% … 但 MC null **8.7%**（0.01% 級）→ density 殺」
- **2026-10-24 更正附註**：遷入腳本 `scripts/05_hadron_dm.py` 實跑 MC null = **43.3%**（域 0.1–10、thr 0.1%）；
  8.7% 對應一個**唔同 protocol**（域 1–10、thr 0.015%），唔可以由現存腳本唯一復現（詳 docs/01 §三）。
  統一 protocol 重判（docs/02）：同域同 N、thr = 命中自身 dev（0.0195%）⟹ **PN = 30.98%**。
  **裁定不變**：density 殺（30.98% > 5%），但數字以 **30.98%** 為準（8.7% 作廢為不可復現）。

- §2.2 寫「m_c/m_s vs γ20/γ5 dev 0.042%——MC null **97.9%**」：97.9% 係 02 腳本嘅**族級**null（質量域）；
  dev 0.042% 嘅**單比/同域族級** PN = **48.07%**（docs/02 §二）。兩者都 ⟹ 隨機級，數字以 docs/02 為準。