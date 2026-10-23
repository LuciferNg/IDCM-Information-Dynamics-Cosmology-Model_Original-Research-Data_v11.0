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
> **2026-10-24 追加**: 🟡 **partial — R1 閘通過、R2 橋精確、R3 lattice 重判、R4 A 攻堅（negative）；未 close（見 §八／docs/13_script13_mdmu_second_instance.md）**
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
- R5: 報告文檔（docs/01–08-03：定案承接／掃描全滅／GUT 錨分層）

## 五、驗證標準

- 每支 script 重跑 PASS 數一致（重建閘）
- 任何新閉合過四硬測（零參數／盲測／樣本外／MC null density 閘）
- u/d/s 誤差帶內命中一律標後驗，唔升格
- 引 O89 γ₁ 階梯前先 RAG 核證 O85 通報後狀態

## 六、狀態

- 📄 **已註冊 — 只開不推**（推導鏈＋報告文檔 = 新 session 第一步）
- **2026-10-24**：🟡 partial — R1 閘 8/8、R2 結構橋（Chebyshev 精確）、R3 lattice 重判（span≈12 升 🟡）、R4 A（🔴 negative）；詳 §八＋`docs/01–13`

## 七、相關文件

- `scripts/01_quark_scan.py` … `scripts/08_span_12_judge.py`（8 支遷入）
- `scripts/09_density_gate.py` … `scripts/12_A_spectral_identity.py`（本輪新增）
- O257 docs/07（R7 定案，上游）
- O89 docs/01_gamma_cascade_closure.md（γ₁ 階梯，O85 通報後更新版）
- O57 docs/16（2.96% 先例，誠實負面）
- O66 docs/02（M_comp 1/76 循環，已降級）

---

## 八、R1–R4 執行結果（2026-10-24 追加；原文 §〇–§七 一字未改）

**推導時序索引**：R1 閘（04:56）→ F1 修 01 [C]（04:58）→ F2 05 MC 判別（05:00）→ R3 lattice（05:03，含 FLAG PDF 實檔）→ R2 橋（05:06）→ R4 A（05:07）→ 文檔 docs/01–13（05:15）

| 輪 | 結果 | 檔 |
|:--|:--|:--|
| R1 閘 | 8/8 跑通；7/8 逐位吻合；01 首跑 crash（TypeError，/tmp 原版 bug）已修；05 MC 43.3% ≠ 記錄 8.7%（不可復現） | docs/01–08 |
| R1-corr | 統一 density protocol（thr = 命中自身 dev；域 = 受測比值域）：**全部 sector PN ≥ 31%** ⟹ 夸克/輕子/hadron 命中全滅 | docs/09_script09_density_gate.md |
| R2 | γ_n = 2−2T_n(c) 精確；遞推 `γ_{n+1}=2cγ_n−γ_{n−1}+2γ₁`；γ₆/γ₁ = deg-20 本原元（ℚ(ζ₈₈)⁺）⟹ **「橋」唔缺常數** | docs/10_script10_gamma1_gamma6_bridge.md |
| R3 | FLAG 2024 lattice：c/s 🔴（PN 78.75%）、s/d=20 🔴（19.8866±0.2452，PN 99.69%） | docs/11_script11_lattice_flag_span.md |
| R3 ★ | **span_Q/span_ν = 12 由 🔴 升 🟡**：R=12.0055±0.0399（0.14σ）；預言 **m_u = 2.1498 MeV**（FLAG 2.14±0.08，+0.12σ） | docs/11_script11_lattice_flag_span.md |
| R3 🟡 | 新命中 **m_d/m_u = L₅/F₅ = 2.19929±0.08605**（dev 0.032%, 0.01σ；密度 P1 1.85%；後驗）— **同日 R3b 已降級，見下行** | docs/11_script11_lattice_flag_span.md |
| R3b 🔴 | m_d/m_u=11/5 **第二實例搜索 = 降級**：K-ratio 族 79 值，密度 P1(0.17%) = 5.69% ⟹ 15 比期望命中 0.85；觀察 2 個 <1% 少於隨機期望 4.34 ⟹ 密度解釋 | docs/11_script11_lattice_flag_span.md §四 |
| R4 | A = 1.64337 eV/γ²（v_SM）；族 13104 組合 **0 命中 <0.1%**，密度 5.21% ⟹ 🔴 未閉合 | docs/12_script12_A_spectral_identity.md |
| R2 🔴 | O89 M_GUT/PS/R 三階梯 **不可引用**（O85 由零重建證階梯唔存在、目標硬編碼） | docs/10_script10_gamma1_gamma6_bridge.md |

**狀態**：M4 = partial（未 close）。剩 🟡 兩項：span≈12（等 m_u <0.5% 精度 ＋ sample-out blind）、m_d/m_u=11/5（等第二實例）。

---

## 九、§2.4 MC 數字更正（append-only 附註，原文不動）

- §2.4 寫「n−p/m_e=γ29/γ16 dev 0.0195% … 但 MC null **8.7%**（0.01% 級）→ density 殺」
- **2026-10-24 更正附註**：遷入腳本 `scripts/05_hadron_dm.py` 實跑 MC null = **43.3%**（域 0.1–10、thr 0.1%）；
  8.7% 對應一個**唔同 protocol**（域 1–10、thr 0.015%），唔可以由現存腳本唯一復現（詳 docs/05 §七）。
  統一 protocol 重判（docs/09_script09_density_gate.md）：同域同 N、thr = 命中自身 dev（0.0195%）⟹ **PN = 30.98%**。
  **裁定不變**：density 殺（30.98% > 5%），但數字以 **30.98%** 為準（8.7% 作廢為不可復現）。

- §2.2 寫「m_c/m_s vs γ20/γ5 dev 0.042%——MC null **97.9%**」：97.9% 係 02 腳本嘅**族級**null（質量域）；
  dev 0.042% 嘅**單比/同域族級** PN = **48.07%**（docs/09_script09_density_gate.md §二）。兩者都 ⟹ 隨機級，數字以 docs/09_script09_density_gate.md 為準。

---

## 十、逐腳本文檔登記（2026-10-24 追加；編號規則 = docs/NN ↔ scripts/NN）

**指令**：所有 Python 一律「假設從未執行過」重跑；每個腳本必須有 markdown 記錄**前設問題 ＋ 結論報告**。
**執行**：13 支全部 fresh-run（`PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3`），exit 0；canonical log ＋ sha256 記錄如下。

| 腳本 | 文檔 | canonical log（sha256 前 12） | 執行 | 結論 |
|:--|:--|:--|:--|:--|
| `01_quark_scan.py` | `docs/01` | `run_01_quark_scan.log` `c7e8b76ec2df` | 102 ms / 76 行 | 🔴 單點譜比全滅（池密度 44%/96%）|
| `02_mc_mcs_judge.py` | `docs/02` | `run_02_mc_mcs_judge.log` `f3d2e2705ec5` | 652 ms / 50 行 | 🔴 97.9% 隨機級 |
| `03_reading_strength.py` | `docs/03` | `run_03_reading_strength.log` `1c7fe47448ce` | 82 ms / 40 行 | 🔴 γ²·A 形式對夸克失敗 |
| `04_lepton_anchor.py` | `docs/04` | `run_04_lepton_anchor.log` `9633b6c19908` | 136 ms / 35 行 | 🔴 sector 因子凑成；✅ anchor 分層 |
| `05_hadron_dm.py` | `docs/05` | `run_05_hadron_dm.log` `86cc1139fb5f` | 116 ms / 36 行 | 🔴 MC 43.3%；🟡 em/strong 反相 |
| `06_hop_probe.py` | `docs/06` | `run_06_hop_probe.log` `772b5ed26c8d` | 91 ms / 38 行 | 🔴 hop 非質量算子 |
| `07_hit_followup.py` | `docs/07` | `run_07_hit_followup.log` `37941c2647ab` | 83 ms / 28 行 | 🔴 s/d=20 引用值巧合；🟡 span 11.995 |
| `08_span_12_judge.py` | `docs/08` | `run_08_span_12_judge.log` `71887104775c` | 632 ms / 36 行 | 🔴（後由 `11` 升 🟡）；MC 0.166% |
| `09_density_gate.py` | `docs/09_script09_density_gate.md` | `run_09_density_gate.log` `06ca6b95974c` | 7 291 ms / 55 行 | 🟢 統一 protocol；全部 sector PN ≥ 30.98% |
| `10_gamma1_gamma6_bridge.py` | `docs/10_script10_gamma1_gamma6_bridge.md` | `run_10_gamma1_gamma6_bridge.log` `7e7b38f94f1a` | 1 010 ms / 51 行 | 🟢 譜代數精確；🔴 階梯族不可辨識 |
| `11_lattice_flag_span.py` | `docs/11_script11_lattice_flag_span.md` | `run_11_lattice_flag_span.log` `24512bc505e7` | 973 ms / 57 行 | 🔴 T1/T2；🟡 T3（span=12, 0.14σ）|
| `12_A_spectral_identity.py` | `docs/12_script12_A_spectral_identity.md` | `run_12_A_spectral_identity.log` `e2d938abba0d` | 3 133 ms / 38 行 | 🔴 0/13 104 命中；密度 5.21% |
| `13_mdmu_second_instance.py` | `docs/13_script13_mdmu_second_instance.md` | `run_13_mdmu_second_instance.log` `60364add23b2` | 259 ms / 40 行 | 🔴 11/5 降級（密度解釋）|

**編號重排說明（append-only）**：本 session 早前建立嘅 6 份輪次報告原為 `docs/01–06`，與「docs/NN ↔ scripts/NN」規則衝突；
2026-10-24 調律者指令「只保留 01–13」⟹ 輪次報告**內容全部併入對應腳本文檔後刪檔**（合併位置見下表），
`docs/` 現時只有 `docs/01–13`。register 舊 event 嘅 ref 由翻譯層（`register.py` `REF_PATH_FIX`）指向合併後文檔，**歷史 event 不改**。

| 已刪輪次檔 | 併入位置 | 併入內容 |
|:--|:--|:--|
| `r01` | `docs/01–08` §六、`docs/01` §七、`docs/05` §七 | 閘 8 支期望 vs 實跑；[C] TypeError 取證；05 MC protocol 不符 |
| `r02` | `docs/09` §六 | protocol 理由、夸克/輕子/hadron 逐命中全表、跨 sector 結論 |
| `r03` | `docs/10` §六–§九 | 代數恆等細節、範疇錯誤判詞、O89 階梯核證表、log 階梯真值、R2 裁定表 |
| `r04` | `docs/11` §六–§七 | FLAG 數據來源表、誤差傳播方法、1σ 帶、T4 density 細節、裁定表 |
| `r04 §四` | `docs/13` §五 | T4b 第二實例全表 |
| `r05` | `docs/12` §六–§七 | v 版本表、γ₆² 差異註、族結構、M_comp 比對、下一步、裁定表 |
| `r06` | `docs/13` §六 | 16 項 derived-vs-not、gap G1–G6、自我約束記錄 |

**合併方式（2026-10-24 追加，調律者糾正「標題全數缺失」）**：併入時**全文照錄**（原 r 檔 H1 標題作 §附錄 標題、原 H2 標題降一級保留原字，
一符一字未改），唔用摘要改寫；逐檔核對 36 個原標題全部存在。`docs/01–08` 另保留逐支閘核對表（一指標對一支）。
