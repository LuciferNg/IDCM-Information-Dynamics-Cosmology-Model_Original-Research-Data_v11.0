# M1 — 質量檔案 schema：帳本編碼（ledger encoding）與判別力

<!-- registry
pid:     M1
title:   質量檔案 schema：帳本編碼（ledger encoding）與判別力（Reality-Architecture; 由 v10 O251 遷移）
state:   partial
verdict: partial
flags:   migrated-from-O251; route-ledger-encoding; format-layer-only; no-flavour-map; first-round-done; r6-nogo
owner:   M1
created: 2026-10-14
updated: 2026-10-28
prereq:  M2 + M4 + M5 + M6 + M8 + M9 + M10 + M11; O211/O212/O204/O237/O239/O140/O216/O144/O120/O209; O115
source:  調律者2026-10-14 O251→v11 M1；2026-10-28 整合 M2–M11 後重寫起手三件套；同日「開始推導M1」首輪完成
-->

## 狀態: 🟡 partial — 首輪＋R6＋下游重推完成（2026-10-29）｜鏈 8 支全 exit 0、JSON 重跑 byte-identical、文檔↔數據 171/171

**一句話**：M2 已證 `m = M₀Σσ_nγ_n`；M1 攻格式三問（正規形／多重性／判別力）。首輪結果：**質量值冇譜判別力（≤0.19 bit）**——此即 M10 W5「dof=0」之格式層證明。

**首輪七條結論**（詳見 `docs/05_prereg_results.md`）

1. 判別力梯度：簡併度 **μ 1 → s 7 → τ 513 → c 4285**；μ 之 ρ=61 對 null p = 0.113（coverage k≤4 = 17.040 %）⟹ 證據 ≈ 0.15 bit；其唯一 k=2 編碼 ⟹ 0.19 bit
2. 二值帳本不完備：9 粒中 **7 粒**在 k ≤ 12 無二值編碼（只 μ k*=2、s k*=3）
3. 輕扇區 = 連續輸入：e／u／d 之編碼係退化解 `S=∅, σ₄₄=t/4`（0.004578693557／0.01935420431／0.04184450655）⟹ 3 個連續 dof
4. 重粒子簡併爆炸：τ 513、c 4285 個等價編碼 ⟹ 重夸克「質量公式」必屬事後選擇
5. 跨量子格式缺口（**唯一新開放**）：b n_q=2、t n_q=69，量子分配律未定義
6. 三鏈 = 一骨架三 preset：O54 0 bit（F_proj 連續）、O239 φ² 域外（√5 conductor 5 ∤ 88）＋428.4σ、M4 Δm²+R 只代內
7. 五層約化：真內容 3 件（γ_n 譜表 × M₀ 單位 × σ）；topology／labels 置換不變 200/200 ⟹ 空層

**R6 跨量子（2026-10-28 同日續推）**：🔴 **no-go（不可識別）**
- E（等分）與 F（滿載＋餘額）兩律皆以零新常數精確重現 b、t ⟹ 質量值層不可識別
- 六個量子含量（h_b 30.184122161227／r_b 59.815877838773／h_t 24.898076881069／r_t 65.101923118931／c_b 74.907938919387／c_t 89.639158306071）k_min 7–39、ρ(k_min) 19–1393 ⟹ 全部無譜身份
- 預註冊 R6-Q1、R6-Q2 **皆 FAIL**（由單一均值外推／預設常數虧損律）
- 新禁區：D1 禁由質量反解量子分配；D2 禁以 γ₁ 格點命中當結構

**下游重推（2026-10-29，`docs/08`）**：M2–M12 依賴點全掃 ⟹ 7 條受影響已 register 推送（M2 守恆式推廣／M4 慣例＋支撐閘／M5 多體合成接口／M6 `n_q ≥`／M10 W5 加强／M11 量子分配裁示／M12 ρ 閘）；4 條不受影響（M3/M7/M8/M9，附證據）。★**17.040 % 係三慣例中最低者**（log-uniform 35.155 %、M4 慣例 24.877 %）⟹ 引用必標慣例；支撐閘 `coverage(k)≤0.1` 只到 k=2。★`M₀` 兩義**解除**：0.033331 MeV 係 v10 O239 `|Δm|^φ²` 之 1 參數 fit（標度被式吸收），非矛盾。

**預註冊對帳**：Q1 🟡 半、Q2 🔴、Q3 🔴、Q4 🔴、Q5 ✅、R6-Q1 🔴、R6-Q2 🔴（原文不改，修正見 `docs/05 §2`、`docs/06 §3`）

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 上游結構 + P1–P5 + R1–R5 + 驗收 + 預註冊 Q1–Q5 + 時序索引 + 禁區 T1–T11 + 界面契約 + **§十一 首輪結果** | 🟡 |
| `SESSION_START.md` | 起手式 + **§七 首輪進度** + 下一輪（R6）接手點 | 🟡 |
| `docs/01_layer_reduction.md` | R1 五層擾動測試 ⟹ 約化定理（3 件 + 2 空層 + 3 外部） | ✅ |
| `docs/02_multiplicity.md` | R2 ρ_bin／ρ_1p／coverage／MC null 10⁷ | ✅ |
| `docs/03_canonical_form.md` | R3 正規形 + 顯式編碼 + 簡併度階梯 + 退化解定理 | ✅ |
| `docs/04_three_chain_presets.md` | R4 三鏈對照 ⟹ 一骨架三 preset | ✅ |
| `docs/05_prereg_results.md` | 預註冊對帳 + 七條結論 + 完整判準 v1 + R6 預註冊 | ✅ |
| `docs/06_multiquantum.md` | R6 跨量子：兩律等價 ⟹ 不可識別（no-go）＋六含量譜測試 | ✅ |
| `docs/07_incoming_handoff_audit.md` | 入站 handoff 審計：入境 1 條判 closed-moot；出境 2 條（§7：M4 已接受改名並修正範圍） | ✅ |
| `docs/08_downstream_impact.md` | 下游影響重推：11 個 M 依賴掃描、P1–P6 重算、7 條推送 | ✅ |
| `scripts/00_layer_perturbation.py` | 擾動測試（4/4） | ✅ |
| `scripts/01_multiplicity_engine.py` | 多重性引擎 MITM（3/3） | ✅ |
| `scripts/02_mc_null.py` | MC null 10⁷（3/3） | ✅ |
| `scripts/03_canonical_form.py` | 正規形／簡併度（3/3） | ✅ |
| `scripts/04_three_chain_presets.py` | 三鏈（3/3） | ✅ |
| `scripts/06_multiquantum.py` | R6 跨量子（8/8） | ✅ |
| `scripts/07_downstream_propagation.py` | 下游重推（13/13） | ✅ |
| `scripts/05_doc_consistency.py` | 文檔↔數據 171/171 | ✅ |
| `scripts/run_chain.sh` | driver（預設 dry-run；`--run` 全跑 ~98 s） | ✅ |
| `data/00–04_*.json`、`data/run_*.log` | 輸出與執行記錄（重跑 byte-identical） | ✅ |
| `backups/*.20261028_1150.bak` | 2026-10-14 原文三件套（史料，不刪） | 📦 |

## 未閉（觸發器）

R6 已判 no-go ⟹ 本 M **唔 close**（按「推盡 ≠ 結論」）。重開觸發器（任一）：
1. M11 通道算子層提供量子分配機制
2. M5 束縛態給出多體合成律
3. 出現可直接量度 `n_q` 或量子含量的獨立觀測
樣本外保留：以 PDG 2022 重跑 `docs/02 §1` ρ 表。

## 依賴

- **骨幹**：M2（`docs/01 §1–§4`：泛函、帳本、partial 被迫）
- **結構界**：M4 `docs/10/19/22`、M6 `docs/07/14`、M9 README
- **語義**：M8 `docs/05/07/08`（M₀ = 單位、功率判準）
- **外部**：M5（束縛態）、M10（合法性）、M11（通道算子）
- **歷史**：v10 O251（`state=merged`）、O237、O239、O140、O212、O204、O216/O144、O120、O209

## 註冊

pid `v11.0/M1`（created 2026-10-14；register 現址 `IDCM-Original/problem_register/`，跨版本共用）