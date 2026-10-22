# M4 — SESSION_START（開題，2026-10-21）

<!-- registry
pid:     M4
title:   質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2, R7 定案承接；只開不推）
state:   registered
verdict: -
flags:   register-only
owner:   M4
created: 2026-10-21
updated: 2026-10-21
prereq:  O257 + O66 + O89 + O57 + P6 + O73 + O212
source:  調律者 2026-10-21: 開v11 新M, 將/tmp全部py移入scripts；新session重建推導鏈及報告文檔
-->

> **狀態**: 📄 已註冊 — 未開始（只開不推）

---

## 〇、調律者目標（原文）

「定案後，應該開新chapter, 這是理論大升級」→「好」→「開v11 新M, 將以上所有/tmp 的py 移入資料夾/scripts 中，要在新session 先重建推導鏈及報告文檔」

**執行狀態（本 session）**：✅ 開 M4 + 8 支 py 遷入 scripts/ + 註冊 register-only
**交接（新 session）**：🔴 重建推導鏈（重跑 8 支）+ 寫報告文檔（docs/）

## 一、接手第一步：重現源 scripts（失敗即停）

**逐支重跑 + 期望輸出**（全部改動路徑：scripts 內冇 hardcoded path 依賴，純計算）：

| Script | 內容 | 期望輸出 |
|:--|:--|:--|
| scripts/01_quark_scan.py | 15 質量比 vs γ比池 + m=γ²v / γv | m_c/m_s γ20/γ5 dev 0.0424%；隨機撞池 97.9% |
| scripts/02_mc_mcs_judge.py | m_c/m_s MC null | 撞池 97.9%＝後驗判死 |
| scripts/03_reading_strength.py | 讀數×譜強度 | m=γ²v 全部 >6%；γv 1.3-1.7% 事後揀 |
| scripts/04_lepton_anchor.py | 輕子 vs γ比池 | μ/e 1.15%、τ/μ 1.04%、τ/e 77%；MC 16.8% |
| scripts/05_hadron_dm.py | hadron Δm vs γ池 | n−p/m_e γ29/γ16 0.0195%（density 殺）；反相 insight |
| scripts/06_hop_probe.py | hop 階梯 vs 質量 | s/d=20 池零偏差（PDG 巧合）；hop 對比 |
| scripts/07_hit_followup.py | s/d 誤差判別 | m_d ±10% 誤差帶 [17.4,22.6]＝巧合 |
| scripts/08_span_12_judge.py | quark/ν span≈12 | R=11.994978 dev 0.042%；MC 0.166%；m_u 誤差 kill |

**閘**：8/8 跑通且輸出吻合 → 開始路線圖；任何一支數值唔對 = 停低報告。

## 二、已確立結構（開題快照，詳細見 01_problem_record §二）

- 🟢 中微子差值閉合（O257 R7）：m=譜讀數差值；m₃=γ₆²A（0.002%）
- 🔴 夸克/輕子/hadron 全滅：γ比池 density（MC null 16.8%-97.9%）
- 🟡 GUT γ₁ 冪階梯（O89）：M_Pl×γ₁ᵏ×{1,1,1,N_d}；M_EW→v 0.80%
- 🟡 quark/ν span≈12=hop_max（MC 0.166% 低但 m_u ±10% 誤差 kill）
- 🔴 γ₁↔γ₆ 結構橋未解（「邊界之間嘅邊界」）

## 三、禁區（唔重推）

1. **「質量比=γ比」直路**（已死：MC null 97.9%/63.4%/16.8%）
2. M_comp 絕對標度閉合（O257 R4-R6：三重 fit 循環、0.051% 降級量級 ±3%）
3. μ/e (γ₇/γ₁)(F_μ/F_e) 0.031%——F 因子後驗，唔當結構
4. s/d=20、span≈12——u/d/s 誤差帶內，後驗
5. O89 γ₁ 階梯直接當閉合證據——**先 RAG 核證 O85 通報後狀態**（scripts/02 已更新，但冇獨立 blind-test 審計）
6. M_comp = M_pl·10^(1/γ₆−1/76) 式任何變體（R4/R5 已否證）

## 四、路線圖（新 session）

1. **重現 8 支 scripts**（閘，§一）
2. **R2 γ₁↔γ₆ 橋**：γ₆/γ₁=35.47≈?、γ₁ᵏγ₆ˡ·A 混合、γ₁ 階梯 log 等距 2.29279 身份（事前定義族＋density 閘）
3. **R3 lattice 高精度值**：FLAG 2024 m_u/m_d/m_s 代入 span/比掃描（誤差收窄 10 倍後重判）
4. **R4 A 譜身份**：A=1.64044 eV/γ² 唔靠 M_comp 嘅譜讀數形式（v²/M_comp→γ 純量重寫）
5. **R5 報告文檔**：docs/01 承接定案／docs/02 掃描全滅判決／docs/03 GUT 錨分層
6. 完成後決定 close 定開子題

## 五、交接/依賴

- 上游：O257（R7）🟢 已定案；O66 M_comp 降級；O57 2.96% 先例；O89 γ₁ 階梯（⚠️ 需核）；P6 §10；O73；O212
- 本 M 唔郁 O257/O66/O57/O89/P6 任何檔案
- 引用 M3（Ω 鎖定）同 M2（同步↔脫同步）——同屬 v11 Reality-Architecture，本 M 質量讀數係 M2 同步邊界嘅具體化

## 六、已知陷阱

- γ比池 density 極高：任何單點 <0.1% 命中 = 後驗（MC null 先行）
- u/d/s 誤差 ±10%：0.01% 級「閉合」無意義
- O89 有污染史：引用前 RAG 核證
- scripts/01-08 由 /tmp 遷入：路徑零依賴但輸出格式要同 §一 期望值對返