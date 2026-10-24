# M4 — 質量 = 譜讀數差值：全費米子譜差閉合

## 狀態: 部分閉合（partial, 2026-10-24）｜審計雙軌：`docs/23`（全量）＋`docs/24`（增量，違規 0）

質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2，R7 定案承接）。
已確立：差值閉合力學（Δm²+R 直出譜，唔經 M_comp）＋ K_chain 譜代數（γ_n = 2−2T_n(c)，全譜 ∈ ℤ[γ₁]）。
攻堅路線：全費米子單點譜比（已全滅，density 閘）、lattice 精度重判（span≈12 升部分閉合）、A 譜身份（開放）。

2026-10-24 追加：**審計**（13+1 支全部 byte-identical 重現、14 檔文檔結構全核）、**開放問題整理**（`01_problem_record.md` §十一：Q1–Q8 ＋ 已判死路線 ＋ 三軸判讀紀律）、
**繼續推**（`docs/14_double_probe.md`：雙 probe 反相證明係**恆等還原**、G5 量化為 `R_em = 3.479111`（密度殺 34.73%）、`E` 無譜形式 0/13104）。

**Q7 追推（`docs/15_exponent8.md`）**：`M_EW = M_Pl γ₁⁸ N_d` 嘅指數 8 **反選確認** —— 8 anchor 只有 v 嘅 best k = 8（Λ／M₀ 皆 9）、`k_mech = 7.2817 ≠ 8`（M1 否證）、M2 唔支持、8 有 17 個 identity（冗餘）、density 5.59% ⟹ 🔴 判死，新增禁區。

**Q6 追推（`docs/16_scale_operator.md`）**：🟢 零參數 α_s 跑動（α_s(M_Z)=0.117985 ⟹ α_s(m_b) 0.236%／α_s(m_c) 1.680%／α_s(2 GeV) 1.182%）；🟡 scale 語義 = **自對偶點**（重夸克 `m_q(m_q)` 恆等、輕夸克無解 ⟹ 2 GeV 共同標度結構強制）；🔴 跑動指數數值檢定阻塞（無乾淨雙標度對）。

**Q8 追推（`docs/17_em_decomposition.md`）**：🔴 最小電荷模型（κ_sΣQ²+κ_xΣQ_iQ_j）**結構否證** —— 三通道係數退化 ⟹ 恆等預言 ΔK=−ΔN，觀測差 2.641 MeV＝125.7σ；🟡 H1（E_K=E_N）被數據偏好：r=1 落 lattice 容許 [0.941,1.252]、X=2.6137（+0.57σ）、隨機偶然率 0.41%，但仍屬 premise（未推導）。

**Q2/Q3 追推（`docs/18_em_constraints.md`）**：8 通道實證譜 ⟹ em operator 四類形式全排除（電荷加性 R²=0.003、sub-count R²=0.089、spectator 獨立 spread 3.578 MeV、單一幂律 worst 117%）；K/N 一致（spread 0.107）再證 H1；operator 必須含 mode-level 資訊 ⟹ **Q2/Q3 blocked-by-O140**（flavour→mode 映射）。

**O140／P6 交付處理（`docs/19_delivery_processing.md`）**：✅ O140 五項全核實（≤1.1e-13）；🔴 ℤ[γ₁] 多項式框架無鑑別力（deg-0 空洞：τ/e 0.0066%；deg≥1 密度 1.97% ⟹ Q5 判死）；硬約束：單模上限 `4/γ₁=784.9647` ⟹ 跨代/跨類必須多模；🟡 A 線索 `87γ₁⁴M₀=1.63673 eV`（0.404%，未過門檻）。

**可辨識性曲線（`docs/20_identifiability.md`）**：密度殺定量化 —— F1 單模比族 0.1% 帶覆蓋率 **20.61%**（1/5 隨機）；所需精度 δ★：F1 0.0004%／F4 0.0081%；質量比精度對各族短欠 **12–1 772×**（精度瓶頸），而 em 精度 0.003% 對 F4 **已足夠**（概念瓶頸）⟹ Q2/Q3 推唔動係 operator 未定義，非量度不足。

**上游交付核實（`docs/21_upstream_verification.md`）**：✅ O142 十項全核（≤9.4e-4）、O124 41.2σ 範疇 no-go、O85 R²=0.310001 精確重現 ⟹ **Q2/Q3/Q8 判死（群論＋Cartan 雙 no-go）、Q4 σ 橋通道關閉**；★真缺口＝**強子束縛態譜／`⟨1/r_ij⟩` 重疊算子**（IDCM 現無任何 O 承擔，見 §11.5）。

**O216／O86 交付核實（`docs/22_o216_o86_verification.md`）**：✅ 標籤層加性被 41.5σ 否證 ⟹ Q2/Q3 **六路線全封**；✅ O86 接手 A 算子身份（O86/N8 open）⟹ Q4 數值側 M4 終結；★P6 鏈 `m_ν=A·γ_n²` 之 m1/m2 反解 n = 4.7283／4.8245 **非整數**（§11.6，發還 P6）。

**全量審計（`docs/23_full_audit.md`）**：22 支腳本重跑 **22/22 byte-identical**、docstring **22/22** 對齊、14 JSON 語義全同、**違規 0**；修復 3 項（`14` log 標題 stale、`12`/`13` §一 節錄、審計腳本自審遞歸）；數字回溯 triage 完成（meta 68／external 264／pending 122 → 90 復算、32 跨文檔或外部，零無源數字）。

**第二輪增量審計（`docs/24_audit_round2.md`）**：審計「審計自身」＋審計後改動 —— C1 `scripts/` 自 `ea59cf9` **零改動**（首輪 22/22 結論可繼承）、C2 `docs/23` 宣稱↔JSON 逐項全對、C3 抽樣 `07`/`15`/`16` 重跑 **byte-identical**、C4 T3 未解釋 **0**、C5 審計後改動之數值**多重集零漂移** ⟹ **違規 0**；期間自我否證兩個審計器假 FAIL（T3 語料缺跨文檔引用值；C5 誤用「改動行數字」而非多重集）。

詳見 `01_problem_record.md`（§八／§十一／§十二）／`SESSION_START.md`／`docs/13`–`docs/23`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 已確立結構 + R1–R4 結果（§八）＋ append-only 更正（§九）＋ 逐腳本文檔登記（§十）＋ **開放問題整理（§十一）** ＋ **審計記錄（§十二）** | 部分閉合 |
| `SESSION_START.md` | 起手式 + 接手起點 + 2026-10-24 執行記錄（§七） | 部分閉合 |
| `docs/` | 23 檔（`01–23`）：逐腳本記錄＝前設問題＋事前判準＋執行＋結論報告（2026-10-24 全文重寫，輪次報告內容已融入正文） | 部分閉合 |
| `scripts/` | 計算腳本（23 檔：01–08 遷入、09–23 新增；含 G5 `14_double_probe`、Q7 `15_exponent8`、Q6 `16_scale_operator`、Q8 `17_em_decomposition`、Q2/Q3 `18_em_constraints`、交付處理 `19_delivery_processing`、可辨識性 `20_identifiability`、上游核實 `21_upstream_verification`、O216/O86 核實 `22_o216_o86_verification`、全量審計 `23_full_audit`） | 部分閉合 |
| `data/` | stdout log ＋ JSON（24 log ＋ 14 JSON，另 `audit_20261024/` 13 檔對照） | 部分閉合 |

## 編號規則

- `docs/NN_<slug>.md`：**一對一**對應 `scripts/NN_<slug>.py`（NN = 01…13），結構固定：一、前設問題（docstring 原文）→ 二、事前判準（掃描前寫死）→ 三、執行記錄／數據／結果 → 四、原始輸出（引自 log）→ 五、結論報告 → 末節 與其他腳本關係
- 輪次／閘報告**已併入**對應腳本文檔，唔再另立檔（2026-10-24 調律者指令：只保留 01–13）：
  - 閘全景（8 支期望 vs 實跑）→ `docs/01` §三；01 [C] bug 取證 → `docs/01` §六
  - 05 MC protocol 不符（8.7% vs 43.3%）→ `docs/05` §五
  - density 統一 protocol ＋ 逐命中全表（夸克／輕子／hadron）→ `docs/09` §二–§六
  - γ₁↔γ₆ 代數身份 ＋ O89 階梯核證 ＋ log 等距 → `docs/10` §三–§六
  - FLAG lattice 數據來源、四項測試、1σ 帶、裁定表 → `docs/11` §二–§五
  - A 歸一化、族 13 104 掃描、M_comp、下一步 → `docs/12` §三–§六
  - T4b 第二實例全表 → `docs/13` §三；本 M derived-vs-not（16 項）／gap G1–G6 → `docs/13` §五–§七
- 所有腳本一律「假設從未執行過」重跑：canonical run log 存 `data/run_NN_<slug>.log`（附 sha256）

## 依賴

O257 (R7 docs/07 定案) + O66 (M_comp 降級) + O89 (γ₁ 冪階梯，2026-10-24 核證：M_GUT/PS/R 三階梯作廢) + O57 (m_ν 盲測) + P6 (§10 差值閉合) + O73 (c=Goldstone) + O212 (dof=1)

## 註冊

pid `M4`（created 2026-10-21；state partial 2026-10-24；handoff → P6 A 身份、O140 selection rule）
