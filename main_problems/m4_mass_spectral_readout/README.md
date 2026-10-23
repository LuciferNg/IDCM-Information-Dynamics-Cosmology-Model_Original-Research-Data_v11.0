# M4 — 質量 = 譜讀數差值：全費米子譜差閉合

## 狀態: 部分閉合（2026-10-24）

質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2，R7 定案承接）。
已確立：差值閉合力學（Δm²+R 直出譜，唔經 M_comp）＋ K_chain 譜代數（γ_n = 2−2T_n(c)，全譜 ∈ ℤ[γ₁]）。
攻堅路線：全費米子單點譜比（已全滅，density 閘）、lattice 精度重判（span≈12 升部分閉合）、A 譜身份（開放）。

2026-10-24 追加：**審計**（13+1 支全部 byte-identical 重現、14 檔文檔結構全核）、**開放問題整理**（`01_problem_record.md` §十一：Q1–Q8 ＋ 已判死路線 ＋ 三軸判讀紀律）、
**繼續推**（`docs/14_double_probe.md`：雙 probe 反相證明係**恆等還原**、G5 量化為 `R_em = 3.479111`（密度殺 34.73%）、`E` 無譜形式 0/13104）。

詳見 `01_problem_record.md`（§八／§十一／§十二）／`SESSION_START.md`／`docs/13_mdmu_second_instance.md`／`docs/14_double_probe.md`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 已確立結構 + R1–R4 結果（§八）＋ append-only 更正（§九）＋ 逐腳本文檔登記（§十）＋ **開放問題整理（§十一）** ＋ **審計記錄（§十二）** | 部分閉合 |
| `SESSION_START.md` | 起手式 + 接手起點 + 2026-10-24 執行記錄（§七） | 部分閉合 |
| `docs/` | 14 檔（`01–14`）：逐腳本記錄＝前設問題＋事前判準＋執行＋結論報告（2026-10-24 全文重寫，輪次報告內容已融入正文） | 部分閉合 |
| `scripts/` | 計算腳本（14 檔：01–08 遷入、09–14 新增；含 G5 `14_double_probe`） | 部分閉合 |
| `data/` | stdout log ＋ JSON（15 log ＋ 5 JSON，另 `audit_20261024/` 13 檔對照） | 部分閉合 |

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
