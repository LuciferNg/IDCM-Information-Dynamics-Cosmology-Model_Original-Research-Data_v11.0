# docs/41 — M12 階段性 close-out（T0–T6，2026-10-29）

**日期**: 2026-10-29 18:05 HKT
**性質**: **記錄層 close-out**（依 skill `idcm-closeout` §2 T0–T6）—— 零新推導、零新判詞。
**授權（調律者原話）**：「然後做一次 M12 全鏈階段性 closeout」。

---

## 〇 判定（一句話）

**階段性 close-out 完成。** R1–R26 之三層驗證（L1／L2／L3）＋全鏈審計＋傳播審計**全綠**；
記錄**維持 `partial`，不 close、不裁決**（理由見 §八）。close 觸發器兩條不變。

---

## 一 T0 — 閉合水位（本輪首次建立）

水位表落 `01_problem_record.md §十`；本輪為**首次建立** ⟹ 依 skill 規則全部檔案**全檔驗**（無「水位內跳過」項）。

| 群 | 件數 | 驗證方式 |
|:--|:-:|:--|
| 根檔（`01_problem_record`／`SESSION_START`／`README`／`final_result`） | 4 | 全文人工 ＋ L3 交叉核對 |
| `docs/01`–`docs/41` | 41 | 全文；數字由 L3 逐條核（55 CLAIMS） |
| `scripts/01`–`scripts/48` ＋ `l2`／`l3` | 50 | `15 --runs 2`（重現性，A3 逐支 ×2）＋ `28` C2（逐支 stdout ↔ log） |

---

## 二 T1 — 獨立重驗（L1／L2／L3 ＋ 兩支審計器）

| 層 | 命令 | 結果 |
|:--|:--|:--|
| **L1** | `scripts/15_full_audit.py --runs 2` | **FAIL = 0 ／ WARN = 0**（log 231 行，×2 穩定）；A3 逐支 ×2 `byte-identical`，**含 `47` 與 `48`** |
| **L2** | `scripts/l2_independent_verify.py`（零 import 主鏈） | **14/14 PASS** |
| **L3** | `scripts/l3_doc_crosscheck.py`（唯讀文檔核對） | **CLAIMS 55／55 PASS，FAIL 0，MISSING 0，未溯源 0**（機器值池 1472） |
| 全鏈 | `scripts/47_allchain_audit.py`（A1–A10） | **FAIL = 0 ／ INFO = 4** |
| 連續性 | `scripts/28_continuity_audit.py`（C1–C5） | **C1–C5 全 0**（`C3d = 4` 白名單披露） |
| 傳播 | `scripts/48_propagation_audit.py`（P1–P7） | **44/44 PASS**（×2 `byte-identical`） |

**首輪 FAIL 揭露（真缺陷 1）**：`47` 之 **A5**（四件套狀態行輪次同步）首跑 FAIL ——
`actual` 由 `final_result.md` 之 `## R(\d+)` 節推得 R25，而 `01_problem_record`／`README` 已宣稱 R26 ⟹
**真缺陷 = `final_result` 缺 R26 節**（非檢查器錯）。已補 `## R26` 節 ⟹ A5 PASS。

**第二輪 FAIL 揭露（真缺陷 2，根因級）**：`15` 之 **A3 對 `28` 判「跑前 canonical 快照漂移」**，且重跑不收敛。
追根因：**`28`／`47`／`48` 三支審計器從不自寫 canonical log**（只 `print` 到 stdout，log 靠 shell `>` 產生）
⟹ `15` 之 A3 假設「腳本會自行覆寫 canonical，避免自證」在此三支**不成立**：快照法退化成
「與上一次手動 redirect 的殘影比對」，任何狀態變更（新增 doc、register note）即漂移。
**修法（改結構，不改判準）**：三支各加 `_Tee`，`stdout` 雙寫（終端 ＋ `data/run_NN.log`），`atexit` 落碟。
修後 `15 --runs 2`：**FAIL = 0／WARN = 0，log 231 行**，連跑兩輪 `byte-identical`（`15`／`28` log 皆穩定）。
備份：`~/.hermes/backups/scripts_{28,47,48}_*.v23pre.20261029_194221`。

★ 分類依 skill `idcm-audit-verify` §10-A：**先懷疑 check，再懷疑被審對象**。本輪三類皆遇：
①`scripts/48` 首 FAIL 係檢查器 bug（快照解析 `M9（negative）`）；②`47` 之 A5 係真缺陷（`final_result` 缺 R26 節）；
③`15` 對 `28` 之 A3 係**根因級真缺陷**（canonical log 非自寫）。

---

## 三 T2 — 文檔完整性

| 項 | 結果 |
|:--|:--|
| `README.md` docs 表 ↔ `docs/` 檔數 | **41 vs 41**（`docs/41` 本檔） |
| 三檔狀態行一致 | `01_problem_record`／`README`／`final_result` **皆宣稱 R1–R26**（`47` A5） |
| L3 溯源 | 未溯源 **0**（1472 值池） |
| 圍欄平衡／U+FFFD | 全部 .md 偶數圍欄、零替換字元（`47` A8） |

**本輪修好的兩處真缺陷**：

1. **`final_result.md` 之「驗證｜水位」停留 R3 期** —— 原文寫「`L1 scripts/01–22 … 21/21`、`L3 55/55＋溯源 0（889 值池）`、`docs 18／scripts 24／data 17`」，與 R25 後之實況（docs 40／scripts 50／值池 1472）不符。已重寫為現況七行。
2. **v11 根 `README.md` 缺 M11／M12 兩列**（`docs/40 §P3`）—— 已補，並修 M1 狀態（「已註冊」→「部分閉合」）與快照（10 → 12 → 13，隨 M13 開）。

---

## 四 T3 — 時序完整性（禁回溯污染）

| 檢查 | 結果 |
|:--|:--|
| 本輪新增推導 | **0**（純記錄層）⟹ 無「doc 早於 script」風險 |
| 產物次序 | `scripts/48`（17:43）→ `docs/40`（17:45）→ `docs/41`（18:05）—— 腳本先於文檔 ✓ |
| 舊 doc 改動 | 5 份（`docs/17`／`24`／`28`／`33`／`34`）**只加檔頭 ⚠️ 指針**，正文原文零改動；改前備份 `~/.hermes/backups/*.v22pre.20261029_173932` |
| 時序索引 | `01_problem_record §三.2` **append** R26 一列；原標題「R4–R25」不改，另加追補註（append-only） |
| 跨版本 | 零跨版本寫入 |

---

## 五 T4 — 登記面（register）

| 動作 | 內容 |
|:--|:--|
| 收據 | 2 條（`【收據 A／M11】` route A；`【收據 B／M1】` `coverage(k,ε)` 閘） |
| 發還 | 1 條 `handoff M12 → M3`（R22 帳本對 M3 `P1` 之部分答覆＋阻礙） |
| note | 1 條（R26 落碟，含兩套慣例披露） |
| progress | 1 條（`verdict=partial`，不 close） |
| 事件數 | M12 由 43 → **48** |
| `check` | **HARD DRIFT 0**；`REF-MISSING`／`PID-NO-FOLDER`／`UNREGISTERED-DIR`／`BLOCK-DRIFT`／`STATUS-DRIFT`／`README-ONLY` **全 0** |
| 基線偏移 | `NO-BLOCK` 632 → **634** = M13 兩檔（無機器塊，與 M11／M12 同慣例）；`INFO-DB-ONLY` 1、`SEED-DISAGREE` 54 不變 |

---

## 六 T5 — `final_result.md`

- 狀態行 → **R1–R26 ＋ 階段性 close-out**（明寫「未 close、未裁決」）。
- 新增 **`## R26`** 節（P1–P7 六列判表 ＋ 兩項 ★）。
- 重寫「驗證｜水位」（原 R3 期數字失效）。
- 行數 **223 → 250**。

---

## 七 T6 — Git（v11 獨立 repo）

- **scoped commit**：只 `git add` 本 M 目錄 ＋ v11 根 `README.md`；`git diff --cached --name-only` 核對 staged 清單。
- message 前綴 `M12 stage closeout:`（可追溯計數沿革）。
- **不 push**：分支上若有他人未 push 之 commit（`git status -sb` 之 `ahead N` 大於本輪 commit 數），只 commit 不 push。
- 驗證：`git log --oneline -1` ＋ `git show --stat HEAD` ＋ `git status --porcelain <本目錄>` 應零殘留。

---

## 八 為何維持 `partial`（三條理由）

1. **唯一未封項唔在推導力之內**：`docs/37` 已證真缺口 = 「令 `⟨w,v_u⟩` 近零之**結構理由**」，而非再找形式族；該理由必須來自**新註冊算子**（建構層已開 **M13**）。在 M13 未回判前，本 M 不得自判。
2. **精度閘未達**：`σ_ρd ÷ (1/41) = 1.4391 > 1`、半格 `0.012195122` 小於 σ ⟹ **任何**有理目標都會落某格點 1σ 內 ⟹ 命中零資訊。此閘只能由外部量測（`σ_ρd ≲ 0.005`）解除。
3. **「推盡 ≠ 結論」**（2026-10-12 教訓）：R1–R26 已窮盡 v11 存貨內可算載體；窮盡嘅意思係「本階段無可再推」，唔係「問題已答」。

---

## 九 相關

- 上游：`docs/39`（全鏈審計）、`docs/40`（傳播審計）
- 產物：`scripts/47`／`scripts/48`；`data/run_15_full_audit.log`、`run_l2.log`、`run_l3.log`、`run_47_allchain_audit.log`、`run_48_propagation_audit.log`
- 承接：**M13 `m13_flavour_orthogonality`**（建構層；判定權仍在本 M）
- 備份：`~/.hermes/backups/*.v22pre.20261029_173932`（6 份 doc）、`*.v23pre.20261029_194221`（3 支腳本）
- 水位：`01_problem_record.md §十`

---

## 十 後記（append，2026-10-29 晚；本節不改上文數字）

本 close-out 之四數係 **R26 檢查點**（`15` log 231 行、L3 值池 1472、`docs/01–41`／`scripts/01–48`）。
其後調律者指示「最後先推 1-2」，執行路線②⟹ 新增 **R27**（`docs/42`／`scripts/49`），鏈隨之右移：

| 指標 | R26 檢查點 | R27 後（現況） |
|:--|:--|:--|
| `docs` | 41 | **42** |
| `scripts` | 50（`01–48`＋`l2`／`l3`） | **51**（`01–49`＋`l2`／`l3`） |
| `15` log | 231 行 | **235 行**（A3 含 `49`） |
| L3 值池 | 1472 | **1545** |
| `47` | docs 41／scripts 46 | docs 42／scripts 47 |

**T0–T6 之判定不變**（記錄仍 `partial`）；R27 只更正 σ 與閘值（`σ_ρd 0.03510 → 0.04478518`、`÷格距 1.4391 → 1.8362`），
並把正式判詞寫入 `final_result.md §判詞`。現況全鏈零漂移（七份 canonical log 兩輪 `sha256` 相同）。
