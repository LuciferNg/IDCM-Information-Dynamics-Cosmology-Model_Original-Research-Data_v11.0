# M7 docs/13 — 檔案刪除審計（2026-10-26 13:33 事件）

> 觸發：調律者「並行 session 檔案已刪除，檢查 scripts／logs／docs 是否有誤刪」
> 狀態：✅ **無誤刪自身 artifact**｜⚠️ 6 個並行 session 檔案已離碟（全未 commit）｜🟢 **6/6 副本完好保存在 scratch**
> 腳本重跑核：`01`／`02b`／`03` 三次重跑 ⟹ JSON **byte-identical**、`19/19`＋`13 ✅`、凍結 hash 不變

## 一、判定（一句）

**本 session 之 21 個 artifact 一個不缺**；離碟之 6 個檔案全屬並行 session、全屬**未 commit**（git 不可復原），但 scratch 有完整副本。

## 二、證據 1：自身 artifact 完整（git 三重核）

| 核 | 結果 |
|:--|:--|
| `git status --porcelain` | **空**（無任何 tracked 檔案被刪） |
| `git show --name-only HEAD` 21 檔逐一 `test -f` | **0 MISSING** |
| 三腳本重跑 vs 已 commit JSON `diff` | `01`／`02b`／`03` **完全一致** ⟹ 零依賴已刪檔 |

盤點（現況）：root 4（`01_problem_record`／`final_result`／`README`／`SESSION_START`）· `docs/` 3 · `scripts/` 4 ＋ `deprecated/` 1 · `data/` 9 = **21**，與 commit `669f28b` stat 逐項吻合。

## 三、證據 2：離碟清單（逐項對照 13:01 盤點）

| 路徑 | 大小 | git | 副本 |
|:--|--:|:--|:--|
| `scripts/01_repro_gate.py` | 15162 | 從未 commit | 🟢 scratch `d577a467dd7f…` |
| `scripts/02_definition_freeze.py` | 8473 | 從未 commit | 🟢 scratch `2cd41c46f7dc…` |
| `data/01_repro_gate.json` | 1848 | 從未 commit | 🟢 scratch `0ad2f06e3d89…` |
| `data/02_definition_freeze.json` | 1048 | 從未 commit | 🟢 scratch `30369fd256a4…` |
| `data/run_01_repro_gate.log` | 2900 | 從未 commit | 🟢 scratch `5aab13273a8c…` |
| `data/run_02_definition_freeze.log` | 6223 | 從未 commit | 🟢 scratch `308f6abe67bc…` |

副本位置：`/home/wsl/.hermes/cache/scratch/m7_wrongdir_backup/`（`scripts/`、`data/` 結構完整）。
另 `~/.hermes/backups/autosave/93139986/01_repro_gate.json.*`（autosave 快照，hash 不同，非權威）。
計數核：`data/` 13→9（−4）、`scripts/` 6→4（−2）⟹ 缺失 **6**，與上表全等 ⟹ **無第七項被刪**。

## 四、證據 3：我方結論不依賴已刪檔

`docs/10` §五 之兩約定重校表（`a = 0.384184 fm`、`t = 1.7711 MeV`、`⟨B⟩ = 0.4812`、`P₀ = 0.3891`）
係由**本 M 自己** `scripts/02b` 內建之 `a1, t1, o1`（`w = v²` 約定）算出，見 `data/run_02b_r5_mode_content_nogo.log` 之 `w=v² 約定重校：a=0.384184 fm, t=1.7711 MeV`；
`docs/12` §3.4 之 181×／829×／2150× 同源。⟹ **weight 約定審計之數字有獨立可重跑來源**，並行檔案刪除不影響任何結論。

## 五、dangling 引用（歷史陳述，本檔為更正註記）

以下 4 處引用現指向不存在之檔案，屬**當時存在**之歷史陳述，**原文不動**：

| 檔 | 位置 | 現況 |
|:--|:--|:--|
| `final_result.md` | 交付層 ⚠️ 並行 session 註記 | 已離碟（副本見 §三） |
| `docs/10` §首 ⚠️ 註 | 「並行 session 另有 `01_repro_gate.py`」 | 已離碟 |
| `docs/10` §五 | 「並行 `01_repro_gate.py` 指出…雙重計數」 | 已離碟；數字由本 M `02b` 獨立重現 |
| `01_problem_record.md` §十 末 | 並行 session 註記 | 已離碟 |

## 六、是否應復原（建議）

- **不建議**自動復原入 `scripts/`：會令 R0 重現閘出現**兩套同題腳本**（`01_reproduction_gate.py` vs `01_repro_gate.py`），屬命名混淆風險源，正是今次審計要消除之物。
- 若並行 session 仍需要：由**該 session** 自行由 scratch 取回，或另立目錄（如 `scripts/parallel/`）明標來源。
- 我方立場：`docs/10`／`12` 之所有數字均可由本 M 三腳本零依賴重跑復現（已驗 byte-identical）。

## 七、與其他檔案關係

- ← `docs/10` §五、`docs/12` §3.4、commit `669f28b`（21 檔 stat）
- → 無（本檔為審計終點）