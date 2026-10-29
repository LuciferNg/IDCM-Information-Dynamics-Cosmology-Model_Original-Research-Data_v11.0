# docs/39 — 全鏈審計（All-chain audit）：六項核對＋四件套＋register＋覆蓋披露

**日期**: 2026-10-29
**腳本**: `scripts/47_allchain_audit.py`（A1–A10）
**log**: `data/run_47_allchain_audit.log`（36 行，`stdout == log`，**×2 byte-identical**）
**依**：skill `idcm-audit-verify` §3.1（①–⑥）、§0.1（四件套同步）、§7.3（L2／L3 覆蓋披露）、§7 陷阱 9／10

---

## 一、四個數（交付格式 §8）

```
L1 46/46（A3 ×2 byte-identical）｜ L2 14/14 ｜ L3 55/55（溯源 0，1468 值池）｜ 全鏈 A10：FAIL 0
```

- `scripts/15_full_audit.py --runs 2`：`PASS 89`、`FAIL 0`、`WARN 0`（log 223 行）
- `scripts/28_continuity_audit.py`：C1–C5 全 0、`C3d = 4`（白名單披露，不計 FAIL）⟹ **全綠**
- `scripts/47_allchain_audit.py`：**FAIL 0 ／ INFO 4**

## 二、首輪發現（5 FAIL）→ 三分類後處置

| # | 項 | 首輪結果 | 分類 | 處置 |
|:-:|:--|:--|:--|:--|
| 1 | A1 log 無 script `[15, 28]` | FAIL | 🟡 **檢查器 bug** | logs 側同樣排除 `EXCL`（15／28 之 log 屬彙總器自身） |
| 2 | **A5 狀態行輪次** | FAIL：`01 = R0`／`README = R3`／`final = R9` vs 實際 **R25** | 🔴 **真缺陷** | 三檔狀態行同步至 R1–R25（見 §三） |
| 3 | **A6 時序索引覆蓋** | FAIL：**15/38 = 39.5%**（缺 `docs/13`–`26` 等 23 份） | 🔴 **真缺陷** | `01_problem_record` 新增 **§三.2** 補錄 R4–R25（26 列） |
| 4 | A7 register ref 缺 `[1, 6, 35]` | FAIL | 🟡 **檢查器 bug** | 原以 `pid == "M12"` 過濾，漏 `target_pid`（handoff 事件）；改以 `register.py log --pid M12` 為準 |
| 5 | A8 圍欄 `docs/12` | FAIL | 🟡 **檢查器 bug** | `docs/12:67` 之 ` ``` ` 係**表格內文字**，非圍欄；檢查改為**只認行首** ` ``` ` |

⟹ **檢查器 bug 3 個、真缺陷 2 個。** 真缺陷皆屬**四件套同步**（skill §0.1 最常見反模式 #4）。

## 三、傳播表（真缺陷）

| # | 位置 | 類型 | 舊 | 新 | 處置 | 狀態 |
|:-:|:--|:--:|:--|:--|:--|:--|
| 1 | `README.md:3` | 狀態行 | `R1/R2/R3 完成` | `R1–R25 完成`＋精度界限 | 改 | ✅ |
| 2 | `final_result.md:2` | 狀態行 | `R1–R9 完成＋三層審計；5 條存活路線` | `R1–R25 完成＋三層＋全鏈；五類 no-go；精度界限` | 改 | ✅ |
| 3 | `01_problem_record.md §八` | 狀態行 | `⚪ registered（未開推）` | append 🟡 `partial — R1–R25`（原文不動） | append | ✅ |
| 4 | `01_problem_record.md §三.2` | 時序索引 | 只至 `docs/12` | +`docs/13`–`38`（26 列，逐輪一句） | append（新節） | ✅ |
| 5 | `SESSION_START.md` | 快照 | `⚪ registered` | **不動**（skill §0：第 1 輪快照，只 append） | 豁免 | ✅ |

★ 傳播後重跑：A5 `{01:25, README:25, final:25}` = 實際 R25 ✓；A6 `38/38 = 100.0%` ✓。

## 四、A1–A10 逐項（終態）

| 項 | 內容 | 結果 |
|:--|:--|:--|
| **A1** | 腳本集對稱 ↔ log 配對：`docs 38 份`／`scripts 45 支（EXCL 15,28）`／`logs 45 份`，三側差集空 | ✅ |
| **A2** | `docs` ↔ README 表格：**38 vs 38**，無漏列、無多列 | ✅ |
| **A3** | 產物綁定：`data/*.json` **無孤兒**（18 份皆有產者） | ✅ |
| **A4** | 路徑衛生：`scripts/*.py` 之 `/tmp` 寫入 **0** | ✅ |
| **A5** | 四件套狀態行：三者皆 **R25** ＝ 實際 | ✅ |
| **A6** | 時序索引覆蓋：**38/38 = 100.0%** | ✅ |
| **A7** | register ref 覆蓋：`M12` 事件覆蓋 **docs/01–38 全** | ✅ |
| **A8** | 圍欄偶數（行首）＋ `U+FFFD = 0`：全部 `.md` 通過 | ✅ |
| **A9** | 更正留痕：`docs/34`（撤回）／`docs/35`（撤回）／`docs/37`（更正）／`docs/38`（定理）皆在 | ✅ |
| **A10** | 覆蓋披露：`15` 之 `*.py` glob **含** L2／L3（49 支）；`28` 之 `[0-9][0-9]_*.py` 不含，但以 `AUX` 收 `run_l*.log` | ✅ INFO |

## 五、覆蓋範圍（誠實披露，skill §7.3）

- **`scripts/15`**：`SCRIPTS.glob("*.py")` ⟹ **含** `l2`／`l3`（重現性覆蓋全 49 支）。
- **`scripts/28`**：`glob("[0-9][0-9]_*.py")` ⟹ **不含** `l2`／`l3`；但已以 `AUX` 收 `run_l*.log` 作溯源池。
- **`scripts/47`**：A7 以 `register.py log` 子進程取 ref（不直接讀 1 MiB+ 之 `register.jsonl`，避開掃描上限）。

⟹ 三者合計：**重現性**（15）＋**溯源**（28 AUX）＋**記錄鏈**（47）三向覆蓋，`L2`／`L3` 無盲區。

## 六、判定

🔴 **無待修項。全鏈審計 PASS（FAIL 0）。**

四件套自此**四份一致**（`SESSION_START` 依 skill 屬快照豁免）；`01_problem_record` 時序索引補齊至 R25，
下一個 session 讀任一份都會得到同一個狀態：**`partial`，R1–R25，卡精度界限**。