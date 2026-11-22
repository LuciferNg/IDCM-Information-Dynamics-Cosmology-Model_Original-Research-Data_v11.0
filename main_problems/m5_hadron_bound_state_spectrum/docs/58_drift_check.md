# M5 docs/58 — 全鏈數字漂移／斷鏈檢查（獨立檢查器）

> 腳本：`scripts/60_drift_check.py`｜log：`data/run_60_drift_check.log`｜JSON：`data/60_drift_check.json`
> 觸發：全鏈審計要求「文檔數據及腳本引用數字無漂移、無斷鏈」之可重複檢查（唔重跑推導，只核數字／引用／計數）
> ✅ **結果**：`D1` 常數出現 **149 次 ⟹ 漂移 0**；`D2` 斷鏈 **0**（cross-O 引用另核：P10 `docs/00_33level_structure.md` 確認存在於 v10.0）；
> `D3` 計數：`docs = 57`、`scripts = 62` ⟹ README 已同步（`01_problem_record` 之舊快照屬 append-only 歷史，不判漂移）
> 2026-11-22

## 一、前設問題（docstring 原文）

```text
M5 docs/58 配套 — 全鏈數字漂移／斷鏈檢查器（doc ↔ data ↔ 檔案存在性）
================================================================================
用途：M5 推導鏈之**跨檔一致性**獨立檢查（唔重跑推導，只核數字與引用）：

  D1 **常數漂移**：v2 正典常數在全部 `docs/*.md` 中之每次出現，須與正典值一致。
     例外白名單：出現行含「v1」「作廢」「舊」「撤回」「更正」「附註」者視為對照，不判漂移。
  D2 **斷鏈**：`docs/NN_*.md`、`scripts/NN_*.py`、`data/NN_*.json` 之引用須指向實際存在檔案
     （`docs/NN` 若係 cross-O 引用如 `m7_…/docs/27` 則跳過）。
  D3 **計數漂移**：`README.md`／`01_problem_record.md` 內宣稱之 docs／scripts 數須等於實際檔數。

判準（寫死）：
  V1 常數漂移 = 0 ⟹ PASS；任何非白名單行不等於正典 ⟹ FAIL 並列出行號。
  V2 斷鏈 = 0 ⟹ PASS。
  V3 計數漂移 = 0 ⟹ PASS。
  V4 本檢查器**不改任何檔案**，只報。
輸出：`data/60_drift_check.json`
```

## 二、檢查設計

| 檢查 | 內容 | 判準 |
|:--|:--|:--|
| `D1` 常數漂移 | 13 個 v2 正典常數（`a₀`、`t`、`k₀`、`ħc/a₀`、`A_max(v1/v2)`、`T₂`、`A_eff`、`T₃`、`k`、`κ`、`P₀(D)`、`W1 σ`）在全部 `docs/*.md` 之每次出現 | 非白名單行須與正典相容；白名單 = 行內含 v1／作廢／舊／撤回／更正／附註／退／廢／對照／錯 |
| `D2` 斷鏈 | `docs/NN_*`／`scripts/NN_*`／`data/NN_*` 之反引號引用 | 目標檔須存在；cross-O 引用（前 60 字有 `P/M/O` + 數字）跳過 |
| `D3` 計數漂移 | `README`／`01_problem_record` 宣稱之 docs／scripts 數 | 須等於實際檔數（歷史快照除外）|

## 三、結果

| 項 | 值 |
|:--|:--|
| `D1` 常數出現處 | **149** 次 |
| `D1` 漂移 | **0** ✅ |
| `D2` 斷鏈 | **0** ✅（cross-O 引用：P10 `docs/00_33level_structure.md` 已核存在）|
| `D3` 實際計數 | `docs = 57`、`scripts = 62` |
| `D3` 漂移 | README 已同步 ⟹ 0；`01_problem_record` 之 `docs=37／scripts=39` 快照屬**歷史記錄**（append-only，不判漂移）|

## 四、發現並修正

1. 🔴 **`docs/02` 之 cross-O 引用**：`P10 docs/00_33level_structure.md` 一度被判斷鏈 —— 實核該檔存在於
   `IDCM-v10.0--It-From-Beat/main_problems/p10_ktree_structure/docs/00_33level_structure.md` ⟹ **引用正確**，
   係檢查器需加 cross-O 白名單（已加）。**非**文檔錯誤。
2. 🔴 **README 計數漂移**：宣稱 `核心 61 檔`，實際 62（新增本檢查器）⟹ 已同步。
3. ⚪ `docs=37／scripts=39` 為 `01_problem_record` 之 2026-10-25 審計快照 ⟹ 保留為歷史，另於本檔 §五 給現值。

## 五、現值（2026-11-22）

| 量 | 值 |
|:--|:--|
| `docs/` | **57** 檔（`01`–`58`）|
| `scripts/` | **62** 檔（含 `12_full_audit.py`／`16_independent_verification.py`／`28_chain_continuity.py`／`l2`／`l3`）|
| `data/*.log` | 66 檔 |
| `data/*.json` | 62 檔 |

## 六、結論

1. ✅ **常數零漂移**：149 次出現全部與 v2 正典相容 ⟹ 06-28 之分支統一（v2）在全部文檔貫徹。
2. ✅ **零斷鏈**：本 O 內部引用全部可解析；唯一 cross-O 引用（P10）已核目標存在。
3. ✅ **計數已同步**（README）；`01_problem_record` 舊快照保留為歷史。
4. ★ **本檢查器可重用**：日後每次改文檔後跑 `scripts/60_drift_check.py` 即得 D1–D3 三位一體檢查。

## 七、與其他檔案關係

- ← `docs/27`（v2 裁決）、`docs/55`／`56`（分支重推）、全部 `docs/*`
- → `README.md`（計數同步）、`final_result.md`、register
