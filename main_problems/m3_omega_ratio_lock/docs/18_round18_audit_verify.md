# M3 docs/18 — 第十八輪：Audit-Verify（全鏈審計，2026-11-03）

> **觸發**：用戶「執行 Audit-verify，檢查上次水位，目標是所有 M3 未驗證腳本與未審計文檔，確保推導全鏈連續性及一致性」
> **水位**：M3 **從無審計基建**（無 L2／L3／driver／audit）⟹ 首次全驗，17 支腳本 ＋ 18 份文檔全覆蓋。

---

## §1 四數（必報）

```
L1    17/17   每支跑兩次 stdout byte-identical 且 == canonical（rc 全 0）
L2    45 keys ＋ 10 identities   raw-input 獨立重算（零 import L1、零讀 docs）
L3    27/27   claims；19/19 份文檔有 claim（含 final_result）
audit 14/14   PASS ｜ FAIL 0 ｜ WARN 0
```

## §2 審計項（14）

| # | 檢查 | 結果 |
|:-:|:--|:--|
| A1 | 每支 L1 有 canonical log | ✅ 17/17 |
| A2 | L1 兩跑 identical 且 == canonical | ✅ 17/17 |
| A3 | 路徑衛生（scripts 無 `/tmp` 寫入）| ✅ 0 命中（跳過審計自身）|
| A4 | driver 腳本集 == L1 宣告集 | ✅ 17/17 |
| A5 | 數值形狀（無 `nan`／`inf`）| ✅ 0 命中 |
| A6 | 每份在碟 doc 有 ≥1 claim | ✅ 19/19 |
| A7 | 被推翻 doc 有前向撤回指針 | ✅ docs/03、docs/11 已補 |
| A8 | README 引用 docs 全存在 | ✅ 0 缺檔 |
| A9 | 四件套齊全 | ✅ 4/4（final_result 本輪補建）|
| A10 | 產物時序（log mtime ≥ script mtime）| ✅ 0 違例 |
| A11 | 負控制：L3 錯值必 FAIL；log 篡改必 FAIL | ✅ 2/2 |
| A12 | 無點號式內部碼（`m14.1` 型）| ✅ 0 命中 |
| A13 | PASS 行數 == 總結 | ✅ |
| A14 | 審計自身兩跑 byte-identical | ✅ |

## §3 本輪捉到之真缺陷（3 項，已修）

| # | 症狀 | 根因 | 處置 |
|:-:|:--|:--|:--|
| **D1** | `README` 引用 `docs/13_round13_index_audit.md` 但檔案不存在（A8 前置）| 第十三輪只寫入 `01_problem_record §二十三`，未落 doc | 補建 `docs/13`（內容同 §二十三，明標 2026-11-03 補建）|
| **D2** | `docs/03`、`docs/11` 仍以撤回前讀法作活文檔（A7）| 只是 §二十三 單向記錄，被推翻 doc 無前向指針（skill §4.7b #2）| 兩檔尾 append 更正塊（含 `3 = rank(P_κ)` 撤回、`3` 改由特徵標供給）|
| **D3** | 四件套缺 `final_result.md`（A9）| 未收官 | 本輪補建（≤ 段落：狀態／結論／關鍵式／依賴／審計／保留／禁區／邊界）|

## §4 新基建（本輪建立）

```
scripts/l2_independent_verify.py    45 keys、10 恆等式（差 < 1e-30）
scripts/l3_doc_crosscheck.py        27 claims；Fraction／e-記法／千分位／上標 全支援
scripts/18_audit_verify.py          14 檢查（含 2 負控制）
scripts/run_chain.sh                17 L1 + L2 + L3 + audit（預設 dry-run，--run 實跑）
```

## §5 一致性結論

- **推導鏈連續**：17 輪腳本 → 17 canonical log → 19 份 doc，逐層有 producer；無孤兒、無 inline 幽靈。
- **數值一致**：L2（獨立路徑）與 L1（原腳本）在全部 45 key 上無分歧；10 條代數恆等差 < 1e-30。
- **文檔一致**：27 條 headline claim 全命中；小整數（`rank_phys=2`）用錨點窗，避免零鑑別力。
- **已知邊界**：`Ω_b` 仍 +2.33σ（R2，等 CMB-S4）；統計上限 3.25σ（整數三分空間 1711 個）已定。

## §6 檔案

```
scripts/l2_independent_verify.py   data/l2_independent_verify.json
scripts/l3_doc_crosscheck.py       data/l3_doc_crosscheck.json
scripts/18_audit_verify.py         data/run_18_audit_verify.log
scripts/run_chain.sh               （driver）
```