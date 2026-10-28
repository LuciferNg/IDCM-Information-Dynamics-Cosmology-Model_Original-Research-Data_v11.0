# M4 docs/41 — 推導鏈審計（Audit-Verify，2026-10-28）

> 腳本：`scripts/41_chain_audit.py`（六項核對）＋ `scripts/run_chain.sh`（L1 全鏈）＋ `scripts/l2_independent_verify.py` ＋ `scripts/l3_doc_crosscheck.py`
> log：`data/run_41_chain_audit.log`、`data/audit_20261028/`、`data/run_l2_independent_verify.log`、`data/run_l3_doc_crosscheck.log`
> JSON：`data/41_chain_audit.json`、`data/l2_independent_verify.json`、`data/l3_doc_crosscheck.json`
> 依據：skill `idcm-audit-verify` §2（三層）＋ §3（六項核對）＋ §4（傳播）

## 一、前設問題

| # | 問題 |
|:-:|:--|
| Q1 | M4 之 41 支 L1 腳本能否在**當前樹**逐支重跑且零非零退出？ |
| Q2 | 重跑後 `data/*.json` 是否 byte-identical（即是否確定性、無隱藏狀態）？ |
| Q3 | 所有 headline 數字能否由**原始輸入**獨立重推（L2），而非只靠 L1 重跑？ |
| Q4 | 文檔所載 headline 數字是否等於 L2 之獨立值（L3），有無 over-precision？ |
| Q5 | 腳本集 ↔ 文檔集、產物 ↔ 生產者是否對稱（無孤兒、無幽靈）？ |
| Q6 | 有無「只存在於 session／文檔」而無落碟 producer 之數字？ |
| Q7 | 四件套（record／SESSION_START／README／final_result）之狀態與數字是否同步？ |

## 二、事前判準（跑之前寫死）

| 項 | 判據 |
|:--|:--|
| L1 | 全部 `rc=0`；每支**跑兩次**且 log／JSON byte-identical |
| L2 | 33 項全 PASS；期望值只可出現於 `IDENTITIES` 段；零 L1 import、零 docs 讀 |
| L3 | 全部 claim PASS 且 `MISS = 0`；期望值只可來自 L2；容差依文內有效位（≥6 位 `1e−4`／≥4 位 `5e−4`／其餘 `5e−3`） |
| audit | 六項核對全 PASS，或差異**已明文申報**（不可靜默放寬） |
| 任何 FAIL | 判「發現」而唔係改容差；發現一律 append-only 記錄 |

## 三、執行

```
L1：40 支（01–40，含 31b/36b）逐支 /usr/bin/python3，OMP/OPENBLAS/MKL_NUM_THREADS=1
    → data/audit_20261028/exit_codes.txt（全 rc=0）
L2：33 checks（pslq 精確最小多項式、Fraction 精確 rank、mpmath dps=50 三角）
L3：25 claims（跨 docs/11、18、29、36、37、40）
audit：12 checks（六項核對 ＋ L2/L3 閘 ＋ L1 覆蓋）
```

確定性：`l2`、`l3`、`41` 三支各跑兩次，JSON md5 與 log **byte-identical**（已 `cmp` 核實）。
L1 之 `data/*.json`：30 個中 **28 個**與審計前 md5 完全一致；差異 2 個 = `23_full_audit.json`、`24_audit_round2.json` —— 兩者**本身就係檔案計數型審計器**，我新增 `docs/41`、`scripts/41`、`run_chain.sh` 後其內部計數改變，屬 skill §7.2 陷阱 17 之預期行為，非缺陷。

## 四、六項核對結果（`scripts/41`，12/12 PASS）

| # | 核對 | 結果 |
|:-:|:--|:--|
| ① | 腳本集對稱 | 41 支 NN 腳本；`script-only = {31b, 35, 36b, 41}`、`doc-only = {32, 37}`，**全部已申報**（共用 doc／宣告式） |
| ② | 逆向綁定 | 5 條 headline 數字 → `data/*.json` → producer，**三重反查零 gap** |
| ③ | 路徑衛生 | `scripts/*.py` 內 `/tmp` **真寫入模式 0 命中**（`docs/11` 之「/tmp 命中」係歷史敘述，非現行寫入） |
| ④ | 產物綁定 | 32 個 `data/*.json` 全有 producer、零孤兒；多產者 3 個全部已申報 |
| ⑤ | inline 幽靈推導 | 本輪全部 headline 數字皆有落碟腳本（`39`、`40`、`l2`、`l3`、`41`）；無「只在 session 算過」者 |
| ⑥ | driver | `scripts/run_chain.sh` 存在（本輪新建），鎖 `OMP/OPENBLAS/MKL_NUM_THREADS=1` |

## 五、發現（F1–F9）

| # | 級 | 發現 | 處置 |
|:-:|:--:|:--|:--|
| F1 | 🟡 | `23_full_audit` 報 **21 違規**，其中 19 條為「`docs/26–40` 缺 §一 前設問題／事前判準／結論報告」＋ docstring 與 docs §一 不一致 → **歷史結構債**，非本輪造成 | 本輪補 `docs/39`、`docs/40`、`docs/41`（新寫者一律帶三節）；`26–38` 列待辦（§九） |
| F2 | 🟡 | `scripts/36`(v1) 與 `scripts/36b`(v2) **同寫** `data/36_em_kaon_external.json`，內容不同 ⟹ 最終檔由**執行次序**決定；`23` 因此報「JSON 不相交」 | 已申報（`41` 之 `DECLARED_MULTI`）；v2 為準（`docs/36` header 已載）。**未改腳本**（需授權） |
| F3 | 🟡 | `σ_R = 0.039896`（`docs/11`／`docs/37`）係 L1 之 2×10⁵ 樣本 MC 值；L2 解析傳播（只計 `σ_mu = 0.08`）= **0.039714**，相差 **0.46%** ⟹ 超出 `5e−4` 容差，**不可逐位重現** | 判決不變（`dev/σ = 0.0461%/0.331% = 0.14σ`）；`docs/11` append 附註；L3 已明文 `EXCLUDED` |
| F4 | 🔴 | `docs/37` C 欄仍寫「`⟨O⟩` 線外部否證 **16σ**」，而正式值（`docs/36 §十一`＋`docs/40`）為 **−9.76σ**（`−15.98σ` 係「`X` 精確」約定值）⟹ live summary 與正文**不一致** | `docs/37` append 更正（§八） |
| F5 | 🟡 | `scripts/__pycache__/23_full_audit.cpython-314.pyc` 殘留（另一解釋器版本） | 已刪；建議 `.gitignore`（未做，需授權） |
| F6 | 🟡 | `README.md` 之檔數宣稱（docs 23／scripts 23／log 24／json 14）與實際（docs 39／scripts 41／json 32／log 41）嚴重不符 | 本輪更新（§八） |
| F7 | 🟡 | `rank_ℚ = 20` 於 L1 **無 producer JSON**（只 doc 宣稱，源自 O137/M6 交叉核） | 現由 **L2 A7** 精確提供（Fraction 消元），逆向綁定改綁 `l2_independent_verify.json` |
| F8 | 🟡 | **L2 覆蓋缺口**：pool 枚舉型數字（`97.9%`、`1.97%`、`5.59%`、`17 identity`、`p = 0.550`、`7.2×`）與上游產權數字（`41.2σ`、`R² = 0.310001`、`a ≤ 0.9729 fm`）**未入 L2**（L2 禁 import L1，重建 pool 枚舉器需獨立實作） | 已在 L3 `EXCLUDED` 逐條明文列出；列下一輪（§九） |
| F9 | 🟡 | `24_audit_round2` 之 C1/C2/C5 屬**恆真 FAIL**（它比較「上次審計之後有無改動」），導致 `23` 報 `24: exit 1` | 非缺陷；建議改用 info 級（未改，需授權） |

### 5.1 易誤判點（記錄用）

`data/audit_20261028/run_27_em_sign_theorem.log` 內有 **4 個「FAIL」字串**，全部係**物理判詞標籤**（「PASS 1／FAIL 2 通道」＝符號預言之成敗），**唔係**腳本自檢失敗。以 `grep FAIL` 判腳本健康會誤報 ⟹ 應以 exit code ＋ `RESULT: n/n PASS` 為準。

## 六、三層驗證（四個數）

| 層 | 結果 |
|:--|:--|
| **L1** | **40/40**（全 `rc=0`；28/30 JSON byte-identical，2 個為檔案計數型審計器之預期變化） |
| **L2** | **33/33 PASS**（`data/l2_independent_verify.json`） |
| **L3** | **25/25 PASS，MISS 0**（`data/l3_doc_crosscheck.json`；另有 7 條 `EXCLUDED` 明文申報） |
| **audit** | **12/12 PASS**（`data/41_chain_audit.json`） |

## 七、水位（2026-10-28 審計基準）

| 檔 | 已驗範圍 | 狀態 |
|:--|:--|:--|
| `scripts/01–41`（41 支） | 全文重跑 ＋ exit code ＋ 確定性 | ✅ 已驗 |
| `scripts/l2`、`l3`、`run_chain.sh` | 全文（新建） | ✅ 已驗 |
| `docs/11, 18, 29, 36, 37, 40` | headline 數字（L3 25 claims） | ✅ 已驗 |
| `docs/01–10, 12–17, 19–28, 30–35, 38, 39` | **未**逐字 L3（只經 23 之結構檢查） | ⏳ 待驗（§九） |
| `data/*.json`（32） | producer 綁定 ＋ L1 重跑 | ✅ 已驗 |

## 八、四件套同步（本輪已做）

| 檔 | 動作 |
|:--|:--|
| `01_problem_record.md` | registry `updated: 2026-10-28`；append §二十八（本審計）；水位見 §七 |
| `SESSION_START.md` | append 一行「審計輪」標記（第 1 輪快照不動） |
| `README.md` | 狀態行改「三層審計」；檔數宣稱更新為實數 |
| `final_result.md` | **本輪新建**（原缺；M4 為 `partial`，但仍需活摘要） |

## 九、未處理（待授權）

1. **F1 後半**：`docs/26–38` 補 §一 三節（13 份，多為歷史輪次；要否補需裁定）。
2. **F2**：`scripts/36b` 改寫 `data/36b_*.json`（避免與 v1 撞檔）。
3. **F8**：建獨立 pool 枚舉器（覆蓋 6 個 pool 型數字）。
4. **F9**：`24_audit_round2` 之 exit code 語義。
5. **F5**：`.gitignore` 加 `__pycache__/`。

## 十、結論報告

1. **推導鏈完整可重現**：L1 40/40、L2 33/33、L3 25/25、audit 12/12，且三支新腳本兩次執行 byte-identical。
2. **本輪發現一個真錯（F4，🔴）**：`docs/37` 之 `16σ` 與正式值 `−9.76σ` 不一致 —— 屬「活摘要未跟正文更新」，已 append 更正；判決（KILL）不變。
3. **一個不可逐位重現之數字（F3）**：`σ_R` 由 MC 而來，與解析傳播差 0.46%；判決不受影響，但文檔精度應降為「MC 值」。
4. **覆蓋缺口已量化（F8）**：6 個 pool 型數字 + 3 個上游產權數字未入 L2，已逐條列出，不假裝全覆蓋。
5. **M4 狀態不變**：`partial`；`Q1 m_u = 2.1498 MeV` 仍為唯一活預言；close 觸發器不變（等獨立 `m_u` <0.5%）。本審計**未**改變任何物理判決。

## 附註三（2026-10-28；F3 結案，append-only）

F3（`σ_R` MC vs 解析，差 0.46%）**已結案**：canonical `σ_R = 0.039833`（精確二階 `0.039833024`；3×10⁷ MC 覆核 0.0398356 ± 0.0000046）。
舊 MC `0.039896` = 2×10⁵ 樣本漲落；L2 舊解析 `0.039714` = 一階近似（低估 0.2825%）。
`scripts/l2` 已改為計算精確值（新增 `R_span_sigma_exact` ＋ check `B4`），`scripts/l3` 之 `σ_R` 例外項已移除並改為 claim。
證據：`scripts/44`（4/5 PASS；A1 判準被否證＝爭議之解釋）、`data/44_sigma_R_resolution.json`。


## 附註四（2026-10-29；新增 F10，append-only）

**F10 `scripts/30` 之 `chi2_sf` 奇 `dof` 分支有 bug**：`sum(... for j in range(0))` ⟹ 第二項恆空，`k = 3` 返回 **1-dof 值**（與 M7 指出 `docs/29` 之同型錯誤相同）。
- **影響**：`docs/30` 之 `p` 欄低估（判準用 `χ² < 3`，**判詞不受影響**）。
- **處置**：`scripts/30` 屬凍結歷史腳本**不改**（避免回溯污染 `data/30`）；正確 `p(3 dof)` 表記於 `docs/47` §六、`docs/30` 附註；新腳本一律用 `mpmath` 正則化上不完全 gamma（`scripts/47` 之 `chi2_sf`）。
- **狀態**：🟡 記錄性（無判詞影響）


## 附註五（2026-10-29；新增 F11，append-only）

**F11 產物綁定之 producer 判定過寬（假 multi-producer）**：舊規則「腳本內有 `json.dump`／`'w'` ⟹ 該檔所有 `.json` 名皆算產者」把**純讀取**誤判為寫入。
- **實例**：`scripts/47` 讀 `data/42_flag2024_x_revision.json`，被報 `multi-producer`（`scripts/42` ＋ `scripts/47`）⟹ `4b FAIL`。
- **根因修法（已落地）**：判**行 ＋ 上一行**之語境 —— 寫模式／`json.dump`／`write_text`（同行），或「路徑常量賦值」（同行或上一行，且該語境無 `json.load`／`read_text`）。不再以「檔案層」為單位。
- **結果**：`4a orphans 0`、`4b multi-producer {}`（原本靠申報掩蓋嘅 `l2_independent_verify.json`／`23_full_audit.json` 假陽性一併消失）、`4c missing 0` ⟹ `12/12 PASS`。
- **副作用**：`DECLARED_MULTI` 內三條申報保留為歷史（現時已無 multi 產者）。
