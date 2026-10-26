# M7 docs/14 — 全 Python 重跑 ＋ 文檔數字逐條核（2026-10-26 第二輪）

> 觸發：調律者「M7 全部 python 重跑，確定和文檔記錄相同」
> 狀態：✅ **4 腳本全 exit 0；JSON 全部 byte-identical；log 全部 identical；文檔數字 88/88 PASS（4 條經 append-only 更正）**
> 新腳本：`scripts/04_doc_consistency.py`（doc-vs-artifact 檢查器，零任意容差）

## 一、重跑（全 4 腳本，`/usr/bin/python3`、`PYTHONDONTWRITEBYTECODE=1`）

| # | 腳本 | exit | log 行 | JSON vs 已 commit | log vs 已 commit |
|:-:|:--|:--:|:--:|:--|:--|
| 1 | `01_reproduction_gate.py` | 0 | 52 | ✅ byte-identical | ✅ identical |
| 2 | `02_r5_mode_content_nogo.py`（v1） | 0 | 51 | ✅ byte-identical | ✅ identical |
| 3 | `02b_r5_mode_content_nogo.py`（v2） | 0 | 69 | ✅ byte-identical | ✅ identical |
| 4 | `03_grep_evidence.py` | 0 | 33 | ✅ byte-identical | ✅ identical |

- 重跑序：`diff` 逐檔比對（`/tmp/m7_pre/*.json` 為重跑前快照）⟹ **4/4 完全一致**（含 mpmath 30 dps 浮點尾部）。
- `01` 之 19/19；`02b` 之 13 條 ✅；`03` 之凍結 hash `883e91d5…92421` 不變。

## 二、文檔↔輸出逐條核（`scripts/04`，88 條）

**判準（無可調參數）**：文檔只印有限位數 ⟹ 判「文檔值 == 腳本值四捨五入到文檔所印精度」——
由文檔數字字串**自身位數**推 ulp，要求 `|doc − art| ≤ 0.5·ulp`；整數 claim 須全等。
（v1 版本曾用 rel-tol 判準，被廢棄：rel-tol 可調 = 可作弊；現判準由文檔精度唯一決定。）

### 覆蓋範圍
`docs/10` §三 19 項重現閘 ＋ §五 兩約定表 10 值 ＋ `docs/11` grep 表 11 詞 × 3 版本（33 格）＋ hash ＋ `docs/12` §二 K1–K8 表 ＋ §3.3–3.4 ＋ `final_result` 11 值 ＋ 跨版本 v1/v2 一致 2 項。

### 結果：**88/88 PASS**

| 類 | 數 |
|:--|--:|
| ✅ 直接相符 | 84 |
| ✅ 經更正行相符（原文不動） | 4 |
| 🔴 不合 | **0** |

## 三、發現並更正之 3 個記錄錯（append-only，不影響結論）

| # | 位置 | 原文 | 正確（腳本） | 性質 |
|:-:|:--|:--|:--|:--|
| 1 | `docs/10` §三 第16行 mpmath `⟨B^{−d}⟩` | `0.472772220017288211`，且註「float64 同、Δ 1.1e−16」 | `0.472772220017289335655`；float64 = `0.4727722200172883`（Δ **1.0e−15**） | 轉錄尾數錯 ＋ Δ 誤記 |
| 2 | `docs/10` §三 第17行 mpmath `P₀` | `0.368921706883373474`，且註「float64 同、Δ 5.6e−17」 | `0.368921706883374606256`；float64 = `0.36892170688337345`（Δ **1.2e−15**） | 同上 |
| 3 | `docs/12` §二 K4 ＋ §3.4 之 Ξ 通道 | `2150×` | `2151.377×`（`shortfall[2]` = 2151.37741233） | 末位四捨五入錯（0.06%） |

**判**：三者皆屬**記錄層**（轉錄／四捨五入），非計算錯；原文保留，更正置 `docs/10` §八、`docs/12` §九。
兩值之 15 位有效數字與全部 19 項判準皆不受影響；Ξ 通道 2151 vs 2150 皆遠超判準 ≥10²。

## 四、附帶核實

- v1 與 v2 之 `K1.max_dev`（`5.773e−15`）與 `K3.χ²`（`18.220727`）**一致** ⟹ v2 新增 K6–K8 未改動既有結果。
- `data/02_r5_mode_content_nogo.json` 與 `.v1.json` 內容相同（v1 執行 + 封存副本，已知重複）。
- 判準 FAIL 仍照實存在：`02b` K4 判準（≥10²）於 N 通道 36.5× **未達**，log 保留 🔴（未因重跑而改判）。

## 五、與其他檔案關係

- ← `docs/10`／`11`／`12`／`final_result`（被核之文檔）
- → `docs/10` §八、`docs/12` §九（更正）、`data/04_doc_consistency.json`（88 條逐項明細）
- 重跑閘 ⚠️：本目錄另有並行 session 之 `m6_scale_dof_theorem/{data,scripts}`（未追蹤），本輪**未觸碰**。
## 六、重跑所揭之可重現性缺陷（已修）：語料定義漂移

**現象**：第二次重跑時 `03_grep_evidence.json` 與 commit 版**不再一致**（1 條 claim FAIL）——
唯一變動 = 「質量泛函非可加」v11 欄 `3f/4h → 6f/10h`，新增命中全來自**未追蹤**之並行 WIP：
`main_problems/m6_scale_dof_theorem/{docs/01_t1_composite_additivity.md, scripts/02_t2_scale_dof_theorem.py, data/02_t2_scale_dof_theorem.json}`。

**根因**：`03` 以 `os.walk` 掃碟 ⟹ 語料 = 「此刻碟上有咩」，其他 M 一動計數即漂 ⟹ **不可重現**。

**修法（v2，不改變任何已記錄數字）**：語料改定義為**各版本 `git HEAD` 之 tracked files**（`git ls-files`），仍排除本 M；
JSON 新增 `corpus_scope` 與 `fingerprint`（每版本 HEAD ＋ 掃描檔數）：

| 版本 | HEAD | 掃描檔數 |
|:--|:--|--:|
| v9.0 | `60d67ae` | 407 |
| v10.0 | `151ec75a` | 10247 |
| v11.0 | `1ee1ee2` | 229 |

**核**：修後 11 詞 × 3 版本全部計數與原記錄**逐格相同**（只多出 fingerprint 欄），`docs/11` 表**無需改動**；`04` 重跑回 **88/88 PASS**。
⚠️ 副作用（正面）：`m6_scale_dof_theorem/docs/01_t1_composite_additivity.md` 之出現，顯示並行 session 正處理本 M 所 register 之 **M2 加性適用範圍** 交接項。

## 七、可重現性定案（03 v3）

再修兩點，令 `03` 之 JSON **跨 commit byte-stable**：

1. 語料身份 = `corpus_sha256` = sha256(sorted 「path\0blob_sha」 from `git ls-tree -r HEAD`) ⟹ 只隨語料**內容**變，不隨 commit 變。
2. **host repo（本 M 所屬 v11）之 head 不入指紋**（每 commit 自變）；v9/v10（外部 repo）保留 head 作 provenance。

最終指紋：

| 版本 | head | files | corpus_sha256（前12） |
|:--|:--|--:|:--|
| v9.0 | `60d67ae` | 407 | `299a36d3848a` |
| v10.0 | `151ec75a` | 10247 | `f71f6d39776e` |
| v11.0（host） | — | 229 | `d5050c59304f` |

**定案核**：commit 後再重跑 ⟹ **5/5 JSON 全部 == HEAD**（`01`／`02`／`02b`／`03`／`04`）；`04` 仍 **88/88 PASS**。
