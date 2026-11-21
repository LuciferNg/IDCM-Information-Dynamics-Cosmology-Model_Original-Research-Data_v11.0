# docs/07 — 審計（機器閘 ＋ 三層重算）

**日期**: 2026-11-14
**腳本**: `scripts/06_audit_verify.py`、`scripts/l1_reproducibility.py`、`scripts/l2_independent_recompute.py`、`scripts/l3_doc_crosscheck.py`

---

## 一、核心腳本水位

| 腳本 | 檢查 | RESULT |
|:--|:--|:--|
| `00_reproduction_gate.py` | G1–G7（含 G6a–e） | **11/11 PASS** |
| `01_frame_census.py` | C1–C7 | **8/8 PASS** |
| `02_bit_decomposition.py` | B1–B7 | **8/8 PASS** |
| `03_admissibility_and_readings.py` | D1–D6 | **8/8 PASS** |
| `04_cbt_scan_density.py` | E1–E5 | **5/5 PASS** |
| `05_registry_build.py` | F1–F5 | **5/5 PASS** |
| `07_slot_scope_theorem.py` | H1–H6 | **6/6 PASS** |

核心合計 = **51/51**。

## 二、三層重算

```
L1  scripts/l1_reproducibility.py   15/15 PASS
    8 支腳本各兩跑 byte-identical（md5 逐支列出）；7 支與入庫 log md5 相同
    （06 自審腳本只比兩跑，唔比 log —— 其 log 內容含自身 docs 狀態）
L2  scripts/l2_independent_recompute.py   7/7 PASS
    零 import 本輪任何腳本；閉式重算 17,280 / 21 / 4.392317 / 14.076816 /
    4320·720·2880 / 鴿籠 6 > 4 / γ₁ 兩式一致（<1e−58）/ 7 奇 ＋ 2 偶 = 9
L3  scripts/l3_doc_crosscheck.py   3/3 PASS
    30 條 claim token 逐字命中產物；docs/01–06 全覆蓋；
    docs 內多位小數／科學記數全部有產物支持（容差 5e-3，0 條未支持）
```

## 三、機器閘 A1–A8（`scripts/06`）

```
A1 腳本 8 支、全含預註冊塊                     PASS
A2 每支有對應 log（data/run_0N_*）              PASS
A3 RESULT n/n 與 chk(" 呼叫數逐支一致            PASS（含多行 chk 呼叫）
A4 產物綁定 6 個 JSON                           PASS
A5 路徑衛生（唯一外部絕對路徑 = 00 之唯讀 grep 根）PASS
A6 引用閘重放 PASS/BLOCK/BLOCK                  PASS
A7 import 白名單（無他 M 模組）                  PASS
A8 docs ≥6 份且含判段                            PASS

RESULT 9/9 PASS
```

總水位 = **51 ＋ 9 ＋ 15 ＋ 7 ＋ 3 = 85/85**。

## 四、零跨 M 寫入（硬閘）

```
所有寫入 = m19_mode_index_registry/{scripts,docs,data}
他 M 目錄：零寫入（只作唯讀 grep：m12/m13/m14 之 final_result 與 docs）
v10 目錄：只唯讀（grep 出處），零寫入
```

## 五、備份（pre-patch，硬編碼）

```
scripts/00_reproduction_gate.py.pre_patch_20261119_200646
scripts/01_frame_census.py.pre_patch_*
scripts/02_bit_decomposition.py.pre_patch_*
scripts/03_admissibility_and_readings.py.pre_patch_*
scripts/04_cbt_scan_density.py.pre_patch_* ×2
scripts/06_audit_verify.py.pre_patch_* ×3
scripts/l3_doc_crosscheck.py.pre_patch_*
```

## 六、禁區遵守（逐條核）

| 禁區 | 本輪做法 |
|:--|:--|
| 由質量值反解 `k_q` | 質量只用於①排序（`docs/01` 定理）②`C_f` 展示；`docs/04` 正控制明標輸入 |
| φ-冪比值 | 零使用 |
| `3 = h11/12` | 零接觸；`3` 只作 `N_c = B−1` 引用 |
| 形狀族 `δ_x` | 零重推 |
| 指數掃描凑 c/b 質量 | 掃描**唔含任何質量目標**；目標為模式集合，密度閘 0.0% |
| 重跑 M12/M13/M14 腳本 | 零執行（只 grep） |

## 七、判

```
✅ 水位 85/85；三層重算全綠；零跨 M 寫入
✅ 預註冊全部先於實跑（8 支 docstring）；兩處 append-only 澄清：
   ① scripts/00 之 Δγ 恆等式常數（π/44 → π/88）—— 推導錯，非閾值放寬
   ② scripts/01 C5 之「6/9」→「5/9」—— 手數錯，實跑為準
✅ 判據零放寬：T2 觸發凍結失敗判準（🔴），T4 觸發診斷更正判準（🔴）
```

## 八、本輪未做（誠實）

```
⚪ 未測「非子集型逐粒子命名規則族」——已由 docs/03 之 rank 定理覆蓋（任何此類命中必經質量 rank）
⚪ 未改任何他 M 檔案（依紀律，只經 register handoff）
⚪ 未處理 S5 載體（屬 M20）
```

---

## 九、收尾審計（2026-11-20）——覆蓋全鏈，取代 §一 之水位

**新增閘（本輪）**：`scripts/run_chain.sh`（driver）＋ `scripts/08_reverse_numeric_coverage.py`（L3b）
＋ `scripts/06` 擴至 A1–A14。§一 之「核心 51／85」為當時值，**現行水位如下**。

| 支 | 檢查 | RESULT |
|:--|:--|:--|
| `00` | G1–G7 | 11/11 |
| `01` | C1–C7 | 8/8 |
| `02` | B1–B7 | 8/8 |
| `03` | D1–D6 | 8/8 |
| `04` | E1–E5 | 5/5 |
| `05` | F1–F5 | 5/5 |
| `07` | H1–H6 | 6/6 |
| `08` | R1–R7（反向覆蓋） | 4/4 |
| `l1` | 兩跑 byte-identical（9 支 ×2） | 17/17 |
| `l2` | 獨立閉式重算 | 7/7 |
| `l3` | 30 claim 逐字命中 | 3/3 |
| `06` | A1–A14 | 15/15 |

**總水位 = 11+8+8+8+5+5+6+4 ＋ 17 ＋ 7 ＋ 3 ＋ 15 = 97/97**（`run_chain.sh` rc=0）。

### 九.1 新增三閘之意義

```
A9  driver 所跑腳本集 == 磁碟腳本集（差集空）      ← 防「支腳本從未入鏈而 audit 照綠」
A10 零 __pycache__／.pyc                          ← 路徑衛生（py_compile 曾遺 .pyc）
A11 四件套狀態同步（四份皆 🟢／已閉合）             ← 修好 01_problem_record／SESSION_START 之 stale 狀態行
A12 register ref 全部 resolve（6 條）              ← handoff／close／note 之 ref 指向實檔
A13 零 /tmp／scratch／cache 讀寫                    ← 排除 06 自身（禁字清單自含）
A14 pre-patch 備份 ≥8 份                           ← 改檔必先備份
```

### 九.2 反向數字覆蓋（L3b）結果

```
pool = 1637 值（data/run_0{0,1,2,3,4,5,7}.log ∪ data/0*.json；排除 audit／L1／L2／L3／L3b 自身）
doc token（4+ dp 或科學記數）= 51（docs/08 加入後重掃：51 → 見 §九.5）
UNTRACEABLE = 0；ALLOW = 0 條（零豁免）
負控制 9.876543 被判 UNTRACEABLE ⟹ 掃描器有鑑別力
```

### 九.2a 附註（append-only，2026-11-20 收尾輪末）

- **ALLOW 由 0 → 1 條**：`docs/08_open_problems_inventory.md` 加入掃描範圍後，`docs/07 §九.2` 自身引述之負控制值 `9.876543` 被掃成 UNTRACEABLE。該值係掃描器自證用，**非物理量**，故 `scripts/08` 增 `ALLOW = {"9.876543": "負控制值"}`（逐條列明理由）。同時修 `R1` 之硬編碼檔數（11 → 動態 `len(docs)+4`）。§九.2 之「ALLOW 0 條」為當時值，現行為 **1 條**。
- 掃描檔數 = **12**（docs 8 ＋ 四件套 4），非 §九.2 當時之 11。

### 九.3 本輪修好之真缺陷（三條，全部為審計自身）

| # | 缺陷 | 修法 |
|:-:|:--|:--|
| 1 | `08` 掃描器把自己輸出之負控制值 `9.876543` 收入 pool ⟹ 負控制自我中毒 | 排除 `08_reverse_numeric_coverage.json` |
| 2 | `06` A13 之禁字清單（`/tmp/` 等）令自審腳本自我命中 | A13 跳過 `06_*` |
| 3 | 四件套中 `01_problem_record`／`SESSION_START` 狀態行仍為「⚪ 已註冊、未開推」 | live field 直接改為 🟢 已閉合（§八／§一 原文保留 ＋ 附註指向） |

### 九.4 driver 覆蓋核（A9 實測）

```
run_chain.sh 跑 12 支：00,01,02,03,04,05,07,08,l1,l2,l3,06
磁碟腳本集（0*.py 9 支 ＋ l*.py 3 支）= 12 ⟹ 差集空 ✅
```

## 十、判（收尾審計）

```
✅ 97/97 全綠；driver rc=0；反向覆蓋 0 UNTRACEABLE（零豁免）
✅ 六項核對：①腳本集對稱 ✅ ②逆向綁定 ✅ ③路徑衛生 ✅ ④產物綁定（7 JSON）✅
             ⑤無 ghost 推導（全部 docs/NN 標明 producer script＋log）✅ ⑥driver ✅
✅ 四件套同步；register ref 6/6 resolve
```

