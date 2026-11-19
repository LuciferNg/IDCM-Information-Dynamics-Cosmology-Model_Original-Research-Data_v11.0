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
