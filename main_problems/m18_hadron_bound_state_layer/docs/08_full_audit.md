# M18 docs/08 — 全量審計
> 腳本：`scripts/08_full_audit.py`｜log：`data/run_08_full_audit.log`｜JSON：`data/08_full_audit.json`
> 狀態：✅ **24 PASS / 0 FAIL**（exit 0）｜2026-11-19

## 一、前設問題

```text
交付層審計：7 支推導腳本 × 7 份 docs 之可重現性、可追溯性、四件套同步。
A1 全部 01–07 重跑 exit 0，且 data/*.json byte-identical（sha256 前後相同）
A2 每支腳本有 docstring，且 docs/NN 存在
A3 data JSON 可解析且含 checks/pass
A4 docs/NN ↔ 腳本編號一一對應
A5 核心數字雙向可追溯（docs ↔ data JSON）
A6 本 M 目錄外零寫入（git status ＋ mtime 歸屬）
```

## 二、事前判準
見上 A1–A6（寫死；審計腳本自身先建立，後跑）。

## 三、執行

41.6 s。`PASS 24 / FAIL 0`，exit 0。

## 四、結果

**A1 重跑 byte-identical（7/7）**

| 腳本 | exit | sha256 前 12 | 判 |
|:--|:-:|:--|:-:|
| `01_reproduction_gate.py` | 0 | `72b48eae603d` | identical |
| `02_support_theorem.py` | 0 | `0e236ac88ee4` | identical |
| `03_composition_law.py` | 0 | `1a839d36d37f` | identical |
| `04_jpc_classification.py` | 0 | `f6d91da893a2` | identical |
| `05_radial_ladder.py` | 0 | `0c78f63fd4ed` | identical |
| `06_shape_parameter_identity.py` | 0 | `da3b3abbcf37` | identical |
| `07_feedback_gate.py` | 0 | `cbf1be3ef1a4` | identical |

**A2/A3/A4**：7/7 docstring 齊、docs 齊、JSON 齊。
**A5**：7/7 核心數字雙向可追溯（docs ↔ JSON）。
**A6**：git status 11 行變更 —— 本 M 3 行（`data/`、`docs/`、`scripts/`）；本 M 以外 8 行
**全部 mtime 早於 session 起點（2026-11-19 20:44:15）** ⟹ 屬其他 session／歷史殘留，**本 session 零越界寫入**。

## 五、判決

1. **交付層合規**（24/24）。
2. 首輪 audit 捉到並已修 3 項：
   - `04` 之 data JSON／log 陳舊（U4b 判詞改寫後未重跑）⟹ 已重跑；
   - `docs/02` 寫「8,388,608」帶千分位 ⟹ A5 追蹤失敗，已改為 `8388608`；
   - `docs/05` 寫「12,424×」帶千分位 ⟹ 同上，已改 `12424×`。
3. A6 判準 v1（「git status 零變更」）**過嚴**：repo 內本來就有其他 session 之未提交變更
   ⟹ 已改為 mtime 歸屬判準（本 M 外變更須全部早於 session 起點）。

## 六、與其他檔案關係

- 覆蓋 `docs/01`–`docs/07` 全部推導；
- 四件套狀態見 `final_result.md`、`01_problem_record.md` §十、`README.md`、`SESSION_START.md` §六。

---

## 七、審計擴充（v2；2026-11-19 append；原文 §一–§六 不動）

隨第二／三輪新增腳本，`DOCMAP` 擴至 **10 支**（加 `09_jpc_geometric_carrier`、`10_residual_census`），
並修 **A6** 判準（v1 過嚴，見下）：

| 判準 | v1 | v2（現行）|
|:--|:--|:--|
| A6 | 「git status 本 M 以外零變更」 | **A6a 硬閘結構證明**（`mode=enforce`、本 session 恰 1 scope、`roots` 只含 m18 ⟹ 越界寫入結構上不可能）＋ **A6b 歸屬**（本 M 外變更全部落喺他 M 目錄 ⟹ 並行 session）|

**為何改**：v1 假設「repo 內只有本 session 改動」，但實際有並行 session（實測 `m20_third_factor_carrier/`
於 11-20 02:53 新增 `scripts/04`），v1 會把他人改動算成自己越界。v2 把判據由「mtime 早晚」改為
「閘結構 ＋ 路徑歸屬」，唔依賴他人工作時間。

**v2 執行結果**

```
重跑 9/9 byte-identical（01–07, 09, 10）
A1–A6 全部 PASS ⟹ 31/31 PASS，exit 0
```

首輪 v2 捉到並修 4 項：`04`/`10` 之 data JSON 陳舊（判詞改寫後未重跑）→ 已重跑；
`docs/02`／`docs/05` 千分位數字令 A5 追蹤失敗 → 已改無千分位；
A6 v1 判準過嚴 → 已改 v2（本節）。
