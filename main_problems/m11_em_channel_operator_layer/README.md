# M11 — em 通道算子層（U5 攻堅）


**pid** `M11` ｜ **created** 2026-10-28 ｜ **version** v11.0

**狀態**：🟡 partial（em 通道八處封閉＋C 路線死；route A 已移交 M12）

## 問題（v2，R1 收窄後）

是否存在**零新參數**結構算子 `Ô`，其 `⟨O⟩_ch ≡ ∂M_h/∂(m_d−m_u)` 同時滿足 5 錨值（π 0／N 1.1263／Σ 1.9090／Ξ 2.1719／K 2.4213）並定出 3 條 out-of-sample 通道 `{Σ⁰−Σ⁺, Σ⁻−Σ⁺, D⁺−D⁰}` 之 `(⟨O⟩, em)` 分割？全滅則寫 scoped no-go（對 inventory `ℐ`），不得寫「不存在」。

原陳述「重現 8 條 `E_ch`」不可證偽（`E_ch` 係 em 與 `⟨O⟩≠1` 之混合量；8 條中 5 條在正核兩體算子類內不可行）。

## 判詞摘要（細節見 `01_problem_record` §五–§十五）

- 🔴 **C9**：5/8 通道不可行；`K⁺−K⁰`、`Σ⁰−Σ⁺` 屬**符號封鎖**（與 `k₀`、`κ` 無關）
- 🔴 **formula-first 終結**：`1040`（R2）＋`18468`（R3 後驗凹形族）候選，**0 PASS**
- 🟡 **κ 係自由參數**：由 M5 自身 K 兩數釘死 `κ=0.300`（帶 [0.2975,0.305]）⟹ `a₀=0.1878 fm`、`k₀=7.67 MeV`
- 🔴 **R4**：A／B／C 三類別模型 OUT；**D 飽和 `dof=0` 唯一存活但不可鑑別**
- 🔴 **G0（HF）**：`∂H/∂X = −(h/m_l)|0⟩⟨0|`，幅度差 **128–141×**，D/K 0.903 無分辨
- 🔴 **路線 C（總量 Δm）**：`χ²/dof = 3154`（`σ_X` 版 87.8）、`max|z| = 128`；死因係**模型空間退化**（26.7σ／9.0σ）
- ★ **`n_s` 依賴係真結構**：脫離 em/strong 分解後，單一 `d→u` 取代隨旁觀者給 1.293／4.807／6.850 MeV（與 M7 BMW 三錨獨立收斂）
- ⚫ **八處封閉**（mass 導數／sector pin／chain index／flavour→mode／常數族／算子槽位／質量參數／2 向量總量）＋**M7 `docs/28` 味盲定理＝第九處**

## 結論與交接

`⟨O⟩_ch` 不可由 v11 存貨推出（scoped）。**解鎖唯一條件**：lattice isovector scalar 矩陣元。route A（CY₃ 投影層 flavour index）→ **M12 `m12_cy3_flavour_operator`**（register `handed_off` 2026-10-28）。

## 檔案

| 類 | 數量 | 內容 |
|:--|:--|:--|
| 四件套 | 4 | `01_problem_record`（主記錄，append-only，§一–§十五）／`SESSION_START`／`README`／`final_result` |
| docs | 13 | `01` R1 判準｜`02` M4 交接裁決｜`03` R2 掃描｜`04` R3 凹形族｜`05` R4 必要性｜`06` R5 `k₀` 身份｜`07` M5 `docs/52` 對接｜`08` κ 決策面｜`09` em 錨溯源｜`10` G0 判死｜`11` 路線 C 判死｜`12` 入站 handoff 全量裁定｜`13` v2 支全鏈重算 |
| scripts | 14 | `01`–`11`（R1–R6＋審計鏈）＋`12`（全檔審計）＋`l2_independent_verify`＋`l3_doc_crosscheck`＋`run_chain.sh` |
| data | 13 JSON＋13 log | `01`–`12` 各一 JSON／log，`l2`／`l3` 各一 JSON／log |

## 驗證（四個數 + 全檔）

`L1 11/11 byte-identical`｜`L2 16/16`｜`L3 42/42`｜`audit 24/24`（FAIL 0、WARN 1＝文檔引 M5 `scripts/54` 跨記錄）｜`12_full_file_audit` 全檔（4 根檔＋11 docs＋14 scripts＋13 JSON＋13 log）

重跑入口：`bash scripts/run_chain.sh --full`（唯一標準；`--fast` 只作快速回歸）。

## 上游／下游

- **上游**：M4 `docs/19`（em 攻堅點＝矩陣元）、M4 `docs/37 §一 B4`（版本紀律下懸空 ⟹ 本記錄承接）
- **下游**：`handoff M11→M12`（route A）、`handoff M11→M5`（`scripts/54` 之 v2 支 `A0`）、`note M11→M12`（σ_O 傳播與達標分析）