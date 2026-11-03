> ✅ **承接已完成（2026-10-30，另一 session）**：`🔴 closed / verdict negative`（四問全滅，scoped no-go）。
> 判決見 `final_result.md`；逐問見 `docs/02`–`docs/06`；發還見 `docs/07`。
> **原文以下全部保留不動**（`SESSION_START` 係起手式記錄，不作改寫）；「只開不推」之閘已於 `docs/01 §一` 記錄解閘。

---

# M14 SESSION_START — u,d,s 之味標籤與站指派：結構來源與可固定性

## 狀態

⚪ **registered（只開不推）**，2026-10-30 開。`docs/`／`scripts/`／`data/` **零檔案**。
調律者明示：「**開 M14 推 u,d,s 但在其他 session 推**」 ⟹ 推導**不在**本 session 進行。
解閘條件 = 承接 session 明示（或調律者指示）；未解閘前**不得**寫任何推導腳本。

## 起點（一句話）

M13 已證 `ρ_x` **只依賴 `α`、與站指派無關**（站無關性定理）⟹「邊個味住邊個站」在比值層**不可觀測**。
本 M 之問題**不是**再算 `ρ`，而是：**u,d,s 之（站指派、味標籤、序）三件事可否由既有結構固定；固定不了就要判 gauge，並講清換什麼觀測量才可固定。**

## 第一步（解閘後依序讀，缺一即停）

1. `../m13_flavour_orthogonality/SESSION_START.md`、`docs/09_completeness_audit.md`（殘餘 = `x=2` 之來歷；四項未盡）
2. `../m13_flavour_orthogonality/docs/05_universality.md`（**站無關性定理**：44 三元組，偏差 ≤1e−34）
3. `../m13_flavour_orthogonality/docs/10_sector_extension.md`（**槽數 `= rank(A₂)+1`**；第 4 味不可外推；`δ_c≈0`）
4. `../m13_flavour_orthogonality/scripts/01_reproduction_gate.py`、`scripts/10_family_x_scan.py`（須獨立重現，禁引用未核數值）
5. `../m12_cy3_flavour_operator/docs/14_r5_isospin_scan.md`（**二體荷 `δ_u/δ_d/δ_s`** 三錨唯一解；`δ_c`）
6. `../m12_cy3_flavour_operator/docs/28_r15_feasibility_theorem.md`（`cone=ℝ³`；**閘 3：判死須標範圍**）
7. `../m12_cy3_flavour_operator/docs/40_propagation_audit.md`（§六 四條恆久約束：`coverage`／慣例／配對／版本閘）
8. `../m4_mass_spectral_readout/docs/32_flavour_structure_scout.md`（**味標籤住 `K_chain`**：`k∈{1,4,5}`、`d∈{3,2,0}`）
9. `../m7_hadron_model_basis/docs/28_*.md`（**味盲定理 `F̂ ∉ 𝒜`**）
10. `../m1_mass_file_schema/docs/08_downstream_impact.md`（`coverage(k,ε)`：`ε=1e−3` ⟹ `k ≤ 2`）

## 承接須知（踩過的坑，勿重犯）

- **`σ_ρd` 正典 = 0.044785178**（M12 `docs/42` R27）。`docs/38 §E` 之 `0.03510` 已 SUPERSEDED —— 唔好再用。
- **`docs/06` 之「樹層核排除」範圍只及 `{4 個具體 w}`**，且同族 `x=2` 反而**通過**（M13 `docs/09 §四`）。族語言禁用。
- **`ρ` 命名衝突**：M12 `docs/14` 之 `ρ_s=1.148613`、`ρ_c=17.09` 係**介子/重子比值**；M13 之 `ρ_x` 係**重子電荷比**。跨記錄引用必寫定義。
- **M12 `docs/15 L61` 之 `Ξc` 預言 `−3.1430` 疑符號錯**（三種一致讀法 `|z| ≤ 1.9`；只有取負才 `11.85σ`）—— 已 `register note` 回 M12，裁決權不在本 M。若本 M 用到帶粲數值，先確認 M12 是否已更正。
- **只開不推**之定義：可以讀、可以 register、**唔可以**寫 `docs/`／`scripts/`／`data/`。

## 邊界

- 唔寫上游任何檔案（M1/M4/M7/M12/M13）；跨問題只經 `register`（`handoff`／`note`）。
- 判死必標範圍（單一／有限集／族）；`N_trials` 全報；禁後驗凑數。
- `dof = 0` 為硬閘；任何「指派權重」連續參數即違規。
- 本 M 不宣稱 M12／M13 閉合。

## 產出（解閘後）

`docs/NN_*.md` ＋ `scripts/NN_*.py` ＋ `data/run_NN_*.log`，一發現一 doc；`final_result.md` 只在全滅或結構閉合時寫。

> ⚠ **2026-10-31 附註（append-only）**：verdict 已由 `negative` 改為 **`partial`**（Q1／Q2②／Q3② 之死判撤回；理由見 `docs/09`、`final_result §附註二`）。原文（第 1 輪快照）保留不動。
