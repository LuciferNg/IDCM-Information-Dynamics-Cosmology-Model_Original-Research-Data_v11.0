# M4 docs/49 — canonical 引用 registry ＋ 全鏈漂移／斷鏈檢查（2026-11-22）

> 觸發：調律者 8 點指令之 1–3（全鏈腳本齊備／未審計項審計／數字無漂移斷鏈）｜腳本 `scripts/49_citation_registry_check.py`（**5/5 PASS**，`docs/49` 落地後）｜數據 `data/49_citation_registry.json`

## 一、前設問題

M4 現有 **51 支腳本**（NN 49 支 ＋ `l2`／`l3`）、**45 份 docs**、**41 個 data JSON**。三問：
1. 每支推導腳本是否都落地為檔案、且其過程與結果有文檔記錄？
2. 文檔內每個被當作**現行**使用之數值，是否有 canonical 來源、且無被取代之舊值冒充？
3. 有無「腳本 → data → doc」斷鏈（無 producer／無 doc／無 data）？

## 二、事前判準

1. **R2**：每個 canonical 值必須在其來源 `data/*.json`（或 log）內以數值容差 `2×10⁻³` 追到 ⟹ 否則 = 斷鏈。
2. **R3**：舊值（已被取代）只可出現於**該文檔內已標明取代關係**者（含 `附註`／`禁區`／`legacy`／`撤回`／`→` 等標記）；否則 = **漂移違規**。
3. **R4**：每支 NN 腳本有 doc（或已於 `scripts/41` 申報共用）；每個 `data/*.json` 有 producer；無腳本之 docs 只可屬敘述型且 ≤ 2 份。
4. **R5**：未追到之 token 只作 **ledger**（不作硬 FAIL）—— 因為 O-era 上游數值、計數整數、衍生算術皆屬合法但不可由本 O data 追。

## 三、canonical registry（20 量）

| ID | 量 | canonical | 來源 | 已被取代之舊值 |
|:--|:--|--:|:--|:--|
| `R_span` | lattice span 比 R | **12.0055** | `data/11_lattice_flag.json` | 12.0049 |
| `sigma_R` | R 之 σ（精確二階） | **0.039833** | `data/44_sigma_R_resolution.json` | 0.039896、0.039714 |
| `X_ratio` | X 比值形式 | **2.5145** | `data/42_flag2024_x_revision.json` | 2.5600、2.56 ± 0.094 |
| `sigma_X` | σ_X 比值形式 | **0.1159** | `data/42_flag2024_x_revision.json` | 0.0943、0.094 |
| `epsilon` | FLAG ε | **0.79** | `data/42_flag2024_x_revision.json` | 0.78 |
| `em_K` | em(K⁺−K⁰) 錨 | **2.26459** | `data/36b_em_kaon_external.json` | — |
| `z_canon` | R5 canonical z（含 σ_X） | **9.76** | `data/46_handoff_adjudication_round4.json` | 15.98 |
| `O_N` | ⟨O⟩_N | **0.8906** | `data/46_handoff_adjudication_round4.json` | — |
| `O_Sigma` | ⟨O⟩_Σ（total−ΔQED） | **1.5469** | `data/46_handoff_adjudication_round4.json` | 1.4980 |
| `O_Xi` | ⟨O⟩_Ξ | **2.293** | `data/46_handoff_adjudication_round4.json` | — |
| `chi2_O1` | χ²(⟨O⟩≡1) 新錨 | **18.9723** | `data/47_downstream_rerun_round4.json` | 18.2207 |
| `p_O1` | p(3 dof) 新錨 | **0.000277** | `data/47_downstream_rerun_round4.json` | 3.96e-4、1.97e-5 |
| `window_a` | em 承載窗 a_max | **0.9729** | `data/45_pool_enumerator.json` | 0.3357 |
| `cap_gamma` | 4/γ₁ | **784.965** | `data/45_pool_enumerator.json` | — |
| `A_dev` | 87γ₁⁴M₀ dev% | **0.404** | `data/45_pool_enumerator.json` | — |
| `M0_eV` | M₀（由 A 反解） | **2.79009e+07** | `data/45_pool_enumerator.json` | — |
| `density_F1` | F1 比值池 coverage（log-uniform） | **20.54** | `data/48_handoff_adjudication_round5.json` | 20.61 |
| `density_legacy` | legacy family density（不可復現） | **98.66** | `data/45_pool_enumerator.json` | 97.9 |
| `pnull_flavour` | 108 模型掃描 p_null（新錨） | **0.565** | `data/47_downstream_rerun_round4.json` | 0.550 |
| `KN_sigma` | K/N σ_Δ（協方差） | **0.2327** | `data/47_downstream_rerun_round4.json` | 0.164 |

**分類規則（R5 ledger 用）**：
- **A 類 — 上游進口**（O120／O140／O142／M1／M10／FLAG／BMW／PDG）：不可由本 O data 追，引用須標來源。
- **B 類 — 計數整數**（`108`、`504`、`1892`、`1936` 等）：由列舉定義決定，不需 data 檔。
- **C 類 — 衍生算術**（同一行內由其他數字即時算出，如 `√18.9723 = 4.356`）：以公式可復現者合法。
- **D 類 — 真斷鏈**：以上皆非 ⟹ 必須補 producer 或改寫（本輪 **0 個**）。

## 四、執行結果

| 檢查 | 結果 |
|:--|:--|
| R2 canonical 存在性 | **0 缺失**（20/20 追到） |
| R3 舊值漂移 | **0 違規**（18 組舊值共出現 100 處，全部在已標記取代關係之文檔內） |
| R4a 腳本↔doc | **0 缺口**（`31b→31`、`35→34`、`36b→36`、`43→42`、`44→11`、`46/47/48/49` 自對） |
| R4b 敘述型 docs | `32`／`37` 兩份（無腳本，內容為分析與收尾清點；已於 41 申報） |
| R4c data↔producer | **41/41** 有 producer |
| R5 ledger | 225 token 未由本 O data 追到（全屬 A／B／C 類；D 類 0） |

## 五、結論報告

1. ✅ **全鏈齊備**：49 支 NN 腳本 ↔ 45 份 docs ↔ 41 個 data JSON 三向對齊，無孤兒。
2. ✅ **零漂移**：20 個 canonical 量全部可追；18 組舊值（`Σ=1.498`、`σ_R=0.039896/0.039714`、`X=2.56±0.094`、`z=−15.98`、`97.9%`、`1.97%`、`5.59%`、`17 identity`、`a≤0.3357`、`χ²=18.2207`、`p=3.97e−4`、`26.68`、`p_null=0.550`、`K/N σ=0.164`、`20.61%`、`0.78`、`12.0049`、`2.5600`）**全部**已落在標記取代關係之文檔內，無一被當現行用。
3. ✅ **引用規則（本 doc 立為硬規則）**：任何引用上述量，一律用 registry 之 canonical 值；要引舊值必須**同一文檔內**附「已被取代」之說明（現有 18 組皆已如此）。
4. ⚠️ **R5 ledger 之 225 token** 不作 FAIL，但 §三 之分類規則已寫死：日後若出現 D 類（既非上游、非計數、非衍生）才算斷鏈。

## 六、腳本↔文檔↔產物清單（本輪實測）

- **腳本 51**：NN `01–49`（含 `31b`、`36b`）＋ `l2_independent_verify`、`l3_doc_crosscheck`
- **docs 45**：`01–49` 中 43 份 ＋ `32`／`37` 兩份敘述型
- **data 41**：每份皆有 producer（`R4c` 0 孤兒）；`41_chain_audit.json` 為審計自身輸出（已申報）


## 附註（2026-11-23；引用規則加入 `k`，append-only）

`docs/51` 立：任何密度／覆蓋率引用一律 **`(k, 測度, 域, N, δ)`** 五元組；因為 `coverage(8, uniform) ≈ coverage(4, log-uniform 池內)`（34.229% vs 35.155%，差 0.926pp）⟹ 單報支撐數或單報測度皆可製造「等效」錯覺。
本 registry 之 `density_F1`／`density_legacy` 兩項亦按此格式引用（其 `k` 為 `0`／`族級`，測度為 log-uniform／legacy）。
