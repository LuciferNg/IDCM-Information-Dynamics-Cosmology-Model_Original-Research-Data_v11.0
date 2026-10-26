# M6 — 尺度自由度定理（兩界）：`dof_scale = 1` **恰**（條件性定理）

## 狀態: 🟢 closed（2026-10-26；第一輪 40/40 ＋ 第二輪 13/13 ＋ 自審 4/4 ＋ 四路 ①②④ 完成 ＋ 全鏈審計 37/37 ＋ M8 T3 回交已收訖）

> **地位演進**：原 T2 = 定義性歸約（`docs/08` 自審）→ **T2′ = 定理**（`docs/10` 上界 ＋ `docs/12` §2 下界 ＋ `docs/11` β no-go）。
> 假設明列：H1 Galois 協變、H2 質量繼承晶格關係（`docs/10` §6）；界外：離散 RG 類比（`docs/11` §4）。

問「絕對標度係唔係獨立結構對象」。命題：原子以上皆複合投影（質量加性）⟹ 中間層標度只係計數／插值；全域尺度自由度 = 1（O212 錨）；剩餘唯一非平凡數 = 純計數 `ω_M0/ω_site = 4.095e38`（T3 → M8）。

**第一輪結果**：T1 加性 ✅ · T2 `dof_scale = 1` ✅ · P4 反例 0 ✅（40/40）。

**第二輪（邊界 item 2–4）**：
- item 2 ｜p（交換機率）= **真獨立 dimensionless dof**（dim p = 0、Jacobian det = 1）；速率預言覆蓋 100% ⟹ 零鑑別力；唯一可證偽點 `r > 1.0397344 Gyr⁻¹`。
- item 3 ｜κ 之 λ：**晶格層不可識別**（單 mode 0 命中、2-mode 10 個、全支撐 1,817 個）。
- item 4 ｜中間層隱藏 **機制已定**：支撐定理 `λ_max(k) = 4k − Σ_{j<k}γ_j` ⟹ 身份清晰區止於 λ < 8；達 M_Pl 需 4.86e18 量子（乘法只需 10.56 層 ⟹ 與 T1 不相容）。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + T1–T3 + 攻擊路線 R1–R4 + 推導時序索引（§八–§十） | 🟢 |
| `SESSION_START.md` | 起手式 + 接手起點 + 禁區 | 🟢 |
| `final_result.md` | 結論摘要（RESULT_TEMPLATE） | 🟢 |
| `docs/01_t1_composite_additivity.md` | R1：複合投影算子 ＋ 加性（ℚ(ζ₈₈) 精確） | 🟢 20/20 |
| `docs/02_t2_scale_dof_theorem.md` | R2：dof 定理 ＋ 中間層 ＋ 第二標度封死 | 🟢 14/14 |
| `docs/03_r4_counterexample_hunt.md` | R4：10 候選反例搜索 ＋ 子集和鑑別力 | 🟢 6/6 |
| `docs/04_t3_handoff_and_verdict.md` | T3 交接 ＋ M5→M6 答覆 ＋ 總判決 | 🟡 |
| `docs/05_p_dof_ledger.md` | item 2：p 之自由度總帳（量綱／Jacobian／預言力） | 🟢 5/5 |
| `docs/06_kappa_lambda_identifiability.md` | item 3：κ 之 λ 可識別性（精確 2⁴⁴ 計數） | 🟢 4/4 |
| `docs/07_hidden_layer_mechanism.md` | item 4：中間層隱藏機制（支撐定理＋加法量子化） | 🟢 4/4 |
| `docs/08_framework_self_audit.md` | 框架自審：T1/T2 鑑別力測試（4/4 = 零鑑別力之證） | 🟡 |
| `docs/09_upgrade_routes.md` | 四條升級路線 ＋ Galois 軌道（44 → 7） | 🟢 5/5 |
| `docs/10_forced_single_scale.md` | **T2′ 上界**：7 軌道標度被迫全等（精確代數） | 🟢 9/9 |
| `docs/11_beta_no_go.md` | **標度生成 no-go**：無連續 β ⟹ 無 transmutation | 🟢 7/7 |
| `docs/12_economy_and_homogeneity.md` | 參數經濟（速率軸 ≤6）＋ **齊次性 no-go**（dof=0 不可能） | 🟢 4/4 |
| `scripts/01–07_*.py` | 推導腳本（精確算術） | 🟢 40/40 ＋ 13/13 |
| `scripts/08`–`11_*.py` | 自審／對稱／強制單一標度／β＋齊次性 | 🟢 25/25 |
| `scripts/12_audit_chain.py` | 全鏈審計（十項，唯讀） | 🟢 30/30 |
| `scripts/04_doc_consistency.py` | 文檔↔數據一致性（32 claims） | 🟢 32/32 |
| `scripts/run_chain.sh` | 全鏈 driver（預設 dry-run，`--run` 才跑） | 🟢 |
| `data/01–12_*.json`、`data/run_*.log` | 輸出數據與執行記錄（重跑 byte-identical） | 🟢 |

## 依賴

O38（docs/15 單位制測試；docs/28 穹頂判據；docs/38 G1 關閉／G2 保留）· O212（錨定理 dof=1）· M2（標度盲引理、質量泛函唯一性）· O211（Σγ=90）· O124（β=0 ⟹ 無 transmutation）· O213/O121（宏觀端錨）

## 註冊

pid `M6`（created 2026-10-25；closed 2026-10-26，verdict = closed）

## 框架地位（2026-10-26 自審後裁定）

> 詳見 `docs/08_framework_self_audit.md`（`scripts/08`，4/4 PASS = 零鑑別力之證明）

| 部份 | 地位 |
|:--|:--|
| T1 之加性檢查（A6–A8） | 🔴 **定義性**：換隨機數據（含正負數）同樣 0 違反 ⟹ 只驗分配律 |
| T2 `dof_scale = 1`（原版） | 🔴 定義性歸約（`docs/08`）→ **已由 T2′ 取代**（下表） |
| **T2′ 上界（`docs/10`）** | ✅ **定理**：Galois 協變 ＋ 關係保持 ⟹ 7 軌道標度全等（500×7 系統解空間 = 1）；反例案（場值標度）6 自由 |
| **T2′ 下界（`docs/12` §2）** | ✅ **齊次性定理**：任何晶格泛涵皆 degree-1 齊次 ⟹ 整體倍率恆自由 ⟹ `dof ≥ 1`；**升格 O38 G1**（偽命題 → 結構不可能） |
| **標度生成 no-go（`docs/11`）** | ✅ **定理**：`3g = 0` ＋ 繞數 running ⟹ 無連續 β ⟹ 無 transmutation ⟹ 無誘發第二標度（界外：離散 RG） |
| 88 邊形代數事實（Σγ=90、配對和、rank 20、**7 個 Galois 軌道**、支撐定理、tan 族度數） | ✅ **真數學**（但係數學事實，未係物理證據） |
| 不可識別性結果（R4 子集和、item 3 κ、item 2 p） | ✅ **真認識論內容**（負面結果有實效） |
| **速率軸自由度** | 🟡 **修正**：≤ 6（依 7 個 Galois 軌道，正規化後）；均勻假設下 1（`docs/05` 附註更正） |
| 「質量 ∈ 20 維 ℚ 空間 ⟹ 21 質量整數關係」 | 🟡 **唯一非定義性物理預言**，暫不可判（需 21 個 @ δ≤1e−10；現有 4 個） |
| 加性 vs M5 束縛能 | 🔴 **未對賬**（H2 之檢驗，仍未做） |

⟹ M6 由「分類」升級為**證明為何恰一個標度**：上界殺 7 軌道方案、下界殺 0、β no-go 殺後天生成。
仍未閉：**H2（強子層加性）未檢驗**；速率軸 ≤6 未定。