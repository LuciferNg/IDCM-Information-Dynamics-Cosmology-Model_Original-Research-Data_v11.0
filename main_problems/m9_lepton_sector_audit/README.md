<!-- registry
pid:     M9
title:   輕子扇區獨立攻：比值算子身份唯一化（譜讀數 vs φ-冪）
state:   closed
verdict: negative
flags:   sigma-not-percent; phi-lock-gate; primary-source-audit; density-gate; no-lepton-ratio-structure
owner:   M9
created: 2026-10-27
updated: 2026-10-27
prereq:  M4(docs/04,19), O34, O35, O120, O209, M6
source:  調律者 2026-10-27「開M9 + M10，但本session不直接推」
-->

# M9 — 輕子扇區獨立攻：比值算子身份唯一化

## 狀態: 🔴 closed（negative，2026-10-27）｜P1–P5 完成｜scripts 01–04 **29/29 PASS** ＋ `scripts/05` 交叉核 **43/43 PASS**

**問題**：同一個觀測量 `m_μ/m_e = 206.7682827`，三個互斥聲稱同時掛住 🟢。

## 三聲稱判決（統一 σ 門檻）

| 公式 | 值 | σ（u = 4.6e−6） | 判 |
|:--|--:|--:|:--|
| `φ^(11+7/88)`（O34/O35） | 206.770253211 | **428.4σ** | 🔴 後驗 fit ＋ 數值否證 ＋ K_φ FAIL |
| `(γ₇/γ₁)(F_μ/F_e)`（memory／M4 L58） | 206.705 記載 | 恆等 0σ | 🔴 **1 dof tautology**（`F_μ/F_e ≡ 目標/(γ₇/γ₁)`）；記載值不可重現 |
| 整數 `207`／「界內可達」（M4 docs/19） | 207 | **50 373σ** | 🔴 deg-0 空洞；「可達」非預言 |

## 病因

三個 🟢 能並存**唔係理論問題，係門檻問題**：全部用 % 級（0.0005–0.11 %）門檻。
輕子權威相對精度 `2.225e−8` ⟹ **% 門檻粗 4–6 個數量級**。

## 判準（本 M 唯一新東西）

同一量**唔可能**有兩個合法身份；判準唔係「邊個準」：

1. **事前性**：參數須在見目標前凍結（φ-lock）。
2. **可證偽性／族密度**：族必須有能力落空。φ-ladder 係**插值族**（528 半格 0.0456 % ⟹ 任何目標都「命中」⟹ 每次檢定 0 bit）；γ-ratio 係**離散族**（1 892 值，可落空）。
3. **dof 成本**：每反解因子／每取整 = 1 dof（L2 雙 γ 比 P_MC(1e−5) = 52 %）。

## 主要結果

- `CAP = γ₄₄/γ₁ = 784.9646644`；`m_τ/m_e = 3477.23 = 4.43× CAP` ⟹ 跨代輕子比**結構上**非單模。
- L0 譜語法（零 dof）：μ/e 最近 1.1524 %、τ/μ 1.0396 %、τ/e 77.4 % ⟹ 全落空。
- O209 唯一分配之純 γ-ratio 讀值錯 46×–6.1e4× ⟹ 輕子質量唔可能係 γ 比值。
- holdout `m_p/m_e`、`m_n/m_e` 同樣 −57 % ⟹ 譜比語言整體無法達 σ 級。
- 上游 `o34/scripts/21` primary source：先算目標後硬編 `(11,7)`、1800 點網格反搜、反解打印、`m_e = 0.511` —— 記載 `0.0003 %`／`0.0015 %` 由此圓整重現，真值 `0.000953 %`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:--|:--|:--:|
| `01_problem_record.md` | 問題 ＋ fence ＋ 禁區 ＋ 判據 ＋ 執行記錄 | ✅ |
| `SESSION_START.md` | 接手起點（已更新為已執行） | ✅ |
| `docs/01_sigma_table.md` | P1/P2 權威精度 ＋ σ 表 | ✅ |
| `docs/02_phi_lock_gate.md` | R2 K_φ ＋ 密度 ＋ primary-source 審計 | ✅ |
| `docs/03_language_consistency.md` | R3 可證偽性判準（含 Q3 更正） | ✅ |
| `docs/04_spectral_uniqueness.md` | R4 語法 ＋ dof 階梯 ＋ holdout | ✅ |
| `docs/05_preregistered_predictions.md` | R5 預註冊預言 | ✅ |
| `scripts/01…05` | 零 import 獨立重建（29/29）＋ 文檔×數據交叉核（43/43） | ✅ |
| `data/run_01…04.log` ＋ 4 JSON | 執行記錄 | ✅ |
| `final_result.md` | 結論摘要 | ✅ |

## 依賴

M4（docs/04,19）｜O34（scripts/21 審計）｜O35（docs/02 階梯）｜O209（Sel(𝒮)、分配）｜O120（C_f dof=0）｜M6（dof=1）｜O38 scripts/33（φ-lock 慣例）｜M8／O38 G1（絕對標度，未碰）

## 註冊

register pid **M9**｜2026-10-27｜`close --verdict negative`｜本 M 未寫入任何其他 M/O 目錄