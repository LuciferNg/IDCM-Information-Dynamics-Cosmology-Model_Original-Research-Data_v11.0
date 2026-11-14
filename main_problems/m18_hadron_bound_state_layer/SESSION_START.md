# M18 — SESSION_START（接手文件）

> **建立**: 2026-11-14｜**狀態**: ⚪ 已註冊、未開推
> **本 session 性質**: 只開不推 —— 零 scripts、零 docs 內容、零數據。

---

## §一 接手起點（第一步讀邊幾份）

依序讀，唔好跳：

1. 本檔（章程 + 禁區 + 重現閘）
2. `01_problem_record.md`（§〇 歸屬 fence、§一 Q1–Q4、§二 T1–T4、§五 R1–R5）
3. `README.md`（依賴與註冊）
4. M4 `main_problems/m4_mass_spectral_readout/§11.5`（缺口原文）
5. M5 `main_problems/m5_hadron_bound_state_spectrum/final_result.md` §未閉合
6. M7 `main_problems/m7_hadron_model_basis/final_result.md`（模基退化）
7. M6 `main_problems/m6_scale_dof_theorem/final_result.md`（dof = 1 定理）

---

## §二 重現閘（硬性）

本 M 係新開，**冇自家 scripts 可重現**。第一輪工作嘅閘係：

```
① 先由「零 import 上游腳本」重建 §三 所列輸入值（γ_n 譜、γ_n+γ_{44−n}=4、c=cos(π/44)）
② 逐項對照 M4／M5／M6 文檔記載值；數唔對 ⟹ 停，唔好開始 R1
③ 唔准重跑 M4/M5/M6/M7 嘅腳本（會寫入他 M 目錄）
```

**期望對照（上游記載，唔係本 M 產物）**：
- `γ₁ = 0.005095770779492917`（M4 `docs/12`）
- `γ₄₄ = 4`，`1/γ₁ = 196.24116611059782`
- `γ₂ = 4γ₁ − γ₁²` EXACT（ℤ[γ₁] 封閉性例證）

---

## §三 禁區（承接 session 不得重推）

```
✗ ⟨O⟩ = w₁ − w₄                 M7 closed(negative)；⟨O⟩ ≡ 1 係歸一化恆等
✗ contact 算子標度／支撐／衰減律  M5 三層掃描全敗（缺口 62–68×）
✗ 以觀測強子質量反解 mode 佔用    後驗，違零參數鐵律
✗ φ-冪比值                      M9 域分離定理
```

**唔可以由「掃描」起步**：R1–R4 每條都要先寫死判準，再算數。掃完才定判準 = 後驗。

---

## §四 與相鄰 M 之分工（防重複）

| M | 做咩 | M18 唔做咩 |
|:-:|:--|:--|
| M5 | 已有譜下之 `⟨1/r⟩` 矩陣元、質量差 | 唔重算 `A_max`／`A_req` |
| M7 | 模基 `w_k(h)`（已判退化） | 唔重用退化模基作起點 |
| M4 | 譜代數、差值力學、`Q1` 活預言 | 唔重推 ℤ[γ₁] |
| M6 | 尺度 dof 定理 | 唔重推 dof 計數 |
| M12 | CY₃ 投影層 | 唔碰味 sector |

---

## §五 交接／依賴

- **上游（提供輸入，唔受本 M 改動）**：M4、M5、M6、M7、M12
- **下游（等本 M 解除阻塞）**：M4 `Q2/Q3/Q8`、M5 質量差 em 閉合
- **跨 M 交接只限 register**（`register.py handoff --target-pid ...`）；禁在其他 M 寫入任何檔案。
- **版本隔離**：本 M 屬 v11.0。禁引用 v10.0 數字作依據（歷史記錄可提及，唔可引用）。

---

## §六 本 session 交付記錄

```
2026-11-14  開檔：01_problem_record.md + SESSION_START.md + README.md
            資料夾 m18_hadron_bound_state_layer/{docs,scripts,data} 建空
            零推導、零 scripts、零 data
```
