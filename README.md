# IDCM-v11.0 -- Reality-Architecture

> **Project**: 公式驗證 —— 所有 IDCM 公式必須通過「第一性原理推導」審計
> **鐵律**: 公式不能以現有數值來後驗擬合 (no post-hoc fitting to measured values)
> **創建**: 2026-09-03 | **狀態**: 🔴 框架階段, 逐 q 處理
> **位置**: `/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/`

---

> **命名記錄 (2026-10-14, 調律者)**: 本 repo 名稱由 `IDCM-v11.0--Equation-Verification` 改為
> `IDCM-v11.0--Reality-Architecture`。方向以 Reality-Architecture 為準；公式驗證（緣起目標）繼續為
> 支柱之一。O115 docs/03 曾訂名「It-From-Selection」——與實際開檔名唔一致，以 Reality-Architecture 為最終名。

---

## 緣起

v10.0 嘅公式 (Ω_Λ=4/(φ⁴−1)、Ω_b、ε=φ⁻¹/4 等) 大量以「MERA 分支閉合」「轉移矩陣不變量」等標籤呈現, 但追源發現多處只有**數值巧合接軌**, 冇結構推導。本 project 逐條公式重新推導, 唔接受後驗擬合。

## 背景判決 (2026-09-03, 前 session 完成)

針對 Ω_Λ=4/(φ⁴−1) 分子 4 嘅追源:

| 項目 | 狀態 |
|:---|:---:|
| Ω_Λ=1/(1+3ε) 精確代數恆等式 (ktree_structure/08 引理3) | ✅ 結構 |
| 分子 4 = 通分產物 (唔係 B 空降) | ✅ 事實 |
| **分子 4 ≡ B** 嘅推導連接 | 🔴 **空缺, 只係數值巧合+標籤** |
| **ε=φ⁻¹/4 嘅「4」來源** | 🔴 **從未推導 (全 corpus 0 命中)** |

→ 觸發 Q1。

---

## Q-series 索引

| ID | 問題 | 狀態 |
|:---|:-----|:----:|
| ~~Q1~~ | ~~真空能 4/(φ⁴−1) 分子 4 第一性原理~~ | ~~🔴 框架已寫, 待下 session~~ |

> Q1 已刪除（2026-10-14, 調律者指示）；備份 `/tmp/q1_vacuum_energy_backup_20261014.tar.gz`（session 用, WSL 重啟會清空）。

---

## main_problems 索引

| ID | 問題 | 狀態 |
|:---|:-----|:----:|
| [M1](main_problems/m1_mass_file_schema/) | 質量檔案 schema：譜點/標度/算子/拓撲/標籤五層格式規範（Reality-Architecture; 由 v10 O251 遷移, 2026-10-14） | 📄 已註冊 — 未開始 |
| [M2](main_problems/m2_sync_desync_exchange/) | 同步態 ↔ 脫同步態互換公式（質能公式；只開不推, 2026-10-15） | 📄 已註冊 — 未開始 |
| [M3](main_problems/m3_omega_ratio_lock/) | Ω 比例硬性鎖定：Ω_b=3/60、Ω_dm=16/60、Ω_Λ=41/60（只開不推, 2026-10-15） | 📄 已註冊 — 未開始 |
| [M4](main_problems/m4_mass_spectral_readout/) | 質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2, O257 R7 定案承接；8 支遷入 scripts；只開不推, 2026-10-21） | 📄 已註冊 — 未開始 |