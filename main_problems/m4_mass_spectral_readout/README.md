# M4 — 質量 = 譜讀數差值：全費米子譜差閉合

## 狀態: 🟡 partial（2026-10-24：R1 閘通過、R2 橋精確、R3 lattice 重判、R4 A negative；未 close）

質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2, R7 定案承接）。詳見 `01_problem_record.md`／`SESSION_START.md`／`docs/06_r_summary.md`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 已確立結構 + §八 R1–R4 結果 + §九 append-only 更正 | 🟡 partial |
| `SESSION_START.md` | 起手式 + 接手起點 + §七 2026-10-24 執行記錄 | 🟡 partial |
| `docs/` | 推導記錄（6 檔） | 🟢 |
| `scripts/` | 計算腳本（12 檔：01–08 遷入、09–12 本輪新增） | 🟢 |
| `data/` | stdout log ＋ JSON（17 檔） | 🟢 |

## 文檔索引

| 檔 | 內容 |
|:--|:--|
| `docs/01_r1_chain_rebuild.md` | R1 閘 8/8；01 [C] TypeError 修復；05 MC protocol 不一致（8.7% 不可復現） |
| `docs/02_density_unified.md` | 統一 density protocol（thr = 命中自身 dev）＋ 夸克/輕子/hadron 全滅判決 |
| `docs/03_gamma1_gamma6_bridge.md` | R2：γ_n = 2−2T_n(c)、遞推、γ₆/γ₁ deg-20 本原元；O89 三階梯作廢核證 |
| `docs/04_lattice_r3.md` | R3：FLAG 2024 實檔值；span≈12 升 🟡（預言 m_u=2.1498 MeV）；m_d/m_u=11/5 候選 |
| `docs/05_A_identity_r4.md` | R4：A = 1.64337 eV/γ²（v_SM）；族 13104 組合 0 命中；M_comp 反解無身份 |
| `docs/06_r_summary.md` | 16 項 derived-vs-not ＋ 結構化 gap（G1–G6） |

## 依賴

O257 (R7 docs/07 定案) + O66 (M_comp 降級) + O89 (γ₁ 冪階梯, **2026-10-24 核證：M_GUT/PS/R 三階梯作廢**) + O57 (m_ν 盲測) + P6 (§10 差值閉合) + O73 (c=Goldstone) + O212 (dof=1)

## 註冊

pid `M4`（created 2026-10-21；state → partial 2026-10-24）
