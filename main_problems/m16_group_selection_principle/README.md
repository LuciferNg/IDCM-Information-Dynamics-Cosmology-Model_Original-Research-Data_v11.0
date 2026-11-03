# M16 — 群選擇原理（Group Selection Principle）

**狀態**: ⚪ 已註冊、未開推 ｜ **version** v11.0 ｜ **性質** 公理層缺口

---

## 定位

`M15` 已把宇宙預算分母由「凑數」升級為「A₅ 之 irrep 恆等式」：

```
Σd² = Σd + Σd(d−1)   ⟺   60 = 16 + 44      （d = [1,3,3,4,5]）
```

但下一步即刻撞牆：**為何係 A₅？** 呢個唔係算術問題，係**群選擇原理**問題 ——
本 O 專門處理呢一條缺口，唔做其他。

| 已閉（勿重推） | 出處 |
|:--|:--|
| `60 = 16 + 44`（A₅ irrep 恆等式） | `M15 docs/25` |
| 容量 `= (B−1)/60`、溢出不可能 | `M15 docs/17`、`docs/23`、`docs/24` |
| `3 = B − 1`（階 4 群非平凡特徵標） | `M15 docs/17` |
| A₅ ↔ 88-site 配分（含舊群選擇準則） | v10 `O171`（只讀） |

## 本 O 之五個缺口

| # | 缺口 |
|:-:|:--|
| **R1** | 準則兩條件（`∏ = |G|`、`Σ = 12`）各自之結構理由 |
| **R2** | 全域唯一性：階 `≤ 500` 之群，同時滿足者是否只有 A₅ |
| **R3** | 「相異」「非平凡」兩修飾詞之必要性（自由度計數） |
| **R4** | 準則之地位：導出 / 事後（look-elsewhere 校正） |
| **R5** | `60 = |A₅|` 與 `Σd² = |G|` 是否同一陳述（若同一，`M3` 之 🟢 須降級） |

## 三件套

```
01_problem_record.md   ★ 主檔（問題、已立、缺口、判準、範圍、參考集）
SESSION_START.md       承接指引（已立事實、陷阱、首步）
README.md              本檔（定位）
```

## 目錄

```
m16_group_selection_principle/
├── 01_problem_record.md
├── SESSION_START.md
├── README.md
├── docs/        ⬜ 待承接 session 建
└── scripts/     ⬜ 待承接 session 建
```

## 參考集（摘要；全文見 `01_problem_record.md §6`）

**v10 主先前工作**

- `o171_a5_partition_theory/` —— A₅ ↔ 88-site 配分；`final_result.md:225` 群選擇準則；
  `docs/04_r1_sector_selection_law.md:67` 準則原文
- `o137_cy3_topological_resonance/` —— `docs/24:259` 明寫「為何 A₅ 集被選」為下一步
- `o19_lowl_dodecahedron/`、`o102_12_pentagons_su12/`、`o238_mass_as_cost/`、
  `o167_cmb_h0_tension/`、`o213_world_dynamics/`、`o220_spin_origin/`、`p4_cmb_cold_spot/`

**v11**

- `m15/docs/25_a5_sixty_uniqueness.md`（★ 來源）、`m15/docs/02`＋`scripts/02`（A₅ 諧波重數）
- `m3/01_problem_record.md:37`（`60 = |A₅|` 🟢）、`m12/SESSION_START.md:36-37`（A₅ 頂點標籤 🔴）
- `m4/docs/{30,33,34,37,47}`、`m7/docs/{17,19}`（A₅ 引用，待逐條核）

## 使用方式

1. 讀 `SESSION_START.md` → 已立事實與陷阱；
2. 讀 `01_problem_record.md §3`（五缺口）與 `§4`（判準）；
3. 自行判斷從哪一缺口入手（建議 `R5 → R2 → R3 → R1`）；
4. **判死亦算結果** —— 若群選擇不可由結構決定，明寫「A₅ = 輸入」即可結案。