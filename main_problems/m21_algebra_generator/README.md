# M21 — 代數產生器層（algebra generator）

## 狀態

🟡 **partial — 部分閉合**（2026-11-20 第一輪）｜本 O 自身殘項 **3**（橋選定、h11/h21 算子關係、`N₂ = B²` 身份）
📎 交付：`scripts/00–06`（7 支）｜`docs/00–05`（6 份）｜`data/` 7 JSON ＋ 7 log｜`final_result.md`

## 一句話

`N_d = 44` **唔係天上掉下來** —— 佢係 A₅ 三條群不變量之零參數後果
`Σ_χ d(d−1) = |G| − #{g: g²=e} = 60 − 16 = 44`，而度多重集 `{1,3,3,4,5}` 被
`(Σd²=60, Σd=16, k=5)` **唯一強迫**（新引理 `CA1`）。未閉的只剩
「鏈 44 模式 ↔ 群 44 非對合」之典範映射 —— 即 `O171` 早已列為 residual 的「為何乘 2 ＋ 為何環」。

## 範圍收窄（照章程）

```
✓ 本 M 只攻：N_d = 44 之 producer ＋ 10^38 = f_M0·f_R/N_d 之兩因子缺口
✗ 唔重推 B = 4        M16 已閉合（dof = 0）
✗ 唔重推 Σγ = 90      M6 環核 12/12
✗ 唔重做 FORCE_44 審計  M15 out_of_scope/docs/14 已完成（只引用）
✗ 唔碰 v10 p9 內部推導  唯讀引用
```

## 判詞

| 項 | 判 | 一句 |
|:-:|:-:|:--|
| Q1／T1 | 🟡 | 群側 ✅（`E4`，dof=0）／晶格側橋 🟡（`BR1` 只到存在，249 種分解） |
| Q2／T3 | 🔴 | 「兩因子」係 gauge 分解；乘積 = `N_O38` = `M₀/ħω_ring` 單一 object |
| Q3／P5 | ✅ | 雙態登記制度（`data/05_registry_Nd.json` 機械生成） |
| Q4／T4 | 🔴 | 13 邊：`F7/T2/M3/A0/D1`；缺口 3 M ＋ 2 T；無孤立節點 |
| T2 | 🔴 | CY₃ 池：k-free 命中 6 ≥ 5、`q_obs` 53.5% ⟹ 零判別力 |

## 現況表

| 量 | 值 | 現況 | 出處 |
|:--|:-:|:--|:--|
| `B` | 4 | 🟢 閉合 | M16 |
| `Σγ` | 90 | 🟢 環核（**T 型**：恆等反演，非 producer） | M6 ＋ M21 `CG1` |
| `h11`/`h21` | 36/98 | 🟢 上游（`N_d` 關係為 **M 型**，無算子身份） | v10 O211 / M21 `E7,E8` |
| **`N_d`** | **44** | 🟡 群側導出（`F`）＋ 橋未閉 | M21 `E4`、`CA1` |
| `K_tree` | 33 | 🟢 `= Σ(d−1)² = 3·11`（與 v10 O83 一致） | M21 `E10` |
| `10^38` | 4.095e38 | 🔴 「兩因子」偽命題 → 單一 object | M21 `TF1` |

## 缺口數字（不變）

```
圖卡跨距 2.894296 dex   vs   T₁ 需 38.612274 dex
log₁₀T₁ = log₁₀(M₀/ħω_ring) − log₁₀(N_d) = 40.2557 − 1.6435 = 38.6123
gap = 38.612274 − 2.894296 = 35.717977 dex = 5.223686e35×   （「缺口」無新內容）
```

## 檔案

| 檔 | 內容 |
|:--|:--|
| `SESSION_START.md` | 接手起點、重現閘、範圍收窄、禁區、分工 |
| `01_problem_record.md` | §〇 fence／Q1–Q4／T1–T4／P1–P5／R1–R5／禁區／**§十 推導結果** |
| `docs/00_preregistration.md` | 判準凍結（跑前寫死）＋ POSTHOC 披露 |
| `docs/01_a5_census_identity.md` | R1／T1：A₅ 普查恆等式（26/26） |
| `docs/02_cy3_route_and_pool_gate.md` | R2／T2：池密度閘（🔴） |
| `docs/03_two_factor_audit.md` | R3／T3：兩因子 gauge（🔴） |
| `docs/04_derivation_graph.md` | R4／T4：互導圖 `CG1/GR1` |
| `docs/05_registry_and_dedup.md` | P5／Q3：查重 ＋ 雙態登記 ＋ 總判詞 |
| `final_result.md` | 最終摘要 |
| `scripts/00–06` | 重現閘 → R1–R4 → 登記查重 → 自查水位 |

## 依賴

```
上游：M16（B=4）· M6（Σγ=90）· M8（docs/04 §4、docs/06 E2/E3）· M15（§三、docs/13 Q3、
      out_of_scope/docs/14）· v10 O171（N_d = 2Σd(d−1)）· O211（h11/h21）· O83／O218（K_tree=3·11）
下游：M15 Q3 · M6 dof 記帳（分支 (a)：O212 dof=1 成立）· M3（Ω 分母 60）
```

## 註冊

`register.py` pid = `M21`｜state = `partial`｜verdict = `partial`｜
note → M15（Q3 提供 `K_tree = Σ(d−1)²`）／M6（殘餘措辭、dof=1）
