# M21 docs/05 — 雙態登記、產物查重與總判詞（P5／Q3 ＋ NDR1）

**腳本**: `scripts/05_registry_and_dedup.py`（log `data/run_05.log`，**18/18 PASS**）
**判準**: `docs/00` §3–§4（`NDR1` 於本輪新增，屬登記制度層，非假設層）
**日期**: 2026-11-20

---

## §1 產物查重（跨 O／跨版本，唯讀）

掃 `v10 ∪ v11` 之 `.md/.py`（排除本 O）。
**快照 = 2026-11-20 第一輪：9,576 檔**（下表數字屬該快照）；
**第三輪重跑（2026-11-20 晚）：9,707 檔** —— corpus 隨他人新增檔增長，命中數隨之（K1 23→26、K6 92→93、K7 151→157、K8 1→2），
**結論方向不變且更強**（各關鍵詞族仍 ≥1 命中）。活值一律以 `data/05_registry_Nd.json` 為準。

| 族 | 命中檔數 | 判定 |
|:--|:-:|:--|
| `K1  Σ_χ d(d−1)` | 23 | 已註冊（O171 docs/33、docs/38；M15 §三；M16 記賬） |
| `K2  N_d = B·L₅ = L₃·L₅` | 167 | 已註冊（M15 docs/13 Q3；O171 `2B·L₅ = 88`） |
| `K3  Σγ = 2N+2` | 77 | 已註冊（O171「三恆等式 ⟹ 鑑別力 0」） |
| `K4  K_tree = 33 = 3·11` | 29 | 已註冊（v10 O83／O218） |
| `K6  f_M0·f_R / N_O38` | 92 | 已註冊（M8 docs/06） |
| `K7  60 = N_d + B²` | 151 | 已註冊（M16） |
| `K8  度集唯一性引理` | 1（標題字樣） | **未見同型** ⟹ 本輪 `CA1` 視為新 |

⟹ **本 O 不宣稱任何新恆等式**。本輪新增限於：
`CA1` 唯一性引理 · `CA2` 普查形式（加權均值 = B）· `CA5/CA6` 密度閘數字 ·
`BR1` 定量（249 種分解）· `GR1` 邊分類（T/M/A 三型）· dof 分支判定。

## §2 雙態登記（機械生成：`data/05_registry_Nd.json`）

```json
{
  "quantity": "N_d = 44",
  "formal_state": "structure-derived (group side) / bridge-open (lattice side)",
  "state_zh": "結構導出（群側）／橋未閉（晶格側）",
  "dof_continuous_added": 0,
  "producers": ["E4（Σd(d−1) = |G| − N₂，dof=0，F）", "E9（B·L₅ = 4×11，F，給定 CA2）"],
  "tautologies_never_producers": ["E5", "E6"],
  "matches_no_operator_identity": ["E7 h11", "E8 h21", "E13 N₂ = B²"],
  "refuted_routes": ["p9-FORCE44（φ-級聯，M15 docs/14 判 🔴）"],
  "open_bridge": "A₅ 44 非對合 ↔ 鏈 44 模式之典範映射（≙ O171「為何乘 2 ＋ 為何環」）"
}
```

**登記制度（Q3／P5 之答覆）**：

```
N_d = 44 之後不得再寫「錨／輸入」—— 應寫：
  「N_d = 44 結構導出（群側，dof = 0）；身份橋未閉（晶格側，BR1 只到存在）」
任何引用必須帶 (producer 碼, 橋狀態) 二欄；只引數值不引橋者，一律標 🟡。
```

⟹ M15 `docs/14` 之「錨」判定**降級**為「群側已導出、晶格側橋未閉」（**非**推翻：其 FORCE_44 判詞不變）。

## §3 總判詞表（Q1–Q4／T1–T4）

| 項 | 問題 | 判 | 依據 |
|:-:|:--|:-:|:--|
| **Q1／T1** | `N_d = 44` 有無 producer？抑或輸入？ | 🟡 | 群側 ✅（`E4`：`60 − 16`，dof=0；度集被 `CA1` 唯一強迫）；晶格側橋 🟡（`BR1` 只到存在） |
| **Q2／T3** | `f_M0`、`f_R` 之結構身份 | 🔴 | `TF1`：乘積 = `N_O38` = `M₀/ħω_ring` 單一 object；分解係 gauge（`docs/03`） |
| **Q3／P5** | 若為輸入如何標示 | ✅ | 雙態登記制度（§2）；`data/05_registry_Nd.json` 機械生成 |
| **Q4／T4** | 互導圖是否閉合 | 🔴（不閉） | `F7 / T2 / M3 / A0 / D1`；無孤立節點；缺口 = 3 條 M ＋ 2 條 T（`docs/04`） |
| T2 | CY₃ 路線 | 🔴 | k-free 命中 6 ≥ 5；`q_obs` 53.5%（`docs/02`） |

## §4 對上游之影響（只 register，禁跨版本推送）

```
v11 M15   Q3（K_tree = 33 歸約）→ 本輪提供 K_tree = Σ(d−1)² = 33（E10）
          與已註冊 K_tree = N_gen·L₅ = 3·11（v10 O83）之關係 = 一致性，非新預言
v11 M6    殘餘措辭「2（f_M0·f_R）」→ 應讀作「乘積 1 object ＋ N_IDCM」；dof = 1 不變（分支 (a)）
v10 O171  residual「為何乘 2 ＋ 為何環」= 本輪 BR1 所定之橋 ⟹ **v11 讀 v10 結論，不反向推送**
```

## §5 交付清單

```
scripts/00_reproduction_gate.py      22/22 PASS  重現閘
scripts/01_a5_census_identity.py     26/26 PASS  R1／T1
scripts/02_cy3_pool_gate.py          4/6   PASS  R2／T2（2 FAIL = 結論）
scripts/03_two_factor_audit.py       8/8   PASS  R3／T3
scripts/04_derivation_graph.py       9/11  PASS  R4／T4（2 FAIL = 結論）
scripts/05_registry_and_dedup.py     18/18 PASS  P5／Q3
data/*.json（6 檔）＋ data/run_0*.log（6 檔）
docs/00_preregistration.md ～ docs/05（本檔）
```
