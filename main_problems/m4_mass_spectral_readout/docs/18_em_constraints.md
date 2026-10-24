# M4 docs/18 — em 讀數 operator 約束集 v1：8 通道實證譜，四類形式全排除 ⟹ Q2/Q3 依賴 O140

> 腳本：`scripts/18_em_constraints.py`｜log：`data/run_18_em_constraints.log`（sha256 `018f044d0972b9e8…`，49 行）｜JSON：`data/18_em_constraints.json`
> Canonical run（2026-10-24）：exit 0，49 行
> 狀態：🔴 電荷加性排除（R² 0.003／0.024）；🔴 substitution-count 線性排除（R² 0.089）；🔴 spectator 獨立部分排除（baryon spread 3.578 MeV）；🔴 標度幂律排除（worst dev 117%）；🟡 K/N 一致（spread 0.107 ⟹ H1 再證）；⟹ **Q2/Q3 依賴 O140（flavour→mode 映射）**

## 一、前設問題

```
M4 Q2/Q3 — em 讀數 operator 約束集 v1：7 通道實證譜，形式限制與依賴圖（2026-10-24）
======================================================================================
上游：docs/14（R_em = 3.479、反相=恆等還原）；docs/17（最小電荷模型結構否證 ⟹ 必須非電荷加性）；
      01_problem_record §十一 Q2/Q3；Q5（selection rule）已 handoff O140。
目標：唔再「揀模型」，而係**由數據反推 operator 必須滿足嘅形式約束**，並定出依賴圖
      （若 operator 需 flavour→mode 映射，則 Q2/Q3 形式上下游於 O140）。

Protocol（事前寫死）：
  X ≡ m_d − m_u = 2.560 MeV（FLAG24，外部固定；唔由本輪擬合）
  C1 通道 em 讀數：E_ch ≡ Δm_ch − (substitution count)·X ，substitution = u↔d 替換數
     通道（PDG 2024）：n−p、K+−K0、π+−π0、Σ−−Σ0、Σ0−Σ+、Σ−−Σ+、Ξ−−Ξ0
  C2 形式約束測試（每項事前指定判準）：
     (a) 電荷加性：E 對 ΔQ_had 或 Δ(ΣQ²) 線性 ⟹ 若 R² > 0.9 且殘差 < 10% ⟹ 未排除
     (b) substitution-count 線性：E ≈ k·(sub count) ⟹ 同上判準
     (c) spectator 獨立：同 substitution、唔同 spectator 嘅通道 E 應相等（|ΔE| < 2σ）
     (d) 標度假設：E ∝ M_ch^p（單一 p）⟹ 由 (n−p, D+−D0) 定 p，預測其餘 ⟹ dev < 10% 才未排除
  C3 依賴圖：若 (a)–(d) 全排除 ⟹ operator 必須含 mode-level（通道內 mode 佔用）資訊
     ⟹ 需 flavour→mode 映射 = Q5 selection rule = O140 ⟹ Q2/Q3 標 blocked-by-O140
輸出：data/18_em_constraints.json ＋ stdout
```

## 二、事前判準

| 項 | 事前寫死判準 |
|:--|:--|
| C1 | `E_ch = Δm_ch − coeff·X`，`X = 2.56 ± 0.094`（FLAG24，**外部固定，唔由本輪擬合**）；coeff = Σ(最終態 − 初始態) 嘅 u↔d 質量項 |
| C2a | 電荷加性：E 對 `Δ(ΣQ²)` 或 `ΔQ_had` 線性 ⟹ `R² > 0.9` 且殘差 < 10% 才「未排除」 |
| C2b | substitution-count 線性：`E ≈ k·|sub|` ⟹ 同上 |
| C2c | spectator 獨立：同方向、唔同 spectator ⟹ `|ΔE| < 2σ`（spread < 0.30 MeV） |
| C2d | 標度幂律：由兩點定 `p`，預測其餘 ⟹ worst dev < 10% 才「未排除」 |
| C3 | 若 (a)–(d) 全排除 ⟹ operator 必須含 mode-level 資訊 ⟹ 需 flavour→mode 映射（屬 O140）⟹ Q2/Q3 標 blocked-by-O140 |
| 全域 | 禁見結果改通道集；禁以「擬合 2 參數」當檢定；電荷值由 O124（`Y = T₃^R + (B−L)/2`）固定 |

## 三、執行記錄

| 項 | 結果 |
|:--|:--|
| Canonical run | exit 0，49 行；重跑 byte-identical（sha `018f044d0972b9e8`） |
| 通道集（8） | n−p、K+−K0、π+−π0、Σ−−Σ0、Σ0−Σ+、Σ−−Σ+、Ξ−−Ξ0、D+−D0（PDG 2024） |
| 外部固定 | `X = m_d − m_u = 2.56 ± 0.094`（FLAG24） |
| **符號公約（首版 bug 已修）** | `Δm = Σ(m_final − m_initial) + E`；`ΔK = m_K+ − m_K0`（u s̄ vs d s̄）⟹ X 係數 **−1**。首版誤用 +1 令 `E_K = −6.494`（假）⟹ 修正後 `E_K = −1.374`，與 `E_N = −1.267` 一致（H1 再證）。此係真 bug，已修並重跑 |

## 四、原始輸出（引自 log）

**C1 通道 em 讀數（MeV）**

| 通道 | Δm | coeff | E_ch | σ_E | 類型 |
|:--|--:|--:|--:|--:|:--|
| n−p | +1.2933 | +1 | **−1.2667** | 0.0940 | baryon u↔d，spectators {u,d} |
| K+−K0 | −3.9340 | −1 | **−1.3740** | 0.0963 | meson u→d，spectator s̄ |
| π+−π0 | +4.5936 | +1 | **+2.0336** | 0.0940 | annihilation |
| Σ−−Σ0 | +4.8070 | +1 | **+2.2470** | 0.1010 | baryon，spectators {d,s} |
| Σ0−Σ+ | +3.2720 | +1 | **+0.7120** | 0.1196 | baryon，spectators {u,s} |
| Σ−−Σ+ | +8.0790 | +2 | **+2.9590** | 0.2028 | baryon 2×u↔d |
| Ξ−−Ξ0 | +6.8500 | +1 | **+4.2900** | 0.2301 | baryon，spectators {s,s} |
| D+−D0 | +4.8200 | +1 | **+2.2600** | 0.1172 | meson，spectator c |

**C2a**：`Δ(ΣQ²)` slope +0.343，**R² = 0.0026**；`ΔQ_had` slope −0.874，**R² = 0.0240** ⟹ 🔴 電荷加性排除
**C2b**：`E = −0.416 + 1.687·|sub|`，**R² = 0.0894**，殘差 −2.646…+3.018 ⟹ 🔴 排除
**C2c**：A 組（n−p −1.267 vs K −1.374）spread **0.107** ⟹ 一致（H1）；B 組（Σ−Σ0 +2.247／Σ0−Σ+ +0.712／Ξ +4.290／D +2.260）spread **3.578** ⟹ 🔴 spectator 獨立排除
**C2d**：由 (n−p, D+−D0) 得 `p = 0.8422`；預言 K 0.740 vs 1.374（46.2%）、π 0.251 vs 2.034（87.7%）、Σ−Σ0 1.552 vs 2.247（30.9%）、Σ0−Σ+ 1.548 vs 0.712（117.4%）、Ξ 1.686 vs 4.290（60.7%）⟹ 🔴 幂律排除

## 五、結論報告

1. **四類「簡單」形式全部排除**：電荷加性（含 docs/17 已否證嘅 κ_sΣQ²+κ_xΣQ_iQ_j，本輪再做線性回歸確認 R²≈0）、substitution-count 線性、spectator 獨立、單一標度幂律 —— 冇一類可以描述 8 通道實證譜。
2. **唯一穩定結構：K/N 通道一致**（−1.374 vs −1.267，spread 0.107 ＝ 1.1σ）⟹ docs/17 嘅 H1（E_K = E_N）在此**獨立再證**（用完全唔同嘅通道集與外部固定 X）。
3. **方向不對稱（新事實）**：+X 方向嘅輕通道 E ≈ −1.3（K、N），而 d 取代 u 嘅重子通道 E = +0.71…+4.29（Σ0−Σ+、Σ−−Σ0、Ξ、D）⟹ **em 讀數對 spectator 內容極敏感**，唔係「替換對」嘅局部性質。
4. **operator 必須含 mode-level 資訊**：E 嘅跨度（−1.37 → +4.29，5.66 MeV）遠超任何電荷／計數／標度嘅單調函數 ⟹ 唯一可能形式 = 通道內**佔用 mode 集合**與邊界類型嘅函數。
5. **依賴圖（關鍵）**：mode 佔用 = quark flavour → K_chain mode 嘅映射 = **Q5 selection rule = O140**（已 handoff）。因此 Q2/Q3 **形式上下游於 O140**；本 O 依 10-16 裁決**唔可自行補此映射** ⟹ **Q2/Q3 標 blocked-by-O140**（由「等數據」改為「等他 O 產出」）。
6. **對 docs/17 嘅影響**：docs/17 嘅 H1 結論加強（非推翻）；R_em = 3.479 應重述為「π(annihilation) 通道 vs substitution 通道嘅 em 讀數比」，屬 operator 待定後才可解釋。

## 六、Q2/Q3 狀態更新

| 子項 | 之前 | 現在 |
|:--|:--|:--|
| operator 形式 | 「需 channel-dependent」 | **四類形式排除**：電荷加性／sub-count／spectator 獨立／標度幂律 |
| 通道實證譜 | 只有 π/K/N 三通道 | **8 通道**（含 Σ、Ξ、D）＋ 精確 E 值與誤差 |
| 方向不對稱 | 未知 | E(u→d 輕) ≈ −1.3；E(d→u 重子) = +0.71…+4.29 |
| 依賴 | 「需 em operator」 | **blocked-by-O140**（flavour→mode 映射 = selection rule）｜**handoff 已登記**：`2026-10-24T10:14:04` M4 → O140（register `handed_off`，ref = 本檔），非轉交產權 |
| 下一步 | — | 等 O140 的 flavour→mode 映射；屆時 M4 可用 8 通道約束集**直接檢定** operator |

## 七、與其他腳本／O 的關係

- ← `docs/14`（R_em、反相）、`docs/17`（最小電荷模型否證、H1）
- ← O124（電荷結構）、O140（selection rule／flavour→mode 映射）
- → `01_problem_record.md` §八（加 18 行）、§十一 Q2/Q3（改判 blocked-by-O140）、§十一.2（加四類形式禁區）
- 誠實標記：C2c 為「部分排除」（A 組一致、B 組排除）；K 通道 X 係數首版寫錯（+1）已修正重跑
