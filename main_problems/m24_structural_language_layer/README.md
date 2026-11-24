# M24 — 結構語言層（structural language layer）

## 狀態

🟡 **partial — 兩輪推完**（2026-11-25）｜L0 ✅ 建成 · L1/L2 🔴 否證 · L3 🟡 收窄（已量化）
重現閘 17/17｜scripts 99 PASS / 0 FAIL（00–08）｜9 docs｜9 data
第二輪：R-a／R-b 復活失敗（統一索引族 0 命中）；R-c 得 ℚ(√5) ⊄ ℚ(ζ₈₈) 障礙；T3 收緊至定義級（內容 ≤3 bits）

## 一句話

M21 已把四殘項嘅缺口講到盡：**唔缺數字、唔缺算子，缺嘅係三層語法 ＋ 一層元語言**。
M24 就係要建（或舉證不可能）呢三層：`L2` 分次標籤代數 · `L1` 上同調-Hodge 函子 · `L3` 測度傳輸層 · `L0` count-identity 判準。

## 三層 ＋ 一層元語言

```
L2  分次標籤代數：每模式帶 B-值標籤 ⟹ 計數相乘
    目標：由同一標籤系統逼出 N₂ = B·B = 16、h11 = 36、h21 = 98
L1  上同調／Hodge 函子：dim H^{1,1} = 36、dim H^{2,1} = 98 由晶格／表示資料**計出**
    現況：h21 = Σ_(i<j) d_i d_j = 98 有維數讀法但係代數恆等；h11 = 36 池內零 canonical 命中
L3  測度傳輸層：非等變之 pushforward（群普查 measure ↔ 模譜 density）
    ⚠ 必先對接 v10 O210 docs/11「**測度層 🔴 永久 no-go**」—— 唔接即重推已否證路線
L0  元語言：結構身份由「交換圖存在」定義，而唔係數字相等
    現況：所有「M 型」判詞靠密度閘（5%／k-free ≤5 條）＝啟發式
```

## 殺手測試（唯一硬判準）→ **FAIL**（第一輪）

```
predict:  (N₂, h11, h21) = (16, 36, 98) 同源 ＋ 密度 ≤5% ＋ 零自由參數
判定：K1 同源 FAIL（三數跨兩系統：repr {16,98}｜chain {36}；交集 = {16}）
      K4 零參 FAIL（36 於模側須逐數選 isotypic pattern {3,3,4}）
      K6 同級 FAIL（L0 grade：16 = F｜98 = M｜36 = M）
⟹ 非同一標籤系統之產物；F 型升格不成立
```

## 範圍收窄

```
✓ 建／否證 L0–L3 四層語法
✗ 唔重推 M21 AB1–AB4 判詞（所有權在 M21，只作測試用例）
✗ 唔重做 v10 O210 測度層 no-go（要挑戰須先指其假設缺口）
✗ 唔碰 h11／h21 數值來源（v10 O211 零輸入）
```

## 檔案

| 檔 | 內容 |
|:--|:--|
| `01_problem_record.md` | §〇 fence／Q1–Q4／T1–T4／P1–P5／R1–R5／禁區／判準／殺手測試／**§十一 第一輪結果** |
| `SESSION_START.md` | 接手起點、重現閘六項、範圍收窄、禁區、分工／**§九 第一輪 session 記錄** |
| `README.md` | 本檔 |
| `docs/00_preregistration.md` | 跑前凍結：𝕃 3 系統／𝔉 5 函子／H1–H4／M1–M3／S1–S4／校準集 |
| `docs/01_l2_label_algebra_result.md` | T1 🔴：三系統庫存、交集定理、MC null 診斷 |
| `docs/02_l1_hodge_functor_result.md` | T2 🔴：F1 (4,0)/(4,4)、F2 59、Galois 障礙定理 |
| `docs/03_l3_measure_transport_result.md` | T3 🟡：O210 範圍閘、零參數傳輸三路、條件式定理 |
| `docs/04_l0_meta_language_result.md` | T4 ✅：S1–S4 判準、10/10 校準、3 盲測、grade 分裂 |
| `docs/05_interdependence_and_ledger.md` | P5 🔴：分裂定理、總帳、復活條件 R-a–R-d |
| `docs/06_preregistration_round2.md` | 第二輪跑前凍結：R-a／R-b 規則族、C1–C3、D1–D3 |
| `docs/07_r_a_r_b_uniform_rule_result.md` | R-a／R-b 🔴：117 格 × 13 濾子 0 命中；後驗 36 讀法（M 型） |
| `docs/08_arithmetic_field_obstruction.md` | R-c 🟡：ℚ(√5) ⊄ ℚ(ζ₈₈)（定理＋機核）⟹ 像 ⊆ ℚ |
| `docs/09_t3_dof_result.md` | T3 收緊：零參數⟹Aut-不變（定義級）；內容 ≤3 bits = {20,24} |
| `scripts/00–08` `data/*.json` `logs/*.txt` | 九輪機核（99 PASS / 0 FAIL，全 exit 0） |

## 依賴

```
上游：M21（四殘項測試用例）· M22（BR1′–BR1‴）· v10 O210（測度層 no-go）· v10 O211（h11/h21）
      M12（Hodge↔譜 無先驗對應）· M16／M6（dof）
下游：M21 AB1–AB4 之最終判 · IDCM 新描述語言（框架長期方向）
```

## 註冊

`register.py` pid = `M24`｜state = `partial`｜verdict = `partial`（2026-11-25）｜
source = 調律者 2026-11-20「欠缺什麼機制和什麼語言」→「開新M，待下一session建構」
