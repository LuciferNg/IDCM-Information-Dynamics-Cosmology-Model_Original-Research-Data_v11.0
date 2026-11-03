# M1 SESSION_START — 帳本編碼（ledger encoding）與判別力

<!-- registry
pid:     M1
title:   質量檔案 schema：帳本編碼（ledger encoding）與判別力（Reality-Architecture; 由 v10 O251 遷移）
state:   partial
verdict: partial
flags:   migrated-from-O251; route-ledger-encoding; format-layer-only; no-flavour-map; first-round-done; r6-nogo
owner:   M1
created: 2026-10-14
updated: 2026-10-28
prereq:  M2 + M4 + M5 + M6 + M8 + M9 + M10 + M11; O211/O212/O204/O237/O239/O140/O216/O144/O120/O209; O115
source:  調律者2026-10-14 O251→v11 M1；2026-10-28 整合 M2–M11 後重寫（原文存 backups/）
-->

> **Session**: 2026-10-28（第二次起手；路線重整，未開推）
> **狀態**: ⚪ registered — 路線已定、零 docs／零 scripts
> **原文**: `backups/{README,01_problem_record,SESSION_START}.md.20261028_1150.bak`（2026-10-14 版，不刪）

---

## 一、接手起點

M1 自 v10 O251 遷移後一直「只開不推」，原因已消失：**骨幹公式已由 M2 證成**，唔需要 M1 自己立質量式。

```
M2 §2（11/11）  : m = M₀·Σ σ_n γ_n        （σ ∈ [0,1]^44，零新常數）
M2 §1（16/16）  : m/M₀ + E_sync = Σγ = 90 （E_sync = 61.5，E_proj = 28.5，K_tree = 33）
M2 §4①（結構）  : Σ_{n≤k}γ_n = 28.5 無解  ⟹ σ 必含 ≥1 partial（f_proj = 0.771960195963）
```

⟹ M1 剩下嘅**全部係格式問題**，三問：

1. 一個質量檔係邊個 σ？（正規形／等價類）
2. 同一質量值有幾個合法 σ？（多重性 ρ）
3. σ 幾時有判別力？（判別力判準 + dof 成本 = partial 數 p）

**唔准做**：另立質量公式、預言質量值、由家族分組反推譜身份、推 M₀ 絕對值、碰束縛態／合法性／通道算子。詳見 `01_problem_record.md` §八 禁區 T1–T11。

---

## 二、先做嘅第一步（開推時，順序不可換）

```
1. R1 層約化表 → docs/01_layer_reduction.md
   五層逐格：vertices（γ_n 譜點）/ transform（M₀ 單位＋M2 泛函）/ operators（σ 選擇）/
            topology（等價類結構）/ labels（只可標註）
   每格標：可運算 / 只可標註 / 外部（M5、M10、M11）＋ 出處（O／M／docs／腳本）
   擾動測試：改該層 ⟹ 預言改變？不變者標「空層」（v11 L3 自由度閘）

2. R2 多重性引擎 → scripts/01_multiplicity_engine.py ＋ docs/02_multiplicity.md
   輸入：γ_n = 2 − 2T_n(cos π/44)（等價 4sin²(nπ/88)，max|Δ| 1.07e−29，mp 50 dps）
   計數：ρ(m; k, ε) = #{ supp(σ) ≤ k : |M₀Σσγ − m| ≤ ε·m }
         全二值支撐 ≤4 可全枚舉 C(44,≤4) = 149,986；k>4 用 meet-in-the-middle 拆分
   輸出：9 費米子（PDG 2024）× k=1..4 × ε∈{1e−3, 1e−4} 之 ρ 表
   註：partial 存在時 ρ 需按「≤p 個 partial」另計一欄（p=0 二值 / p=1 / p≤2）

3. MC null 10⁷ → scripts/02_mc_null.py
   m ~ U(0,90) 隨機目標之 ρ 分佈 ⟹ 判 e/μ 之 ρ 是否顯著低（Q1/Q2）
   同時報 M8 之功率判準（null 證據需 ρ_族 ≥ 3；命中需 ρ_族 ≤ 0.1）

4. 正規形／等價類 → scripts/03_canonical_form.py ＋ docs/03_canonical_form.md
   canonical 全序：n_q ↑ → partial 數 p ↑ → 支撐 k ↑ → lex(σ 索引序列) ↑
   證明存在／唯一；反例搜索 0 命中；若歧義 ⟹ format no-go（記入 docs/03）

5. 三鏈對照 → scripts/04_three_chain_presets.py ＋ docs/04_three_chain_presets.md
   O54 投影鏈／O239 |Δm|^{φ²}（M9 已判 428.4σ 死，只引用不重推）／M4 Δm²+R
   判「同一格式三個 preset」或「異格式」；φ² preset 封存，Δm²+R 判為代內差值算子
```

**收尾（每次必做）**：`scripts/05_doc_consistency.py`（文檔↔data 逐條核，ulp 由文檔精度推）＋ `scripts/run_chain.sh`（預設 dry-run，`--run` 才跑）＋ 全鏈審計（重跑 JSON byte-identical）。

---

## 三、已知陷阱（2026-10-28 更新）

```
⚠️ M₀ 兩義（T1）：錨值 27.9009145112 MeV（O204/O212/M4/M6/M8）vs O239 之 0.033331 MeV（1 參數 fit）。
   本 M 一律用錨值；引 O239 數字必寫「O239-M₀」。
⚠️ 命名碰撞（T2）：M4 docs/13/15 之「機制 M1／機制 M2」= 機制代號，唔係問題 M1/M2。引用寫「機制 M1（M4）」。
   改名請求只經 register handoff，**不在 M4 目錄改字**。
⚠️ Σγ = 90 係**每一量子**之上限，唔係質量上限：b/M₀ = 149.8159（n_q ≥ 2）、t/M₀ = 6185.1019（n_q ≥ 69）。
   任何「質量 > 90M₀」即出錯，或漏咗 n_q 欄。
⚠️ 完備性唔係成就：σ 連續（或二值 2^43 組合）落地極密 ⟹ 任何質量都「可達」⟹ 完備性**證明零預言力**。
   唔准把「找到編碼」當命中；必須報 ρ 與 dof（partial 數）。
⚠️ partial ≠ 免費：每個 partial = 1 dof。M2 已證最少需 1 個（f_proj）；若需 p ≥ 2 ⟹ 發還 M2（Q4）。
⚠️ 單模界：任何單模比值 ≤ CAP = γ₄₄/γ₁ = 784.9646644；τ/e = 3477.23 = 4.43×CAP ⟹ 跨代必多模（M4/M9）。
⚠️ 密度／判別力閘：ρ_族 ≤ 0.1 才算命中；ρ_族 ≥ 3 才可以「落空」當證據；單目標證據上限 ≈ 1 bit（M8 R7′）。
⚠️ 免死金牌（O239 R36 遺留）：schema 未完唔可以當 O239 60% 殘差之解釋；殘差歸 O239/O136。
⚠️ γ_n 數值：n 大時 4sin²(nπ/88) 接近 4；用 2−2T_n(cos π/44) 或 mpmath 50 dps，禁 float64 判等。
⚠️ 版本紀律：v11 內消化，**禁 v11→v10 反向推送**（register 硬拒）；跨版本引用必先核歸一化。
⚠️ register 已遷：`IDCM-Original/problem_register/`（跨版本共用），pid 帶版本 = `v11.0/M1`；舊路徑 v10 `problem_register/` 已失效。
⚠️ 檔案衛生：改任何檔先 backup（`backups/<name>.<TS>.bak`）；只寫本 M 目錄；跨 M 只經 register handoff／note。
```

---

## 四、關鍵數字（本 M 直接使用）

| 量 | 值 | 出處 |
|:---|:---|:-----|
| Σγ | **90**（sympy 精確 Σsin² = 45/2） | O211／M2 §1 |
| E_sync / E_proj / K_tree | 61.5 / 28.5 / 33 | M2 §1 |
| f_proj | 0.771960195963 | M2 §1 |
| ‖P²−P‖ | 0.176037651812（收縮非冪等） | M2 §2 |
| M₀（錨） | 27.9009145112 MeV | O204／O212 |
| γ_n | 2 − 2T_n(cos π/44) ≡ 4sin²(nπ/88) | M4 docs/10；M8 R8′ |
| γ₁ / γ₄₄ | 0.00509577077949 / 4 | — |
| γ_n + γ_{44−n} | 4 | M8 R8′ |
| CAP = γ₄₄/γ₁ | 784.9646644 | M9／M4 docs/19 |
| λ_max(k) | 4k − Σ_{j<k}γ_j：4 / 7.994904 / 11.974547 / 15.928841 / 31.296348（k=1,2,3,4,8） | M6 docs/07 |
| 身份清晰界 | λ < 8（λ_max(2) = 7.994904） | M6 docs/07 |
| m/M₀（PDG 2024） | e 0.0183148・μ 3.7869144・τ 63.6846509・p 33.6287217・c 45.6257446・b 149.8158778・t 6185.101923 | 本檔 §二 計算（mp 40 dps） |
| k_min（λ ≤ 4k 界） | e 1・μ 1・τ 19・p 9・c 13；b/t 超單量子 | 同上 |
| n_q ≥ ⌈λ/90⌉ | p 1・b 2・t 69・v_SM 99 | M6 docs/07 §3；本檔核 |
| τ/e | 3477.23 = 4.43×CAP | M9 |
| 質量層 dof | 0（9! 精確解）⟹ 零預言力 | O120／M10 W5 |
| 判別力門檻 | 命中 ρ_族 ≤ 0.1；落空證據 ρ_族 ≥ 3 | M8 R7′／M9 |

---

## 五、交接 / 依賴

```
- 取（只讀）: M2（泛函/帳本）、M4（譜代數/單模界/密度閘）、M6（λ_max/量子數）、
              M8（M₀=單位/功率判準）、M9（CAP/四判準）、M10（W5 牆）、M5（束縛態存在）、M11（通道算子缺口）
- 回（僅經 register）:
    M2 ← 若 Q4 判 p ≥ 2（泛函唯一性受挑戰）
    M6 ← 若 ρ 計數與 λ_max 界不一致（先對帳再續推）
    M4 ← 機制代號改名建議（M1/M2 → K1/K2），不代改
    M5/M10/M11 ← 介面欄位需求（`bound_state` / `legality` / `channel`），只提契約不提推導
- 範圍警告:
    唔推夸克時空存在（O247 Q2 已裁）；唔重推 O248（相位橋 no-go）；
    唔重推 O239 殘差分類；唔碰 M5/M10/M11 之攻堅；唔預言任何質量值
```

---

## 六、首輪腳本清單（命名先定，避免日後改名）

```
scripts/01_multiplicity_engine.py   ρ(m;k,ε) 精確計數（mp 50 dps；k≤4 全枚舉 + MITM）
scripts/02_mc_null.py               10⁷ 隨機目標 null；功率判準對照
scripts/03_canonical_form.py        canonical 全序＋唯一性＋反例搜索
scripts/04_three_chain_presets.py   O54／O239(φ², 封存)／M4 Δm²+R 三鏈對照
scripts/05_doc_consistency.py       文檔↔data 逐條核（ulp 由文檔精度推，零可調容差）
scripts/run_chain.sh                driver（預設 dry-run；--run 逐支執行並寫 log）
```

環境：`/usr/bin/python3`、`PYTHONDONTWRITEBYTECODE=1`；精確算術用 `Fraction`／`mpmath dps≥50`；浮點只准最終顯示。

---

## 七、首輪進度（2026-10-28 完成）

```
R1 層約化  ✅ docs/01（擾動：鏡射恆等 4.2169e−81；N=86/88/90 令 9/9 改變；L4 置換 200/200 不變）
R2 多重性  ✅ docs/02（ρ_bin：7/9 為零；μ k*=2 deg=1；s k*=3 deg=7；coverage k≤4 = 17.040 %；MC 10⁷ p(μ)=0.113）
R3 正規形  ✅ docs/03（τ k=19 deg=513；c k=13 deg=4285；退化解 σ₄₄=t/4；b n_q=2、t n_q=69）
R4 三鏈    ✅ docs/04（一骨架三 preset；φ² 域外：√5 conductor 5 ∤ 88）
對帳       ✅ docs/05（Q1 🟡、Q2/Q3/Q4 🔴、Q5 ✅）
鏈         ✅ 6 支全 exit 0、JSON 重跑 byte-identical、文檔↔數據 89/89（~98 s）
```

**接手起點（下一輪，唯一開放）**：**R6 跨量子組合律**

```
1. 目標：b（n_q=2，餘額 59.8158778 = 149.8158778 − 90）
             t（n_q=69，餘額 24.8980769 = 69×90 − 6185.10192312）
2. 預註冊（docs/05 §5，先寫死）：
   R6-Q1 每量子 ≥ 86 = 90 − γ₄₄（否則與 t 之均值 89.63915831 衝突）
   R6-Q2 兩餘額必須同源（同一算子作用 n_q 次）
3. 判死線：候選律需 ≥2 個新參數 ⟹ 直接否證（dof 預算）
4. 禁止：由質量值反解量子分配（即 ρ 爆炸陷阱）；必須先立分配律再檢定
5. 若 R6 判死 ⟹ 本 M 收於 partial/negative：格式層「跨量子不可唯一化」本身即結論
```

**樣本外保留**：本輪質量用 PDG 2024 單快照；下一輪須以 PDG 2022 重跑 `docs/02 §1` 之 ρ 表作樣本外核。

---

## 八、R6 完成（2026-10-28 同日）—— no-go

```
結果：跨量子組合律不可識別
  E 等分      b: 2×74.90793892   t: 69×89.63915831     ✓ 精確
  F 滿載+餘額  b: 90+59.815877839 t: 68×90+65.10192312  ✓ 精確
  ⟹ 兩律等價（σ 連續 + n_q 不可觀測）⟹ 質量值層不可判定
六個量子含量（h_b/r_b/h_t/r_t/c_b/c_t）：k_min 7–39、ρ(k_min) 19–1393、k≤2 全零命中 ⟹ 全部無譜身份
預註冊：R6-Q1 🔴（≥86 係由 t 均值外推）、R6-Q2 🔴（常數虧損律假設）
新禁區：D1 禁由質量反解量子分配｜D2 禁以 γ₁ 格點命中當結構
文檔：docs/06_multiquantum.md；script 06（8/8）；一致性 149/149
```

**現時唯一未閉事項**：無（內部攻堅窮盡）。重開觸發器見 README「未閉（觸發器）」三條。

---

## 九、下游重推（2026-10-29，docs/08）

```
掃 11 個 M 之依賴（n_q/⌈λ/90⌉/17.040/20.61/CAP/Σγ=90/Σσ_nγ_n/密度閘/coverage）
⟹ 7 受影響（已 push：M2 M4 M5 M6 M10 M11 M12）；4 不變（M3 M7 M8 M9，附證據）
★ P3 慣例：k≤4 coverage = uniform 17.040% / log-uniform 35.155% / M4 慣例 24.877% ⟹ 引用必標慣例
★ P4 支撐閘：coverage(k,1e−3)≤0.1 ⟹ uniform 只到 k=2；log-uniform 只 k=1
★ P1 M₀ 兩義解除（0.033331 = v10 O239 1 參數 fit；非矛盾）
scripts/07 13/13 PASS；全鏈 8 支 exit0；docs↔data 171/171
```

---

> ⚠️ **本檔係起手式快照，唔係狀態權威。**
> 當前狀態以 `README.md` 為準；生命週期事件以 `IDCM-Original/problem_register/`（pid `v11.0/M1`）為準。