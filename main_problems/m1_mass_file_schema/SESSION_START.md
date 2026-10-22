# M1 SESSION_START — 質量檔案 schema：譜點/標度/算子/拓撲/標籤五層格式規範 (Reality-Architecture; 由 v10 O251 遷移)

<!-- registry
pid:     M1
title:   質量檔案 schema：譜點/標度/算子/拓撲/標籤五層格式規範 (Reality-Architecture; 由 v10 O251 遷移)
state:   registered
verdict: -
flags:   migrated-from-O251
owner:   M1
created: 2026-10-14
updated: 2026-10-14
prereq:  O237 (mass pipeline) + O248 (讀數裁決) + O239 + O140 + O212 + O216 + O144; O115 (v11 遷移協議)
source:  調律者2026-10-14: O251 遷移至 v11 Reality-Architecture main_problems/m1 (承接 2026-10-13 讀數reframe)
-->

> **Session**: 2026-10-14（遷移至 v11 Reality-Architecture）
> **狀態**: 📄 **已註冊 — 未開始（推導用新 session）**
> **來源**: v10 O251 遷移——調律者 2026-10-13「開新o, 但不直接推」+ 2026-10-14「O251 由 v10 移動到 v11 main_problems, 編號 M1」
> **前置**: O237（mass pipeline roadmap）；O248（讀數裁決）；O239（閉合式/殘差）；O140（階梯正典）；O212（dof=1）；O115（遷移協議）
> **位置**: `main_problems/m1_mass_file_schema/`

---

## 一、接手起點

未開始——v10 O251 已完整遷移（01_problem_record + SESSION_START + 起手目錄）；推導留新 session。
原 M1 前身 O251 喺 v10 只開未推（2026-10-14 註冊, register state 已轉 merged → v11 M1）。
承接: 調律者 2026-10-13「質量檔案 schema」reframe——唔推「夸克在時空如何存在」，推「質量檔案嘅函數結構同變量」。
O248 已裁（11/11 PASS）：相位模→質量污染無載體，質量 = Object−Frame 脫相位讀數（O237）——本 O 係呢個裁決嘅建設性後續。

## 二、先做嘅第一步（新 session 開推時）

```
1. 讀 O237 docs/01_mass_pipeline_roadmap.md —— 現有 mass pipeline（S0-S4 譜層到公式層）,
   將佢重劃為 schema 五層: vertices（譜點）→ transform（M₀）→ operators（γ_n/|Δm|/F_proj/S_f）
   → topology（9 點 family）→ labels（色/flavor）——出「五層表」底稿（每層標 已閉合/部分/缺口 + 來源）
2. 逐層 audit:
   vertices: ε_n = 4sin²(nπ/88)（O204/O140）——補 O237 roadmap 嘅 S0 表
   transform: dof=1, M₀=27.9009145112（O212/O204/O38）
   operators: m_f = M₀γ_nF_proj（O54/O64, C_f≡1 under O137）vs m = M₀|Δm|^{φ²}（O239, α=φ²）——兩鏈對照
   topology: O239 docs/35 c(n) 表 → 9 點 family U/D/L 分組有冇譜身份（過 O206 池密度閘先可宣稱）
   labels: O216（色=嵌入標籤）+ O144（禁閉定理）——呢層大概率已閉合, 只需收錄
3. topology 攻堅（成本 medium）: c(n) vs 譜池 {γ_n, |Δm|, F_proj, C_f, S_f} 掃描
   —— 條件: 池密度閘 P=0.4175 先例; 任何「命中」要 MC null + 樣本外
4. modifier 鏈統一: O54 投影鏈 vs O239 |Δm|^{φ²} 鏈——逐項對照輸入/輸出, 尋同源算子身份
   —— no-go 都係答案: 兩條 chain = 兩種 export serialization（寫明）
```

## 三、已知陷阱

```
⚠️ 免死金牌警示（O239 R36）: schema 未完唔可以當 O239 60% 殘差嘅解釋——殘差歸屬仍屬 O239/O136
⚠️ 池子密度（O206 P=0.4175）: 任何「family/topology 命中」要過池閘先算數; φ-power trap 照舊
⚠️ 版本隔離: 跨版本質量公式不相通（v5/v1/v3 先例係方向證據唔係 v10 公式）——M1 屬 v11 Reality-Architecture,
   引用 v10 數字必先核歸一化; v11 VERIFICATION 密度基準落地後本 O 屬 P 類須重驗（O115 docs/03）
⚠️ frame 配對: 質量公式用 O206 CSSG-v3 assignment（e(1,3.6) μ(7,15.5) τ(20,37) u(1,15.5) d(1,33) s(14,3.6) c(34,13) b(33,44) t(44,1554)）
⚠️ freeze v2 質量集規則（O234 docs/03）: m_b = 4990.9（v2-3GeV）, 唔好再用 ×3.291 = 3581.6（混標度偽影）
⚠️ 唔重複: O247 Q2（核子/原子絕對質量 no-go 已裁）; O248（相位橋已閉合）; O239 已做嘅殘差 audit
⚠️ v11 現狀: 未有 CHARTER/VERIFICATION（O115 憲章未落地）——本 O 先以 v10 結構開推, 渡航時再重驗
```

## 四、關鍵數字

| 量 | 值 | 出處 |
|:---|:---|:-----|
| 唯一有量綱錨 | dof=1（G≡ω_ring≡M_Pl） | O212 |
| M₀ | 27.9009145112 MeV | O204/O38 |
| 譜點 ε_n | 4sin²(nπ/88), n=1..44 | O204 T1/O140 |
| n/p BBN | ε₇/ε₂₀ = 0.14261761 | O140 docs/16 |
| 閉合式殘差 | 0.5876 → RGI 0.6567 (1.117×) | O239 R35 |
| 污染因子 c(n) | b0.574/μ0.623/d0.624/e0.648/t0.978/τ0.957/u0.935/c2.508/s3.146 | O239 docs/35 R1 |
| 色荷 | 嵌入層標籤（非譜量） | O216 docs/03 |
| 禁閉 | Z₄ 子格存在、Z₃ 不存在（divisors(176) 無 3 倍數） | O144 docs/24 |

## 五、交接 / 依賴

```
- 前置: O237（mass pipeline 現狀）; O248（讀數裁決）; O239（閉合式/殘差/topology 輸入）;
         O140（階梯正典）; O212（dof=1）; O216/O144（標籤層）; O115（v11 遷移協議）
- 移交出: topology 規則若成立 → O206 池閘 / O136 逐點公式; 殘差 → O239（照舊）;
           export 翻譯層 → O139/O122（本 O 只 reference 唔重推）
- 範圍警告: 唔推夸克時空存在（O247 Q2 已裁: 複合層唔喺 K_chain 譜）;
           唔重複 O248（相位橋已閉合 verdict=negative）;
           唔重複 O239 已做（scheme 污染 audit / 殘差分類）
```

---

> ⚠️ **本檔係起手式快照，唔係狀態權威。**
> 當前狀態以 `README.md` 為準；生命週期事件以 v10 `problem_register/`（M1 遷移記錄）為準。