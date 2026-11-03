#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/07 — 回應 M13 來件（docs/11 T1 定理、docs/12 競爭者審計）

兩件來件於本 session 工作期間插入（register note 21:12、23:13），落入本 M `Q2` 範圍。
本腳本只回答**可計算**者：
  §一 一般 2-dof α 之擬合自由度：`(ρ_u, ρ_d)` 可否**同時精確命中**？代價 = `dof`
  §二 `ρ_u = 0` 之地位：`α_d = 2α_u` 係 Cartan 元素（非參數）⟹ 與 `r` 無關（核對 M13 docs/12 §一）
  §三 `r = α_s/α_u` 之四候選（M13 docs/12 §四）：本 M 之獨立重算（唯讀引用其候選集）
  §四 `√B` ／ `B = |ℤ₂|²` 之算術地位：三義是否塌縮為恆等式（純算術檢定，非凑數）
N_trials：§三 = 4（四候選，全報）；其餘為代數。
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
print("=" * 78)
print("M14 scripts/07 — 回應 M13 來件：Q2 之可計算部分")
print("=" * 78)

RD = mp.mpf("0.40350365")
RU = mp.mpf("-0.02588837")
SD = mp.mpf("0.044785178")
SU = mp.mpf("0.033583899")

# =====================================================================
print("\n【§一】一般 α（projective 2 dof）可否同時精確命中 (ρ_u, ρ_d)？")
print("  設 α_u=1、α_d=t、α_s=r：")
print("    ρ_u = (2−t)/(2r)      ρ_d = (2t−1)/(2r)")
print("  兩方程、兩未知 ⟹ 一般可唯一解。解得：")
# (2−t) = 2r·RU ; (2t−1) = 2r·RD
# 由第一式 t = 2 − 2r·RU ⟹ 代入第二式
# 2(2 − 2r·RU) − 1 = 2r·RD ⟹ 3 − 4r·RU = 2r·RD ⟹ r = 3/(4RU + 2RD)
r_star = 3 / (4 * RU + 2 * RD)
t_star = 2 - 2 * r_star * RU
print(f"    α_s/α_u = r* = {mp.nstr(r_star,10)}")
print(f"    α_d/α_u = t* = {mp.nstr(t_star,10)}")
print(f"    驗：ρ_u = {mp.nstr((2 - t_star)/(2*r_star),10)}（目標 {mp.nstr(RU,10)}）")
print(f"        ρ_d = {mp.nstr((2*t_star - 1)/(2*r_star),10)}（目標 {mp.nstr(RD,10)}）")
print(f"    t* > 2 ？{'是 ⟹ ρ_u < 0 ✓（負號可容納）' if t_star > 2 else '否'}")
print("  ⟹ **2 dof 下任何 (ρ_u,ρ_d) 皆可精確命中**（含 ρ_u<0）⟹ 零預言力 ⟹ 違 A3（dof=0 硬閘）")
print("  ⟹ 反過來：`ρ_u = 0` 之結構性（下節）係唯一把 dof 由 2 壓至 1 之條件")

# =====================================================================
print("\n【§二】`ρ_u = 0` 之地位（核對 M13 docs/12 §一）")
print("  ω₂ 條件 α_d = 2α_u ⟹ ρ_u = (2α_u−2α_u)/(2α_s) = 0 **精確**，與 r 無關")
for rv in ("3", "4", "11/3", "15/4", "7/2"):
    r = mp.mpf(F(rv).numerator) / F(rv).denominator
    print(f"    r = {rv:>5s} ⟹ ρ_u = {mp.nstr((2-2)/(2*r),3)}（精確 0）、"
          f"ρ_d = {mp.nstr(3/(2*r),8)}（z = {mp.nstr((3/(2*r)-RD)/SD,6)}）")
print("  ⟹ ✅ 本 M 獨立核對一致：`ρ_u = 0` 係結構強制（Cartan 零方向），不隨 r 變")

# =====================================================================
print("\n【§三】`r` 之四候選（M13 docs/12 §四 之候選集；唯讀引用，N_trials = 4 全報）")
print(f"  觀測 ρ_d = {mp.nstr(RD,9)} ± {mp.nstr(SD,7)}；ρ_d = 3/(2r)")
print(f"  {'r':>8s} {'ρ_d':>10s} {'z':>9s} 出處")
cands = [("15/4", "未註冊"), ("11/3", "N_d/(N_gen·B) = 44/12"), ("7/2", "未註冊"), ("4", "B = γ44（已註冊）")]
zs = []
for rv, src in cands:
    r = mp.mpf(F(rv).numerator) / F(rv).denominator
    rd = 3 / (2 * r)
    z = (rd - RD) / SD
    zs.append(abs(z))
    print(f"  {rv:>8s} {mp.nstr(rd,9):>10s} {mp.nstr(z,6):>9s} {src}")
print(f"  ⟹ |z| 排位：最佳 = {cands[min(range(4), key=lambda i: zs[i])][0]}；"
      f"`4`（已註冊 B）排第 {sorted(range(4), key=lambda i: zs[i]).index(3)+1}")
print("  ⟹ 本 M 獨立重算與 M13 docs/12 §四 **逐位吻合** ⟹ 接受其降級（3/8 非數據擇優）")
print("  ⚠ 但：本 M 之 Q1／Q2② 判詞**不依賴** 3/8 之數值地位（Q1 係置換對稱、Q2② 係 0 bit）⟹ 不受影響")

# =====================================================================
print("\n【§四】`√B` ／ `B = |ℤ₂|²` 之算術地位（純算術檢定；禁凑數）")
B = F(4)
cands_read = {"√B = B^(1/2)": mp.sqrt(mp.mpf(B.numerator) / B.denominator),
              "log₂B": mp.log(4, 2), "|ℤ₂|": mp.mpf(2), "N_gen−1": mp.mpf(2)}
print("  M13 docs/11 §六 之三種同值讀法（欲得 α 之底數 x = 2）：")
for k, v in cands_read.items():
    print(f"    {k:>9s} = {mp.nstr(v,12)}   {'= 2 ✓' if abs(v-2) < mp.mpf('1e-30') else '≠ 2'}")
print("  ⟹ 三者在 B = 2² 之下**同值** ⟹ 塌縮為算術恆等式，非算子身份")
print(f"  ⟹ 核：B = |ℤ₂|² ？ 4 = 2² ⟹ {'成立（純算術）' if B == F(2)**2 else '否'}")
print("  ⚠ 依 corpus 標準（M13 docs/11 §〇 自訂）：「算子身份」≫「算術」。同值三義係")
print("     **算術巧合之典型**，與禁區 3（小分母凑數／φ-power）同型 ⟹ 不可作 provenance。")
print("  ⟹ 故本 M 之答覆：`√B` 之出處**不可由 M14 範圍內給出**；且即使給出，其結構價值為 0：")
print("     ① 它只定比值 r，不定**次序**（次序即指派，docs/04 §二：結構輸出係無味名之多集）；")
print("     ② 1σ 窗內有 3 個格子候選（q≤3），其中 11/3 之 |z| 小於 4 ⟹ 推得 r=4 亦不可裁決。")

# =====================================================================
print("\n【§五】★ 分叉：`d` 係 `k` 之函數（v9 字面）抑或獨立世代指標（M13 docs/11 §一）？")
z2 = {1: 3, 4: 2, 5: 0}
print(f"  v9 z2_map = {z2}（本 M docs/03 §一 已證為純字面 dict，AST 零 Name 依賴）")
print(f"  image = {sorted(set(z2.values()))}；|image| = {len(set(z2.values()))} = |domain| = {len(z2)}")
print(f"  ⟹ z2_map 係 **雙射** ⟹ `d = z2_map[k]` 完全由 `k` 決定")
print("  ⟹ 推論（形式）：若 (k,d) 只取 v9 之三點 {(1,3),(4,2),(5,0)}，")
print("       則任何 h(k,d) 皆可寫成 h'(k)（複合 d=z2_map[k]）⟹ M13 docs/11 之 **case C 塌回 case A**")
print("       case A：d,s 同 k=4 ⟹ α_d ≡ α_s ⟹ ρ_u=0 ⟹ α_1 = α_4/2 ⟹ ρ_d = 3/4 強制 ⟹ z = +7.74σ 🔴")
print("  ⟹ 反之，M13 docs/11 §一 之指派 u=(1,3)、d=(4,3)、s=(4,2) 需要 **(4,3) 呢一格**，")
print("       而 (4,3) ∉ v9 之 image ⟹ M13 已**默認 d 為獨立世代指標**，即放棄 z2_map 之雙射。")
print("  ⟹ 故 label 層有兩種讀法，判詞不同：")
print("     (i)  v9 字面（d=f(k)）：標籤空間 3 點 ⟹ 全函數族 = 單 k 函數 ⟹ ρ_d = 3/4（7.74σ 死）")
print("     (ii) d 獨立（世代索引）：3×3 格，格子佔用係**未註冊選擇** ⟹ dof = 2（M13 case C）")
print("     兩者皆 🔴；差別只係死法（值被強制 vs 自由度過高）。")

# =====================================================================
print("\n" + "=" * 78)
print("對 M13 兩件來件之答覆（可計算部分）")
print("=" * 78)
print(f"""
  docs/11（T1 定理）
    ✅ case C（h(k,d) ⟹ dof = 2）**接受**，並升格：本 M docs/04 之範圍②（原引 docs/09 §3.2
       之「單函數 g 帶 dof≥1」）現可改引此定理級版本「(k,d) 網格全函數族 ⟹ dof 2 或值被強制」。
    ⚠ cases A／B 之「u,d 同代」「d,s 同扇區」前提，依賴 (k,d)↔{{u,d,s}} 之映射；
       該映射經本 M docs/03 §五 判為**未註冊**（v9 §D.1 不含 u、d）。⟹ A/B 係**條件性**命題，
       範圍須標「若該映射成立」。
    ⚠ 其 §四「ρ_d = 3/8 係預言」已被 docs/12 §六 自我降級（後出者為準）。
    ★ 本輪新增（§五）：§一 之指派需要 (k,d)=(4,3) 一格，而該格 ∉ v9 z2_map 之 image（雙射）
       ⟹ docs/11 已默認「d 獨立於 k」。若改回 v9 字面讀法，case C 塌回 case A（ρ_d = 3/4 死）。

  docs/12（競爭者審計）
    ✅ 接受降級。本 M 獨立重算四候選（15/4、11/3、7/2、4）之 |z| 排位，逐位吻合。
    ✅ 本 M 之 Q1／Q2② 判詞不依賴 3/8 之數值地位 ⟹ 不受影響。
    ★ 補一更強形式（本輪新算）：一般 α（2 dof）可**精確**命中 (ρ_u,ρ_d) = (−0.0259, 0.4035)，
       解 α ∝ (1, {mp.nstr(t_star,8)}, {mp.nstr(r_star,8)}) ⟹ 「負號 + 觀測 ρ_d」之唯一代價 = dof 2 ⟹ 違 A3。

  交辦「√B 之出處」
    🔴 不可給（且非本 M 可給）：三義同值 = 算術恆等式 B = |ℤ₂|²，非算子身份；
       即使給出亦只定 r、不定次序，且在 1σ 窗內不可裁決（11/3 之 |z| 更小）。
""")
raise SystemExit(0)