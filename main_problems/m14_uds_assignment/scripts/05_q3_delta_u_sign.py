#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/05 — Q3：`δ_u < 0` 可否由結構產生？

問題（record §二 Q3）：觀測 δ_u = −0.0887 ± σ；M13 之 ρ_u = 0 結構解給 δ_u = 0 精確。
問：負號可否由結構產生？抑或本語言（A₂ 方向 ω₂）結構上只能給 0，負號屬語言外殘餘？

三段檢定：
  §一  𝒜 層（一體＋兩體）：δ_x ≡ 0 是否精確（承 scripts/02 L3）
  §二  α-模型層：ρ_d = 3/8 與 ρ_u = 0 是否**同一條件的兩面**（族內代數）
  §三  觀測層：負號是否已達可裁決標準（3σ）；以 χ²(x) profile 定裁決力
判準：A3（dof=0）／A5（範圍）／A6（味盲）。禁加負號自由度（= dof>0）。
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
print("=" * 78)
print("M14 scripts/05 — Q3：δ_u < 0 之結構來源")
print("=" * 78)

RD = mp.mpf("0.40350365")
RU = mp.mpf("-0.02588837")
SD = mp.mpf("0.044785178")     # M12 docs/42 R27 正典
SU = mp.mpf("0.033583899")

# 觀測量（scripts/01 R1/R2 已重現；此處只列 σ）
DU = mp.mpf("0.088667660")     # |δ_u| MeV
S_DU = mp.mpf("0.11244")       # σ(δ_u) = √(σ(ΔΣ)² + (σ(ΔΞ)/2)²)，scripts/01 §R2
print(f"\n【§〇】觀測：δ_u = −{mp.nstr(DU,6)} ± {mp.nstr(S_DU,5)} MeV  ⟹ |δ_u|/σ = {mp.nstr(DU/S_DU,4)}")
print(f"      ρ_u = {mp.nstr(RU,8)} ± {mp.nstr(SU,6)}  ⟹ z(ρ_u vs 0) = "
      f"{mp.nstr(-RU/SU,4)}")
print(f"      3σ 判別負號需 |δ_u| ≥ 3σ = {mp.nstr(3*S_DU,4)} MeV（現值之 "
      f"{mp.nstr(3*S_DU/DU,3)} 倍）")

# =====================================================================
print("\n【§一】𝒜 層：δ_x 是否精確為 0")
print("  scripts/02 §四 已證：𝒜 之任一算子只依『佔用模之多集』與『模距之多集』")
print("  ⟹ 任何味對差 Δ(f,f') 之期望值 ≡ 0（精確，10⁷ null bit-identical）")
print("  ⟹ 於 𝒜 內，ρ_u 之結構值 = 0 對**任何** α 皆成立（不止 ω₂ 方向）")
print("  ⟹ 故『δ_u < 0』在 𝒜 內 **不可生成**（0 無符號可言）")

# =====================================================================
print("\n【§二】α-模型層：ρ_d = 3/8 與 ρ_u = 0 是否同一條件之兩面")
print("  族 α = (1, x, x²)（docs/10 之家族；x=2 唯一雙通）：")
print("    ρ_u(x) = (2−x)/(2x²)        ρ_d(x) = (2x−1)/(2x²)")
print("  解 ρ_d(x) = 3/8（精確 Fraction）：")
print("    (2x−1)/(2x²) = 3/8 ⟹ 16x − 8 = 6x² ⟹ 3x² − 8x + 4 = 0")
disc = mp.sqrt(64 - 4 * 3 * 4)
roots = [F(8 + 4, 6), F(8 - 4, 6)]
for r in roots:
    ru = (2 - r) / (2 * r * r)
    print(f"    x = {r}  ⟹  ρ_u = {ru}  ({'+' if ru > 0 else ''}{float(ru):.4f})")
print("  ⟹ 族內**只有兩個根**：x=2（ρ_u = 0）與 x=2/3（ρ_u = +3/2）")
print(f"     x=2/3 對應之 z_u = {mp.nstr((mp.mpf(3)/2 - RU)/SU,4)}σ ⟹ 已死（docs/10 x=4 之外另一死點）")
print("  ⟹ 因此：在**唯一存活之族**內，ρ_d = 3/8 ⟺ ρ_u = 0（精確等價）")
print("  ⟹ 『δ_u < 0』與『ρ_d = 3/8』在族內 **精確互斥**（無同時解）")

print("\n  代數：ρ_u < 0 ⟺ x > 2（因 2−x < 0、分母恆正）⟹ 此時 ρ_d(x) < 3/8")
for xv in (2.1, 2.2, 2.3, 2.5, 3.0):
    x = mp.mpf(xv)
    ru = (2 - x) / (2 * x * x)
    rd = (2 * x - 1) / (2 * x * x)
    print(f"    x = {xv:.2f}：ρ_u = {mp.nstr(ru,6):>10s}（z_u = {mp.nstr((ru-RU)/SU,6):>8s}）"
          f"   ρ_d = {mp.nstr(rd,6):>9s}（z_d = {mp.nstr((rd-RD)/SD,6):>8s}）")
x_star = mp.findroot(lambda t: (2 - t) / (2 * t * t) - RU, 2.2)
rd_star = (2 * x_star - 1) / (2 * x_star ** 2)
print(f"  使 ρ_u 恰等於觀測之根 x* = {mp.nstr(x_star,10)} ⟹ ρ_d(x*) = {mp.nstr(rd_star,10)}、"
      f"z_d = {mp.nstr((rd_star-RD)/SD,6)}")
print("  ⟹ 欲得觀測之負 ρ_u，須 ρ_d 偏離 3/8 至 −0.85σ（1.34σ 拉低）")

# =====================================================================
print("\n【§三】觀測層：χ²(x) profile（判別力）")
def chi2_at(x):
    ru = (2 - x) / (2 * x * x)
    rd = (2 * x - 1) / (2 * x * x)
    return ((ru - RU) / SU) ** 2 + ((rd - RD) / SD) ** 2


xs = [mp.mpf(1) + mp.mpf(i) / 1000 for i in range(0, 6001)]     # x ∈ [1, 7]
vals = [(chi2_at(x), x) for x in xs]
cmin, xbest = min(vals)
print(f"  最佳 x̂ = {mp.nstr(xbest,8)}（χ²_min = {mp.nstr(cmin,8)}，2 dof ⟹ p = "
      f"{mp.nstr(mp.e**(-cmin/2),4)}）")
print(f"  x̂ > 2 ？{'是（數據偏負 ρ_u）' if xbest > 2 else '否（數據偏正 ρ_u）'}")
print("  1σ 區（Δχ² = 1）與 2σ 區（Δχ² = 4）之邊界：")
for lvl in (1, 4):
    inside = [x for c, x in vals if c - cmin <= lvl]
    print(f"    Δχ² ≤ {lvl}：x ∈ [{mp.nstr(min(inside),6)}, {mp.nstr(max(inside),6)}]")
c_at2 = chi2_at(mp.mpf(2))
print(f"  x=2（ρ_u ≡ 0）之 χ² = {mp.nstr(c_at2,8)} ⟹ 距最小僅 Δχ² = {mp.nstr(c_at2-cmin,6)}")
print(f"  強制 ρ_u = 觀測（x = x*）之 χ² = {mp.nstr(chi2_at(x_star),8)} ⟹ "
      f"Δχ² = {mp.nstr(chi2_at(x_star)-cmin,6)}（< 1 ⟹ 不可裁決）")
print(f"  ⟹ 數據對『ρ_u 之符號』之判別力：**不足**（兩分支 Δχ² < 1）")
print(f"  ⟹ 負號若要達 3σ，需 σ(δ_u) ≲ {mp.nstr(DU/3,4)} MeV（現 {mp.nstr(S_DU,4)}，需改善 "
      f"{mp.nstr(S_DU/(DU/3),3)}×）")

# =====================================================================
print("\n" + "=" * 78)
print("Q3 判決")
print("=" * 78)
print(f"""
  §一 𝒜 層：δ_x ≡ 0 精確 ⟹ **負號不可由 𝒜 生成**（0 無符號）
  §二 α-模型層：族內 ρ_d = 3/8 ⟺ ρ_u = 0（x=2 為唯一存活根；另一根 x=2/3 給 +3/2，已死）
        ⟹ 『ρ_u < 0』與『ρ_d = 3/8』**精確互斥**；欲得負號須 ρ_d 降至 0.344（z_d = −1.34σ）
  §三 觀測層：|δ_u|/σ = {mp.nstr(DU/S_DU,3)} ⟹ 負號**從未達 3σ**；χ² 對符號無判別力（Δχ² < 1）

  ⟹ **Q3 = ①「現精度不可裁決」＋ ②「結構生成 = 死」**
    ① 負號是否存在：不可裁決（{mp.nstr(DU/S_DU,3)}σ；G1 式『嚴格負』從非可裁決要求）
    ② 若存在，可否由結構產生：**死**。範圍（A5）：{{族：𝒜 全部算子}}（§一）＋
       {{有限集：α=(1,x,x²) 族內之兩個 ρ_d=3/8 根}}（§二）
  ⟹ 語言層判定：本語言（A₂ 方向 ω₂ ⊕ 單態）**結構上只能給 0**；
     負號（若有）屬**語言外殘餘（Type-3：語言不適用）**，其地位與 M7 之 A_eff = 4.4492 MeV 同級
     —— 即外部輸入，非可由現有算子生成。
  ⟹ 反過來警告（禁區 4）：**不得**為容納負號而加 `c·sgn` 型自由參數（dof +1，違 A3），
     亦不得改 ω₂ 條件以求命中（ρ_u = 0.77σ 相容零，係結構成果，不可丟）。
""")
raise SystemExit(0)