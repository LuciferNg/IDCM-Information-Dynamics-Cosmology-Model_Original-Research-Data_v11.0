#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/10_family_x_scan.py — 家族 α=(1,x,x²) 之**註冊值**掃描（完備性審計產物）

由審計發現：`docs/04` 之候選 `α=(1,2,4)` 並非孤立 —— 佢係家族
    α = (1, x, x²)   （即 α_ℓ ∝ x^{ℓ−1}，ℓ = 0,1,2）
之 **x=2** 成員；而 `docs/06` 只測過 x=B=4 一個點（樹層分枝），從未掃過家族本身。

x 之**註冊**候選（全部來自既有常數，非自由參數）：
  x=2  ：`log₂B = 2` ＝ `|ℤ₂| = 2`（grading 群階）＝ `N_gen − 1`
  x=3  ：`N_gen = 3`
  x=4  ：`B = 4` ＝ `γ44 = 4`
（`x=1` 作控制：均勻權重）

判準（跑前凍結，承 docs/02 §T5 與 docs/01 §Q2）：
  z_x = (ρ_x − ρ_x^obs)/σ_x ；χ² = z_u² + z_d² ；σ_ρd = 0.044785178 、σ_ρu = 0.033583899
  N_trials = 4（x=1,2,3,4），全報，不得只報命中。
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
RD = mp.mpf("0.40350365")
RU = mp.mpf("-0.02588837")
SD = mp.mpf("0.044785178")
SU = mp.mpf("0.033583899")


def rho(x):
    """α=(1,x,x²) ⟹ 由 A₂ Cartan 方向：ρ_u=(2α_u−α_d)/(2α_s)、ρ_d=(2α_d−α_u)/(2α_s)。"""
    au, ad, as_ = F(1), F(x), F(x*x)
    return (2*au-ad)/(2*as_), (2*ad-au)/(2*as_)


print("=== 0. 家族 α=(1,x,x²) ；x 之註冊候選 ===")
print("  x=2 = log₂B = |ℤ₂|；  x=3 = N_gen；  x=4 = B = γ44；  x=1 = 控制（均勻）")
print(f"  觀測：ρ_d = {mp.nstr(RD,9)} ± {mp.nstr(SD,7)} ； ρ_u = {mp.nstr(RU,9)} ± {mp.nstr(SU,7)}")

print("\n=== 1. 全報（N_trials = 4）===")
print(f"  {'x':>3s} {'α_s/α_u':>9s} {'ρ_u(精確)':>14s} {'ρ_d(精確)':>12s} {'z_u':>8s} {'z_d':>8s} {'χ²':>8s}  判")
rows = []
for x in (1, 2, 3, 4):
    ru, rd = rho(x)
    zu = (mp.mpf(ru.numerator)/ru.denominator - RU)/SU
    zd = (mp.mpf(rd.numerator)/rd.denominator - RD)/SD
    chi2 = zu*zu + zd*zd
    verdict = "✅ <1σ 雙通" if chi2 < 4 else ("🟡 2σ" if chi2 < 9 else "🔴 排除")
    rows.append((x, F(x*x, 1), ru, rd, zu, zd, chi2))
    print(f"  {x:>3d} {x*x:>9d} {str(ru):>14s} {str(rd):>12s} "
          f"{mp.nstr(zu,5):>8s} {mp.nstr(zd,5):>8s} {mp.nstr(chi2,5):>8s}  {verdict}")

print("\n=== 2. 相對判別（以最佳者為基準）===")
best = min(rows, key=lambda r: r[6])
_p = mp.e**(-best[6]/2)          # χ² 2 dof 之尾機率
print(f"  最佳 = x={best[0]}（χ² = {mp.nstr(best[6],6)}，2 dof ⟹ p = {mp.nstr(_p,4)}）")
for (x, _, _, _, _, _, chi2) in rows:
    if x == best[0]:
        continue
    dz = mp.sqrt(abs(chi2-best[6]))
    print(f"  x={best[0]} vs x={x}：Δχ² = {mp.nstr(chi2-best[6],6)} ⟹ 約 {mp.nstr(dz,4)}σ；"
          f"ρ_d 分離 = {mp.nstr(abs(mp.mpf(rows[[r[0] for r in rows].index(x)][3].numerator)/rows[[r[0] for r in rows].index(x)][3].denominator - mp.mpf(best[3].numerator)/best[3].denominator)/SD,5)}σ")

print("\n=== 3. look-elsewhere 與門檻 ===")
print(f"  註冊 x 值共 4 個（含控制）⟹ naive p ≈ 1/4 = 0.25；"
      f"若只計非平凡 x∈{{2,3,4}} ⟹ 1/3 = 0.333")
rd2, rd3 = mp.mpf(rho(2)[1].numerator)/rho(2)[1].denominator, mp.mpf(rho(3)[1].numerator)/rho(3)[1].denominator
gap = abs(rd2-rd3)
print(f"  x=2 vs x=3 之 ρ_d 分離 = {mp.nstr(gap/SD,4)}σ（現況）；3σ 判別需 σ_ρd ≤ {mp.nstr(gap/3,6)} "
      f"⟹ 改善 {mp.nstr(SD/(gap/3),4)}×")

print("\n=== 4. 判 ===")
print("  · x=2 (α=(1,2,4))：ρ=(0, 3/8)，z=(+0.77,−0.64)，χ²≈1.0/2dof ⟹ **雙通（<1σ）**")
print("  · x=4 (α=(1,4,16))：ρ=(−1/16, 7/32)，z_d=−4.13 ⟹ **4.1σ 排除**（此結論與 docs/06 一致，但範圍只限此點）")
print("  · x=3 (α=(1,3,9))：ρ=(−1/18, 5/18)，z_d=−2.81 ⟹ 2.8σ，**未排除亦未通過**")
print("  · x=1（控制）：ρ=(1/2,1/2) ⟹ ρ_u 偏 +15.7σ ⟹ 排除（同 docs/03 之 φ(3)φ(4) 結果）")
print("\n  ⟹ 家族**可裁決**：x=2 為唯一雙通成員；x=4 已死；x=3 於 2.8σ 邊緣（3σ 需 σ_ρd ≤ 0.0324，即 1.38×）。")
print("  ⟹ 殘餘缺口由『α_s/α_u 係邊個整數』改寫為『x=2 之結構來歷（|ℤ₂| vs log₂B vs N_gen−1）』—— 三義同值，需一條裁決。")
sys.exit(0)