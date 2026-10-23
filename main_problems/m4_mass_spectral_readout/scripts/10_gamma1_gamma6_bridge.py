#!/usr/bin/env python3
"""
M4 R2 — γ₁↔γ₆「同步邊界之間嘅邊界」結構橋 (2026-10-24)
==========================================================
問題（SESSION_START §四 R2）：γ₆/γ₁ = 35.47 ≈ ?、γ₁^k γ₆^l·A 混合、γ₁ 階梯「log 等距 2.29279」身份。
禁區：M_GUT/M_PS/M_R = M_Pl γ₁^k 階梯 —— O85 由零重建已證「GUT/PS/R 階梯唔存在，舊污染已刪」
      （o85_su12_breaking/docs/01 §零硬編碼標度; docs/06「Λ 係唯一自由標度」）。

本輪做三件事（全部可精確驗證，唔用 fit）：
  A. 譜代數：γ_n = 4sin²(nπ/88) = 2−2T_n(c),  c = cos(π/44)
     ⟹ 三項遞推 γ_{n+1} = 2c·γ_n − γ_{n−1} + 2γ₁（精確, 用 sympy 證）
     ⟹ 全譜係單一生成元嘅多項式（ℤ[c] = ℤ[γ₁/2]）—— 即「橋」= Chebyshev 恆等，冇新常數
  B. γ₆/γ₁ 精確代數身份：寫成 c 嘅多項式 + 最小多項式 (degree)
     ⟹ 「35.47 ≈ ?」唔係巧合問題, 係代數數身份
  C. 指數階梯族密度：family = {M_Pl·γ₁^k·S}；算 multiplicity（幾多組合可以「解釋」同一目標）
     ⟹ 判 M_EW = M_Pl γ₁⁸ N_d（dev 0.80%）係結構定揀出嚟
  D. 「log 等距 2.29279」= log10(1/γ₁) —— 定義式，唔含額外結構
"""
import json
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 50
L = 88
print("="*80)
print("M4 R2: γ₁↔γ₆ 結構橋 — Chebyshev 譜代數（精確）")
print("="*80)

# ── A. 精確譜恆等 ──
c = sp.symbols('c')
x = sp.symbols('x')
gamma1_expr = sp.simplify(4*sp.sin(sp.pi/L)**2)
print(f"γ₁ = 4sin²(π/88) 精確式 = {gamma1_expr}")
g1 = mp.mpf(4)*mp.sin(mp.pi/L)**2
print(f"γ₁ (50dps) = {mp.nstr(g1, 25)}")

# T_n(cos θ) = cos(nθ) ⟹ 4sin²(nπ/88) = 2 − 2cos(nπ/44) = 2 − 2T_n(c), c=cos(π/44)
Tn = {}
Tn[0] = sp.Integer(1); Tn[1] = c
for k in range(2, 46):
    Tn[k] = sp.expand(2*c*Tn[k-1] - Tn[k-2])

gam_sym = {n: sp.expand(2 - 2*Tn[n]) for n in range(0, 45)}
gam_sym[0] = sp.Integer(0)  # mode 0 = 0

# 數值核對 1..44（50 dps）
cm = mp.cos(mp.pi/44)
maxerr = mp.mpf(0)
for n in range(1, 45):
    lhs = 4*mp.sin(n*mp.pi/L)**2
    rhs = mp.mpf(str(sp.N(gam_sym[n].subs(c, sp.Float(str(cm), 50)), 40)))
    maxerr = max(maxerr, abs(lhs-rhs)/lhs)
print(f"\n[A1] γ_n = 2−2T_n(c) 核對 n=1..44: 最大相對誤差 {mp.nstr(maxerr, 5)} "
      + ("✅ (50 dps 精度內)" if maxerr < mp.mpf('1e-30') else "❌"))

# 遞推 γ_{n+1} = 2c γ_n − γ_{n−1} + 2γ₁
rec_ok = True
g1s = gam_sym[1]
for n in range(1, 44):
    lhs = gam_sym[n+1]
    rhs = sp.expand(2*c*gam_sym[n] - gam_sym[n-1] + 2*g1s)
    if sp.simplify(lhs-rhs) != 0:
        rec_ok = False
        print(f"  ❌ n={n}: 差 {sp.simplify(lhs-rhs)}")
        break
print(f"[A2] 遞推 γ_(n+1) = 2c·γ_n − γ_(n−1) + 2γ₁ : 精確成立 "
      + ("✅ (sympy 對 n=1..43 全部 simplify 至 0)" if rec_ok else "❌"))
print("     ⟹ 全譜由 (c, γ₁) 生成，且 c = 1 − γ₁/2 ⟹ 單一生成元 γ₁ ⟹ 譜內冇獨立自由度")
print("     ⟹ 「γ₁↔γ₆ 橋」= 上述 Chebyshev/Toeplitz 恆等（唔需要、亦唔存在額外結構常數）")

# ── B. γ₆/γ₁ 精確代數身份 ──
ratio6 = sp.simplify(gam_sym[6]/g1s)
print(f"\n[B1] γ₆ = {gam_sym[6]}")
poly, rem = sp.div(sp.Poly(sp.expand(gam_sym[6]), c), sp.Poly(sp.expand(g1s), c))
print(f"     γ₆/γ₁ = 多項式 {poly.as_expr()} + 餘式 {rem.as_expr()}")
r6 = mp.mpf(0)
r6 = 4*mp.sin(6*mp.pi/L)**2/g1
print(f"     γ₆/γ₁ (50dps) = {mp.nstr(r6, 25)}")
print(f"     mode 結構: γ₆/γ₁ = [sin(6π/88)/sin(π/88)]² = (U₅(cosπ/88))², U₅ = 2nd-kind Chebyshev")
# 最小多項式（ℚ 上）: 用 resultant 消除 c —— γ₆/γ₁ = P₅(c), c = cos(π/44), minpoly(c) | T₄₄(c)−1
cc, yy = sp.symbols('c y')
P5 = 32*cc**5 + 32*cc**4 - 16*cc**3 - 16*cc**2 + 2*cc + 2
mc_c = sp.expand(sp.chebyshevt(44, cc) + 1)          # c = cos(π/44) ⟹ T₄₄(c) = cos(π) = −1 ⟹ T₄₄(c)+1 = 0
res = sp.factor(sp.resultant(sp.expand(yy - P5), mc_c, cc))
print(f"\n     resultant((y − P₅(c)), T₄₄(c)+1; c) 分解 = {res}")
# 取數值根 y = r6 所屬嘅不可約因子（用高精度 + 係數尺度正規化）
fac = sp.factor_list(res)[1]
best = None
for f, _e in fac:
    d = int(sp.degree(f, yy))
    if d == 0:
        continue
    val = mp.mpf(str(sp.N(f.subs(yy, sp.Float(str(r6), 46)), 50)))
    cf = sp.Poly(f, yy).all_coeffs()
    scale = sum(abs(float(k))*(float(r6)**(d-i)) for i, k in enumerate(cf))
    rel = abs(val)/mp.mpf(str(scale))
    print(f"     candidate factor deg {d:2d}: |f(r6)|/coeffscale = {mp.nstr(rel,4)}")
    if best is None or rel < best[0]:
        best = (rel, d, f, val)
rel, deg6, mp6, val6 = best
print(f"\n     γ₆/γ₁ 最小多項式 (ℚ 上) degree = {deg6}: {mp6}")
print(f"     核對 |f(r6)|/scale = {mp.nstr(rel,4)}  (值 {mp.nstr(val6,4)})")
print(f"     ⟹ 35.4678… 係 degree-{deg6} 代數數（ℚ(ζ₈₈)⁺ 內嘅 {deg6}-次子域）；唔係「≈35.47 巧合」")

# ── C. 指數階梯族密度 ──
print("\n" + "-"*80)
print("[C] 指數階梯族 density: family = {M_Pl · γ₁^k · S}")
print("-"*80)
M_Pl = mp.mpf('1.2209e19')   # GeV (O89 用值，同 register 一致)
v_SM = mp.mpf('246.2196')    # GeV (IDCM 錨)
S_pool = {
    "1": 1, "B=4": 4, "Nd=44": 44, "L=88": 88, "Ngen=3": 3, "F5=5": 5, "F7=13": 13,
    "F8=21": 21, "L5=11": 11, "L7=29": 29, "L9=76": 76, "2B=8": 8, "Ngen²=9": 9,
    "L6=18": 18, "Nd+B=48": 48, "L+Nd=132": 132, "√5": float(np.sqrt(5)), "φ=1.618": (1+np.sqrt(5))/2,
}
ks = list(range(-45, 46))
combos = []
for k in ks:
    for sn, S in S_pool.items():
        combos.append((k, sn, M_Pl*(g1**k)*mp.mpf(str(S))))
print(f"  族大小: k∈[-45,45]×{len(S_pool)} 結構常數 = {len(combos)} 組合")

def mult(target, thr_pct):
    hits = [(k, sn, val) for k, sn, val in combos if abs(val/target - 1)*100 < thr_pct]
    return hits

for tname, tval in [("v_SM", v_SM), ("M_EW/O89", mp.mpf('244.24'))]:
    for thr in (0.5, 0.8, 1.0):
        h = mult(tval, thr)
        print(f"  目標 {tname} = {mp.nstr(tval,8)} GeV, thr {thr}%: multiplicity = {len(h)}")
    h = mult(tval, 1.0)
    h.sort(key=lambda z: abs(z[2]/tval-1))
    for k, sn, val in h[:8]:
        print(f"      γ₁^{k:<3}×{sn:<10} = {mp.nstr(val,8):>14} GeV  dev {float(abs(val/tval-1)*100):.4f}%")
    print(f"      ⟹ {'🔴 多解（唔可辨識）' if len(h) > 3 else '🟡 族內單解；但 k 由目標反選 ⟹ 唔算預言'}")
# O89 原式
o89 = M_Pl*(g1**8)*44
print(f"\n  O89 原式 M_Pl·γ₁⁸·N_d = {mp.nstr(o89, 10)} GeV;  v_SM dev = {float(abs(o89/v_SM-1)*100):.4f}%;"
      f" O89 自報 244.24/dev 1.05%")
print(f"  ⟹ 指數 8 係事後揀（Λ²(12)=8 由目標反選）; 族內 {len(mult(o89,0.80))} 個組合同樣 hit ⟹ 唔可辨識")

# ── D. log 等距 2.29279 ──
print("\n" + "-"*80)
print("[D] 「γ₁ 階梯 log 等距 2.29279」身份")
print("-"*80)
l10 = mp.log10(1/g1)
print(f"  log10(1/γ₁) = {mp.nstr(l10, 20)}")
print(f"  記錄值 2.29279 → 差 {float(abs(l10-mp.mpf('2.29279'))/mp.mpf('2.29279')*100):.2e}%")
print("  ⟹ 恆等式：2.29279 = log10(1/γ₁) = −log10(4sin²(π/88)) —— 定義式，零額外結構")
# 階梯間距實際上係 log10(γ_k/γ_1)，逐 k 列
print("  各 mode log10(γ_k/γ₁): " + ", ".join(
    f"{k}:{float(mp.log10(4*mp.sin(k*mp.pi/L)**2/g1)):.4f}" for k in (2, 3, 4, 5, 6, 11, 12, 20)))

out = {
    "gamma1": mp.nstr(g1, 30),
    "chebyshev_identity": "γ_n = 2 − 2T_n(c), c = cos(π/44); recurrence γ_(n+1)=2cγ_n−γ_(n−1)+2γ₁",
    "max_rel_err_1_44": mp.nstr(maxerr, 6),
    "recurrence_exact": bool(rec_ok),
    "gamma6_over_gamma1": mp.nstr(r6, 25),
    "minpoly_gamma6_over_gamma1": str(mp6),
    "ladder_family_size": len(combos),
    "v_SM_multiplicity_0.8pct": len(mult(v_SM, 0.80)),
    "log10_inv_gamma1": mp.nstr(l10, 20),
}
with open("data/10_gamma1_gamma6_bridge.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\n→ data/10_gamma1_gamma6_bridge.json")
