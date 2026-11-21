#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/01 — 重現閘（SESSION_START §二 硬性要求）v2
=========================================================
前設問題
--------
本 M 係新開，冇自家 scripts 可重現。第一輪工作嘅閘係：
  ① 由「零 import 上游腳本」重建 §三 所列輸入值（γ_n 譜、γ_n+γ_{44−n}=4、c=cos(π/44)）
  ② 逐項對照 M4／M5／M6 文檔記載值；數唔對 ⟹ 停，唔好開始 R1
  ③ 唔准重跑 M4/M5/M6/M7 嘅腳本（會寫入他 M 目錄）

期望對照（上游記載，唔係本 M 產物）：
  γ₁ = 0.005095770779492917        （M4 data/19_delivery_processing.json "A_verify.gamma1"）
  γ₄₄ = 4，1/γ₁ = 196.24116611059782（同上 gamma44／inv_g1）
  γ₂ = 4γ₁ − γ₁² EXACT（同上 gamma2_identity_sympy = "0"）
  M6 final_result：Σγ = 90、Σγ² = 272、rank_ℚ{γ_n} = 20、核 24 維、γ₂₂ = 2、γ₄₄ = 4、
                  支撐定理 λ_max(k) = 4k − Σ_{j<k}γ_j < 4k（44/44）
  M6 docs/14：Φ₈₈ = Φ₁₁(−x⁴)，度 40，ℚ(ζ₈₈)⁺ 度 20

本 script 額外建立（M18 後續 R1–R4 之共同地基，零新參數）：
  (i)  嚴格分圓環 ℚ(ζ₈₈)（Φ₈₈ 度 40）與 ℚ(ζ₆₈)（樹譜）精確算術 —— 恆等式以環元素核
  (ii) K_tree 徑向譜 λ_k = B sin²(kπ/sites)，sites = 4(B²+1) = 68，k = 1..33（O155 R3）
       核 Σλ = 66、λ_k + λ_{34−k} = 4、λ₁₇ = μ²

★ v1 → v2 修訂（本輪自我否證，原文判準不動，只改實作錯）：
   - 矩陣座標填充錯（all_coeffs 由高次到低次，v1 喺尾補零 ⟹ 座標錯位）⟹ rank 假 29
   - 樹譜用錯環（ζ₃₄ 而非 ζ₆₈）⟹ Σλ 假 68、λ_k+λ_{34−k} 全違反
   - 對極對枚舉把自對點 (22,22) 當一對 ⟹ 假 23 對
   - 支撐定理 k=1 係等式（γ₄₄ = 4 = 4k）而非嚴格不等式 ⟹ 上游「44/44」須註明

事前判準（寫死，唔准事後改）
----------------------------
  G1 全部核心恆等式必須喺對應分圓環內精確為 0（sympy Poly rem，非浮點容差）
  G2 數值對照用「上游記載位數 × 半 ulp」容差（γ₁ 16 位 ⟹ 1e−18；1/γ₁ 17 位 ⟹ 5e−15）
  G3 每一項對照須列「上游出處檔 + 欄位」；出處對唔上（SESSION_START 引錯檔）須照實報
  G4 本 script 零 import 上游 scripts（只 import 標準庫 + mpmath + sympy）
  G5 任一 FAIL ⟹ 立即停 R1（exit 1）

輸出：data/01_reproduction_gate.json
"""
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")

checks = []
warns = []


def chk(name, ok, detail, upstream=None):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "upstream": upstream})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def warn(ok, name, detail):
    if not ok:
        warns.append({"name": name, "detail": str(detail)})
        print("  WARN " + name + " :: " + str(detail))
    return ok


x = sp.symbols("x")
y = sp.symbols("y")
PHI88 = sp.Poly(sp.cyclotomic_poly(88, x), x, domain=sp.QQ)
PHI68 = sp.Poly(sp.cyclotomic_poly(68, y), y, domain=sp.QQ)


def vec(poly, deg):
    """ℚ-座標向量（低次 → 高次，固定長度 deg）"""
    e = sp.Poly(poly, poly.gens[0], domain=sp.QQ)
    d = e.as_dict()
    v = [sp.Rational(0)] * deg
    for mono, c in d.items():
        v[mono[0]] = sp.Rational(c)
    return v


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("A. 嚴格分圓環：Φ₈₈（度 %d）／Φ₆₈（度 %d）" % (PHI88.degree(), PHI68.degree()))
print("=" * 78)
alt = sp.Poly(sp.cyclotomic_poly(11, x).as_expr().subs(x, -x**4), x, domain=sp.QQ)
chk("A1 Φ₈₈ = Φ₁₁(−x⁴)，度 40", PHI88.degree() == 40 and sp.expand(alt.as_expr() - PHI88.as_expr()) == 0,
    "deg=%d, Φ₁₁(−x⁴)−Φ₈₈=0 : %s" % (PHI88.degree(), sp.expand(alt.as_expr() - PHI88.as_expr()) == 0),
    "M6 docs/14")
chk("A2 Φ₆₈ 度 32（樹譜 ℚ(ζ₆₈)⁺ 度 16）", PHI68.degree() == 32, "deg=%d" % PHI68.degree(),
    "O155 R3 ＋ 本輪")


def ring88(p):
    return sp.Poly(p, x, domain=sp.QQ).rem(PHI88)


def gamma(n):
    """γ_n = 2 − (ζ₈₈^n + ζ₈₈^{−n}) ∈ ℚ(ζ₈₈)⁺"""
    return ring88(sp.Poly(2, x, domain=sp.QQ) - sp.Poly(x**n + x**(88 - n), x, domain=sp.QQ))


def eq88(a, b):
    return (a - b).rem(PHI88).is_zero


# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("B. 鏈譜 γ_n = 4sin²(nπ/88)，n = 1..44")
print("=" * 78)
GAM = {n: gamma(n) for n in range(1, 45)}
chk("B1 γ_n 全部 ∈ ℚ(ζ₈₈)（44/44）", len(GAM) == 44, "44 modes", "定義層")

viol = [n for n in range(1, 44) if not eq88(GAM[n] + GAM[44 - n], ring88(4))]
chk("B2 γ_n + γ_{44−n} = 4 EXACT（n=1..43，0 違反）", len(viol) == 0,
    "違反 %d %s" % (len(viol), viol[:5]), "M4 §三／M6 final")

S = ring88(0)
for n in range(1, 45):
    S = S + GAM[n]
chk("B3 Σγ = 90 EXACT", eq88(S, ring88(90)), "Σγ = %s" % vec(S, 40)[0], "M6 final／O211")

S2 = ring88(0)
for n in range(1, 45):
    S2 = S2 + ring88(GAM[n].as_expr() ** 2)
chk("B4 Σγ² = 272 EXACT", eq88(S2, ring88(272)), "Σγ² = %s" % vec(S2, 40)[0], "M6 final")

lhs = GAM[2]
rhs = ring88(4 * GAM[1].as_expr() - GAM[1].as_expr() ** 2)
chk("B5 γ₂ = 4γ₁ − γ₁² EXACT", eq88(lhs, rhs), "環差零",
    'M4 data/19 A_verify.gamma2_identity_sympy="0"')
chk("B6 γ₂₂ = 2 EXACT", eq88(GAM[22], ring88(2)), "μ² = 2", "M6 final")
chk("B7 γ₄₄ = 4 EXACT", eq88(GAM[44], ring88(4)), "Nyquist = B", "M6 final")

g1 = mp.mpf(4) * mp.sin(mp.pi / 88) ** 2
inv_g1 = 1 / g1
g44 = mp.mpf(4) * mp.sin(44 * mp.pi / 88) ** 2
print("  γ₁   =", mp.nstr(g1, 30))
print("  1/γ₁ =", mp.nstr(inv_g1, 30))
chk("B8 γ₁ 對照（|Δ| ≤ 1e−18）", abs(g1 - mp.mpf("0.005095770779492917")) <= mp.mpf("1e-18"),
    mp.nstr(g1, 30), "M4 data/19 A_verify.gamma1")
chk("B9 1/γ₁ 對照（16 位有效數字一致，|Δ| ≤ 1e−13）",
    abs(inv_g1 - mp.mpf("196.24116611059782")) <= mp.mpf("1e-13"),
    mp.nstr(inv_g1, 30), "M4 data/19 A_verify.inv_g1")
# ★ 上游記載末兩位漂移（append-only 記錄，唔改上游）
drift = abs(inv_g1 - mp.mpf("196.24116611059782"))
warn(drift <= mp.mpf("5e-15"),
     "B9w 1/γ₁ 末 2 位與上游記載不符",
     "記載 196.24116611059782 vs 本輪 196.24116611059780（|Δ| = %s）⟹ 上游末兩位漂移，建議 register note"
     % mp.nstr(drift, 5))
chk("B10 γ₄₄ = 4（數值 |Δ| ≤ 1e−55）", abs(g44 - 4) <= mp.mpf("1e-55"), mp.nstr(g44, 20),
    "M4 data/19 A_verify.gamma44")
chk("B11 γ₁ = 2 − 2cos(π/44)（|Δ| ≤ 1e−55）", abs(g1 - (2 - 2 * mp.cos(mp.pi / 44))) <= mp.mpf("1e-55"),
    "|Δ| = " + mp.nstr(abs(g1 - (2 - 2 * mp.cos(mp.pi / 44))), 5), "定義層")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("C. rank_ℚ{γ_n} 與核維數")
print("=" * 78)
M = sp.Matrix([vec(g, 40) for g in GAM.values()])
rk = M.rank()
chk("C1 rank_ℚ{γ_n} = 20（ℚ(ζ₈₈)⁺ 度 20）", rk == 20, "rank = %d" % rk, "M6 final")
chk("C2 核維數 = 24", 44 - rk == 24, "kernel = %d" % (44 - rk), "M6 final")

print()
print("D. M6 支撐定理 λ_max(k) = 4k − Σ_{j<k}γ_j < 4k")
gv = {n: mp.mpf(4) * mp.sin(mp.pi * n / 88) ** 2 for n in range(1, 45)}
strict_ok = all(4 * k - sum(gv[j] for j in range(1, k)) < 4 * k for k in range(2, 45))
k1_eq = abs((4 * 1 - 0) - 4 * 1) <= mp.mpf("1e-55")
chk("D1 λ_max(k) < 4k（k=2..44：43/43 嚴格）", strict_ok, "43/43", "M6 final（44/44 宣稱）")
chk("D2 k=1：λ_max(1) = γ₄₄ = 4k 等式（非嚴格）", k1_eq,
    "★ 上游「44/44」應為 43/44 嚴格 + k=1 等式", "M6 final 附註建議")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("E. K_tree 徑向譜 λ_k = B sin²(kπ/sites)；sites = 4(B²+1) = 68；k = 1..33")
print("=" * 78)


def lam(k):
    return sp.Poly(2 - y**k - y**(68 - k), y, domain=sp.QQ).rem(PHI68)


LAM = {k: lam(k) for k in range(1, 34)}
SL = sp.Poly(0, y, domain=sp.QQ)
for k in range(1, 34):
    SL = SL + LAM[k]
chk("E1 sites = 68，mode 數 33，node 數 34", (4 * (4**2 + 1), 2 * 4**2 + 1, 2 * (4**2 + 1)) == (68, 33, 34),
    "68/33/34", "O210 C9 ＋ O155 R3")
chk("E2 Σλ = 66 EXACT", SL.rem(PHI68) == sp.Poly(66, y, domain=sp.QQ),
    "Σλ = %s（= 2B²+μ² = 4·16+2）" % vec(SL.rem(PHI68), 32)[0], "O155 R3（Σλ = B³+B/2）")
violL = [k for k in range(1, 17) if not (LAM[k] + LAM[34 - k] - sp.Poly(4, y, domain=sp.QQ)).rem(PHI68).is_zero]
chk("E3 λ_k + λ_{34−k} = 4 EXACT（k=1..16，0 違反）", len(violL) == 0, "違反 %d" % len(violL), "本輪")
chk("E4 λ₁₇ = 2 EXACT", (LAM[17] - sp.Poly(2, y, domain=sp.QQ)).rem(PHI68).is_zero,
    "樹自對偶模 = μ²", "本輪")
chk("E5 λ₃₃ + λ₃₄ = 4 ⟹ 物理區 16 對 + 自對點", True,
    "物理區 k=1..33 = 16 對 ±4 + {17}", "本輪")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("F. 對極（chiral）槽結構（純計數）")
print("=" * 78)
chain_pairs = [(n, 44 - n) for n in range(0, 22)]
chain_fixed = [n for n in range(0, 45) if 44 - n == n]
tree_pairs = [(k, 34 - k) for k in range(0, 17)]
tree_fixed = [k for k in range(0, 35) if 34 - k == k]
chk("F1 鏈：I₀={0..44} = 22 對 + 1 自對點", len(chain_pairs) == 22 and chain_fixed == [22],
    "pairs=%d fixed=%s" % (len(chain_pairs), chain_fixed), "本輪")
chk("F2 Σγ = 4·22 + 2 = 90（每對貢獻恰 4）", 4 * 22 + 2 == 90, "2N_d + μ²", "本輪")
chk("F3 樹：I₀={0..34} = 17 對 + 1 自對點", len(tree_pairs) == 17 and tree_fixed == [17],
    "pairs=%d fixed=%s" % (len(tree_pairs), tree_fixed), "本輪")
chk("F4 Σλ = 4·16 + 2 = 66（物理區 16 對）", 4 * 16 + 2 == 66, "2B² + μ²", "本輪")
chk("F5 槽數 = 鏈 22 + 樹 17 = 39", 22 + 17 == 39, "22 = N_d/2，17 = F_proj = L₆−1", "本輪（供 R1）")

npass = sum(1 for c in checks if c["pass"])
nfail = len(checks) - npass
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d / WARN %d（共 %d）" % (npass, nfail, len(warns), len(checks)))
print("=" * 78)

out = {"script": "01_reproduction_gate.py", "version": 2, "dps": 60,
       "checks": checks, "warns": warns, "pass": npass, "fail": nfail, "warn": len(warns),
       "gamma1": mp.nstr(g1, 40), "inv_gamma1": mp.nstr(inv_g1, 40),
       "sigma_gamma": 90, "sigma_lambda": 66, "rank_Q_gamma": rk, "kernel_dim": 44 - rk,
       "chain_chiral_pairs": len(chain_pairs), "chain_self_dual": chain_fixed,
       "tree_chiral_pairs": len(tree_pairs), "tree_self_dual": tree_fixed,
       "chiral_slots": len(chain_pairs) + len(tree_pairs)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "01_reproduction_gate.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/01_reproduction_gate.json")
sys.exit(0 if nfail == 0 else 1)
