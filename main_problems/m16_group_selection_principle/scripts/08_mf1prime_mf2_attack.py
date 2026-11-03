#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/08_mf1prime_mf2_attack.py — M16 追加輪 3：攻 MF1′（等距共振之依據）與 MF2（perfect 之獨立痕跡）

§A  MF2：perfect（k=1）是否被觀測測定？
    IDCM：Ω_b = (B−1)/(B³−B) = 3/60 = 1/20；若群階被 abelian 膨脹 k 倍 ⟹ Ω_b = 1/(20k)
    ⟹ 由 Planck Ω_b 反解 k，判 k=1 vs k≥2 之顯著度
§B  MF1′：γ 等差 ⟺ (對極互補 ∧ 中項半最大值) —— 精確等價
§C  MF1′：光譜刻劃族（10 條）在 56 個帶三元組上之檢定 ⟹ {3,4,5} 是否被過度決定
"""
import os, sys
from fractions import Fraction
import sympy as sp

NPASS = NFAIL = 0
OUT = []


def chk(n, c, d=""):
    global NPASS, NFAIL
    if c:
        NPASS += 1
        OUT.append(f"PASS {n}" + (f"  [{d}]" if d else ""))
    else:
        NFAIL += 1
        OUT.append(f"FAIL {n}" + (f"  [{d}]" if d else ""))


# ─────────────────────────────────────────────────────────────
OUT.append("§A MF2：k = |G:G'| 可否被 Ω_b 測定？")
# 精確結構值
OmB_struct = Fraction(3, 60)
OUT.append(f"  IDCM 結構值 Ω_b = 3/60 = {OmB_struct} = {float(OmB_struct):.6f}；若膨脹 k 倍 ⟹ Ω_b = 1/(20k)")
# Planck 2018（TT,TE,EE+lowE+lensing）
om_b, s_om = 0.02237, 0.00015        # Ω_b h²
h, s_h = 0.6736, 0.0054
Om = om_b / h ** 2
s_Om = Om * (((s_om / om_b) ** 2 + (2 * s_h / h) ** 2) ** 0.5)
OUT.append(f"  Planck 2018：Ω_b h² = {om_b} ± {s_om}，h = {h} ± {s_h}")
OUT.append(f"  ⟹ Ω_b = {Om:.6f} ± {s_Om:.6f}")
k_hat = 1 / (20 * Om)
s_k = k_hat * s_Om / Om
OUT.append(f"  ⟹ 反解 k = 1/(20Ω_b) = {k_hat:.4f} ± {s_k:.4f}")
for k in [1, 2, 3, 4]:
    z = (k_hat - k) / s_k
    OUT.append(f"     k={k}: z = {z:+.2f}σ" + ("   ← 唯一相容" if abs(z) < 2 else ""))
chk("A:k=1-compatible", abs((k_hat - 1) / s_k) < 2, f"{(k_hat-1)/s_k:+.2f}σ")
chk("A:k=2-excluded", abs((k_hat - 2) / s_k) > 20, f"{(k_hat-2)/s_k:+.1f}σ")
# IDCM 自身 H0 版本
h_idcm = 0.674071
Om_idcm = om_b / h_idcm ** 2
k_idcm = 1 / (20 * Om_idcm)
OUT.append(f"  用 IDCM 結構 H₀ = {h_idcm}：Ω_b = {Om_idcm:.6f} ⟹ k = {k_idcm:.4f}（與 Planck 版差 {abs(k_idcm-k_hat)/s_k:.2f}σ）")
chk("A:two-h-agree", abs(k_idcm - k_hat) / s_k < 0.5, f"{abs(k_idcm-k_hat)/s_k:.2f}σ")

# ─────────────────────────────────────────────────────────────
OUT.append("§B MF1′：γ-等差 ⟺ (對極互補 ∧ 中項半最大值)？")
n = sp.symbols("n", integer=True, positive=True)
gam = lambda m: sp.simplify(2 * (1 - sp.cos(m * sp.pi / 8)))
# 對極互補：γ_n + γ_(8−n) = γ_max = 4
comp = sp.simplify(gam(n) + gam(8 - n))
OUT.append(f"  對極互補：γ_n + γ_(8−n) = {comp}  （= γ_max = 4，恆等）")
chk("B:complement-identity", comp == 4)
# 中項半最大：γ_4 = 2
chk("B:mid-half", gam(4) == 2)
# 等差 ⟺ 互補 ∧ 中項（在三連續帶之情形）
for B in range(2, 8):
    lhs = sp.simplify(gam(B - 1) + gam(B + 1) - 2 * gam(B))
    ap = (lhs == 0)
    rhs = sp.simplify(gam(B - 1) + gam(B + 1) - (gam(4) + gam(4)))
    OUT.append(f"  B={B}: 帶 {{{B-1},{B},{B+1}}} 等差={ap}；γ_(B−1)+γ_(B+1)−2γ_max/2 = {sp.N(rhs,8)}")
OUT.append("  ⟹ 三連續帶中，等差 ⟺ 中項 = 帶 4（半最大值）⟺ 兩端恰為對極互補對 {3,5}")

# ─────────────────────────────────────────────────────────────
OUT.append("§C MF1′：光譜刻劃族之 look-elsewhere（56 個帶三元組）")
import itertools
gv = {m: sp.N(gam(m), 30) for m in range(1, 9)}
gmax = sp.Integer(4)
triples = list(itertools.combinations(range(1, 9), 3))
TARGET = (3, 4, 5)


def f(x):
    return float(x)


CONDS = {
    "C1 γ 等差（γ1+γ3=2γ2）": lambda a, b, c: abs(f(gv[a] + gv[c] - 2 * gv[b])) < 1e-20,
    "C2 對極互補（γ1+γ3=γ_max）": lambda a, b, c: abs(f(gv[a] + gv[c] - 4)) < 1e-20,
    "C3 中項半最大（γ2=γ_max/2）": lambda a, b, c: abs(f(gv[b] - 2)) < 1e-20,
    "C4 ω 等差（√γ1+√γ3=2√γ2）": lambda a, b, c: abs(f(sp.sqrt(gv[a]) + sp.sqrt(gv[c]) - 2 * sp.sqrt(gv[b]))) < 1e-20,
    "C5 γ 幾何（γ2²=γ1γ3）": lambda a, b, c: abs(f(gv[b] ** 2 - gv[a] * gv[c])) < 1e-20,
    "C6 γ 倒數等差": lambda a, b, c: abs(f(1 / gv[a] + 1 / gv[c] - 2 / gv[b])) < 1e-20,
    "C7 三帶索引連續": lambda a, b, c: (b - a == 1 and c - b == 1),
    "C8 γ 和為整數": lambda a, b, c: abs(f(gv[a] + gv[b] + gv[c]) - round(f(gv[a] + gv[b] + gv[c]))) < 1e-12,
    "C9 兩端和為整數": lambda a, b, c: abs(f(gv[a] + gv[c]) - round(f(gv[a] + gv[c]))) < 1e-12,
    "C10 全為 sin² 之精確角（n=4,8 型）": lambda a, b, c: all(m in (4, 8) for m in (a, b, c)),
}
for name, fcn in CONDS.items():
    hits = [t for t in triples if fcn(*t)]
    has_target = TARGET in hits
    OUT.append(f"  {name:26s} 命中 {len(hits):2d} 組 {'∋目標' if has_target else '∌目標'}   "
               + (f"{hits}" if len(hits) <= 6 else ""))
score = {t: sum(1 for fcn in CONDS.values() if fcn(*t)) for t in triples}
best = max(score.values())
argbest = [t for t in triples if score[t] == best]
OUT.append(f"  每組滿足之條件數：{TARGET} → {score[TARGET]}（最大 = {best}，達到者 {argbest}）")
chk("C:target-maximal", score[TARGET] == best and argbest == [TARGET],
    f"score={score[TARGET]}, argmax={argbest}")
n_ap_all = [t for t in triples if CONDS["C1 γ 等差（γ1+γ3=2γ2）"](*t)]
n_ap_cons = [t for t in triples if CONDS["C1 γ 等差（γ1+γ3=2γ2）"](*t) and CONDS["C7 三帶索引連續"](*t)]
OUT.append(f"  ★ 更正（相對 docs/06）：γ-等差單獨有 {len(n_ap_all)} 組 = {n_ap_all}（全部形如 (4−k,4,4+k)）")
OUT.append(f"     惟「γ-等差 ∧ 三帶連續」唯一 = {n_ap_cons}")
chk("C:ap-alone-3", sorted(n_ap_all) == [(1, 4, 7), (2, 4, 6), (3, 4, 5)])
chk("C:ap-and-consecutive-unique", n_ap_cons == [TARGET])
OUT.append(f"  ★ ω-等差（√γ）在帶 1..8 內**零解** ⟹ 「等距」只可能在 γ 幣制中存在（γ 方為預算幣制，Σγ=90）")
OUT.append(f"  ⟹ 粗 look-elsewhere：56 個三元組中，恰 1 組同時滿足最多之自然刻劃 = {TARGET}")
OUT.append(f"     若盲選三元組：命中 C1 之概率 = 3/56 = {3/56:.4f}；命中「連續 ∧ C1」= 1/6 = {1/6:.4f}")

print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)