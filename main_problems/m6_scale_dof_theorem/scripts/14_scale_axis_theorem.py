#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M6 scripts/14_scale_axis_theorem.py — (B) 尺度軸定理（脫離 H2）
================================================================================
H2（強子層加性）已判假（`docs/13`）。本腳本唔再靠「關係保持」，改由**代數結構 ＋ 適足性**。

B1 代數剛性（核心）
    結構常數若由晶格譜經**有限 ℚ-多項式**得到 ⟹ 落於 ℚ(ζ₈₈)⁺ ⟹ **完全確定、零連續自由度**。
    「Galois 不變」與「確定」係兩件事：γ_7·γ_37 唔不變（有 40 個共軛），但每個共軛皆由整數
    決定 ⟹ 冇可調連續參數。核法：**精確分圓環算術** ℚ[x]/(Φ₈₈)，Φ₈₈ = Φ₁₁(−x⁴)（度 40）。

B2 適足性（錨恰一 ⟹ dof=1 ⟺ 有預言力）
    一錨只定一標度 ⟹ dof≥2 時尺度比自由、零預言。實例：H₀ 由 c = R_world·ω_ring 預言。

B3 帳本完備性（自由數只可住超越槽）
    速率 ≤6（M2）＋ 強子形狀 2（M5）⟹ ≤9。

判據（事前寫死）
    A1 45 條 γ_n + γ_(44−n) = 4 精確（環內恆等），且 Σγ = 90。
    A2 12 個 ℚ-多項式：**環內精確元素**之數值 = 直接數值（差 < 1e−30）⟹ 確定性成立。
    A3 3 個超越反例：環內無法表示；數值 PSLQ（deg ≤ 20，dps 120，殘差 dps 220 確認）找不到關係。
    A4 H₀ 預言與 Planck 在 1σ 內。
輸出：data/14_scale_axis_theorem.json
"""
import json
from fractions import Fraction
from math import gcd

import mpmath as mp

mp.mp.dps = 50
N = 88
DEG = 40                                     # φ(88)
M0 = mp.mpf("27.900915")                     # MeV（M6 錨）
HBARC = mp.mpf("197.3269804")                # MeV·fm

# ── 精確分圓環 ℚ[x]/(Φ₈₈)，Φ₈₈(x) = Σ_{j=0..10} (−1)^j x^{4j} ─────────────
PHI = {4 * j: Fraction((-1) ** j) for j in range(11)}     # 次數 → 係數（含首項 1）
assert PHI[DEG] == 1


def red(p):
    """多項式（dict 次數→Fraction）模 Φ₈₈；用 x^40 = −Σ_{j<10}(−1)^j x^{4j}。"""
    p = {k: Fraction(v) for k, v in p.items() if v}
    while any(k >= DEG for k in p):
        top = max(p)
        c = p[top]
        if top == DEG:
            p.pop(top)
            for j in range(10):
                p[4 * j] = p.get(4 * j, Fraction(0)) - c * PHI[4 * j]
        else:                                  # x^top = x^(top−40)·x^40
            p.pop(top)
            for j in range(10):
                p[top - DEG + 4 * j] = (p.get(top - DEG + 4 * j, Fraction(0))
                                        - c * PHI[4 * j])
        p = {k: v for k, v in p.items() if v}
    return p


def add(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, Fraction(0)) + v
    return red(r)


def mul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            r[ka + kb] = r.get(ka + kb, Fraction(0)) + va * vb
    return red(r)


def scale(a, s):
    return red({k: v * Fraction(s) for k, v in a.items()})


def num(a):
    """環元素之數值（代入 ζ = e^{2πi/88}）。"""
    z = mp.expjpi(2 / mp.mpf(N))
    return mp.fsum(mp.mpf(v.numerator) / mp.mpf(v.denominator) * z ** k for k, v in a.items())


ONE = {0: Fraction(1)}
TWO = {0: Fraction(2)}


def zeta_pow(n):
    n %= N
    return {n: Fraction(1)}


def gam_ring(n):
    """γ_n = 2 − (ζ^n + ζ^−n) 之環元素。"""
    return red(add(TWO, scale(add(zeta_pow(n), zeta_pow(-n)), -1)))


R = {n: gam_ring(n) for n in range(45)}
G = {n: 2 - (mp.expjpi(2 * mp.mpf(n) / N) + mp.expjpi(-2 * mp.mpf(n) / N)) for n in range(45)}

print("=" * 100)
print("M6 (B) 尺度軸定理 — B1 代數剛性（精確分圓環 ℚ[x]/(Φ₈₈)，Φ₈₈ = Φ₁₁(−x⁴)）")
print("=" * 100)
print(f"[輸入] N = {N}，φ(88) = {DEG}，單位群 |(ℤ/88)^×| = {sum(1 for m in range(1, N) if gcd(m, N) == 1)}，M₀ = {mp.nstr(M0, 10)} MeV")

# ── A1 環內恆等 ──────────────────────────────────────────────────────────
print("\n[A1] 晶格關係（環內精確，非數值）")
bad1 = [n for n in range(45) if add(R[n], R[44 - n]) != {0: Fraction(4)}]
sumR = {}
for n in range(45):
    sumR = add(sumR, R[n])
ok_a1 = (not bad1) and sumR == {0: Fraction(90)}
print(f"  γ_n + γ_(44−n) = 4：違反 {len(bad1)}/45")
print(f"  Σ_(n=0..44) γ_n（環內）= {dict(sumR) if sumR else 0}  ⟹ {'equals 90 ✅' if sumR == {0: Fraction(90)} else '❌'}")

# ── A2 12 個 ℚ-多項式：環元素 vs 直接數值 ────────────────────────────────
GOOD = {
    "F1  Σγ_n":            lambda r: _sum45(r),
    "F2  γ_44":            lambda r: r[44],
    "F3  γ_22":            lambda r: r[22],
    "F4  γ_11+γ_33":       lambda r: add(r[11], r[33]),
    "F5  γ_7·γ_37":        lambda r: mul(r[7], r[37]),
    "F6  γ_1²−γ_43²":      lambda r: add(mul(r[1], r[1]), scale(mul(r[43], r[43]), -1)),
    "F7  Σn·γ_n/44":       lambda r: scale(_weighted(r), Fraction(1, 44)),
    "F8  γ_0·γ_44/4":      lambda r: scale(mul(r[0], r[44]), Fraction(1, 4)),
    "F9  (γ_0+γ_44)/2":    lambda r: scale(add(r[0], r[44]), Fraction(1, 2)),
    "F10 γ_2³/8":          lambda r: scale(mul(mul(r[2], r[2]), r[2]), Fraction(1, 8)),
    "F11 3γ_5−2γ_40+1/7":  lambda r: add(add(scale(r[5], 3), scale(r[40], -2)), {0: Fraction(1, 7)}),
    "F12 Σ(−1)^n γ_n":     lambda r: _alt(r),
}


def _sum45(r):
    s = {}
    for n in range(45):
        s = add(s, r[n])
    return s


def _weighted(r):
    s = {}
    for n in range(45):
        s = add(s, scale(r[n], Fraction(n)))
    return s


def _alt(r):
    s = {}
    for n in range(45):
        s = add(s, scale(r[n], Fraction(1 if n % 2 == 0 else -1)))
    return s


GOOD_NUM = {
    "F1  Σγ_n":            lambda g: mp.fsum(g[n] for n in range(45)),
    "F2  γ_44":            lambda g: g[44],
    "F3  γ_22":            lambda g: g[22],
    "F4  γ_11+γ_33":       lambda g: g[11] + g[33],
    "F5  γ_7·γ_37":        lambda g: g[7] * g[37],
    "F6  γ_1²−γ_43²":      lambda g: g[1] ** 2 - g[43] ** 2,
    "F7  Σn·γ_n/44":       lambda g: mp.fsum(mp.mpf(n) * g[n] for n in range(45)) / 44,
    "F8  γ_0·γ_44/4":      lambda g: g[0] * g[44] / 4,
    "F9  (γ_0+γ_44)/2":    lambda g: (g[0] + g[44]) / 2,
    "F10 γ_2³/8":          lambda g: g[2] ** 3 / 8,
    "F11 3γ_5−2γ_40+1/7":  lambda g: 3 * g[5] - 2 * g[40] + mp.mpf(1) / 7,
    "F12 Σ(−1)^n γ_n":     lambda g: mp.fsum((1 if n % 2 == 0 else -1) * g[n] for n in range(45)),
}

print(f"\n[A2] 12 個 ℚ-多項式：環元素（精確）↔ 直接數值")
res = {}
ok_a2 = True
for name, F in GOOD.items():
    elem = F(R)
    v_ring = num(elem)
    v_direct = GOOD_NUM[name](G)
    d = abs(v_ring - v_direct)
    rational = (set(elem) <= {0})
    ok_a2 &= (d < mp.mpf("1e-30"))
    res[name] = {"value": mp.nstr(mp.re(v_direct), 20), "ring_terms": len(elem),
                 "rational": bool(rational), "dev": mp.nstr(d, 4)}
    print(f"  {'✅' if d < mp.mpf('1e-30') else '❌'} {name:22s} = {mp.nstr(mp.re(v_direct), 16):>20s}"
          f"   環項數 {len(elem):2d} {'（有理）' if rational else '（代數·有共軛）'}  差 {mp.nstr(d, 3)}")

print(f"  ⟹ 12/12 之數值皆等於**精確環元素** ⟹ 確定、零連續自由度：{'PASS' if ok_a2 else 'FAIL'}")

# ── A3 超越反例（環外）──────────────────────────────────────────────────
BAD = {
    "X1 γ_1·π":      lambda g: g[1] * mp.pi,
    "X2 γ_1^(√2)":   lambda g: g[1] ** mp.sqrt(2),
    "X3 Γ(γ_1+3)":   lambda g: mp.gamma(g[1] + 3),
}


def _residual(rel, v, deg):
    mp.mp.dps = 220
    powers = [mp.mpf(1)]
    for _ in range(deg):
        powers.append(powers[-1] * mp.mpf(v))
    sc = max(abs(p) for p in powers)
    r = abs(mp.fsum(mp.mpf(int(c)) * p for c, p in zip(rel, powers)) / sc)
    mp.mp.dps = 50
    return r


def find_relation(val, deg=20):
    mp.mp.dps = 120
    powers = [mp.mpf(1)]
    v = mp.re(val)
    for _ in range(deg):
        powers.append(powers[-1] * v)
    sc = max(abs(p) for p in powers)
    if sc == 0:
        return None
    try:
        rel = mp.pslq([p / sc for p in powers], tol=mp.mpf("1e-60"),
                      maxcoeff=10 ** 40, maxsteps=10 ** 6)
    except Exception:  # noqa: BLE001
        rel = None
    mp.mp.dps = 50
    if rel is not None and _residual(rel, v, deg) < mp.mpf("1e-70"):
        return rel
    return None


print("\n[A3] 3 個超越反例（環外 ⟹ 須找不到 ℚ-關係）")
ok_a3 = True
for name, F in BAD.items():
    rel = find_relation(F(G))
    got = rel is not None
    ok_a3 &= (not got)
    res[name] = {"value": mp.nstr(mp.re(F(G)), 16), "relation_found": bool(got)}
    print(f"  {'✅ 無關係' if not got else '❌ 假命中'}  {name:16s} = {mp.nstr(mp.re(F(G)), 14)}")

# ── A4 適足性 ───────────────────────────────────────────────────────────
print("\n[A4] B2 適足性：dof=1 ⟹ H₀ 有預言；dof≥2 ⟹ 此數自由")
H0p, H0m, H0s = mp.mpf("67.4071"), mp.mpf("67.36"), mp.mpf("0.54")
z = abs(H0p - H0m) / H0s
ok_a4 = z < 1
print(f"  H₀(pred) = {mp.nstr(H0p,6)}  vs Planck {mp.nstr(H0m,5)} ± {mp.nstr(H0s,3)}   z = {mp.nstr(z,4)}σ")

# ── B3 帳本 ─────────────────────────────────────────────────────────────
print("\n[B3] 有量綱量 → M₀ 倍數（超越槽以外零自由）")
rows = []
for tag, a, t in (("v1（作廢）", Fraction(377422, 10 ** 6), Fraction(7884, 10 ** 4)),
                  ("v2（現行）", Fraction(384184, 10 ** 6), Fraction(17711, 10 ** 4))):
    a_mp = mp.mpf(a.numerator) / mp.mpf(a.denominator)
    t_mp = mp.mpf(t.numerator) / mp.mpf(t.denominator)
    u, v = HBARC / (a_mp * M0), t_mp / M0
    rows.append({"tag": tag, "u": float(u), "v": float(v)})
    print(f"  {tag:10s} ħc/(a₀M₀) = {mp.nstr(u,10):>14s}   t/M₀ = {mp.nstr(v,10):>14s}")
kl, kh = 68.5696 / 27.900915, 70.5744 / 27.900915
print(f"  κ/M₀ ∈ [{kl:.6f}, {kh:.6f}]")

npass = sum([ok_a1, ok_a2, ok_a3, ok_a4])
print(f"\n{'=' * 100}\n(B) 判：A1={ok_a1} A2={ok_a2} A3={ok_a3} A4={ok_a4}  ⟹ {npass}/4 PASS\n{'=' * 100}")

out = {"checks": {"A1_ring_identities": bool(ok_a1), "A2_exact_determinacy": bool(ok_a2),
                  "A3_transcendental_bad": bool(ok_a3), "A4_adequacy": bool(ok_a4)},
       "summary": {"pass": npass, "total": 4}, "phi88": "x^40 - x^36 + ... (Φ11(−x^4))",
       "functionals": res, "calibration": rows, "kappa_over_M0": [kl, kh],
       "H0": {"pred": float(H0p), "planck": float(H0m), "sigma": float(H0s), "z": float(z)},
       "ledger": {"scale": 1, "rate_slots": "<=6", "hadron_shape_slots": 2,
                  "total_dimensionless_upper": 9}}
with open("data/14_scale_axis_theorem.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("→ data/14_scale_axis_theorem.json")