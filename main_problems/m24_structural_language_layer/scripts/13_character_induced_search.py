#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/13 — R-a″：非平凡特徵標誘導模（1 維／2 維）
=====================================================
判據（`docs/15_preregistration_round4.md`，跑前凍結）：

  χ ∈ Irr(H)，dim χ ∈ {1,2}；Ind_H^G χ(g) = (1/|H|) Σ_{x: x⁻¹gx∈H} χ(x⁻¹gx)
  分解 m_i = (1/60) Σ_類 |c| χ_G conj(χ_i)（須非負整數）
  構作 D1(k)／D2(k)／D3／D4／D5；統一族 = (固定構作, |H| 或字符索引)
  通過 ⟺ 三值集合 = {16,36,98}

方法誠實標記：字符值用 mpmath dps60（複數）；多重度以 |dev|<1e-30 判整數（明示數值步驟）。
零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations as _perms
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


try:
    import mpmath as mp
    mp.mp.dps = 60
except ImportError:
    raise SystemExit("mpmath required")

E5 = tuple(range(5))


def parity(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2


G = [p for p in _perms(range(5)) if parity(p) == 0]
Gset = set(G)


def mul(a, b):
    return tuple(a[b[i]] for i in range(5))


def invp(a):
    r = [0] * 5
    for i, v in enumerate(a):
        r[v] = i
    return tuple(r)


def order(g):
    x, k = E5, 1
    while True:
        x = mul(x, g)
        if x == E5:
            return k
        k += 1


ORD = {g: order(g) for g in G}
CLS, _seen = [], set()
for g in G:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    CLS.append(c)
    _seen |= c
CLS.sort(key=lambda c: (order(min(c)), min(c)))
CLASS_SIZES = [len(c) for c in CLS]
CLASS_REP = [min(c) for c in CLS]
D = [1, 3, 3, 4, 5]
N_d, B = 44, 4
TARGETS = (16, 36, 98)

# A₅ 特徵標（ℚ(√5) → mpc）
A_ = (Fraction(1, 2), Fraction(1, 2))
B_ = (Fraction(1, 2), Fraction(-1, 2))


def q5mp(x):
    """ℚ(√5) 精確 → mpc（Fraction 逐位轉 mpf，避免 float 截斷）。"""
    a = mp.mpf(x[0].numerator) / x[0].denominator
    b = mp.mpf(x[1].numerator) / x[1].denominator
    return mp.mpc(a + b * mp.sqrt(5), 0)


CHI_G = {
    "χ1": [mp.mpc(1)] * 5,
    "χ3": [mp.mpc(3), mp.mpc(-1), mp.mpc(0), q5mp(A_), q5mp(B_)],
    "χ3'": [mp.mpc(3), mp.mpc(-1), mp.mpc(0), q5mp(B_), q5mp(A_)],
    "χ4": [mp.mpc(4), mp.mpc(0), mp.mpc(1), mp.mpc(-1), mp.mpc(-1)],
    "χ5": [mp.mpc(5), mp.mpc(1), mp.mpc(-1), mp.mpc(0), mp.mpc(0)],
}
DEG = {"χ1": 1, "χ3": 3, "χ3'": 3, "χ4": 4, "χ5": 5}

print("=" * 78)
print("M24 scripts/13 — R-a″：非平凡特徵標誘導模")
print("=" * 78)

# ── 子群（一代表／共軛類）────────────────────────────────────
def closure(gens):
    S, frontier = {E5}, {E5}
    gens = list(gens)
    while frontier:
        new = set()
        for x in frontier:
            for g in gens:
                for z in (mul(x, g), mul(g, x)):
                    if z not in S:
                        S.add(z)
                        new.add(z)
        frontier = new
    return frozenset(S)


subgroups = {frozenset({E5})}
for a in G:
    subgroups.add(closure([a]))
for a, b in combinations(G, 2):
    subgroups.add(closure([a, b]))
subgroups = [S for S in subgroups if S <= Gset]

# 以 Ind_H(1) 之字符去重（共軛子群給同一 Ind）
def perm_char(H):
    Hl = list(H)
    seen, CS = set(), []
    for g in G:
        if g in seen:
            continue
        C = tuple(sorted(mul(g, h) for h in Hl))
        CS.append(set(C))
        seen |= set(C)
    chi = {}
    for g in G:
        chi[g] = sum(1 for c in CS if {mul(g, x) for x in c} == c)
    return chi


REPS = {}
for H in subgroups:
    chi = perm_char(H)
    key = tuple(chi[CLS_REP[k] if False else min(c)] for c in CLS)   # 類代表值
    key = tuple(chi[min(c)] for c in CLS)
    REPS.setdefault(key, (len(H), H))
print(f"  子群類 = {len(REPS)}｜階 = {sorted(v[0] for v in REPS.values())}")

# ── 1 維／2 維特徵標 ─────────────────────────────────────────
def _gens_of(H):
    gens, cur = [], {E5}
    for g in sorted(H):
        if g not in cur:
            gens.append(g)
            cur = set(closure(gens))
    return gens


def _lcm_exp(H):
    from math import lcm
    e = 1
    for g in H:
        e = lcm(e, order(g))
    return e


def one_dim_chars(H):
    """H 之全部一維（可能複）特徵標：以生成元賦值為 N 次單位根，驗乘性。

    N = H 之指數（元素階之 lcm）。字符值用 mpmath dps60（明示數值步驟）。"""
    from itertools import product
    gens = _gens_of(H)
    N = _lcm_exp(H)
    roots = [mp.expjpi(mp.mpf(2) * k / N) for k in range(N)]
    out = {}
    for tup in product(range(N), repeat=len(gens)):
        vals = {E5: mp.mpc(1)}
        ok = True
        for i, g in enumerate(gens):
            vals[g] = roots[tup[i]]
        frontier = list(gens)
        while frontier and ok:
            nxt = []
            for x in frontier:
                for g in gens:
                    z = mul(x, g)
                    v = vals[x] * vals[g]
                    if z in vals:
                        if abs(vals[z] - v) > mp.mpf("1e-40"):
                            ok = False
                            break
                    else:
                        vals[z] = v
                        nxt.append(z)
                if not ok:
                    break
            frontier = nxt
        if ok and len(vals) == len(H):
            key = tuple(mp.nint(mp.re(vals[g]) * 10 ** 12) for g in sorted(H)) + \
                  tuple(mp.nint(mp.im(vals[g]) * 10 ** 12) for g in sorted(H))
            out.setdefault(key, (vals, f"1d_N{N}_{tup}"))
    return list(out.values())


def induce_within(H, K, chiK):
    """H 內部誘導（造 2 維 irrep 用）。"""
    Hl = list(H)
    chi = {}
    for g in Hl:
        s = mp.mpc(0)
        for x in Hl:
            y = mul(mul(invp(x), g), x)
            if y in K:
                s += chiK[y]
        chi[g] = s / len(K)
    return chi


def two_dim_chars(H):
    """2 維 irrep：S₃ = Ind_{Z₂⊂H}(1) − 1 − sign；D₁₀ = χ₃|_H − χ_sign(r↦1, s↦−1)。

    兩者皆以「已知 character 之整數組合」構造，並在 inducer 後由分解檢定其正確性
    （多重度須為非負整數）。"""
    n = len(H)
    out = []
    ones = one_dim_chars(H)
    triv = next((v for v, nm in ones if all(abs(v[g] - 1) < mp.mpf("1e-40") for g in H)), None)
    if n == 6:
        # 取 H 之一個 order-2 子群 Z₂（S₃ 內全共軛）
        z2 = None
        for g in H:
            if order(g) == 2:
                z2 = frozenset({E5, g})
                break
        if z2 is None:
            return out
        chiK = {E5: mp.mpc(1), [g for g in z2 if g != E5][0]: mp.mpc(1)}
        perm = induce_within(H, z2, chiK)          # = 1 + sign + 2dim
        sign = None
        for v, nm in ones:
            if any(abs(v[g] + 1) < mp.mpf("1e-40") for g in H) and all(
                    abs(v[g] - 1) < mp.mpf("1e-40") or abs(v[g] + 1) < mp.mpf("1e-40") for g in H):
                if not all(abs(v[g] - 1) < mp.mpf("1e-40") for g in H):
                    sign = v
                    break
        if sign is not None and triv is not None:
            out.append(({g: perm[g] - triv[g] - sign[g] for g in H}, "S3_std2"))
    if n == 10:
        # χ_sign：r ↦ 1、s ↦ −1（即階2 元素取 −1、階5 取 +1）
        chi_sign = {g: mp.mpc(1 if order(g) == 5 else -1) if g != E5 else mp.mpc(1) for g in H}
        # 檢定乘性（若不乘性則不用）
        okmul = all(abs(chi_sign[mul(a, b)] - chi_sign[a] * chi_sign[b]) < mp.mpf("1e-40")
                    for a in H for b in H)
        if okmul:
            vals = {}
            for g in H:
                idx = [k for k, c in enumerate(CLS) if g in c][0]
                vals[g] = CHI_G["χ3"][idx] - chi_sign[g]
            out.append((vals, "D10_2dim"))
    return out


# ── 誘導 ＋ 分解 ────────────────────────────────────────────
def induce(H, chiH):
    Ns = len(G)
    NH = len(H)
    chi = {}
    for g in G:
        s = mp.mpc(0)
        for x in G:
            y = mul(mul(invp(x), g), x)
            if y in H:
                s += chiH[y]
        chi[g] = s / NH
    return chi


def decompose(chi):
    mult = {}
    for name, cv in CHI_G.items():
        s = mp.mpc(0)
        for k in range(5):
            s += CLASS_SIZES[k] * chi[CLASS_REP[k]] * mp.conj(cv[k])
        s /= 60
        if abs(s.imag) > mp.mpf("1e-30") or abs(s.real - mp.nint(s.real)) > mp.mpf("1e-30"):
            return None, f"non-int mult for {name}: {mp.nstr(s, 20)}"
        mult[name] = int(mp.nint(s.real))
    return mult, None


rows = []
print("\n  誘導模分解表（只列成功者）：")
for key, (nH, H) in sorted(REPS.items(), key=lambda kv: kv[1][0]):
    chars = []
    if nH <= 60:
        chars += one_dim_chars(H)
    chars += two_dim_chars(H)
    for chiH, cname in chars:
        chiG = induce(H, chiH)
        dim = int(mp.nint(chiG[E5].real))
        mult, err = decompose(chiG)
        if err:
            print(f"    |H|={nH:>2} {cname:<12} dim={dim:>2}  ⚠ {err}")
            continue
        if any(v < 0 for v in mult.values()):
            print(f"    |H|={nH:>2} {cname:<12} dim={dim:>2}  ⚠ negative mult {mult}")
            continue
        rows.append({"order": nH, "char": cname, "dim": dim, "mult": mult})

for r in rows:
    print(f"    |H|={r['order']:>2} {r['char']:<12} dim={r['dim']:>2}｜m = {r['mult']}")

# N1 校驗
def D1(r, k):
    return sum(r["mult"].get(nm, 0) * comb(DEG[nm], k) for nm in CHI_G)


def D2(r, k):
    return sum(r["mult"].get(nm, 0) * comb(DEG[nm] + k - 1, k) for nm in CHI_G)


def D3(r):
    it = [(DEG[nm], r["mult"].get(nm, 0)) for nm in CHI_G]
    return sum(it[i][1] * it[j][1] * it[i][0] * it[j][0]
               for i in range(5) for j in range(i + 1, 5))


def D4(r):
    return sum(r["mult"].get(nm, 0) ** 2 * DEG[nm] for nm in CHI_G)


def D5(r):
    return D1(r, 2) + D3(r)


bad = [r for r in rows if D1(r, 1) != r["dim"]]
ck("R-a″-1  全部誘導模 Σ m_i d_i = dim（0 違例）", len(bad) == 0,
   f"violations = {len(bad)}｜表大小 = {len(rows)}")

# 掃描
print("\n  命中（16／36／98）：")
hits = []
for r in rows:
    vals = {"D1": [D1(r, k) for k in (1, 2, 3)], "D2": [D2(r, k) for k in (1, 2, 3)],
            "D3": [D3(r)], "D4": [D4(r)], "D5": [D5(r)]}
    flat = [x for v in vals.values() for x in v]
    got = sorted(set(flat) & set(TARGETS))
    if got:
        hits.append((r["order"], r["char"], r["dim"], got, vals))
        print(f"    |H|={r['order']:>2} {r['char']:<12} dim={r['dim']:>2} → 命中 {got}｜{vals}")
if not hits:
    print("    （無）")

allvals = []
for r in rows:
    allvals += [D1(r, 1), D1(r, 2), D1(r, 3), D2(r, 1), D2(r, 2), D2(r, 3), D3(r), D4(r), D5(r)]
print(f"\n  全表評估值 = {len(allvals)}｜" +
      "｜".join(f"命中 {t} = {sum(1 for x in allvals if x == t)}" for t in TARGETS))

# 統一族：(固定構作, |H| 索引) 或 (固定構作, 字符索引)
triple = []
for cname in ("D1", "D2", "D3", "D4", "D5"):
    for kk in (0, 1, 2):
        fam = []
        for r in rows:
            vals = {"D1": [D1(r, k) for k in (1, 2, 3)], "D2": [D2(r, k) for k in (1, 2, 3)],
                    "D3": [D3(r)], "D4": [D4(r)], "D5": [D5(r)]}[cname]
            fam.append(vals[min(kk, len(vals) - 1)])
        if set(TARGETS) <= set(fam):
            triple.append((cname, kk))
ck("R-a″-2  **無任何** (構作, 索引) 族含齊 {16,36,98}", len(triple) == 0, f"triple = {triple}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/13 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "13_character_induced_search.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "n_modules": len(rows),
               "table": [{"order": r["order"], "char": r["char"], "dim": r["dim"],
                          "mult": r["mult"],
                          "D1": [D1(r, k) for k in (1, 2, 3)],
                          "D2": [D2(r, k) for k in (1, 2, 3)],
                          "D3": D3(r), "D4": D4(r), "D5": D5(r)} for r in rows],
               "hits": [[h[0], h[1], h[2], h[3], {k: v for k, v in h[4].items()}] for h in hits],
               "triple_families": triple,
               "verdict": "R-a″ 見 docs/16"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
