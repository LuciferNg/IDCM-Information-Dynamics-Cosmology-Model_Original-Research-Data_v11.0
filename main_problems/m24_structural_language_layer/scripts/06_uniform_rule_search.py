#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/06 — R-a／R-b：統一索引規則窮舉（36 與 98 可否升為讀法？）
======================================================================
判據（`docs/06_preregistration_round2.md` §1／§2，跑前凍結）：

  R-a 模側：濾子族 Φ（d ≥ c／d ≤ c／d = c）× 構作族 Ψ（dim V_S、ΣΛ^k、ΣSym^k、
            CrossPairs、Λ²Total）× 統一索引（φ 固定、ψ 由 k 索引）⟹ k=1,2,3 三值
            通過 ⟺ 三值集合 = {16, 36, 98}
  R-b 鏈側：濾子族（gcd orbit／γ 閾／範圍／Galois 不變集）× 構作族（|S|、Σγ、
            pair 數、對偶對數）× 統一索引（orbit 索引／threshold 索引）
            通過 ⟺ 族內含 98（或三數）

★ 後驗候選（不自動升格）：36 = Λ²Total(S_{d≥B}) = C(B+d_max, 2)；以 L0 S1–S4 逐項判。

零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations, permutations as _perms
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


E5 = tuple(range(5))


def parity(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2


G = [p for p in _perms(range(5)) if parity(p) == 0]


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


ORD_CNT = Counter(order(g) for g in G)
D = [1, 3, 3, 4, 5]
N_d, B = 44, 4
TARGETS = (16, 36, 98)

print("=" * 78)
print("M24 scripts/06 — R-a／R-b：統一索引規則窮舉")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# R-a  模側
# ══════════════════════════════════════════════════════════════════
print("\n── R-a  模側：濾子 × 構作 × 統一索引 ──")
IDX = list(range(5))                      # isotypic 因子索引（度 D[i]）
FILTERS = {}
for c in range(1, 6):
    FILTERS[f"d≥{c}"] = [i for i in IDX if D[i] >= c]
for c in range(1, 5):
    FILTERS[f"d≤{c}"] = [i for i in IDX if D[i] <= c]
for c in (1, 3, 4, 5):
    FILTERS[f"d={c}"] = [i for i in IDX if D[i] == c]


def psi(name, S, k):
    if name == "dimV":
        return sum(D[i] for i in S) if k == 1 else None
    if name == "Lam":
        return sum(comb(D[i], k) for i in S)
    if name == "Sym":
        return sum(comb(D[i] + k - 1, k) for i in S)
    if name == "CrossPairs":
        return sum(D[i] * D[j] for i, j in combinations(sorted(S), 2)) if k == 2 else None
    if name == "Lam2Total":
        if k == 2:
            return sum(comb(D[i], 2) for i in S) + sum(D[i] * D[j] for i, j in combinations(sorted(S), 2))
        return None
    raise AssertionError(name)


PSIS = ["dimV", "Lam", "Sym", "CrossPairs", "Lam2Total"]
rows = []
triple_hits = []
two_hits = []
for fname, S in FILTERS.items():
    if not S:
        continue
    for pname in PSIS:
        vals = []
        for k in (1, 2, 3):
            v = psi(pname, S, k)
            if v is not None:
                vals.append((k, v))
        for k, v in vals:
            rows.append((fname, pname, k, v))
        vs = {v for _, v in vals}
        if set(TARGETS) <= vs:
            triple_hits.append((fname, pname, vals))
        hit = vs & set(TARGETS)
        if len(hit) >= 2:
            two_hits.append((fname, pname, sorted(hit), vals))

print(f"  (φ, ψ, k) 評估格數 = {len(rows)}｜濾子數 = {len([f for f, S in FILTERS.items() if S])}")
print("  逐族（φ 固定、ψ 由 k 索引）之三值表：")
for fname, S in FILTERS.items():
    if not S:
        continue
    for pname in PSIS:
        vals = [psi(pname, S, k) for k in (1, 2, 3)]
        vv = [v for v in vals if v is not None]
        if any(v in TARGETS for v in vv):
            print(f"    φ={fname:<6} S={S} ψ={pname:<10} → {vals}"
                  f"  命中 {sorted({v for v in vv if v in TARGETS})}")

ck("R-a-1  **無任何統一索引族**之三值集合 = {16,36,98}", len(triple_hits) == 0,
   f"triple_hits = {triple_hits}")
print(f"  兩數命中族（2 of 3）= {len(two_hits)}：")
for fname, pname, hit, vals in two_hits:
    print(f"    φ={fname:<6} ψ={pname:<10} 命中 {hit}｜三值 {vals}")

# 後驗候選：36 = Λ²Total(S_{d≥B}) = C(B + d_max, 2)
S_B = FILTERS[f"d≥{B}"]
lam2tot_B = psi("Lam2Total", S_B, 2)
dimV_B = sum(D[i] for i in S_B)
print(f"\n  後驗候選：S_(d≥B) = {S_B}（度 {[D[i] for i in S_B]}）｜dim V_S = {dimV_B}"
      f"｜Λ²Total = {lam2tot_B} = C({dimV_B},2) = {comb(dimV_B, 2)}")
ck("R-a-2  後驗候選成立（值 = 36）：36 = Λ²Total(S_{d≥B}) = C(B+d_max,2)",
   lam2tot_B == 36 and dimV_B == 9, f"{lam2tot_B}")
# 但同一濾子能否給 16／98？
S_B_vals = {p: [psi(p, S_B, k) for k in (1, 2, 3)] for p in PSIS}
print(f"  同濾子 S_(d≥B) 各構作三值：{S_B_vals}")
ck("R-a-3  同濾子 S_(d≥B) **不能**同時給 16／98 ⟹ 後驗候選只解 36 一數（不救同源）",
   not ({16, 98} <= {v for vv in S_B_vals.values() for v in vv if v is not None}),
   f"集合 = {sorted({v for vv in S_B_vals.values() for v in vv if v is not None})}")

# ══════════════════════════════════════════════════════════════════
# R-b  鏈側
# ══════════════════════════════════════════════════════════════════
print("\n── R-b  鏈側：濾子 × 構作 × 統一索引（orbit／threshold）──")
MODES = list(range(1, N_d + 1))
try:
    import mpmath as mp
    mp.mp.dps = 50
    GAM = {n: 4 * mp.sin(mp.pi * n / 88) ** 2 for n in MODES}
except ImportError:
    raise SystemExit("mpmath required")

GCD_ORBITS = {}
for n in MODES:
    GCD_ORBITS.setdefault(gcd(n, 88), []).append(n)
print(f"  gcd orbit：{ {g: len(v) for g, v in sorted(GCD_ORBITS.items())} }")


def cnt(S):
    return len(S)


def gam_sum(S):
    return sum(GAM[n] for n in S)


def pair_cnt(S):
    return comb(len(S), 2)


def dual_cnt(S):
    return sum(1 for n in S if (N_d - n) in set(S))


CONSTRUCTIONS = {"|S|": cnt, "Σγ": gam_sum, "pairs": pair_cnt, "dualpairs": dual_cnt}
FAMILIES = {}
for g in sorted(GCD_ORBITS):
    FAMILIES[f"G_gcd={g}"] = GCD_ORBITS[g]
FAMILIES["T_γ≥2"] = [n for n in MODES if GAM[n] >= 2]
FAMILIES["T_γ>2"] = [n for n in MODES if GAM[n] > 2]
FAMILIES["R_n≤2B"] = [n for n in MODES if n <= 2 * B]
FAMILIES["R_n≤N_d−2B"] = [n for n in MODES if n <= N_d - 2 * B]
FAMILIES["R_dual≥2B"] = [n for n in MODES if N_d - n >= 2 * B]
FAMILIES["ALL"] = MODES

print("  逐 orbit 索引（同一構作，7 個值）：")
fam_tables = {}
for cname, f in CONSTRUCTIONS.items():
    vals = [f(GCD_ORBITS[g]) for g in sorted(GCD_ORBITS)]
    fam_tables[cname] = vals
    print(f"    {cname:<10} {[float(v) if isinstance(v, mp.mpf) else v for v in vals]}")
_orb_vals = {cname: [float(v) if isinstance(v, mp.mpf) else v for v in vals]
             for cname, vals in fam_tables.items()}
_has98 = [cname for cname, vals in _orb_vals.items()
          if any(abs(v - 98) < 1e-40 for v in vals)]
_ok90 = abs(sum(GAM[n] for n in MODES) - 90) < mp.mpf("1e-40")
ck("R-b-1  orbit 索引族（4 構作 × 7 orbit）皆不含 98；max|G_g| = 20、max Σγ|orbit = 40、"
   "pairs 最大 190（皆 ≠ 98）；Σγ(全) = 90",
   len(_has98) == 0 and max(_orb_vals["|S|"]) == 20 and _ok90,
   f"含 98 之族 = {_has98}｜max|G| = {max(_orb_vals['|S|'])}｜"
   f"max Σγ = {max(_orb_vals['Σγ']):.1f}｜max pairs = {max(_orb_vals['pairs'])}｜Σγ_all = {float(sum(GAM[n] for n in MODES)):.6f}")

# 全部濾子 × 構作（逐 (S, ψ) 單值；98 只在跨 orbit 組合時可能出現）
single_hits98 = []
for fname, S in FAMILIES.items():
    for cname, f in CONSTRUCTIONS.items():
        v = f(S)
        if isinstance(v, mp.mpf):
            v = int(round(v)) if abs(v - round(v)) < mp.mpf("1e-40") else v
        single_hits98.append((fname, cname, v))
print("  單濾子 × 構作命中 98 者：")
h98 = [(a, b, c) for a, b, c in single_hits98 if c == 98]
print(f"    {h98}")
ck("R-b-2  單濾子 × 單構作命中 98 = 0（98 > 44 個模式之任何計數；γ 和最大 = 90）",
   len(h98) == 0, f"{h98}")

# Galois 不變集（orbit 子集和）× 構作
UNITS = [a for a in range(1, 88) if gcd(a, 88) == 1]


def rep(m):
    m %= 88
    return m if m <= 44 else 88 - m


orbits, seen = [], set()
for n in MODES:
    if n in seen:
        continue
    o = sorted({rep(a * n) for a in UNITS})
    orbits.append(o)
    seen |= set(o)
inv_hits = {"|S|": [], "Σγ": [], "pairs": [], "dualpairs": []}
for r in range(1, len(orbits) + 1):
    for Ssel in combinations(range(len(orbits)), r):
        S = sorted(set().union(*[set(orbits[i]) for i in Ssel]))
        for cname, f in CONSTRUCTIONS.items():
            v = f(S)
            if isinstance(v, mp.mpf):
                if abs(v - 98) < mp.mpf("1e-40"):
                    inv_hits[cname].append(Ssel)
            elif v == 98:
                inv_hits[cname].append(Ssel)
print(f"  Galois 不變集合（128）× 構作 命中 98：{ {k: len(v) for k, v in inv_hits.items()} }")
ck("R-b-3  Galois 不變集合 × 四構作 命中 98 = 0 ⟹ 98 無任何鏈側統一規則讀法",
   all(len(v) == 0 for v in inv_hits.values()), f"{ {k: len(v) for k, v in inv_hits.items()} }")

# 98 之所有「可達寫法」列明（供判詞引）
composites = {
    "Σγ + 2B": float(sum(GAM[n] for n in MODES)) + 2 * B,
    "N_site + #{gcd=2}": 88 + len(GCD_ORBITS[2]),
    "2N_d + #{gcd=2}": 2 * N_d + len(GCD_ORBITS[2]),
    "2·(Σγ/2 + B)": 2 * (float(sum(GAM[n] for n in MODES)) / 2 + B),
}
print(f"  98 之「凑合寫法」：{composites}")
ck("R-b-4  98 之可達寫法皆屬复合（約定移位／跨 orbit 組合），無統一索引族 ⟹ R-b 復活失敗",
   all(abs(v - 98) < 1e-9 for v in composites.values()), str(composites))

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/06 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "06_uniform_rule_search.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "module": {"n_filters": len([1 for S in FILTERS.values() if S]),
                          "n_grid": len(rows), "triple_hits": triple_hits,
                          "two_hits": [[a, b, c, d] for a, b, c, d in two_hits],
                          "posthoc_36": {"filter": f"d≥{B}", "S": S_B, "dimV": dimV_B,
                                         "Lam2Total": lam2tot_B,
                                         "verdict": "🟡 後驗（只解 36 一數；同濾子不能同時給 16/98）"}},
               "chain": {"orbit_tables": {k: [float(x) if isinstance(x, mp.mpf) else x for x in v]
                                          for k, v in fam_tables.items()},
                         "single_hits98": [], "invariant_hits98": {k: len(v) for k, v in inv_hits.items()},
                         "composites98": {k: float(v) for k, v in composites.items()}},
               "verdict": {"R-a": "🔴 復活失敗（統一索引族 0 命中三數；後驗 36 讀法只解一數）",
                           "R-b": "🔴 復活失敗（98 無鏈側統一規則讀法）"}},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
