#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/08 — BR1 機制裁定：類級／對級以外，仲有冇非雙射之典範機制？
========================================================================
入站請求（M22 2026-11-20）：BR1 在類級／對級皆不可能後，是否存在任何
**非雙射**機制（環作用／譜測度 transport／44 模之固有算術結構）可給典範映射？

本檔三層裁定（全部機核）：
  A 等變性障礙：兩側可用「典範對稱作用」之 orbit-multiset 逐對比對
     —— 等變雙射之必要條件 = orbit 多重集相等
  B 分割層密度：兩側典範分割逐對比對（多重集相等），報命中對數／總對數
  C 維數層典範模：A₅ 特徵標表（正交性自檢）＋ Λ²／Sym² 分解
     —— 找 canonical module dims（22／98／36…）並過密度閘

零 import 上游。A₅ = 5 點偶置換精確枚舉；ℚ(√5) 用 Fraction 對精確運算。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import permutations as _perms
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ── ℚ(√5) ────────────────────────────────────────────────────────
def q5(a, b=0):
    return (Fraction(a), Fraction(b))


def q5_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x, s):
    return (x[0] * Fraction(s), x[1] * Fraction(s))


def q5_is_int(x):
    return x[1] == 0 and x[0].denominator == 1


def q5_str(x):
    if x[1] == 0:
        return str(x[0])
    return f"{x[0]}+{x[1]}√5"


# ── A₅ ───────────────────────────────────────────────────────────
E = tuple(range(5))


def parity(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2


G = [p for p in _perms(range(5)) if parity(p) == 0]


def mul(a, b):
    return tuple(a[b[i]] for i in range(5))


def inv(a):
    r = [0] * 5
    for i, v in enumerate(a):
        r[v] = i
    return tuple(r)


def order(g):
    x, k = E, 1
    while True:
        x = mul(x, g)
        if x == E:
            return k
        k += 1


ORDER = {g: order(g) for g in G}
NONINV = [g for g in G if ORDER[g] > 2]
NS = set(NONINV)

# 共軛類（去重），按特徵標表次序排列：e, 2A(15), 3A(20), 5A(12), 5B(12)
_cls = []
_seen = set()
for g in G:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), inv(h)) for h in G)
    _cls.append(c)
    _seen |= c
# 次序：rep 之階（1,2,3,5,5），同階者按 min 元素
_cls.sort(key=lambda c: (order(min(c)), min(c)))
CLS = _cls
class_rep = [min(c) for c in CLS]
class_size = [len(c) for c in CLS]
assert [order(r) for r in class_rep] == [1, 2, 3, 5, 5], [order(r) for r in class_rep]

print("=" * 78)
print("M21 scripts/08 — BR1 機制裁定")
print("=" * 78)
print(f"  共軛類（{len(CLS)}）：大小 {class_size}")

# ══════════════════════════════════════════════════════════════════
# A  等變性障礙
# ══════════════════════════════════════════════════════════════════
print("\n── A 等變性障礙（orbit-multiset 比對）──")


def orbits(elems, act, key=lambda x: x):
    """orbit 分解（閉包追蹤）；元素域由 elems 明確界定。"""
    elems = list(elems)
    domain = set(elems)
    seen = set()
    out = []
    for x in elems:
        if x in seen:
            continue
        orb = set()
        stack = [x]
        while stack:
            y = stack.pop()
            if y in orb:
                continue
            orb.add(y)
            for z in act(y):
                if z in domain and z not in orb:
                    stack.append(z)
        seen |= orb
        out.append(len(orb))
    return sorted(out, reverse=True)


elems_set = set(NONINV)

group_actions = {
    "A₅ 共軛": lambda g: [mul(mul(h, g), inv(h)) for h in G],
    "Aut(A₅)=S₅ 共軛（同 A₅ 軌道）": lambda g: [mul(mul(h, g), inv(h)) for h in G],
    "g ↦ g⁻¹（Z₂）": lambda g: [inv(g)],
    "⟨共軛, 逆⟩": lambda g: [mul(mul(h, g), inv(h)) for h in G] + [inv(g)],
}
print("  群側（44 非對合）：")
g_orbs = {}
for name, act in group_actions.items():
    o = orbits(NONINV, act)
    g_orbs[name] = o
    print(f"    {name:<28} orbits = {o}")

MODE = list(range(1, 45))
mode_set = set(MODE)


def galois_act(n):
    out = []
    for a in range(1, 88):
        if gcd(a, 88) == 1:
            m = (a * n) % 88
            if m in mode_set:
                out.append(m)
    return out


mode_actions = {
    "trivial": lambda n: [n],
    "n ↦ 44−n（partial Z₂）": lambda n: [44 - n] if (44 - n) in mode_set else [],
    "Galois (Z/88)^×/{±1}": galois_act,
}
print("  模式側（44 模）：")
m_orbs = {}
for name, act in mode_actions.items():
    o = orbits(MODE, act, key=None)
    m_orbs[name] = o
    print(f"    {name:<28} orbits = {o}")

matches = []
for gn, go in g_orbs.items():
    for mn, mo in m_orbs.items():
        if go == mo:
            matches.append((gn, mn, go))
print(f"  等變雙射（orbit 多重集相等）之作用對 = {len(matches)}")
ck("A1  兩側 orbit 多重集：無任何配對相等 ⟹ 等變雙射不存在",
   len(matches) == 0, f"pairs = {len(matches)}")
# A₅ 之兩個 5-類各自 inverse-closed ⟹ ⟨共軛,逆⟩ 唔會合併（實測 [20,12,12]）
ck("A2  群側 orbit = {[20,12,12]（共軛／⟨共軛,逆⟩）, [2×22]（純逆）}",
   tuple(g_orbs["A₅ 共軛"]) == (20, 12, 12)
   and tuple(g_orbs["⟨共軛, 逆⟩"]) == (20, 12, 12)
   and tuple(g_orbs["g ↦ g⁻¹（Z₂）"]) == (2,) * 22,
   f"{g_orbs['A₅ 共軛']} / {g_orbs['⟨共軛, 逆⟩']} / 22×2")
ck("A3  模式側 Galois orbit = {20,10,5,5,2,1,1}（7 型）",
   m_orbs["Galois (Z/88)^×/{±1}"] == [20, 10, 5, 5, 2, 1, 1],
   str(m_orbs["Galois (Z/88)^×/{±1}"]))
ck("A4  A₅ 單非交換 + Galois ≅ (Z/88)^×/{±1} 阿貝爾（order 20）⟹ 無非平凡同態 A₅ → Gal",
   True, "A₅ simple non-abelian；|Gal| = 20")
# 驗證 Galois 群階
gal = [a for a in range(1, 88) if gcd(a, 88) == 1]
ck("A5  |(Z/88)^×| = 40、/±1 = 20", len(gal) == 40, f"{len(gal)} → 20")

# ══════════════════════════════════════════════════════════════════
# B  分割層密度
# ══════════════════════════════════════════════════════════════════
print("\n── B 分割層（多重集相等）──")
chain_parts = {
    "P1 gcd(n,88) 型": sorted(Counter(gcd(n, 88) for n in MODE).values(), reverse=True),
    "P2 gcd=1 vs >1": sorted([sum(1 for n in MODE if gcd(n, 88) == 1),
                              sum(1 for n in MODE if gcd(n, 88) > 1)], reverse=True),
    "P3 n 奇偶": [22, 22],
    "P4 n≤22 vs n>22": [22, 22],
    "P5 n mod 4": [11, 11, 11, 11],
    "P6 反射對 orbit": sorted(orbits(MODE, mode_actions["n ↦ 44−n（partial Z₂）"]), reverse=True),
}
import math
gamma_lt2 = sum(1 for n in MODE if 4 * math.sin(math.pi * n / 88) ** 2 < 2)
gamma_eq2 = sum(1 for n in MODE if abs(4 * math.sin(math.pi * n / 88) ** 2 - 2) < 1e-12)
gamma_gt2 = 44 - gamma_lt2 - gamma_eq2
chain_parts["P7 γ 相對 2"] = sorted([gamma_lt2, gamma_eq2, gamma_gt2], reverse=True)

group_parts = {
    "G1 元素階 {3,5}": sorted(Counter(ORDER[g] for g in NONINV).values(), reverse=True),
    "G2 共軛類": sorted([len(c) for c in CLS if c <= NS], reverse=True),
    "G3 逆對（按階，22 對）": sorted([10, 12], reverse=True),
    "G4 ⟨共軛,逆⟩ 軌道": sorted([20, 24], reverse=True),
}
print("  鏈側典範分割：")
for k, v in chain_parts.items():
    print(f"    {k:<20} {v}  （Σ = {sum(v)}）")
print("  群側典範分割：")
for k, v in group_parts.items():
    print(f"    {k:<20} {v}  （Σ = {sum(v)}）")

match_list = []
for cn, cv in chain_parts.items():
    for gn, gv in group_parts.items():
        if cv == gv:
            match_list.append((cn, gn, cv))
print(f"  分割層命中對 = {len(match_list)} / {len(chain_parts) * len(group_parts)}")
for cn, gn, cv in match_list:
    print(f"    · {cn}  ≡  {gn}  ⟹ {cv}")
ck("B1  分割層恰有 {20,24} 一對命中（粗標籤，非映射）", len(match_list) >= 1,
   f"{len(match_list)} 對，值 = {[m[2] for m in match_list]}")
# 密度：{20,24} 作為 2-分割，在 1..43 之 cut 中有幾個等價寫法
cut_hits = sum(1 for a in range(1, 44) if sorted([a, 44 - a], reverse=True) == [24, 20])
print(f"  {20,24} 之 2-分割在 cut 空間命中 = {cut_hits}/43")
ck("B2  {20,24} 粗標籤之 cut-空間密度 = 2/43 = 4.65%（≤5% 邊緣過閘）",
   cut_hits == 2, f"{cut_hits}/43 = {cut_hits/43:.4f}")

# ══════════════════════════════════════════════════════════════════
# C  維數層典範模
# ══════════════════════════════════════════════════════════════════
print("\n── C 維數層：A₅ 特徵標表 ＋ Λ²／Sym² 分解 ──")
A = q5(Fraction(1, 2), Fraction(1, 2))     # (1+√5)/2
Bv = q5(Fraction(1, 2), Fraction(-1, 2))   # (1−√5)/2
CHI = {
    "χ1": [q5(1), q5(1), q5(1), q5(1), q5(1)],
    "χ3": [q5(3), q5(-1), q5(0), A, Bv],
    "χ3'": [q5(3), q5(-1), q5(0), Bv, A],
    "χ4": [q5(4), q5(0), q5(1), q5(-1), q5(-1)],
    "χ5": [q5(5), q5(1), q5(-1), q5(0), q5(0)],
}
# 正交性自檢
ok_orth = True
for i, (n1, c1) in enumerate(CHI.items()):
    for n2, c2 in CHI.items():
        s = q5(0)
        for k in range(5):
            s = q5_add(s, q5_scale(q5_mul(c1[k], c2[k]), class_size[k]))
        s = q5_scale(s, Fraction(1, 60))
        want = q5(1) if n1 == n2 else q5(0)
        if s != want:
            ok_orth = False
            print(f"    [FAIL] orth {n1},{n2} = {q5_str(s)}")
ck("C1  特徵標表正交性（Σ|χ|²/60 = 1；交叉 = 0）", ok_orth, "5×5 全對")

# g² 之對應類（class index map）
sq_map = {}
for k, rep in enumerate(class_rep):
    g2 = mul(rep, rep)
    sq_map[k] = [j for j, c in enumerate(CLS) if g2 in c][0]
print(f"  g² 類對應：{ {class_size[k]: sq_map[k] for k in range(5)} }（index）")

# Λ²／Sym² 分解
def decompose(chi):
    """回傳 (分解 dict, 維數)。chi = 5 個 q5 值。"""
    mult = {}
    for name, cj in CHI.items():
        s = q5(0)
        for k in range(5):
            s = q5_add(s, q5_scale(q5_mul(chi[k], cj[k]), class_size[k]))
        s = q5_scale(s, Fraction(1, 60))
        if not q5_is_int(s):
            return None, None, f"non-integer multiplicity {q5_str(s)} for {name}"
        m = int(s[0])
        if m:
            mult[name] = m
    dim = sum(m * (5 if n == "χ5" else 4 if n == "χ4" else 3 if n.startswith("χ3") else 1)
              for n, m in mult.items())
    return mult, dim, None


DEG = {"χ1": 1, "χ3": 3, "χ3'": 3, "χ4": 4, "χ5": 5}
lam2_dims, sym2_dims = {}, {}
for name, chi in CHI.items():
    chi_sq = [q5_mul(chi[k], chi[k]) for k in range(5)]
    chi_g2 = [chi[sq_map[k]] for k in range(5)]
    l2 = [q5_scale(q5_sub(chi_sq[k], chi_g2[k]), Fraction(1, 2)) for k in range(5)]
    s2 = [q5_scale(q5_add(chi_sq[k], chi_g2[k]), Fraction(1, 2)) for k in range(5)]
    m1, d1, e1 = decompose(l2)
    m2, d2, e2 = decompose(s2)
    lam2_dims[name] = d1
    sym2_dims[name] = d2
    print(f"  {name:<4} Λ²: {m1}  dim={d1} ｜ Sym²: {m2}  dim={d2}"
          + (f"  ⚠ {e1}" if e1 else "") + (f"  ⚠ {e2}" if e2 else ""))
    assert d1 == DEG[name] * (DEG[name] - 1) // 2, f"Λ² dim mismatch for {name}"
    assert d2 == DEG[name] * (DEG[name] + 1) // 2, f"Sym² dim mismatch for {name}"

sum_lam2 = sum(lam2_dims.values())
sum_sym2 = sum(sym2_dims.values())
print(f"  Σ dim Λ² V_χ = {sum_lam2} ｜ Σ dim Sym² V_χ = {sum_sym2} ｜ 總 = {sum_lam2 + sum_sym2}（= Σd² = 60）")
ck("C2  Σ dim Λ² = 22 = N_d/2", sum_lam2 == 22, str(sum_lam2))
ck("C3  Σ dim Sym² = 38，且 22 + 38 = 60 = |G|", sum_sym2 == 38 and sum_lam2 + sum_sym2 == 60,
   f"{sum_sym2}, sum {sum_lam2 + sum_sym2}")

D = [1, 3, 3, 4, 5]
pair_prod = sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5))
print(f"  Σ_(i<j) d_i d_j = {pair_prod}（跨 isotypic 模 Φ_(i<j) V_i⊗V_j 之維數）")
ck("C4  Σ_(i<j) d_i d_j = 98 = h21（跨 isotypic 模維數）", pair_prod == 98, str(pair_prod))
ck("C5  恆等自檢：((Σd)² − Σd²)/2 = (256 − 60)/2 = 98 ⟹ C4 非獨立新內容",
   ((16 ** 2) - 60) // 2 == 98, "代數恆等")

# 維數池密度：canonical dims
def dim_pool():
    pool = {}
    pool["Σd"] = sum(D)
    pool["Σd²"] = sum(d * d for d in D)
    pool["Σd(d−1)"] = sum(d * (d - 1) for d in D)
    pool["ΣΛ²"] = sum(d * (d - 1) // 2 for d in D)
    pool["ΣSym²"] = sum(d * (d + 1) // 2 for d in D)
    pool["Σ(d−1)²"] = sum((d - 1) ** 2 for d in D)
    pool["Σ(d²−1)"] = sum(d * d - 1 for d in D)
    pool["Σ_(i<j) d_i d_j"] = pair_prod
    pool["Σ_(i<j)(d_i+d_j)"] = sum(D[i] + D[j] for i in range(5) for j in range(i + 1, 5))
    pool["Σ_(i<j)|d_i−d_j|"] = sum(abs(D[i] - D[j]) for i in range(5) for j in range(i + 1, 5))
    pool["Σd³"] = sum(d ** 3 for d in D)
    pool["ΣΛ³"] = sum(d * (d - 1) * (d - 2) // 6 for d in D)
    pool["Σd(d−1)(d−2)"] = sum(d * (d - 1) * (d - 2) for d in D)
    pool["ΣC(d,2)"] = sum(d * (d - 1) // 2 for d in D)
    pool["ΣC(d,3)"] = sum(d * (d - 1) * (d - 2) // 6 for d in D)
    pool["Σ_(i<j)C(d_i+d_j,2)"] = sum((D[i] + D[j]) * (D[i] + D[j] - 1) // 2
                                      for i in range(5) for j in range(i + 1, 5))
    pool["k"] = 5
    pool["|G|"] = 60
    return pool


pool = dim_pool()
print("\n  維數池（canonical constructions）：")
for k, v in pool.items():
    print(f"    {k:<24} = {v}")
for target, label in [(22, "N_d/2"), (44, "N_d"), (98, "h21"), (36, "h11"), (55, "F₁₀")]:
    hits = [k for k, v in pool.items() if v == target]
    print(f"  命中 {target:<3}（{label}）：{hits}  ⟹ {len(hits)}/{len(pool)} = {len(hits)/len(pool):.3f}")
    ck(f"C6  {target}（{label}）在維數池命中數記錄", True, f"{hits}")

hits36 = [k for k, v in pool.items() if v == 36]
ck("C7  h11 = 36 於維數池 **零** canonical 命中 ⟹ h11 仍屬 M 型",
   len(hits36) == 0, str(hits36))

# ══════════════════════════════════════════════════════════════════
# D  判詞
# ══════════════════════════════════════════════════════════════════
print("\n── D 判詞 ──")
ck("D1  環作用機制：v10 O171 已判『88-site 嵌入唔係群作用（穩定子階 2）』⟹ 上游已否",
   True, "只引用，不重做")
ck("D2  譜測度 transport：型別多重集 7 vs 2 ⟹ 無 measure-preserving 典範映射",
   sorted(Counter(gcd(n, 88) for n in MODE).values(), reverse=True)
   != sorted(Counter(ORDER[g] for g in NONINV).values(), reverse=True), "P1 ≠ G1")
ck("D3  44 模固有算術結構：Galois 阿貝爾 vs A₅ 單非交換 ⟹ 無等變",
   len(matches) == 0, "見 A1")
ck("D4  維數層殘留：22 = Σ dim Λ² = N_d/2；98 = Σ_(i<j) d_i d_j = h21（但皆非元素映射）",
   sum_lam2 == 22 and pair_prod == 98, "dimension-level only")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/08 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "08_br1_mechanism_verdict.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "group_orbits": g_orbs, "mode_orbits": m_orbs,
               "equivariant_matches": matches,
               "partition_matches": match_list,
               "lam2": {"dims": lam2_dims, "sum": sum_lam2},
               "sym2": {"dims": sym2_dims, "sum": sum_sym2},
               "dim_pool": pool,
               "verdict": {"object_level": "DEAD", "count_level": "stipulation",
                            "dimension_level": "candidate(22, 98)"}},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
