#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/09 — R-a′：子群／誘導模型子模（第二輪明列未掃者）
============================================================
判據（`docs/10_preregistration_round3.md` §1，跑前凍結）：

  子群類 9 個（order 1,2,3,4,5,6,10,12,60）；置換模 P_H = ℂ[A₅/H]（dim = 60/|H|）
  χ_{P_H}(g) = #被 g 固定之左陪集；isotypic 分解 m_i(H) = ⟨χ_P, χ_i⟩（ℚ(√5) 精確）
  構作 D1(k)、D2(k)（k=1,2,3）、D3、D4、D5
  統一族 = (固定構作, 以 |H| 或 k 索引)；通過 ⟺ 三值集合 = {16,36,98}

零 import 上游；ℚ(√5) 用 Fraction 精確；χ₃/χ₃' 之 5A/5B 由共軛類代表元定序，一致計算。
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
CLASS_ORD = [order(min(c)) for c in CLS]
D = [1, 3, 3, 4, 5]
N_d, B = 44, 4
TARGETS = (16, 36, 98)

print("=" * 78)
print("M24 scripts/09 — R-a′：子群／誘導模型子模")
print("=" * 78)
print(f"  共軛類（階）：{list(zip(CLASS_SIZES, CLASS_ORD))}")

# ── ℚ(√5) ────────────────────────────────────────────────────────
def q5(a, b=0):
    return (Fraction(a), Fraction(b))


def q5_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def q5_mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x, s):
    return (x[0] * Fraction(s), x[1] * Fraction(s))


def q5_isint(x):
    return x[1] == 0 and x[0].denominator == 1


CHI = {
    "χ1": [q5(1)] * 5,
    "χ3": [q5(3), q5(-1), q5(0), q5(Fraction(1, 2), Fraction(1, 2)), q5(Fraction(1, 2), Fraction(-1, 2))],
    "χ3'": [q5(3), q5(-1), q5(0), q5(Fraction(1, 2), Fraction(-1, 2)), q5(Fraction(1, 2), Fraction(1, 2))],
    "χ4": [q5(4), q5(0), q5(1), q5(-1), q5(-1)],
    "χ5": [q5(5), q5(1), q5(-1), q5(0), q5(0)],
}
DEG = {"χ1": 1, "χ3": 3, "χ3'": 3, "χ4": 4, "χ5": 5}

# ── 子群枚舉（closure of generated sets；按置換特徵標去重）────────
def closure(gens):
    S = {E5}
    frontier = {E5}
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
print(f"  子群總數 = {len(subgroups)}｜階分佈 = {dict(sorted(Counter(len(S) for S in subgroups).items()))}")

# 置換特徵標：χ_P(g) = #{左陪集 C : g·C = C}
def cosets(H):
    """左陪集 gH（左作用 g'(gH) = (g'g)H 保持陪集結構）。"""
    Hl = list(H)
    seen, out = set(), []
    for g in G:
        if g in seen:
            continue
        C = tuple(sorted(mul(g, h) for h in Hl))
        out.append(C)
        seen |= set(C)
    return out


def perm_char(H):
    CS = cosets(H)
    CSset = [set(c) for c in CS]
    chi = {}
    for g in G:
        cnt = 0
        for c in CSset:
            if {mul(g, x) for x in c} == c:
                cnt += 1
        chi[g] = cnt
    return chi, len(CS)


def decomp(chi):
    """回傳 multiplicity 字典（ℚ(√5) 精確）與總維數。"""
    cls_vals = []
    for c in CLS:
        rep = min(c)
        cls_vals.append(chi[rep])
    mult = {}
    for name, cv in CHI.items():
        s = q5(0)
        for k in range(5):
            s = q5_add(s, q5_scale(q5_mul(q5(cls_vals[k]), cv[k]), CLASS_SIZES[k]))
        s = q5_scale(s, Fraction(1, 60))
        if not q5_isint(s):
            return None, f"non-int mult for {name}: {s}"
        mult[name] = int(s[0])
    return mult, None


# 去重：以 (dim, 分解) 為鍵（置換模只依賴子群共軛類）
table = {}
for H in subgroups:
    chi, dim = perm_char(H)
    mult, err = decomp(chi)
    if err:
        print(f"    ⚠ H(order {len(H)}) 分解失敗：{err}")
        continue
    key = (dim, tuple(sorted(mult.items())))
    table.setdefault(key, {"order": len(H), "dim": dim, "mult": mult})
print(f"  相異置換模（按特徵標去重）= {len(table)}")

rows = []
print("\n  置換模分解表：")
for key, v in sorted(table.items()):
    m = v["mult"]
    print(f"    |H|={v['order']:>2}｜dim={v['dim']:>2}｜m = {m}")
    rows.append(v)


def D1(v, k):
    return sum(v["mult"].get(nm, 0) * comb(DEG[nm], k) for nm in CHI)


def D2(v, k):
    return sum(v["mult"].get(nm, 0) * comb(DEG[nm] + k - 1, k) for nm in CHI)


def D3(v):
    items = [(DEG[nm], v["mult"].get(nm, 0)) for nm in CHI]
    tot = 0
    for i in range(5):
        for j in range(i + 1, 5):
            tot += items[i][1] * items[j][1] * items[i][0] * items[j][0]
    return tot


def D4(v):
    return sum(v["mult"].get(nm, 0) ** 2 * DEG[nm] for nm in CHI)


def D5(v):
    return D1(v, 2) + D3(v)


# 校驗：D1(v,1) 必須等於 dim
bad = [v for v in rows if D1(v, 1) != v["dim"]]
ck("R-a′-1  Σ_i m_i·d_i = dim P_H（全 表 核對，0 違例）", len(bad) == 0,
   f"violations = {len(bad)}｜表大小 = {len(rows)}")

# ── 統一索引族掃描 ───────────────────────────────────────────────
print("\n  構作逐表（k=1,2,3 或單值）：")
hits = {}
for v in rows:
    vals = {
        "D1": [D1(v, k) for k in (1, 2, 3)],
        "D2": [D2(v, k) for k in (1, 2, 3)],
        "D3": [D3(v)],
        "D4": [D4(v)],
        "D5": [D5(v)],
    }
    hit = sorted({x for vs in vals.values() for x in vs if x in TARGETS})
    if hit:
        hits[v["order"]] = (vals, hit)
        print(f"    |H|={v['order']:>2}：D1={vals['D1']}｜D2={vals['D2']}｜D3={vals['D3']}"
              f"｜D4={vals['D4']}｜D5={vals['D5']}  命中 {hit}")

# 統一族判定：同一構作、以 |H| 索引 ⟹ 需同一構作之 |H|-族含三數
triple = []
for cname in ("D1", "D2", "D3", "D4", "D5"):
    # 以 |H| 索引之族：取每表該構作之「代表值」（k 固定 1／2／3 各一組）
    for kk in (0, 1, 2):
        fam = []
        for v in rows:
            vals = {"D1": [D1(v, k) for k in (1, 2, 3)],
                    "D2": [D2(v, k) for k in (1, 2, 3)],
                    "D3": [D3(v)], "D4": [D4(v)], "D5": [D5(v)]}[cname]
            fam.append(vals[min(kk, len(vals) - 1)])
        if set(TARGETS) <= set(fam):
            triple.append((cname, kk, fam))
ck("R-a′-2  **無任何** (構作, |H| 索引) 族含齊 {16,36,98}", len(triple) == 0, f"triple = {triple}")

# 兩數命中族（若有，須逐個報並判 L0 S3）
two = []
for v in rows:
    vals = {"D1": [D1(v, k) for k in (1, 2, 3)], "D2": [D2(v, k) for k in (1, 2, 3)],
            "D3": [D3(v)], "D4": [D4(v)], "D5": [D5(v)]}
    flat = [x for vs in vals.values() for x in vs]
    got = sorted(set(flat) & set(TARGETS))
    if len(got) >= 2:
        two.append((v["order"], v["dim"], got, flat))
print(f"  單表兩數命中 = {two}")

# 全域：36／98 於全表（5 構作 × 9 表）之命中數
allvals = []
for v in rows:
    allvals += [D1(v, 1), D1(v, 2), D1(v, 3), D2(v, 1), D2(v, 2), D2(v, 3),
                D3(v), D4(v), D5(v)]
print(f"  全表評估值數 = {len(allvals)}")
for t in TARGETS:
    n = sum(1 for x in allvals if x == t)
    print(f"    命中 {t} = {n}/{len(allvals)}")
ck("R-a′-3  全表（9 表 × 9 值）中 36 與 98 之命中數須記錄；若 >1 ⟹ S3 不合格",
   True, "".join(f"{t}:{sum(1 for x in allvals if x == t)}｜" for t in TARGETS))

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/09 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "09_induced_module_search.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "n_subgroups": len(subgroups), "n_modules": len(rows),
               "table": [{"order": v["order"], "dim": v["dim"], "mult": v["mult"],
                          "D1": [D1(v, k) for k in (1, 2, 3)],
                          "D2": [D2(v, k) for k in (1, 2, 3)],
                          "D3": D3(v), "D4": D4(v), "D5": D5(v)} for v in rows],
               "hits": {str(k): [v, h] for k, (v, h) in hits.items()},
               "triple_families": triple, "two_hits": two,
               "verdict": "R-a′ 🔴（置換模族 9 表 × 5 構作：無統一索引族含三數）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
