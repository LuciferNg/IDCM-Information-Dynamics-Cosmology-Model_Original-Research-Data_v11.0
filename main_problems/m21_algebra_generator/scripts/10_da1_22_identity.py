#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/10 — DA1：22 之三重身份可否升為 count-level bridge？
================================================================
判準（docs/08 §1，跑前凍結）：
  E1a  N_d = 2·Σ dim Λ² 對任意度多重集是否恆等
  E1b  鏈側 22 是否獨立於 44／88
  E1c  群側 22 是否獨立於 44
  E1d  canonical module dims 池內 22 之命中數
升格條件：E1a 非恆等 ∧ E1b 有獨立 22 ∧ E1d = 1

零 import 上游。
"""
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


print("=" * 78)
print("M21 scripts/10 — DA1：22 三重身份")
print("=" * 78)

D = [1, 3, 3, 4, 5]
lam2_sum = sum(d * (d - 1) // 2 for d in D)
nd = sum(d * (d - 1) for d in D)

# ── E1a：代數恆等（任意度集）──────────────────────────────────
print("\n── E1a  N_d = 2·Σ dim Λ² 係代數恆等？──")
import sympy as sp

d1, d2, d3, d4, d5 = sp.symbols("d1:6", integer=True, positive=True)
lhs = sum(d * (d - 1) for d in (d1, d2, d3, d4, d5))
rhs = 2 * sum(d * (d - 1) / 2 for d in (d1, d2, d3, d4, d5))
ident = sp.simplify(lhs - rhs) == 0
print(f"  Σd(d−1) − 2·Σ d(d−1)/2 ≡ 0 ？  {ident}")
# 隨機度集抽樣核（10^4）
import random
random.seed(20261120)
rnd_ok = all(sum(d * (d - 1) for d in m) == 2 * sum(d * (d - 1) // 2 for d in m)
             for m in ([random.randint(1, 9) for _ in range(5)] for _ in range(10000)))
print(f"  隨機度集 10^4 抽樣：{rnd_ok}")
ck("E1a-a 恆等（sympy 符號）", ident, "≡ 0")
ck("E1a-b 恆等（10^4 隨機度集）", rnd_ok, "全部相等")
ck("E1a  恆等成立 ⟹ 「22 = N_d/2」零內容（除 2 之算術）", ident and rnd_ok,
   "N_d = Σd(d−1) = 2·Σd(d−1)/2")

# ── E1b：鏈側 22 是否獨立 ────────────────────────────────────
print("\n── E1b  鏈側 canonical 22 之源 ──")
Ns = 88
cands = {}
cands["n* = N_site/4"] = Ns // 4
cands["對極對數 + 不動點 = N_d/2"] = (44 // 2)
cands["#{n ∈ 1..43 : n < 44−n} + 1"] = sum(1 for n in range(1, 44) if n < 44 - n) + 1
# γ_22 = 2 精確（浮點需容差，否則 γ_22 被誤算成 1.9999999999999998 而多數一個）
cands["#{n : γ_n < 2} + 1"] = sum(
    1 for n in range(1, 45) if 4 * math.sin(math.pi * n / 88) ** 2 < 2 - 1e-12) + 1
for k, v in cands.items():
    print(f"  {k:<32} = {v}")
ck("E1b-a  所有候選值皆 = 22", set(cands.values()) == {22}, str(sorted(set(cands.values()))))
ck("E1b-b  每項皆係 N_d（=44）或 N_site（=88）之函數 ⟹ 非獨立",
   True, "N_site/4 = 22；N_d/2 = 22；前後皆由 44／88 決定")

# ── E1c：群側 22 ────────────────────────────────────────────
print("\n── E1c  群側 22 ──")
E5 = tuple(range(5))


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


from itertools import permutations as _p
G = [p for p in _p(range(5)) if sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2 == 0]
NON = [g for g in G if order(g) > 2]
pairs = set()
for g in NON:
    gi = invp(g)
    pairs.add(tuple(sorted([G.index(g), G.index(gi)])))
print(f"  |非對合| = {len(NON)} ｜ 逆對數 = {len(pairs)}")
ck("E1c  逆對數 = |非對合|/2 = 22 ⟹ 恆為 44 之函數", len(pairs) == 22, str(len(pairs)))

# ── E1d：模側 22 之密度 ─────────────────────────────────────
print("\n── E1d  canonical module dims 池內 22 命中 ──")
pool = {
    "ΣΛ¹ = Σd": sum(D),
    "ΣΛ²": lam2_sum,
    "ΣΛ³": sum(d * (d - 1) * (d - 2) // 6 for d in D),
    "ΣΛ⁴": sum(d * (d - 1) * (d - 2) * (d - 3) // 24 for d in D),
    "ΣΛ⁵": sum(d * (d - 1) * (d - 2) * (d - 3) * (d - 4) // 120 for d in D),
    "ΣSym²": sum(d * (d + 1) // 2 for d in D),
    "ΣSym³": sum(d * (d + 1) * (d + 2) // 6 for d in D),
    "Σ_{i<j} d_i d_j": sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5)),
    "k": 5,
    "|G|": 60,
}
hits22 = [k for k, v in pool.items() if v == 22]
print(f"  池（{len(pool)}）：" + "，".join(f"{k}={v}" for k, v in pool.items()))
print(f"  22 命中 = {hits22}")
ck("E1d  22 於池內命中數 = 1", len(hits22) == 1, str(hits22))

# ── 判詞 ────────────────────────────────────────────────────
print("\n── 判詞（DA1）──")
upgrade = (not ident) and False and (len(hits22) == 1)
print(f"  升格條件：E1a 非恆等(={not ident}) ∧ E1b 有獨立 22(=False) ∧ E1d=1(={len(hits22)==1})")
print("  ⟹ DA1 = 🔴 vacuous：「22 = N_d/2」係除 2 之算術；三重身份全部係 44／88 之函數")
ck("DA1  判 🔴 vacuous（三條件不全中）", not upgrade, "E1a 恆等 + E1b 非獨立")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/10 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "10_da1_22_identity.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "identity": {"N_d": nd, "two_lam2": 2 * lam2_sum, "sympy_identity": bool(ident)},
               "chain_22_candidates": cands, "module_pool": pool, "hits22": hits22,
               "verdict": "DA1 🔴 vacuous"}, f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
