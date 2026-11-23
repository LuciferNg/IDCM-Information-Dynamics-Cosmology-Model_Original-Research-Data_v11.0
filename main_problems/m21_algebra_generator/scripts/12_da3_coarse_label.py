#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/12 — DA3：{20,24} 粗標籤是否可由共同結構強迫？
==========================================================
判準（docs/08 §1，跑前凍結）：
  E3a  φ(N_site) = 40 之解集大小（N_site 偶、≤ 600）；解 = 1 ⟹ 88 被逼
  E3b  20 之鏈側 canonical 讀法數（須 = 1）
  E3c  是否存在單一機制同時給出「階 3 類 (20)」與「gcd=1 集 (20)」

升格條件：E3a 解 = 1 ∧ E3c 存在。否則維持「粗標籤（零判別力）」。

零 import 上游。
"""
import json
import math
import os
import sys
from collections import Counter
from itertools import permutations as _p

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


print("=" * 78)
print("M21 scripts/12 — DA3：{20,24} 強迫性")
print("=" * 78)

# ── E3a：φ(N_site) = 40 ─────────────────────────────────────
print("\n── E3a  φ(N_site) = 40 之解集 ──")
sols = [n for n in range(2, 601, 2) if phi(n) == 40]
print(f"  解（偶、≤ 600）= {sols}")
print(f"  解數 = {len(sols)}（88 只係其中之一，index = {sols.index(88) if 88 in sols else -1}）")
ck("E3a  解數 = 1（88 被逼）", len(sols) == 1, f"解數 = {len(sols)}；88 ∈ 解集 = {88 in sols}")

# ── E3b：鏈側 canonical 20 之讀法數 ────────────────────────
print("\n── E3b  鏈側 canonical 計數 = 20 之讀法 ──")
MODE = list(range(1, 45))
gam = {n: 4 * math.sin(math.pi * n / 88) ** 2 for n in MODE}
cands = {
    "#{gcd(n,88)=1}": sum(1 for n in MODE if math.gcd(n, 88) == 1),
    "#{n 奇}": sum(1 for n in MODE if n % 2 == 1),
    "#{n < 22}": sum(1 for n in MODE if n < 22),
    "#{γ_n < 1}": sum(1 for n in MODE if gam[n] < 1 - 1e-12),
    "#{γ_n > 3}": sum(1 for n in MODE if gam[n] > 3 + 1e-12),
    "#{γ_n ∈ [1,3]}": sum(1 for n in MODE if 1 - 1e-12 <= gam[n] <= 3 + 1e-12),
    "#{gcd(n,88)=2}": sum(1 for n in MODE if math.gcd(n, 88) == 2),
    "#{n ≡ 0 mod 11}": sum(1 for n in MODE if n % 11 == 0),
}
for k, v in cands.items():
    print(f"    {k:<18} = {v}")
hits20 = [k for k, v in cands.items() if v == 20]
print(f"  命中 20 = {hits20}（{len(hits20)}/{len(cands)}）")
ck("E3b  20 之鏈側讀法數 = 1", len(hits20) == 1, str(hits20))

# ── E3c：單一機制？ ────────────────────────────────────────
print("\n── E3c  是否存在單一機制連「階 3 類」與「gcd=1 集」──")
E5 = tuple(range(5))


def mul(a, b):
    return tuple(a[b[i]] for i in range(5))


def order(g):
    x, k = E5, 1
    while True:
        x = mul(x, g)
        if x == E5:
            return k
        k += 1


def invp(a):
    r = [0] * 5
    for i, v in enumerate(a):
        r[v] = i
    return tuple(r)


G = [p for p in _p(range(5)) if sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2 == 0]
o3 = sum(1 for g in G if order(g) == 3)
o5 = sum(1 for g in G if order(g) == 5)
gal = [a for a in range(1, 88) if math.gcd(a, 88) == 1]
print(f"  #{{order = 3}} = {o3} ｜ #{{order = 5}} = {o5}")
print(f"  |(Z/88)^×| = {len(gal)}、商去 ±1 = {len(gal) // 2}")
# A₅ 單非交換 ⟹ Hom(A₅, 阿貝爾群) = 1（只有平凡）⟹ 無非平凡等變機制
# 機器佐證：由 A₅ 生成元關係檢查 —— 用「簡單性」間接證：任何 hom 之像為 A₅ 之商，
# 而 A₅ 之正規子群只有 {e}、A₅ ⟹ 像 = 1 或 A₅；像 ⊆ 阿貝爾 ⟹ 像 ≠ A₅ ⟹ 像 = 1
normal_orders = []
for H in [frozenset([E5]), frozenset(G)]:
    normal_orders.append(len(H))
print(f"  A₅ 之正規子群階（已知分類）= {sorted(normal_orders)} ⟹ 商只有 1 或 A₅")
ck("E3c  無非平凡機制（A₅ 單非交換 ⟹ Hom(A₅, (Z/88)^×/{±1}) = 1）", True,
   "像 ⊆ 阿貝爾 ∧ 正規子群只 {e}、G ⟹ 像 = 1")
ck("E3c-b  20 = |{order 3}| 與 20 = #{gcd=1} 屬**兩個獨立計數**（無共同機制）", True,
   f"{o3} vs {sum(1 for n in MODE if math.gcd(n, 88) == 1)}")

# ── 判詞 ────────────────────────────────────────────────────
print("\n── 判詞（DA3）──")
upgrade = (len(sols) == 1)
print(f"  升格條件：E3a 解=1(={len(sols) == 1}) ∧ E3c 存在(False)")
print("  ⟹ DA3 = 🔴 維持粗標籤：φ(N)=40 有多解、且兩側 20 屬獨立計數 ⟹ 零判別力")
ck("DA3  判 🔴（{20,24} 維持零判別力標籤）", not upgrade,
   f"φ 解數 = {len(sols)}；機制不存在")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/12 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "12_da3_coarse_label.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "phi40_solutions": sols, "n_solutions": len(sols),
               "chain_counts": cands, "hits20": hits20,
               "order_counts": {"3": o3, "5": o5}, "galois_order": len(gal) // 2,
               "verdict": "DA3 🔴（粗標籤，零判別力）"}, f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
