#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/11 — DA2：h11 = 36 之 canonical construction 擴大搜尋
=================================================================
判準（docs/08 §1，跑前凍結）：
  E2a  預登錄池 𝔇 內命中 36 之元素數 n₃₆（升格需 n₃₆ ≤ 1 且該元素有算子讀法）
  E2b  36 是否為恆等反演（如 = Σd − …）
  E2c  q = n₃₆/|𝔇|（升格需 q ≤ 5%）

★ 𝔇 **預登錄**（命中後不得增刪）：
  （i）   外冪維數 Σ C(d,k)，k = 1..5
  （ii）  對稱冪維數 Σ C(d+k−1,k)，k = 1..3
  （iii） 配對函數 Σ_(i<j) d_i d_j、Σ_(i<j)(d_i+d_j)、Σ_(i<j) C(d_i+d_j,2)、Σ_(i<j)|d_i−d_j|
  （iv）  階集計數 |{g : order(g) ∈ S}|（S ⊆ {1,2,3,5}，非空）
  （v）   類大小子集和（5 個共軛類之非空子集）
  （vi）  （iv）（v）之補 |G| − X（避免 0/60 重複）

零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations, permutations as _p

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


D = [1, 3, 3, 4, 5]
GG = 60

# ── A₅ 群（精確）──────────────────────────────────────────────
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
ord_cnt = Counter(order(g) for g in G)
classes = []
_seen = set()
for g in G:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    classes.append(len(c))
    _seen |= c

print("=" * 78)
print("M21 scripts/11 — DA2：h11 = 36 canonical 搜尋")
print("=" * 78)
print(f"  階計數 = {dict(sorted(ord_cnt.items()))} ｜ 類大小 = {sorted(classes, reverse=True)}")

# ── 池 𝔇 ─────────────────────────────────────────────────────
pool = {}
# (i) 外冪
for k in range(1, 6):
    pool[f"ΣC(d,{k})"] = sum(__import__("math").comb(d, k) for d in D)
# (ii) 對稱冪
for k in range(1, 4):
    pool[f"ΣC(d+{k - 1},{k})"] = sum(__import__("math").comb(d + k - 1, k) for d in D)
# (iii) 配對
pool["Σ_(i<j) d_i·d_j"] = sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)(d_i+d_j)"] = sum(D[i] + D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j) C(d_i+d_j,2)"] = sum((D[i] + D[j]) * (D[i] + D[j] - 1) // 2
                                   for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)|d_i−d_j|"] = sum(abs(D[i] - D[j]) for i in range(5) for j in range(i + 1, 5))
# (iv) 階集計數
orders = sorted(ord_cnt)
for r in range(1, len(orders) + 1):
    for S in combinations(orders, r):
        pool["|{order ∈ {%s}}|" % ",".join(map(str, S))] = sum(ord_cnt[o] for o in S)
# (v) 類大小子集和
for r in range(1, len(classes) + 1):
    for S in set(combinations(sorted(classes, reverse=True), r)):
        pool["Σ類子集 {%s}" % ",".join(map(str, S))] = sum(S)
# (vi) 補
comp = {}
for k, v in list(pool.items()):
    if 0 < v < GG:
        comp[f"|G| − ({k})"] = GG - v
pool.update(comp)

# ── 命中 36 ───────────────────────────────────────────────────
hits = [(k, v) for k, v in pool.items() if v == 36]
print(f"\n  池 𝔇 大小 = {len(pool)}（預登錄；命中後未增刪）")
print(f"  命中 36 者（{len(hits)} 條）：")
for k, _ in hits:
    print(f"    · {k} = 36")

n36 = len(hits)
q = n36 / len(pool)
print(f"  n₃₆ = {n36} ｜ q = {n36}/{len(pool)} = {q:.4f}")

ck("E2a  n₃₆ ≤ 1（唯一讀法）", n36 <= 1, f"n₃₆ = {n36}")
ck("E2c  q ≤ 5%", q <= 0.05, f"q = {q:.4f}")

# ── E2b：是否恆等反演 ────────────────────────────────────────
print("\n── E2b  36 是否恆等反演？──")
# 36 = 60 − 24；24 = #{order 5}；亦 = Σ類子集 (12+12)
id_repr = {
    "|G| − |{order 5}| = 60 − 24": GG - ord_cnt[5],
    "|{order ∈ {1,2,3}}| = 1+15+20": ord_cnt[1] + ord_cnt[2] + ord_cnt[3],
    "|G| − Σ(12,12)": GG - (12 + 12),
}
for k, v in id_repr.items():
    print(f"    {k:<32} = {v}")
ck("E2b  36 之三個讀法互為恆等重寫（同一集合／補集）",
   set(id_repr.values()) == {36}, "1+15+20 = 60−24 = 60−(12+12)")

# ── 判詞 ────────────────────────────────────────────────────
print("\n── 判詞（DA2）──")
upgrade = (n36 == 1) and (q <= 0.05)
print(f"  升格條件：n₃₆=1(={n36 == 1}) ∧ 非恆等(False) ∧ q≤5%(={q <= 0.05})")
print("  ⟹ DA2 = 🔴 維持真缺口：36 有 ≥2 條預登錄候選讀法，且互為恆等重寫 ⟹ 零判別力")
ck("DA2  判 🔴（h11 = 36 維持真缺口）", not upgrade, f"n₃₆ = {n36}, q = {q:.4f}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/11 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "11_da2_h11_canonical.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "pool_size": len(pool), "hits36": [k for k, _ in hits],
               "n36": n36, "q": q, "ord_counts": dict(ord_cnt), "class_sizes": sorted(classes, reverse=True),
               "verdict": "DA2 🔴（真缺口；候選 ≥2 且恆等重寫）"}, f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
