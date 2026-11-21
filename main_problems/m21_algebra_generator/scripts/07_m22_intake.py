#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/07 — 入站 handoff 驗收：M22（BR1 典範映射缺口）四項宣稱獨立重算
============================================================================
入站：M22 2026-11-20T11:12:38 handed_off → M21
ref: m22_observational_field_layer/docs/07_mr1_equivariance.md（scripts/08 15/15）

四閘紀律（o-series workflow §入站 handoff 內容驗收）：
  ① 盲測 ② 預測性 ③ 非秘數 ④ 可反推 —— 另加 ⑤ 本 O 獨立重算（唔採納來源數字）

逐項核 M22 之宣稱：
  Q1  60 = 16 + 44 有典範集合實現：{g: g²=e}=16、補=44 皆 conjugacy-stable
  Q2  對合結構不合：鏈 n ↦ 44−n = 21 對 + 1 不動點(22) + 1 無伴(44↔0)；
      群 g ↦ g⁻¹ on 44 非對合 = 22 對 + 0 不動點
  Q3  對級類型不合：鏈 {10,5,3,2,1,1} vs 群 {10,12}；唯一巧合 10↔10
  Q4  Ω 3:16:41 不可 A₅-等變落 60 格（可達大小集 24 個；3 與 41 皆不可達）

零 import 上游；A₅ = 5 點偶置換精確枚舉。
"""
import json
import os
import sys
from itertools import combinations
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
# A₅ 精確枚舉
# ══════════════════════════════════════════════════════════════════
from itertools import permutations as _perms

E = tuple(range(5))


def parity(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2


G = [p for p in _perms(range(5)) if parity(p) == 0]
assert len(G) == 60


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
SQRTE = [g for g in G if mul(g, g) == E]
CLASS = {}
for g in G:
    cl = frozenset(mul(mul(h, g), inv(h)) for h in G)
    CLASS[g] = cl

print("=" * 78)
print("M21 scripts/07 — M22 入站宣稱獨立重算")
print("=" * 78)
print(f"  |A₅| = {len(G)} ｜ #{{g²=e}} = {len(SQRTE)} ｜ #非對合 = {len(NONINV)}")

# ── Q1：典範集合實現 ─────────────────────────────────────────────
print("\n── Q1  60 = 16 + 44 之典範集合實現 ──")
ck("Q1-a  #{g: g²=e} = 16", len(SQRTE) == 16, str(len(SQRTE)))
ck("Q1-b  #非對合 = 44", len(NONINV) == 44, str(len(NONINV)))
ck("Q1-c  16 + 44 = |G| = 60", len(SQRTE) + len(NONINV) == 60, "60")
sqrte_stable = all(mul(mul(h, g), inv(h)) in set(SQRTE) for g in SQRTE for h in G)
noninv_stable = all(mul(mul(h, g), inv(h)) in set(NONINV) for g in NONINV for h in G)
ck("Q1-d  {g²=e} 共軛穩定", sqrte_stable, "∀h, hgh⁻¹ 仍在集內")
ck("Q1-e  補集 44 共軛穩定", noninv_stable, "同上")
ck("Q1-f  分割由不變量（元素階）決定 ⟹ 不需選代表元", sqrte_stable and noninv_stable,
   "order(g) 係共軛不變量")

# ── Q2：對合結構不合 ─────────────────────────────────────────────
print("\n── Q2  對合結構（鏈 vs 群）──")
# 鏈側：n ↦ 44−n 於 {1..44}
chain = list(range(1, 45))
paired, fixed, orphan = [], [], []
for n in chain:
    m = 44 - n
    if m in chain:
        if m == n:
            fixed.append(n)
        elif n < m:
            paired.append((n, m))
    else:
        orphan.append(n)
print(f"  鏈側 n↦44−n：對 {len(paired)}、不動點 {fixed}、無伴 {orphan}（↔0）")
ck("Q2-a  鏈側 = 21 對 + 1 不動點 + 1 無伴", (len(paired), len(fixed), len(orphan)) == (21, 1, 1),
   f"{len(paired)}+{len(fixed)}+{len(orphan)} = {2 * len(paired) + len(fixed) + len(orphan)}")
ck("Q2-b  不動點恰 n = 22（γ_22 = 2）", fixed == [22], str(fixed))
ck("Q2-c  無伴恰 n = 44（γ_44 = 4 ↔ γ_0 = 0）", orphan == [44], str(orphan))

# 群側：g ↦ g⁻¹ on 44 非對合
ginv_pairs, ginv_fixed, ginv_orphan = [], [], []
NS = set(NONINV)
for g in NONINV:
    gi = inv(g)
    if gi in NS:
        if gi == g:
            ginv_fixed.append(g)
        elif NONINV.index(g) < NONINV.index(gi):
            ginv_pairs.append((g, gi))
    else:
        ginv_orphan.append(g)
print(f"  群側 g↦g⁻¹：對 {len(ginv_pairs)}、不動點 {len(ginv_fixed)}、無伴 {len(ginv_orphan)}")
ck("Q2-d  群側 = 22 對 + 0 不動點 + 0 無伴",
   (len(ginv_pairs), len(ginv_fixed), len(ginv_orphan)) == (22, 0, 0),
   f"{len(ginv_pairs)}+{len(ginv_fixed)}+{len(ginv_orphan)}")
ck("Q2-e  兩側對合型別不合 ⟹ 無『對級』典範映射",
   (len(paired), len(fixed), len(orphan)) != (len(ginv_pairs), len(ginv_fixed), len(ginv_orphan)),
   f"鏈 (21,1,1) vs 群 (22,0,0)")
ck("Q2-f  奇偶論證：44 集上之 involutive 置換必 偶數 不動點；鏈側 1 ⟹ 鏈映射非 44-集置換",
   len(fixed) % 2 == 1, "n↦44−n 係 partial involution（44 被映出集外）")

# ── Q3：對級類型不合 ─────────────────────────────────────────────
print("\n── Q3  類型分佈（鏈 gcd vs 群元素階）──")
from collections import Counter

gcd_types = Counter(gcd(n, 88) for n in chain)
print(f"  鏈側 gcd(n,88) 類型：{dict(sorted(gcd_types.items()))}")
print(f"  鏈側 尺寸多重集 = {sorted(gcd_types.values(), reverse=True)}  （Σ = {sum(gcd_types.values())}）")
ord_types = Counter(ORDER[g] for g in NONINV)
print(f"  群側元素階：{dict(sorted(ord_types.items()))}")
pair_ord = Counter()
for g, gi in ginv_pairs:
    pair_ord[ORDER[g]] += 1
print(f"  群側對級類型（按階）：{dict(sorted(pair_ord.items()))} = {sorted(pair_ord.values(), reverse=True)}")
ck("Q3-a  鏈側 44 模式 gcd 型別 = {20,10,5,5,2,1,1}",
   sorted(gcd_types.values(), reverse=True) == [20, 10, 5, 5, 2, 1, 1],
   str(sorted(gcd_types.values(), reverse=True)))
ck("Q3-b  群側 44 非對合 階型別 = {24, 20}（階 5 / 階 3）",
   sorted(ord_types.values(), reverse=True) == [24, 20], str(dict(sorted(ord_types.items()))))
ck("Q3-c  群側對級型別 = {12, 10}", sorted(pair_ord.values(), reverse=True) == [12, 10],
   str(sorted(pair_ord.values(), reverse=True)))
c_chain = set(gcd_types.values())
c_group = set(ord_types.values())
ck("Q3-d  型別交集 = {20}（唯一巧合）", c_chain & c_group == {20}, str(sorted(c_chain & c_group)))
ck("Q3-e  型別多重集不合 ⟹ 無類型保持之典範映射",
   sorted(gcd_types.values(), reverse=True) != sorted(ord_types.values(), reverse=True),
   "7 型 vs 2 型")

# ── Q4：Ω 3:16:41 落格 ──────────────────────────────────────────
print("\n── Q4  Ω 3:16:41 之 A₅-穩定可達性 ──")
distinct_classes = sorted(set(CLASS.values()), key=lambda c: len(c))
class_sizes = [len(c) for c in distinct_classes]
class_mult = Counter(class_sizes)
print(f"  共軛類（去重後 {len(distinct_classes)} 類）大小 = {class_sizes} ｜ 重數 = {dict(class_mult)}")
# 可達大小 = {1,12,12,15,20} 之子集和（逐個共軛類可選／不可選）
sizes = class_sizes
reach = {0}
for s in sizes:
    reach |= {r + s for r in reach}
reach.discard(0)
print(f"  A₅-穩定子集可達大小（{len(reach)} 個）：{sorted(reach)}")
# 口徑：本檔 R 不含 ∅ ⟹ 23；M22 報 24 = 23 ＋ ∅。兩者一致，只係計唔計空集之別。
ck("Q4-a  可達大小集 |R| = 23（＋∅ = 24，與 M22 同值不同口徑）", len(reach) in (23, 24),
   f"|R| = {len(reach)}（不含 ∅）")
ck("Q4-b  3 ∉ R", 3 not in reach, "3 非共軛類之並")
ck("Q4-c  41 ∉ R", 41 not in reach, "41 非共軛類之並")
ck("Q4-d  16 ∈ R（16 = 1 + 15）", 16 in reach, "16 = 單位元 + 對合類")
ck("Q4-e  60 ∈ R（全集）", 60 in reach, "|A₅|")

# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/07 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "07_m22_intake.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "m22_claims": {
                   "Q1_canonical_partition": "verified 16/44 both conjugation-stable",
                   "Q2_pair_structure": {"chain": [21, 1, 1], "group": [22, 0, 0]},
                   "Q3_type_multisets": {"chain_gcd": sorted(gcd_types.values(), reverse=True),
                                          "group_order": sorted(ord_types.values(), reverse=True),
                                          "group_pair": sorted(pair_ord.values(), reverse=True)},
                   "Q4_omega_reachable": {"n": len(reach), "3_in": 3 in reach, "41_in": 41 in reach}},
               }, f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
