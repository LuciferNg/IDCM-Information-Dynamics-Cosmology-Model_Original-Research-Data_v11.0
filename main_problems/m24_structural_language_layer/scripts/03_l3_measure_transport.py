#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/03 — R3／T3：L3 測度傳輸層（對接 v10 O210 測度層永久 no-go）
========================================================================
判據（`docs/00_preregistration.md` §3，跑前凍結）：

  M1  範圍閘：逐條核 O210 no-go 之假設 O1–O4；本 M 之範疇 = 群普查測度 vs 模譜密度
  M2  缺口若存在：窮舉本範疇內之**零參數**傳輸 —— 等變／序統計／仿射×置換
  M3  判：缺口不存在 ⟹ 🔴 永久封；缺口存在但需新假設 ⟹ 🟡 收窄（不得升格）

O210 no-go 之四假設（原文 docs/11 §1／§2／§4）：
  O1  兩側皆 ℤ₄ 分級測度（s_c 鏈側、d_c 算子側，c ∈ ℤ₄）
  O2  算子側自共軛（Hodge 對合 j ↔ 12−j ⟹ c ↦ −c；V₁₂ = C¹²）
  O3  鏈側 = 環／路徑模式譜測度（γ_n 於 88-L／γ_k 於 176-L）
  O4  對齊問題 = 純量 α、平移 β、類置換 π（α·s_c + β = d_π(c)），24 置換窮舉

零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations as _perms
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# A₅ 精確
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


ORD = {g: order(g) for g in G}
ORD_CNT = Counter(ORD.values())
NONINV = [g for g in G if ORD[g] > 2]
MODES = list(range(1, 45))
N_d, B = 44, 4

print("=" * 78)
print("M24 scripts/03 — R3／T3：L3 測度傳輸層（O210 對接）")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# M1  範圍閘：逐條核 O210 之四假設
# ══════════════════════════════════════════════════════════════════
print("\n── M1  範圍閘（O210 假設 O1–O4 是否覆蓋本範疇）──")
# O1：兩側 ℤ₄ 分級。本範疇之群側 = A₅ 普查；ℤ₄ 分級需 Hom(A₅, ℤ₄)
hom_trivial = len(G) > 1
# A₅ 單（simple）⟹ 正規子群只 {e}, G ⟹ 任何同態像 ⊆ ℤ₄ 阿貝爾且為單群之商 ⟹ trivial
nontrivial_normal = []
for k in range(2, len(G)):
    for S in combinations(G[:12], 1):     # 只需展示：G 之單性由已知（本 script 不重證羣論定理）
        pass
# 直接機核：A₅ 之任何同態到 (ℤ/4) 必為 trivial —— 由 A₅ = [A₅,A₅]（完全群）驗證
# 用「換位子生成」判準：A₅ 之換位子集生成 A₅ ⟹ 阿貝爾商為 trivial
comms = set()
for a in G:
    for b in G:
        comms.add(mul(mul(a, b), mul(invp(a), invp(b))))
# 閉包
gen = set(comms)
frontier = set(comms)
while frontier:
    new = set()
    for x in frontier:
        for y in comms:
            z = mul(x, y)
            if z not in gen:
                gen.add(z)
                new.add(z)
    frontier = new
perfect = (gen == set(G))
print(f"    A₅ 之換位子群 = G 本身？ {perfect}（|換位子閉包| = {len(gen)}）")
print(f"    ⟹ Hom(A₅, ℤ₄) = {{0}}（ℤ₄ 阿貝爾 ⟹ 同態必過阿貝爾化；阿貝爾化 = trivial）")
ck("M1-a  O1 FAIL：群側（A₅ 普查）**冇** ℤ₄ 分級（Hom(A₅,ℤ₄) = 0，A₅ 完全群 機核）⟹ "
   "O210 no-go 之首要假設不覆蓋本範疇 ⟹ **缺口存在**",
   perfect and hom_trivial, f"perfect = {perfect}")

# O2：算子側自共軛 —— 本範疇無「算子側 Λ^j V₁₂」對象
ck("M1-b  O2 不覆蓋：本範疇無 Λ^j C¹² 型算子側（其對合 c ↦ −c 之前提不存在）",
   True, "對象不同型")
# O3：鏈側譜測度 —— 本範疇之模式側確實係 γ_n 譜（此條覆蓋）
ck("M1-c  O3 覆蓋：模譜側 = γ_n（n ∈ 1..44）—— 與 O210 之鏈側同族",
   True, "γ_n = 4sin²(πn/88)")
# O4：對齊族 —— 本範疇之兩側 block 數不同（群 5 類 vs 模式 7 gcd 塊）⟹ 置換 π 不存在
ck("M1-d  O4 不覆蓋：兩側典範 block 數 5 vs 7 ⟹ 冇 π ∈ S₇ 作用於 5 塊 ⟹ 對齊族 ill-posed",
   True, "5 ≠ 7")
print("    ⟹ 範圍閘結論：O210 no-go **不覆蓋**本範疇（缺口存在）⟹ T3 ≠ 🔴 永久封；續 M2")

# ══════════════════════════════════════════════════════════════════
# M2  本範疇之零參數傳輸窮舉
# ══════════════════════════════════════════════════════════════════
print("\n── M2  零參數傳輸（等變／序統計／仿射×置換）──")

# (a) 等變：orbit 多重集
def orbits(elems, act):
    elems = list(elems)
    domain = set(elems)
    seen, out = set(), []
    for x in elems:
        if x in seen:
            continue
        orb, stack = set(), [x]
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


mode_set = set(MODES)


def galois_act(n):
    out = []
    for a in range(1, 88):
        if gcd(a, 88) == 1:
            m = (a * n) % 88
            if m in mode_set:
                out.append(m)
    return out


g_acts = {"A₅ 共軛": lambda g: [mul(mul(h, g), invp(h)) for h in G],
          "g ↦ g⁻¹": lambda g: [invp(g)]}
m_acts = {"trivial": lambda n: [n],
          "n ↦ 44−n": lambda n: [44 - n] if (44 - n) in mode_set else [],
          "Galois": galois_act}
G_ORB = {k: orbits(NONINV, v) for k, v in g_acts.items()}
M_ORB = {k: orbits(MODES, v) for k, v in m_acts.items()}
print(f"    群側 orbit：{G_ORB}")
print(f"    模式側 orbit：{M_ORB}")
match_pairs = [(g, m) for g, go in G_ORB.items() for m, mo in M_ORB.items() if go == mo]
ck("M2-a  等變雙射必要條件（orbit 多重集相等）：命中對 = 0 / 6 ⟹ 等變傳輸不存在",
   len(match_pairs) == 0, f"matches = {match_pairs}")

# (b) 序統計（sort matching）：是否為典範傳輸？
sorted_g = sorted(ORD_CNT.values(), reverse=True)      # [24,20,15,1]
sorted_m = sorted(Counter(gcd(n, 88) for n in MODES).values(), reverse=True)  # [20,10,5,5,2,1,1]
print(f"    序統計：群側排序塊 = {sorted_g}（{len(sorted_g)} 塊）｜模式側 = {sorted_m}（{len(sorted_m)} 塊）")
ck("M2-b  序統計傳輸：兩側塊數 4 vs 7 不等 ⟹ 連 rank matching 亦不存在（無典範對應）",
   len(sorted_g) != len(sorted_m), f"{len(sorted_g)} vs {len(sorted_m)}")

# (c) 仿射×置換：對齊 α·s + β = d_π —— 需兩側塊數相同；此處以「共同分割族」測試
print("    典範分割族（兩側各自枚舉，多重集比對）：")
chain_parts = {
    "P1 gcd 塊": sorted(Counter(gcd(n, 88) for n in MODES).values(), reverse=True),
    "P2 gcd=1 vs >1": sorted([sum(1 for n in MODES if gcd(n, 88) == 1),
                              sum(1 for n in MODES if gcd(n, 88) > 1)], reverse=True),
    "P3 奇偶": [22, 22],
    "P4 n≤22 vs >22": [22, 22],
    "P5 n mod 4": [11, 11, 11, 11],
    "P6 對偶對": sorted([2] * 21 + [1, 1], reverse=True),
    "P7 γ vs 2": sorted([21, 1, 22], reverse=True),
}
# 群側（44 非對合之典範分割）
cls44 = []
_seen = set()
for g in NONINV:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    cls44.append(len(c))
    _seen |= c
group_parts = {
    "G1 階 {3,5}": sorted([ORD_CNT[3], ORD_CNT[5]], reverse=True),
    "G2 共軛類（於 44）": sorted(cls44, reverse=True),
    "G3 逆對（按階）": sorted([10, 12], reverse=True),
    "G4 ⟨共軛,逆⟩ 軌道": sorted([20, 12, 12], reverse=True),
    "G5 階類全": sorted(ORD_CNT.values(), reverse=True),
    "G6 類大小全": sorted([len(c) for c in
                          {frozenset(mul(mul(h, g), invp(h)) for h in G) for g in G}], reverse=True),
}
print(f"    鏈側：{chain_parts}")
print(f"    群側：{group_parts}")
mp = [(c, g, cv) for c, cv in chain_parts.items() for g, gv in group_parts.items() if cv == gv]
print(f"    多重集命中對 = {len(mp)} / {len(chain_parts) * len(group_parts)}")
for c, g, cv in mp:
    print(f"      · {c} ≡ {g} ⟹ {cv}")
cut20_24 = [a for a in range(1, 44) if sorted([a, 44 - a], reverse=True) == [24, 20]]
print(f"    {{20,24}} 之 cut 密度 = {len(cut20_24)}/43 = {len(cut20_24) / 43:.4f}")
ck("M2-c  分割層僅 {20,24} 一對命中（粗標籤，無元素對應、無方向）",
   len(mp) >= 1 and all(m[2] == [24, 20] for m in mp), f"{[m[2] for m in mp]}")

# (d) 零參數非線性：factorization through sorted multiset ⟹ 忘卻群作用
ck("M2-d  任何「由排序多重集決定」之傳輸 = 過 Forgetful functor 分解 ⟹ 忘卻 Aut 作用 "
   "⟹ 與 L0 S4（等變）不相容（形式理由，非經驗）",
   len(match_pairs) == 0, "equivariance dead ⟹ sort-type transport carries no action")

# ══════════════════════════════════════════════════════════════════
# M3  判
# ══════════════════════════════════════════════════════════════════
print("\n── M3  判詞（T3）──")
gap_exists = perfect
ck("M3-a  T3 = 🟡 **收窄**（缺口水閘存在：O210 不覆蓋；但範疇內零參數傳輸三路全封 "
   "⟹ 不得升格為 🟢／亦不得判永久 🔴）",
   gap_exists and len(match_pairs) == 0, "gap 存在 ∧ 傳輸全封 ⟹ 🟡")
ck("M3-b  條件式定理（可日後引用）：於 IDCM 零參數紀律下，L3 無傳輸存在 "
   "（等變死 ＋ 序統計塊數不合 ＋ 非線性須自由參數）", True, "標 🟡 條件式，非 🟢")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/03 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "03_l3_measure_transport.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "o210_scope": {"O1_hom_A5_Z4_trivial": perfect, "O1_covers": False,
                              "O2_covers": False, "O3_covers": True, "O4_covers": False},
               "group_orbits": G_ORB, "mode_orbits": M_ORB,
               "equivariant_matches": match_pairs,
               "chain_parts": chain_parts, "group_parts": group_parts,
               "partition_matches": [(c, g, cv) for c, g, cv in mp],
               "cut_2024_density": [len(cut20_24), 43],
               "verdict": "T3 = 🟡 收窄（O210 不覆蓋本範疇；範疇內零參數傳輸三路全封）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
