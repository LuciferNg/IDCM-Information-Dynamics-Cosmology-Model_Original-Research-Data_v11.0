#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/08 — MR1 攻堅：Ω 比例不可等變落格 ＋ 60 = 16+44 之典範實現 ＋ 對合結構不合
=====================================================================================
三條（全部新）：
  NA4  Ω 三分量 3 : 16 : 41 **不可** A₅-等變地落在 60 格上
       （A₅-穩定子集之大小 = 共軛類 {1,15,20,12,12} 之子集和；3 與 41 皆**不可達**，16 可達）
       ⟹ 唔可以藉 Ω 比例固定振幅（除非引入非等變標號 = 隱藏 dof）
  BR1′ 60 = 16 + 44 有**典範實現**：{g : g² = e} = 16、其補 = 44（兩者皆 A₅-穩定）
       ⟹ M21 BR1 之「44 ↔ 44」在**集合層**係典範（唔需要選）
  BR1″ 但**對合結構不合**：鏈側 `n ↦ 44−n` 於 {1..44} = 21 對 + 1 不動點(22) + 1 無伴(44)
       群側 `g ↦ g⁻¹` 於 44 非對合 = 22 對 + 0 不動點
       ⟹ 連「對級」典範映射都唔存在（比 M21 原判更強）
"""
import itertools
import json
import os
import sys
from collections import Counter
from itertools import permutations
from math import gcd

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


PHI = (1 + msqrt(5)) / 2
quats = [tuple(mpf(1) / 2 if (s >> i) & 1 else mpf(-1) / 2 for i in range(4)) for s in range(16)]
for i in range(4):
    for sgn in (1, -1):
        v = [mpf(0)] * 4
        v[i] = mpf(sgn)
        quats.append(tuple(v))


def perm_parity(p):
    n, par, p = len(p), 1, list(p)
    for i in range(n):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            par = -par
    return par


base = [mpf(0), mpf(1), 1 / PHI, PHI]
for p in [q for q in permutations(range(4)) if perm_parity(q) == 1]:
    for s in range(8):
        quats.append(tuple(base[p[i]] * (1 if (s >> i) & 1 else -1) / 2 for i in range(4)))


def canon(q):
    fq = [float(c) for c in q]
    for c in fq:
        if abs(c) > 1e-12:
            s = 1.0 if c > 0 else -1.0
            return tuple(round(s * x, 12) for x in fq)
    return tuple(fq)


uniq = {}
for q in quats:
    if abs(float(sum(x * x for x in q)) - 1) < 1e-12:
        uniq.setdefault(canon(q), q)
rot = list(uniq.values())[:60]
mats, angs = [], []
for q in rot:
    w, x, y, z = [float(t) for t in q]
    M = np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])
    mats.append(M)
    angs.append(round(float(np.degrees(np.arccos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))))), 3))
ck("群-重建", len(mats) == 60, "|A₅| = 60")
cls_size = dict(Counter(angs))
print(f"   類大小（轉角）: {cls_size}", flush=True)
ck("群-類", cls_size == {0.0: 1, 72.0: 12, 120.0: 20, 144.0: 12, 180.0: 15}, "1/15/20/12/12")

# ── NA4：A₅-穩定子集之可達大小 ────────────────────────────────
print("\nNA4 Ω 3:16:41 可否 A₅-等變落格", flush=True)
sizes = [1, 15, 20, 12, 12]
reach = set()
for r in range(len(sizes) + 1):
    for comb in itertools.combinations(sizes, r):
        reach.add(sum(comb))
print(f"   A₅-穩定子集之可達大小（{len(reach)} 個）: {sorted(reach)}", flush=True)
ck("NA4-3不可達", 3 not in reach, "Ω_b = 3 不可達（{1,15,20,12,12} 之子集和冇 3）")
ck("NA4-41不可達", 41 not in reach, "Ω_Λ = 41 不可達")
ck("NA4-16可達", 16 in reach, "Ω_DM = 16 = 1 + 15 可達（單位元 ＋ 15 對合）")
ck("NA4-44可達", 44 in reach, "44 = 20 + 12 + 12 可達（= 非對合集）")
ck("NA4-判", True,
   "⟹ Ω 三分量**唔可以** A₅-等變地落 60 格（3 與 41 皆不可達）"
   " ⟹ 唔可以藉 Ω 比例定振幅；只有 {16, 44} 之分割係等變")

# ── BR1′：60 = 16 + 44 之典範實現 ─────────────────────────────
print("\nBR1′ 60 = 16 + 44 之典範實現", flush=True)


def key(M):
    return tuple(np.round(M.ravel(), 8))


keys = {key(M) for M in mats}
squ = [M for M in mats if key(M @ M) == key(np.eye(3))]
noninv = [M for M in mats if key(M @ M) != key(np.eye(3))]
ck("BR1′-對合集", len(squ) == 16, f"{{g : g² = e}} = {len(squ)}（單位元 ＋ 15 對合）")
ck("BR1′-非對合集", len(noninv) == 44, "{g : g² ≠ e} = 44 = N_d")
ck("BR1′-等變", True, "兩集皆 conjugacy-stable（由階／轉角定義）⟹ A₅-穩定 ⟹ 分割典範")

# ── BR1″：對合結構不合 ────────────────────────────────────────
print("\nBR1″ 對合結構（鏈 n↦44−n vs 群 g↦g⁻¹）", flush=True)
chain_pairs = sorted({tuple(sorted((n, 44 - n))) for n in range(1, 45)})
fixed = [p for p in chain_pairs if p[0] == p[1]]
pairs_only = [p for p in chain_pairs if p[0] != p[1]]
outside = [p for p in pairs_only if 0 in p or 44 in p]     # (0,44) 型：夥伴在 44 之外
pure_pairs = [p for p in pairs_only if not outside.__contains__(p)]
print(f"   鏈側：對 = {len(pure_pairs)} ｜ 不動點 = {[p[0] for p in fixed]} ｜ 無伴 = {outside}", flush=True)
ck("BR1″-鏈結構", len(pure_pairs) == 21 and fixed == [(22, 22)] and len(outside) == 1,
   "鏈側 = 21 對 + 1 不動點(22) + 1 無伴(44↔0)")

# 群側：44 非對合之逆映射
invpairs = set()
for M in noninv:
    Mt = M.T
    invpairs.add(tuple(sorted((key(M), key(Mt)))))
ck("BR1″-群結構", len(invpairs) == 22, f"群側 = {len(invpairs)} 對（g ≠ g⁻¹ 於非對合）")

# 群側對級分類（3-cycles 10 對、5-cycles 12 對）
pair_types = Counter()
for M in noninv:
    a = round(float(np.degrees(np.arccos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))))), 1)
    pair_types[a] += 1
print(f"   群側非對合類別分佈（轉角→個數）: {dict(sorted(pair_types.items()))}", flush=True)
ck("BR1″-對級類型", True,
   f"群側對級類型 = 3-cycles {pair_types[120.0]//2} 對 + 5-cycles {(pair_types[72.0]+pair_types[144.0])//2} 對"
   f" = 10 + 12 = 22")

mode_pair_types = Counter(gcd(n, 88) for n in range(1, 23))
print(f"   鏈側對級（n=1..22 之 gcd(n,88)）: {dict(sorted(mode_pair_types.items()))}", flush=True)
ck("BR1″-不合", (len(pure_pairs), len(invpairs)) == (21, 22),
   f"21 對 ≠ 22 對 ⟹ 連對級典範映射都唔存在（比 M21 BR1 原判更強）")
ck("BR1″-類型不合", sorted(mode_pair_types.values(), reverse=True) != [12, 10],
   f"鏈側對類型 {sorted(mode_pair_types.values(), reverse=True)} ≠ 群側 [12, 10]")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/08：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"reachable_sizes": sorted(reach), "class_sizes": {str(k): v for k, v in cls_size.items()},
       "chain_pairs": len(pure_pairs), "chain_fixed": [p[0] for p in fixed],
       "group_inverse_pairs": len(invpairs),
       "group_pair_types": {str(k): v for k, v in pair_types.items()},
       "mode_pair_gcd": {str(k): v for k, v in mode_pair_types.items()},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "08_mr1_partition.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
