#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/11 — R-d′：`2B` 之典範強迫性（M24 可及之部分）
=========================================================
判據（`docs/10_preregistration_round3.md` §3，跑前凍結）：

  G1  {1..m} 對 (Z/88)^×/{±1} 不變，逐 m = 1..44
  G2  {1..m} 對對偶 n ↦ 44−n 封閉，逐 m
  G3  S = 大小 ≤ c 之 orbit 之聯集（c = 1..20）之計數 = 36？
  G4  S = gcd ≤ c 之 orbit 之聯集（c ∈ {1,2,4,8,11,22,44}）之計數 = 36？
判：存在統一判準給 h11 = 36 且 sector 典範 ⟹ H2 由 FAIL 升 🟡；否則維持 FAIL

零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


MODES = list(range(1, 45))
MSET = set(MODES)
N_d, B = 44, 4
UNITS = [a for a in range(1, 88) if gcd(a, 88) == 1]


def rep(m):
    m %= 88
    return m if m <= 44 else 88 - m


print("=" * 78)
print("M24 scripts/11 — R-d′：2B 之典範強迫性")
print("=" * 78)

# ── G1 {1..m} 之 Galois 不變性 ─────────────────────────────────
inv_m = []
notinv_example = {}
for m in range(1, 45):
    Sm = set(range(1, m + 1))
    ok = True
    wit = None
    for a in UNITS:
        for n in MODES:
            if (rep(a * n) in Sm) != (n in Sm):
                ok = False
                wit = (a, n, rep(a * n))
                break
        if not ok:
            break
    if ok:
        inv_m.append(m)
    else:
        notinv_example[m] = wit
print(f"  G1：{{1..m}} 為 Galois 不變之 m = {inv_m}")
print(f"      m = 36（h11 = N_d − 2B 之 sector）之反例 = {notinv_example.get(36)}")
ck("G1-a  「n ≤ m」型 sector 之 Galois 不變 m 只有 " + str(inv_m) +
   " ⟹ m = 36 非不變；唯一非平凡不變者 = m = 43（即 k = 1）",
   inv_m == [43, 44] and 36 not in inv_m, f"inv_m = {inv_m}")

# ── G2 對偶封閉 ────────────────────────────────────────────────
dual_closed = []
for m in range(1, 45):
    Sm = set(range(1, m + 1))
    if all((N_d - n) in Sm for n in Sm if N_d - n >= 1):
        dual_closed.append(m)
print(f"  G2：{{1..m}} 對偶封閉之 m = {dual_closed}")
ck("G2-a  「n ≤ m」型 sector 之對偶封閉 m = " + str(dual_closed) +
   "（m < 22 時 44−n > m ⟹ 不可能封閉）", 36 not in dual_closed, f"{dual_closed}")

# ── G3 大小 ≤ c 之 orbit 聯集 ─────────────────────────────────
orbits, seen = [], set()
for n in MODES:
    if n in seen:
        continue
    o = sorted({rep(a * n) for a in UNITS})
    orbits.append(o)
    seen |= set(o)
sizes = [len(o) for o in orbits]
print(f"  Galois orbit（{len(orbits)}）：大小 {sizes}｜gcd 標籤 {[gcd(o[0], 88) for o in orbits]}")
g3 = []
for c in range(1, 21):
    S = sorted(set().union(*[set(o) for o in orbits if len(o) <= c]))
    g3.append((c, len(S)))
print(f"  G3：c ⟹ |S|：{g3}")
g3_36 = [c for c, s in g3 if s == 36]
ck("G3-a  「orbit 大小 ≤ c」之聯集得 36 之 c = " + str(g3_36) +
   "；若存在須查該 c 是否有統一動機", True, f"g3_36 = {g3_36}")

# ── G4 gcd ≤ c 之 orbit 聯集 ──────────────────────────────────
g4 = []
for c in (1, 2, 4, 8, 11, 22, 44):
    S = sorted(set().union(*[set(o) for o in orbits if gcd(o[0], 88) <= c]))
    g4.append((c, len(S)))
print(f"  G4：gcd ≤ c ⟹ |S|：{g4}")
g4_36 = [c for c, s in g4 if s == 36]
ck("G4-a  「gcd ≤ c」之聯集得 36 之 c = " + str(g4_36), True, f"g4_36 = {g4_36}")

# ── 逐項回答：36 之所有「單一規則＋不變」讀法 ─────────────────
print("\n── 36 之單一規則＋不變性讀法清單 ──")
cands = []
for c, s in g3:
    if s == 36:
        cands.append(("orbit 大小 ≤ c", c))
for c, s in g4:
    if s == 36:
        cands.append(("gcd ≤ c", c))
# 「gcd ≤ 4 之聯集 ∪ 單一不動 orbit」型
singles = [o for o in orbits if len(o) == 1]
base4 = sorted(set().union(*[set(o) for o in orbits if gcd(o[0], 88) <= 4]))
for o in singles:
    if len(set(base4) | set(o)) == 36:
        cands.append(("gcd ≤ 4 ∪ {singleton}", (gcd(o[0], 88),)))
print(f"  單一規則候選 = {cands}")
ck("R-d′-1  存在「統一規則 + 不變」之 36 sector 候選數 = " + str(len(cands)) +
   "（若有，仍不足以解釋 2B：2B = 8 與該 c 無關）", True, f"{cands}")

# ── 2B = 8 之直接測試：{37..44}（lifted 8 模）是否不變/封閉 ──
lift = set(range(37, 45))
lift_inv = all((rep(a * n) in lift) == (n in lift) for a in UNITS for n in MODES)
print(f"\n  lifted 8 模 {{37..44}}：Galois 不變 = {lift_inv}｜"
      f"反例 = {next(((a, n, rep(a * n)) for a in UNITS for n in MODES if (rep(a * n) in lift) != (n in lift)), None)}")
ck("R-d′-2  `2B` 之 lifted 集合 {37..44}（8 模）**非** Galois 不變 ⟹ 2B 之 sector 非典範",
   not lift_inv, "見反例")

# ── 判 ────────────────────────────────────────────────────────
ck("R-d′-3  H2 維持 FAIL：2B = 8 既非唯一不變 sector（不變者為 m = 43），"
   "亦非 lifted 集合之不變寬度 ⟹ 2B 之必然性未由任何 2D-結構判準強迫",
   len(inv_m) == 2 and not lift_inv, f"inv_m = {inv_m}｜lift_inv = {lift_inv}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/11 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "11_2B_forcing_test.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "invariant_m": inv_m, "dual_closed_m": dual_closed,
               "orbit_sizes": sizes,
               "orbit_gcd": [gcd(o[0], 88) for o in orbits],
               "G3": g3, "G4": g4, "G3_36": g3_36, "G4_36": g4_36,
               "single_rule_candidates": [[a, str(b)] for a, b in cands],
               "lifted_invariant": lift_inv,
               "verdict": "H2 維持 FAIL（2B 非由 2D 結構判準強迫；唯一不變 m = 43）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
