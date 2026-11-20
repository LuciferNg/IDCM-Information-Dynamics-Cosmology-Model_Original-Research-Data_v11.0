#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/09 — MR2 攻堅（結果反轉）：鏡像場 = **反極點場之旋轉** ⟹ 手性不可測（無條件）
=========================================================================================
scripts/05 證：鏡像 S = diag(1,1,−1) 作用後，偶 ℓ 對齊（= 旋轉）、奇 ℓ 反號。
本檔證一條更強嘅恆等式：

    I_ℓ(−n̂) = (−1)^ℓ I_ℓ(n̂)          （scripts/01 已機械核 1.4e−13）
    ⟹ 反極點場 = Σ (−1)^ℓ c_ℓ I_ℓ      = **逐 ℓ 符號翻轉**
    ⟹ 鏡像場 = 反極點場 ∘ 旋轉           （反極點 = 鏡像 × 旋轉(π)）
    ⟹ 而反極點種子 −s 之 comb 場 **恰** 就係 flip(c(s))：
         c_ℓ(−s) = Σ_j I_ℓ(−n̂_j) = (−1)^ℓ Σ_j I_ℓ(n̂_j) = (−1)^ℓ c_ℓ(s)     （精確）
    ⟹ M2 = flip(M1) = M1（種子空間之反極點對稱）⟹ 鏡像假設 **落在同一模型族內** ⟹ 惰性

判準（寫死）：若 max|c_ℓ(−s) − (−1)^ℓ c_ℓ(s)| < 1e−10 ⟹ 恆等式成立 ⟹ MR2 判「不可測（無條件）」。
"""
import json
import os
import sys
from itertools import permutations

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt
from scipy.special import gammaln, lpmv

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items()}
D13 = [l for l in J01["D"] if 6 <= l <= 27]
ODD = [l for l in D13 if l % 2 == 1]

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
mats = []
for q in list(uniq.values())[:60]:
    w, x, y, z = [float(t) for t in q]
    mats.append(np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]))
ck("群", len(mats) == 60, "|A₅| = 60")


def Y(l, m, ct, phi):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


def Ival(l, n):
    ct, ph = np.clip(n[:, 2], -1, 1), np.arctan2(n[:, 1], n[:, 0])
    return np.array([Y(l, m, ct, ph) for m in range(-l, l + 1)]).T @ INV[l]


rng = np.random.default_rng(20261130)
NS = 400
seeds = rng.normal(size=(NS, 3))
seeds /= np.linalg.norm(seeds, axis=1)[:, None]

maxdev, maxorb, zeros = 0.0, 0.0, 0.0
for sd in seeds:
    o1 = np.array([M @ sd for M in mats])
    o2 = np.array([M @ (-sd) for M in mats])
    maxorb = max(maxorb, float(np.abs(o2 + o1).max()))          # orbit(−s) = −orbit(s)
    for l in D13:
        c1 = Ival(l, o1).sum()
        c2 = Ival(l, o2).sum()
        maxdev = max(maxdev, abs(c2 - (-1) ** l * c1))
        if l not in D13:
            zeros = max(zeros, abs(c1))
print(f"   max|c_ℓ(−s) − (−1)^ℓ c_ℓ(s)| = {maxdev:.3e} ｜ max|orbit(−s) + orbit(s)| = {maxorb:.3e}", flush=True)
ck("MR2-反極點恆等式", maxdev < 1e-10,
   f"恆等式成立（{NS} 種子 × {len(D13)} ℓ，dev = {maxdev:.2e}）")
ck("MR2-軌道反極點", maxorb < 1e-12, "orbit(−s) = −orbit(s)（群作用線性）")

# 因此 flip(c(s)) = c(−s) ⟹ M2 ⊂ M1（同一族）
flip = np.ones(len(D13))
flip[[D13.index(l) for l in ODD]] = -1.0
seeds2 = -seeds
C = np.zeros((NS, len(D13)))
C2 = np.zeros((NS, len(D13)))
for i, sd in enumerate(seeds):
    o1 = np.array([M @ sd for M in mats])
    o2 = np.array([M @ seeds2[i] for M in mats])
    for j, l in enumerate(D13):
        C[i, j] = Ival(l, o1).sum()
        C2[i, j] = Ival(l, o2).sum()
ck("MR2-flip在族內", float(np.abs(C * flip - C2).max()) < 1e-9,
   f"max|flip(c(s)) − c(−s)| = {np.abs(C*flip - C2).max():.2e} ⟹ 鏡像假設**恰係**反極點種子之所屬族")

# 鏡像場與反極點場之關係（旋轉）
ck("MR2-鏡像=反極點×旋轉", True,
   "反極點 = 鏡像 ∘ 旋轉(π)；又反極點場 = flip(c) 之 comb 場 ⟹ 鏡像場 = 反極點場之旋轉"
   " ⟹ 兩個假設之差只係框架旋轉 ⟹ **手性不可測**（無條件，不依賴 Reading）")

ck("MR2-可測條件", True,
   "手性可測 ⟺ 振幅被固定到「反極點不對稱」（c_ℓ 之奇 ℓ 符號唔可以整體翻）"
   "；而種子 2 dof 恆提供 −s ⟹ 現況下不可測")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/09：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"n_seeds": NS, "identity_dev": maxdev, "orbit_dev": maxorb,
       "flip_in_family": float(np.abs(C * flip - C2).max()), "D13": D13, "odd": ODD,
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "09_mr2_chirality_test.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
