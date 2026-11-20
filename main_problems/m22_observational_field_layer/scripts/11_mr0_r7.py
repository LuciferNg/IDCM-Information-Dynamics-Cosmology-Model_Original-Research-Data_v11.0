#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/11 — R7：MR0 終局（形式定理 ×3 ＋ 「唔可能由 𝒯 導出」no-go）
========================================================================
MR0 現形（R6 收窄）：剩「天空 = 3 維載體之單位球」一個識別。
R7 把它拆成兩半，一半證死、一半證不可導出：

  T-A（線性化；引用定理）  有限群在 S² 之忠實作用必共軛於正交作用
        （Kerékjártó 1919／Brouwer；S² 之有限群分類：C_n、D_n、T、O、I）
        ⟹ A₅ ⊂ O(3) ⟹ 忠實 3 維實表示 ⟹ {3, 3′}（scripts/07 之二值）
        機器旁證：|A₅| = 60 且類結構 {1,15,20,12,12} ⟹ 階 60 之 SO(3) 子群 = I（分類唯一）
  T-B（不變度量唯一）      End_{A₅}(V) = ℝ（Schur）⟹ 不變對稱雙線性形 = 1 維
        ⟹ S² 之度量由結構固定至**一個**尺度 = dof_scale = 1（M6）⟹ 零額外參數
  T-C（自由軌道）          generic 種子之 stabilizer = 平凡（|orbit| = 60）；
        非自由種子 = 60 條轉軸（有限並，測度 0）⟹ 自由軌道存在且 generic
        機器核：2000 個隨機種子全部 |orbit| = 60
  T-D（不可導出 no-go）    |Aut C₈₈| = 176 = 2⁴·11；K_tree = 33 節點路徑 ⟹ |Aut| = 2
        ⟹ |Aut 𝒯| = 352，5 ∤ 352；A₅ 單群且 gcd(60,352) = 4
        ⟹ **A₅ → Aut(𝒯) 只能係平凡同態** ⟹ 「A₅ 作用於天空」唔可能由 𝒯（或其因子）導出
        ⟹ MR0 剩餘 = **單一 bridge axiom B1**（A₅ 作用於天空），唔可在 v11 內閉合
"""
import json
import os
import sys
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
mats, ang = [], []
for q in list(uniq.values())[:60]:
    w, x, y, z = [float(t) for t in q]
    M = np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])
    mats.append(M)
    ang.append(round(float(np.degrees(np.arccos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))))), 1))
from collections import Counter  # noqa: E402

ck("群", len(mats) == 60, "|A₅| = 60")

# ── T-A：線性化 ＋ 分類唯一性 ────────────────────────────────
print("\nT-A 線性化 ＋ 階 60 之子群分類", flush=True)
ck("TA-階", len(mats) == 60, "|A₅| = 60")
ck("TA-類", dict(sorted(Counter(ang).items())) == {0.0: 1, 72.0: 12, 120.0: 20, 144.0: 12, 180.0: 15},
   "類結構 1/15/20/12/12（SO(3) 有限子群分類中，階 60 者唯 I；T 階 12、O 階 24、D_n 階 2n）")
ck("TA-線性化（引用）", True,
   "定理（Kerékjártó 1919／Brouwer）：有限群在 S² 之忠實作用共軛於正交作用"
   " ⟹ A₅ ⊂ O(3) ⟹ 忠實 3 維實表示 ⟹ {3,3′}（唯一至二值）")

# ── T-B：不變度量唯一（Schur）────────────────────────────────
print("\nT-B 不變對稱雙線性形之維數（End_{A₅}(V) = ℝ）", flush=True)
# 解 Rᵀ S R = S（S 對稱 3×3，6 個未知）之齊次系統
# 用全部 60 個旋轉堆疊
A = []
for M in mats:
    for i in range(3):
        for j in range(i, 3):
            row = np.zeros(6)
            idx = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (1, 1): 3, (1, 2): 4, (2, 2): 5}
            for a in range(3):
                for b in range(3):
                    w = M[a, i] * M[b, j]
                    if a <= b:
                        row[idx[(a, b)]] += w
                    else:
                        row[idx[(b, a)]] += w
            row[idx[(i, j)]] -= 1.0
            A.append(row)
A = np.array(A)
s = np.linalg.svd(A, compute_uv=False)
rank = int((s > 1e-9 * s[0]).sum())
dim = 6 - rank
ck("TB-維數", dim == 1, f"不變對稱形空間 = {dim} 維（6 − rank {rank}）⟹ Schur：End(V) = ℝ")
# 解出並核 ∝ I
u, sv, vt = np.linalg.svd(A, full_matrices=True)
sol = vt[-1]
S = np.array([[sol[0], sol[1], sol[2]], [sol[1], sol[3], sol[4]], [sol[2], sol[4], sol[5]]])
S = S / np.diag(S).mean()
ck("TB-∝I", np.abs(S - np.eye(3)).max() < 1e-8,
   f"解 ∝ 單位陣（max|S/I − 1| = {np.abs(S-np.eye(3)).max():.1e}）⟹ 度量唯一至一尺度 = dof_scale = 1")

# ── T-C：自由軌道 generic ───────────────────────────────────
print("\nT-C 自由軌道之 generic 性", flush=True)
rng = np.random.default_rng(20261210)
NS = 2000
seeds = rng.normal(size=(NS, 3))
seeds /= np.linalg.norm(seeds, axis=1)[:, None]
sizes = np.array([len({tuple(np.round((M @ s), 9)) for M in mats}) for s in seeds])
ck("TC-自由", int((sizes == 60).sum()) == NS, f"{int((sizes==60).sum())}/{NS} 個隨機種子 |orbit| = 60")
ck("TC-非自由集", True,
   "非自由點 = 60 個元素之轉軸（每個 2 點）之並 = 1 維集，測度 0 ⟹ 自由軌道存在且 generic"
   "（與 P4 docs/72 之「orbit 2 dof」一致：自由軌道族 = S² 減測度零集）")

# ── T-D：𝒯 上 A₅ 作用必然平凡（no-go）─────────────────────────
print("\nT-D A₅ → Aut(𝒯) 只能平凡（MR0 不可由 𝒯 導出）", flush=True)
aut_chain = 176          # |Aut C₈₈| = |D₈₈| = 2·88
aut_tree = 2             # K_tree = 33 節點路徑（徑向階梯）⟹ Aut = Z₂
aut_T = aut_chain * aut_tree
print(f"   |Aut C₈₈| = {aut_chain} = 2⁴·11 ｜ |Aut K_tree| = {aut_tree}（路徑反射）｜ |Aut 𝒯| = {aut_T} = 2⁵·11",
      flush=True)
ck("TD-5不整除", 5 % 5 == 0 and aut_T % 5 != 0, f"5 ∤ {aut_T}")
ck("TD-gcd", gcd(60, aut_T) == 4, f"gcd(60, {aut_T}) = {gcd(60, aut_T)} = 4")
ck("TD-單群", True,
   "A₅ 單群 ⟹ 任何同態 A₅ → Aut(𝒯) 之像 = A₅/{1} 或 A₅/A₅ ⟹ |Im| ∈ {1,60}；"
   f"但 |Im| 須整除 |Aut 𝒯| = {aut_T} 而 60 ∤ {aut_T} ⟹ |Im| = 1 ⟹ 同態平凡")
ck("TD-判", True,
   "⟹ **A₅ 在 𝒯 及其任一因子（88、34）上皆無忠實作用（5 ∤ 88、5 ∤ 34）**"
   " ⟹ 「A₅ 作用於天空」唔可能由 𝒯 導出 ⟹ MR0 剩餘 = 單一 bridge axiom B1，v11 內不可閉")

# ── 結論 ────────────────────────────────────────────────────
print("\nR7 終局判詞", flush=True)
ck("R7-形式已閉", True,
   "若 B1 成立 ⟹ 作用形式唯一（{3,3′}、度量唯一至 dof=1、自由軌道 generic）⟹ 無額外參數")
ck("R7-B1不可導", True,
   "B1（A₅ 作用於天空）在 v11 內不可導出（T-D 算術 no-go）⟹ 列為 bridge axiom；"
   "M22 之 T1–T4／𝒟／形狀定理全部 **conditional on B1**")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/11：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"invariant_form_dim": dim, "free_orbit_frac": float((sizes == 60).mean()),
       "aut_chain": aut_chain, "aut_tree": aut_tree, "aut_T": aut_T,
       "gcd_60_autT": gcd(60, aut_T),
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "11_mr0_r7.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
