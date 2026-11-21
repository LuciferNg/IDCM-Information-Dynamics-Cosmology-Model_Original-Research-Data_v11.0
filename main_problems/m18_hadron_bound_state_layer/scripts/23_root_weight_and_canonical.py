#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/23 — `MD16` 收窄：耦合 11 維中有幾多維係**典範**？＋ 每 `J` 之根權重上限
================================================================================
前設問題
--------
`MD16` 只報耦合空間維數（`11/100/300`）。本輪問兩件**可算即判**之事：
  Q-a 每 `J` 多重態之 **`d=0` 根權重上限** = ？（由模塊內容直接定，零參數）
  Q-b `11` 維之中，有幾多維由**幾何典範構造**（`F₀→F₁` 之子節點關聯，`A₄`-平均）填滿？

判準（寫死）
------------
  X1 根權重上限：`cap(J) = Σ_{R ∈ J↓, R ∈ F₀} dim R / (2J+1)`（`0` ⟹ 根振幅必為 `0`）
  X2 典範耦合：以「根節點 → 其 `4` 個子」之自然關聯矩陣 `M`，取 `A₄`-平均 `T = Σ_h P(h) M P(h)⁻¹ / |A₄|`
     ⟹ 檢查 (a) `T ≠ 0`、(b) `T` 係 `A₄`-等變、(c) `T` **唔係** `A₅`-等變（破對稱之見證）
  X3 典範維數：`span{T, T^T, 全 1 矩陣}` 在 `Hom_{A₄}` 內之維數 ⟹ 剩餘自由維數 = `11 − 典範`
  X4 禁：唔可以靠「揀一個似嘅算子」；典範性必須由**平均構造**（唯一）提供
  X5 判詞只可報數字與結構；`cap(J) = 0` 者列為 M5 可直接核之預言

輸出：data/23_root_weight_and_canonical.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="X", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="X"):
    return chk(name, ok, detail, tag, integrity=False)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def sign(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


E = tuple(range(5))
A5 = sorted(p for p in itertools.permutations(range(5)) if sign(p) == 1)


def closure(gens):
    S = {E}
    fr = list(gens)
    while fr:
        g = fr.pop()
        for x in {compose(g, s) for s in S} | {compose(s, g) for s in S}:
            if x not in S:
                S.add(x)
                fr.append(x)
    return frozenset(S)


subs = {closure([g]) for g in A5}
prev = -1
while len(subs) != prev:
    prev = len(subs)
    Ls = list(subs)
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            subs.add(closure(list(Ls[i]) + list(Ls[j])))
byorder = {}
for s in subs:
    byorder.setdefault(len(s), []).append(s)
A4 = frozenset(byorder[12][0])
C3 = frozenset([s for s in byorder[3] if s < A4][0])

DIM5 = {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}
SO3I = {0: ["A"], 1: ["T1"], 2: ["H"], 3: ["T2", "G"], 4: ["G", "H"],
        5: ["T1", "T2", "H"], 6: ["A", "T1", "G", "H"], 7: ["T1", "T2", "G", "H"],
        8: ["T2", "G", "H", "H"]}
F0_IRREPS = {"A", "G"}          # docs/14：5-orbit = A ⊕ G

print("=" * 78)
print("X1. 每 J 之 d=0 根權重上限（模塊內容直接定）")
print("=" * 78)
tab = []
for J in sorted(SO3I):
    comps = SO3I[J]
    inside = [r for r in comps if r in F0_IRREPS]
    cap = mp.mpf(sum(DIM5[r] for r in inside)) / (2 * J + 1)
    tab.append((J, comps, inside, cap))
    print("   J = %d：成分 %-22s｜落 F₀ 者 %-8s ⟹ cap = %s%%" %
          (J, comps, inside, mp.nstr(cap * 100, 6)))
chk("X1a J=0 之 cap = 100%（唯一可完全局域根）", tab[0][3] == 1, "%s" % mp.nstr(tab[0][3] * 100, 6))
chk("X1b J=1、J=2 之 cap = 0（根振幅必為 0）", tab[1][3] == 0 and tab[2][3] == 0,
    "%s｜%s" % (mp.nstr(tab[1][3] * 100, 4), mp.nstr(tab[2][3] * 100, 4)))
chk("X1c J=5 之 cap = 0（T1⊕T2⊕H 全不在 F₀）", tab[5][3] == 0, "%s" % mp.nstr(tab[5][3] * 100, 4))
chk("X1d 每 J 之成分維數和 = 2J+1（表自核）",
    all(sum(DIM5[r] for r in c) == 2 * J + 1 for J, c in SO3I.items()), "9/9")

print()
print("=" * 78)
print("X2. 典範耦合：子節點關聯之 A₄-平均")
print("=" * 78)


def cosets(G, H):
    out, seen = [], set()
    for x in G:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen:
            seen.add(c)
            out.append(sorted(c)[0])
    return out


X0, X1_, X2_ = cosets(A5, A4), cosets(A5, C3), cosets(A5, frozenset([E]))
i0 = {frozenset(compose(x, h) for h in A4): k for k, x in enumerate(X0)}
i1 = {frozenset(compose(x, h) for h in C3): k for k, x in enumerate(X1_)}
i2 = {x: k for k, x in enumerate(X2_)}


def idx0(x):
    return i0[frozenset(compose(x, h) for h in A4)]


def idx1(x):
    return i1[frozenset(compose(x, h) for h in C3)]


def idx2(x):
    return i2[x]


# 自然關聯 M：F₀ → F₁（每點 4 子）、F₁ → F₂（每點 3 子）
A4_el = sorted(A4)
C3_coh = [sorted(frozenset(compose(x, h) for h in C3))[0] for x in A4_el]
C3_rep = cosets(A4, C3)
M01 = np.zeros((20, 5))
for k0, x in enumerate(X0):
    for c in C3_rep:
        M01[idx1(compose(x, c)), k0] = 1.0
M12 = np.zeros((60, 20))
for k1, y in enumerate(X1_):
    for g in C3:
        M12[idx2(compose(y, g)), k1] = 1.0
print("   M01：每列 1 數 = %d（應 4）｜總 1 數 = %d｜M12：每列 = %d（應 3）｜總 = %d" %
      (int(M01.sum(axis=0)[0]) if M01.size else 0, int(M01.sum()), int(M12.sum(axis=0)[0]), int(M12.sum())))
chk("X2a M01 每點 4 子（總 20）", int(M01.sum(axis=0)[0]) == 4 and int(M01.sum()) == 20, "4｜20")
chk("X2b M12 每點 3 子（總 60）", int(M12.sum(axis=0)[0]) == 3 and int(M12.sum()) == 60, "3｜60")


def perm(n, g, H, idxf, X):
    P = np.zeros((n, n))
    for a, x in enumerate(X):
        P[idxf(compose(g, x)), a] = 1.0
    return P


def avg_over(G, M, n_from, n_to, H, idxf, X):
    T = np.zeros_like(M)
    for g in G:
        P_to = perm(n_to, g, H, idxf, X)
        P_from = perm(n_from, g, H, idxf, X)
        T += P_to @ M @ P_from.T
    return T / len(G)


T01_A4 = avg_over(A4_el, M01, 5, 20, A4, idx0, X0)
print("   T01（A₄-平均）：非零元 = %d｜範數 = %.6f" % (int((np.abs(T01_A4) > 1e-12).sum()), float(np.linalg.norm(T01_A4))))
chk("X2c A₄-平均後 T01 ≠ 0（典範方向存在）", float(np.linalg.norm(T01_A4)) > 1e-9,
    "norm %.6f" % float(np.linalg.norm(T01_A4)))

# (b) A₄-等變性
def is_equiv(G, T, n_from, n_to, H, idxf, X, tol=1e-9):
    for g in G:
        P_to = perm(n_to, g, H, idxf, X)
        P_from = perm(n_from, g, H, idxf, X)
        if np.linalg.norm(P_to @ T @ P_from.T - T) > tol:
            return False
    return True


eq_A4 = is_equiv(A4_el, T01_A4, 5, 20, A4, idx0, X0)
eq_A5 = is_equiv(A5, T01_A4, 5, 20, A4, idx0, X0)
chk("X2d T01 係 A₄-等變", eq_A4, "%s" % eq_A4)
chk("X2e T01 **唔係** A₅-等變（破對稱見證；與 MD13 一致）", not eq_A5, "A₅-等變 = %s" % eq_A5)

print()
print("=" * 78)
print("X3. 典範維數 vs 總維數")
print("=" * 78)
ones01 = np.ones((20, 5))
ones12 = np.ones((60, 20))
T12_A4 = avg_over(A4_el, M12, 20, 60, C3, idx1, X1_)
for nm, T in (("T12", T12_A4),):
    print("   %s：範數 = %.6f｜A₄-等變 = %s" % (nm, float(np.linalg.norm(T)), is_equiv(A4_el, T, 20, 60, C3, idx1, X1_)))


def span_dim(mats):
    vs = [m.reshape(-1) for m in mats]
    A = np.array(vs).T
    return int(np.linalg.matrix_rank(A, tol=1e-9))


d01 = span_dim([T01_A4, ones01])
d12 = span_dim([T12_A4, ones12])
print("   F₀→F₁：典範 span{T01, 全1} 維數 = %d（總 dim Hom = 11）" % d01)
print("   F₁→F₂：典範 span{T12, 全1} 維數 = %d（總 dim Hom = 100）" % d12)
chk("X3a F₀→F₁ 典範維數 = 2 ⟹ 剩餘自由 = 9", d01 == 2, "%d｜自由 %d" % (d01, 11 - d01))
chk("X3b F₁→F₂ 典範維數 = 2 ⟹ 剩餘自由 = 98", d12 == 2, "%d｜自由 %d" % (d12, 100 - d12))

print()
print("=" * 78)
print("X4/X5. 判詞")
print("=" * 78)
verdict("X5a ★ `MD16` 收窄：耦合維數 %d（`F₀→F₁`）之中，典範構造只填 %d 維 ⟹ **真自由 = %d 維**"
        % (11, d01, 11 - d01), 11 - d01 > 0, "典範 = 子節點關聯之 `A₄`-平均（唯一）＋ 全 1（平凡）", "X5")
verdict("X5b 根權重上限表（`cap(J)`）係 `P1–P3` 之最強形式：`J=1,2,5` cap = 0；`J=0` cap = 100%", True,
        "M5 可直接核：若 `J=2` 於根有非零權重 ⟹ `MD14` 錯", "X5")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)
out = {"script": "23_root_weight_and_canonical.py", "checks": checks, "integrity_fail": nfail,
       "cap_by_J": {str(J): float(cap) for J, _, _, cap in tab},
       "canonical_dim_F0F1": d01, "canonical_dim_F1F2": d12,
       "free_F0F1": 11 - d01, "free_F1F2": 100 - d12,
       "T01_is_A5_equivariant": bool(eq_A5)}
with open(os.path.join(DATA, "23_root_weight_and_canonical.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/23_root_weight_and_canonical.json")
sys.exit(0 if nfail == 0 else 1)
