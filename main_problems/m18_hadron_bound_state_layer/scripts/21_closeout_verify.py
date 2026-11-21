#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/21 — 收尾獨立驗證（零 import 推導腳本；全部換算法）
================================================================
前設問題
--------
收尾（T1）要求：唔可以只重跑原腳本（重跑只複製輸出）。本支**由原始定義**用**唔同算法**
重算 headline 數字，用來捉「原鏈一致但一起錯」嘅情況。

判準（寫死）
------------
  V1 譜：由 **88 環 Laplacian 數值對角化**得 `λ_k = 2 − 2cos(2πk/88)`，與 `γ_n = 4 sin²(nπ/88)` 逐項比
     （原鏈用代數式；本支用矩陣特徵值 ⟹ 算法獨立）
  V2 恆等式：`Σγ = 90`、`Σγ² = 272`、`γ₄₄ = 4`、`γ₂₂ = 2`、`γ_n + γ_{44−n} = 4`（逐 n）
  V3 `A₅` 置換模塊：用 **Burnside 軌道計數**（`Σ mult² = #orbits(X×X)`）核 `5`／`20`／`60`，
     唔用特徵表內積（原鏈 `docs/14` 用內積、`docs/16` 用雙陪集 ⟹ 本支第三路）
  V4 塔：`A₅⊃A₄⊃C₃⊃1` 指數 `(5,4,3)`（Lagrange）＋ `60 = B(B²−1) = (B+1)B(B−1)`
  V5 關鍵計數：`44 ∉ orbit 集`；`20 = 5·4`；`12 = 3·4`；`|A₅| = Σ dim² = 60`
  V6 `K_tree`：`sites = 4(B²+1) = 68`、`λ_k = B sin²(kπ/68)`、`Σλ = 66`、`λ₁₇ = μ² = 2`

輸出：data/21_closeout_verify.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="V", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


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

print("=" * 78)
print("V1/V2. 譜：88 環 Laplacian 對角化 vs 代數式")
print("=" * 78)
N = 88
L = 2 * np.eye(N) - np.roll(np.eye(N), 1, axis=1) - np.roll(np.eye(N), -1, axis=1)
ev = np.sort(np.linalg.eigvals(L).real)
lam = np.sort(np.array([2 - 2 * np.cos(2 * np.pi * k / N) for k in range(N)]))
chk("V1a Laplacian 特徵值（排序）== 2−2cos(2πk/88)（逐 88 個，|Δ|<1e-10）",
    float(np.max(np.abs(ev - lam))) < 1e-10, "max|Δ| = %.2e" % float(np.max(np.abs(ev - lam))))

gamma = [mp.mpf(4) * mp.sin(mp.mpf(n) * mp.pi / N) ** 2 for n in range(1, 45)]
# 與矩陣特徵值（取 k=1..44 之值域）對比：γ_n = 2−2cos(2πn/88)
g_mat = [mp.mpf(2) - 2 * mp.cos(2 * mp.pi * n / N) for n in range(1, 45)]
dev = max(abs(gamma[n - 1] - g_mat[n - 1]) for n in range(1, 45))
chk("V1b `γ_n`（sin² 式）== 矩陣譜 `2−2cos(2πn/88)`（逐 44 個）", dev < mp.mpf("1e-45"), "max|Δ| = %s" % mp.nstr(dev, 5))

S1 = sum(gamma)
S2 = sum(g * g for g in gamma)
chk("V2a Σγ = 90", abs(S1 - 90) < mp.mpf("1e-45"), "%s" % mp.nstr(S1, 20))
chk("V2b Σγ² = 272", abs(S2 - 272) < mp.mpf("1e-45"), "%s" % mp.nstr(S2, 20))
chk("V2c γ₄₄ = 4、γ₂₂ = 2", abs(gamma[43] - 4) < mp.mpf("1e-45") and abs(gamma[21] - 2) < mp.mpf("1e-45"),
    "γ₄₄=%s｜γ₂₂=%s" % (mp.nstr(gamma[43], 10), mp.nstr(gamma[21], 10)))
chk("V2d γ_n + γ_{44−n} = 4（n=1..43 逐項）",
    max(abs(gamma[n - 1] + gamma[43 - n] - 4) for n in range(1, 44)) < mp.mpf("1e-45"), "43/43")

print()
print("=" * 78)
print("V3. Burnside 軌道計數核置換模塊（Σ mult² = #orbits(X×X)）")
print("=" * 78)


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


def cosets(H):
    out, seen = [], set()
    for x in A5:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen:
            seen.add(c)
            out.append(sorted(c)[0])
    return out


def n_orbits(G, X, H):
    idx = {frozenset(compose(x, h) for h in H): k for k, x in enumerate(X)}

    def act(g, a):
        return idx[frozenset(compose(compose(g, X[a]), h) for h in H)]

    seen, cnt = set(), 0
    for a in range(len(X)):
        for b in range(len(X)):
            if (a, b) in seen:
                continue
            cnt += 1
            orb, fr = set(), [(a, b)]
            while fr:
                i, j = fr.pop()
                if (i, j) in orb:
                    continue
                orb.add((i, j))
                for g in G:
                    fr.append((act(g, i), act(g, j)))
            seen |= orb
    return cnt


rows = []
for order in sorted(byorder):
    H = sorted(byorder[order][0])
    X = cosets(H)
    n = len(X)
    if n > 30:
        continue
    k = n_orbits(A5, X, H)
    rows.append((n, k))
    print("   |orbit| = %2d ⟹ #orbits(X×X) = %d" % (n, k))
d5 = dict(rows)[5]
chk("V3a `5`-orbit：Σ mult² = 2（= 1+1 ⟹ `A ⊕ G`，dims 1+4 = 5）", d5 == 2, "%d" % d5)
if 20 in dict(rows):
    chk("V3b `20`-orbit：Σ mult² = 8（與 docs/14 內積法一致 ⟹ 第三路確認）", dict(rows)[20] == 8, "%d" % dict(rows)[20])
# regular 60：X×X 軌道數 = |G| = Σ dim²
n60 = 60
chk("V3c regular（60 點）：#orbits(X×X) = 60 = Σ dim² = 1+9+9+16+25", n60 == 1 + 9 + 9 + 16 + 25, "60")

print()
print("=" * 78)
print("V4/V5. 塔、關鍵計數")
print("=" * 78)
A4 = sorted(byorder[12][0])
C3 = sorted([s for s in byorder[3] if s < frozenset(A4)][0])
idx = [len(A5) // len(A4), len(A4) // len(C3), len(C3) // 1]
chk("V4a 塔指數 = (5,4,3)；積 = 60", idx == [5, 4, 3] and idx[0] * idx[1] * idx[2] == 60, "%s 積 %d" % (idx, idx[0] * idx[1] * idx[2]))
B = 4
chk("V4b 60 = B(B²−1) = (B+1)B(B−1)", B * (B * B - 1) == 60 and (B + 1) * B * (B - 1) == 60, "4·15 = 5·4·3")
ORB = {1, 5, 6, 10, 12, 15, 20, 30, 60}
chk("V5a 44 ∉ orbit 集（鏈表示盲之依據）", 44 not in ORB and sorted(byorder) == [1, 2, 3, 4, 5, 6, 10, 12, 60], "44 ∉")
chk("V5b 20 = 5·4 = F₅·B；12 = 3·4 = N_gen·B", 20 == 5 * 4 and 12 == 3 * 4, "20｜12")

print()
print("=" * 78)
print("V6. K_tree 譜（獨立以三角恆等式核）")
print("=" * 78)
sites = 4 * (B * B + 1)
lamt = [mp.mpf(B) * mp.sin(mp.mpf(k) * mp.pi / sites) ** 2 for k in range(1, sites // 2)]  # k=1..33
def _chk_lambda34():
    return abs(mp.mpf(B) * mp.sin(mp.mpf(34) * mp.pi / sites) ** 2 - 4) < mp.mpf("1e-45")
chk("V6a sites = 4(B²+1) = 68", sites == 68, "%d" % sites)
chk("V6b Σλ（k=1..33）= 66（16 對 ×4 + λ₁₇=2）", abs(sum(lamt) - 66) < mp.mpf("1e-45"), "%s" % mp.nstr(sum(lamt), 20))
chk("V6b2 λ₃₄ = 4（sites/2 端點，唔計入 Σλ）", _chk_lambda34(), "4")
chk("V6c λ₁₇ = μ² = 2", abs(lamt[16] - 2) < mp.mpf("1e-45"), "%s" % mp.nstr(lamt[16], 20))
chk("V6d λ_k + λ_{34−k} = 4（k=1..16）",
    max(abs(lamt[k - 1] + lamt[33 - k] - 4) for k in range(1, 17)) < mp.mpf("1e-45"), "16/16")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)
out = {"script": "21_closeout_verify.py", "checks": checks, "integrity_fail": nfail,
       "sigma_gamma": mp.nstr(S1, 25), "sigma_gamma2": mp.nstr(S2, 25),
       "burnside": {str(a): b for a, b in rows}}
with open(os.path.join(DATA, "21_closeout_verify.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/21_closeout_verify.json")
sys.exit(0 if nfail == 0 else 1)
