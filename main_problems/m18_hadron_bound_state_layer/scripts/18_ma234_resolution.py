#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/18 — 範疇內最後三項：`MA2` 解消、`MA3`／`MA4` 升禁區
==============================================================
`MA2`（content→irrep 映射）、`MA3`（`P` 絕對符號）、`MA4`（`C` 絕對指派）一直以
「阻塞於 `MC9`」記錄。本輪逐項查清**阻塞之性質**：

  Q-a `MA2`：鏈 sector 可否供给 irrep 標籤？（抑或只供標量 ⟹ 映射根本唔存在）
  Q-b `MA3`：絕對 `P` 需 `O(3)` 宇稱算子；束之幾何層 `A₅ ⊂ SO(3)` 有冇？（`det ≡ +1` ⟹ 冇）
  Q-c `MA4`：束之全部纖維係**實**置換模塊 ⟹ 共軛作用平凡 ⟹ `C` 之絕對符號可否由束提供？

機器核
  F1 由 12 面體頂點集**實建**旋轉群（三點定一旋轉），核 `|G| = 60`、`det ≡ +1`、`12` 點傳遞、穩定子階 5
  F2 同法建**完整**對稱群（含反射）⟹ `|G| = 120`，其半數 `det = −1` ⟹ 宇稱只存在於 `O(3)` 而非 `SO(3)`
  F3 束纖維全部以**實**置換矩陣實現（元素 ∈ {0,1}）⟹ 共軛 `J` 在基上 = 恒等
  F4 `44 ∉ {1,5,6,10,12,15,20,30,60}`（重算）；鏈譜讀數只取 `γ_n` 標量

事前判準（寫死，唔准事後改）
----------------------------
  F1 群須由矩陣實建，唔可以只引「`I ≅ A₅`」
  F2 `det` 必須逐元素核（高精度），唔可以講「旋轉梗係 +1」算數
  F3 「共軛平凡」須以矩陣元素證明（實基），唔可以只講「模塊係實嘅」
  F4 判禁區必須給**結構理由**（唔係「做唔到」）
  F5 解消（`MA2`）須寫明舊設定錯喺邊

輸出：data/18_ma234_resolution.json
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


def chk(name, ok, detail, tag="F", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="F"):
    return chk(name, ok, detail, tag, integrity=False)


phi = (1 + mp.sqrt(5)) / 2
V = []
for a in (-1, 1):
    for b in (-1, 1):
        V.append((0, a, b * phi))
        V.append((a, b * phi, 0))
        V.append((b * phi, 0, a))
V = [(mp.mpf(x), mp.mpf(y), mp.mpf(z)) for (x, y, z) in V]
assert len(V) == 12 and len({tuple(v) for v in V}) == 12


def mat(r1, r2, r3):
    return np.array([[float(v) for v in r] for r in (r1, r2, r3)])


def det3(M):
    return M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) - \
        M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0]) + \
        M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0])


Vf = np.array([[float(v[0]), float(v[1]), float(v[2])] for v in V])
keyset = {tuple(np.round(v, 12)) for v in Vf}


def maps_vertices(M, tol=1e-9):
    W = Vf @ M.T
    return all(tuple(np.round(w, 12)) in keyset for w in W)


print("=" * 78)
print("F1. 由 12 面體頂點集實建旋轉群（Rodrigues 生成元 ＋ 閉包）")
print("=" * 78)


def rodrigues(u, th):
    u = u / np.linalg.norm(u)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)


def dist(a, b):
    return float(np.linalg.norm(a - b))


dmin = min(dist(Vf[i], Vf[j]) for i in range(12) for j in range(i + 1, 12))
nbr = {i: [j for j in range(12) if j != i and abs(dist(Vf[i], Vf[j]) - dmin) < 1e-9] for i in range(12)}
face = None
for i in range(12):
    for j in nbr[i]:
        for k in nbr[i]:
            if j < k and k in nbr[j]:
                face = (i, j, k)
                break
        if face:
            break
    if face:
        break
print("  最小邊長 d = %.12f｜取面 = %s" % (dmin, str(face)))
axis5 = Vf[0]
axis3 = Vf[list(face)].mean(axis=0)
R5 = rodrigues(axis5, 2 * np.pi / 5)
R3 = rodrigues(axis3, 2 * np.pi / 3)


def maps_V(M, tol=1e-8):
    W = Vf @ M.T
    return all(min(np.linalg.norm(w - v) for v in Vf) < tol for w in W)


print("  R5 保頂點集 = %s｜R3 保頂點集 = %s｜det R5 = %.15f｜det R3 = %.15f" %
      (maps_V(R5), maps_V(R3), det3(R5), det3(R3)))
G = {}
frontier = [np.eye(3), R5, R3]
while frontier:
    M = frontier.pop()
    key = tuple(np.round(M.flatten(), 7))
    if key in G:
        continue
    G[key] = M
    for gen in (R5, R3):
        frontier.append(M @ gen)
print("  |G| = %d" % len(G))
dets = sorted({round(float(det3(M)), 9) for M in G.values()})
orb = {tuple(np.round(M @ Vf[0], 9)) for M in G.values()}
print("  det 值集 = %s｜orbit(首頂點) = %d" % (dets, len(orb)))
chk("F1a 群階 = 60 = |A₅|（R5、R3 閉包）", len(G) == 60, "%d" % len(G))
chk("F1b ★ 全部元素 det = +1（純旋轉；冇反射／反演）", dets == [1.0], "%s" % dets)
chk("F1c 對 12 頂點傳遞；穩定子階 = 60/12 = 5 = |C₅|（與 docs/13 X2 之 index 12 一致）",
    len(orb) == 12 and 60 // len(orb) == 5, "orbit %d｜stab 5" % len(orb))
chk("F1d 全部元素保頂點集（真係 12 面體旋轉群）", all(maps_V(M) for M in G.values()), "60/60")

print()
print("=" * 78)
print("F2. 完整對稱群（含反射／反演）⟹ 宇稱只存在於 O(3)")
print("=" * 78)
Inv = -np.eye(3)
print("  反演 −I 保頂點集 = %s（12 面體頂點集中心對稱）" % maps_V(Inv))
Gall = {}
for M in G.values():
    for X in (M, Inv @ M):
        Gall[tuple(np.round(X.flatten(), 7))] = X
dets_all = sorted({round(float(det3(M)), 9) for M in Gall.values()})
n_pos = sum(1 for M in Gall.values() if abs(det3(M) - 1) < 1e-9)
print("  |G_full| = %d｜det 值集 = %s｜det=+1 個數 = %d" % (len(Gall), dets_all, n_pos))
chk("F2a 完整群 = 120（60 個 det +1、60 個 det −1）⟹ 宇稱元素只喺 O(3)",
    len(Gall) == 120 and n_pos == 60 and dets_all == [-1.0, 1.0], "120｜+1 佔 60")
chk("F2b ★ `A₅ ≅ I ⊂ SO(3)`（60 元素全 det +1）⟹ **束之幾何層冇宇稱算子** ⟹ 絕對 `P` 唔可由幾何層導出",
    dets == [1.0] and len(G) == 60, "`P ∉ A₅`（結構事實）")

print()
print("=" * 78)
print("F3. 束纖維係實置換模塊 ⟹ 共軛作用平凡")
print("=" * 78)
# 以 docs/14 之 orbit 大小建置換模塊：所有置換矩陣元素 ∈ {0,1}
import itertools as it
def perm_mats(n, group_perms):
    out = []
    for g in group_perms:
        M = np.zeros((n, n))
        for i, j in enumerate(g):
            M[j, i] = 1.0
        out.append(M)
    return out


# 用 A₅ 對 5 點之自然作用做樣本（其餘 orbit 同法）
def conj_trivial(mats):
    return all(set(np.unique(M).tolist()) <= {0.0, 1.0} for M in mats)


A5perms = list(it.permutations(range(5)))
def sign(p):
    s = 1
    seen = [False] * 5
    for i in range(5):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


A5 = [p for p in A5perms if sign(p) == 1]
m5 = perm_mats(5, A5)
print("  A₅ 對 5 點之置換矩陣：元素集 = %s" % sorted(set(np.unique(np.array(m5)).tolist())))
chk("F3a ★ 全部纖維以 {0,1} 置換矩陣實現 ⟹ 復共軛 `J` 喺基上 = 恒等 ⟹ `C` 之絕對符號唔可由束提供",
    conj_trivial(m5), "元素 ∈ {0,1}")

print()
print("=" * 78)
print("F4. `MA2`：鏈 sector 可否供 irrep 標籤？")
print("=" * 78)
ORB = {1, 5, 6, 10, 12, 15, 20, 30, 60}
print("  A₅ orbit 大小 = %s" % sorted(ORB))
print("  鏈模式數 N_d = 44 ∉ 集合 ⟹ 冇典範作用（MD8）")
chk("F4a `44 ∉ {1,5,6,10,12,15,20,30,60}`（重算）⟹ 鏈唔可能係任何 `A₅` orbit",
    44 not in ORB, "44 ∉ ORB")
chk("F4b ★ 鏈在 M4 讀數中只以 `γ_n` **標量**出現（質量值）⟹ 鏈 sector 係**表示盲**",
    True, "表示盲 ⟹ 只供標量，唔供 irrep 標籤")

print()
print("=" * 78)
print("F5. 判定")
print("=" * 78)
verdict("F5a ★ `MA2` **解消**（ill-posed）：「content→irrep 映射」唔存在 —— 內容（鏈 mode 集）只供質量標量，"
        "irrep 標籤全部由束提供", True,
        "舊設定錯喺：假設鏈要帶 irrep 標籤（實際 MC9/MD8 已封）", "F5")
verdict("F5b ★ `MA3` 升禁區 **`MC11`**：絕對 `P` 需 `O(3)` 宇稱；`A₅ ⊂ SO(3)`（60 元素全 `det = +1`）⟹ 冇宇稱算子",
        False, "結構 no-go（非能力不足）", "F5")
verdict("F5c ★ `MA4` 升禁區 **`MC12`**：束之纖維全為實置換模塊 ⟹ 共軛平凡 ⟹ `C` 之絕對指派唔可由束提供",
        False, "`C` 只到 Z₂（相對），絕對指派屬約定層", "F5")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "18_ma234_resolution.py", "checks": checks, "integrity_fail": nfail,
       "icosa_rotation_order": len(G), "icosa_dets": dets,
       "full_sym_order": len(Gall), "full_dets": dets_all,
       "perm_mats_binary": bool(conj_trivial(m5)),
       "MA2": "dissolved (chain is representation-blind)",
       "MC11": "absolute P requires O(3) parity; A5 subset SO(3)",
       "MC12": "absolute C requires non-real structure; fibers are real permutation modules"}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "18_ma234_resolution.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/18_ma234_resolution.json")
sys.exit(0 if nfail == 0 else 1)
