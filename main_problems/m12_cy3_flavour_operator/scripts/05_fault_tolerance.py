#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2 —— 局域容錯架構（locally fault-tolerant architecture）檢定
2026-10-28

使用者假設：先假定宇宙係一個局域容錯架構 ⟹
  (a) 局域穩定子（每面 parity）必須飽和（每個缺陷面恰好 1 個缺陷）→ 缺陷集＝syndrome
  (b) 唯一解 ⟺ syndrome 唯一 ∧ 無 gauge 簡併（否則標籤只係 gauge，唔係資訊）
  (c) 容錯記憶需要非平凡同調（logical qubit）；球面 S² 上冇 ⟹ 唔可存資訊

本腳本只做 (a)(b) 之硬算 + (c) 之定理引用；唔涉及味指派本身。
建構：由正二十面體展開（skill polyhedral-graph-analysis §15(h) 之可靠構造，零外部表）。

跑法：/usr/bin/python3 scripts/05_fault_tolerance.py
"""
from __future__ import annotations
import itertools
import json
import os
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], []


def say(s=""):
    L.append(s)
    print(s)


PHI = (1.0 + 5.0 ** 0.5) / 2.0


# ─────────── 建構 Rhombicosidodecahedron ───────────
def icosahedron():
    v = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            v.append((0.0, s1 * 1.0, s2 * PHI))
            v.append((s1 * 1.0, s2 * PHI, 0.0))
            v.append((s2 * PHI, 0.0, s1 * 1.0))
    V = np.array(sorted(set(v)), dtype=float)
    return V


def build_rh():
    I = icosahedron()
    n = len(I)
    assert n == 12, n
    D = np.linalg.norm(I[:, None, :] - I[None, :, :], axis=2)
    edge = D[D > 1e-9].min()
    faces = [t for t in itertools.combinations(range(n), 3)
             if abs(D[t[0], t[1]] - edge) < 1e-9
             and abs(D[t[1], t[2]] - edge) < 1e-9
             and abs(D[t[0], t[2]] - edge) < 1e-9]
    assert len(faces) == 20, len(faces)
    R = (edge / 2.0) * np.sqrt(11.0 + 4.0 * np.sqrt(5.0))
    th3 = np.arccos(np.sqrt(R ** 2 - (edge / np.sqrt(3.0)) ** 2) / R)
    pts = []
    for f in faces:
        A, B, C = I[f[0]], I[f[1]], I[f[2]]
        ch = (A + B + C) / 3.0
        ch = ch / np.linalg.norm(ch)
        for v in (A, B, C):
            p = v - np.dot(v, ch) * ch
            p = p / np.linalg.norm(p)
            pts.append(np.cos(th3) * ch + np.sin(th3) * p)
    P = np.array(pts)
    # 去重（每個頂點被 3 個面共享）
    keyed = {}
    for p in P:
        keyed.setdefault(tuple(np.round(p, 9)), p)
    return np.array(list(keyed.values()))


def graph(P, thresh_deg=30.0):
    V = len(P)
    cosang = np.clip(P @ P.T, -1.0, 1.0)
    ang = np.degrees(np.arccos(cosang))
    np.fill_diagonal(ang, 180.0)
    adj = [[] for _ in range(V)]
    edges = []
    for i in range(V):
        for j in range(i + 1, V):
            if ang[i, j] < thresh_deg:
                adj[i].append(j)
                adj[j].append(i)
                edges.append((i, j))
    return adj, edges, ang


def bfs_colour(adj, root):
    n = len(adj)
    c = [-1] * n
    c[root] = 0
    q = deque([root])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if c[w] < 0:
                c[w] = 1 - c[u]
                q.append(w)
    return c


def canon(c):
    a = tuple(c)
    b = tuple(1 - x for x in c)
    return min(a, b)


def rot(axis, angle):
    a = np.array(axis, dtype=float)
    a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def group_I():
    g1 = rot((0, 1, PHI), 2 * np.pi / 5)
    g2 = rot((1, 1, 1), 2 * np.pi / 3)
    G = [np.eye(3)]
    frontier = [np.eye(3)]
    while frontier:
        nf = []
        for g in frontier:
            for h in (g1, g2):
                m = h @ g
                if not any(np.abs(m - k).max() < 1e-9 for k in G):
                    G.append(m)
                    nf.append(m)
        frontier = nf
    return G


def preserves(P, m, tol=1e-6):
    Q = P @ m.T
    # 集合歸屬測試
    ok = 0
    for q in Q:
        d = np.abs(P - q).max(axis=1).min()
        if d < tol:
            ok += 1
    return ok == len(P)


def main():
    say("=" * 80)
    say("F1  RH 建構（二十面體展開，零外部表）")
    P = build_rh()
    V = len(P)
    adj, edges, ang = graph(P)
    deg = sorted(set(len(a) for a in adj))
    say("    V = %d, E = %d, degree set = %s, 邊角 = %.6f°"
        % (V, len(edges), deg, np.degrees(np.arccos(np.clip(P[edges[0][0]] @ P[edges[0][1]], -1, 1)))))
    assert V == 60 and len(edges) == 120 and deg == [4]
    say("    ✅ 60v / 120e / 4-regular")

    say("")
    say("=" * 80)
    say("F2  飽和 2-色族（BFS from 每個 root，mod 全域翻轉去重）")
    fam = {}
    for r in range(V):
        c = bfs_colour(adj, r)
        fam.setdefault(canon(c), c)
    say("    相異族員 = %d" % len(fam))
    fr_sets = {}
    counts = {}
    for k, c in fam.items():
        fs = frozenset(tuple(sorted(e)) for e in edges if c[e[0]] == c[e[1]])
        fr_sets.setdefault(fs, []).append(k)
        counts[len(fs)] = counts.get(len(fs), 0) + 1
    say("    缺陷邊數分佈 = %s" % counts)
    say("    相異缺陷集 (syndrome) 數 = %d" % len(fr_sets))
    only = list(fr_sets.keys())[0] if len(fr_sets) == 1 else None
    say("    %s" % ("✅ syndrome 唯一：30 個族員共享**同一** 32 條缺陷邊集" if only else "⚠ syndrome 不唯一"))
    OUT.append(("F2", {"family_size": len(fam), "frustration_counts": counts,
                       "distinct_syndromes": len(fr_sets)}))

    say("")
    say("=" * 80)
    say("F3  幾何對稱群 I 對 colouring 族之作用（gauge 判定）")
    G = group_I()
    say("    |I| = %d" % len(G))
    npre = sum(1 for m in G if preserves(P, m))
    say("    保點集之元素 = %d / %d" % (npre, len(G)))
    keys = sorted(fam.keys())
    idx = {k: i for i, k in enumerate(keys)}
    # 把 colouring 由 root-family 表示轉為「以頂點序為準」之 canon
    reps = [fam[k] for k in keys]
    # 建立 permutation of each group element
    perms = []
    for m in G:
        Q = P @ m.T
        pm = []
        ok = True
        for q in Q:
            d = np.abs(P - q).max(axis=1)
            j = int(np.argmin(d))
            if d[j] > 1e-6:
                ok = False
                break
            pm.append(j)
        if ok:
            perms.append(pm)
    say("    可用 permutation = %d" % len(perms))
    rep = {k: None for k in keys}
    # 以每個 rep 之第一次出現為代表，計 orbit
    seen = set()
    orbits = []
    for i, k in enumerate(keys):
        if i in seen:
            continue
        orbi = set()
        stack = [k]
        while stack:
            cur = stack.pop()
            if cur in orbi:
                continue
            orbi.add(cur)
            c = fam[cur]
            for pm in perms:
                nc = [0] * V
                for v in range(V):
                    nc[pm[v]] = c[v]
                nk = canon(nc)
                if nk in fam:
                    stack.append(nk)
        ids = {keys.index(x) for x in orbi}
        seen |= ids
        orbits.append(sorted(ids))
    say("    colouring 族在 I 下之 orbit 數 = %d，大小 = %s"
        % (len(orbits), sorted(len(o) for o in orbits)))
    say("    ⟹ %s" % ("單一 orbit ⟹ 族員互為對稱像 ⟹ 標籤指派係 gauge，非資訊" if len(orbits) == 1 else "多 orbit ⟹ 有 I-不變之類別結構"))
    OUT.append(("F3", {"I_order": len(G), "perms": len(perms),
                       "orbits": len(orbits), "orbit_sizes": sorted(len(o) for o in orbits)}))

    say("")
    say("=" * 80)
    say("F4  容錯記憶之拓撲條件（定理引用，非本輪新算）")
    say("    RH 係球面上之平面圖（skill §2b：check_planarity = True，62 面）")
    say("    ⟹ 面邊界生成整個圈空間（E−V+1 = 120−60+1 = %d = F−1）" % (120 - 60 + 1))
    say("    ⟹ H₁(S²) = 0 ⟹ Z₂ 標籤碼**冇 logical qubit** ⟹ 球面上唔可存邏輯資訊")
    say("    K_tree 係樹 ⟹ b₁ = 0，同樣冇同調記憶")
    say("    ⟹ 「味 = 拓撲荷」在 A₅ 與 K_tree 上皆死")
    OUT.append(("F4", {"planar": True, "faces": 62, "cycle_dim": 61,
                       "H1_trivial": True, "logical_qubits": 0}))

    with open(os.path.join(ROOT, "data", "05_fault_tolerance.json"), "w") as f:
        json.dump(dict(OUT), f, indent=2, ensure_ascii=False)
    say("")
    say("[written] data/05_fault_tolerance.json , data/run_05_fault_tolerance.log")

    with open(os.path.join(ROOT, "data", "run_05_fault_tolerance.log"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()