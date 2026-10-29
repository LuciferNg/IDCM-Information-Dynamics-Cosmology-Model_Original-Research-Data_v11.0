#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R3-P2 —— 味源改在 CMB/拓撲層：阻挫應力算子之範疇性檢定

背景（使用者提議）：味既不在譜、不在樹、亦不在 CY₃ ⟹ 可能在 CMB 層由手性／阻挫應力生成。

事前判準（開算前寫死）：
  C1【尖銳】阻挫簽名鄰接 A_s（±1）之 nullity_ℚ == 9 ?  （9 = O137 之 9 模式數 = N_gen²）
       → 若等於 9 ⟹ 「9 模式 = 阻挫算子之零模」成立（命中）；否則落空。
  C2  rank_ℚ(A_s) ≥ 4 ?（承載 4 通道之必要條件；對照 κ 之 rank 2）
  C3【修正，開算前】3 重簡併之地位：30 個基態色為單一 I-orbit（|orbit|=30, 穩定子階 2）
       ⟹ A_s 只與階 2 子群交換，**不**與整個 I 交換 ⟹ 其譜**不**按 I irrep 分塊，
       ⟹ 重數不受對稱強制 ⟹ 若見 3 重簡併係**偶然／非平凡**，可報（但仍需結構預言才可作證據）。
       （初稿誤以為 A_s 為 I-對稱；此處更正。）
  C4  譜之重數分佈：報告全部；若重數皆屬 {1,2,3,4,5} 亦不足為奇（60 維、多巧合）。
  C5  每頂點阻挫度之分佈（若只有 3 個值 ⟹ 範疇性三分，記為數據非命中）
  C6  無符號鄰接 A 之譜作對照

跑法：/usr/bin/python3 scripts/13_frustration_operator.py
"""
from __future__ import annotations
import itertools
import json
import os
from collections import Counter, deque

from fractions import Fraction

import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


def icosahedron():
    v = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            v.append((0.0, s1 * 1.0, s2 * PHI))
            v.append((s1 * 1.0, s2 * PHI, 0.0))
            v.append((s2 * PHI, 0.0, s1 * 1.0))
    return np.array(sorted(set(v)), dtype=float)


def build_rh():
    I = icosahedron()
    n = len(I)
    D = np.linalg.norm(I[:, None, :] - I[None, :, :], axis=2)
    edge = D[D > 1e-9].min()
    faces = [t for t in itertools.combinations(range(n), 3)
             if abs(D[t[0], t[1]] - edge) < 1e-9
             and abs(D[t[1], t[2]] - edge) < 1e-9
             and abs(D[t[0], t[2]] - edge) < 1e-9]
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
    keyed = {}
    for p in np.array(pts):
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
    return adj, edges


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


def exact_rank(M):
    """精確 rank_ℚ：Fraction 高斯消去。M = list of list of int。"""
    rows = [[Fraction(x) for x in r] for r in M]
    n = len(rows)
    m = len(rows[0])
    r = 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if rows[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][c]
        if pv != 1:
            rows[r] = [x / pv for x in rows[r]]
        for i in range(n):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == n:
            break
    return r


def cluster(ev, tol=1e-7):
    """把排序後特徵值聚成 (值, 重數)。"""
    out = []
    for x in ev:
        if out and abs(x - out[-1][0]) < tol:
            out[-1] = (out[-1][0], out[-1][1] + 1)
        else:
            out.append((x, 1))
    return out


def main():
    say("=" * 84)
    say("R3-P2  RH 重建（零外部表）")
    P = build_rh()
    adj, edges = graph(P)
    V = len(P)
    say("  V = %d, E = %d, degree set = %s" % (V, len(edges), sorted(set(len(a) for a in adj))))
    assert V == 60 and len(edges) == 120

    # 阻挫簽名：BFS 2-色 ⟹ 同色邊 = 阻挫
    c = bfs_colour(adj, 0)
    frust = [e for e in edges if c[e[0]] == c[e[1]]]
    sat = [e for e in edges if c[e[0]] != c[e[1]]]
    say("  阻挫邊 = %d ｜滿足邊 = %d ｜Ising H = −%d + %d = %d"
        % (len(frust), len(sat), len(sat), len(frust), len(frust) - len(sat)))

    say("")
    say("=" * 84)
    say("C5  每頂點阻挫度分佈")
    fdeg = Counter()
    for u, v in frust:
        fdeg[u] += 1
        fdeg[v] += 1
    dist = Counter(fdeg.get(i, 0) for i in range(V))
    say("  阻挫度 → 頂點數：%s" % dict(sorted(dist.items())))
    say("  ⟹ 相異值數 = %d ⟹ %s"
        % (len(dist), "★ 範疇性三分（可作 3 味之標籤，但僅為標籤）" if len(dist) == 3
           else "非三分"))
    OUT["C5"] = {"frustration_degree_hist": {str(k): v for k, v in sorted(dist.items())}}

    say("")
    say("=" * 84)
    say("C1/C2  A_s（±1 簽名鄰接）之精確秩與零模數")
    As = [[0] * V for _ in range(V)]
    for u, v in frust:
        As[u][v] = 1
        As[v][u] = 1
    for u, v in sat:
        As[u][v] = -1
        As[v][u] = -1
    rk = exact_rank(As)
    null = V - rk
    say("  rank_ℚ(A_s) = %d / %d ⟹ nullity = %d" % (rk, V, null))
    say("  C1 判準 nullity == 9 ⟹ %s"
        % ("★ 命中：9 模式 = 阻挫零模" if null == 9 else "🔴 落空（實得 %d ≠ 9）" % null))
    say("  C2 判準 rank ≥ 4 ⟹ %s（承載 4 通道之容量）"
        % ("✅ 通過" if rk >= 4 else "🔴 不足"))
    OUT["C1C2"] = {"rank": rk, "nullity": null, "nullity_eq_9": null == 9}

    say("")
    say("=" * 84)
    say("C3/C4  譜與重數（A_s 與 A）")
    Asf = np.array(As, dtype=float)
    ev_s = np.sort(np.linalg.eigvalsh(Asf))
    cl_s = cluster(ev_s)
    say("  A_s 譜（%d 個相異值）：" % len(cl_s))
    tot = 0
    for val, m in cl_s:
        say("     λ = %+9.6f   重數 %d" % (val, m))
        tot += m
    say("  重數總和 = %d（應 = 60）" % tot)
    irrep = {1, 3, 4, 5}
    bad = [(v, m) for v, m in cl_s if m not in irrep]
    say("  重數皆屬 I 之不可約維數 {1,3,4,5}？%s" % ("✅ 是 ⟹ I-對稱（譜按 irrep 分塊）"
                                              if not bad else "🔴 否：%s" % bad))
    say("  C3 之地位：A_s 僅與階 2 穩定子交換 ⟹ 重數非對稱強制 ⟹ 3 重簡併若出現即非平凡")
    OUT["C3C4"] = {"distinct_values": len(cl_s),
                   "multiplicities": [[float(v), m] for v, m in cl_s]}

    E = {tuple(sorted(e)) for e in edges}
    Adj = np.array([[1.0 if (min(i, j), max(i, j)) in E else 0.0 for j in range(V)]
                    for i in range(V)])
    ev_a = np.sort(np.linalg.eigvalsh(Adj))
    cl_a = cluster(ev_a)
    say("")
    say("  A（無符號，4-regular，對照）譜（%d 個相異值）：" % len(cl_a))
    for val, m in cl_a:
        say("     λ = %+9.6f   重數 %d" % (val, m))
    OUT["C6"] = {"unsigned_multiplicities": [[float(v), m] for v, m in cl_a]}

    say("")
    say("=" * 84)
    say("D1  【關鍵】A_s 是否規範等價於 A？  檢 A_s == S A S")
    say("")
    s = [1 if x == 0 else -1 for x in c]        # s_i = (−1)^{c_i}
    AdjI = np.array(As, dtype=float)
    Ab = np.array([[1.0 if (min(i, j), max(i, j)) in E else 0.0 for j in range(V)]
                   for i in range(V)])
    S = np.diag(np.array(s, dtype=float))
    SAS = S @ Ab @ S
    gap = float(np.abs(SAS - AdjI).max())
    say("  max|S A S − A_s| = %.1f ⟹ %s" % (gap, "✅ 恆等（A_s = S A S）" if gap == 0.0 else "⚠ 非恆等"))
    say("  ⟹ 由 2-色所簽之阻挫算子係 A 之么正相似 ⟹ 譜相同 ⟹ **零新資訊**")
    say("  ⟹ 「阻挫應力算子之譜」作味源 = 🔴 作廢（同 A₅ 頂點標籤之 gauge 論證）")
    OUT["D1"] = {"max_abs_gap": gap, "identity": gap == 0.0}

    say("")
    say("=" * 84)
    say("D2  三分 8/40/12 是否為非規範不變量？（掃全部 30 個最優着色）")
    say("")
    fam = {}
    for r in range(V):
        cc = bfs_colour(adj, r)
        fam.setdefault(min(tuple(cc), tuple(1 - x for x in cc)), cc)
    variants = Counter()
    fcounts = set()
    for key, cc in fam.items():
        fr = [e for e in edges if cc[e[0]] == cc[e[1]]]
        fcounts.add(len(fr))
        fd = Counter()
        for u, v in fr:
            fd[u] += 1
            fd[v] += 1
        variants[tuple(sorted(Counter(fd.get(i, 0) for i in range(V)).items()))] += 1
    say("  相異着色 = %d ｜阻挫邊數集合 = %s" % (len(fam), sorted(fcounts)))
    say("  阻挫度分佈變體數 = %d" % len(variants))
    for d, k in variants.items():
        say("     %s  ×%d" % (d, k))
    say("  ⟹ %s"
        % ("✅ 30/30 同一分佈 ⟹ 三分係非規範不變量" if len(variants) == 1
           else "🔴 隨着色改變 ⟹ 非不變量"))
    OUT["D2"] = {"n_colourings": len(fam), "variants": len(variants),
                 "invariant": len(variants) == 1}

    say("")
    say("=" * 84)
    say("R3-P2 判詞")
    say("")
    say("  C1（nullity = 9）：%s" % ("★ 命中" if null == 9 else "🔴 落空，實得 %d" % null))
    say("  C2（rank ≥ 4）：%s" % ("✅" if rk >= 4 else "🔴"))
    say("  C3（3 重簡併）：A_s 僅與階 2 穩定子交換 ⟹ 重數非對稱強制；且 A_s = S A S ⟹ 簡併即 A 之簡併")
    say("  C5（三分）：%d 個阻挫度值，30/30 不變 ⟹ 非規範不變量，可作標籤（無力學內容）" % len(dist))
    say("  🔴 D1：阻挫簽名算子 = A 之么正相似 ⟹ 譜零新資訊 ⟹ 此路線作廢")
    say("  ⟹ 殘存唯一載體 = 三分 8/40/12；由標籤到數值仍缺**入射規則**。")


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/13_frustration_operator.json , data/run_13_frustration_operator.log")
    with open(os.path.join(ROOT, "data", "13_frustration_operator.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False, default=lambda o: str(o))
    with open(os.path.join(ROOT, "data", "run_13_frustration_operator.log"), "w") as f:
        f.write("\n".join(L) + "\n")