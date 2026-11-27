#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R4：22 個 δ 之支撐 4／5 稀疏關係普查
   形式：δ = Σ_{i∈P} γ_i − Σ_{j∈M} γ_j（|P|+|M| = 4 或 5，P∩M = ∅）
   float64 預篩（tol 1e−9）⟹ mpmath dps 50 精確確認
"""
import itertools
import json

import mpmath as mp
import numpy as np

mp.mp.dps = 50
G = np.array([0.0] + [float(4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2) for k in range(1, 45)])
Gmp = [0] + [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]

d = json.load(open("/home/wsl/idcm_ckpt/e8/mrB/offdelta.json"))
deltas = sorted({r["offset"] for r in d["rows"]})
print(f"相異 δ 數 = {len(deltas)}")
DEL = np.array([float(x) for x in deltas])
tolf = 1e-9

hits = {4: {}, 5: {}}

# ---- 支撐 4 ----
for (a, b, c_, dd) in itertools.combinations(range(1, 45), 4):
    quad = [a, b, c_, dd]
    for pidx in itertools.combinations(range(4), 2):
        s = sum(G[quad[i]] for i in pidx) - sum(G[quad[i]] for i in range(4) if i not in pidx)
        near = np.where(np.abs(DEL - s) < tolf)[0]
        for k in near:
            P = [quad[i] for i in pidx]
            M = [quad[i] for i in range(4) if i not in pidx]
            val = sum(Gmp[i] for i in P) - sum(Gmp[i] for i in M)
            if abs(val - mp.mpf(deltas[k])) < mp.mpf(10) ** -22:
                hits[4].setdefault(deltas[k], []).append((P, M))
print(f"支撐 4 命中 {len(hits[4])} / {len(deltas)}")

# ---- 支撐 5 ----
for comb in itertools.combinations(range(1, 45), 5):
    q = list(comb)
    for np_ in (1, 2):  # 正號個數 1 或 2（等價於 4 負／3 負）
        for pidx in itertools.combinations(range(5), np_):
            s = sum(G[q[i]] for i in pidx) - sum(G[q[i]] for i in range(5) if i not in pidx)
            near = np.where(np.abs(DEL - s) < tolf)[0]
            for k in near:
                P = [q[i] for i in pidx]
                M = [q[i] for i in range(5) if i not in pidx]
                val = sum(Gmp[i] for i in P) - sum(Gmp[i] for i in M)
                if abs(val - mp.mpf(deltas[k])) < mp.mpf(10) ** -22:
                    hits[5].setdefault(deltas[k], []).append((P, M))
print(f"支撐 5 命中 {len(hits[5])} / {len(deltas)}")

out = {"n_deltas": len(deltas),
       "support4_hits": {k: [list(map(list, v)) for v in vs[:3]] for k, vs in hits[4].items()},
       "support5_hits": {k: [list(map(list, v)) for v in vs[:3]] for k, vs in hits[5].items()},
       "n_support4": len(hits[4]), "n_support5": len(hits[5])}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/r4_support45.json", "w"), indent=1, ensure_ascii=False)
for k in deltas:
    print(f"δ={k[:20]}｜4: {hits[4].get(k, [])[:1]}｜5: {hits[5].get(k, [])[:1]}")
