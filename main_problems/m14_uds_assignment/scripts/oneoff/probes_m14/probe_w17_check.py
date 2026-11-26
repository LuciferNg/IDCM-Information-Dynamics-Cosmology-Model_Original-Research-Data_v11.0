#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""核 §十一 之「w ∈ ℤ，17 值 [−10,+30]」與頂點整數 W_table（值域 {−1,0,1,2,3}）之關係
   方法：以 3↔5 建構（同 scripts 之正規化）取 amp，反算 w_rec = amp·90·√(NT3·NT5)/(6·emb)
   對照：W_table(v)（16 遮罩版）於同一簽名
"""
import itertools
import json
import os
from collections import Counter

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

N, MOD = 44, 90
G = np.array([4.0 * np.sin(np.pi * k / 90) ** 2 for k in range(0, N + 2)])


def nrm(t):
    return float(np.sqrt(sum(G[a] for a in t)))


SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]


def W_table(v):
    s = 0
    for c, sg in zip(CS, SIG):
        if sum(a * b for a, b in zip(v, sg)) % MOD == 0:
            s += c
    return s


ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))
QP = list(itertools.combinations_with_replacement(range(1, N + 1), 5))
RANK5 = {t: i for i, t in enumerate(QP)}
NT5 = np.array([nrm(t) for t in QP])
from collections import Counter as C  # noqa: E402
import math as _math


def rank5(Qs):
    return np.array([RANK5[tuple(int(x) for x in q)] for q in Qs], dtype=np.int32)


rng = np.random.default_rng(20261126)
recs = []
for _ in range(4000):
    rem = int(rng.integers(1, N + 1))
    tri = [int(x) for x in rng.choice(np.arange(1, N + 1), size=2, replace=True)]  # surv pair
    A = np.column_stack([np.full(len(ADD3), rem, dtype=np.int32), np.array(ADD3, dtype=np.int32)])
    # W4v 等價：用 W_table 於 (rem, add1, add2, add3)
    w4 = np.array([W_table((rem, a, b, c)) for a, b, c in ADD3], dtype=np.int64)
    nz = np.nonzero(w4)[0]
    if len(nz) == 0:
        continue
    Qs = np.sort(np.column_stack([np.repeat(np.array(tri, dtype=np.int32)[None, :], len(nz), axis=0),
                                  np.array(ADD3, dtype=np.int32)[nz]]), axis=1)
    r = rank5(Qs)
    cs = C(tri)
    emb = np.ones(len(nz))
    for m, cnt in cs.items():
        matches = (Qs == m).sum(1)
        emb *= np.array([float(_math.comb(int(x), cnt)) for x in matches])
    NT3 = nrm(tri)
    amp = 6.0 * (w4[nz] / 90.0) * emb / np.sqrt(NT3 * NT5[r])
    wrec = amp * 90.0 * np.sqrt(NT3 * NT5[r]) / (6.0 * emb)
    for k in range(min(3, len(nz))):
        recs.append(dict(rem=rem, tri=tri, add=list(ADD3[nz[k]]), q=int(r[k]),
                         w_table=int(w4[nz[k]]), w_rec=float(wrec[k]), amp=float(amp[k])))

vals = Counter(round(x["w_rec"], 9) for x in recs)
bad = sum(1 for x in recs if abs(x["w_rec"] - x["w_table"]) > 1e-9)
print(f"樣本 {len(recs)} 條路徑")
print(f"w_table 值集合：{sorted(set(x['w_table'] for x in recs))}")
print(f"w_rec   值集合（{len(vals)} 個）：{sorted(vals)[:20]}")
print(f"w_rec == w_table 之不符合數：{bad}")
print("前 6 例：")
for x in recs[:6]:
    print("  ", x)
json.dump({"n_rec": len(recs), "w_table_values": sorted(set(x["w_table"] for x in recs)),
           "w_rec_values": sorted(vals), "n_mismatch": bad, "sample": recs[:50]},
          open("/home/wsl/idcm_ckpt/e8/mrB/w17_check.json", "w"), indent=1, ensure_ascii=False)
print("已寫 /home/wsl/idcm_ckpt/e8/mrB/w17_check.json")
