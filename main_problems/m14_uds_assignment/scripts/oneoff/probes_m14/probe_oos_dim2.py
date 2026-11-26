#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
OOS-2（跑前凍結）：維數 ⟹ δ 巧合數之【單調性】檢定（四次鏈）
  已測：m=90（γ∈ℚ(ζ₉₀)^+，dim=φ(90)/2=12）⟹ distinct δ = 22
        m=88（dim=φ(176)/2=40）              ⟹ distinct δ = 28
  ★ 凍結預言（單調：dim 越大 ⟹ δ 巧合越少 ⟹ distinct δ 越多）：
     m=45（dim=φ(45)/2=6）  ⟹ distinct δ ≤ 22（巧合不少於 m=90）
     m=180（dim=φ(180)/2=24）⟹ distinct δ ≥ 22（巧合不多於 m=90）
  判定：兩條皆須成立；任一 FAIL ⟹ 單調律否。
  次級讀數（不作主判）：EQ distinct、最大重數。
"""
import itertools
import json
import os

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

mp.mp.dps = 45
N = 44
OUT = "/home/wsl/idcm_ckpt/e8/mrB/oos_dim2.json"


def build(m):
    G = np.array([4.0 * np.sin(np.pi * k / m) ** 2 for k in range(0, N + 2)])
    QP = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 5)), dtype=np.int32)
    EQ = np.zeros(len(QP), dtype=float)
    for t in range(5):
        EQ += G[QP[:, t]]
    return G, QP, EQ


def census(m, tag):
    G, QP, EQ = build(m)
    srt = np.sort(EQ)
    d = np.diff(srt)
    cuts = np.nonzero(d > 1e-12)[0]
    ndist = len(np.concatenate([[srt[0]], srt[cuts + 1]]))
    sizes = np.diff(np.concatenate([[-1], cuts, [len(srt) - 1]]))
    offs = []
    for n in range(1, N + 1):
        gn = float(G[n])
        j = int(np.argmin(np.abs(EQ - gn)))
        idx = [int(x) for x in QP[j]]
        e = sum(mp.mpf(4) * mp.sin(mp.pi * mp.mpf(k) / m) ** 2 for k in idx)
        g = mp.mpf(4) * mp.sin(mp.pi * mp.mpf(n) / m) ** 2
        offs.append((n, j, idx, e - g))
    groups = []
    for n, j, idx, off in offs:
        for grp in groups:
            if abs(off - grp["off"]) < mp.mpf(10) ** -30:
                grp["ns"].append(n)
                break
        else:
            groups.append({"off": off, "ns": [n]})
    print(f"[{tag}] 模數 {m}｜EQ distinct {ndist:,}｜最大重數 {int(sizes.max())}｜相異 δ = {len(groups)} / {N}")
    for grp in sorted(groups, key=lambda x: -len(x["ns"]))[:5]:
        print("   δ={:>28s}｜n={}".format(mp.nstr(grp["off"], 22), grp["ns"]))
    return dict(modulus=m, eq_distinct=int(ndist), max_mult=int(sizes.max()), n_distinct_delta=len(groups),
                groups=[{"offset": mp.nstr(g["off"], 26), "ns": g["ns"]} for g in groups])


print("【OOS-2】維數 ⟹ δ 巧合數 單調性（凍結預言：m45 ≤ 22；m180 ≥ 22）")
a = census(45, "m45/dim6")
b = census(180, "m180/dim24")
ok1 = a["n_distinct_delta"] <= 22
ok2 = b["n_distinct_delta"] >= 22
print(f"\n★ m45 = {a['n_distinct_delta']}（預言 ≤22）⟹ {'PASS' if ok1 else 'FAIL'}")
print(f"★ m180 = {b['n_distinct_delta']}（預言 ≥22）⟹ {'PASS' if ok2 else 'FAIL'}")
print("★ 總判：{}".format("PASS" if (ok1 and ok2) else "FAIL"))
print("  對照：m45 dim6 = {}｜m90 dim12 = 22｜m180 dim24 = {}｜m88 dim40 = 28".format(
    a["n_distinct_delta"], b["n_distinct_delta"]))
json.dump({"m45_dim6": a, "m180_dim24": b, "ref_m90_dim12": 22, "ref_m88_dim40": 28,
           "verdict": "PASS" if (ok1 and ok2) else "FAIL",
           "pred_m45": "<=22", "pred_m180": ">=22"}, open(OUT, "w"), indent=1, ensure_ascii=False)
print("已寫", OUT)
