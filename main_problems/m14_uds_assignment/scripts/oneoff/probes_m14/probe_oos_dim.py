#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
OOS 樣本外測試（跑前凍結）：機制「δ 之量子化程度由 ℚ(ζ)^+ 之維數決定」
  鏈 A（/90，已測）：γ_n = 4sin²(πn/90)，γ ∈ ℚ(ζ₉₀)^+，dim = φ(90)/2 = 12
                     ⟹ 44 個值塞 12 維 ⟹ 觀測 distinct δ = 22
  ★ 事前凍結預言（鏈 B，/88）：γ_n = 4sin²(πn/88)，γ ∈ ℚ(ζ₁₇₆)^+，dim = φ(176)/2 = 40
                     ⟹ 44 個值塞 40 維 ⟹ 預言 distinct δ ≥ 40（即幾乎無精確巧合；上限 44）
  判準：B 鏈 distinct δ 之計數；A 鏈重算作對照（須復現 22）
  另計：兩鏈之 EQ distinct 數、最大重數（作次級讀數，不作主判）
"""
import itertools
import json
import os
import sys

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

mp.mp.dps = 45

N = 44
OUT = "/home/wsl/idcm_ckpt/e8/mrB/oos_dim.json"


def build(m):
    """m = 鏈模數（90 或 88）：γ_k = 4 sin²(πk/m)"""
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
    print(f"[{tag}] 模數 {m}｜五模態 {len(EQ):,}｜EQ distinct {ndist:,}｜最大重數 {int(sizes.max())}")
    # 逐 n 之最近模偏移（dps 45 精算）
    offs = []
    for n in range(1, N + 1):
        gn = float(G[n])
        j = int(np.argmin(np.abs(EQ - gn)))
        idx = [int(x) for x in QP[j]]
        e = sum(mp.mpf(4) * mp.sin(mp.pi * mp.mpf(k) / m) ** 2 for k in idx)
        g = mp.mpf(4) * mp.sin(mp.pi * mp.mpf(n) / m) ** 2
        offs.append((n, j, idx, e - g))
    # 精確分組（dps 45 下相鄰差 < 1e-30 視為相等）
    groups = []
    for n, j, idx, off in offs:
        placed = False
        for grp in groups:
            if abs(off - grp["off"]) < mp.mpf(10) ** -30:
                grp["ns"].append(n)
                placed = True
                break
        if not placed:
            groups.append({"off": off, "ns": [n]})
    print(f"[{tag}] 相異 δ 值 = {len(groups)} / {N}")
    for grp in sorted(groups, key=lambda x: -len(x["ns"]))[:6]:
        print("   δ={:>28s}｜n={}".format(mp.nstr(grp["off"], 22), grp["ns"]))
    return dict(modulus=m, n_states=int(len(EQ)), eq_distinct=int(ndist), max_mult=int(sizes.max()),
                n_distinct_delta=len(groups),
                groups=[{"offset": mp.nstr(g["off"], 26), "ns": g["ns"]} for g in groups],
                rows=[{"n": n, "j": j, "set": idx, "offset": mp.nstr(off, 26)} for n, j, idx, off in offs])


print("【OOS】ℚ-維數 ⟹ distinct δ 之樣本外檢定（跑前凍結：預言 B 鏈 ≥ 40）")
res = {"A_mod90": census(90, "A/90"), "B_mod88": census(88, "B/88")}
pred = "B ≥ 40"
got = res["B_mod88"]["n_distinct_delta"]
print(f"\n★ 預言（凍結）：{pred}｜實測 B = {got}｜A 對照 = {res['A_mod90']['n_distinct_delta']}")
print("★ 判定：{}".format("PASS" if got >= 40 else "FAIL"))
res["prediction"] = pred
res["verdict"] = "PASS" if got >= 40 else "FAIL"
with open(OUT, "w") as f:
    json.dump(res, f, indent=1, ensure_ascii=False)
print("已寫", OUT)
