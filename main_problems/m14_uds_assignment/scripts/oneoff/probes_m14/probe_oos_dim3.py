#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
OOS-3（跑前凍結）：機制修正版 —— 控制量係【內部反射關係數】，非 ℚ-維數
  γ_k = 4 sin²(πk/m) 之兩條精確反射關係：
    R1：γ_k = γ_{m−k}          （同值；內部 ⟺ m−k ∈ [1,44]）
    R2：γ_k + γ_{m/2−k} = 4    （和 4；內部 ⟺ m/2−k ∈ [1,44]，m 須偶）
  已測四鏈（distinct δ）與關係數：
    m=45 ：R1 22 對（44 值塌成 22）｜R2 無（m 奇）        ⟹ δ 20
    m=90 ：R1 無｜R2 22 對                                ⟹ δ 22
    m=88 ：R1 僅 k=44 自對｜R2 21 對 + k=22 自對           ⟹ δ 28
    m=180：R1 無｜R2 無（90−k ≥ 46）                       ⟹ δ 42
  ⟹ 關係越多 ⟹ 精確巧合越多 ⟹ distinct δ 越少
  ★ 凍結預言：
    m=60 ：R1 14 對（k=16..29）＋k=30 自對｜R2 29 對（k=1..29）⟹ 關係極多 ⟹ distinct δ ≤ 22
    m=176：R1 無｜R2 僅 k=44 自對（88−44=44）⟹ 關係近乎零 ⟹ distinct δ ≥ 40
  判定：兩條皆須成立。
"""
import itertools
import json
import os

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

mp.mp.dps = 45
N = 44


def pair_counts(m):
    r1 = sum(1 for k in range(1, N + 1) if 1 <= m - k <= N and m - k > k)
    r1self = sum(1 for k in range(1, N + 1) if m - k == k)
    r2 = 0
    r2self = 0
    if m % 2 == 0:
        h = m // 2
        for k in range(1, N + 1):
            p = h - k
            if 1 <= p <= N:
                if p == k:
                    r2self += 1
                elif p > k:
                    r2 += 1
            elif p == 0:
                pass
    return dict(r1_pairs=r1, r1_self=r1self, r2_pairs=r2, r2_self=r2self)


def census(m, tag):
    G = np.array([4.0 * np.sin(np.pi * k / m) ** 2 for k in range(0, N + 2)])
    QP = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 5)), dtype=np.int32)
    EQ = np.zeros(len(QP), dtype=float)
    for t in range(5):
        EQ += G[QP[:, t]]
    offs = []
    for n in range(1, N + 1):
        j = int(np.argmin(np.abs(EQ - float(G[n]))))
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
    pc = pair_counts(m)
    print("[{}] m={}｜R1 對 {}（自 {}）｜R2 對 {}（自 {}）｜distinct δ = {} / {}｜γ_1 = {:.6f}".format(
        tag, m, pc["r1_pairs"], pc["r1_self"], pc["r2_pairs"], pc["r2_self"], len(groups), N, float(G[1])))
    for grp in sorted(groups, key=lambda x: -len(x["ns"]))[:4]:
        print("   δ={:>26s}｜n={}".format(mp.nstr(grp["off"], 20), grp["ns"]))
    return dict(modulus=m, **pc, n_distinct_delta=len(groups),
                groups=[{"offset": mp.nstr(g["off"], 26), "ns": g["ns"]} for g in groups])


print("【OOS-3】內部反射關係數 ⟹ distinct δ（凍結預言：m60 ≤22；m176 ≥40）")
a = census(60, "m60")
b = census(176, "m176")
ok1 = a["n_distinct_delta"] <= 22
ok2 = b["n_distinct_delta"] >= 40
print(f"\n★ m60 = {a['n_distinct_delta']}（預言 ≤22）⟹ {'PASS' if ok1 else 'FAIL'}")
print(f"★ m176 = {b['n_distinct_delta']}（預言 ≥40）⟹ {'PASS' if ok2 else 'FAIL'}")
print("★ 總判：{}".format("PASS" if (ok1 and ok2) else "FAIL"))
print("  五鏈總覽（關係多→少 ⟹ δ 少→多）：m45 R1×22 → 20｜m90 R2×22 → 22｜m88 R2×21+1 → 28｜"
      "m60 = {}｜m176 = {}｜m180 R×0 → 42".format(a["n_distinct_delta"], b["n_distinct_delta"]))
json.dump({"m60": a, "m176": b, "verdict": "PASS" if (ok1 and ok2) else "FAIL",
           "reference": {"m45": 20, "m90": 22, "m88": 28, "m180": 42}},
          open("/home/wsl/idcm_ckpt/e8/mrB/oos_dim3.json", "w"), indent=1, ensure_ascii=False)
print("已寫 /home/wsl/idcm_ckpt/e8/mrB/oos_dim3.json")
