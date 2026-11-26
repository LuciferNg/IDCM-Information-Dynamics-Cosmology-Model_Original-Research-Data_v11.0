#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""§十七 殘留 2：主導不等式之結構判據（極簡；只用 rank1.json 既有量 + EQ）
   恆等：w* = Cu*/den* ⟹ |Cu*| = |w*|·gap*
         ratio := |Cz*|/(|E₂||w*|) = (|Cz*|/|Cu*|)·gap*/|E₂|
   判據：若 gap* < |E₂|/C（C = max|Cz*/Cu*|）則 ratio < 1 ⟹ E₂ 項主導
"""
import itertools
import json

import numpy as np

N, MOD = 44, 90
G = np.array([4.0 * np.sin(np.pi * k / 90) ** 2 for k in range(0, N + 2)])
QP = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 5)), dtype=np.int32)
EQ = np.zeros(len(QP), dtype=float)
for t in range(5):
    EQ += G[QP[:, t]]

d = json.load(open("/home/wsl/idcm_ckpt/e8/mrB/rank1.json"))
rows = d["rows"]
res = []
for r in rows:
    n, j = r["n"], r["j"]
    gap = abs(float(G[n]) - float(EQ[j]))
    cu = abs(r["wj"]) * gap          # |Cu*| = |w*|·gap*
    cz = abs(r["Czj"])
    E2 = abs(r["E2"])
    ratio_obs = cz / (E2 * abs(r["wj"]))
    ratio_pred = (cz / cu) * gap / E2
    res.append(dict(n=n, j=j, gap=gap, cu=cu, cz=cz, E2=E2,
                    cz_over_cu=cz / cu, ratio_pred=ratio_pred, ratio_obs=ratio_obs))
Cmax = max(x["cz_over_cu"] for x in res)
dev = max(abs(x["ratio_pred"] - x["ratio_obs"]) / x["ratio_obs"] for x in res)
ok = sum(1 for x in res if x["gap"] < x["E2"] / Cmax)
print("【§DI】主導不等式之結構判據")
print(f"  恆等 ratio = (|Cz*|/|Cu*|)·gap*/|E₂|：最大相對差 {dev:.3e}")
print(f"  C = max|Cz*/Cu*| = {Cmax:.6f}")
print(f"  判據 gap* < |E₂|/C 命中 {ok}/44")
print(f"  界限 |E₂|/C：min {min(x['E2']/Cmax for x in res):.6e}｜gap* max {max(x['gap'] for x in res):.6e}")
print(f"  ratio_obs max {max(x['ratio_obs'] for x in res):.6f}｜上界 C·gap/E₂ max {max(x['gap']*Cmax/x['E2'] for x in res):.6f}")
bad = [x["n"] for x in res if not (x["gap"] < x["E2"] / Cmax)]
print(f"  未達判據之 n：{bad}")
print(f"  ratio<1 之 n 數：{sum(1 for x in res if x['ratio_obs']<1)}/44")
json.dump({"identity_max_rel_dev": dev, "C_max_cz_over_cu": Cmax, "criterion_hits": ok, "n": len(res),
           "bound_min": min(x["E2"] / Cmax for x in res), "gap_max": max(x["gap"] for x in res),
           "ratio_obs_max": max(x["ratio_obs"] for x in res),
           "bound_max": max(x["gap"] * Cmax / x["E2"] for x in res),
           "criterion_miss": bad, "ratio_lt1": sum(1 for x in res if x["ratio_obs"] < 1),
           "rows": res},
          open("/home/wsl/idcm_ckpt/e8/mrB/dominance.json", "w"), indent=1, ensure_ascii=False)
print("已寫 /home/wsl/idcm_ckpt/e8/mrB/dominance.json")
