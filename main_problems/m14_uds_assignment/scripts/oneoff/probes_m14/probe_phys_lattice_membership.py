#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""物理應用 3.2 之判準（1e−6 相對）：觀測質量比是否落在 {γ_b/γ_a} 之 946 個比值上？
   直接比較 m_hi/m_lo（p=1 讀法）與 γ_b/γ_a；報告最佳相對殘差
"""
import itertools
import json

import mpmath as mp

mp.mp.dps = 40
G = [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]
RAT = {(a, b): G[b - 1] / G[a - 1] for a in range(1, 45) for b in range(a + 1, 45)}
M = {"u": "0.00216", "d": "0.00467", "s": "0.0934", "c": "1.27", "b": "4.18", "t": "172.69",
     "e": "0.00051099895", "mu": "0.1056583755", "tau": "1.77686"}
PAIRS = [("t", "u"), ("tau", "e"), ("mu", "e"), ("tau", "mu"), ("b", "s"), ("c", "u"), ("t", "c"), ("b", "u")]
rows = []
for hi, lo in PAIRS:
    obs = mp.mpf(M[hi]) / mp.mpf(M[lo])
    best = min(((abs(g / obs - 1), (a, b)) for (a, b), g in RAT.items()), key=lambda t: t[0])
    within = obs <= max(RAT.values())
    rows.append({"pair": f"{hi}/{lo}", "obs": mp.nstr(obs, 12),
                 "best_rel_residual": mp.nstr(best[0], 6),
                 "best_stations": list(best[1]),
                 "within_range": bool(within),
                 "meets_1e-6": bool(within and best[0] < mp.mpf("1e-6"))})
    print(f"{hi:>4}/{lo:<4} obs {mp.nstr(obs,10):>16}｜最佳殘差 {mp.nstr(best[0],6):>12}｜站 {best[1]}｜在域內 {within}｜合格 {rows[-1]['meets_1e-6']}")
out = {"criterion": "|γ_b/γ_a − m_hi/m_lo| / (m_hi/m_lo) < 1e−6（p=1 讀法）",
       "n_ratios": len(RAT), "rows": rows,
       "n_pass": sum(r["meets_1e-6"] for r in rows),
       "verdict": ("判死：無任何觀測質量比落在 946 個 γ 比上（1e−6 級）"
                   if all(not r["meets_1e-6"] for r in rows) else "有命中，需逐一核")}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/phys_lattice_membership.json", "w"), indent=1, ensure_ascii=False)
print("判詞:", out["verdict"], "｜合格", out["n_pass"], "/", len(rows))
