#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""S2：④／⑤ 出口之可行性判定 —— 「頂點整數權重 ⇒ SM 4-點耦合」是否可作預言？
   方法：凑數陷阱偵測（而非凑數本身）
   ① 取 corpus 已註冊之結構整數池 P（不含任意整數）
   ② 枚舉 P×P 之比值 p/q，計有多少落在 λ_SM = 0.12907 之 1% 內
   ③ 若命中數 ≫ 1 ⟹ 任何「命中」皆無信息量 ⟹ 該路線不可作可證偽預言
"""
import itertools
import json

import mpmath as mp

mp.mp.dps = 30
POOL = {
    "N_d": 44, "Σγ": 88, "M90": 90, "60": 60, "κinv16": 16, "N_gen": 3, "B": 4,
    "Ω_m19": 19, "Ω_Λ41": 41, "45": 45, "22": 22, "12": 12, "36": 36, "98": 98,
    "421": 421, "3136": 3136, "sumw22230": 22230, "sigs667920": 667920,
    "w3_88": 88, "L6_18": 18, "L7_29": 29, "F5_5": 5, "L8_47": 47, "2Nd": 88 * 1,
    "N12_9": 9, "16b": 16, "17": 17, "11": 11, "8": 8, "24": 24, "4!": 24,
}
LAM_OBS = mp.mpf("0.12907")      # SM Higgs quartic λ = m_h²/(2v²)（v=246.2196 GeV, m_h=125.25）
TOL = mp.mpf("0.01")

hits = []
for (ka, a), (kb, b) in itertools.product(POOL.items(), repeat=2):
    if b == 0:
        continue
    r = mp.mpf(a) / mp.mpf(b)
    if r <= 0:
        continue
    if abs(r / LAM_OBS - 1) < TOL:
        hits.append({"ratio": f"{ka}/{kb}", "value": mp.nstr(r, 8),
                     "dev_pct": mp.nstr(100 * (r / LAM_OBS - 1), 4)})
hits = sorted(hits, key=lambda h: abs(mp.mpf(h["dev_pct"])))
print(f"λ_SM = {mp.nstr(LAM_OBS,7)}（1% 容差）｜結構整數池 {len(POOL)} 個｜比值命中 {len(hits)}")

# 容差掃描（陷阱起爆點）
sweep = {}
ratios = []
for (ka, a), (kb, b) in itertools.product(POOL.items(), repeat=2):
    if b:
        r = mp.mpf(a) / mp.mpf(b)
        if r > 0:
            ratios.append((f"{ka}/{kb}", r))
for tol_pct in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0):
    n = sum(1 for _, r in ratios if abs(r / LAM_OBS - 1) < mp.mpf(tol_pct) / 100)
    sweep[f"{tol_pct}%"] = n
    print(f"  容差 {tol_pct:>4}% ⟹ 命中 {n}")
for h in hits[:12]:
    print(f"  {h['ratio']:<22} = {h['value']:>12}  dev {h['dev_pct']:>7}%")
out = {"lambda_obs": mp.nstr(LAM_OBS, 8), "tol": "0.01", "pool_size": len(POOL),
       "tolerance_sweep": sweep, "n_ratio_pairs_distinct": len(ratios),
       "n_ratio_pairs": len(POOL) ** 2, "n_hits": len(hits), "hits": hits[:40],
       "verdict": ("凑數陷阱：命中數 ≫ 1 ⟹ 任何單一命中皆無信息量 ⟹ 該對照不可作可證偽預言"
                   if len(hits) > 3 else "命中稀少，需獨立算子推導核實")}
with open("/home/wsl/idcm_ckpt/e8/mrB/s2_lambda_trap.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("判詞:", out["verdict"])
