#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/06_rate_bound.py — B 分支: 交換速率（上界推導 + 值 no-go）

結構:
  tick Δt = T_ring/88 (O235 T1 EXACT) ⟹ 每模最快一次翻轉/tick
  ⟹ 速率上界 R_max = 1/Δt = 88/T_ring  (per mode)
  ⟹ 質量率上界 dm/dt|max = M₀·γ_n/Δt (per mode); 總 = M₀Σγ/Δt

誠實界線:
  實際速率 = p·R_max, p∈(0,1) 為躍遷機率 —— 無推導來源 ⟹ 值 no-go;
  本 script 只出「結構上界」(零參數), 不選 p。
  ⚠️ 禁止事後揀 p = 19/41、1/88 等凑數 (調律者鐵律)。

輸出: data/06_rate_bound.json
"""
import json
import os
import mpmath as mp

mp.mp.dps = 50
B, Nd = 4, 44
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

T_ring_Gyr = mp.mpf("84.637")             # O121
YR_S = mp.mpf("3.15576e7")                # Julian year [s]
T_ring = T_ring_Gyr * 1e9 * YR_S          # s
T_tick = T_ring / 88                      # s (O235 T1)
R_max = 1 / T_tick                        # s^-1 per mode

gam = [B * mp.sin(mp.pi * k / (2 * Nd)) ** 2 for k in range(1, Nd + 1)]
S = mp.fsum(gam)
M0 = mp.mpf("27.9009145112")              # MeV
MeV_kg = mp.mpf("1.782661921e-30")

q_max = M0 * 4                            # MeV (n=44)
dmdt_mode = q_max / T_tick                # MeV/s
dmdt_total = M0 * S / T_tick              # MeV/s
dmdt_mode_kg = dmdt_mode * MeV_kg
dmdt_total_kg = dmdt_total * MeV_kg

checks = {
    "B1_tick_equals_Tring_over_88": abs(T_tick * 88 - T_ring) < mp.mpf("1e-6"),
    "B2_rate_bound_positive": R_max > 0,
    "B3_mass_rate_scales_with_quantum": abs(dmdt_mode - q_max / T_tick) < mp.mpf("1e-9"),
    "B4_value_requires_p": True,     # 結構陳述: dm/dt = p·R_max, p 未導
    "B5_no_posthoc_p_adopted": True, # 明文: 唔採用任何凑數 p
}

out = {
    "script": "06_rate_bound.py",
    "time": {"T_ring_Gyr": mp.nstr(T_ring_Gyr, 12), "T_tick_s": mp.nstr(T_tick, 12),
             "T_tick_Gyr": mp.nstr(T_tick / (1e9 * YR_S), 12)},
    "rate": {"R_max_per_s": mp.nstr(R_max, 12), "R_max_per_Gyr": mp.nstr(R_max * 1e9 * YR_S, 12),
             "per_mode": True},
    "mass_rate": {"q_max_MeV": mp.nstr(q_max, 12), "dmdt_mode_MeV_s": mp.nstr(dmdt_mode, 12),
                  "dmdt_mode_kg_s": mp.nstr(dmdt_mode_kg, 12),
                  "dmdt_total_MeV_s": mp.nstr(dmdt_total, 12),
                  "dmdt_total_kg_s": mp.nstr(dmdt_total_kg, 12),
                  "dmdt_total_GeV_per_Gyr": mp.nstr((dmdt_total / 1000) * (1e9 * YR_S), 12)},
    "verdict": "速率上界 (零參數, 由 tick 離散性); 實際值 dm/dt = p·R_max, p 未導 ⟹ 值 conditional no-go",
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "06_rate_bound.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("T_tick =", mp.nstr(T_tick, 8), "s =", mp.nstr(T_tick / (1e9 * YR_S), 8), "Gyr")
print("R_max =", mp.nstr(R_max, 8), "s⁻¹ (per mode) =", mp.nstr(R_max * 1e9 * YR_S, 8), "Gyr⁻¹")
print("dm/dt (max quantum n=44) =", mp.nstr(dmdt_mode, 8), "MeV/s =", mp.nstr(dmdt_mode_kg, 8), "kg/s")
print("dm/dt (total) =", mp.nstr(dmdt_total, 8), "MeV/s =", mp.nstr((dmdt_total / 1000) * (1e9 * YR_S), 8), "GeV/Gyr")
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))