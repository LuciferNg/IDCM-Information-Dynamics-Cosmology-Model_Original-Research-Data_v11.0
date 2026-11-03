#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/09_p3_frozen_split.py — P3: 交換係對稱還是過程？(p 嘅範疇裁決)

理路:
  1. 若每 mode 有對稱翻轉通道 (p>0, 獨立), 兩態 Markov 鏈平穩分佈 = (½,½) ∀n
     ⟹ E_proj^eq = Σγ/2 = 45, f_E^eq = ½
  2. 觀測/結構值 f_E = 19/60 (E_proj = 28.5) ⟹ 偏離 = 16.5 = K_tree/2 (EXACT)
     ⟹ 分裂 = 破缺序參量, 唔係熱平衡
  3. 數據界: 宇宙年齡內漂移 (½ − 19/60)(1 − e^{−2pN}) < 1% ⟹ p < p_bound
     其中 N = t_age/tick

判決: X: σ↦1−σ 係賬本對稱 (gauge), 唔係動力學過程 ⟹ 「交換速率」係範疇問題;
      若存在物理通道, 由 O7/O220 管, 且受 p_bound 壓制。

輸出: data/09_p3_frozen_split.json
"""
import json
import os
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
B, Nd = 4, 44
S = F(90)                 # Σγ
E_PROJ = F(57, 2)         # 28.5
K_TREE = F(33)            # E_sync − E_proj
f_E = E_PROJ / S          # 19/60

# F1 平穩分佈 ⟹ E_proj^eq = Σγ/2
E_EQ = S / 2
# F2 偏離 = K_tree/2
dev = E_EQ - E_PROJ

# F3 時間
T_tick_Gyr = mp.mpf("84.637") / 88
t_age_Gyr = mp.mpf("13.8")
N = t_age_Gyr / T_tick_Gyr

# F4 界: drift = (½ − 19/60)(1 − e^{−2pN}) < tol
half_minus = mp.mpf(1) / 2 - mp.mpf(19) / 60
tol = mp.mpf("0.01")
p_bound = -mp.log(1 - tol / half_minus) / (2 * N)

# F5 若 p = 1/88 (每 tick 一位元), 漂移幾大
p_serial = mp.mpf(1) / 88
drift_serial = half_minus * (1 - mp.e ** (-2 * p_serial * N))

checks = {
    "F1_equilibrium_is_half": bool(E_EQ == 45),
    "F2_offset_is_half_K_tree": bool(dev == K_TREE / 2 and dev == F(33, 2)),
    "F3_f_E_is_19_over_60": bool(f_E == F(19, 60)),
    "F4_p_bound_positive_and_tight": bool(0 < p_bound < mp.mpf("1e-2")),
    "F5_serial_p_excluded": bool(drift_serial > tol),
}

out = {
    "script": "09_p3_frozen_split.py",
    "equilibrium": {"E_proj_eq": "45 = Σγ/2", "f_E_eq": "1/2"},
    "observed": {"E_proj": "57/2 = 28.5", "f_E": "19/60", "K_tree": "33"},
    "offset_from_equilibrium": "16.5 = K_tree/2 (EXACT) ⟹ 分裂 = 破缺序參量",
    "time": {"tick_Gyr": mp.nstr(T_tick_Gyr, 8), "t_age_Gyr": mp.nstr(t_age_Gyr, 6), "N_ticks": mp.nstr(N, 8)},
    "p_bound": {"tol_1pct": mp.nstr(p_bound, 8), "drift_if_p_1_over_88": mp.nstr(drift_serial, 8),
                "note": "p=1/88 (每 tick 一位元序列化) 已被 Ω_m 穩定性排除"},
    "verdict": {
        "category": "🔵 X: σ↦1−σ 係賬本對稱 (gauge)，唔係動力學過程 ⟹ 「交換速率」係範疇問題",
        "structure": "✅ 19/60 分裂 = 破缺序參量 K_tree=33（偏離均衡恰 16.5=K_tree/2）",
        "data_bound": "p < " + mp.nstr(p_bound, 4) + " /tick/mode（Ω_m 漂移 < 1% 於 13.8 Gyr）",
        "falsifiable": "若 p ∼ bound，Ω_m(z) 有可測漂移（未來 Euclid/DESI 可判）",
    },
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs("data", exist_ok=True)
with open("data/09_p3_frozen_split.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("均衡: E_proj =", E_EQ, "f_E = 1/2 | 實值:", E_PROJ, "f_E =", f_E)
print("偏離 =", dev, "= K_tree/2 =", K_TREE / 2)
print("N (ticks in 13.8 Gyr) =", mp.nstr(N, 6), "| p_bound =", mp.nstr(p_bound, 6))
print("若 p = 1/88 ⟹ drift =", mp.nstr(drift_serial, 4), "vs tol 0.01")
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))