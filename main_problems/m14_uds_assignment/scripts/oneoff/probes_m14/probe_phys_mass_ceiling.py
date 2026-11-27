#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""物理應用第一級：質量比可否由 /90 譜 γ 比承載？（天花板測試，零擬合）
   判準：m ∝ γ^p ⟹ 可承載之最大質量比 = (γ_44/γ_1)^p
   若 (γ_44/γ_1)^p < 觀測比 ⟹ 該通道對任何標記指派皆不可能（不需味站指派）
"""
import json

import mpmath as mp

mp.mp.dps = 50
G1 = 4 * mp.sin(mp.pi * mp.mpf(1) / 90) ** 2
G44 = 4 * mp.sin(mp.pi * mp.mpf(44) / 90) ** 2
R = G44 / G1

# PDG 2024（MS̄ 2 GeV for u,d,s,c,b；t 用 pole；帶電輕子 pole）
M = {
    "u": mp.mpf("0.00216"), "d": mp.mpf("0.00467"), "s": mp.mpf("0.0934"),
    "c": mp.mpf("1.27"), "b": mp.mpf("4.18"), "t": mp.mpf("172.69"),
    "e": mp.mpf("0.00051099895000"), "mu": mp.mpf("0.1056583755"), "tau": mp.mpf("1.77686"),
}
pairs = [("t", "u"), ("tau", "e"), ("mu", "e"), ("tau", "mu"), ("b", "s"), ("c", "u"), ("t", "c"), ("b", "u")]
rows = []
for a, b in pairs:
    obs = M[a] / M[b]
    need_p = mp.log(obs) / mp.log(R)
    rows.append({"a": a, "b": b, "obs_ratio": mp.nstr(obs, 12),
                 "need_power_p": mp.nstr(need_p, 10),
                 "p_half_ceiling": mp.nstr(mp.sqrt(R), 10),
                 "p_one_ceiling": mp.nstr(R, 10),
                 "feasible_p_half": bool(mp.sqrt(R) >= obs),
                 "feasible_p_one": bool(R >= obs)})

out = {"gamma_1": mp.nstr(G1, 20), "gamma_44": mp.nstr(G44, 20),
       "gamma_ratio": mp.nstr(R, 20),
       "mass_ceiling_p_half": mp.nstr(mp.sqrt(R), 12),
       "mass_ceiling_p_one": mp.nstr(R, 12),
       "required_power_min_over_all_pairs": mp.nstr(max(mp.log(M[a] / M[b]) / mp.log(R) for a, b in pairs), 10),
       "rows": rows,
       "verdict": "純 γ 比不可能承載觀測質量比（任何標記指派皆不可能）" if not any(r["feasible_p_one"] for r in rows) else "部分可能"}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/phys_mass_ceiling.json", "w"), indent=1, ensure_ascii=False)
print(f"γ_1 = {mp.nstr(G1,12)}｜γ_44 = {mp.nstr(G44,12)}｜γ_44/γ_1 = {mp.nstr(R,12)}")
print(f"m ∝ γ^(1/2) 之質量比天花板 = {mp.nstr(mp.sqrt(R),8)}")
print(f"m ∝ γ^1   之質量比天花板 = {mp.nstr(R,8)}")
print(f"全部配對所需最小指數 p_min = {out['required_power_min_over_all_pairs']}")
for r in rows:
    print(f"{r['a']:>4}/{r['b']:<4} 觀測 {r['obs_ratio']:>16}  需 p {r['need_power_p']:>10}  p=1 可行 {r['feasible_p_one']}")
print("判詞:", out["verdict"])
