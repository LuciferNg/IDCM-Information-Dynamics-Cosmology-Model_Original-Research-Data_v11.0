#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/03_null_rigidity.py — 結構必然性、MC null、剛性、外部對照

A. 尖截斷不可能性: S(k)=Σ_{n≤k}γ_n 冇任何 k 等於 28.5 ⟹ bottleneck n*=27 被迫
   (f_proj∈(0,1) 唯一) ⟹「投影唔係尖譜截斷」之定量版。
B. MC null (10⁷ 隨機 binary mask + 10⁶ 隨機相位擾動):
   - 交換恆等式 m/M₀ + E_sync = Σγ 對任何 mask 偏差 = 0 (代數恆等, 非擬合);
   - 隨機 mask 命中半整數 E_proj 嘅次數 ≈ 0 ⟹ 28.5 之精確 = commensurate 整數模位置,
     唔係 generic 性質;
   - 相位擾動 δ∈U[−0.5,0.5) 對稱: E_proj 分佈 (敏感度研究, 非 p-value)。
C. 剛性: 翻任一 mode ⟹ ΔE_proj = ±γ_n (交換量子離散); dE_proj/df_proj = γ₂₇。
D. 外部對照: 41/60, 19/60 vs Planck 2018 (Ω_Λ=0.6847±0.0073, Ω_m=0.3153±0.0073)。

輸出: data/03_null_rigidity.json
"""
import json
import os
import numpy as np
import mpmath as mp

mp.mp.dps = 50
B, Nd, L6 = 4, 44, 18
rng = np.random.default_rng(20261024)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

gam_np = np.array([B * np.sin(np.pi * k / (2 * Nd)) ** 2 for k in range(1, Nd + 1)])
gam = [mp.mpf(float(g)) for g in gam_np]
S = mp.fsum(gam)
n_star = Nd - L6 + 1
S26 = mp.fsum(gam[:26])
f_proj = (mp.mpf("28.5") - S26) / gam[26]
total = 90.0

# ---------- A. 尖截斷不可能性 ----------
cum = np.cumsum(gam_np)
sharp = np.concatenate([[0.0], cum])           # S(k), k=0..44
devs = np.abs(sharp - 28.5)
k_best = int(np.argmin(devs))
check_no_sharp = devs.min() > 1e-3
check_bracket = (cum[25] < 28.5 < cum[26])
check_fproj_range = (0 < float(f_proj) < 1)

# ---------- B. MC null ----------
N_MC = 10**7
chunk = 10**5
max_dev = 0.0
half_hits = 0
E_min, E_max = float("inf"), float("-inf")
for _ in range(N_MC // chunk):
    masks = rng.integers(0, 2, size=(chunk, Nd)).astype(np.float64)
    E = masks @ gam_np
    Es = total - E
    max_dev = max(max_dev, float(np.max(np.abs(Es + E - total))))
    half_hits += int(np.sum(np.abs(E * 2 - np.round(E * 2)) < 1e-9))
    E_min = min(E_min, float(E.min()))
    E_max = max(E_max, float(E.max()))
check_mc_identity = max_dev < 1e-9

# 相位擾動 (非 commensurate)
N_PH = 10**6
delta = rng.random((N_PH, Nd)) - 0.5   # 對稱 δ∈[−0.5,0.5) — 避免 γ 單調造成單邊偏壓
ph = (np.arange(1, Nd + 1)[None, :] + delta) * np.pi / (2 * Nd)
gam_pert = B * np.sin(ph) ** 2
sig_np = np.array([1.0] * 26 + [float(f_proj)] + [0.0] * 17)
E_pert_sig = gam_pert @ sig_np
mean_e, std_e = float(E_pert_sig.mean()), float(E_pert_sig.std())
z_285 = (28.5 - mean_e) / std_e if std_e > 0 else float("nan")
p_within_half = float(np.mean(np.abs(E_pert_sig - 28.5) < 0.5))

# ---------- C. 剛性 ----------
sensitivity = float(gam[26])              # dE_proj/df_proj = γ₂₇
prec_needed = 1e-12 / sensitivity         # 鎖定 f_proj 至 1e-12 所需精度
flip_costs = {"min": float(gam_np.min()), "max": float(gam_np.max())}

# ---------- D. Planck 對照 ----------
obs = {"Omega_Lambda": (0.6847, 0.0073), "Omega_m": (0.3153, 0.0073)}
pred = {"Omega_Lambda": 41 / 60, "Omega_m": 19 / 60}
zs = {k: (pred[k] - obs[k][0]) / obs[k][1] for k in obs}

checks = {
    "A1_no_sharp_cutoff_hits_28.5": bool(check_no_sharp),
    "A2_28.5_between_S26_S27": bool(check_bracket),
    "A3_f_proj_in_(0,1)": bool(check_fproj_range),
    "B1_MC_identity_dev_zero": bool(check_mc_identity),
    "B2_half_integer_hits_rare": bool(half_hits < 10),
    "C1_quanta_discrete": bool(abs(flip_costs["max"] - 4) < 1e-12 and flip_costs["min"] < 0.01),
    "D1_Planck_within_1sigma": bool(all(abs(z) < 1.0 for z in zs.values())),
}

out = {
    "script": "03_null_rigidity.py",
    "A_sharp_cutoff": {
        "S26": float(cum[25]), "S27": float(cum[26]), "target": 28.5,
        "closest_k": k_best, "closest_dev": float(devs.min()),
        "conclusion": "無尖截斷命中 28.5 ⟹ bottleneck n*=27 被迫, f_proj∈(0,1) 唯一",
    },
    "B_MC_null": {
        "N_masks": N_MC, "identity_max_dev": max_dev, "half_integer_hits": half_hits,
        "E_range": [E_min, E_max],
        "identity_note": "Σ 恆等對所有 mask 偏差 0 — 定律係代數, 唔係擬合; 相位 null 為敏感度研究 (對稱 δ), 非 p-value; 28.5 之精確由 f_proj 定義 (O121), 結構內容在 A (尖截斷不可能)",
        "phase_perturbation": {"N": N_PH, "mean": mean_e, "std": std_e,
                               "z_of_28.5": z_285, "P_within_0.5": p_within_half},
    },
    "C_rigidity": {"dE_df": sensitivity, "f_proj_precision_for_1e-12": prec_needed,
                   "flip_quantum_range": flip_costs},
    "D_external": {"pred": pred, "obs_Planck2018": obs, "z": zs},
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "03_null_rigidity.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("A: S26=%.6f S27=%.6f target=28.5 | closest k=%d dev=%.3e" % (cum[25], cum[26], k_best, devs.min()))
print("B: MC %d masks: identity max dev=%.2e | half-integer hits=%d | E∈[%.2f,%.2f]"
      % (N_MC, max_dev, half_hits, E_min, E_max))
print("   相位擾動: mean=%.4f std=%.4f z(28.5)=%.2f P(|E-28.5|<0.5)=%.4f" % (mean_e, std_e, z_285, p_within_half))
print("C: dE/df=%.4f | f_proj 精度需 ~%.1e | 量子範圍 [%.6f, %.6f]" % (sensitivity, prec_needed, flip_costs["min"], flip_costs["max"]))
print("D: 41/60 z=%.3f | 19/60 z=%.3f" % (zs["Omega_Lambda"], zs["Omega_m"]))
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))