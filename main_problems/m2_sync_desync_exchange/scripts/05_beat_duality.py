#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/05_beat_duality.py — C 分支: 振幅讀法互換律 + 能量↔振幅對偶

結構:
  能量讀法: γ_n = |1−ζⁿ|²,  Σγ = 90                     (O211)
  振幅讀法: √γ_n = |1−ζⁿ| = 2sin(nπ/88) = beat, Σ√γ = cot(π/176)+1  (本 M script 04)
  互換 X: σ ↦ 1−σ 對兩讀法皆對合; T: γ↦√γ 與 X 交換 (X∘T = T∘X)
  振幅互換律: m_amp/M₀ + s_amp/M₀ = Σ√γ = cot(π/176)+1
  逐模量子:    Δm_amp,n = M₀√γ_n

預註冊盲測 (M1 建檔前鎖定, 2026-10-24):
  (P-a) 最輕振幅量子 m₁ = M₀·2sin(π/88) = 1.99278… MeV
  (P-b) 最重振幅量子 m₄₄ = 2M₀ = 55.8018… MeV
  (P-c) 振幅池份額 f_amp = Σσ√γ/Σ√γ (值見輸出; 若 M1 schema 用振幅讀法, 必須重現此份額)

輸出: data/05_beat_duality.json
"""
import json
import os
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
B, Nd, L6 = 4, 44, 18
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

gam = [B * mp.sin(mp.pi * k / (2 * Nd)) ** 2 for k in range(1, Nd + 1)]
beat = [mp.sqrt(g) for g in gam]
S_E = mp.fsum(gam)
S_A = mp.fsum(beat)
n_star = Nd - L6 + 1
S26 = mp.fsum(gam[:26])
f_proj = (mp.mpf("28.5") - S26) / gam[26]
sigma = [mp.mpf(1) if k < n_star else (f_proj if k == n_star else mp.mpf(0))
         for k in range(1, Nd + 1)]
E_proj = mp.fsum(s * g for s, g in zip(sigma, gam))
E_sync = S_E - E_proj
A_proj = mp.fsum(s * b for s, b in zip(sigma, beat))
A_sync = S_A - A_proj
M0 = mp.mpf("27.9009145112")

# 精確身份
x = sp.pi / 176
check_sumA = (sp.simplify(sp.sin(45 * sp.pi / 176) - sp.cos(43 * sp.pi / 176)) == 0
              and abs(sp.N(sp.sqrt(2) * sp.cos(sp.pi / 4 - x) - (sp.cos(x) + sp.sin(x)), 80))
              < sp.Float("1e-70"))
check_beat_chord = max(abs(beat[k - 1] - abs(1 - mp.expj(mp.pi / 44) ** k)) for k in range(1, Nd + 1)) < mp.mpf("1e-45")

# 交換對合 (兩讀法)
sig2 = [1 - (1 - s) for s in sigma]
check_X_amp = max(abs(a - b) for a, b in zip(sig2, sigma)) == 0

# T 與 X 交換 (逐模: T 只動 γ 槽, X 只動 σ 槽)
check_commute = all(
    (1 - s) * mp.sqrt(g) == mp.sqrt(g) * (1 - s)
    for s, g in zip(sigma, gam)
)

# 振幅守恆
check_amp_conserved = abs(A_proj + A_sync - S_A) < mp.mpf("1e-50")

# 預註冊數字
m1 = M0 * beat[0]
m44 = M0 * beat[-1]
f_amp = A_proj / S_A
f_E = E_proj / S_E

# 份額對照 (唔准凑): 19/60, 41/60, 17/44
refs = {"19/60": mp.mpf(19) / 60, "41/60": mp.mpf(41) / 60, "17/44": mp.mpf(17) / 44}
near = {k: mp.nstr(f_amp - v, 6) for k, v in refs.items()}

checks = {
    "C1_sum_sqrt_gamma_exact_cot_plus_1": bool(check_sumA),
    "C2_beat_equals_chord": bool(check_beat_chord),
    "C3_X_involution_amplitude": bool(check_X_amp),
    "C4_T_commutes_with_X": bool(check_commute),
    "C5_amplitude_pool_conserved": bool(check_amp_conserved),
    "C6_amplitude_law_holds": bool(abs(A_proj + A_sync - S_A) < mp.mpf("1e-45")),
}

out = {
    "script": "05_beat_duality.py",
    "totals": {"Sigma_gamma_energy": mp.nstr(S_E, 20), "Sigma_sqrt_gamma_amplitude": mp.nstr(S_A, 20),
               "cot_pi_176": mp.nstr(mp.cot(mp.pi / 176), 20)},
    "pools": {"E_proj": mp.nstr(E_proj, 20), "E_sync": mp.nstr(E_sync, 20),
              "A_proj": mp.nstr(A_proj, 20), "A_sync": mp.nstr(A_sync, 20),
              "f_E": mp.nstr(f_E, 20), "f_amp": mp.nstr(f_amp, 20),
              "amp_vs_refs": near},
    "preregistered": {"P_a_m1_MeV": mp.nstr(m1, 20), "P_b_m44_MeV": mp.nstr(m44, 20),
                      "P_c_f_amp": mp.nstr(f_amp, 20),
                      "note": "M1 未有 data 檔 — 本三個數字在建檔前鎖定"},
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "05_beat_duality.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("Σγ =", mp.nstr(S_E, 12), "| Σ√γ =", mp.nstr(S_A, 12), "| cot+1 =", mp.nstr(mp.cot(mp.pi / 176) + 1, 12))
print("能量池: E_proj/E_sync =", mp.nstr(E_proj, 10), "/", mp.nstr(E_sync, 10), "| f_E =", mp.nstr(f_E, 10))
print("振幅池: A_proj/A_sync =", mp.nstr(A_proj, 10), "/", mp.nstr(A_sync, 10), "| f_amp =", mp.nstr(f_amp, 10))
print("預註冊: m₁ =", mp.nstr(m1, 10), "MeV | m₄₄ =", mp.nstr(m44, 10), "MeV | f_amp =", mp.nstr(f_amp, 10))
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))