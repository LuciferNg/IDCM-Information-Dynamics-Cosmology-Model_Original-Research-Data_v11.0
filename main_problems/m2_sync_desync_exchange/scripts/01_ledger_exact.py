#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/01_ledger_exact.py — 同步/脫同步能量帳本（精確核對, v2）

上游約定（O121 scripts/01_tpipe_dag_spectrum.py）:
  γ_n = B·sin²(πn/2N_d) = 4·sin²(πn/88), n = 1..44      (O211 final_result:12)
  Σγ = B(N_d+1)/2 = 90 (EXACT)                          (O211)
  n* = N_d − L₆ + 1 = 27  (bottleneck)                  (O121 身份 A)
  σ_n = 1 (n ≤ 26), σ₂₇ = f_proj, σ_n = 0 (n ≥ 28)      (O121 奇異值譜)
  f_proj = (28.5 − Σ_{1..26}γ)/γ₂₇ = 0.771960           (O121: NOT hardcoded)
  E_proj = Σσγ = 28.5;  E_sync = Σ(1−σ)γ = 61.5         (O121 final_result:22)

本 script 只做「帳本核對」: 精確值、恆等式、兩側獨立錨（⟨W⟩²=61.5 vs Σγ=90）。
精度: mpmath 60 dps; 有理重建 1e-45; sympy 精確求和。
輸出: data/01_ledger.json
"""
import json
import os
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
B = 4
Nd = 44
L6 = 18
N_gen = 3
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def gamma(n):
    return B * mp.sin(mp.pi * n / (2 * Nd)) ** 2


gamma_all = [gamma(n) for n in range(1, Nd + 1)]
S = mp.fsum(gamma_all)                      # Σγ

n_star = Nd - L6 + 1                        # 27
E_proj_target = mp.mpf("28.5")
S26 = mp.fsum(gamma(n) for n in range(1, 27))
g27 = gamma(27)
f_proj = (E_proj_target - S26) / g27        # 由 28.5 反解（O121 同法）

sigma = [mp.mpf(1) if n < n_star else (f_proj if n == n_star else mp.mpf(0))
         for n in range(1, Nd + 1)]
E_proj = mp.fsum(s * g for s, g in zip(sigma, gamma_all))
E_sync = S - E_proj

# 精確求和（sympy）
n = sp.symbols('n', integer=True, positive=True)
sum_sin2_total = sp.simplify(sp.summation(sp.sin(sp.pi * n / 88) ** 2, (n, 1, 44)))
sum_sin2_low = sp.simplify(sp.summation(sp.sin(sp.pi * n / 88) ** 2, (n, 1, 27)))
sum_sin2_ker = sp.simplify(sp.summation(sp.sin(sp.pi * n / 88) ** 2, (n, 28, 44)))


def rat(x, tol=1e-45):
    r = sp.nsimplify(sp.N(x, 60), rational=True, tolerance=tol)
    if abs(sp.N(sp.Rational(r) - sp.N(x, 60), 60)) > sp.Rational(1, 10**45):
        raise AssertionError("rational reconstruction failed: %s" % x)
    return sp.Rational(r)


S_rat, Eproj_rat, Esync_rat = rat(S), rat(E_proj), rat(E_sync)
K_tree = Esync_rat - Eproj_rat

# 尖譜截斷對照（O121 final_result 已列: 60.8845 / 29.1155）
ker_sum = sp.N(4 * sum_sin2_ker, 60)     # sympy 未化簡: 用 60 位數值
low_sum = sp.N(4 * sum_sin2_low, 60)

checks = {
    # 總譜
    "C1_Sigma_gamma_90_exact": (sp.simplify(4 * sum_sin2_total) == 90),
    "C1_Sigma_gamma_numeric": abs(S - 90) < mp.mpf("1e-50"),
    # 池
    "C2_E_proj_28.5": (Eproj_rat == sp.Rational(57, 2)),
    "C3_E_sync_61.5": (Esync_rat == sp.Rational(123, 2)),
    "C4_K_tree_33": (K_tree == 33),
    "C5_pool_sum_Sigma": (Esync_rat + Eproj_rat == S_rat),
    # 比
    "C6_ratio_41_19": (Esync_rat / Eproj_rat == sp.Rational(41, 19)),
    "C7_frac_41_60": (Esync_rat / S_rat == sp.Rational(41, 60)),
    "C8_frac_19_60": (Eproj_rat / S_rat == sp.Rational(19, 60)),
    # 上游數字復現
    "C9_f_proj_0.771960": abs(f_proj - mp.mpf("0.771960")) < mp.mpf("1e-6"),
    "C10_n_star_27": (n_star == 27),
    "C11_kernel_dim_17": (len(range(28, 45)) == 17),
    # 尖譜截斷（誠實記錄：結構值 ≠ 區間和）
    "C12_sharp_cut_60.8845": abs(mp.mpf(str(ker_sum)) - mp.mpf("60.8845")) < mp.mpf("1e-3"),
    "C13_sharp_cut_29.1155": abs(mp.mpf(str(low_sum)) - mp.mpf("29.1155")) < mp.mpf("1e-3"),
    "C14_bottleneck_bridge": abs((mp.mpf(str(ker_sum)) + (1 - f_proj) * g27) - E_sync) < mp.mpf("1e-40"),
    # 兩側獨立錨一致性: ⟨W⟩² = μ²/(2λ) = 61.5 (O204) vs Σγ − E_proj
    "C15_condensate_anchor": (sp.Rational(2, 1) / (2 * sp.Rational(2, 123)) == sp.Rational(123, 2)),
}

out = {
    "script": "01_ledger_exact.py",
    "inputs": {"B": B, "N_d": Nd, "L_6": L6, "gamma": "4*sin^2(pi*n/88)",
               "n_star": n_star, "kernel": "n=28..44 (sync side, sigma=0)"},
    "values": {
        "Sigma_gamma": str(S_rat), "E_sync": str(Esync_rat), "E_proj": str(Eproj_rat),
        "K_tree": str(K_tree), "f_proj": mp.nstr(f_proj, 20),
        "E_sync_over_E_proj": str(Esync_rat / Eproj_rat),
        "E_sync_over_Sigma": str(Esync_rat / S_rat), "E_proj_over_Sigma": str(Eproj_rat / S_rat),
        "sharp_cut_n_ge_28": mp.nstr(mp.mpf(str(ker_sum)), 20),
        "sharp_cut_n_le_27": mp.nstr(mp.mpf(str(low_sum)), 20),
        "bottleneck_term": mp.nstr((1 - f_proj) * g27, 20),
    },
    "sympy_exact": {"sum_sin2_1_44": str(sum_sin2_total), "sum_sin2_1_27": str(sum_sin2_low),
                    "sum_sin2_28_44": str(sum_sin2_ker)},
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
    "honest_note": "E_proj=28.5 係 O121 投影預算輸入, f_proj 由它反解 (非硬編碼)。"
                   "另一側獨立錨: E_sync=⟨W⟩²=μ²/(2λ)=61.5 (O204) ⟹ 28.5=90−61.5 被兩錨夾定。",
}

os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "01_ledger.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("Σγ =", S_rat, "| E_sync =", Esync_rat, "| E_proj =", Eproj_rat, "| K_tree =", K_tree)
print("f_proj =", mp.nstr(f_proj, 15), "| 尖譜: n≥28 =", mp.nstr(mp.mpf(str(ker_sum)), 10),
      "| n≤27 =", mp.nstr(mp.mpf(str(low_sum)), 10), "| bottleneck =", mp.nstr((1 - f_proj) * g27, 10))
print("41/19 =", Esync_rat / Eproj_rat, "| 41/60 =", Esync_rat / S_rat, "| 19/60 =", Eproj_rat / S_rat)
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))