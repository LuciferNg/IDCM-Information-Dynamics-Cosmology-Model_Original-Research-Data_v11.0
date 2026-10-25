#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/07_phase_slip_link.py — P3 缺口: 互換量子 vs O7 寫入成本

問句 (01_problem_record §三 P3): ΔE = 2|h_i| (O7 寫入) 係唔係每 mode 互換嘅能量量子?

O7 事實 (docs/03 G1-G2, docs/04 §5):
  flip ⟺ 相位滑移 π ⟺ U₈₈ = e^{iπ} = −1
  h_i = A·cosθ_i, A = 134.2 μK (冷斑振幅); θ_i = arccos(T_i/A) = CMB 輸入 (5 假設否證)
  ΔE_write = B(n)·γ·|cosθ_i|, B(n) = 2E_min = 123/88, 44 寫 = 61.5
  ⚠ O7 docs/03 G2 明令: 兩個「2」(譜側 / 圖側淨+2) 屬平行層, 不強行統一

本 script 只測可證偽的結構身份, 零新常數, 全部精確算術。

輸出: data/07_phase_slip_link.json
"""
import json
import os
from fractions import Fraction as F
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
B, Nd = 4, 44
S_E = F(90)            # Σγ = B(N_d+1)/2 = 90 (O211)
E_SYNC = F(123, 2)     # 61.5 (O121)
E_PROJ = S_E - E_SYNC  # 28.5
E_MIN = E_SYNC / 88    # 61.5/88 = 123/176 (O41 每 bit 能量)
BN = 2 * E_MIN         # 123/88 (O41 每寫成本單位)
M0 = mp.mpf("27.9009145112")
kB_MeV_K = mp.mpf("8.617333262e-11")
A_uK = mp.mpf("134.2")

checks = {}

# D1 — Nyquist link ↔ 最高 rung: |1 − U₈₈|² = γ₄₄ = 4 (EXACT)
U88 = sp.exp(sp.I * sp.pi)
gam44 = 4 * sp.sin(sp.pi * 44 / 88) ** 2
checks["D1_U88_gives_top_rung_gamma44"] = bool(sp.simplify(sp.Abs(1 - U88) ** 2 - gam44) == 0 and gam44 == 4)

# D2 — O41 寫入預算 = M2 同步池 (EXACT)
checks["D2_44_writes_equal_E_sync"] = bool(44 * BN == E_SYNC)

# D3 — 總帳: E_sync + E_proj = Σγ (EXACT)
checks["D3_pools_sum_to_90"] = bool(E_SYNC + E_PROJ == S_E and E_PROJ == F(57, 2))

# D4 — 🔴 no-go: 123/88 唔等於任何 γ_n (代數整數論證, 非數值比較)
#   γ_n = 2 − 2cos(πn/44)。若 γ_n = 123/88 ⟹ 2cos(πn/44) = 2 − 123/88 = 53/88 (有理數)
#   但 2cos(π·有理) 係代數整數; 有理代數整數必為整數; 53/88 ∉ ℤ ⟹ 矛盾。
two_cos = 2 - BN
gamma_is_alg_int = True   # 2cos(2πk/m) 為代數整數 (標準定理)
rational_alg_int_is_integer = (two_cos.denominator == 1)
checks["D4_no_rung_equals_123_over_88"] = bool(not rational_alg_int_is_integer and two_cos == F(53, 88))
# 數值旁證: 解 4sin²(πn/88) = 123/88
n_sol = mp.asin(mp.sqrt(BN / 4)) * 88 / mp.pi
checks["D4b_solution_noninteger"] = bool(abs(n_sol - round(float(n_sol))) > mp.mpf("1e-3"))

# D5 — 均值比 = 90/61.5 = 60/41 = 1/Ω_Λ (EXACT)
ratio = (S_E / 44) / (E_SYNC / 44)
checks["D5_mean_ratio_60_over_41"] = bool(ratio == F(60, 41))

# D6 — 單位縫: ΔE = 2A = 268.4 μK 對應能量 vs 最輕 rung m₁
E_h = 2 * A_uK * mp.mpf("1e-6") * kB_MeV_K     # MeV
gam1 = 4 * mp.sin(mp.pi / 88) ** 2
m1 = M0 * gam1
unit_gap = m1 / E_h
checks["D6_unit_gap_reported"] = bool(unit_gap > mp.mpf("1e12"))

out = {
    "script": "07_phase_slip_link.py",
    "exact": {
        "gamma44": 4, "U88": "-1", "E_min": "123/176", "B(n)=2E_min": "123/88",
        "E_sync": "123/2", "E_proj": "57/2", "Sigma_gamma": "90",
        "mean_ratio": "60/41 = 1/Omega_Lambda",
        "two_cos_if_equal": "53/88 (非整數 ⟹ 矛盾)",
    },
    "unit_bridge": {
        "2A_uK": str(2 * A_uK), "2A_in_MeV": mp.nstr(E_h, 12),
        "m1_MeV": mp.nstr(m1, 12), "unit_gap_m1_over_2A": mp.nstr(unit_gap, 6),
        "n_solution_noninteger": mp.nstr(n_sol, 12),
    },
    "verdict": {
        "Z2_phase_slip_link": "✅ EXACT — γ₄₄ = |1−U₈₈|² = 4（相位滑移 π 坐喺質量階梯最高 rung）",
        "write_budget_link": "✅ EXACT — 44 寫 × 2E_min = 61.5 = E_sync（O7/O41 寫入預算 = M2 同步池）",
        "per_mode_quantum_identity": "🔴 NO-GO — 123/88 唔係任何 γ_n（代數整數論證），且 2|h_i| 屬 μK 讀數層（差 1e13）",
        "blocking_gap": "60 節點 ↔ 88-site 鏈嵌入（O7 docs/04 §5 唯一剩餘 🔴）；θ_i 係 CMB 輸入非圖輸出",
    },
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs("data", exist_ok=True)
with open("data/07_phase_slip_link.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("γ₄₄ = |1−U₈₈|² =", 4, "| B(n) =", BN, "| 44·B =", 44 * BN, "= E_sync")
print("γ_n = 123/88 ⟹ 2cos =", two_cos, "∉ ℤ ⟹ 矛盾（無任何 rung 等於 123/88）")
print("n 解 (數值旁證) =", mp.nstr(n_sol, 10), "（非整數）")
print("均值比 =", ratio, "= 1/Ω_Λ")
print("單位縫: 2A =", mp.nstr(E_h, 6), "MeV vs m₁ =", mp.nstr(m1, 6), "MeV ⟹", mp.nstr(unit_gap, 4), "×")
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))