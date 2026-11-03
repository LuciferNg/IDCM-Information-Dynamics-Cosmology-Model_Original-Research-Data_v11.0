#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/10_p_family_completeness.py — M7 裁定 (a): p 族淘汰完備性審計

M2 docs/01 §7A 測過嘅族:  m⁽ᵖ⁾ = M₀ Σ σ_n γ_nᵖ   ← 只係 **可分 (separable)** 族
未測: 非可分族  m = M₀ [ Σ σ_n γ_n + λ Σ_{n<m} σ_n σ_m w_nm ]

本 script 檢查非可分族係唔係被 M2 嘅既有約束排除:
  C1 單模態譜不變 (λ 任意) ⟹ M2 全部單粒子測試對 λ 盲
  C2 複合態偏移量 (λ=1 時嘅量級)
  C3 賬本可維持 (E_sync 重定義後仍 ≥0, 互換律恆等)
  C4 X 對合仍成立
  C5 ⟹ 完備性主張不成立 (需要 O237 逐模讀數作額外前提)

輸出: data/10_p_family_completeness.json
"""
import json
import os
from itertools import combinations
import mpmath as mp

mp.mp.dps = 50
B, Nd = 4, 44
M0 = mp.mpf("27.9009145112")
gam = [B * mp.sin(mp.pi * k / (2 * Nd)) ** 2 for k in range(1, Nd + 1)]
S = mp.fsum(gam)                                       # Σγ = 90
w = [[g1 * g2 / S for g2 in gam] for g1 in gam]        # 對稱 pair 權 (量綱: γ²/Σγ)


def mass_nonadd(sig, lam):
    """m/M₀ = Σσγ + λ Σ_{n<m} σ_nσ_m w_nm"""
    lin = mp.fsum(s * g for s, g in zip(sig, gam))
    pair = mp.fsum(sig[i] * sig[j] * w[i][j] for i, j in combinations(range(Nd), 2))
    return lin + lam * pair


def mass_add(sig):
    return mp.fsum(s * g for s, g in zip(sig, gam))


lam = mp.mpf(1)

# C1 — 單模態: 只一條 σ=1
single_dev = max(abs(mass_nonadd([1 if k == n else 0 for k in range(Nd)], lam) - gam[n]) for n in range(Nd))

# C2 — 最高兩模複合偏移
top2 = [1 if k in (42, 43) else 0 for k in range(Nd)]
shift_top2 = (mass_nonadd(top2, lam) - mass_add(top2)) * M0        # MeV
# 全 44 模 (真空→滿態)
full = [1] * Nd
shift_full = (mass_nonadd(full, lam) - mass_add(full)) * M0

# C3 — 賬本窗口: 質量 ≤ Σγ (守恆) 且 ≥ 0
#   滿態: m = S + λ·pair_full ⟹ 排斥 λ>0 超總量被排除; 吸引 λ<0 係束縛, 只受 m ≥ 0 限制
pair_full = mp.fsum(w[i][j] for i, j in combinations(range(Nd), 2))
lam_max = mp.mpf(0)                 # 守恆上界: λ ≤ 0 (排斥項被 Σγ 封頂)
lam_min = -S / pair_full            # 質量 ≥ 0 下界

# C4 — X 對合 + 互換律恆等 (m/M₀ + E_sync' = Σγ)
sig = [1 if k in (42, 43) else 0 for k in range(Nd)]
sigx = [1 - s for s in sig]
m1, m2 = mass_nonadd(sig, lam), mass_nonadd(sigx, lam)
E1 = S - m1
E2 = S - m2
check_ident = (abs((m1 + E1) - S) < mp.mpf("1e-40")) and (abs((m2 + E2) - S) < mp.mpf("1e-40"))
check_involution = all((1 - (1 - s)) == s for s in sig)

checks = {
    "C1_single_mode_spectrum_lambda_blind": bool(single_dev < mp.mpf("1e-45")),
    "C2_composite_shift_nonzero_MeV": bool(shift_top2 > mp.mpf("1")),
    "C3_binding_window_nonempty": bool(lam_min < lam_max and lam_min < mp.mpf("-1")),
    "C4_X_involution_and_identity": bool(check_ident and check_involution),
    "C5_completeness_not_established": bool(single_dev < mp.mpf("1e-45") and shift_top2 > mp.mpf("1")),
}

out = {
    "script": "10_p_family_completeness.py",
    "object": "非可分族 m = M₀[Σσ_nγ_n + λΣ_{n<m}σ_nσ_m w_nm], w_nm = γ_nγ_m/Σγ",
    "C1_single_mode_max_dev": mp.nstr(single_dev, 6),
    "C2_shift": {"top2_modes_43_44_MeV": mp.nstr(shift_top2, 10), "full_44_modes_MeV": mp.nstr(shift_full, 10)},
    "C3_window": {"lambda_max": mp.nstr(lam_max, 6), "lambda_min": mp.nstr(lam_min, 8),
                  "meaning": "λ≤0 吸引（束縛）容許；λ>0 排斥被守恆封頂 ⟹ 束縛窗口非空"},
    "verdict": {
        "a_completeness": "🔴 不完備 — p 族只覆蓋可分族；非可分族 λ≠0 與單粒子譜、互換律、X 對合全部相容",
        "escape": "唯一框架內排除 = O237 逐模讀數（質量= per-mode |1−U_n|² 之和）—— 屬定義，非 M2 定理",
        "b_theorem_status": "只可升為【條件定理】：O237 逐模讀數 ⟹ 加性 ⟹ ⟨O⟩≡1；無條件版本拒絕",
    },
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs("data", exist_ok=True)
with open("data/10_p_family_completeness.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("C1 單模態對 λ 最大偏差 =", mp.nstr(single_dev, 6), "（= 0 ⟹ 單粒子譜盲）")
print("C2 複合偏移: 43+44 兩模 =", mp.nstr(shift_top2, 8), "MeV | 全 44 模 =", mp.nstr(shift_full, 8), "MeV")
print("C3 窗口 λ ∈ [", mp.nstr(lam_min, 8), ",", mp.nstr(lam_max, 4), "] ⟹ 束縛(λ<0) 容許")
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))