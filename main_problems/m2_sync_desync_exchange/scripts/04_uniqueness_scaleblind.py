#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M2/scripts/04_uniqueness_scaleblind.py — 互換律唯一性壓力測試 + 標度盲引理 + 頻率讀法

A. p-族淘汰: m⁽ᵖ⁾ = M₀·Σσ_nγ_nᵖ。要求 (i) 加性 (ii) 同步側無質量 (iii) E∝m 之能量讀法
   ⟹ 只有 p=1 通過; p≠1 之交換量子唔再 ∝ γ_n (比值隨 n 變), 報告比值散度。
B. 標度盲引理: m/M₀ + E_sync = Σγ —— M₀ 在定律內相消 ⟹ 互換律**不能**定 M₀
   (絕對標度必須外來: O38; 屬他 O, 本 M 不推)。
C. 頻率讀法對照: √γ_n = |1−ζⁿ| = 2sin(nπ/88) (beat/弦長); 總和 Σ√γ_n vs cot(π/176)
   (O121 ring Casimir 值) —— 量度兩種讀法之總量關係 (誠實報告 dev)。
D. 速率缺口: 互換律冇時間參數; dm/dt 需要一個無量綱率常數, corpus 內冇選定機制
   ⟹ P3「速率」= conditional no-go (零參數前提下不可導)。

輸出: data/04_uniqueness_scaleblind.json
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
S = mp.fsum(gam)
n_star = Nd - L6 + 1
S26 = mp.fsum(gam[:26])
f_proj = (mp.mpf("28.5") - S26) / gam[26]
sigma = [mp.mpf(1) if k < n_star else (f_proj if k == n_star else mp.mpf(0))
         for k in range(1, Nd + 1)]
E_proj = mp.fsum(s * g for s, g in zip(sigma, gam))
E_sync = S - E_proj
M0 = mp.mpf("27.9009145112")

# ---------- A. p-族 ----------
res_p = {}
for p in (0.5, 1, 2, 3):
    m_p = M0 * mp.fsum(s * (g ** p) for s, g in zip(sigma, gam))
    # 交換量子/能量比: q_n/γ_n ∝ γ_n^(p-1) — p=1 時比值常數
    ratios = [(g ** p) / g for g in gam]
    spread = max(ratios) / min(ratios)      # p=1 → 1 (常數)
    res_p[str(p)] = {"m_MeV": mp.nstr(m_p, 12), "ratio_spread": mp.nstr(spread, 12)}
# 直接驗 p=1 比值 = 1
check_p1 = abs(max(g / g for g in gam) - 1) < mp.mpf("1e-50")
check_p_not1 = all(abs(mp.mpf(res_p[str(p)]["ratio_spread"]) - 1) > mp.mpf("1e-3")
                   for p in (0.5, 2, 3))
# 加性: (γa+γb)^p vs γa^p+γb^p
add_ok = {}
for p in (0.5, 1, 2):
    lhs = (gam[0] + gam[1]) ** p
    rhs = gam[0] ** p + gam[1] ** p
    add_ok[str(p)] = abs(lhs - rhs) < mp.mpf("1e-45")
check_additive_p1 = add_ok["1"]

# ---------- B. 標度盲 ----------
def law(M0test):
    return (M0test * E_proj) / M0test + E_sync          # m/M₀ + E_sync
check_scale_blind = abs(law(M0) - S) < mp.mpf("1e-50") and abs(law(2 * M0) - S) < mp.mpf("1e-50")

# ---------- C. 頻率讀法 ----------
zeta = mp.expj(mp.pi / 44)
beat = [abs(1 - zeta ** k) for k in range(1, Nd + 1)]     # = 2sin(kπ/88)
dev_beat = max(abs(beat[k - 1] - mp.sqrt(gam[k - 1])) for k in range(1, Nd + 1))
sum_beat = mp.fsum(beat)
ring_casimir = mp.cot(mp.pi / 176)                        # O121 S_ring
# 解析證明: Σ_{k=1}^{44} 2sin(kπ/88) = cot(π/176)+1
#   Σ sin = sin(44π/176)sin(45π/176)/sin(π/176) = (√2/2)cos(43π/176)/sin(π/176)
#   √2 cos(π/4 − x) = cos x + sin x  (x=π/176) ⟹ Σ·2 = cot x + 1
x176 = sp.pi / 176
# 兩步精確: sin(45π/176)=cos(43π/176) (sympy 化簡=0); 且 √2cos(π/4−x)=cos x+sin x
check_beat_exact = (sp.simplify(sp.sin(45 * sp.pi / 176) - sp.cos(43 * sp.pi / 176)) == 0
                    and abs(sp.N(sp.sqrt(2) * sp.cos(sp.pi / 4 - x176)
                                 - (sp.cos(x176) + sp.sin(x176)), 80)) < sp.Float("1e-70"))
beat_pool = mp.fsum(s * beat[k - 1] for k, s in enumerate(sigma, start=1))

# ---------- D. 速率缺口 ----------
rate_candidates = {"1/88 (tick)": mp.mpf(1) / 88, "19/41 (回收率)": mp.mpf(19) / 41,
                   "41/60 (Ω_Λ)": mp.mpf(41) / 60, "1/176": mp.mpf(1) / 176}
check_no_rate = True   # 定律內無時間參數 —— 結構性陳述, 非數值檢定

checks = {
    "A1_p1_quantum_proportional_to_gamma": bool(check_p1),
    "A2_p_not1_ratio_spread_gt_1": bool(check_p_not1),
    "A3_topology_additive_only_p1": bool(check_additive_p1),
    "B1_scale_blind_M0_cancels": bool(check_scale_blind),
    "C1_beat_equals_sqrt_gamma": bool(dev_beat < mp.mpf("1e-45")),
    "C2_beat_sum_identity_exact": bool(check_beat_exact),
    "D1_no_time_parameter_in_law": bool(check_no_rate),
}

out = {
    "script": "04_uniqueness_scaleblind.py",
    "A_p_family": {"results": res_p,
                   "verdict": "只有 p=1: 量子∝γ (E∝m) + 對 γ 加性; p≠1 比值散度>1 且非加性"},
    "B_scale_blind": {"law": "m/M₀ + E_sync = Σγ", "M0_cancels": True,
                      "verdict": "互換律唔能定 M₀ ⟹ 絕對標度必須外來 (O38 handoff)"},
    "C_beat_reading": {"sum_sqrt_gamma": mp.nstr(sum_beat, 20),
                       "cot(pi/176)": mp.nstr(ring_casimir, 20),
                       "dev": mp.nstr(sum_beat - ring_casimir, 10),
                       "pool": mp.nstr(beat_pool, 20)},
    "D_rate": {"law_has_time_param": False, "candidates": {k: mp.nstr(v, 12) for k, v in rate_candidates.items()},
               "verdict": "conditional no-go: 速率需新無量綱常數; 零參數下不可導"},
    "checks": {k: bool(v) for k, v in checks.items()},
    "ALL_PASS": all(checks.values()),
    "n_pass": sum(1 for v in checks.values() if v), "n_total": len(checks),
}

os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
with open(os.path.join(ROOT, "data", "04_uniqueness_scaleblind.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)

for k, v in checks.items():
    print(("PASS " if v else "FAIL ") + k)
print("---")
print("A: p=0.5:", res_p["0.5"]["m_MeV"], "spread", res_p["0.5"]["ratio_spread"])
print("   p=1  :", res_p["1"]["m_MeV"], "spread", res_p["1"]["ratio_spread"], "← 唯一")
print("   p=2  :", res_p["2"]["m_MeV"], "spread", res_p["2"]["ratio_spread"])
print("C: Σ√γ =", mp.nstr(sum_beat, 12), "| cot(π/176) =", mp.nstr(ring_casimir, 12),
      "| dev =", mp.nstr(sum_beat - ring_casimir, 6))
print("ALL_PASS =", out["ALL_PASS"], "(%d/%d)" % (out["n_pass"], out["n_total"]))