#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/08_framework_self_audit.py — 框架自審：M6 嘅定理有冇鑑別力？
================================================================================
調律者問：「M6 框架是否錯誤」。
本腳本唔推新結果，只做**鑑別力測試**（discriminating-power test）：

  A. T1（加性）到底測咗咩？—— 用隨機數據重跑同一套檢查
  B. T2（dof_scale = 1）到底係定理還是定義？—— 用任意玩具模型重跑
  C. M6 框架唯一**非定義性**嘅預言係咩？要咩數據才可判死？

判準（事前）：
  - 若「真 γ」與「隨機數據」通過率相同 ⟹ 該檢查**零鑑別力**（屬定義性檢查）
  - 若玩具模型亦得 dof = 1 ⟹ T2 係「單一有量綱輸入」之重述，唔係 IDCM 之發現

執行：/usr/bin/python3 scripts/08_framework_self_audit.py
"""
import json
import os
import random
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 40
OUT = {"checks": {}}
NP = NF = 0


def ck(ok, msg, detail=""):
    global NP, NF
    NP, NF = (NP + 1, NF) if ok else (NP, NF + 1)
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}" + (f"  {detail}" if detail else ""))
    return bool(ok)


ND = 44
GAMMA = np.array([float(4 * mp.sin(mp.pi * n / 88) ** 2) for n in range(1, ND + 1)])

print("=" * 92)
print("M6 框架自審 — 鑑別力測試（勿視為新結果）")
print("=" * 92)

# ══════════════════════════════════════════════════════════════
print("\nA — T1（加性）：同一套檢查，換數據")
rng = np.random.default_rng(20261026)
M0 = 27.900915


def additivity_violations(weights, n_trials=9100, seed=7):
    """複製 scripts/01 A6 之檢查：隨機分割 ⟹ Σ M(S_k) 是否 = M(full)。"""
    r = random.Random(seed)
    viol = 0
    n = len(weights)
    for _ in range(n_trials):
        k = r.randint(2, 8)
        cuts = sorted(r.sample(range(1, n), k - 1))
        parts, prev = [], 0
        for c in cuts + [n]:
            parts.append(weights[prev:c])
            prev = c
        total = sum(float(w.sum()) for w in parts)
        full = float(weights.sum())
        if abs(total - full) > 1e-9:          # 相對容差極嚴（這些量 O(90)）
            viol += 1
    return viol


v_real = additivity_violations(GAMMA)
print(f"      真 γ 譜（Σ = {GAMMA.sum():.3f}）：違反 {v_real}/9100")
rand_pos = [rng.uniform(0.001, 4.0, ND) for _ in range(3)]
v_rand = [additivity_violations(w) for w in rand_pos]
print(f"      隨機正數（3 seeds）：違反 {v_rand}/9100")
rand_signed = rng.uniform(-4.0, 4.0, ND)
v_signed = additivity_violations(rand_signed)
print(f"      隨機正負數：違反 {v_signed}/9100")
ck(v_real == 0 and all(v == 0 for v in v_rand) and v_signed == 0,
   "真譜／隨機正數／隨機正負數：**違反數全部 = 0** ⟹ T1 檢查零鑑別力",
   "⟹ T1 只驗證程式（分配律），唔驗證物理：任何 44 個數都「通過」")
OUT["checks"]["A"] = {"real": v_real, "random_positive": v_rand,
                      "random_signed": v_signed,
                      "verdict": "zero discriminating power"}

# ══════════════════════════════════════════════════════════════
print("\nB — T2（dof_scale = 1）：用任意玩具模型重跑")
# M6 docs/02 之論證：q = M₀^a × (無量綱數) ⟹ Jacobian 單列 ⟹ rank 1
def jacobian_rank(masses, m0):
    """(∂ln m_i/∂ln M₀) 之「矩陣」只有一列；rank 恆 = 1（若任何 m_i ≠ 0）。"""
    col = np.array([1.0 for _ in masses])     # m_i = M₀·u_i ⟹ 指數恆為 1
    return int(np.linalg.matrix_rank(col.reshape(-1, 1)))


ranks = []
for _ in range(100):
    u = rng.uniform(0.01, 100.0, 12)          # 12 個任意「質量」
    ranks.append(jacobian_rank(u * M0, M0))
print(f"      100 個玩具模型（m_i = M₀·u_i，u_i 任意）：rank 值集合 = {sorted(set(ranks))}")
ck(set(ranks) == {1},
   "全部 100 個玩具模型都得 dof_scale = 1 ⟹ T2 對**任何**「單一有量綱輸入」模型成立",
   "⟹ T2 係『IDCM 只有一個有量綱輸入』之重述，唔係對 IDCM 結構嘅發現")
OUT["checks"]["B"] = {"ranks": sorted(set(ranks)), "n_trials": 100,
                      "verdict": "definitional (single dimensionful input)"}

# ══════════════════════════════════════════════════════════════
print("\nC — M6 框架唯一非定義性預言：ℚ-span 關係")
# rank_ℚ{γ_n} = 20 ⟹ 所有 λ 落在 20 維 ℚ-向量空間 V
# ⟹ 任何 21 個質量（單量子或多量子）都滿足整數關係 Σ c_i m_i = 0
# 判死門檻：係數上界 H 之偽關係能以精度 δ 擬合 ⟹ δ ≳ H^{-(dim)}
dimV = 20
print(f"      dim_ℚ V = {dimV}（M6 docs/01）⟹ 21 個質量必有整數關係 Σc_i m_i = 0")
print(f"      M₀ 因子自動消去（λ_i = m_i/M₀，共同標量）⟹ 只須量質量，不須知 M₀")
print(f"\n      {'相對精度 δ':>14}{'可偽擬合之最大係數 H':>24}{'判死力':>12}")
rows = []
for d in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
    H = d ** (-1.0 / (dimV + 1))
    rows.append({"delta": d, "H_max": H})
    tag = "無（係數要 ≤ 1）" if H < 2 else ("弱" if H < 10 else "中")
    print(f"      {d:>14.0e}{H:>24.2f}{tag:>12}")
ck(rows[-1]["H_max"] < 10,
   f"要判死（H ≲ 3 ⟹ 偽關係幾乎不可能），需 21 個質量同時達到 δ ≲ {1e-10:.0e}",
   "⟹ **現時數據未夠格**：高精度質量（≤1e-8）只有 m_e / m_p / m_n / m_μ 等少數，湊唔夠 21 個")
OUT["checks"]["C"] = {"dimV": dimV, "thresholds": rows,
                      "verdict": "test exists but not yet decisive (need 21 masses @ 1e-10)"}

# 現時可用高精度質量數目（PDG 級，相對精度 ≤ 1e-8）
hi_prec = {"m_e": 1e-10, "m_p": 1e-10, "m_n": 1e-9, "m_μ": 1e-8}
print(f"\n      現時 ≤1e-8 精度之質量：{len(hi_prec)} 個（{', '.join(hi_prec)}）⟹ 距 21 個甚遠")
ck(len(hi_prec) < 21, f"高精度質量只有 {len(hi_prec)} 個 < 21 ⟹ ℚ-關係測試**暫不可判**")

print("\n" + "=" * 92)
print(f"自審：{NP}/{NP + NF} PASS")
print("=" * 92)
OUT["pass"], OUT["fail"] = NP, NF
os.makedirs("data", exist_ok=True)
with open("data/08_framework_self_audit.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/08_framework_self_audit.json")
sys.exit(0 if NF == 0 else 1)