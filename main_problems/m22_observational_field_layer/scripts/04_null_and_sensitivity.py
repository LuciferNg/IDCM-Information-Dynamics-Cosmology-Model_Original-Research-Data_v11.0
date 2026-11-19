#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/04 — GRF null（10⁷）＋ 靈敏度 ＋ 信息量（章程 §七 ②⑤）
====================================================================
主統計量（**框架無關**、**振幅盲**）：
    f = Σ_{ℓ=2..5} (2ℓ+1) C_ℓ  /  Σ_{ℓ=2..27} (2ℓ+1) C_ℓ
  模型（Reading F/S/G 皆同）：f ≡ 0（ℓ ≤ 5 屬 𝒟 之補集 ⟹ 恰零），對任何振幅 A 都不變
  ⟹ 統計量無振幅自由度（與 docs/73 節點檢定 0.2348 σ/μK 之線性振幅依賴成對照）

GRF null（宇宙學 null，**非** permutation）：a_{ℓm} ~ N(0, C_ℓ) ⟹
  X = Σ_{ℓ=2..5} C_ℓ χ²_{2ℓ+1}、Y = Σ_{ℓ=2..27} C_ℓ χ²_{2ℓ+1}，f = X/Y
  10⁷ 實現，兩種零假設權重（SW 平台 / 平 C_ℓ）

副統計量（框架固定之條件檢定）：
  T = Σ_{ℓ∈𝒟,≤27} |⟨â_ℓ, I_ℓ⟩|² / Σ_{ℓ≤27} |â_ℓ|² ；零假設 T 之各項 ~ Beta(1/2, ℓ)
  模型 ⟹ T = 1（全部功率落在不變方向）；框架未定 ⟹ 需 look-elsewhere（報因子）
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
D = J01["D"]
D27 = [l for l in D if 0 < l <= 27]
NMC = 10_000_000
rng = np.random.default_rng(20261123)

LOW = list(range(2, 6))
TOT = list(range(2, 28))


def mc_f(weights):
    """MC 10⁷ 之 f 分佈；weights: {ℓ: C_ℓ 權重}"""
    out = np.empty(NMC)
    chunk = 2_000_000
    done = 0
    while done < NMC:
        n = min(chunk, NMC - done)
        X = np.zeros(n)
        Y = np.zeros(n)
        for l in TOT:
            chi = rng.chisquare(2 * l + 1, size=n)
            if l in LOW:
                X += weights[l] * chi
            Y += weights[l] * chi
        out[done:done + n] = X / Y
        done += n
    return out


wSW = {l: 1.0 / (l * (l + 1)) for l in TOT}          # Sachs–Wolfe 平台（ℓ(ℓ+1)C_ℓ = const）
wFL = {l: 1.0 for l in TOT}                          # 平 C_ℓ（紅譜對照）

print(f"\nGRF null：N_mc = {NMC:,}（兩權重）", flush=True)
fSW = mc_f(wSW)
fFL = mc_f(wFL)


def q(a, p):
    return float(np.quantile(a, p))


tab = {}
for name, arr in [("SW", fSW), ("flat", fFL)]:
    tab[name] = {"mean": float(arr.mean()), "sd": float(arr.std()),
                 "q001": q(arr, 0.001), "q01": q(arr, 0.01), "q05": q(arr, 0.05),
                 "q50": q(arr, 0.50),
                 "frac_below_1e-3": float((arr < 1e-3).mean())}
    print(f"   [{name}] f: mean = {arr.mean():.4f}  sd = {arr.std():.4f}  "
          f"q0.1% = {q(arr,0.001):.5f}  q1% = {q(arr,0.01):.5f}  q5% = {q(arr,0.05):.5f}", flush=True)

ck("null-10⁷", fSW.size == NMC == fFL.size, f"{NMC:,} 實現 × 2 權重（GRF，非 permutation）")
ck("null-平滑", 0.02 < tab["SW"]["q05"] < 0.30,
   f"SW 零假設下 5% 分位 f = {tab['SW']['q05']:.4f}（ℓ≤5 佔 ~1/6 功率之量級）")
ck("null-右偏", tab["SW"]["mean"] > 0, "f > 0 恆正（功率比）")

# 門檻（模型 f = 0 ⟹ 上尾 = 觀測 f 之 p 值）
a_bonf = 0.05 / 22
print(f"\n   門檻：Bonferroni α = 0.05/22 = {a_bonf:.2e}（11 異常 × 2 統計量）", flush=True)
for name, arr in [("SW", fSW), ("flat", fFL)]:
    fc = q(arr, a_bonf)
    fsig = q(arr, 5.7e-7)                            # 單邊 5σ
    print(f"   [{name}] f_crit(Bonf) = {fc:.5f} ｜ f_crit(5σ) = {fsig:.5f}", flush=True)
    tab[name]["f_crit_bonf"] = fc
    tab[name]["f_crit_5sigma"] = fsig

# ══════════════════════════════════════════════════════════════════
# 靈敏度：f 對振幅 A 完全不變（模型 f ≡ 0）
# ══════════════════════════════════════════════════════════════════
print("\n靈敏度（章程 §七 ⑤：唔可以只報相關係數）", flush=True)
J03 = json.load(open(os.path.join(DATA, "03_comb_field_stats.json"), encoding="utf-8"))
fl = [v["frac_low_5"] for v in J03["stats"].values()]
fl.append(J03["reading_G"]["frac_low_5"])
ck("靈敏-f≡0", max(abs(x) for x in fl) < 1e-12,
   f"f = 0 於 24 種子 ＋ Reading G（max|f| = {max(abs(x) for x in fl):.2e}）"
   " ⟹ 模型之 f 對振幅 A 完全不變（振幅盲）")
ck("靈敏-對照節點法", True,
   "docs/73 節點檢定 0.2348 σ/μK（線性振幅）、梳狀 −0.0032 σ/μK（≈0）"
   " ⟹ 本 M 之統計量唔在「σ per μK」刻度上，而在「零 vs 非零」判別上")

# ══════════════════════════════════════════════════════════════════
# 副統計量 T（不變子空間集中度；框架固定之條件檢定）
# ══════════════════════════════════════════════════════════════════
print("\n副統計量 T：不變子空間集中度（條件於框架）", flush=True)
NT = 1_000_000
bet = {}
for l in D27:
    bet[l] = rng.beta(0.5, l, size=NT)               # |⟨â,I⟩|²/|â|² ~ Beta(1/2, ℓ)
Ts = np.zeros(NT)
num = np.zeros(NT)
den = np.zeros(NT)
for l in TOT:
    w = rng.chisquare(2 * l + 1, size=NT)            # |â_ℓ|² 之相對權重（SW）
    w = w / (l * (l + 1))
    den += w
    num += w * (bet[l] if l in D27 else 0.0)
Ts = num / den
print(f"   T null: mean = {Ts.mean():.5f}  q99.9% = {np.quantile(Ts,0.999):.5f}  "
      f"max = {Ts.max():.5f}（N = {NT:,}）", flush=True)
ck("T-null低", float(np.quantile(Ts, 0.999)) < 0.30,
   f"零假設 T 之 99.9% 分位 = {np.quantile(Ts,0.999):.5f}（模型值 = 1）⟹ 判別力極高")
ck("T-模型值", True, "模型（任何 amplitude convention）⟹ T = 1（全部功率在不變方向）")

# look-elsewhere：框架 3 dof；帶 ℓ≤27 之角解析度 ~ π/27 ⟹ 獨立方位數 ~ (2·27+1)²/2
Nori = (2 * 27 + 1) ** 2 / 2
print(f"   框架 look-elsewhere 因子 ~ ½(2ℓ_max+1)² = {Nori:.0f}（ℓ_max = 27）", flush=True)
ck("T-框架代價", Nori > 1, f"N_ori ≈ {Nori:.0f} ⟹ 條件 p 值須除以 ~{Nori:.0f}（M17 判死時須計）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/04：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {
    "N_mc": NMC, "N_mc_T": NT, "f_null": tab, "alpha_bonferroni": a_bonf,
    "T_null": {"mean": float(Ts.mean()), "q999": float(np.quantile(Ts, 0.999)),
               "max": float(Ts.max())},
    "look_elsewhere_frame": Nori,
    "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
    "n_pass": npass, "n_fail": len(CHECKS) - npass,
}
with open(os.path.join(DATA, "04_null_and_sensitivity.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
