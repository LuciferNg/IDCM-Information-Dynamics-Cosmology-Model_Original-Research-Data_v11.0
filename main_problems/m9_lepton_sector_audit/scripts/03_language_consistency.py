#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/03 — R3 扇區語言一致性（Q3 / T3）

(A) φ-ladder 加法一致性是否真證據？→ 格點整數算術 ⟹ 自動成立（構造 artifact）
(B) 兩個語言嘅可表示性 ＋ 解像度對照（φ-ladder 格 vs γ-ratio 譜集 vs 觀測精度）
(C) 判準結論：插值族（覆蓋任何目標）≠ 預言族（離散、可證偽）
零 import 上游。
"""
import json
import os

from mpmath import mp, mpf, log, sqrt

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

PHI = (1 + sqrt(5)) / 2
LN_PHI = log(PHI)

R_MUE, U_MUE = mpf("206.7682827"), mpf("0.0000046")
R_TAUMU, U_TAUMU = mpf("16.8170"), mpf("0.0011")
R_TAUE, U_TAUE = mpf("3477.23"), mpf("0.23")

checks = []


def chk(name, cond, detail=""):
    checks.append({"check": name, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


def gamma(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print("=" * 74)
print("M9 scripts/03 — R3 語言一致性")
print("=" * 74)

X = {"mu/e": (R_MUE, U_MUE), "tau/mu": (R_TAUMU, U_TAUMU), "tau/e": (R_TAUE, U_TAUE)}
XL = {k: log(v[0]) / LN_PHI for k, v in X.items()}

# ════════════════════════════════════════════════════════════════════
# (A) 加法一致性：格點算術自動成立
# ════════════════════════════════════════════════════════════════════
print("\n[A] φ-ladder 加法一致性檢查")
D = 528
b_best = {}
for k, xo in XL.items():
    frac = xo - mp.floor(xo)
    b_star = D * frac
    b_best[k] = int(round(b_star))
    print(f"  {k:7s} x_obs = {mp.nstr(xo, 12)}   b* = {mp.nstr(b_star, 8)}  → 取整 {b_best[k]:4d}"
          f"  (殘差 {mp.nstr(abs(b_star - b_best[k]), 4)} 格)")
sum_check = b_best["mu/e"] + b_best["tau/mu"]
print(f"  b(μ/e) + b(τ/μ) = {b_best['mu/e']} + {b_best['tau/mu']} = {sum_check}"
      f"   vs b(τ/e) = {b_best['tau/e']}")
chk("A1 格點加法恆等 b(μ/e)+b(τ/μ) = b(τ/e)（整數算術 ⟹ 自動）",
    sum_check == b_best["tau/e"],
    f"{sum_check} vs {b_best['tau/e']}")
# 只要三個比值都係同一格族嘅取整，加法一致性必然成立（或差 1 格）
# ⟹ O34 docs/07「加法一致性 ✅」係構造 artifact，非獨立證據
print("  ⟹ 加法一致性係取整算術之必然結果，對任何三個比值都成立 ⟹ 零資訊（非證據）")

# ════════════════════════════════════════════════════════════════════
# (B) 可表示性 ＋ 解像度
# ════════════════════════════════════════════════════════════════════
print("\n[B] 兩語言可表示性／解像度")
G = {n: gamma(n) for n in range(1, 45)}
CAP = G[44] / G[1]
print(f"  γ-ratio 譜集上限 CAP = γ₄₄/γ₁ = {mp.nstr(CAP, 10)}")
for k, (v, u) in X.items():
    reach = v <= CAP
    print(f"  {k:7s} obs = {mp.nstr(v, 9):>12s}   ≤ CAP : {reach}"
          f"   (需要 {mp.nstr(v / CAP, 5)}× CAP)" if not reach else
          f"  {k:7s} obs = {mp.nstr(v, 9):>12s}   ≤ CAP : {reach}")
chk("B1 τ/e 超出單模 γ-ratio 上限（結構上不可表示）", R_TAUE > CAP,
    f"{mp.nstr(R_TAUE / CAP, 5)}× CAP")

# γ-ratio 最近鄰（L0；詳細掃描喺 scripts/04）
vals = []
for a in range(1, 45):
    for b in range(1, 45):
        if a != b:
            vals.append((G[a] / G[b], a, b))
best = {}
for k, (v, u) in X.items():
    cand = min(vals, key=lambda t: abs(t[0] / v - 1))
    dev = abs(cand[0] / v - 1)
    best[k] = {"pair": (cand[1], cand[2]), "value": mp.nstr(cand[0], 10),
               "dev_rel": mp.nstr(dev, 6), "dev_pct": mp.nstr(dev * 100, 6),
               "sigma": mp.nstr(dev * v / u, 6)}
    print(f"  {k:7s} L0 最近 γ 比 = γ{cand[1]}/γ{cand[2]} = {mp.nstr(cand[0], 10)}"
          f"  dev = {mp.nstr(dev * 100, 6)}%  ({best[k]['sigma']}σ)")

# 鄰居間距（log 空間）統計 — 離散族詳情
lv = sorted(float(log(t[0])) for t in vals)
sp = [lv[i + 1] - lv[i] for i in range(len(lv) - 1) if lv[i + 1] - lv[i] > 1e-12]
sp.sort()
med = sp[len(sp) // 2]
print(f"  γ-ratio 族：{len(vals)} 值，log 空間鄰居間距 中位 = {med:.3e}"
      f"（= 相對 {med * 100:.3f}%），最細 = {sp[0]:.2e}")

# 解像度對照表
print("\n  解像度對照（『呢個語言可以講到幾細』）")
rows = []
for name, res, note in [
    ("φ-ladder 88 格（半格）", float(PHI ** (mpf(1) / 176) - 1), "插值族；任何目標 ≤ 0.274%"),
    ("φ-ladder 132 格（半格）", float(PHI ** (mpf(1) / 264) - 1), "插值族"),
    ("φ-ladder 528 格（半格）", float(PHI ** (mpf(1) / 1056) - 1), "插值族；任何目標 ≤ 0.0456%"),
    ("γ-ratio 族（鄰居中位）", med, "離散族；可完全錯開（可證偽）"),
    ("觀測 μ/e", float(U_MUE / R_MUE), "權威相對 u"),
    ("觀測 τ/μ", float(U_TAUMU / R_TAUMU), "權威相對 u"),
]:
    rows.append({"language": name, "resolution_rel": res, "note": note})
    print(f"    {name:26s} {res:.3e}   {note}")
chk("B2 φ-ladder 528 半格粗於 μ/e 精度 >1e4×",
    (float(PHI ** (mpf(1) / 1056) - 1)) / float(U_MUE / R_MUE) > 1e4,
    f"{(float(PHI ** (mpf(1) / 1056) - 1)) / float(U_MUE / R_MUE):.3e}×")

# ════════════════════════════════════════════════════════════════════
# (C) 判準
# ════════════════════════════════════════════════════════════════════
crit = {
    "插值族 (φ-ladder)": {
        "coverage": "任何目標都在半格內 ⟹ 命中率 100%（τ ≥ 半格）",
        "falsifiability": "零 —— 冇任何觀測值可以否證（除精確度 < 半格時全部否證）",
        "info_per_test_bits": "0（構造上必然命中）",
    },
    "離散族 (γ-ratio L0)": {
        "coverage": f"{len(vals)} 值，log 空間中位間距 {med:.3e}",
        "falsifiability": "有 —— 目標可以完全落空（實測：三個比值全部落空 ≥1%）",
        "info_per_test_bits": f"{float(log(mpf('784.96') / mpf('1') ) / LN_PHI):.2f} 之 log₂ 折減（見 scripts/04）",
    },
    "verdict": (
        "T3 判：φ-ladder 對三個輕子比值『全部命中』係插值族之必然（P=1），不構成同語言證據；"
        "γ-ratio 離散族對三個比值全部落空（≥1%）⟹ 輕子扇區在 σ 級無結構比值。"
        "『同一扇區兩種語言』之真正判準唔係一致性，而係可證偽性：插值族無資格做身份。"
    ),
}
print("\n[C] 判準")
for k, v in crit.items():
    print(f"  {k}: {v}")

out = os.path.join(DATA, "03_language.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump({"additivity_b": b_best, "cap": mp.nstr(CAP, 10), "l0_best": best,
               "gamma_ratio_n": len(vals), "neighbour_median_log": med,
               "resolution_table": rows, "criteria": crit, "checks": checks,
               "pass": sum(1 for c in checks if c["pass"]), "total": len(checks)},
              f, ensure_ascii=False, indent=2)
print("\n" + "=" * 74)
print(f"PASS {sum(1 for c in checks if c['pass'])}/{len(checks)}   -> {out}")
print("=" * 74)