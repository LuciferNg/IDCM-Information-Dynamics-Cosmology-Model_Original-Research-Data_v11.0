#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/17 — R13：κ₃ 估計量之公約穩健性自檢（交 M17 前之最後一關）
========================================================================
風險：M17 若用不同公約（實/複基底、Condon–Shortley 相位、3j 呼叫次序、網格 vs 求積），
      κ₃ 之值可能唔同 ⟹ 判死會誤。

本檔用**公約無關**之路徑交叉核：
  (i)  由實基底（cos/sin 係數）重建 f(θ,φ)，與由複基底 a_m 重建之 f 比較
  (ii) 直接球面求積 ∫ f³ dΩ，與 3j 收縮 I₃ 比較
       標準恆等式：∫Y_{ℓm₁}Y_{ℓm₂}Y_{ℓm₃}dΩ = √((2ℓ+1)³/4π)·(ℓℓℓ|000)·(ℓℓℓ|m₁m₂m₃)
       ⟹ ∫f³ = √((2ℓ+1)³/4π)(ℓℓℓ|000) · I₃
  (iii) 求積解析度獨立性（NG = 120 vs 200）
  (iv) 奇 ℓ：∫f³ ≡ 0（與 (ℓℓℓ|000) = 0 一致）
"""
import json
import os
import sys

import numpy as np
from scipy.special import gammaln, lpmv
from sympy.physics.wigner import wigner_3j

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
J13 = json.load(open(os.path.join(DATA, "13_shape_invariants.json"), encoding="utf-8"))
D13 = [l for l in J01["D"] if 6 <= l <= 27]
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items()}


def real_to_complex(coef, l):
    a = np.zeros(2 * l + 1, dtype=complex)
    a[l] = coef[l]
    for m in range(1, l + 1):
        am = (coef[l + m] - 1j * coef[l - m]) / np.sqrt(2)
        a[l + m] = am
        a[l - m] = (-1) ** m * np.conj(am)
    return a


def norm_lm(l, m):
    am = abs(m)
    return np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))


def f_real(coef, l, ct, ph):
    """實基底：m>0 → cos、m<0 → sin、m=0 → P"""
    out = np.zeros_like(ct, dtype=float)
    for m in range(-l, l + 1):
        c = coef[l + m]
        if abs(c) < 1e-16:
            continue
        am = abs(m)
        n = norm_lm(l, m) * (np.sqrt(2) if m != 0 else 1.0)
        P = lpmv(am, l, ct)
        out += c * n * P * (1.0 if m == 0 else (np.cos(am * ph) if m > 0 else np.sin(am * ph)))
    return out


def f_complex(a, l, ct, ph):
    out = np.zeros_like(ct, dtype=complex)
    for m in range(-l, l + 1):
        out += a[l + m] * norm_lm(l, m) * lpmv(abs(m), l, ct) * np.exp(1j * m * ph)
    return out


print("(i) 實基底 vs 複基底之 f(θ,φ) 重建", flush=True)
worst = 0.0
for l in D13:
    ct = np.cos(np.linspace(0.02, np.pi - 0.02, 137))
    ph = np.linspace(0, 2 * np.pi, 151)[:-1]
    CT, PH = np.meshgrid(ct, ph, indexing="ij")
    fr = f_real(INV[l], l, CT, PH)
    fc = np.real(f_complex(real_to_complex(INV[l], l), l, CT, PH))
    worst = max(worst, float(np.max(np.abs(fr - fc))))
ck("R13-基底一致", worst < 1e-10, f"max|f_real − Re f_complex| = {worst:.2e}（13 個 ℓ）")

print("\n(ii) 求積 ∫f³dΩ vs 3j 收縮 I₃", flush=True)


def quad(l, NG=200):
    xq, wq = np.polynomial.legendre.leggauss(NG)
    phq = 2 * np.pi * np.arange(2 * NG) / (2 * NG)
    CT, PH = np.meshgrid(xq, phq, indexing="ij")
    Wq = np.outer(wq, np.full(2 * NG, 2 * np.pi / (2 * NG))).ravel()
    f = f_real(INV[l], l, CT, PH)
    return float(np.sum((f ** 3).ravel() * Wq))


def I3_exact(l):
    a = real_to_complex(INV[l], l)
    s = 0j
    for m1 in range(-l, l + 1):
        for m2 in range(-l, l + 1):
            m3 = -m1 - m2
            if -l <= m3 <= l:
                s += float(wigner_3j(l, l, l, m1, m2, m3)) * a[l + m1] * a[l + m2] * a[l + m3]
    return s


ratios = {}
for l in D13:
    q = quad(l)
    I3 = I3_exact(l)
    pred_const = np.sqrt((2 * l + 1) ** 3 / (4 * np.pi)) * float(wigner_3j(l, l, l, 0, 0, 0))
    if abs(I3) > 1e-14:
        r = q / float(np.real(I3)) if l % 2 == 0 else float("nan")
        ratios[l] = {"quad": q, "I3": float(np.real(I3)), "ratio": r, "pred": pred_const}
        print(f"   ℓ = {l:2d}：∫f³ = {q:+.8f}｜I₃ = {np.real(I3):+.8f}"
              f"｜比 = {r:.8f}｜預期常數 = {pred_const:.8f}"
              f"｜{'✅' if abs(r - pred_const) < 1e-6 else '❌'}", flush=True)
    else:
        ratios[l] = {"quad": q, "I3": 0.0, "ratio": None, "pred": pred_const}
        print(f"   ℓ = {l:2d}（奇）：∫f³ = {q:+.2e}｜I₃ = 0｜(ℓℓℓ|000) = {pred_const:.2e}", flush=True)

evenr = {l: v for l, v in ratios.items() if v["ratio"] is not None}
ck("R13-求積=收縮", all(abs(v["ratio"] - v["pred"]) < 1e-6 for v in evenr.values()),
   f"偶 ℓ（{len(evenr)} 個）：∫f³/I₃ 全等於 √((2ℓ+1)³/4π)(ℓℓℓ|000) ⟹ **公約無關路徑一致**")
ck("R13-奇ℓ恆零", all(abs(ratios[l]["quad"]) < 1e-8 for l in D13 if l % 2 == 1),
   "奇 ℓ：∫f³ ≈ 0（4 個；與 (ℓℓℓ|000)=0 及 I₃≡0 三重一致）")

print("\n(iii) 求積解析度獨立性（NG = 120 vs 200）", flush=True)
devs = {}
for l in [6, 10, 12, 26]:
    q1, q2 = quad(l, 120), quad(l, 200)
    ref = max(abs(q2), 1e-30)
    devs[l] = abs(q1 - q2) / ref
    print(f"   ℓ = {l:2d}：NG=120 → {q1:+.9f}｜NG=200 → {q2:+.9f}｜相對差 = {devs[l]:.2e}", flush=True)
ck("R13-解析度獨立", all(d < 1e-6 for d in devs.values()),
   f"NG 120 vs 200：max 相對差 = {max(devs.values()):.2e}")

print("\n(iv) 由求積反推 κ₃（M17 可直接用之路徑）", flush=True)
ok = {}
for l in D13:
    v = ratios[l]
    if v["ratio"] is None:
        continue
    k3_quad = (v["quad"] / v["pred"]) / 1.0        # I₁ = 1（I_ℓ 已歸一）
    k3_13 = J13["model"][str(l)]["kappa3_real"]
    ok[l] = abs(k3_quad - k3_13)
    print(f"   ℓ = {l:2d}：κ₃（求積路徑）= {k3_quad:+.6f}｜κ₃（3j 路徑）= {k3_13:+.6f}"
          f"｜差 = {abs(k3_quad-k3_13):.2e}", flush=True)
ck("R13-κ₃兩路一致", all(d < 1e-6 for d in ok.values()),
   f"偶 ℓ：兩路徑 κ₃ 最大差 = {max(ok.values()):.2e} ⟹ M17 可用「∫f³ 求積」實作，無公約風險")

npass = sum(1 for _, okk, _ in CHECKS if okk)
print(f"\n=== scripts/17：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"basis_dev": worst, "ratios": {str(k): {kk: vv for kk, vv in v.items()} for k, v in ratios.items()},
       "res_dev": {str(k): v for k, v in devs.items()}, "k3_two_paths": {str(k): v for k, v in ok.items()},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "17_r13_estimator_robustness.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
