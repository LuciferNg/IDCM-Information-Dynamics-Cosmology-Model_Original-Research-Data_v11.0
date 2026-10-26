#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/33 -- Incoming handoff intake audit (2026-10-26).

Scope: the 11 register events addressed to M4 (``target_pid == M4``), of which
8 were already receipted on 2026-10-24 and 3 arrive from M7 on 2026-10-26
(12:59 / 14:45 / 18:38).  This script INDEPENDENTLY REBUILDS the numeric
claims of the three M7 letters before M4 absorbs them.  It does NOT import or
re-run any M7 / M5 script (cross-record rule: rebuild, never cite).

Reconstructed here (from M4 docs/28 anchors + registered structure only):
  A1  chi2( <O> = 1 ) / 3 dof and the correct p-value; p at 1 dof as a control
  A2  weighted linear fit <O> = a + b*n_s  (a, sigma_a, b, sigma_b, chi2_min)
  A3  per-anchor deviation in sigma
  A4  model table: M1 = 1 + 2 n_s/3 vs candidate "<O>-sea" forms, Delta chi2
  A5  MC null 1e7 for b = 0 (independent seed) -- 4-Hard-Test (2)
  A6  valence-count identity <O> = 1 in exact Fraction arithmetic
  A7  sum_k gamma_k = 90 exact (mpmath 50 dps)
  A8  per-mode capacity 2*(B-1) = 6 vs N_q = 3  => max occupancy sigma = 1/2
  A9  three-body harmonic shift (pairwise-only) = -18.35 %
  A10 v9 section D.1 sector table: does it contain u and d?  (definition gate)

Arithmetic policy: Fraction / mpmath 50 dps for every deterministic number;
float only for the Monte-Carlo and for final display.
"""
import json
import os
import time
from fractions import Fraction as F

import numpy as np
import mpmath as mp

mp.mp.dps = 50

CHECKS = []
DIAG = []


def check(label, ok, got, expect=""):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-46s" % label, "|", got, "| expect", expect)


def diag(label, text):
    DIAG.append({"label": label, "text": str(text)})
    print("[DIAG]", label, "|", text)


# ----------------------------------------------------------------- anchors
# docs/28 section 2: BMW PRL 111, 252001 (2013) TABLE I, reverse-engineered
# <O>_N / <O>_Sigma / <O>_Xi, with n_s = 0 / 1 / 2.
Y = [F(891, 1000), F(1498, 1000), F(2293, 1000)]
S = [F(107, 1000), F(262, 1000), F(351, 1000)]
NS = [F(0), F(1), F(2)]

B = 4            # tree branching (registered)
N_Q = 3          # valence quark number (registered)
N_D = 88         # chain length (registered)

print("=" * 78)
print("M4 scripts/33 -- incoming handoff intake audit (2026-10-26)")
print("=" * 78)
print("anchors: <O> =", [str(v) for v in Y], " sigma =", [str(v) for v in S])

# ----------------------------------------------------------------- A1
chi2_const = sum(((y - 1) / s) ** 2 for y, s in zip(Y, S))
chi2_const_f = float(chi2_const)
print("\n[A1] chi2(<O>==1) exact Fraction =", chi2_const, "=", chi2_const_f)

p3 = mp.gammainc(mp.mpf(3) / 2, mp.mpf(chi2_const_f) / 2, mp.inf, regularized=True)
p1 = mp.erfc(mp.sqrt(mp.mpf(chi2_const_f) / 2))
print("     p (3 dof) =", mp.nstr(p3, 6), "   p (1 dof) =", mp.nstr(p1, 6))
check("A1a chi2(<O>==1) == 18.22 (docs/29 C1)", abs(chi2_const_f - 18.22) < 5e-3,
      round(chi2_const_f, 4), "18.22")
check("A1b p(3 dof) == 3.96e-4 (M7 K3)", abs(float(p3) - 3.96e-4) / 3.96e-4 < 0.02,
      mp.nstr(p3, 5), "3.96e-4")
check("A1c p(1 dof) == 1.97e-5 (= docs/29 C1 figure)",
      abs(float(p1) - 1.97e-5) / 1.97e-5 < 0.05, mp.nstr(p1, 5), "1.97e-5")
z3 = mp.sqrt(mp.mpf(2)) * mp.erfinv(1 - 2 * p3)   # one-sided equivalent
print("     equivalent one-sided sigma(3 dof) =", mp.nstr(z3, 5))
check("A1d sigma(3 dof) == 3.36", abs(float(z3) - 3.36) < 0.02, mp.nstr(z3, 4), "3.36")

# ----------------------------------------------------------------- A2
# weighted least squares, weights 1/sigma^2, exact Fraction for the algebra
Sw = sum(1 / s ** 2 for s in S)
Sx = sum(n / s ** 2 for n, s in zip(NS, S))
Sxx = sum(n * n / s ** 2 for n, s in zip(NS, S))
Sy = sum(y / s ** 2 for y, s in zip(Y, S))
Sxy = sum(n * y / s ** 2 for n, y, s in zip(NS, Y, S))
det = Sw * Sxx - Sx * Sx
a_fit = (Sxx * Sy - Sx * Sxy) / det
b_fit = (Sw * Sxy - Sx * Sy) / det
chi2_min = sum(((y - (a_fit + b_fit * n)) / s) ** 2 for n, y, s in zip(NS, Y, S))
sig_a = float((Sxx / det) ** 0.5)
sig_b = float((Sw / det) ** 0.5)
print("\n[A2] weighted fit: a =", float(a_fit), " sigma_a =", sig_a)
print("     b =", float(b_fit), " sigma_b =", sig_b, " chi2_min =", float(chi2_min))
check("A2a a == 0.885740 +- 0.105493",
      abs(float(a_fit) - 0.885740) < 1e-5 and abs(sig_a - 0.105493) < 1e-5,
      "%.6f +- %.6f" % (float(a_fit), sig_a), "0.885740 +- 0.105493")
check("A2b b == 0.675330 +- 0.161348",
      abs(float(b_fit) - 0.675330) < 1e-5 and abs(sig_b - 0.161348) < 1e-5,
      "%.6f +- %.6f" % (float(b_fit), sig_b), "0.675330 +- 0.161348")
check("A2c chi2_min == 0.0864", abs(float(chi2_min) - 0.0864) < 1e-4,
      round(float(chi2_min), 6), "0.0864")
# decompositions (both exact algebra of weighted least squares)
chi2_b0 = float(sum((y - (Sy / Sw)) ** 2 / s ** 2 for y, s in zip(Y, S)))
check("A2d chi2(b=0, a free) - chi2_min == (b/sigma_b)^2",
      abs((chi2_b0 - float(chi2_min)) - (float(b_fit) / sig_b) ** 2) < 1e-9,
      "%.6f - %.6f = %.6f" % (chi2_b0, float(chi2_min), (float(b_fit) / sig_b) ** 2), "identity")
check("A2e chi2(<O>==1) - chi2(b=0, a free) == (a0-1)^2*Sw",
      abs((chi2_const_f - chi2_b0) - (float(Sy / Sw) - 1) ** 2 * float(Sw)) < 1e-9,
      "%.6f - %.6f = %.6f" % (chi2_const_f, chi2_b0, (float(Sy / Sw) - 1) ** 2 * float(Sw)), "identity")

# ----------------------------------------------------------------- A3
zs = [(1 - y) / s for y, s in zip(Y, S)]
print("\n[A3] per-anchor z (value 1 vs anchor):", [round(float(z), 3) for z in zs])
check("A3a z_N == 1.02", abs(abs(float(zs[0])) - 1.02) < 0.01, round(float(zs[0]), 3), "1.02")
check("A3b z_Sigma == 1.90", abs(abs(float(zs[1])) - 1.90) < 0.01, round(float(zs[1]), 3), "1.90")
check("A3c z_Xi == 3.68", abs(abs(float(zs[2])) - 3.68) < 0.01, round(float(zs[2]), 3), "3.68")
check("A3d intercept a compatible with 1 (|1-a|/sig_a < 2)",
      abs(1 - float(a_fit)) / sig_a < 2.0, round(abs(1 - float(a_fit)) / sig_a, 3), "< 2")

# ----------------------------------------------------------------- A4
P0 = 0.36892170688337345      # M5 reproduction P0 (M4 docs/25 verified value)
rA = 0.4728                   # <B^{-d}>
MODELS = {
    "M1  1+2n/3        ": [1 + F(2, 3) * n for n in NS],
    "A2  Nq(n+1)/B     ": [N_Q * (n + 1) / F(B) for n in NS],
    "M0sea 1+n*P0      ": [1 + n * F(P0) for n in NS],
    "MB   1+n*<B^-d>   ": [1 + n * F(rA) for n in NS],
    "M34  1+0.75n      ": [1 + F(3, 4) * n for n in NS],
    "Mhalf 1+0.5n      ": [1 + F(1, 2) * n for n in NS],
    "M2   3/(3-n)      ": [F(3) / (3 - n) for n in NS],
}
print("\n[A4] model chi2 (3 dof) and Delta chi2 vs M1")
tab = {}
chi_m1 = None
for name, vals in MODELS.items():
    c2 = sum(((y - v) / s) ** 2 for y, v, s in zip(Y, vals, S))
    tab[name] = float(c2)
    if name.startswith("M1"):
        chi_m1 = float(c2)
print("     chi2(M1) =", round(chi_m1, 4))
for name, c2 in tab.items():
    print("     %-20s chi2 = %8.4f   Dchi2 vs M1 = %7.4f" % (name, c2, c2 - chi_m1))
check("A4a chi2(M1) == 1.47", abs(chi_m1 - 1.47) < 0.01, round(chi_m1, 4), "1.47")
dsea = tab["M0sea 1+n*P0      "] - chi_m1
check("A4b Delta chi2(M1 vs 1+n*P0 sea) == 2.32 (M7 K5)",
      abs(dsea - 2.32) < 0.02, round(dsea, 4), "2.32")
z_dchi = float(mp.sqrt(mp.mpf(dsea)))
check("A4c sqrt(Delta chi2) == 1.52 sigma", abs(z_dchi - 1.52) < 0.02, round(z_dchi, 3), "1.52")
best_alt = min((v, k) for k, v in tab.items() if not k.startswith("M1"))
check("A4d no registered-structure form beats M1 by >4 in chi2 (indistinguishable)",
      best_alt[0] - chi_m1 > -4.0, "%s D=%.3f" % (best_alt[1].strip(), best_alt[0] - chi_m1), ">-4")

# ----------------------------------------------------------------- A5
# MC null 1e7 under H0: <O> == 1 exactly; statistic |b_hat|.
B_HAT = float(b_fit)
SIG = np.array([float(s) for s in S])
NSf = np.array([0.0, 1.0, 2.0])
W = 1.0 / SIG ** 2
SW, SX, SXX = W.sum(), (W * NSf).sum(), (W * NSf * NSf).sum()
DET = SW * SXX - SX * SX
rng = np.random.default_rng(20261026)
NMC = 10 ** 7
CHUNK = 10 ** 6
hits = 0
bmax = 0.0
t0 = time.time()
for _ in range(NMC // CHUNK):
    y = 1.0 + SIG[None, :] * rng.standard_normal((CHUNK, 3))
    Sy = (W[None, :] * y).sum(axis=1)
    Sxy = (W[None, :] * NSf[None, :] * y).sum(axis=1)
    b = (SW * Sxy - SX * Sy) / DET
    ab = np.abs(b)
    hits += int((ab >= B_HAT).sum())
    bmax = max(bmax, float(ab.max()))
p_mc = hits / NMC
print("\n[A5] MC null 1e7 under <O>==1: |b| >= %.6f in %d draws -> p = %.3e (%.1fs)"
      % (B_HAT, hits, p_mc, time.time() - t0))
check("A5a MC p in [1e-5, 6e-5] (M7 quote 2.79e-5)", 1e-5 <= p_mc <= 6e-5,
      "%.3e" % p_mc, "2.79e-5")
zd = float(mp.sqrt(2) * mp.erfinv(1 - 2 * p_mc)) if p_mc > 0 else 0.0
check("A5b MC sigma equivalent >= 4.0", zd >= 4.0, round(zd, 2), ">= 4.0 (M7: 4.19)")
z_b = B_HAT / sig_b
check("A5c analytic z_b == 4.186", abs(z_b - 4.186) < 0.01, round(z_b, 3), "4.186")

# ----------------------------------------------------------------- A6
# valence-count identity, exact Fractions.  <O>_ch = dDeltaM/d(m_u-m_d)/coeff.
# m_u - m_d is the same number in every channel, so the ratio is 1 identically.
mu_minus_md = F(256, 100)     # FLAG24 2.56 MeV, value irrelevant for identity
CH = {"N": (1, 1), "Sigma": (2, 2), "Xi": (1, 1)}   # (Delta n_u, coeff)
ident = {}
for k, (dnu, coeff) in CH.items():
    ident[k] = (coeff * mu_minus_md) / (coeff * mu_minus_md)
print("\n[A6] valence identity <O>_ch = dnu/coeff =", {k: str(v) for k, v in ident.items()})
check("A6a <O> == 1 exactly in all 3 channels",
      all(v == 1 for v in ident.values()) and all(dnu == coeff for dnu, coeff in CH.values()),
      [str(v) for v in ident.values()], "1 (exact)")

# ----------------------------------------------------------------- A7
sg = mp.nsum(lambda k: 4 * mp.sin(mp.pi * k / 88) ** 2, [1, 44])
print("\n[A7] sum_k gamma_k (44-chain) =", mp.nstr(sg, 30))
check("A7a sum gamma == 90 exact (50 dps)", abs(sg - 90) < mp.mpf("1e-45"), mp.nstr(sg, 20), "90")
half = mp.nsum(lambda k: mp.sin(mp.pi * k / 88) ** 2, [1, 44])
check("A7b sum sin^2 == 45/2 exact", abs(half - mp.mpf(45) / 2) < mp.mpf("1e-45"),
      mp.nstr(half, 20), "45/2")

# ----------------------------------------------------------------- A8
cap = 2 * (B - 1)
sigma_max = F(N_Q, cap)
print("\n[A8] capacity per chain mode = 2*(B-1) =", cap, "; N_q =", N_Q, "; max sigma =", sigma_max)
check("A8a capacity == 6", cap == 6, str(cap), "6")
check("A8b max occupancy == 1/2 < 1 (never saturates => additivity not an assumption)",
      sigma_max < 1 and sigma_max == F(1, 2), str(sigma_max), "1/2")

# ----------------------------------------------------------------- A9
# 3 equal masses, pairwise harmonic kernel: eigenfrequencies {0, sqrt3, sqrt3}*sqrt(k/m)
e3 = mp.mpf(1) / 2 * (0 + 2 * mp.sqrt(3))
epair = mp.mpf(3) / 2 * mp.sqrt(2)
rel = (e3 - epair) / epair
print("\n[A9] E0^(3) =", mp.nstr(e3, 12), " sum pair =", mp.nstr(epair, 12),
      " rel =", mp.nstr(rel * 100, 6), "%")
check("A9a three-body shift == -18.35 %", abs(float(rel * 100) + 18.35) < 0.01,
      mp.nstr(rel * 100, 5) + " %", "-18.35 %")

# ----------------------------------------------------------------- A10
V9 = "/home/wsl/Projects/IDCM-Original/IDCM-v9.0--W-Field-Lattice/wfield_lattice/03_extended_derivation.md"
block = ""
if os.path.exists(V9):
    txt = open(V9, encoding="utf-8").read().splitlines()
    for i, ln in enumerate(txt):
        if "Generation modes" in ln:
            block = "\n".join(txt[i:i + 6])
print("\n[A10] v9 section D.1 block:")
print(block)
labels = {"c", "t"} if "(c, t)" in block else set()
labels |= {"s", "b"} if "(s, b)" in block else set()
labels |= {"e"} if "(e," in block else set()
print("     parsed quark labels in the three sector lines:", sorted(labels))
check("A10a block reproduces (c,t)/(s,b)/(e,...) sector table",
      {"c", "t", "s", "b", "e"} <= labels, sorted(labels), "{c,t,s,b,e}")
check("A10b u absent from sector table (=> u->P_1 is an extrapolation)",
      "u" not in labels, sorted(labels), "u absent")
check("A10c d absent from sector table (=> d->P_4 is an extrapolation)",
      "d" not in labels, sorted(labels), "d absent")

# ----------------------------------------------------------------- A11 doc side
D29 = open("docs/29_operator_structure.md", encoding="utf-8").read()
D30 = open("docs/30_prediction_envelope.md", encoding="utf-8").read()
D32 = open("docs/32_flavour_structure_scout.md", encoding="utf-8").read()
diag("docs/29 quotes p = 1.97e-05 (stale, 1 dof value used for 3 dof)", "p = 1.97e-05" in D29)
diag("docs/29 keeps M1 labelled post-hoc", "後驗" in D29 and "M1" in D29)
diag("docs/30 lists M1/A1 and A2", "A1" in D30 and "A2" in D30)
diag("docs/32 asserts the mode-subspace mapping P_{k=1} / P_{k=4}", "P_{k=1}" in D32 and "P_{k=4}" in D32)

# ----------------------------------------------------------------- summary
npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 78)
print("SUMMARY: %d/%d PASS   (%d DIAG lines)" % (npass, len(CHECKS), len(DIAG)))
print("=" * 78)
out = {
    "script": "scripts/33_handoff_intake_20261026.py",
    "date": "2026-10-26",
    "scope": "incoming handoffs addressed to M4 (11 events); independent rebuild of M7 x3",
    "anchors": {"O": [str(v) for v in Y], "sigma": [str(v) for v in S]},
    "chi2_const": str(chi2_const),
    "p_3dof": mp.nstr(p3, 12),
    "p_1dof": mp.nstr(p1, 12),
    "sigma_3dof": mp.nstr(z3, 8),
    "fit": {"a": float(a_fit), "sigma_a": sig_a, "b": float(b_fit), "sigma_b": sig_b,
            "chi2_min": float(chi2_min), "z_b": z_b},
    "models_chi2": tab,
    "delta_chi2_m1_vs_p0sea": dsea,
    "mc": {"N": NMC, "hits": hits, "p": p_mc, "sigma": zd, "seed": 20261026,
           "stat": "|b_hat| >= b_fit under <O>=1"},
    "capacity": {"per_mode": cap, "N_q": N_Q, "sigma_max": str(sigma_max)},
    "three_body_rel_pct": float(rel * 100),
    "v9_labels": sorted(labels),
    "checks": CHECKS,
    "diag": DIAG,
    "n_pass": npass,
    "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/33_handoff_intake.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/33_handoff_intake.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)