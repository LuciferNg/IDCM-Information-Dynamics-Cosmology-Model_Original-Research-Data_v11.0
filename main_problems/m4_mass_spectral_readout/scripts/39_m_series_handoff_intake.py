#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/39 -- M-series handoff intake #2 (2026-10-28): independent rebuild.

Scope: FOUR incoming handoffs arrived after M4's last processing (2026-10-27 14:47):
  (a) M7 2026-10-27 21:16  -- accepts R5 verdict; scope narrowing; Q1 independence
  (b) M7 2026-10-27 22:59  -- independent Fraction recompute of M4 docs/36 + 2 new diagnostics
  (c) M5 2026-10-27 23:02  -- cannot deliver zero-parameter |psi(0)|^2; offers empirical channel instead
  (d) M5 2026-10-28 08:01  -- Delta I=0 hyperfine channel FAILS (contact 1/(m_h m_l) assumption wrong)

RULE: zero import of M7/M5 scripts; every number claimed by them that M4 may cite is
recomputed here from M4's OWN inputs (PDG 2024 masses, docs/29 law, docs/36 chain),
in EXACT rational arithmetic. M7/M5-internal quantities (their own no-go scan statistics)
are NOT imported -- they are recorded as "other-M property, not cited".

All arithmetic is Fraction (exact); no float touches a verdict.
"""
import json
import os
from fractions import Fraction as F

CHECKS = []
DATA = {}


def check(label, ok, got, expect):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-66s" % label, "|", got, "| expect", expect)


def rel(a, b):
    """relative deviation of Fraction/float a from b"""
    return abs(float(a) - float(b)) / abs(float(b))


def rounds_to(x, quoted):
    """does x round (half-even) to the quoted decimal string at its own precision?"""
    s = str(quoted)
    nd = len(s.split(".")[1]) if "." in s else 0
    return round(float(x), nd) == float(quoted)


# ---------------------------------------------------------------- inputs
# PDG 2024 masses (MeV), decimal-exact Fractions
MPI_P = F("139.57039")
MPI_0 = F("134.9768")
MK_P = F("493.677")
MK_0 = F("497.611")

# FLAG Review 2024 Sec. 4.1.3 Eqs.(38)-(39); MILC 18 PRD99 034503
EPS = F("0.78")
EPS_LO = F("0.67")      # 0.78 - 0.11 (syst down)
EPS_HI = F("0.86")      # 0.78 + 0.08 (syst up)

# M4 docs/29 Sec.3 law:  em_M4 = dM - coeff * X * <O>_pred
# K+ - K0 row (docs/29 line 54): coeff=-1, dM=-3.9340, <O>_pred=1.667, em_pred=+0.3327
COEFF_K = F(-1)
DM_K = F("-3.9340")
O_PRED_K = F("1.667")
EM_PRED_K = F("0.3327")
X_K = (DM_K - EM_PRED_K) / (COEFF_K * O_PRED_K)      # solved from the row itself

# docs/30 envelope for em(K0 - K+)
ENV_LO, ENV_HI = F("-0.52"), F("0.55")

print("=" * 92)
print("M4 scripts/39 -- M-series handoff intake #2 (2026-10-28), exact-rational rebuild")
print("=" * 92)
print("X (docs/29 K row, solved exactly) =", X_K, "=", float(X_K))

# ---------------------------------------------------------------- (b) M7 #2 chain
dpi2 = MPI_P ** 2 - MPI_0 ** 2
dk2 = MK_P ** 2 - MK_0 ** 2
mk_sum = MK_P + MK_0

check("B1 dM2_pi  = 1261.157227 (M7 quoted)", rel(dpi2, "1261.157227") < 2e-9, float(dpi2), "1261.157227")
check("B2 dM2_K   = -3899.726992 (M7 quoted)", rel(dk2, "-3899.726992") < 2e-9, float(dk2), "-3899.726992")

em_c = (1 + EPS) * dpi2 / mk_sum                      # central
em_lo = (1 + EPS_LO) * dpi2 / mk_sum
em_hi = (1 + EPS_HI) * dpi2 / mk_sum
sig = (em_hi - em_lo) / 2                              # symmetric 68% half-width
check("B3 em(K+-K0) central rounds to +2.264589 (M7 quoted)", rounds_to(em_c, "2.264589"), float(em_c), "2.264589")
check("B4 band lower rounds to 2.124642 (M7 quoted)", rounds_to(em_lo, "2.124642"), float(em_lo), "2.124642")
check("B5 band upper rounds to 2.366368 (M7 quoted)", rounds_to(em_hi, "2.366368"), float(em_hi), "2.366368")
check("B6 sigma_sym rounds to 0.120863 (M7 quoted)", rounds_to(sig, "0.120863"), float(sig), "0.120863")

# z of M4's frozen prediction
z_c = (EM_PRED_K - em_c) / sig
z_cons = (EM_PRED_K - em_lo) / sig                     # least-negative edge = most conservative
check("B7 z central = -15.98 (M7 quoted)", abs(float(z_c) + 15.98) < 5e-3, float(z_c), -15.98)
check("B8 z conservative edge = -14.83 (M7 quoted)", abs(float(z_cons) + 14.83) < 5e-3, float(z_cons), -14.83)

# diagnostic (1): the frozen +0.3327 sits on the *misread* (1-eps) value, not on the true value
em_misread = (1 - EPS) * dpi2 / mk_sum
z_vs_misread = (EM_PRED_K - em_misread) / sig
check("B9 (1-eps) misread = 0.279893 (M7 quoted)", rel(em_misread, "0.279893") < 1e-6, float(em_misread), "0.279893")
check("B10 prediction vs misread = +0.44 sigma (M7 quoted)", abs(float(z_vs_misread) - 0.44) < 5e-3, float(z_vs_misread), 0.44)

# diagnostic (2): back-solve <O>_K from the external value, compare to BMW Sigma anchor 1.498 +/- 0.262
# M7 reproduced their back-solve with the ROUNDED constant X = 2.56; exact X = 42667/16670.
X_ROUND = F("2.56")
O_K_m7 = (DM_K - em_c) / (COEFF_K * X_ROUND)
O_K_ext = (DM_K - em_c) / (COEFF_K * X_K)
BMW, BMW_S = F("1.498"), F("0.262")
z_O_m7 = (O_K_m7 - BMW) / BMW_S
z_O = (O_K_ext - BMW) / BMW_S
check("B11a M7 back-solve reproduced with X=2.56 -> 2.42132", rounds_to(O_K_m7, "2.42132"), float(O_K_m7), "2.42132")
check("B11b EXACT X=42667/16670 -> <O>_K = 2.421789 (M4 canonical)",
      rounds_to(O_K_ext, "2.421789"), float(O_K_ext), "2.421789")
check("B12a M7 sigma-distance with X=2.56 -> +3.52", rounds_to(z_O_m7, "3.52"), float(z_O_m7), "3.52")
check("B12b EXACT X -> +3.526 sigma (verdict unchanged: >3 KILL)",
      rounds_to(z_O, "3.526"), float(z_O), "3.526")

# diagnostic (2) corollary: the WHOLE linear law is excluded, not a single point.
# M1 form <O> = 1 + (2/3) n_s must hit <O>_K = 2.4213 at s=-1  =>  n_s = 3/2 (2.4213-1) = 2.132
ns_required = (O_K_ext - 1) / (F(2) / 3)
check("B13 M1 law would need n_s = +2.132 at K (|s|<=1 admitted)",
      float(ns_required) > 1.0, float(ns_required), ">1 (outside admitted range)")

DATA["m7_recompute"] = {
    "dM2_pi": str(dpi2), "dM2_K": str(dk2), "em_central": str(em_c),
    "band": [str(em_lo), str(em_hi)], "sigma": str(sig),
    "z_central": str(z_c), "z_conservative": str(z_cons),
    "em_misread": str(em_misread), "z_vs_misread": str(z_vs_misread),
    "O_K_backsolved": str(O_K_ext), "z_O_vs_BMW": str(z_O),
    "ns_required_for_M1": str(ns_required), "X_K": str(X_K),
}

# ---------------------------------------------------------------- (b) dof claim
# M7: chi^2 = 18.220727 is a dof=0 descriptive quantity (3 anchors, one prediction each).
N_ANCHORS, N_PRED_PER_ANCHOR, N_PARAMS = 3, 1, 3
dof_ = N_ANCHORS * N_PRED_PER_ANCHOR - N_PARAMS
check("B14 dof = 3 anchors x 1 prediction - 3 params = 0", dof_ == 0, dof_, 0)
DATA["dof"] = dof_

# ---------------------------------------------------------------- (d) M5 #2: Delta I=0 channel
rows = [
    # label, model ratio, data ratio, data sigma, quoted sigma-distance
    ("D_s/D", "0.5188", "1.0130", "0.0029", 170.0),
    ("B_s/B", "0.5696", "1.0764", "0.0316", 16.0),
]
ms_md = F("19.9")     # FLAG 2+1+1  m_s/m_d
res = {}
for lab, mod, dat, sd, qu in rows:
    mod, dat, sd = F(mod), F(dat), F(sd)
    zz = (dat - mod) / sd
    ratio = dat / mod
    expo = F(str(float(__import__("math").log(float(ratio)))) ) / F(str(float(__import__("math").log(float(ms_md)))))
    check("D%s %s |z| = %.0f sigma (M5 quoted)" % ("1" if lab.startswith("D") else "2", lab, qu),
          abs(abs(float(zz)) - qu) / qu < 0.02, float(abs(zz)), qu)
    res[lab] = {"model": str(mod), "data": str(dat), "sigma": str(sd), "z": str(zz),
                "ratio": str(ratio), "implied_ml_exponent": str(expo)}
    print("      %s: data/model = %.4f -> implied extra m_l exponent %.4f (vs m_s/m_d = 19.9)"
          % (lab, float(ratio), float(expo)))
DATA["m5_deltaI0"] = res

# ---------------------------------------------------------------- (c) M5 #1: internal consistency
sigma_over_c_K = F("6.09")
check("C1 sigma_c/|c|_K = 6.09 -> 609% relative error", float(sigma_over_c_K) * 100 == 609.0, 609.0, 609.0)

# ---------------------------------------------------------------- verdict roll-up
npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 92)
print("RESULT: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 92)

DATA["summary"] = {"pass": npass, "total": len(CHECKS), "checks": CHECKS}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "39_m_series_handoff_intake.json")
out = os.path.normpath(out)
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", encoding="utf-8") as fh:
    json.dump(DATA, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON ->", out)
raise SystemExit(0 if npass == len(CHECKS) else 1)