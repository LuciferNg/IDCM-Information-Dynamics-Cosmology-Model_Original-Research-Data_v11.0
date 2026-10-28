#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/40 -- recheck of three new incoming handoffs (2026-10-28 09:23 / 09:36 / 09:36).

  (1) M7 09:23  : (a) X provenance challenge  (b) challenge to M4's "|n_s|<=1" wording
                  (c) cleaner kill: n_s(K)=1 fixed by strangeness => M1 gives 5/3 exactly
  (2) M10 09:36 : W4 representation dependence of the sign-precision requirement (734x)
  (3) M10 09:36 : W6 pool-density 20.61% is channel incompleteness, not slot density

CRITICAL SELF-CHECK TRIGGERED BY M7:
  M4 docs/18 Sec.2 states  X = m_d - m_u = 2.560 MeV (FLAG24, EXTERNAL, sigma_X = 0.094).
  M4's docs/39 quoted 42667/16670 = 2.5595080984 as "exact canonical".  That value was
  BACK-SOLVED from the 4-decimal rounded row of docs/29 (dM=-3.9340, em=0.3327, <O>=1.667)
  => it is a ROUND-TRIP ARTIFACT, not a precision improvement.  X = 2.560 is canonical.
  FOLLOW-UP FOUND HERE (not raised by M7): sigma_X = 0.094 propagates to the PREDICTION
  with |coeff|*<O> = 5/3 => a 0.157 MeV error bar that docs/36's z IGNORED.
"""
import json
import os
from fractions import Fraction as F

CHECKS = []
DATA = {}


def check(label, ok, got, expect):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-70s" % label, "|", got, "| expect", expect)


# ---------------- inputs (sources pinned)
X_FLAG = F("2.560")        # docs/18 Sec.2 line 18: X = m_d - m_u = 2.560 MeV (FLAG24, external)
SIG_X = F("0.094")         # docs/18 line 35/49: X = 2.56 +- 0.094 (FLAG24)
DM_K = F("-3.9340")        # docs/29 Sec.3 K row
COEFF_K = F(-1)
O_ROW = F("1.667")         # docs/29 K row <O>_pred (rounded)
EM_ROW = F("0.3327")       # docs/29 K row em_pred (rounded, frozen 2026-10-25)
EM_EXT = F("2.26458896223049")   # docs/36 lattice chain (scripts/36b/39)
SIG_EM = F("0.12086289405162727")  # docs/36: sigma from eps band only
BMW, BMW_S = F("1.498"), F("0.262")

print("=" * 96)
print("M4 scripts/40 -- M7/M10 recheck + self-correction of X provenance and error budget")
print("=" * 96)

# ---------------- (1a) X canonical vs back-solve
X_BACK = (DM_K - EM_ROW) / (COEFF_K * O_ROW)          # = 42667/16670
check("E1 X back-solved from docs/29 rounded row = 42667/16670",
      X_BACK == F(42667, 16670), str(X_BACK), "42667/16670")
check("E2 |X_back - X_FLAG|/X_FLAG < 0.05% => back-solve is a rounding artifact, NOT canonical",
      abs(float(X_BACK - X_FLAG)) / float(X_FLAG) < 5e-4,
      float(abs(X_BACK - X_FLAG) / X_FLAG), "<5e-4")
check("E3 X_BACK < X_FLAG (artifact sits below the FLAG value)", X_BACK < X_FLAG,
      float(X_BACK), float(X_FLAG))

# ---------------- (1a') <O>_K under both X
O_K_flag = (DM_K - EM_EXT) / (COEFF_K * X_FLAG)
O_K_back = (DM_K - EM_EXT) / (COEFF_K * X_BACK)
z_flag = (O_K_flag - BMW) / BMW_S
z_back = (O_K_back - BMW) / BMW_S
check("E4 <O>_K (X=2.560, canonical) = 2.421324 (M7 quoted 2.421324)",
      round(float(O_K_flag), 6) == 2.421324, float(O_K_flag), 2.421324)
check("E5 z (X=2.560) = +3.524 (M7 quoted +3.524)",
      round(float(z_flag), 3) == 3.524, float(z_flag), 3.524)
check("E6 <O>_K (X back) = 2.421789 (M4 docs/39 value; now DEMOTED to consistency check)",
      round(float(O_K_back), 6) == 2.421789, float(O_K_back), 2.421789)
check("E7 relative difference between the two X <O>_K values = 0.019% (M7 quoted 0.019%)",
      round(float(abs(O_K_back - O_K_flag) / O_K_flag) * 100, 3) == 0.019,
      float(abs(O_K_back - O_K_flag) / O_K_flag) * 100, 0.019)

# ---------------- (1b) required n_s  vs the REAL domain {0,1,2}
ns_flag = (O_K_flag - 1) / (F(2) / 3)
ns_back = (O_K_back - 1) / (F(2) / 3)
check("E8 n_s required (X=2.560) = 2.131986 (M7 quoted 2.131986)",
      round(float(ns_flag), 6) == 2.131986, float(ns_flag), 2.131986)
check("E9 n_s required (X back) = 2.132684 (M7 quoted 2.132684)",
      round(float(ns_back), 6) == 2.132684, float(ns_back), 2.132684)
NS_DOMAIN = (0, 1, 2)     # docs/29 Sec.5 V2: n_s = 0,1,2 must fit SIMULTANEOUSLY; Xi anchor = 2
check("E10 both required n_s exceed max reachable n_s = 2 (domain {0,1,2})",
      float(ns_flag) > max(NS_DOMAIN) and float(ns_back) > max(NS_DOMAIN),
      "%.4f / %.4f vs max 2" % (float(ns_flag), float(ns_back)), "both > 2")

# ---------------- (1c) cleaner kill: K+ = u s-bar => n_s(K) = 1 exactly
N_S_K = F(1)              # strangeness of K+ (u s-bar) = +1
O_M1_K = 1 + (F(2) / 3) * N_S_K
check("E11 M1 at n_s(K)=1 gives <O>_K = 5/3 exactly", O_M1_K == F(5, 3), str(O_M1_K), "5/3")
check("E12 docs/29 K row <O>_pred = 1.667 == 5/3 rounded (no interpolation used)",
      round(float(O_M1_K), 3) == float(O_ROW), round(float(O_M1_K), 3), float(O_ROW))
em_M1_K = DM_K - COEFF_K * X_FLAG * O_M1_K
check("E13 M1@n_s=1 reproduces the frozen prediction: em = 0.332667 -> 0.3327",
      round(float(em_M1_K), 4) == float(EM_ROW), float(em_M1_K), float(EM_ROW))
z_direct = (em_M1_K - EM_EXT) / SIG_EM
check("E14 direct z at n_s=1 = -15.98 (M7: no n_s extrapolation or bound needed)",
      round(float(z_direct), 2) == -15.98, float(z_direct), -15.98)

# ---------------- (1d) ★ M4 FOUND: the error budget of docs/36 z ignored sigma_X
sig_pred = abs(COEFF_K * O_M1_K) * SIG_X                      # 5/3 * 0.094
sig_tot = (sig_pred ** 2 + SIG_EM ** 2) ** F(1, 2) if False else None
# exact sqrt via floats for the width only; verdict uses |z| >> 3
sig_tot_f = (float(sig_pred) ** 2 + float(SIG_EM) ** 2) ** 0.5
z_witherr = (float(em_M1_K) - float(EM_EXT)) / sig_tot_f
check("E15 sigma_pred from sigma_X: (5/3)*0.094 = 0.156667 MeV",
      round(float(sig_pred), 6) == 0.156667, float(sig_pred), 0.156667)
check("E16 combined sigma = 0.197869 MeV (pred + eps band)", round(sig_tot_f, 6) == 0.197869,
      round(sig_tot_f, 6), 0.197869)
check("E17 z WITH X error budget = -9.76 (docs/36 quoted -15.98 ignored sigma_X)",
      round(abs(z_witherr), 2) == 9.76, z_witherr, -9.76)
check("E18 verdict unchanged: |z| = 9.76 >> 3 => still KILL",
      abs(z_witherr) > 3, abs(z_witherr), ">3")

DATA["X"] = {"flag": str(X_FLAG), "sigma": str(SIG_X), "back": str(X_BACK),
             "back_minus_flag_rel": str(abs(X_BACK - X_FLAG) / X_FLAG),
             "verdict": "X=2.560 (FLAG24) canonical; 42667/16670 = rounding artifact"}
DATA["O_K"] = {"flag": str(O_K_flag), "back": str(O_K_back),
               "z_flag": str(z_flag), "z_back": str(z_back)}
DATA["ns"] = {"required_flag": str(ns_flag), "required_back": str(ns_back),
              "domain": list(NS_DOMAIN), "n_s_K": str(N_S_K), "O_M1_K": str(O_M1_K)}
DATA["z"] = {"direct_ns1": str(z_direct), "sigma_pred": str(sig_pred),
             "sigma_tot": str(sig_tot_f), "z_with_X_err": str(z_witherr)}

# ---------------- (2) M10 W4 : representation dependence
M_REP, X_REP = F("4.17e-4"), F("3.06e-1")
ratio = X_REP / M_REP
check("E19 M10 W4 ratio X/M representation = 733.8x (M10 quoted 734x)",
      round(float(ratio), 1) == 733.8, round(float(ratio), 1), 733.8)
em_req, x_obs = F("1.2667"), F("0.7823")
check("E20 M10 W4 em/|x| = 161.9% (M10 quoted 162%)",
      round(float(em_req / x_obs) * 100, 1) == 161.9, round(float(em_req / x_obs) * 100, 1), 161.9)
DATA["m10_w4"] = {"ratio": str(ratio), "em_over_x": str(em_req / x_obs)}

# ---------------- (3) M10 W6 : pool density
DATA["m10_w6"] = {"claim": "legal set 90 = 20.61% of F1 pool; channel incompleteness not slot density",
                  "verifiable_by_M4": False,
                  "note": "M10-internal (F1 pool size not in M4 corpus) -> recorded, NOT cited"}

npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 96)
print("RESULT: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 96)
DATA["summary"] = {"pass": npass, "total": len(CHECKS), "checks": CHECKS}

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "data", "40_recheck_m7_m10.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump(DATA, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON ->", out)
raise SystemExit(0 if npass == len(CHECKS) else 1)