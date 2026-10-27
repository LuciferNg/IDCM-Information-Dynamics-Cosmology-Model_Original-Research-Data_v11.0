#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/36b -- R5 (v2): K-channel EXTERNAL JUDGEMENT of the pre-registered
prediction (2026-10-27).  External test only: no refitting, no parameter touching.

WHY v2 (honest history -- scripts/36 v1 is KEPT unchanged, do not delete)
  v1 tried two GUESSED readings of the Dashen-violation parameter (linear-in-mass,
  and a "1 - eps" squared-mass reading) because the definition had not been pinned.
  Both guesses were WRONG.  The definition is explicit in FLAG Review 2024
  (arXiv:2411.04268), Sec. 4.1.3, Eqs. (38)-(39):

      eps = (dM2_K - dM2_pi)_gamma / dM2_pi ,
      dM2_pi = M_pi+^2 - M_pi0^2 ,  dM2_K = M_K+^2 - M_K0^2 ,
  and since G-parity forbids O(m_u-m_d) corrections to dM2_pi,
      eps = (dM2_K)_gamma / dM2_pi - 1   =>   (dM2_K)_gamma = (1 + eps) dM2_pi .

  So the kaon squared EM splitting is LARGER than the pion's by (1+eps) ~ 1.78:
  em(K+ - K0) ~ 2.26 MeV -- not ~0.3 MeV and not ~1.0 MeV.  v1's "1 - eps"
  reading had the sign inverted, and the verdict inverts with it.

External anchor
  MILC 18 (Phys.Rev.D 99, 034503 (2019), arXiv:1807.05556):
      eps = 0.78 (1)_stat (+8 / -11)_syst  ->  eps in [0.67, 0.86]

Pre-registered prediction (frozen 2026-10-25 in M4 docs/29, before this comparison)
      em(K+ - K0) = +0.3327 MeV
  and the type-independent envelope of docs/30
      em(K0 - K+ ) in [-0.52, +0.55] MeV .

PRE-REGISTERED criteria (unchanged from v1, fixed before numbers):
  |z| < 1 PASS ; 1-2 TENSION ; 2-3 STRAIN ; >=3 KILL
"""
import json
import os

CHECKS = []


def check(label, ok, got, expect=None):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-64s" % label, "|", got, "| expect", expect)


MPI_P, MPI_0 = 139.57039, 134.9768       # MeV, PDG 2024
MK_P, MK_0 = 493.677, 497.611            # MeV, PDG 2024
EPS, EPS_STAT, EPS_UP, EPS_DN = 0.78, 0.01, 0.08, 0.11

PRED = 0.3327                            # docs/29 (2026-10-25), frozen
ENV_LO, ENV_HI = -0.52, 0.55             # docs/30, quoted for em(K0 - K+)

print("=" * 90)
print("M4 scripts/36b -- R5 K-channel external judgement (FLAG Eq.39 definition pinned)")
print("=" * 90)

dpi2 = MPI_P ** 2 - MPI_0 ** 2
dk2_total = MK_P ** 2 - MK_0 ** 2
mk_sum = MK_P + MK_0

check("S1 dM2_pi = 1261.15 MeV^2", abs(dpi2 - 1261.15) < 0.05, round(dpi2, 3), 1261.15)
check("S2 dM2_K total = -3899.7 MeV^2 (K0 heavier)", abs(dk2_total + 3899.7) < 0.5,
      round(dk2_total, 2), -3899.7)
check("S3 M_K+ - M_K0 = -3.934 MeV", abs((MK_P - MK_0) + 3.934) < 5e-3,
      round(MK_P - MK_0, 4), -3.934)

eps_lo, eps_hi = EPS - EPS_DN, EPS + EPS_UP          # 0.67 .. 0.86


def em_kaon_from_eps(eps):
    return (1.0 + eps) * dpi2 / mk_sum


c, lo, hi = em_kaon_from_eps(EPS), em_kaon_from_eps(eps_lo), em_kaon_from_eps(eps_hi)
sig = (hi - lo) / 2.0
print("\n(dM2_K)_gamma = (1+eps) dM2_pi = 1.78 x %.2f = %.2f MeV^2" % (dpi2, (1 + EPS) * dpi2))
print("em(K+ - K0)   = %.4f MeV   68%% range [%.4f, %.4f]   sigma_sym = %.4f" % (c, lo, hi, sig))
print("reference: Dashen (eps=0) -> em(K+ - K0) = %.4f MeV" % em_kaon_from_eps(0.0))

z_pred = (PRED - c) / sig
inside_env = ENV_LO <= -c <= ENV_HI
z_env = 0.0 if inside_env else (max(-c - ENV_HI, ENV_LO - (-c)) / sig)
print("\npre-registered prediction %+0.4f MeV -> z = %+.2f sigma ; ratio pred/em = %.3f"
      % (PRED, z_pred, PRED / c))
print("envelope: -em = %+.4f vs [%+.2f, %+.2f] -> %s (z_edge = %.2f)"
      % (-c, ENV_LO, ENV_HI, "inside" if inside_env else "OUTSIDE", z_env))


def cat(z):
    az = abs(z)
    return "PASS" if az < 1 else ("TENSION" if az < 2 else ("STRAIN" if az < 3 else "KILL"))


v_pred, v_env = cat(z_pred), cat(z_env)
print("\nCATEGORY  prediction: %-8s (z=%+.2f)   envelope: %-8s (z_edge=%.2f)"
      % (v_pred, z_pred, v_env, z_env))

check("V1 em(K+ - K0) from MILC eps = 2.2646 MeV", abs(c - 2.2646) < 2e-3, round(c, 4), 2.2646)
check("V2 pre-registered prediction KILLED (|z| >= 3)", v_pred == "KILL", v_pred, "KILL")
check("V3 docs/30 envelope KILLED as well", v_env == "KILL", v_env, "KILL")
check("V4 prediction is ~6.8x too small", abs(PRED / c - 0.1469) < 5e-3, round(PRED / c, 4), 0.1469)
check("V5 Dashen-scale value 1.2722 MeV reproduced", abs(em_kaon_from_eps(0.0) - 1.2722) < 2e-3,
      round(em_kaon_from_eps(0.0), 4), 1.2722)

npass = sum(1 for x in CHECKS if x["ok"])
print("\n" + "=" * 90)
print("SUMMARY: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 90)
print("Reading: the <O>-corrected kaon prediction (+0.333 MeV) is excluded by lattice")
print("QCD+QED at ~%.0f sigma, and so is the type-independent envelope of docs/30." % abs(z_pred))
print("The lattice confirms instead the plain Dashen-scale splitting 2.26 MeV -- exactly")
print("the number the <O> line was introduced to remove.")

out = {
    "script": "scripts/36b_em_kaon_external_judgement.py (v2)",
    "date": "2026-10-27",
    "round": "M4 R5 (external judgement, no fitting)",
    "definition_source": "FLAG Review 2024, arXiv:2411.04268, Sec.4.1.3 Eqs.(38)-(39)",
    "definition": "eps = (dM2_K)_gamma/dM2_pi - 1 => (dM2_K)_gamma = (1+eps) dM2_pi",
    "v1_note": "scripts/36 v1 guessed two readings (linear; 1-eps) -- both wrong, eps sign inverted",
    "external_anchor": {"paper": "MILC 18, PRD 99 034503 (arXiv:1807.05556)",
                        "eps": EPS, "stat": EPS_STAT, "syst_up": EPS_UP, "syst_dn": EPS_DN,
                        "range": [eps_lo, eps_hi]},
    "masses_MeV": {"mpi_p": MPI_P, "mpi_0": MPI_0, "mK_p": MK_P, "mK_0": MK_0},
    "dM2_pi": dpi2, "dM2_K_total": dk2_total,
    "em_K": {"central": c, "lo": lo, "hi": hi, "sigma_sym": sig,
             "dashen_eps0": em_kaon_from_eps(0.0)},
    "pre_registered": {"prediction": PRED, "envelope": [ENV_LO, ENV_HI],
                       "source": "docs/29, docs/30 (2026-10-25)"},
    "verdict": {"prediction": v_pred, "z_pred": z_pred, "envelope": v_env, "z_env": z_env},
    "checks": CHECKS, "n_pass": npass, "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/36_em_kaon_external.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/36_em_kaon_external.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)