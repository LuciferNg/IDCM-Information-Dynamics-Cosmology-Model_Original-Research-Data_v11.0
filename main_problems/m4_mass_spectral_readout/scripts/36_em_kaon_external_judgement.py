#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/36 -- R5: K-channel external judgement of the pre-registered prediction
(2026-10-27).  EXTERNAL TEST ONLY -- no refitting, no parameter touching.

Situation
  M4 `docs/29` (2026-10-25) PRE-REGISTERED the sharpest prediction of the <O> line:
      em(K+ - K0) = +0.3327 MeV        (orientation: dM = M(high Q) - M(low Q))
  and `docs/30` recorded the type-independent envelope
      em(K0 - K+ ) in [-0.52, +0.55] MeV .
  Up to now this was marked "waiting for an independent lattice QED kaon value".
  Such a value exists: MILC (Phys.Rev.D 99, 034503 (2019), arXiv:1807.05556)
  quotes the Dashen-violation parameter for the K+ - K0 EM splitting,
      epsilon = 0.78 (1)_stat (+8 / -11)_syst   -> epsilon in [0.67, 0.86].

The comparison hinges on ONE convention question: Dashen's theorem is normally
stated for SQUARED masses (chiPT: M^2 linear in quark mass).  Two readings:

  Branch A (linear masses):   em(K+ - K0) = (1 - eps) (m_pi+ - m_pi0)
  Branch B (squared masses):  em(K+ - K0) = (1 - eps) (M_pi+^2 - M_pi0^2) / (M_K+ + M_K0)

Both are computed.  A verdict is only issued as PASS/FAIL when BOTH branches agree;
otherwise the honest verdict is "convention not pinned" (CONV).

PRE-REGISTERED criteria (fixed before any number was computed, 2026-10-27):
  |z| < 1        -> PASS      (support)
  1 <= |z| < 2   -> TENSION   (not decisive)
  2 <= |z| < 3   -> STRAIN
  |z| >= 3       -> KILL      (falsified)
  branches disagree on category -> CONV (no claim allowed)

Inputs are the published PDG 2024 central values and the published MILC epsilon;
nothing here is fitted.
"""
import json
import math
import os

CHECKS = []


def check(label, ok, got, expect=None):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-62s" % label, "|", got, "| expect", expect)


# ------------------------------------------------------------------ inputs
MPI_P, MPI_0 = 139.57039, 134.9768      # MeV, PDG 2024
MK_P, MK_0 = 493.677, 497.611           # MeV, PDG 2024
EPS, EPS_STAT, EPS_SYST_UP, EPS_SYST_DN = 0.78, 0.01, 0.08, 0.11  # MILC 2019

PRED = 0.3327                            # docs/29, pre-registered 2026-10-25
ENV_LO, ENV_HI = -0.52, 0.55             # docs/30 envelope for em(K0 - K+)

print("=" * 88)
print("M4 scripts/36 -- R5 K-channel external judgement (MILC epsilon) ")
print("=" * 88)

# ------------------------------------------------------------ sanity of the anchors
dpi = MPI_P - MPI_0
dpi2 = MPI_P ** 2 - MPI_0 ** 2
dk2 = MK_P ** 2 - MK_0 ** 2
check("S1 pi+ - pi0 = 4.59359 MeV (PDG 2024)", abs(dpi - 4.59359) < 5e-4, round(dpi, 5), 4.59359)
# NOTE (self-correction 2026-10-27): 1261.2, not 1266.6 -- my arithmetic slip, corrected
# against 139.57039^2 - 134.97680^2 = 1261.15 MeV^2.
check("S2 M_pi+^2 - M_pi0^2 = 1261.2 MeV^2", abs(dpi2 - 1261.2) < 0.5, round(dpi2, 2), 1261.2)
check("S3 M_K+^2 - M_K0^2 = -3899.8 MeV^2 (negative: K0 heavier)",
      dk2 < 0 and abs(dk2 + 3899.8) < 5.0, round(dk2, 1), -3899.8)

eps_lo, eps_hi = EPS - EPS_SYST_DN, EPS + EPS_SYST_UP   # 0.67 .. 0.86


def em_branch_A(eps):
    """linear masses"""
    return (1.0 - eps) * dpi


def em_branch_B(eps):
    """squared masses (standard chiPT reading of Dashen)"""
    return (1.0 - eps) * dpi2 / (MK_P + MK_0)


for name, fn in (("A linear", em_branch_A), ("B squared", em_branch_B)):
    lo, c, hi = fn(eps_hi), fn(EPS), fn(eps_lo)     # eps up -> em down
    print("\n[%s] em(K+ - K0) = %.4f MeV   68%% range [%.4f, %.4f]   sigma_sym = %.4f"
          % (name, c, lo, hi, (hi - lo) / 2.0))
    print("      naive Dashen (eps=0): %.4f MeV" % fn(0.0))

res = {}
for tag, fn in (("A_linear", em_branch_A), ("B_squared", em_branch_B)):
    lo, c, hi = fn(eps_hi), fn(EPS), fn(eps_lo)
    sig = (hi - lo) / 2.0
    z_pred = (c - PRED) / sig
    # envelope: em(K0-K+) = -em(K+K0)  -> compare -c against [-0.52, +0.55]
    z_env = (max(-c - ENV_HI, ENV_LO - (-c), 0.0)) / sig
    z_dashen = (c - fn(0.0)) / sig
    res[tag] = {"em": c, "lo": lo, "hi": hi, "sig": sig,
                "z_pred": z_pred, "z_env": z_env, "z_dashen_vs_lattice": z_dashen,
                "em_dashen": fn(0.0)}
    print("\n%s: em = %.4f +- %.4f MeV | pre-registered pred %+.4f -> z = %+.3f"
          % (tag, c, sig, PRED, z_pred))
    print("    envelope test: -em = %+.4f vs [%+.2f, %+.2f] -> %s (z_edge = %.3f)"
          % (-c, ENV_LO, ENV_HI,
             "inside" if ENV_LO <= -c <= ENV_HI else "OUTSIDE", z_env))


def cat(z):
    z = abs(z)
    if z < 1.0:
        return "PASS"
    if z < 2.0:
        return "TENSION"
    if z < 3.0:
        return "STRAIN"
    return "KILL"


print("\n" + "-" * 88)
print("PRE-REGISTERED category per branch:")
cats = {}
for tag in ("A_linear", "B_squared"):
    c1 = cat(res[tag]["z_pred"])
    c2 = cat(res[tag]["z_env"])
    cats[tag] = (c1, c2)
    print("   %-10s prediction: %-8s (z=%+.3f) | envelope: %-8s (z_edge=%.3f)"
          % (tag, c1, res[tag]["z_pred"], c2, res[tag]["z_env"]))

agree = (cats["A_linear"] == cats["B_squared"])
verdict = cats["A_linear"][0] if agree else "CONV"
print("\nBranches agree: %s -> joint verdict on the pre-registered single-form prediction: %s"
      % (agree, verdict))
print("Interpretation guide:")
print("  PASS  -> the <O> line's sharpest prediction is CONFIRMED by external lattice QED")
print("  CONV  -> the answer turns on the epsilon definition; must be pinned from the")
print("           primary source before any claim (no claim allowed now)")
print("  TENSION/KILL -> see magnitude")

check("V1 branch A gives em(K+-K0) near 1.01 MeV", abs(res["A_linear"]["em"] - 1.0106) < 5e-3,
      round(res["A_linear"]["em"], 4), 1.0106)
check("V2 branch B gives em(K+-K0) near 0.281 MeV", abs(res["B_squared"]["em"] - 0.2806) < 5e-3,
      round(res["B_squared"]["em"], 4), 0.2806)
check("V3 prediction 0.3327 lies INSIDE branch B 68% range",
      res["B_squared"]["lo"] <= PRED <= res["B_squared"]["hi"],
      (round(res["B_squared"]["lo"], 4), round(res["B_squared"]["hi"], 4)), "0.179 .. 0.422")
check("V4 prediction 0.3327 lies OUTSIDE branch A 68% range",
      not (res["A_linear"]["lo"] <= PRED <= res["A_linear"]["hi"]),
      (round(res["A_linear"]["lo"], 4), round(res["A_linear"]["hi"], 4)), "0.643 .. 1.516")
check("V5 branches disagree -> verdict must be CONV", agree is False, verdict, "CONV")
check("V6 naive Dashen is excluded by the lattice anchor in BOTH branches",
      abs(res["A_linear"]["z_dashen_vs_lattice"]) > 3 and abs(res["B_squared"]["z_dashen_vs_lattice"]) > 3,
      (round(res["A_linear"]["z_dashen_vs_lattice"], 2), round(res["B_squared"]["z_dashen_vs_lattice"], 2)),
      ">3 both")

npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 88)
print("SUMMARY: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 88)

out = {
    "script": "scripts/36_em_kaon_external_judgement.py",
    "date": "2026-10-27",
    "round": "M4 R5 (external judgement; no fitting)",
    "pre_registered": {"prediction_MeV": PRED, "envelope_MeV": [ENV_LO, ENV_HI],
                       "source": "docs/29 (2026-10-25), docs/30 (2026-10-25)",
                       "criteria": {"<1": "PASS", "1-2": "TENSION", "2-3": "STRAIN", ">=3": "KILL",
                                    "disagree": "CONV"}},
    "external_anchor": {"paper": "MILC, Phys.Rev.D 99, 034503 (2019), arXiv:1807.05556",
                        "epsilon": EPS, "stat": EPS_STAT,
                        "syst_up": EPS_SYST_UP, "syst_dn": EPS_SYST_DN,
                        "range": [eps_lo, eps_hi]},
    "masses_MeV": {"mpi_p": MPI_P, "mpi_0": MPI_0, "mK_p": MK_P, "mK_0": MK_0},
    "branches": res, "categories": {k: list(v) for k, v in cats.items()},
    "verdict": verdict, "checks": CHECKS, "n_pass": npass, "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/36_em_kaon_external.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/36_em_kaon_external.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)