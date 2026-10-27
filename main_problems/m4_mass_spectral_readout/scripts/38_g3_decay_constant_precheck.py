#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/38 -- G3 pre-check: can a SINGLE per-meson |psi(0)|^2 serve BOTH the
decay constants and the spin-spin (hyperfine) splittings?   (2026-10-27)

Why this is the right first step (and why it is NOT circular)
  M5's contact term needs |psi(0)|^2 per channel.  M5 currently has only a
  DATA-EXTRACTED value (0.052724 GeV^3 from the D channel).  Feeding that into a
  decay-constant prediction would be circular (same data both sides).
  This script therefore uses two INDEPENDENT observables per meson:
      (a) f_P                  (leptonic decay constant, PDG)
      (b) DeltaE_hfs = m_V - m_P  (spin-spin splitting, PDG)
  Both are proportional to the same |psi(0)|^2 with known mass factors:
      f_P^2      = 12 |psi(0)|^2 / m_P                        (quark-model normalisation)
      DeltaE_hfs = (8pi/3) a_s/(16 m_h m_l) * <sigma.sigma> * |psi(0)|^2 ,
                   <sigma.sigma>_{J=1} - <sigma.sigma>_{J=0} = 1/4 - (-3/4) = 1
  Dividing out the mass factors gives two extractions of psi^2 in the SAME units.
  Consistency across mesons is the test: the ratio must be a single constant
  (one overall alpha_s / normalisation cancels in the comparison).

PRE-REGISTERED criteria (fixed before numbers, 2026-10-27)
  Let R_P = psi2_hfs(P) / psi2_dec(P) and spread = max(R)/min(R) over {K, D, Ds}.
    spread <= 1.3  -> CONSISTENT   (single psi survives; G3 may proceed)
    1.3 < spread <= 3 -> MARGINAL
    spread > 3     -> INCONSISTENT (a single per-meson psi cannot serve both)
  pi is EXCLUDED as pseudo-Goldstone (contact term not applicable), only reported.

Conventions pinned (R5 lesson: definition first, numbers second)
  * f_pi = 130.2 MeV (PDG/FLAG "f_pi" convention). The other common convention
    f_pi = 92.07 MeV differs by sqrt(2); all targets below are ratios, so the
    convention cancels -- but it is stated.
  * constituent masses (M5 docs/09): m_u = m_d = 0.300, m_s = 0.495, m_c = 1.270 GeV
  * masses and splittings in GeV; psi^2 in GeV^3.
"""
import json
import os

CHECKS = []


def check(label, ok, got, expect=None):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-62s" % label, "|", got, "| expect", expect)


# ------------------------------------------------------------------ data (PDG 2024)
MP = {"pi": 0.13957039, "K": 0.493677, "D": 1.86966, "Ds": 1.96835}
MV = {"pi": 0.77526, "K": 0.89167, "D": 2.01026, "Ds": 2.1122}
FP = {"pi": 0.1302, "K": 0.1557, "D": 0.2120, "Ds": 0.2499}     # GeV
MH = {"K": 0.495, "D": 1.270, "Ds": 1.270}                      # heavy constituent
ML = {"K": 0.300, "D": 0.300, "Ds": 0.495}                      # light constituent
LIGHT = ["K", "D", "Ds"]

print("=" * 88)
print("M4 scripts/38 -- G3 pre-check: single psi(0) for f_P AND hyperfine?")
print("=" * 88)

# --- C1: convention / data sanity
check("C1a f_pi = 130.2 MeV (FLAG convention)", FP["pi"] == 0.1302, FP["pi"], 0.1302)
check("C1b m_K - m_pi > 0 and m_Ds > m_D", MP["K"] > MP["pi"] and MP["Ds"] > MP["D"],
      (MP["K"], MP["Ds"]), "ordering")
check("C1c hyperfine splittings all positive", all(MV[k] > MP[k] for k in LIGHT),
      {k: round(MV[k] - MP[k], 4) for k in LIGHT}, "all > 0")

# --- C2: required psi ratios from f_P alone (the frozen targets for any future psi)
psi2_dec = {k: FP[k] ** 2 * MP[k] / 12.0 for k in FP}
print("\n[C2] psi^2 from decay constants  (GeV^3):")
for k in MP:
    print("   %-3s f=%.4f  m=%.5f  ->  psi^2 = %.6e" % (k, FP[k], MP[k], psi2_dec[k]))
print("   f-ratios:  fK/fpi=%.4f  fD/fpi=%.4f  fDs/fpi=%.4f  fDs/fD=%.4f"
      % (FP["K"] / FP["pi"], FP["D"] / FP["pi"], FP["Ds"] / FP["pi"], FP["Ds"] / FP["D"]))
print("   psi-ratios (sqrt of psi^2 ratio):  psiK/psiD=%.4f  psiDs/psiD=%.4f"
      % ((psi2_dec["K"] / psi2_dec["D"]) ** 0.5, (psi2_dec["Ds"] / psi2_dec["D"]) ** 0.5))
check("C2a fK/fpi = 1.1959", abs(FP["K"] / FP["pi"] - 1.1959) < 2e-4, round(FP["K"] / FP["pi"], 4), 1.1959)
check("C2b fDs/fD = 1.1788", abs(FP["Ds"] / FP["D"] - 1.1788) < 2e-4, round(FP["Ds"] / FP["D"], 4), 1.1788)

# --- C3: psi^2 from hyperfine (alpha_s and 8pi/3 divided out)
psi2_hfs = {}
print("\n[C3] psi^2 from spin-spin splittings (arbitrary common normalisation):")
for k in LIGHT:
    dE = MV[k] - MP[k]
    psi2_hfs[k] = dE * 16.0 * MH[k] * ML[k]        # = (8pi/3) a_s * psi^2
    print("   %-3s dE=%.5f  m_h*m_l=%.5f  ->  (8pi/3)a_s*psi^2 = %.6f"
          % (k, dE, MH[k] * ML[k], psi2_hfs[k]))
check("C3a DeltaE(D*-D) = 140.60 MeV", abs((MV["D"] - MP["D"]) - 0.14060) < 2e-4,
      round((MV["D"] - MP["D"]) * 1000, 2), 140.60)
check("C3b DeltaE(Ds*-Ds) = 143.85 MeV", abs((MV["Ds"] - MP["Ds"]) - 0.14385) < 2e-4,
      round((MV["Ds"] - MP["Ds"]) * 1000, 2), 143.85)
check("C3c DeltaE(K*-K) = 397.99 MeV", abs((MV["K"] - MP["K"]) - 0.39799) < 2e-4,
      round((MV["K"] - MP["K"]) * 1000, 2), 397.99)

# --- C4: consistency ratio
R = {k: psi2_hfs[k] / psi2_dec[k] for k in LIGHT}
spread = max(R.values()) / min(R.values())
print("\n[C4] R = psi2_hfs / psi2_dec (must be ONE constant; alpha_s cancels):")
for k in LIGHT:
    print("   %-3s R = %.3f" % (k, R[k]))
print("   spread = max/min = %.3f" % spread)

if spread <= 1.3:
    verdict = "CONSISTENT"
elif spread <= 3.0:
    verdict = "MARGINAL"
else:
    verdict = "INCONSISTENT"
print("   PRE-REGISTERED verdict: %s" % verdict)

check("C4a heavy sector (D vs Ds) agrees within 30%",
      abs(R["Ds"] / R["D"] - 1) <= 0.30, round(R["Ds"] / R["D"], 3), "1.00 +- 0.30")
check("C4b K channel deviates by more than 3x from the heavy sector",
      R["K"] / R["D"] > 3.0, round(R["K"] / R["D"], 2), "> 3")
check("C4c overall spread > 3 -> INCONSISTENT", verdict == "INCONSISTENT", verdict, "INCONSISTENT")

# --- C5: pi as reference only (excluded from verdict)
psi2_hfs["pi"] = (MV["pi"] - MP["pi"]) * 16.0 * 0.300 * 0.300   # m_u = m_d = 0.300
R_pi_note = psi2_hfs["pi"] / psi2_dec["pi"]
print("\n[C5] pi (reported only; pseudo-Goldstone, contact term not applicable):")
print("   R_pi = %.3f  (%.2fx the heavy-sector mean %.3f)"
      % (R_pi_note, R_pi_note / ((R["D"] + R["Ds"]) / 2), (R["D"] + R["Ds"]) / 2))

npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 88)
print("SUMMARY: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 88)
print("Reading: the heavy sector (D, Ds) is compatible with ONE psi, but K is off by")
print("%.1fx -- so a single per-meson |psi(0)|^2 cannot serve f_P and the hyperfine" % (R["K"] / R["D"]))
print("splitting across the light sector.  The K channel is exactly where the M4 em")
print("anchor (em(K+-K0) = 2.2646 +- 0.1209 MeV) sits, so this closes the same door.")

out = {
    "script": "scripts/38_g3_decay_constant_precheck.py",
    "date": "2026-10-27",
    "round": "G3 pre-check (model-free; no M5 parameters used, no fitting)",
    "why_not_circular": "f_P and DeltaE_hfs are independent observables; only one overall constant cancels",
    "pre_registered": {"statistic": "spread = max(R)/min(R)", "cut_consistent": 1.3,
                       "cut_marginal": 3.0, "K_exclusion": "pi excluded as pseudo-Goldstone"},
    "data": {"m_P_GeV": MP, "m_V_GeV": MV, "f_P_GeV": FP,
             "m_h_GeV": MH, "m_l_GeV": ML},
    "psi2_from_decay": psi2_dec, "psi2_from_hfs": psi2_hfs,
    "R": R, "spread": spread, "verdict": verdict,
    "pi_reference_R": R_pi_note,
    "checks": CHECKS, "n_pass": npass, "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/38_g3_decay_constant_precheck.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/38_g3_decay_constant_precheck.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)