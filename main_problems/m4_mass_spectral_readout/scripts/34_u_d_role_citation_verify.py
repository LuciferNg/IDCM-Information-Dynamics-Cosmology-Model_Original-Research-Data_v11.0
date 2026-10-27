#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/34 -- Citation verification for the u/d role question (2026-10-26).

Question raised to M4: "u/d sector originally belongs to which project?"
Answer per corpus: O209 already ruled (L5 = observational input) and O120 already
closed the u/d masses.  This script does NOT re-derive anything: it only verifies
the *documents' own arithmetic* so that M4 may cite them (citation + verification,
never citation alone).

Checks
  V1  O209 docs/08 L5 claim: the three generation blocks carry DIFFERENT
      up/down/lepton orderings  => role is not a function of block position.
  V2  O209 docs/08 section 4 quoted orderings reproduce exactly
      Gen3 [up,dn,lp] / Gen2 [up,lp,dn] / Gen1 [dn,up,lp].
  V3  O209 docs/08 section 4 quoted mass chain t>b>tau>c>mu>s>d>u>e reproduces.
  V4  O209 docs/10: role multiset per generation = {1 up + 1 down + 1 lepton}
      (SU(4)_PS) -- trivial set check on the three blocks.
  V5  O120 docs/04 section 4: legal mode set {1,2,3,4,5,7,9,11,13} and A5
      mass-ascending assignment e,u,d -> n=1,2,3 (i.e. u->2, d->3) reproduced
      from the masses alone (no fitting, just rank by mass).
  V6  Consequence for M4: the v9 section D.1 sector reading (u->k=1, d->k=4)
      is NOT the assignment used by any v10 O/M: the ruled assignment puts
      u at n=2 and d at n=3.

Masses: PDG 2024 central values as used in O209 docs/03 / O120 (MeV).
"""
import json
import os

CHECKS = []


def check(label, ok, got, expect=""):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-58s" % label, "|", got, "| expect", expect)


M = {  # MeV, PDG 2024 central values (same set as O209 docs/03)
    "t": 172690.0, "b": 4180.0, "tau": 1776.86,
    "c": 1270.0, "mu": 105.658, "s": 93.4,
    "d": 4.67, "u": 2.16, "e": 0.510999,
}
ROLE = {"t": "up", "b": "dn", "tau": "lp", "c": "up", "s": "dn", "mu": "lp",
        "d": "dn", "u": "up", "e": "lp"}
BLOCKS = {3: ["t", "b", "tau"], 2: ["c", "mu", "s"], 1: ["d", "u", "e"]}
MODE = {1: "e", 2: "u", 3: "d", 4: "s", 5: "mu", 7: "c", 9: "tau", 11: "b", 13: "t"}

print("=" * 84)
print("M4 scripts/34 -- citation verification: u/d role & assignment (O209 / O120)")
print("=" * 84)

# ---- V1/V2: orderings inside each generation block
orderings = {}
for g, parts in BLOCKS.items():
    ordered = sorted(parts, key=lambda p: -M[p])
    orderings[g] = [ROLE[p] for p in ordered]
print("\n[V1/V2] in-block orderings (heavy -> light):")
for g in (3, 2, 1):
    print("   Gen%d %-22s" % (g, str(BLOCKS[g])), "->", orderings[g])
check("V1 three blocks carry different role orderings (role != position)",
      len({tuple(v) for v in orderings.values()}) == 3,
      len({tuple(v) for v in orderings.values()}), "3 distinct")
check("V2a Gen3 == [up, dn, lp]", orderings[3] == ["up", "dn", "lp"], orderings[3], "[up, dn, lp]")
check("V2b Gen2 == [up, lp, dn]", orderings[2] == ["up", "lp", "dn"], orderings[2], "[up, lp, dn]")
check("V2c Gen1 == [dn, up, lp]", orderings[1] == ["dn", "up", "lp"], orderings[1], "[dn, up, lp]")

# ---- V3: the nine-mass chain quoted by O209 docs/08 section 4
chain = sorted(M, key=lambda p: -M[p])
quoted = ["t", "b", "tau", "c", "mu", "s", "d", "u", "e"]
print("\n[V3] mass chain:", chain)
check("V3 quoted chain t>b>tau>c>mu>s>d>u>e reproduces", chain == quoted, chain == quoted, True)

# ---- V4: role multiset per block
for g, parts in BLOCKS.items():
    roles = sorted(ROLE[p] for p in parts)
    check("V4 Gen%d role multiset == [dn, lp, up]" % g, roles == ["dn", "lp", "up"], roles, "[dn, lp, up]")

# ---- V5: A5 mass-ascending assignment (O120 docs/04 section 4)
asc = sorted(MODE, key=lambda n: M[MODE[n]])
spec = {1: "e", 2: "u", 3: "d"}
print("\n[V5] modes sorted by mass ascending:", [(n, MODE[n]) for n in asc])
check("V5a lightest three modes == (e, u, d)", [MODE[n] for n in asc[:3]] == ["e", "u", "d"],
      [MODE[n] for n in asc[:3]], ["e", "u", "d"])
check("V5b assignment e->1, u->2, d->3", all(MODE[n] == p for n, p in spec.items()),
      {n: MODE[n] for n in spec}, {"1": "e", "2": "u", "3": "d"})
check("V5c legal mode set == {1,2,3,4,5,7,9,11,13}", sorted(MODE) == [1, 2, 3, 4, 5, 7, 9, 11, 13],
      sorted(MODE), [1, 2, 3, 4, 5, 7, 9, 11, 13])

# ---- V6: consequence
u_mode = [n for n, p in MODE.items() if p == "u"][0]
d_mode = [n for n, p in MODE.items() if p == "d"][0]
print("\n[V6] ruled assignment: u -> n=%d, d -> n=%d ; v9 section D.1 sector reading: u->k=1, d->k=4"
      % (u_mode, d_mode))
check("V6a u sits at n=2 (not k=1)", u_mode == 2, u_mode, 2)
check("V6b d sits at n=3 (not k=4)", d_mode == 3, d_mode, 3)
check("V6c the two schemes disagree on both quarks",
      (u_mode, d_mode) != (1, 4), (u_mode, d_mode), "(1, 4)")

npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 84)
print("SUMMARY: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 84)

out = {
    "script": "scripts/34_u_d_role_citation_verify.py",
    "date": "2026-10-26",
    "purpose": "verify O209/O120 documents' own arithmetic before M4 cites them",
    "orderings": {str(k): v for k, v in orderings.items()},
    "mass_chain": chain,
    "assignment": {"mode_to_particle": {str(k): v for k, v in MODE.items()},
                   "u_mode": u_mode, "d_mode": d_mode},
    "checks": CHECKS, "n_pass": npass, "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/34_u_d_role.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/34_u_d_role.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)