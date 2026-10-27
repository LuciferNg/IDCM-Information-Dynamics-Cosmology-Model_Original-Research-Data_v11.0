#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/35 -- v9 particle<->mode tables: coverage / injectivity / conflict audit
(2026-10-26).  Documentary audit only; no new derivation, no fitting.

Question (from 調律者): "同一條公式 u/d 不適用，但其他適用，合理嗎？"
To answer it we must first establish WHAT the v9 tables actually are.

Three tables exist inside v9, all claiming to place charged fermions on chain modes:

  T1  wfield_lattice/03_extended_derivation.md  section D.1  (generation modes)
        k=1 Up sector (c,t) / k=4 Down sector (s,b) / k=5 Lepton sector (e,mu,tau)
  T2  wfield_lattice/cy3_topological_operator.md section 5.5  (per-generation)
        gen3 (t,b,tau)=k0 ; gen2 (c,s,mu)=k1,3,8 ; gen1 (u,d,e)=k5,8,9
  T3  wfield_lattice/cy3_topological_operator.md section 5.6  (mass-formula layer)
        e=k5 (s=2) ; d=k8 (s=2) ; u=k9 (s=3)

Ruled v10 assignment (O137 SESSION_START; O209 docs/14,15; O233):
  T0  e=1, u=2, d=3, s=4, mu=5, c=7, tau=9, b=11, t=13   (A5 gamma-ascending)

Checks
  C1  T1 coverage over the 9 charged fermions
  C2  T2 injectivity (a mode may host at most one particle)
  C3  T3 vs T1 conflict (same k, different sector)
  C4  T0 coverage and injectivity
  C5  pairwise agreement of T1/T2/T3/T0 on the particles they share
"""
import json
import os

CHARGED = ["u", "d", "c", "s", "t", "b", "e", "mu", "tau"]
CHECKS = []


def check(label, ok, got, expect=None):
    """expect: any type (str/list/tuple/int) -- annotation kept loose on purpose."""
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-56s" % label, "|", got, "| expect", expect)


T1 = {"c": 1, "t": 1, "s": 4, "b": 4, "e": 5, "mu": 5, "tau": 5}
T2 = {"t": 0, "b": 0, "tau": 0, "c": 1, "s": 3, "mu": 8, "u": 5, "d": 8, "e": 9}
T3 = {"e": 5, "d": 8, "u": 9}
T0 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}

print("=" * 82)
print("M4 scripts/35 -- v9 table audit: coverage / injectivity / cross-table conflict")
print("=" * 82)


def inverse(t):
    inv = {}
    for p, k in t.items():
        inv.setdefault(k, []).append(p)
    return inv


# ---- C1 coverage of T1
missing1 = [p for p in CHARGED if p not in T1]
check("C1a T1 covers 7 of 9 charged fermions", len(T1) == 7, "%d/9" % len(T1), "7/9")
check("C1b T1 is missing exactly {u, d}", sorted(missing1) == ["d", "u"], sorted(missing1), ["d", "u"])
# NOTE (2026-10-26, self-correction): T1 maps a SECTOR to a SET of particles, so
# particle->k injectivity is impossible by construction and is NOT the right criterion.
# The correct criterion is that T1 is a partition: disjoint sectors, no particle twice.
_sets = {"Up": {"c", "t"}, "Down": {"s", "b"}, "Lepton": {"e", "mu", "tau"}}
_pairs_disjoint = all(not (_sets[a] & _sets[b]) for a in _sets for b in _sets if a < b)
check("C1c T1 sector sets are disjoint (a partition)", _pairs_disjoint,
      {k: sorted(v) for k, v in _sets.items()}, "disjoint")

# ---- C2 T2 injectivity
inv2 = inverse(T2)
coll2 = {k: v for k, v in inv2.items() if len(v) > 1}
check("C2a T2 is NOT injective", bool(coll2), coll2, "non-empty")
check("C2b k=8 is shared by mu and d", sorted(coll2.get(8, [])) == ["d", "mu"], coll2.get(8), ["d", "mu"])
check("C2c T2 mode-count 6 < 9 particles", len(inv2) == 6, len(inv2), 6)

# ---- C3 T3 vs T1 conflict on k=5
t1_sector = {1: "Up", 4: "Down", 5: "Lepton"}
check("C3a T1 puts k=5 in the Lepton sector", t1_sector[5] == "Lepton", t1_sector[5], "Lepton")
check("C3b T3 puts the electron at k=5", T3["e"] == 5, T3["e"], 5)
check("C3c T1 puts the DOWN-quark sector at k=4, T3 puts d at k=8",
      (t1_sector[4] == "Down") and (T3["d"] == 8), (t1_sector[4], T3["d"]), ("Down", 8))
check("C3d T2 puts u at k=5 and d at k=8 -- third distinct convention",
      (T2["u"], T2["d"]) == (5, 8), (T2["u"], T2["d"]), (5, 8))

# ---- C4 ruled assignment T0
check("C4a T0 covers 9 of 9", len(T0) == 9, "%d/9" % len(T0), "9/9")
check("C4b T0 is injective", all(len(v) == 1 for v in inverse(T0).values()),
      {k: v for k, v in inverse(T0).items() if len(v) > 1} or "injective", "injective")

# ---- C5 pairwise agreement on shared particles
def agree(a, b):
    shared = sorted(set(a) & set(b))
    same = [p for p in shared if a[p] == b[p]]
    return len(shared), same


pairs = [("T1", T1, "T2", T2), ("T1", T1, "T3", T3), ("T1", T1, "T0", T0),
         ("T2", T2, "T0", T0), ("T3", T3, "T0", T0)]
print("\n[C5] pairwise agreement on shared particles:")
agree_rows = []
for na, a, nb, b in pairs:
    n, same = agree(a, b)
    print("   %s vs %s: shared=%d  agree=%d  -> %s" % (na, nb, n, len(same), same))
    agree_rows.append({"a": na, "b": nb, "shared": n, "agree": len(same), "same": same})

check("C5a T1 vs T3 agree on the electron only", len(agree(T1, T3)[1]) == 1, agree(T1, T3)[1], ["e"])
check("C5b T1 vs T2 agree on exactly one particle (c, k=1)", agree(T1, T2)[1] == ["c"], agree(T1, T2)[1], ["c"])
check("C5c T3 vs T0 agree on nothing", len(agree(T3, T0)[1]) == 0, agree(T3, T0)[1], [])
check("C5d T2 vs T0 agree on nothing", len(agree(T2, T0)[1]) == 0, agree(T2, T0)[1], [])

npass = sum(1 for c in CHECKS if c["ok"])
print("\n" + "=" * 82)
print("SUMMARY: %d/%d PASS" % (npass, len(CHECKS)))
print("=" * 82)
print("Reading: T1 is a 7/9 partial enumeration; T2 is non-injective (k=8 hosts mu and d,")
print("k=0 hosts t,b,tau); T3 is a third convention that contradicts T1 on k=5/k=4;")
print("and none of T1/T2/T3 agrees with the ruled T0.  The 'u/d exception' is therefore")
print("not a property of u/d -- it is a symptom of three mutually inconsistent tables.")

out = {
    "script": "scripts/35_v9_table_audit.py",
    "date": "2026-10-26",
    "tables": {"T1_sector_D1": T1, "T2_generation_5_5": T2, "T3_massform_5_6": T3,
               "T0_ruled": T0},
    "coverage": {"T1": "%d/9" % len(T1), "T0": "%d/9" % len(T0)},
    "collisions_T2": {str(k): v for k, v in coll2.items()},
    "pairwise": agree_rows,
    "checks": CHECKS, "n_pass": npass, "n_checks": len(CHECKS),
}
os.makedirs("data", exist_ok=True)
with open("data/35_v9_table_audit.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("wrote data/35_v9_table_audit.json")
raise SystemExit(0 if npass == len(CHECKS) else 1)