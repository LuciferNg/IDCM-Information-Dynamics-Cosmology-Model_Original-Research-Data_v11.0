#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""S3 支線：味對標籤身份核實（修正：946 = C(44,2)，群作用 mod 88）
   命題（corpus）：946 價對在 D 作用下得 43 軌、Δ=b−a 完備、Σ=a+b 為軌內座標
   本探針：獨立枚舉 (k,s) ∈ Z_88 × {0,1} 作用 x ↦ ±x + k (mod 88)，核軌道數／大小／不變量
"""
import json
from itertools import combinations

MOD = 88
ST = list(range(1, 45))                       # 44 站（1..44）
pairs = list(combinations(ST, 2))
pset = set(pairs)


def act(k, s, p):
    o = [((-x) % MOD if s else x) for x in p]
    o = [(y + k) % MOD for y in o]
    a, b = sorted(o)
    return (a, b)


seen, orbits = set(), []
for p in pairs:
    if p in seen:
        continue
    st, orb = [p], []
    seen.add(p)
    while st:
        q = st.pop()
        orb.append(q)
        for k in range(MOD):
            for s in (0, 1):
                r = act(k, s, q)
                if r in pset and r not in seen:
                    seen.add(r)
                    st.append(r)
    orbits.append(sorted(orb))

sizes = sorted(len(o) for o in orbits)
delta_of = lambda p: (p[1] - p[0]) % MOD
delta_complete = len({delta_of(o[0]) for o in orbits}) == len(orbits)
delta_in_orbit = all(len({delta_of(p) for p in o}) == 1 for o in orbits)
size_eq_delta = all(len(o) == 44 - delta_of(o[0]) for o in orbits)
sigma_sets = sorted({len({(p[0] + p[1]) % MOD for p in o}) for o in orbits})

print("價對 946｜軌道", len(orbits), "｜大小", sizes)
print("Σ_{Δ=1..43} Δ =", sum(range(1, 44)))
print("Δ 完備:", delta_complete, "｜Δ 軌內恆定:", delta_in_orbit, "｜軌大小 == 44−Δ:", size_eq_delta)
print("Σ mod 88 之軌內取值數:", sigma_sets)

out = {"n_pairs": len(pairs), "n_orbits": len(orbits), "orbit_sizes": sizes,
       "sum_sizes": sum(sizes), "sum_delta_1_43": sum(range(1, 44)),
       "delta_complete": delta_complete, "delta_constant_in_orbit": delta_in_orbit,
       "orbit_size_equals_44_minus_delta": size_eq_delta, "sigma_within_orbit_card": sigma_sets,
       "orbit_table": {str(delta_of(o[0])): len(o) for o in orbits},
       "verdict": ("946 = Σ_{Δ=1}^{43}(44−Δ)；每 Δ 恰一軌、軌大小 = 44−Δ、Σ 軌內取 44−Δ 值 ⟹ 味對標籤 = Δ ∈ 1..43（站 = 44）"
                   "⟹ 味對標籤集（43）與站集（44）基數不同 ⟹ 味站指派 ill-posed"
                   if len(orbits) == 43 and size_eq_delta and delta_complete else "需重核")}
with open("/home/wsl/idcm_ckpt/e8/mrB/s3_pair_orbits.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("判詞:", out["verdict"])
