#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""算子化 III：律乙（符號同質）之推導核實
   引理：對 v ∈ [1,44]⁴，σ_S = v1 + s2v2 + s3v3 + s4v4 ∈ [−131,132]，
        σ_S ≡ 0 (mod 90) ⟹ σ_S ∈ {−90, 0, 90}
   推導：任兩條同時成立之條件，其 σ 差 = 2·(±v_a ± v_b) ∈ {0,±90,±180}
        ⟹ (±v_a ± v_b) ∈ {0,±45,±90}，而範圍 [−88,88] ⟹ ∈ {0,±45}
   本探針：全枚舉核實 ①σ 值域 ②兩條件同時成立時 σ 差集合 ③同號性 ④|mask| ≤ 3
"""
import itertools
import json
from collections import Counter

MOD = 90
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
TPROD = [s[1] * s[2] * s[3] for s in SIG]        # = CS

sigma_vals, diffs, mixed, maxmask = Counter(), Counter(), 0, 0
for rem in range(1, 45):
    for a3 in itertools.combinations_with_replacement(range(1, 45), 3):
        v = (rem,) + a3
        sig = [v[0] * s[0] + v[1] * s[1] + v[2] * s[2] + v[3] * s[3] for s in SIG]
        sat = [i for i in range(8) if sig[i] % MOD == 0]
        maxmask = max(maxmask, len(sat))
        for i in sat:
            sigma_vals[sig[i]] += 1
        for i, j in itertools.combinations(sat, 2):
            diffs[sig[i] - sig[j]] += 1
            if TPROD[i] != TPROD[j]:
                mixed += 1

print("① 成立條件之 σ 值集合:", dict(sigma_vals))
print("② 兩條件同時成立之 σ 差集合:", dict(diffs))
print("③ 混號對（符號積不同）計數:", mixed)
print("④ max |mask| =", maxmask)

out = {"sigma_value_set": {str(k): v for k, v in sorted(sigma_vals.items())},
       "sigma_diff_set": {str(k): v for k, v in sorted(diffs.items())},
       "n_mixed_pairs": mixed, "max_mask_card": maxmask,
       "lemma": "σ_S ∈ {−90,0,90}；σ_S−σ_T = 2(±v_a±v_b) ∈ {0,±90,±180} ⟹ ±v_a±v_b ∈ {0,±45}",
       "verdict": ("律乙為定理：σ 差只取 {0,90,−90} ⟹ 同號；全枚舉混號對 0" if mixed == 0 else "有混號對，律乙非定理")}
with open("/home/wsl/idcm_ckpt/e8/mrB/op_sign_homogeneity_proof.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("判詞:", out["verdict"])
