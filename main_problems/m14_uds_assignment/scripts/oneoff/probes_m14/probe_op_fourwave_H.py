#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""算子化 IV：四波頂點 H = (1/4!)[:A⁴: + 6Σ:A²: + 3N_d] 之算符核實
   A = Σ_{n=1}^{44} (S^n − 1)  ∈ End(ℓ²(Z_90))
   核實：⟨A²⟩、⟨A⁴⟩、Wick 分解、以及係數 6／3N_d=132 是否為結構常數
   （全部整數精確：90×90 整數矩陣）
"""
import json

import numpy as np

MOD = 90


def shift(n):
    S = np.zeros((MOD, MOD), dtype=np.int64)
    for m in range(MOD):
        S[(m + n) % MOD, m] = 1
    return S


I90 = np.eye(MOD, dtype=np.int64)
A = np.zeros((MOD, MOD), dtype=np.int64)
for n in range(1, 45):
    A += shift(n) - I90

e0 = np.zeros(MOD, dtype=np.int64)
e0[0] = 1


def vev(M):
    return int(e0 @ (M @ e0))


# 注意 A 非自伴；⟨0|A²|0⟩ 與 ⟨0|A⁴|0⟩ 之定義取 A（非 A†）之真空期望
A2 = vev(A @ A)
A4 = vev(A @ A @ A @ A)
# 亦算 A† 版本以備
Ah = A.T
A2h = vev(Ah @ Ah)
A4h = vev(Ah @ Ah @ Ah @ Ah)

print("⟨0|A²|0⟩ =", A2, "｜⟨0|A⁴|0⟩ =", A4)
print("⟨0|A†²|0⟩ =", A2h, "｜⟨0|A†⁴|0⟩ =", A4h)
print("Wick 3·(⟨A²⟩)² =", 3 * A2 ** 2, "⟹ connected = ⟨A⁴⟩ − Wick =", A4 - 3 * A2 ** 2)
print("3N_d = 132｜6·⟨A²⟩ =", 6 * A2)

# 對照：四腿 correlator 之有序總和
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
# 有序 ⟨Φ_{v1}..Φ_{v4}⟩ 總和 = Σ_{ordered v} Σ_S (Πs)·1[Σ s_j v_j ≡ 0]
tot = 0
for s in SIG:
    pr = s[1] * s[2] * s[3]
    # 數 (v1..v4) ∈ [1,44]⁴ 滿足 Σ s_j v_j ≡ 0 mod 90
    cnt = 0
    for v1 in range(1, 45):
        for v2 in range(1, 45):
            for v3 in range(1, 45):
                for v4 in range(1, 45):
                    if (s[0] * v1 + s[1] * v2 + s[2] * v3 + s[3] * v4) % MOD == 0:
                        cnt += 1
    tot += pr * cnt
print("Σ_ordered ⟨Φ^4⟩ =", tot)

out = {"A2": A2, "A4": A4, "A2h": A2h, "A4h": A4h,
       "wick_3A2sq": 3 * A2 ** 2, "connected": A4 - 3 * A2 ** 2,
       "six_times_A2": 6 * A2, "three_Nd": 132,
       "ordered_4pt_sum": tot,
       "note": "A = Σ_{n=1..44}(S^n − 1)，S = 平移；全部整數精確（90×90 int64）"}
with open("/home/wsl/idcm_ckpt/e8/mrB/op_fourwave_H.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("→ 落 mrB/op_fourwave_H.json")
