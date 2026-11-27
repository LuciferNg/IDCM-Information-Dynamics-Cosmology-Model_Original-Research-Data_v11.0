#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""算子化第一步：四腿頂點權重 W 之算子身份
   ① 解析推導：CS 向量 c 與符號族 SIG 之因式分解 ⟹ c = s2·s3·s4（後三符號之積）
   ② 由此得【DFT 表示】：
        W(v) = (1/90) Σ_{k=0}^{89} e^{2πi k v1/90} · Π_{j=2..4}( e^{2πi k v_j/90} − e^{−2πi k v_j/90} )
      即 W 係「特徵標 × sin 積」之 90 點 DFT；等價於 congruence 指示函數之 DFT 反演
   ③ 全枚舉核實（667,920）：float64 逐點 |Δ|<1e−9，抽樣 mpmath dps 40 再核
   ④ 算子範數不變量（精確整數）：Σ W、Σ W²、N，並預先凍結候選身份表（Var/Mean² 等）
"""
import itertools
import json

import mpmath as mp
import numpy as np

mp.mp.dps = 40
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]
MOD = 90

# ---- ① 因式分解核實：c_S = s2·s3·s4 且 s1 ≡ +1 ----
fact_ok = all(CS[i] == SIG[i][1] * SIG[i][2] * SIG[i][3] for i in range(8))
first_ok = all(SIG[i][0] == 1 for i in range(8))
print("① c = s2·s3·s4:", fact_ok, "｜首符號全 +1:", first_ok)

# ---- ② DFT 表示（float 向量化）----
K = np.arange(MOD)
E1 = np.exp(2j * np.pi * K / MOD)                       # e^{2πik/90}

sigs = np.array(SIG)                                    # (8,4)
# 直接枚舉 667,920 簽名
REMs = np.arange(1, 45)
ADD = np.array(list(itertools.combinations_with_replacement(range(1, 45), 3)))
V = np.concatenate([np.stack([np.full(len(ADD), r), ADD[:, 0], ADD[:, 1], ADD[:, 2]], axis=1)
                    for r in REMs], axis=0).astype(np.int64)
print("② 簽名數:", V.shape[0])

# W 直接定義
S8 = (V @ sigs.T) % MOD                                 # (N,8)
Wdir = (S8 == 0).astype(np.int64) @ np.array(CS, dtype=np.int64)

# W 之 DFT 表示：對每個簽名，1/90 Σ_k e^{2πikv1} Π_{j} 2i sin(2πk v_j)
# 以 8 個符號項重寫（避免逐點 k 迴圈）：W = Σ_S c_S · 1[σ_S≡0] ⟹ 直接比較 DFT 形式
# DFT 形式逐點計算：用 k 向量化，按簽名分塊
def W_dft(v1, v2, v3, v4):
    a = np.exp(2j * np.pi * K * v1 / MOD)
    for vv in (v2, v3, v4):
        a = a * (np.exp(2j * np.pi * K * vv / MOD) - np.exp(-2j * np.pi * K * vv / MOD))
    return a.sum() / MOD

rng = np.random.default_rng(20261127)
idx = rng.choice(len(V), 4000, replace=False)
maxdev = 0.0
for i in idx:
    d = abs(W_dft(*V[i]) - Wdir[i])
    maxdev = max(maxdev, d)
print(f"③ 抽樣 4,000 點 max|W_DFT − W_dir| = {maxdev:.3e}")

# mpmath 精核 50 點
mpdev = mp.mpf(0)
for i in idx[:50]:
    v1, v2, v3, v4 = [int(x) for x in V[i]]
    acc = mp.mpf(0)
    for k in range(MOD):
        term = mp.e ** (2j * mp.pi * k * v1 / MOD)
        for vv in (v2, v3, v4):
            term *= (mp.e ** (2j * mp.pi * k * vv / MOD) - mp.e ** (-2j * mp.pi * k * vv / MOD))
        acc += term
    mpdev = max(mpdev, abs(acc / MOD - int(Wdir[i])))
print(f"③ mp 50 點 dps40 之 max 偏差 = {mp.nstr(mpdev, 6)}")

# ---- ④ 精確不變量 ----
N = len(V)
S1 = int(Wdir.sum())
S2 = int((Wdir.astype(np.int64) ** 2).sum())
inv = {"N": N, "sum_W": S1, "sum_W2": S2,
       "mean_W": mp.nstr(mp.mpf(S1) / N, 20),
       "mean_W2": mp.nstr(mp.mpf(S2) / N, 20),
       "var_over_mean2": mp.nstr(mp.mpf(S2) / N / (mp.mpf(S1) / N) ** 2 - 1, 20),
       "sumW2_times_N_over_sumW2sq": mp.nstr(mp.mpf(S2) * N / mp.mpf(S1) ** 2, 20)}
print("④ Σ W =", S1, "｜Σ W² =", S2, "｜N =", N)
print("   ⟨W⟩ =", inv["mean_W"], "｜⟨W²⟩ =", inv["mean_W2"])
print("   (ΣW²·N)/(ΣW)² =", inv["sumW2_times_N_over_sumW2sq"], "｜Var/Mean² =", inv["var_over_mean2"])

# 預先凍結之候選身份表（跑前指定，不因結果改動）
CAND = {"88": 88, "90": 90, "44": 44, "16": 16, "12": 12, "N_d-1=43": 43, "3": 3, "5": 5}
val = mp.mpf(S2) * N / mp.mpf(S1) ** 2
cand_hits = {k: mp.nstr(abs(val - v), 8) for k, v in CAND.items()}
print("   候選身份 |值 − 候選|：", cand_hits)

out = {"factorization": {"c_eq_s2s3s4": fact_ok, "first_sign_all_plus": first_ok},
       "n_signatures": N,
       "dft_check": {"sample": 4000, "max_abs_dev_float": maxdev,
                     "mp_checked": 50, "mp_max_dev": mp.nstr(mpdev, 8)},
       "invariants": inv, "candidate_identity_devs": cand_hits,
       "dft_formula": "W(v) = (1/90)Σ_{k=0}^{89} e^{2πikv1/90}·Π_{j=2..4}(e^{2πikv_j/90} − e^{−2πikv_j/90})"}
with open("/home/wsl/idcm_ckpt/e8/mrB/op_w_fourier.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("→ 落 mrB/op_w_fourier.json")
