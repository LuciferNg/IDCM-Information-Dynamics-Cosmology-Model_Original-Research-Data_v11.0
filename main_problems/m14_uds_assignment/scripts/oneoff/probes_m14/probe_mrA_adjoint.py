#!/usr/bin/python3
"""MR-A 探針：5→3 與 3→5 之伴隨一致性（引擎原生雙向），並用重複模案例判別
∏C(c_Q,c_S)  vs  ∏ c_Q^{c_S} 兩個多元公約。
只用 modes {1,2,3} 子空間（3 模 = C(5,3)=10；5 模 = C(7,5)=21）⟹ 210 元素，精確可枚舉。
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import itertools
import e7chain as C

N = C.N
print("N =", N, "｜NRK3 =", C.NRK[3], "｜NRK5 =", C.NRK.get(5), flush=True)

MODES = (1, 2, 3)
T3 = np.array([t for t in C.MULTI[3] if max(t) <= 3], dtype=np.int32)
T5 = np.array([t for t in C.MULTI[5] if max(t) <= 3], dtype=np.int32)
print("3 模子集 =", len(T3), "｜5 模子集 =", len(T5), flush=True)

# ---- 引擎：3→5（去 1 添 3）----
A35 = C.op(T3.astype(np.int32), np.ones(len(T3)), na=1, ncre=3, kt=None, out_dim=C.NRK[5])
# ---- 引擎：5→3（去 3 添 1）----
A53 = C.op(T5.astype(np.int32), np.ones(len(T5)), na=3, ncre=1, kt=None, out_dim=C.NRK[3])

# 索引對照
r3 = {tuple(t): C._rank(np.array([t]))[0] for t in T3}
r5 = {tuple(t): C._rank(np.array([t]))[0] for t in T5}
i3 = [r3[tuple(t)] for t in T3]
i5 = [r5[tuple(t)] for t in T5]
sub53 = A53[i3][:, i5]            # 5→3 限制（列=3模目標，欄=5模源）
sub35 = A35[i5][:, i3]            # 3→5 限制（列=5模目標，欄=3模源）
d = sub53 - sub35.T
print("\n[MR-A1] 引擎雙向伴隨差 max|A53 − A35ᵀ| =", np.max(np.abs(d)), flush=True)
print("          非零元素 =", int((d != 0).sum()), "/", d.size, flush=True)
print("          ‖A53‖ =", np.linalg.norm(sub53), "‖A35ᵀ‖ =", np.linalg.norm(sub35.T), flush=True)
nonzero = int((sub53 != 0).sum())
print("          5→3 非零 =", nonzero, "/", sub53.size, flush=True)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_sub53.npy", sub53)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_sub35.npy", sub35)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_T3.npy", T3)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_T5.npy", T5)
print("\n已存 probe_mrA_*.npy", flush=True)
