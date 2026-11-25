#!/usr/bin/python3
"""MR-A 探針：5→3 與 3→5 之伴隨一致性（引擎原生雙向）。
op() 一次只做「一個源扇區 → 一個目標扇區」之整批映射（out 為 dense 向量、累加）
⟹ 取矩陣元素須逐源呼叫（源 = 1×k，φ=[1]）。
子空間 modes {1,2,3}：3 模 10 個、5 模 21 個 ⟹ 210 元素。
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import e7chain as C

N = C.N
MODES = (1, 2, 3)
T3 = [t for t in C.MULTI[3] if max(t) <= 3]
T5 = [t for t in C.MULTI[5] if max(t) <= 3]
print("N =", N, "｜3 模子集 =", len(T3), "｜5 模子集 =", len(T5), flush=True)

rank = C._rank
def col(src, na, ncre, out_dim):
    v = C.op(np.array([src], dtype=np.int32), np.ones(1), na=na, ncre=ncre, kt=None, out_dim=out_dim)
    return v

# 3→5 : 去 1 添 3
A35 = np.zeros((C.NRK[5], len(T3)))
for j, t in enumerate(T3):
    A35[:, j] = col(t, 1, 3, C.NRK[5])
# 5→3 : 去 3 添 1
A53 = np.zeros((C.NRK[3], len(T5)))
for j, t in enumerate(T5):
    A53[:, j] = col(t, 3, 1, C.NRK[3])

i3 = [int(rank(np.array([t]))[0]) for t in T3]
i5 = [int(rank(np.array([t]))[0]) for t in T5]
sub53 = A53[i3, :]               # 列 = rank(3 模目標)；欄 = T5 之位置（源）
sub35 = A35[i5, :]               # 列 = rank(5 模目標)；欄 = T3 之位置（源）

d = sub53 - sub35.T
print("\n[MR-A1] adjoint diff max = {:.6e}".format(np.max(np.abs(d))), flush=True)
print("        相異元素 = {}/{}".format(int((np.abs(d) > 0).sum()), d.size), flush=True)
print("        5->3 非零 = {}/{:.0f}".format(int((sub53 != 0).sum()), sub53.size), flush=True)
print("        3->5 非零 = {}/{:.0f}".format(int((sub35 != 0).sum()), sub35.size), flush=True)
print("        norm(A53) = {:.15e} | norm(A35^T) = {:.15e}".format(np.linalg.norm(sub53), np.linalg.norm(sub35.T)), flush=True)

print("\n[MR-A2] nonzero (T3 | T5 | A53 | A35T)", flush=True)
for a in range(len(T3)):
    for b in range(len(T5)):
        if sub53[a, b] != 0.0 or sub35[b, a] != 0.0:
            print("  {} | {} | {:+.15e} | {:+.15e}".format(
                tuple(int(x) for x in T3[a]), tuple(int(x) for x in T5[b]), sub53[a, b], sub35[b, a]), flush=True)

np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_sub53.npy", sub53)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA_sub35.npy", sub35)
print("\n已存 probe_mrA_sub53/sub35.npy", flush=True)
