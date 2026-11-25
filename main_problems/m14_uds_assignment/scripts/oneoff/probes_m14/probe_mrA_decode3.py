#!/usr/bin/python3
"""MR-A 診斷 4：直接檢 stored 矩陣 A53[0, :] 與 T5 對照（找出索引錯）。"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import e7chain as C

T5 = [t for t in C.MULTI[5] if max(t) <= 3]
print("T5 =", [tuple(int(x) for x in t) for t in T5])
A53 = np.zeros((C.NRK[3], len(T5)))
for j, t in enumerate(T5):
    A53[:, j] = C.op(np.array([t], dtype=np.int32), np.ones(1), na=3, ncre=1, kt=None, out_dim=C.NRK[3])
print("\nA53[0, :] =", ["{:+.6e}".format(x) for x in A53[0, :]])
j = [i for i, t in enumerate(T5) if tuple(int(x) for x in t) == (1, 1, 1, 1, 3)][0]
print("index of (1,1,1,1,3) in T5 =", j, "| A53[0, j] =", A53[0, j])
print("rank(1,1,1,1,3) =", int(C._rank(np.array([(1, 1, 1, 1, 3)]))[0]))
