#!/usr/bin/python3
"""MR-A 診斷 2：單源 (1,1,1,1,3) 之 5→3 全非零表；並核 WT31/WT13 同秩值與 chi4。"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import e7chain as C

def decode(k, r):
    return tuple(int(x) for x in C.COLEX[k][r])

tri = (1, 1, 3)
rk = int(C._rank(np.array([tri]))[0])
print("rank(tri) =", rk, "｜COLEX[3][rk] =", decode(3, rk), "｜MULTI[3][rk] =", tuple(int(x) for x in C.MULTI[3][rk]))
print("chi4([[1,1,1,3]]) =", C.chi4(np.array([[1, 1, 1, 3]])))
print("WT31[0, rk] =", int(C.WT31[0, rk]), "｜WT13[0, rk] =", int(C.WT13[0, rk]))

src = [1, 1, 1, 1, 3]
v = C.op(np.array([src], dtype=np.int32), np.ones(1), na=3, ncre=1, kt=None, out_dim=C.NRK[3])
nz = np.nonzero(v)[0]
print("\nsrc =", tuple(src), "（5→3）非零", len(nz), "個")
for r in nz:
    print("   ->", decode(3, int(r)), " = {:+.15e}".format(v[r]))
