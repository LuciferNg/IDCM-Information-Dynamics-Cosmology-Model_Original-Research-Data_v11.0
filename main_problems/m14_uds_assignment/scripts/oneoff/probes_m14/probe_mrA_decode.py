#!/usr/bin/python3
"""MR-A 診斷：解碼 op() 之非零目標多重集，核 5→3 與 3→5 之目標是否合乎「去 n 添 m」記賬。"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import e7chain as C

def decode(k, r):
    return tuple(int(x) for x in C.COLEX[k][r])

def show(src, na, ncre, kout):
    v = C.op(np.array([src], dtype=np.int32), np.ones(1), na=na, ncre=ncre, kt=None, out_dim=C.NRK[kout])
    nz = np.nonzero(v)[0]
    print(f"\nsrc={tuple(src)} na={na} ncre={ncre} -> |T|={kout}；非零 {len(nz)} 個")
    for r in nz[:12]:
        print("   ", decode(kout, int(r)), "({:+.15e})".format(v[r]))

show([2, 2, 2, 2, 2], 3, 1, 3)
show([1, 1, 1], 1, 3, 5)
show([1, 1, 1, 1, 1], 3, 1, 3)
show([1, 1, 2], 1, 3, 5)
