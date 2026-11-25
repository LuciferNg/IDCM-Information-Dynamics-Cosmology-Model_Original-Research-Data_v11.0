#!/usr/bin/python3
"""MR-A 診斷 3：核 rank/COLEX 往返一致性（k=3,5），找出矩陣建構之索引錯。"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
import numpy as np
import e7chain as C

for k in (3, 5):
    r = C._rank(C.MULTI[k])
    inv = C.COLEX[k][np.argsort(r)]
    ok = (inv == C.MULTI[k]).all()
    print(f"k={k}: rank 往返 == MULTI（lex）? {ok}｜rank 範圍 [{r.min()},{r.max()}]｜NRK={C.NRK[k]}")

t5 = (1, 1, 1, 1, 3)
rk5 = int(C._rank(np.array([t5]))[0])
print("\nrank", t5, "=", rk5, "｜COLEX[5][rk] =", tuple(int(x) for x in C.COLEX[5][rk5]))
print("RANK[5] 中該值是否存在:", bool((C.RANK[5] == rk5).any()), "｜其位置 =", int(np.argmax(C.RANK[5] == rk5)))
print("MULTI[5] 該位置 =", tuple(int(x) for x in C.MULTI[5][int(np.argmax(C.RANK[5] == rk5))]))

t3 = (1, 1, 1)
print("\nrank", t3, "=", int(C._rank(np.array([t3]))[0]))
