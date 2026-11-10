#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPU 全鏈：路徑 A／B，E_1…E_7（7 模態）
E_1 = γ_n（單模態能量，定義值）
輸出一式兩份：data/gpu_chain_<path>.json ＋ log
"""
import os, sys, json, time
import numpy as np
import cupy as cp
from math import comb

sys.path.insert(0, "/home/wsl/.hermes/cache/scratch")
import e7chain as C
import gpu_op as G

CK = C.CK
OUTD = "/home/wsl/.hermes/cache/scratch/gpuout"
os.makedirs(OUTD, exist_ok=True)
NQ3, NQ5, NQ7 = comb(46, 3), comb(48, 5), comb(50, 7)
G_ = C.G
TRI, QUI, S1 = C.COLEX[3], C.COLEX[5], C.COLEX[1]
T7 = np.load("/home/wsl/.hermes/cache/scratch/e7ckptB/table7.npy", mmap_mode="r")
PATH = sys.argv[1] if len(sys.argv) > 1 else "A"
SH = 5.5
LOG = lambda m: print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
T0 = time.time()

if PATH == "A":
    DEN3 = G_[1] - G_[TRI].sum(1)
    DEN5 = G_[1] - G_[QUI].sum(1)
    DEN7 = G_[1] - np.load(f"{CK}/E_7.npy") if os.path.exists(f"{CK}/E_7.npy") else None
else:
    DEN3 = G_[1] - G_[TRI].sum(1) - (1.5 + SH)
    DEN5 = G_[1] - G_[QUI].sum(1) - (2.5 + SH)
    DEN7 = None

if DEN7 is None:
    f7 = f"{OUTD}/E7e.npy"
    if os.path.exists(f7):
        e7e = np.load(f7)
    else:
        e7e = np.empty(len(T7))
        for s in range(0, len(T7), 2_000_000):
            e = min(s + 2_000_000, len(T7))
            e7e[s:e] = G_[np.asarray(T7[s:e])].sum(1)
        np.save(f7, e7e)
    DEN7 = G_[1] - e7e - (3.5 + SH if PATH == "B" else 0.0)


def ck(nm, fn):
    f = f"{OUTD}/{PATH}_{nm}.npy"
    if os.path.exists(f):
        return np.load(f)
    r = fn(); np.save(f, r); LOG(f"{nm} ✓"); return r


def sm(src, phi, kt, create, od):
    return G.gpu_op(src, phi, 0, 2 if create else -2, kt, od, C.WT22)


def main():
    e1 = np.zeros(44); e1[0] = 1.0
    if PATH == "A":
        V33 = lambda p, tag=None: G.gpu_op(TRI, p, 2, 2, 3, NQ3, C.WT22)
        V35 = lambda p, tag=None: G.gpu_op(TRI, p, 1, 3, 5, NQ5, C.WT13)
        V53 = lambda p, tag=None: G.gpu_op(QUI, p, 3, 1, 3, NQ3, C.WT31)
        V55 = lambda p, tag=None: G.gpu_op(QUI, p, 2, 2, 5, NQ5, C.WT22)
        V57 = lambda p, tag=None: G.gpu_op(QUI, p, 1, 3, 7, NQ7, C.WT13)
        V75 = lambda p, tag=None: G.gpu_op(T7, p, 3, 1, 5, NQ5, C.WT31)
        V77 = lambda p, tag=None: G.gpu_op(T7, p, 2, 2, 7, NQ7, C.WT22)
        b = ck("b", lambda: G.gpu_op(S1, e1, 1, 3, 3, NQ3, C.WT13))
    else:
        V33 = lambda p, tag=None: G.gpu_op(TRI, p, 2, 2, 3, NQ3, C.WT22)
        V35 = lambda p, tag=None: G.gpu_op(TRI, p, 1, 3, 5, NQ5, C.WT13) + sm(TRI, p, 5, True, NQ5)
        V53 = lambda p, tag=None: G.gpu_op(QUI, p, 3, 1, 3, NQ3, C.WT31) + sm(QUI, p, 3, False, NQ3)
        V55 = lambda p, tag=None: G.gpu_op(QUI, p, 2, 2, 5, NQ5, C.WT22)
        V57 = lambda p, tag=None: G.gpu_op(QUI, p, 1, 3, 7, NQ7, C.WT13) + sm(QUI, p, 7, True, NQ7)
        V75 = lambda p, tag=None: G.gpu_op(T7, p, 3, 1, 5, NQ5, C.WT31) + sm(T7, p, 5, False, NQ5)
        V77 = lambda p, tag=None: G.gpu_op(T7, p, 2, 2, 7, NQ7, C.WT22)
        b = ck("b", lambda: G.gpu_op(S1, e1, 1, 3, 3, NQ3, C.WT13) + sm(S1, e1, 3, True, NQ3))

    E1 = float(G_[1])
    p1 = ck("phi1", lambda: b / DEN3)
    E2 = float(np.sum(b * p1))
    p2_3 = ck("phi2_3", lambda: V33(p1) / DEN3)
    p2_5 = ck("phi2_5", lambda: V35(p1) / DEN5)
    E3 = float(np.sum(b * p2_3))
    LOG(f"E_2 {E2:+.12f}  E_3 {E3:+.12f}")
    p3_3 = ck("phi3_3", lambda: (V33(p2_3) + V53(p2_5) - E2 * p1) / DEN3)
    E4 = float(np.sum(b * p3_3))
    p3_5 = ck("phi3_5", lambda: (V55(p2_5) + V35(p2_3)) / DEN5)
    p3_7 = ck("phi3_7", lambda: V57(p2_5) / DEN7)
    LOG(f"E_4 {E4:+.12f}  ‖φ_3^7‖ {np.linalg.norm(p3_7):.6e}")
    p4_3 = ck("phi4_3", lambda: (V33(p3_3) + V53(p3_5) - E2 * p2_3 - E3 * p1) / DEN3)
    E5 = float(np.sum(b * p4_3))
    p4_5 = ck("phi4_5", lambda: (V55(p3_5) + V75(p3_7) + V35(p3_3) - E2 * p2_5) / DEN5)
    p4_7 = ck("phi4_7", lambda: V77(p3_7) / DEN7)
    LOG(f"E_5 {E5:+.12f}  ‖φ_4^7‖ {np.linalg.norm(p4_7):.6e}")
    p5_3 = ck("phi5_3", lambda: (V33(p4_3) + V53(p4_5) - E2 * p3_3 - E3 * p2_3 - E4 * p1) / DEN3)
    E6 = float(np.sum(b * p5_3))
    LOG(f"E_6 {E6:+.12f}")
    p5_5 = ck("phi5_5", lambda: (V55(p4_5) + V75(p4_7) + V35(p4_3) - E2 * p3_5 - E3 * p2_5) / DEN5)
    p6_3 = ck("phi6_3", lambda: (V33(p5_3) + V53(p5_5) - E2 * p4_3 - E3 * p3_3 - E4 * p2_3 - E5 * p1) / DEN3)
    E7 = float(np.sum(b * p6_3))
    res = {"path": PATH, "engine": "GPU GTX1060-3GB, float64 accumulate, float32 amp",
           "E1": E1, "E2": E2, "E3": E3, "E4": E4, "E5": E5, "E6": E6, "E7": E7,
           "norms": {"phi3_7": float(np.linalg.norm(p3_7)), "phi4_7": float(np.linalg.norm(p4_7)),
                     "phi6_3": float(np.linalg.norm(p6_3))},
           "sec": round(time.time() - T0, 1)}
    json.dump(res, open(f"{OUTD}/gpu_chain_{PATH}.json", "w"), indent=1)
    LOG(f"E_7 {E7:+.12f}  GPU_{PATH}_DONE")


if __name__ == "__main__":
    main()