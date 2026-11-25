#!/usr/bin/python3
"""MR-A5：raw 公約（∏ c_Q^{c_R} 取代 ∏C(c_Q,c_R)）之伴隨破缺量（判別閘）。
只改移除側多重度因子，其餘（fpA/mulpos/noC/fpC、W 表）與引擎全同。
"""
import os, sys, itertools
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
from math import comb
import numpy as np
import e7chain as C

RAW = os.environ.get("MR_A_RAW", "0") == "1"   # False = ∏C（引擎公約）；True = ∏c_Q^{c_R}

def op_weighted(src, na, ncre, out_dim, raw):
    """單源版 op()，多重度因子可切換 ∏C / raw。"""
    S = np.array(src, dtype=np.int32)
    ksrc = len(S)
    add = C.MULTI[ncre]; addrank = C.RANK[ncre]
    out = np.zeros(out_dim, dtype=np.float64)
    for rm in itertools.combinations(range(ksrc), na):
        keep = [j for j in range(ksrc) if j not in rm]
        A = S[list(rm)]; A_sort = np.sort(A)
        Ssur = S[keep]
        if na == 1:
            fpA = 1.0
        elif na == 2:
            fpA = 2.0 if A_sort[0] == A_sort[1] else 1.0
        else:
            e01 = A_sort[0] == A_sort[1]; e12 = A_sort[1] == A_sort[2]
            fpA = 6.0 if (e01 and e12) else (2.0 if (e01 or e12) else 1.0)
        noA = 1.0; mulpos = 1.0
        for j in range(na):
            v = A_sort[j]
            cT = float((S == v).sum())
            cb = float((A_sort[:j] == v).sum())
            noA *= (cT - cb); mulpos *= (cT - cb) / (cb + 1.0)
        if raw:
            # ∏_m c_Q(m)^{c_R(m)}
            num = 1.0
            for v in set(A_sort.tolist()):
                cQv = float((S == v).sum()); cRv = float((A_sort == v).sum())
                num *= cQv ** cRv
            noA = num * fpA
        if na == 3:
            rkt = int(C._rank(np.array([A_sort]))[0])
        for ai in range(len(add)):
            addv = add[ai]
            if na == 1 and ncre == 3:
                W = int(C.WT13[int(A_sort[0]) - 1, int(addrank[ai])])
            elif na == 3 and ncre == 1:
                W = int(C.WT31[ai, rkt])
            else:
                raise ValueError("僅支援 1↔3")
            if W == 0:
                continue
            noC = 1.0; fpC = 1.0
            as_ = np.sort(addv)
            for j in range(ncre):
                v = as_[j]
                cS = float((Ssur == v).sum())
                cb = float((as_[:j] == v).sum())
                noC *= (cS + cb + 1.0); fpC *= (cb + 1.0)
            T = np.sort(np.concatenate([Ssur, addv]))
            amp = (W / 90.0) * np.sqrt(noA * noC) / (fpA * fpC * mulpos)
            out[int(C._rank(np.array([T]))[0])] += amp
    return out

MODES = (1, 2, 3, 4)
T3 = [t for t in C.MULTI[3] if max(t) <= 4]
T5 = [t for t in C.MULTI[5] if max(t) <= 4]
i3 = [int(C._rank(np.array([t]))[0]) for t in T3]
i5 = [int(C._rank(np.array([t]))[0]) for t in T5]

for raw in (False, True):
    A35 = np.zeros((C.NRK[5], len(T3)))
    for j, t in enumerate(T3):
        A35[:, j] = op_weighted(t, 1, 3, C.NRK[5], raw)
    A53 = np.zeros((C.NRK[3], len(T5)))
    for j, t in enumerate(T5):
        A53[:, j] = op_weighted(t, 3, 1, C.NRK[3], raw)
    s53 = A53[i3, :]; s35 = A35[i5, :]
    d = s53 - s35.T
    tag = "raw ∏c_Q^{c_R}" if raw else "∏C（引擎）"
    print("{:16s} max|A(5→3) − A(3→5)ᵀ| = {:.6e}｜>1e−12 元素 = {}/{}｜max|A| = {:.6e}".format(
        tag, np.max(np.abs(d)), int((np.abs(d) > 1e-12).sum()), d.size, np.max(np.abs(s53))), flush=True)
    if not raw:
        ref53 = s53.copy()
    else:
        print("   兩公約之差 max|A_∏C − A_raw| = {:.6e}".format(np.max(np.abs(ref53 - s53))), flush=True)
