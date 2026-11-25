#!/usr/bin/python3
"""MR-A 最終：① 引擎多重度因子之恆等式（noA/fpA ≡ ∏C(c_Q,c_R)、noC/fpC ≡ ∏C(c_S+c_A,c_A)、mulpos）
② 引擎 5↔3 伴隨性（更大窗口）
判準：恆等式逐項精確；伴隨差 ≤ 1e−15 為 float64 噪底。
"""
import os, sys, itertools
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
from math import comb
import numpy as np
import e7chain as C

def counts(t):
    return {m: t.count(m) for m in set(t)}

def binom_pair(cQ, cR, cS, cA):
    """∏C(c_Q,c_R) 與 ∏C(c_S+c_A,c_A)（未出現之模視為 0）"""
    BR = 1
    for m, cr in cR.items():
        BR *= comb(cQ.get(m, 0), cr)
    BA = 1
    for m, ca in cA.items():
        BA *= comb(cS.get(m, 0) + ca, ca)
    return BR, BA

def engine_factors(S, rm, addv):
    """重現 e7chain.op 之多重度因子（單一源、單一 rm、單一 add）"""
    S = list(S); rm = list(rm); addv = list(addv)
    keep = [j for j in range(len(S)) if j not in rm]
    A = [S[j] for j in rm]; A_sort = sorted(A); Ssur = [S[j] for j in keep]
    na, ncre = len(rm), len(addv)
    if na == 1: fpA = 1
    elif na == 2: fpA = 2 if A_sort[0] == A_sort[1] else 1
    else:
        e01 = A_sort[0] == A_sort[1]; e12 = A_sort[1] == A_sort[2]
        fpA = 6 if (e01 and e12) else (2 if (e01 or e12) else 1)
    noA, mulpos = 1, 1
    for j in range(na):
        v = A_sort[j]
        cT = S.count(v)
        cb = sum(1 for jj in range(j) if A_sort[jj] == v)
        noA *= (cT - cb); mulpos *= (cT - cb) / (cb + 1.0)
    noC, fpC = 1, 1
    as_ = sorted(addv)
    for j in range(ncre):
        v = as_[j]
        cS = Ssur.count(v)
        cb = sum(1 for jj in range(j) if as_[jj] == v)
        noC *= (cS + cb + 1); fpC *= (cb + 1)
    return dict(fpA=fpA, noA=noA, mulpos=mulpos, noC=noC, fpC=fpC)

# ---------- ① 恆等式 ----------
print("=" * 70); print("[MR-A3] 多重度因子恆等式（窗口 modes 1..4，全 3/5 模多重集）"); print("=" * 70)
MODE = (1, 2, 3, 4)
bad_R = bad_A = bad_mul = 0; n_R = n_A = 0
maxdev_R = maxdev_A = 0.0
for t in C.MULTI[3]:
    S = [int(x) for x in t]
    if max(S) > 4: continue
    for rm in itertools.combinations(range(3), 3):   # 全去
        pass
maxdev = 0.0
for k, na in ((3, 1), (3, 3), (5, 3), (5, 1)):
    for t in C.MULTI[k]:
        S = [int(x) for x in t]
        if max(S) > 4: continue
        cQ = counts(S)
        for rm in itertools.combinations(range(k), na):
            keep = [j for j in range(k) if j not in rm]
            Ssur = [S[j] for j in keep]; cS = counts(Ssur)
            cR = counts([S[j] for j in rm])
            f = engine_factors(S, rm, [])
            BR, BA = binom_pair(cQ, cR, cS, {})
            # noA/fpA ≡ ∏C(c_Q,c_R)
            dev = abs(f["noA"] / f["fpA"] - BR); maxdev = max(maxdev, dev)
            if dev != 0: bad_R += 1
            n_R += 1
            # mulpos ≡ B_R / ∏ c_R!  ⟹ 檢 ∏C(c_Q,c_R)/mulpos == ∏ c_R!
            fact = 1
            for cr in cR.values():
                fact *= np.math.factorial(cr) if hasattr(np, "math") else 1
            n_A += 0
print("  noA/fpA vs ∏C(c_Q,c_R)：樣本 {}｜不符 {}｜max|dev| = {:.3e}".format(n_R, bad_R, maxdev), flush=True)

# 創建側（另計，避免上面 fact 佔位碼）
import math
maxdevA = 0.0; bad_A = 0; n_A = 0
for k, ncre in ((3, 3), (5, 3), (5, 1)):
    for t in C.MULTI[k]:
        S = [int(x) for x in t]
        if max(S) > 4: continue
        for rm in itertools.combinations(range(k), k - 2):   # 去至剩 2
            keep = [j for j in range(k) if j not in rm]
            Ssur = [S[j] for j in keep]; cS = counts(Ssur)
            for addv in itertools.combinations_with_replacement(MODE, ncre):
                cA = counts(list(addv))
                f = engine_factors(S, rm, list(addv))
                BR, BA = binom_pair(counts(S), counts([S[j] for j in rm]), cS, cA)
                dev = abs(f["noC"] / f["fpC"] - BA); maxdevA = max(maxdevA, dev)
                if dev != 0: bad_A += 1
                n_A += 1
                # 移除側與 mulpos 之關係：∏C(c_Q,c_R) / mulpos == ∏ c_R!
                fr = 1
                for cr in counts([S[j] for j in rm]).values():
                    fr *= math.factorial(cr)
                devm = abs(BR / f["mulpos"] - fr)
                maxdev = max(maxdev, devm)
print("  noC/fpC vs ∏C(c_S+c_A,c_A)：樣本 {}｜不符 {}｜max|dev| = {:.3e}".format(n_A, bad_A, maxdevA), flush=True)
print("  ∏C(c_Q,c_R)/mulpos vs ∏c_R!：max|dev| = {:.3e}".format(maxdev), flush=True)

# ---------- ② 伴隨性（窗口 modes 1..4）----------
print("\n" + "=" * 70); print("[MR-A4] 引擎 5↔3 伴隨性（窗口 modes 1..4）"); print("=" * 70)
T3 = [t for t in C.MULTI[3] if max(t) <= 4]
T5 = [t for t in C.MULTI[5] if max(t) <= 4]
print("  3 模 = {}｜5 模 = {}".format(len(T3), len(T5)), flush=True)
A35 = np.zeros((C.NRK[5], len(T3)))
for j, t in enumerate(T3):
    A35[:, j] = C.op(np.array([t], dtype=np.int32), np.ones(1), na=1, ncre=3, kt=None, out_dim=C.NRK[5])
A53 = np.zeros((C.NRK[3], len(T5)))
for j, t in enumerate(T5):
    A53[:, j] = C.op(np.array([t], dtype=np.int32), np.ones(1), na=3, ncre=1, kt=None, out_dim=C.NRK[3])
i3 = [int(C._rank(np.array([t]))[0]) for t in T3]
i5 = [int(C._rank(np.array([t]))[0]) for t in T5]
s53 = A53[i3, :]; s35 = A35[i5, :]
d = s53 - s35.T
print("  max|A(5→3) − A(3→5)ᵀ| = {:.6e}｜相異元素 = {}/{}".format(np.max(np.abs(d)), int((np.abs(d) > 1e-15).sum()), d.size), flush=True)
print("  非零元素 = {}/{}".format(int((np.abs(s53) > 0).sum()), s53.size), flush=True)
print("  ‖A‖ = {:.15e}".format(np.linalg.norm(s53)), flush=True)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA4_s53.npy", s53)
np.save("/home/wsl/idcm_ckpt/e8/probe_mrA4_s35.npy", s35)
print("\n已存 probe_mrA4_s53/s35.npy", flush=True)
