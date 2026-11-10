#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""路徑 B 之 E_6／E_7（全 CPU）：位移分母 N_T/2+5.5 ＋ 同模 Δ=±2 通道
φ_k = R_B( V_B φ_{k−1} − Σ_q E_q φ_{k−1−q} ) ,  R_B = (γ_n − γ − (k/2+5.5))^{-1}
校驗：E_2…E_5 須復現舊 GPU 值（−0.633720033859 / +0.040240117635 / −0.045644034981 / +0.015000708912）
"""
import os, sys, json, time
import numpy as np
import multiprocessing as mp

sys.path.insert(0, "/home/wsl/.hermes/cache/scratch")
import e7chain as C
import e7run as R

CK = C.CK
CKB = "/home/wsl/.hermes/cache/scratch/e7ckptB"
os.makedirs(CKB, exist_ok=True)
NW = int(os.environ.get("NW", "12"))
T7FILE = f"{CKB}/table7.npy"
TRI, QUI, S1 = C.COLEX[3], C.COLEX[5], C.COLEX[1]
NQ3, NQ5, NQ7 = R.NQ3, R.NQ5, R.NQ7
SH = 5.5
G = C.G
T7 = np.load(T7FILE, mmap_mode="r")


def energies7():
    f = f"{CKB}/E7e.npy"
    if os.path.exists(f):
        return np.load(f)
    out = np.empty(len(T7), dtype=np.float64)
    for s in range(0, len(T7), 2_000_000):
        e = min(s + 2_000_000, len(T7))
        out[s:e] = G[np.asarray(T7[s:e])].sum(1)
    np.save(f, out); return out


E3e = G[TRI].sum(1); E5e = G[QUI].sum(1); E7e = energies7()
DEN3 = G[1] - E3e - (1.5 + SH)
DEN5 = G[1] - E5e - (2.5 + SH)
DEN7 = G[1] - E7e - (3.5 + SH)

LOG = lambda m: print(m, flush=True)


def ck(nm, fn):
    f = f"{CKB}/{nm}.npy"
    if os.path.exists(f):
        return np.load(f)
    r = fn(); np.save(f, r); LOG(f"  {nm} ✓"); return r


def V33(p): return C.op(TRI, p, 2, 2, 3, NQ3)
def V35(p): return C.op(TRI, p, 1, 3, 5, NQ5) + C.sm_op(TRI, p, 5, True, NQ5)
def V53(p): return C.op(QUI, p, 3, 1, 3, NQ3) + C.sm_op(QUI, p, 3, False, NQ3)
def V55(p, tag): return R.par_op(QUI, p, 2, 2, 5, NQ5, tag)
def V57(p, tag): return R.par_op(QUI, p, 1, 3, 7, NQ7, tag) + C.sm_op(QUI, p, 7, True, NQ7)
def V75(p7): return R.par_op(T7, p7, 3, 1, 5, NQ5, "B_t75") + C.sm_op(T7, p7, 5, False, NQ5)


def _w77(w, lo, hi):
    led = f"{CKB}/v77_w{w}.json"; part = f"{CKB}/v77_w{w}.npy"
    st = json.load(open(led)) if os.path.exists(led) else {"off": lo}
    if st.get("done"):
        return
    out = np.load(part) if os.path.exists(part) else np.zeros(NQ7)
    phi = np.load(f"{CKB}/phi3_7.npy")
    off = st["off"]; t0 = time.time()
    while off < hi:
        end = min(off + 2_000_000, hi)
        out += C.op(np.array(T7[off:end]), phi[off:end], 2, 2, 7, NQ7)
        off = end
        np.save(part, out); json.dump({"off": off, "done": off >= hi}, open(led, "w"))
        LOG(f"  B77 w{w}: {off}/{hi}  {time.time()-t0:.0f}s")
    st["done"] = True; json.dump(st, open(led, "w"))


def V77B(phi):
    if os.path.exists(f"{CKB}/op_t77.npy"):
        return np.load(f"{CKB}/op_t77.npy")
    np.save(f"{CKB}/phi3_7.npy", phi)
    step = (NQ7 + NW - 1) // NW
    ps = [mp.Process(target=_w77, args=(w, w * step, min((w + 1) * step, NQ7))) for w in range(NW)]
    [p.start() for p in ps]; [p.join() for p in ps]
    tot = np.zeros(NQ7)
    for w in range(NW):
        tot += np.load(f"{CKB}/v77_w{w}.npy")
    np.save(f"{CKB}/op_t77.npy", tot)
    return tot


def main():
    t0 = time.time()
    e1 = np.zeros(44); e1[0] = 1.0
    b = ck("b", lambda: C.op(S1, e1, 1, 3, 3, NQ3) + C.sm_op(S1, e1, 3, True, NQ3))
    p1 = ck("phi1", lambda: b / DEN3)
    E2 = float(np.sum(b * p1))
    p2_3 = ck("phi2_3", lambda: V33(p1) / DEN3)
    p2_5 = ck("phi2_5", lambda: V35(p1) / DEN5)
    E3 = float(np.sum(b * p2_3))
    p3_3 = ck("phi3_3", lambda: (V33(p2_3) + V53(p2_5) - E2 * p1) / DEN3)
    E4 = float(np.sum(b * p3_3))
    p3_5 = ck("phi3_5", lambda: (V55(p2_5, "B_t55_p2_5") + V35(p2_3)) / DEN5)
    p3_7 = ck("phi3_7", lambda: V57(p2_5, "B_t57_p2_5") / DEN7)
    p4_3 = ck("phi4_3", lambda: (V33(p3_3) + V53(p3_5) - E2 * p2_3 - E3 * p1) / DEN3)
    E5 = float(np.sum(b * p4_3))
    LOG(f"[{time.time()-t0:.0f}s] E_2 {E2:+.12f}  E_3 {E3:+.12f}  E_4 {E4:+.12f}  E_5 {E5:+.12f}")
    LOG(f"  舊值:  −0.633720033859 / +0.040240117635 / −0.045644034981 / +0.015000708912")
    json.dump({"E2": E2, "E3": E3, "E4": E4, "E5": E5}, open(f"{CKB}/stage1.json", "w"), indent=1)

    p4_5 = ck("phi4_5", lambda: (V55(p3_5, "B_t55_p3_5") + V75(p3_7) + V35(p3_3) - E2 * p2_5) / DEN5)
    LOG(f"[{time.time()-t0:.0f}s] φ_4^5 ✓")
    p5_3 = ck("phi5_3", lambda: (V33(p4_3) + V53(p4_5) - E2 * p3_3 - E3 * p2_3 - E4 * p1) / DEN3)
    E6 = float(np.sum(b * p5_3))
    LOG(f"[{time.time()-t0:.0f}s] E_6(B，7 模未計) = {E6:+.12f}")

    n77 = V77B(p3_7)
    LOG(f"[{time.time()-t0:.0f}s] V^{{7→7}}φ_3^7 完成 ‖·‖={np.linalg.norm(n77):.6f}")
    p4_7 = ck("phi4_7", lambda: n77 / DEN7)
    p5_5 = ck("phi5_5", lambda: (V55(p4_5, "B_t55_p4_5") + V75(p4_7) + V35(p4_3) - E2 * p3_5 - E3 * p2_5) / DEN5)
    p6_3 = ck("phi6_3", lambda: (V33(p5_3) + V53(p5_5) - E2 * p4_3 - E3 * p3_3 - E4 * p2_3 - E5 * p1) / DEN3)
    E7 = float(np.sum(b * p6_3))
    res = {"n": 1, "path": "B", "E2": E2, "E3": E3, "E4": E4, "E5": E5, "E6": E6, "E7": E7,
           "sec": round(time.time() - t0, 1)}
    json.dump(res, open(f"{CKB}/chain_resultB.json", "w"), indent=1)
    LOG(f"[{time.time()-t0:.0f}s] E_7(B) = {E7:+.12f}")
    LOG("B_DONE")


if __name__ == "__main__":
    main()