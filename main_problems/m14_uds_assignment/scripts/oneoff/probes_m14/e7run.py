#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_7 全 CPU 鏈（重建版，路徑 A）
φ_k = R( V φ_{k−1} − Σ_{q=2}^{k−1} E_q φ_{k−1−q} )
所有算子對 scripts/51 vertex_apply 逐點過閘（≤1e−15）
"""
import os, sys, json, time, itertools
from math import comb
import numpy as np
import multiprocessing as mp

sys.path.insert(0, "/home/wsl/.hermes/cache/scratch")
import e7chain as C

N = 44
CK = C.CK
NW = 12
G = C.G
TRI, D2, S1 = C.COLEX[3], C.COLEX[2], C.COLEX[1]
QUI = C.COLEX[5]
NQ3, NQ5, NQ7 = comb(N + 2, 3), comb(N + 4, 5), comb(N + 6, 7)
_T7 = os.environ.get("TABLE7")
if _T7 is None or not os.path.exists(_T7):
    for _c in ("/home/wsl/idcm_ckpt/e7/table7.npy",
               "/home/wsl/.hermes/cache/scratch/e7ckptB/table7.npy",
               "/home/wsl/.hermes/cache/scratch/e7ckpt/table7.npy"):
        if os.path.exists(_c):
            _T7 = _c; break
D2TAB = np.load(_T7, mmap_mode="r")


def _w(w, nw, src, phi, na, ncre, kt, out_dim, tag):
    sp = len(src); step = (sp + nw - 1) // nw
    lo = w * step; hi = min(lo + step, sp)
    if lo >= hi:
        np.save(f"{CK}/tmp_{tag}_{w}.npy", np.zeros(out_dim)); return
    np.save(f"{CK}/tmp_{tag}_{w}.npy", C.op(src[lo:hi], phi[lo:hi], na, ncre, kt, out_dim))


def par_op(src, phi, na, ncre, kt, out_dim, tag, nw=NW):
    final = f"{CK}/op_{tag}.npy"
    if os.path.exists(final):
        return np.load(final)
    ps = [mp.Process(target=_w, args=(w, nw, src, phi, na, ncre, kt, out_dim, tag)) for w in range(nw)]
    [p.start() for p in ps]; [p.join() for p in ps]
    tot = np.zeros(out_dim)
    for w in range(nw):
        f = f"{CK}/tmp_{tag}_{w}.npy"
        tot += np.load(f); os.remove(f)
    np.save(final, tot)
    return tot


def seq_op(src, phi, na, ncre, kt, out_dim, tag):
    final = f"{CK}/op_{tag}.npy"
    if os.path.exists(final):
        return np.load(final)
    r = C.op(src, phi, na, ncre, kt, out_dim)
    np.save(final, r); return r


def energies(NQ, tabfile, tag):
    f = f"{CK}/E_{tag}.npy"
    if os.path.exists(f):
        return np.load(f)
    tab = np.load(f"/home/wsl/.hermes/cache/scratch/e7ckpt/{tabfile}", mmap_mode="r") if tabfile else None
    if tab is not None:
        out = np.empty(len(tab), dtype=np.float64)
        for s in range(0, len(tab), 2_000_000):
            e = min(s + 2_000_000, len(tab))
            out[s:e] = G[np.asarray(tab[s:e])].sum(1)
        del tab
    else:
        out = G[NQ].sum(1) if NQ.ndim == 2 else G[NQ]
    np.save(f, out); return out


def main():
    t0 = time.time()
    log = lambda m: print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)

    def ck(nm, fn):
        f = f"{CK}/{nm}.npy"
        if os.path.exists(f):
            return np.load(f)
        r = fn(); np.save(f, r); log(f"{nm} ✓"); return r

    E3e = energies(TRI, None, "3"); E5e = energies(QUI, None, "5"); E7e = energies(None, "table7.npy", "7")
    DEN3 = G[1] - E3e; DEN5 = G[1] - E5e; DEN7 = G[1] - E7e
    log("能量表 ✓")

    e1 = np.zeros(N); e1[0] = 1.0        # |n=1⟩（S1 有 44 行）
    b = ck("b", lambda: C.op(S1, e1, 1, 3, 3, NQ3))
    phi1 = ck("phi1", lambda: b / DEN3)
    E2 = float(np.sum(b * phi1))
    phi2_3 = ck("phi2_3", lambda: seq_op(TRI, phi1, 2, 2, 3, NQ3, "t33_p1") / DEN3)
    phi2_5 = ck("phi2_5", lambda: seq_op(TRI, phi1, 1, 3, 5, NQ5, "t35_p1") / DEN5)
    E3 = float(np.sum(b * phi2_3))
    log(f"E_2 {E2:+.12f}  E_3 {E3:+.12f}")

    phi3_3 = ck("phi3_3", lambda: (seq_op(TRI, phi2_3, 2, 2, 3, NQ3, "t33_p2_3")
                                   + seq_op(QUI, phi2_5, 3, 1, 3, NQ3, "t53_p2_5") - E2 * phi1) / DEN3)
    E4 = float(np.sum(b * phi3_3))
    phi3_5 = ck("phi3_5", lambda: (seq_op(QUI, phi2_5, 2, 2, 5, NQ5, "t55_p2_5")
                                   + seq_op(TRI, phi2_3, 1, 3, 5, NQ5, "t35_p2_3")) / DEN5)
    phi3_7 = ck("phi3_7", lambda: par_op(QUI, phi2_5, 1, 3, 7, NQ7, "t57_p2_5") / DEN7)
    log(f"E_4 {E4:+.12f}  ‖φ_3^7‖ {np.linalg.norm(phi3_7):.6e}")
    json.dump({"E2": E2, "E3": E3, "E4": E4, "phi3_7_norm": float(np.linalg.norm(phi3_7))},
              open(f"{CK}/partial1.json", "w"), indent=1)

    # φ_4
    phi4_3 = ck("phi4_3", lambda: (seq_op(TRI, phi3_3, 2, 2, 3, NQ3, "t33_p3_3")
                                   + seq_op(QUI, phi3_5, 3, 1, 3, NQ3, "t53_p3_5")
                                   - E2 * phi2_3 - E3 * phi1) / DEN3)
    E5 = float(np.sum(b * phi4_3))
    phi4_5 = ck("phi4_5", lambda: (seq_op(QUI, phi3_5, 2, 2, 5, NQ5, "t55_p3_5")
                                   + par_op(D2TAB, phi3_7, 3, 1, 5, NQ5, "t75_p3_7")
                                   + seq_op(TRI, phi3_3, 1, 3, 5, NQ5, "t35_p3_3")
                                   - E2 * phi2_5) / DEN5)
    n7 = float(np.linalg.norm(phi3_7)); n7_old = 39.832939935398855
    log(f"φ_3^7 範數核對：{n7:.15f} vs 舊 {n7_old:.15f}  相對差 {abs(n7-n7_old)/n7_old:.3e}")
    if abs(n7 - n7_old) / n7_old > 1e-8:   # 舊值來自 GPU float32，容差 1e−8
        raise SystemExit("🔴 φ_3^7 與舊版差 >1e−8 ⟹ v77.npy 不可重用")
    phi4_7 = ck("phi4_7", lambda: (np.load("/home/wsl/.hermes/cache/scratch/e7ckpt/v77.npy") / DEN7))
    log(f"E_5 {E5:+.12f}  ‖φ_4^7‖ {np.linalg.norm(phi4_7):.6e}")
    json.dump({"E2": E2, "E3": E3, "E4": E4, "E5": E5}, open(f"{CK}/partial2.json", "w"), indent=1)

    # φ_5
    phi5_3 = ck("phi5_3", lambda: (seq_op(TRI, phi4_3, 2, 2, 3, NQ3, "t33_p4_3")
                                   + seq_op(QUI, phi4_5, 3, 1, 3, NQ3, "t53_p4_5")
                                   - E2 * phi3_3 - E3 * phi2_3 - E4 * phi1) / DEN3)
    E6 = float(np.sum(b * phi5_3))
    phi5_5 = ck("phi5_5", lambda: (seq_op(QUI, phi4_5, 2, 2, 5, NQ5, "t55_p4_5")
                                   + par_op(D2TAB, phi4_7, 3, 1, 5, NQ5, "t75_p4_7")
                                   + seq_op(TRI, phi4_3, 1, 3, 5, NQ5, "t35_p4_3")
                                   - E2 * phi3_5 - E3 * phi2_5) / DEN5)
    log(f"E_6 {E6:+.12f}")
    json.dump({"E2": E2, "E3": E3, "E4": E4, "E5": E5, "E6": E6},
              open(f"{CK}/partial3.json", "w"), indent=1)

    phi6_3 = ck("phi6_3", lambda: (seq_op(TRI, phi5_3, 2, 2, 3, NQ3, "t33_p5_3")
                                   + seq_op(QUI, phi5_5, 3, 1, 3, NQ3, "t53_p5_5")
                                   - E2 * phi4_3 - E3 * phi3_3 - E4 * phi2_3 - E5 * phi1) / DEN3)
    E7 = float(np.sum(b * phi6_3))
    res = {"n": 1, "path": "A", "E2": E2, "E3": E3, "E4": E4, "E5": E5, "E6": E6, "E7": E7,
           "phi3_7_norm": float(np.linalg.norm(phi3_7)), "phi4_7_norm": float(np.linalg.norm(phi4_7)),
           "phi6_3_norm": float(np.linalg.norm(phi6_3)), "sec": round(time.time() - t0, 1)}
    json.dump(res, open(f"{CK}/chain_result.json", "w"), indent=1)
    log(f"E_7 {E7:+.12f}")
    log("CHAIN_DONE")


if __name__ == "__main__":
    main()