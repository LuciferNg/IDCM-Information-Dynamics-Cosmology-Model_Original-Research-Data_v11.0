#!/usr/bin/python3
"""MR-B 探針：n=35／39／41 符號例外之逐項分解（先以 n=30／33／43 復現 docs/49 §二 作錨）。
定義（docs/49 §二）：
  e1 = |n⟩；b = V^{1→3} e1；φ1 = b/DEN3, DEN3 = G[n] − E3e
  E_2 = ⟨b|φ1⟩
  X   = V^{3→5} φ1；p₂⁵ = X/DEN5, DEN5 = G[n] − E5e
  量：‖p₂⁵‖、⟨p₂⁵|X⟩、E_2‖p₂⁵‖²、⟨p₂⁵|B p₂⁵⟩（B = V^{5→5}）、Δ₂ = ⟨p₂⁵|B p₂⁵⟩ − E_2‖p₂⁵‖²
用法：python3 probe_mrB_decomp.py <n> [<n> ...]
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, "/home/wsl/idcm_ckpt/e8")
from math import comb
import numpy as np
import e7chain as C

N = 44
G = C.G
S1, TRI, QUI = C.COLEX[1], C.COLEX[3], C.COLEX[5]
NQ3, NQ5 = comb(N + 2, 3), comb(N + 4, 5)
OUT = "/home/wsl/idcm_ckpt/e8/mrB"
os.makedirs(OUT, exist_ok=True)

# 能量表（快取）
def energies(kind, tag, nq):
    f = f"{OUT}/E_{tag}.npy"
    if os.path.exists(f):
        return np.load(f)
    out = G[nq].sum(1)
    np.save(f, out); return out

E3e = energies("3", "3", TRI)
E5e = energies("5", "5", QUI)
print("能量表 ✓  E3e", E3e.shape, "E5e", E5e.shape, flush=True)

def run(n, need_B=True):
    t0 = time.time()
    ln = lambda m: print(f"  [n={n} {time.time()-t0:6.1f}s] {m}", flush=True)
    G1 = float(G[n])
    DEN3 = G1 - E3e
    DEN5 = G1 - E5e
    e1 = np.zeros(N); e1[n - 1] = 1.0
    b = C.op(S1, e1, 1, 3, 3, NQ3)
    phi1 = b / DEN3
    E2 = float(np.sum(b * phi1))
    ln(f"E_2 = {E2:+.12f}｜min|DEN3| = {np.min(np.abs(DEN3)):.6e}")
    X = C.op(TRI, phi1, 1, 3, 5, NQ5)
    ln(f"X 非零 = {int((X != 0).sum())}")
    p25 = X / DEN5
    nrm = float(np.linalg.norm(p25))
    dot = float(np.sum(p25 * X))
    tE2 = E2 * nrm * nrm
    ln(f"‖p₂⁵‖ = {nrm:.6e}｜⟨p₂⁵|X⟩ = {dot:+.6e}｜E_2‖p₂⁵‖² = {tE2:+.6e}｜min|DEN5| = {np.min(np.abs(DEN5)):.6e}")
    res = dict(n=n, E2=E2, den3_min=float(np.min(np.abs(DEN3))), den5_min=float(np.min(np.abs(DEN5))),
               p25_norm=nrm, p25_dot_X=dot, term_E2=tE2)
    if need_B:
        Bp = C.op(QUI, p25, 2, 2, 5, NQ5)
        Bq = float(np.sum(p25 * Bp))
        ln(f"⟨p₂⁵|B p₂⁵⟩ = {Bq:+.6e}｜Δ₂ = {Bq - tE2:+.6e}｜sign(Δ₂)·sign(E_2) = {np.sign(Bq - tE2)*np.sign(E2):+.0f}")
        res.update(p25_Bp=Bq, Delta2=Bq - tE2, prod_sign=float(np.sign(Bq - tE2) * np.sign(E2)),
                   sec=round(time.time() - t0, 1))
        del Bp
    else:
        res["sec"] = round(time.time() - t0, 1)
    return res

if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] or [30, 33, 43]
    import json
    allr = []
    for n in ns:
        print(f"\n=== n = {n} ===", flush=True)
        r = run(n, need_B=os.environ.get("MR_B_NOB", "0") != "1")
        allr.append(r)
        print("  ", json.dumps(r, ensure_ascii=False), flush=True)
    f = f"{OUT}/decomp.json"
    old = json.load(open(f)) if os.path.exists(f) else []
    old = [x for x in old if x["n"] not in ns] + allr
    old.sort(key=lambda x: x["n"])
    json.dump(old, open(f, "w"), indent=1, ensure_ascii=False)
    print("\n已寫", f, "（累計", len(old), "個 n）", flush=True)
