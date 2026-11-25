#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_8 尾鏈（路徑 A 慣例：四波通道 only，與 e7run.py 一致）

φ_5^7 = (V77 φ_4^7 + V97 φ_4^9 + V57 φ_4^5 − E2 φ_3^7)/DEN7
φ_6^5 = (V55 φ_5^5 + V75 φ_5^7 + V35 φ_5^3 − E2 φ_4^5 − E3 φ_3^5 − E4 φ_2^5)/DEN5
φ_7^3 = (V33 φ_6^3 + V53 φ_6^5 − E2 φ_5^3 − E3 φ_4^3 − E4 φ_3^3 − E5 φ_2^3 − E6 φ_1)/DEN3
E_8   = ⟨b, φ_7^3⟩
"""
import os, sys, json, time
import numpy as np
from math import comb

sys.path.insert(0, "/home/wsl/.hermes/cache/scratch")
import e7chain as C
import e7run as R

CK = "/home/wsl/.hermes/cache/scratch/e8ckpt"
CK7 = "/home/wsl/.hermes/cache/scratch/e7ckpt2"
NQ3, NQ5, NQ7 = comb(46, 3), comb(48, 5), comb(50, 7)
T0 = time.time()
log = lambda m: print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
L = lambda nm: np.load(f"{CK7}/{nm}.npy")


def main():
    res = json.load(open(f"{CK7}/chain_result.json"))
    E = {k: res[k] for k in ("E2", "E3", "E4", "E5", "E6", "E7")}
    log("上游 E2..E7 = " + " ".join(f"{E[k]:+.9f}" for k in E))

    parts = [f"{CK}/p97_w{k}.npy" for k in range(12)] + [f"{CK}/p97x_{k}.npy" for k in range(12)] + [f"{CK}/p97z_{k}.npy" for k in range(12)]
    have = [p for p in parts if os.path.exists(p)]
    log(f"V97 分片 {len(have)}（主 8 ＋ 續 12）")
    V97 = np.zeros(NQ7)
    for p in have:
        V97 += np.load(p)
    log(f"‖V97 φ_4^9‖ = {np.linalg.norm(V97):.9e}")

    E3e = C.energies(C.TRI, None, "3"); E5e = C.energies(C.QUI, None, "5")
    E7e = C.energies(None, "table7.npy", "7")
    DEN3 = C.G[1] - E3e; DEN5 = C.G[1] - E5e; DEN7 = C.G[1] - E7e

    phi4_7 = L("phi4_7"); phi4_5 = L("phi4_5"); phi3_7 = L("phi3_7")
    phi5_5 = L("phi5_5"); phi5_3 = L("phi5_3"); phi4_3 = L("phi4_3")
    phi3_5 = L("phi3_5"); phi3_3 = L("phi3_3"); phi2_5 = L("phi2_5")
    phi2_3 = L("phi2_3"); phi6_3 = L("phi6_3"); phi1 = L("phi1"); b = L("b")

    t = time.time()
    V77 = R.par_op(R.D2TAB, phi4_7, 2, 2, 7, NQ7, "e8_t77")
    log(f"V77 φ_4^7 ✓ ‖·‖={np.linalg.norm(V77):.9e}（{time.time()-t:.0f}s）")
    V57 = R.par_op(C.QUI, phi4_5, 1, 3, 7, NQ7, "e8_t57_p4_5")
    phi5_7 = (V77 + V97 + V57 - E["E2"] * phi3_7) / DEN7
    np.save(f"{CK}/phi5_7.npy", phi5_7)
    log(f"φ_5^7 ✓ ‖·‖={np.linalg.norm(phi5_7):.9e}")

    V55 = R.par_op(C.QUI, phi5_5, 2, 2, 5, NQ5, "e8_t55_p5_5")
    V75 = R.par_op(R.D2TAB, phi5_7, 3, 1, 5, NQ5, "e8_t75_p5_7")
    V35 = R.par_op(C.TRI, phi5_3, 1, 3, 5, NQ5, "e8_t35_p5_3")
    phi6_5 = (V55 + V75 + V35 - E["E2"] * phi4_5 - E["E3"] * phi3_5 - E["E4"] * phi2_5) / DEN5
    np.save(f"{CK}/phi6_5.npy", phi6_5)
    log(f"φ_6^5 ✓ ‖·‖={np.linalg.norm(phi6_5):.9e}")

    V33 = R.par_op(C.TRI, phi6_3, 2, 2, 3, NQ3, "e8_t33_p6_3")
    V53 = R.par_op(C.QUI, phi6_5, 3, 1, 3, NQ3, "e8_t53_p6_5")
    phi7_3 = (V33 + V53 - E["E2"] * phi5_3 - E["E3"] * phi4_3 - E["E4"] * phi3_3
              - E["E5"] * phi2_3 - E["E6"] * phi1) / DEN3
    np.save(f"{CK}/phi7_3.npy", phi7_3)
    E8 = float(np.sum(b * phi7_3))
    log(f"φ_7^3 ✓ ‖·‖={np.linalg.norm(phi7_3):.9e}")
    log(f"★ E_8 = {E8:+.12f}")
    json.dump({"E2": E["E2"], "E3": E["E3"], "E4": E["E4"], "E5": E["E5"],
               "E6": E["E6"], "E7": E["E7"], "E8": E8,
               "phi5_7_norm": float(np.linalg.norm(phi5_7)),
               "phi6_5_norm": float(np.linalg.norm(phi6_5)),
               "phi7_3_norm": float(np.linalg.norm(phi7_3)),
               "sec": round(time.time() - T0, 1)},
              open(f"{CK}/e8_result.json", "w"), indent=1)
    log("E8_DONE")


if __name__ == "__main__":
    main()