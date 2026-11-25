#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_8 尾鏈（路徑 A 慣例：四波通道 only）—— 乾淨版 V97（Σ p97C_0..7）

φ_5^7 = (V77 φ_4^7 + V97 φ_4^9 + V57 φ_4^5 − E2 φ_3^7)/DEN7
φ_6^5 = (V55 φ_5^5 + V75 φ_5^7 + V35 φ_5^3 − E2 φ_4^5 − E3 φ_3^5 − E4 φ_2^5)/DEN5
φ_7^3 = (V33 φ_6^3 + V53 φ_6^5 − E2 φ_5^3 − E3 φ_4^3 − E4 φ_3^3 − E5 φ_2^3 − E6 φ_1)/DEN3
E_8   = ⟨b, φ_7^3⟩

與舊 e8_tail.py 之差異：
  ① V97 只取 p97C_0..7（乾淨重算，每檔硬核 799,075,328 B）
  ② 上游 E_3/E_5/E_7 直接 np.load（舊版誤用 C.energies）
  ③ TRI／QUI 取 R.TRI／R.QUI（舊版誤用 R.TRI／R.QUI）
  ④ 輸出一律加 _C 後綴，不覆寫舊（污染）結果
"""
import os, sys, json, time
import numpy as np
from math import comb

for _p in ("/home/wsl/idcm_ckpt/e8", "/home/wsl/.hermes/cache/scratch"):
    if os.path.exists(_p + "/e7chain.py"):
        sys.path.insert(0, _p); break
import e7chain as C
import e7run as R

CK = os.environ.get("CK", "/home/wsl/idcm_ckpt/e8")
CK7 = os.environ.get("CK7", "/home/wsl/idcm_ckpt/e7/e7ckpt2")
TAB7 = os.environ.get("TAB7", "/home/wsl/idcm_ckpt/e7/table7.npy")
R.CK = CK7
C.CK = CK7
NQ3, NQ5, NQ7 = comb(46, 3), comb(48, 5), comb(50, 7)
T0 = time.time()
log = lambda m: print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
L = lambda nm: np.load(f"{CK7}/{nm}.npy")
SIZE7 = NQ7 * 8 + 128


def main():
    res = json.load(open(f"{CK7}/chain_result.json"))
    E = {k: res[k] for k in ("E2", "E3", "E4", "E5", "E6", "E7")}
    log("上游 E2..E7 = " + " ".join(f"{E[k]:+.9f}" for k in E))

    parts = [f"{CK}/p97C_{k}.npy" for k in range(8)]
    bad = [p for p in parts if not (os.path.exists(p) and os.path.getsize(p) == SIZE7)]
    if bad:
        raise SystemExit(f"FATAL p97C 分片缺／截斷: {bad}")
    V97 = np.zeros(NQ7)
    for p in parts:
        V97 += np.load(p)
    log(f"V97 分片 8（全部 {SIZE7} B）｜‖V97 φ_4^9‖ = {np.linalg.norm(V97):.9e}")

    E3e = np.load(f"{CK7}/E_3.npy"); E5e = np.load(f"{CK7}/E_5.npy"); E7e = np.load(f"{CK7}/E_7.npy")
    log(f"E_3/E_5/E_7 形狀 {E3e.shape} / {E5e.shape} / {E7e.shape}")
    DEN3 = C.G[1] - E3e; DEN5 = C.G[1] - E5e; DEN7 = C.G[1] - E7e

    phi4_7 = L("phi4_7"); phi4_5 = L("phi4_5"); phi3_7 = L("phi3_7")
    phi5_5 = L("phi5_5"); phi5_3 = L("phi5_3"); phi4_3 = L("phi4_3")
    phi3_5 = L("phi3_5"); phi3_3 = L("phi3_3"); phi2_5 = L("phi2_5")
    phi2_3 = L("phi2_3"); phi6_3 = L("phi6_3"); phi1 = L("phi1"); b = L("b")

    t = time.time()
    V77 = R.par_op(R.D2TAB, phi4_7, 2, 2, 7, NQ7, "e8_t77")
    log(f"V77 φ_4^7 ✓ ‖·‖={np.linalg.norm(V77):.9e}（{time.time()-t:.0f}s，快取）")
    V57 = R.par_op(R.QUI, phi4_5, 1, 3, 7, NQ7, "e8_t57_p4_5")
    phi5_7 = (V77 + V97 + V57 - E["E2"] * phi3_7) / DEN7
    np.save(f"{CK}/phi5_7_C.npy", phi5_7)
    log(f"φ_5^7 ✓ ‖·‖={np.linalg.norm(phi5_7):.9e}")

    V55 = R.par_op(R.QUI, phi5_5, 2, 2, 5, NQ5, "e8_t55_p5_5")
    V75 = R.par_op(R.D2TAB, phi5_7, 3, 1, 5, NQ5, "e8_t75_p5_7")
    V35 = R.par_op(R.TRI, phi5_3, 1, 3, 5, NQ5, "e8_t35_p5_3")
    phi6_5 = (V55 + V75 + V35 - E["E2"] * phi4_5 - E["E3"] * phi3_5 - E["E4"] * phi2_5) / DEN5
    np.save(f"{CK}/phi6_5_C.npy", phi6_5)
    log(f"φ_6^5 ✓ ‖·‖={np.linalg.norm(phi6_5):.9e}")

    V33 = R.par_op(R.TRI, phi6_3, 2, 2, 3, NQ3, "e8_t33_p6_3")
    V53 = R.par_op(R.QUI, phi6_5, 3, 1, 3, NQ3, "e8_t53_p6_5")
    phi7_3 = (V33 + V53 - E["E2"] * phi5_3 - E["E3"] * phi4_3 - E["E4"] * phi3_3
              - E["E5"] * phi2_3 - E["E6"] * phi1) / DEN3
    np.save(f"{CK}/phi7_3_C.npy", phi7_3)
    E8 = float(np.sum(b * phi7_3))
    log(f"φ_7^3 ✓ ‖·‖={np.linalg.norm(phi7_3):.9e}")
    log(f"★ E_8 = {E8:+.15f}")
    json.dump({"chain": "A-clean(p97C)", "E2": E["E2"], "E3": E["E3"], "E4": E["E4"],
               "E5": E["E5"], "E6": E["E6"], "E7": E["E7"], "E8": E8,
               "V97_norm": float(np.linalg.norm(V97)),
               "phi5_7_norm": float(np.linalg.norm(phi5_7)),
               "phi6_5_norm": float(np.linalg.norm(phi6_5)),
               "phi7_3_norm": float(np.linalg.norm(phi7_3)),
               "sec": round(time.time() - T0, 1)},
              open(f"{CK}/e8_result_C.json", "w"), indent=1)
    log("E8C_DONE")


if __name__ == "__main__":
    main()
