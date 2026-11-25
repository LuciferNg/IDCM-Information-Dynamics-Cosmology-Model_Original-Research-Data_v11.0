#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_8 尾鏈（路徑 B 慣例：四波 ＋ 同模 (0,2)/(2,0)；DEN 位移 k/2+5.5）

φ_5^7 = (V77 φ_4^7 + V97 φ_4^9 + V57 φ_4^5 − E2 φ_3^7)/DEN7_B
φ_6^5 = (V55 φ_5^5 + V75 φ_5^7 + V35 φ_5^3 − E2 φ_4^5 − E3 φ_3^5 − E4 φ_2^5)/DEN5_B
φ_7^3 = (V33 φ_6^3 + V53 φ_6^5 − E2 φ_5^3 − E3 φ_4^3 − E4 φ_3^3 − E5 φ_2^3 − E6 φ_1)/DEN3_B
E_8   = ⟨b, φ_7^3⟩
DEN_k^B = γ₁ − E_k − (k/2 + 5.5)   ⟹ k=3: 7.0｜k=5: 8.0｜k=7: 9.0

算子慣例與 e7runB.py 一致（sm_op 只掛於 k→k±2 通道）：
  V33(2,2) 無 sm｜V55(2,2) 無 sm｜V35(1,3)+sm｜V53(3,1)+sm｜V57(1,3)+sm｜V75(3,1)+sm
輸入：
  V97 = Σ p97D_0..7（乾淨重算，每檔 799,075,328 B）
  V77 = e7ckptB/op_t77.npy（‖·‖=1.370869e−01，B 之 φ_3^7 重算值；非 A 之 397.45）
快取 tag 全部加前綴 Bc_ ⟹ 保證不命中舊污染快取
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
CKB = os.environ.get("CKB", "/home/wsl/idcm_ckpt/e7/e7ckptB")
TAB7 = os.environ.get("TABLE7", "/home/wsl/idcm_ckpt/e7/table7.npy")
R.CK = CKB
C.CK = CKB
NQ3, NQ5, NQ7 = comb(46, 3), comb(48, 5), comb(50, 7)
SIZE7 = NQ7 * 8 + 128
SH = 5.5
T0 = time.time()
log = lambda m: print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)
L = lambda nm: np.load(f"{CKB}/{nm}.npy")


def main():
    res = json.load(open(f"{CKB}/chain_resultB.json"))
    E = {k: res[k] for k in ("E2", "E3", "E4", "E5", "E6")}
    log("上游 E2..E6 = " + " ".join(f"{k}={E[k]:+.12f}" for k in E))
    log(f"記錄 E7(B) = {res.get('E7'):+.18f}（sec={res.get('sec')}）")

    parts = [f"{CK}/p97D_{k}.npy" for k in range(8)]
    bad = [p for p in parts if not (os.path.exists(p) and os.path.getsize(p) == SIZE7)]
    if bad:
        raise SystemExit(f"FATAL p97D 分片缺／截斷: {bad}")
    V97 = np.zeros(NQ7)
    for p in parts:
        V97 += np.load(p)
    log(f"V97(B) 8 片全 {SIZE7} B｜‖V^{{9→7}}φ_4^9(B)‖ = {np.linalg.norm(V97):.6e}")

    E3e = np.load("/home/wsl/idcm_ckpt/e7/e7ckpt2/E_3.npy")
    E5e = np.load("/home/wsl/idcm_ckpt/e7/e7ckpt2/E_5.npy")
    E7e = np.load("/home/wsl/idcm_ckpt/e7/e7ckpt2/E_7.npy")
    log(f"E_3/E_5/E_7 形狀 {E3e.shape} / {E5e.shape} / {E7e.shape}")
    DEN3 = C.G[1] - E3e - (1.5 + SH)
    DEN5 = C.G[1] - E5e - (2.5 + SH)
    DEN7 = C.G[1] - E7e - (3.5 + SH)
    log(f"DEN9 位移 10.0｜DEN7 {3.5+SH}｜DEN5 {2.5+SH}｜DEN3 {1.5+SH}")

    phi4_7 = L("phi4_7"); phi4_5 = L("phi4_5"); phi3_7 = L("phi3_7")
    phi5_5 = L("phi5_5"); phi5_3 = L("phi5_3"); phi4_3 = L("phi4_3")
    phi3_5 = L("phi3_5"); phi3_3 = L("phi3_3"); phi2_5 = L("phi2_5")
    phi2_3 = L("phi2_3"); phi6_3 = L("phi6_3"); phi1 = L("phi1"); b = L("b")

    V77 = np.load(f"{CKB}/op_t77.npy")
    if abs(np.linalg.norm(V77) - 3.974540e+02) < 1e-3:
        raise SystemExit("FATAL op_t77 係 A 之結果（‖397.45‖）⟹ 拒絕使用")
    log(f"V77(B) ✓ ‖·‖={np.linalg.norm(V77):.6e}（e7ckptB 快取）")
    V57 = R.par_op(R.QUI, phi4_5, 1, 3, 7, NQ7, "Bc_t57_p4_5") + C.sm_op(R.QUI, phi4_5, 7, True, NQ7)
    phi5_7 = (V77 + V97 + V57 - E["E2"] * phi3_7) / DEN7
    np.save(f"{CK}/phi5_7_B.npy", phi5_7)
    log(f"φ_5^7(B) ✓ ‖·‖={np.linalg.norm(phi5_7):.6e}")

    V55 = R.par_op(R.QUI, phi5_5, 2, 2, 5, NQ5, "Bc_t55_p5_5")
    V75 = R.par_op(R.D2TAB, phi5_7, 3, 1, 5, NQ5, "Bc_t75_p5_7") + C.sm_op(R.D2TAB, phi5_7, 5, False, NQ5)
    V35 = C.op(R.TRI, phi5_3, 1, 3, 5, NQ5) + C.sm_op(R.TRI, phi5_3, 5, True, NQ5)
    phi6_5 = (V55 + V75 + V35 - E["E2"] * phi4_5 - E["E3"] * phi3_5 - E["E4"] * phi2_5) / DEN5
    np.save(f"{CK}/phi6_5_B.npy", phi6_5)
    log(f"φ_6^5(B) ✓ ‖·‖={np.linalg.norm(phi6_5):.6e}")

    V33 = C.op(R.TRI, phi6_3, 2, 2, 3, NQ3)
    V53 = C.op(R.QUI, phi6_5, 3, 1, 3, NQ3) + C.sm_op(R.QUI, phi6_5, 3, False, NQ3)
    phi7_3 = (V33 + V53 - E["E2"] * phi5_3 - E["E3"] * phi4_3 - E["E4"] * phi3_3
              - E["E5"] * phi2_3 - E["E6"] * phi1) / DEN3
    np.save(f"{CK}/phi7_3_B.npy", phi7_3)
    E8 = float(np.sum(b * phi7_3))
    log(f"φ_7^3(B) ✓ ‖·‖={np.linalg.norm(phi7_3):.6e}")
    log(f"★ E_8(B) = {E8:+.15f}")
    json.dump({"chain": "B-clean(p97D)", "E2": E["E2"], "E3": E["E3"], "E4": E["E4"],
               "E5": E["E5"], "E6": E["E6"], "E7": res.get("E7"), "E8": E8,
               "V97_norm": float(np.linalg.norm(V97)),
               "phi5_7_norm": float(np.linalg.norm(phi5_7)),
               "phi6_5_norm": float(np.linalg.norm(phi6_5)),
               "phi7_3_norm": float(np.linalg.norm(phi7_3)),
               "DEN_shifts": {"9": 9/2+SH, "7": 3.5+SH, "5": 2.5+SH, "3": 1.5+SH},
               "sec": round(time.time() - T0, 1)},
              open(f"{CK}/e8_result_B.json", "w"), indent=1)
    log("E8B_DONE")


if __name__ == "__main__":
    main()
