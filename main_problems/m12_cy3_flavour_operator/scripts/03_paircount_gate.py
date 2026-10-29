#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / G6–G7 —— D9 ansatz 之配對計數自洽性檢查（反向法之真正入口）
2026-10-28

G6：以 D6 之配對計數（重子 3 價 → 所觸對 2；介子 2 價 → 所觸對 1）重建最小系統：
      N(udd) = a_uu + a_ud        （u↔d 同構：a_dd = a_uu）
      Σ(uds) = a_ud + a_us
      Ξ(uss) = a_us + a_ss        ← D9 寫 2·a_us，等價於額外假設 a_ss = a_us
      K(u s̄) = a_usbar
    未知 5（a_uu, a_ud, a_us, a_ss, a_usbar）、錨 4 ⟹ **1 自由參數**（D9 之 dof=0 係
    由未推導之 a_ss = a_us 造成）。

G7：model-free 之不變量（唔含任何 ansatz）與其 σ。

跑法：/usr/bin/python3 scripts/03_paircount_gate.py
"""
from __future__ import annotations
import json
import os
from mpmath import mp, mpf, sqrt

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


# 錨：⟨O⟩ 與 σ（同 scripts/01／02）
O = {"N": mpf("1.1263"), "S": mpf("1.9090"), "X": mpf("2.1719"), "K": mpf("2.4213")}
SD = {"N": mpf("0.1848"), "S": mpf("0.1575"), "X": mpf("0.1038"), "K": mpf("0.1007")}
D = {k: O[k] - 1 for k in O}   # 增量


def g6():
    say("=" * 80)
    say("G6  配對計數系統（D6 一致）：5 未知 / 4 錨 ⟹ 1 自由參數")
    say("")
    say("  解（以 a_ss 為自由參數 t）:")
    say("    a_usbar = ⟨O⟩_K − 1                    = %+0.6f" % D["K"])
    say("    a_us    = (⟨O⟩_Ξ − 1) − t              = %+0.6f − t" % D["X"])
    say("    a_ud    = (⟨O⟩_Σ − 1) − a_us           = %+0.6f + t" % (D["S"] - D["X"]))
    say("    a_uu    = (⟨O⟩_N − 1) − a_ud           = %+0.6f − t" % (D["N"] - (D["S"] - D["X"])))
    say("")
    # D9 等價條件 a_ss = a_us
    t9 = D["X"] / 2
    say("  ★ D9 之隱含假設：a_ss = a_us ⟹ t = (⟨O⟩_Ξ−1)/2 = %0.6f" % t9)
    t = t9
    vals = {"usbar": D["K"], "us": D["X"] - t, "ud": (D["S"] - D["X"]) + t,
            "uu": (D["N"] - (D["S"] - D["X"])) - t}
    say("      ⟹ a_uu = %+0.6f  a_ud = %+0.6f  a_us = %+0.6f  a_ss = %+0.6f"
        % (vals["uu"], vals["ud"], vals["us"], t))
    say("      ⟹ 與 D9 之 (−0.19675, +0.32305, +0.58595) 逐位相同 ✅")
    say("      ⟹ 即 D9 之 a_uu<0 **完全由** a_ss = a_us 造成，非由數據決定。")
    say("")
    say("  若 t 自由：a_uu = %+0.6f − t ⟹ a_uu < 0 ⟺ t > %0.6f"
        % (D["N"] - (D["S"] - D["X"]), D["N"] - (D["S"] - D["X"])))
    say("  ⟹ 「a_uu < 0」**不是約束**，係 t 之選取（後驗）。")
    say("")
    say("  反過來，一條**真結構**輸入可切成 t：同味對短程正核（g(d)>0 於 d < 18.5933079682）")
    say("  ⟹ a_uu ≥ 0 ⟹ t ≤ %0.6f ⟹ a_us ≥ %0.6f（vs D9 之 %0.6f）"
        % (D["N"] - (D["S"] - D["X"]), D["X"] - (D["N"] - (D["S"] - D["X"])), D["X"] - t9))
    OUT["G6"] = {
        "family": {"usbar": str(D["K"]), "us": str(D["X"]) + " - t",
                   "ud": str(D["S"] - D["X"]) + " + t",
                   "uu": str(D["N"] - (D["S"] - D["X"])) + " - t"},
        "d9_t": str(t9), "d9_reproduces": True,
        "d9_equivalent_to": "a_ss = a_us (unforced)",
        "positivity_bound_t_max": str(D["N"] - (D["S"] - D["X"])),
        "positivity_bound_a_us_min": str(D["X"] - (D["N"] - (D["S"] - D["X"]))),
    }


def g7():
    say("")
    say("=" * 80)
    say("G7  model-free 不變量（唯一可入反向法之量）")
    say("")
    rows = [("n−p", D["N"], SD["N"]), ("Σ⁻−Σ⁰", D["S"], SD["S"]),
            ("Ξ⁻−Ξ⁰", D["X"], SD["X"]), ("K⁺−K⁰", D["K"], SD["K"])]
    for n, v, s in rows:
        say("    %-8s ⟨O⟩−1 = %+0.4f ± %0.4f   z(=0) = %+0.2f" % (n, v, s, v / s))
    say("")
    say("  (i)  全正：4/4（最弱 +0.68σ，最強 +14.11σ）")
    d1, s1 = D["S"] - D["N"], sqrt(SD["S"] ** 2 + SD["N"] ** 2)
    d2, s2 = D["X"] - D["S"], sqrt(SD["X"] ** 2 + SD["S"] ** 2)
    say("  (ii) s 數單調：Δ1 = Σ−N = %+0.4f ± %0.4f (%.2fσ)；Δ2 = Ξ−Σ = %+0.4f ± %0.4f (%.2fσ)"
        % (d1, s1, d1 / s1, d2, s2, d2 / s2))
    say("       ⟹ 單調成立但 Δ2 僅 %.2fσ ⟹ 非獨立穩健。" % (d2 / s2))
    dK, sK = D["K"] - D["X"], sqrt(SD["K"] ** 2 + SD["X"] ** 2)
    say("  (iii) s̄ 通道 > 2 個 s：(⟨O⟩_K−1) − (⟨O⟩_Ξ−1) = %+0.4f ± %0.4f (%.2fσ)"
        % (dK, sK, dK / sK))
    say("")
    say("  ⟹ 可入反向法之獨立約束 ≈ 3 個（正性方向、單調方向、s̄>2s 方向），")
    say("     全部只係**方向性**（0.7σ–14σ 不等），冇一個釘得住幅度。")
    OUT["G7"] = {"delta1": [str(d1), str(s1), str(d1 / s1)],
                 "delta2": [str(d2), str(s2), str(d2 / s2)],
                 "dK_minus_dX": [str(dK), str(sK), str(dK / sK)],
                 "independent_directional_constraints": 3}


def main():
    g6()
    g7()
    with open(os.path.join(ROOT, "data", "03_paircount_gate.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    say("")
    say("[written] data/03_paircount_gate.json , data/run_03_paircount_gate.log")

    with open(os.path.join(ROOT, "data", "run_03_paircount_gate.log"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()