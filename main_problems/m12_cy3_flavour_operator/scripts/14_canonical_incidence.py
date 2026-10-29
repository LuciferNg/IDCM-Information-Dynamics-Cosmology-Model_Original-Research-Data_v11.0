#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R3-P3 —— 唯一典範入射規則（O137 切片）之預註冊檢定

入射規則（來源 O137 `docs/30 §4`，非本輪自創）：
  9 物質模態 = {1,2,3,4,5,7,9,11,13}（n 指標）
  世代 = 按 n 排序每 N_gen=3 個切一塊：g1{1,2,3} g2{4,5,7} g3{9,11,13}
  味   = 區塊內序位（slot）：S1{1,4,9} S2{2,5,11} S3{3,7,13}
  ⟹ 3×3 格，**零自由度**（無需人手指定）。

事前判準（開算前寫死；不事後改）：
  K1  味荷 = slot 欄之 Σn：W = (14, 18, 23) → 歸一 → 與 (1, r_Σ, r_Ξ) 比
  K2  味荷 = slot 欄之 Σγ：γ_n = 4sin²(nπ/88) → 歸一 → 同上比
  K3【範疇性】γ > 0 於 9 模式 ⟹ 任何欄和／欄差皆 > 0 ⟹ **K 通道之負號不可由
      「加性 slot 荷」產生**（categorical kill of the K channel）
  K4【範疇性】slot 欄和是否嚴格單調（若否，則連「序」都無）
  只看 slot 序（1,2,3）↔ (N, Σ, Ξ)，**不作任何置換掃描**（置換 = 自由度）。
  判死門檻：|z| ≥ 3（σ 見 docs/08：r_Σ 0.029705、r_Ξ 0.163837）。

跑法：/usr/bin/python3 scripts/14_canonical_incidence.py
"""
from __future__ import annotations
import json
import os

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


MODES = [1, 2, 3, 4, 5, 7, 9, 11, 13]
NG = 3
BLOCKS = [MODES[i:i + NG] for i in range(0, len(MODES), NG)]
SLOTS = [[BLOCKS[g][f] for g in range(NG)] for f in range(NG)]

# 目標（docs/08 §二；純質量分割比，σ 由 PDG 質量精度支配）
R_SIG = mp.mpf("3.716756")
R_XI = mp.mpf("5.296397")
S_SIG = mp.mpf("0.029705")
S_XI = mp.mpf("0.163837")


def main():
    say("=" * 84)
    say("R3-P3  典範入射規則（O137 切片，零自由度）")
    say("")
    say("  9 模式 = %s ｜Σn = %d（= F₁₀，與 O137 一致）" % (MODES, sum(MODES)))
    say("  世代分塊：%s" % BLOCKS)
    say("  味 slot 欄：%s" % SLOTS)
    say("  欄之和（檢查）：%s" % [sum(s) for s in SLOTS])
    say("")
    say("  註：O137 之『Σ = 55 = F₁₀』係 Σn，非 Σγ。Σγ(9 模式) = %s"
        % mp.nstr(sum(4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2 for n in MODES), 12))
    OUT["setup"] = {"modes": MODES, "blocks": BLOCKS, "slots": SLOTS,
                    "sum_n": sum(MODES)}

    say("")
    say("=" * 84)
    say("K1  味荷 = slot 欄 Σn")
    Wn = [sum(s) for s in SLOTS]
    wn = [mp.mpf(w) / Wn[0] for w in Wn]
    say("  W = %s ⟹ 歸一 = %s" % (Wn, [mp.nstr(x, 8) for x in wn]))
    say("  目標（slot1→N, slot2→Σ, slot3→Ξ）：1, %.6f, %.6f" % (float(R_SIG), float(R_XI)))
    z2 = (wn[1] - R_SIG) / S_SIG
    z3 = (wn[2] - R_XI) / S_XI
    say("  slot2：%s vs %.6f ⟹ z = %+.1f" % (mp.nstr(wn[1], 8), float(R_SIG), float(z2)))
    say("  slot3：%s vs %.6f ⟹ z = %+.1f" % (mp.nstr(wn[2], 8), float(R_XI), float(z3)))
    say("  ⟹ %s" % ("★ 未判死" if abs(z2) < 3 and abs(z3) < 3 else "🔴 判死（|z| ≥ 3）"))
    OUT["K1"] = {"W_n": Wn, "norm": [mp.nstr(x, 10) for x in wn],
                 "z_sigma": float(z2), "z_xi": float(z3)}

    say("")
    say("=" * 84)
    say("K2  味荷 = slot 欄 Σγ（γ_n = 4sin²(nπ/88)，精確於 ℚ(ζ₈₈)）")
    g = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
    Wg = [sum(g(n) for n in s) for s in SLOTS]
    wg = [w / Wg[0] for w in Wg]
    say("  γ 值：%s" % {n: mp.nstr(g(n), 8) for n in MODES})
    say("  W = %s ⟹ 歸一 = %s" % ([mp.nstr(x, 8) for x in Wg], [mp.nstr(x, 8) for x in wg]))
    z2 = (wg[1] - R_SIG) / S_SIG
    z3 = (wg[2] - R_XI) / S_XI
    say("  slot2：%s vs %.6f ⟹ z = %+.1f" % (mp.nstr(wg[1], 8), float(R_SIG), float(z2)))
    say("  slot3：%s vs %.6f ⟹ z = %+.1f" % (mp.nstr(wg[2], 8), float(R_XI), float(z3)))
    say("  ⟹ %s" % ("★ 未判死" if abs(z2) < 3 and abs(z3) < 3 else "🔴 判死（|z| ≥ 3）"))
    OUT["K2"] = {"W_gamma": [mp.nstr(x, 10) for x in Wg],
                 "norm": [mp.nstr(x, 10) for x in wg],
                 "z_sigma": float(z2), "z_xi": float(z3)}

    say("")
    say("=" * 84)
    say("K3  範疇性：K 通道之負號可否由加性 slot 荷產生？")
    say("")
    say("  9 模式下 γ_n 全部 > 0（min γ_1 = %s）；Σn 全部 > 0" % mp.nstr(g(1), 8))
    say("  ⟹ 任何『欄和』或『欄差（後−前）』皆 > 0 ⟹ 無負值可用")
    say("  目標 r_K = −3.041756 ± 0.015940（**負**）")
    say("  ⟹ 🔴 **加性 slot 荷範疇上不可能產生 K 通道之負號**（categorical kill）")
    OUT["K3"] = {"all_positive": True, "K_sign_impossible": True}

    say("")
    say("=" * 84)
    say("K4  範疇性：slot 欄和之單調性")
    for name, W in (("Σn", Wn), ("Σγ", Wg)):
        mono = all(W[i] < W[i + 1] for i in range(len(W) - 1))
        say("  %s：%s ⟹ %s" % (name, [mp.nstr(mp.mpf(x), 8) for x in W],
                              "嚴格遞增" if mono else "非嚴格遞增"))
    OUT["K4"] = {"Wn": Wn, "Wg": [mp.nstr(x, 10) for x in Wg]}

    say("")
    say("=" * 84)
    say("R3-P3 判詞")
    say("")
    say("  K1（Σn 欄和）：z = %+.1f / %+.1f ⟹ 🔴 判死" % (float(OUT["K1"]["z_sigma"]),
                                                     float(OUT["K1"]["z_xi"])))
    say("  K2（Σγ 欄和）：z = %+.1f / %+.1f ⟹ %s"
        % (float(OUT["K2"]["z_sigma"]), float(OUT["K2"]["z_xi"]),
           "🔴 判死" if max(abs(OUT["K2"]["z_sigma"]), abs(OUT["K2"]["z_xi"])) >= 3 else "★ 未判死"))
    say("  K3：加性 slot 荷不可能出負號 ⟹ 🔴 K 通道範疇上不可表示")
    say("  ⟹ 唯一典範入射規則（O137 切片）已檢；**加性（對角）路線於此入射下亦死**。")
    say("  ⟹ 剩下唯一未試者：**非對角**（slot 間相互作用）——但唯一可用核 g 已於 R2-A 判死。")
    say("  ⟹ M12 至此三層 + 入射規則 + 目標層全部關閉，且無任何未檢之零參數構造。")


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/14_canonical_incidence.json , data/run_14_canonical_incidence.log")
    with open(os.path.join(ROOT, "data", "14_canonical_incidence.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False, default=lambda o: str(o))
    with open(os.path.join(ROOT, "data", "run_14_canonical_incidence.log"), "w") as f:
        f.write("\n".join(L) + "\n")