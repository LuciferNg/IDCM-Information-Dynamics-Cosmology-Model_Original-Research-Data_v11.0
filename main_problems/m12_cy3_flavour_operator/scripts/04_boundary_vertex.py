#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2 —— 崩潰邊界交集（collapse-boundary intersection）
2026-10-28

使用者方法論（釘死版）：
  「假設宇宙只有一個解 ⟹ 每個崩潰邊界 = 可行域的一條邊 ⟹ 解 = 所有不崩潰邊界之交集（點）」

形式化：
  參數空間 ℝ^m（m = 味型幅度數）
  錨通道 ⟹ 等式（超平面）；結構裁決 ⟹ 半空間（不等式，邊界 = 崩潰面）
  唯一解 ⟺ (i) 獨立等式數 = m − 1（收窄成線） ∧ (ii) 再有一條結構邊界與該線相交成點

本腳本做三件事（全部確定性、可重跑）：
  H1 錨系統之秩與零空間（精確 Fraction）→ 確認係「線」而非「點」
  H2 半空間邊界清單 → 線被切成段；段寬之 σ；D9 之點在段外否
  H3 t-無關之崩潰邊界：註冊模式指派 ＋ 衰減核 ⟹ 排序約束 ⟹ 與實測比較

跑法：/usr/bin/python3 scripts/04_boundary_vertex.py
"""
from __future__ import annotations
import json
import os
from fractions import Fraction as F
from mpmath import mp, mpf, sqrt

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


# ── 錨（⟨O⟩, σ）── 與 scripts/01/02 一致
ON = mpf("1.1263"); SN = mpf("0.1848")
OS = mpf("1.9090"); SS = mpf("0.1575")
OX = mpf("2.1719"); SX = mpf("0.1038")
OK = mpf("2.4213"); SK = mpf("0.1007")

dN, dS, dX, dK = ON - 1, OS - 1, OX - 1, OK - 1


# ═══════════ H1：秩與零空間（精確）═══════════
def h1():
    say("=" * 80)
    say("H1  錨系統之秩（精確 Fraction）—— 確認交集係「線」唔係「點」")
    say("")
    say("  變數序 = (a_uu, a_ud, a_us, a_ss, a_usbar)")
    say("    N(udd) : a_uu + a_ud                       = ⟨O⟩_N − 1")
    say("    Σ(uds) :         a_ud + a_us               = ⟨O⟩_Σ − 1")
    say("    Ξ(uss) :                a_us + a_ss        = ⟨O⟩_Ξ − 1")
    say("    K(u s̄) :                       a_usbar     = ⟨O⟩_K − 1")
    A = [
        [F(1), F(1), F(0), F(0), F(0)],
        [F(0), F(1), F(1), F(0), F(0)],
        [F(0), F(0), F(1), F(1), F(0)],
        [F(0), F(0), F(0), F(0), F(1)],
    ]
    # 精確高斯消元求秩
    M = [row[:] for row in A]
    rank = 0
    piv = []
    for c in range(5):
        p = None
        for r in range(rank, len(M)):
            if M[r][c] != 0:
                p = r
                break
        if p is None:
            continue
        M[rank], M[p] = M[p], M[rank]
        pv = M[rank][c]
        M[rank] = [x / pv for x in M[rank]]
        for r in range(len(M)):
            if r != rank and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[rank])]
        piv.append(c)
        rank += 1
    say("  秩 = %d / 5 ⟹ 零空間維度 = %d ⟹ 解集係 **%d 維仿射子空間（線）**"
        % (rank, 5 - rank, 5 - rank))
    say("  自由參數 t ≡ a_ss：")
    say("    a_usbar = %+0.6f" % dK)
    say("    a_us    = %+0.6f − t" % dX)
    say("    a_ud    = %+0.6f + t" % (dS - dX))
    say("    a_uu    = %+0.6f − t" % (dN - (dS - dX)))
    say("  ⟹ 單靠錨（= 觀測）永遠得線，唔會得點。點必須由**結構邊界**提供。")
    OUT["H1"] = {"rank": rank, "nullity": 5 - rank,
                 "line": {"usbar": str(dK), "us": "%s - t" % dX,
                          "ud": "%s + t" % (dS - dX), "uu": "%s - t" % (dN - (dS - dX))}}


# ═══════════ H2：半空間邊界 → 段；D9 之點在段外否 ═══════════
def h2():
    say("")
    say("=" * 80)
    say("H2  結構半空間邊界（正性）對該線之切割")
    say("")
    # 正性邊界：a_uu>=0, a_ud>=0, a_us>=0, a_ss>=0
    lo_ud = -(dS - dX)          # a_ud>=0 ⟹ t >= −(dS−dX)
    hi_uu = (dN - (dS - dX))    # a_uu>=0 ⟹ t <= dN−dS+dX
    hi_us = dX                  # a_us>=0 ⟹ t <= dX
    say("    B1  a_ud ≥ 0 ⟹ t ≥ %+0.6f" % lo_ud)
    say("    B2  a_uu ≥ 0 ⟹ t ≤ %+0.6f" % hi_uu)
    say("    B3  a_us ≥ 0 ⟹ t ≤ %+0.6f（較鬆）" % hi_us)
    say("    B4  a_ss ≥ 0 ⟹ t ≥ 0")
    lo = max(lo_ud, mpf(0))
    hi = min(hi_uu, hi_us)
    width = hi - lo
    say("  ⟹ 交集 = 段 [%+0.6f, %+0.6f]，寬 = %+0.6f" % (lo, hi, width))
    sw = sqrt(SS ** 2 + SX ** 2 + SN ** 2)
    say("     段寬 = ⟨O⟩_N − 1 = %+0.6f ± %0.6f ⟹ z(=0) = %+0.3f"
        % (width, SN, width / SN))
    say("     ⟹ 段存在但**未分辨成一點**（0.68σ）⟹ 現有邊界不足，須加結構邊界。")
    say("")
    t9 = dX / 2
    say("  D9 之點（a_ss = a_us ⟹ t = %0.6f）：" % t9)
    viol = t9 - hi
    say("     t − t_max = %+0.6f ⟹ **落在段外**（a_uu = %+0.6f < 0）"
        % (viol, (dN - (dS - dX)) - t9))
    say("     ⟹ 若接受「同味對正核（g(d)>0 於 d<18.5933）」，D9 之點崩潰。")
    say("     ⟹ 但注意：此判決係 **0.79σ 級**（見 docs/02 §十），不可作強判。")
    OUT["H2"] = {"lo": str(lo), "hi": str(hi), "width": str(width),
                 "width_z": str(width / SN), "d9_t": str(t9), "d9_outside": True,
                 "d9_violation": str(viol)}


# ═══════════ H3：t-無關之崩潰邊界（衰減核 vs 註冊模式）═══════════
def h3():
    say("")
    say("=" * 80)
    say("H3  t-無關崩潰邊界：衰減徑向核 ＋ 註冊模式指派 ⟹ 排序約束")
    say("")
    say("  已註冊：味↔模式 u→2、s∈{4,5,7}（O137 裁決 → O209 d14 → O233；見 docs/02 §4.5）")
    say("  任何**衰減**徑向核 k(d)（含 M5 docs/03 §六.3 之 B^(−d_tree)，B = 4）：")
    say("    d(u,u) = 0 < d(u,s) ⟹ k_uu ≥ k_us ⟹ a_uu ≥ a_us")
    say("")
    d_ = (dN - dS)                       # a_uu − a_us = ⟨O⟩_N − ⟨O⟩_Σ
    sd = sqrt(SN ** 2 + SS ** 2)
    say("  代數：a_uu − a_us = (a_N − a_ud) − (a_S − a_ud) = a_N − a_S")
    say("      ⟹ 與 t **完全無關**：a_uu − a_us ≡ ⟨O⟩_N − ⟨O⟩_Σ")
    say("  實測：⟨O⟩_N − ⟨O⟩_Σ = %+0.6f ± %0.6f ⟹ z = %+0.3f" % (d_, sd, d_ / sd))
    say("  ⟹ 衰減核要求 ≥ 0，實測 = %+0.4f（%.2fσ 負）⟹ **衰減核全族崩潰**"
        % (d_, abs(d_ / sd)))
    say("     受影響：K_tree 徑向族 1/b^ℓ、M5 docs/02 之 B^(−d_tree)、" )
    say("             M5 docs/03 §六.3 之 s_A s_B B^(−|ℓ_A−ℓ_B|)（T4）")
    say("  ⟹ 此判決 **唔依賴 a_ss = a_us**，亦唔依賴段內位置 ⟹ 真·t-無關。")
    say("")
    say("  唯一存活方向：味耦合必須**隨味距遞增**（k_us > k_uu）")
    say("      ⟹ 味部分係「荷」型（label magnitude），唔係「距離衰減」型")
    say("      ⟹ 度規只能定標度，唔能定味序 ⟹ 與 R1 §2.3「標籤⊗度規」相容，且進一步指定：")
    say("         標籤部分須**單調遞增**，而 M5 之 Z₂ 標籤（±1）只給號、唔給遞增幅度 ⟹ 不足。")
    say("")
    # 遞增方向之可數邊界
    say("  遞增方向之可數邊界（供 R2 飽和計數）：")
    say("    M1 s 數單調：Δ1 = Σ−N = %+0.4f ± %+0.4f (%.2fσ)"
        % (dS - dN, sd, (dS - dN) / sd))
    d2 = dX - dS; s2 = sqrt(SX ** 2 + SS ** 2)
    say("    M2 s 數單調：Δ2 = Ξ−Σ = %+0.4f ± %+0.4f (%.2fσ)" % (d2, s2, d2 / s2))
    dKX = dK - dX; sKX = sqrt(SK ** 2 + SX ** 2)
    say("    M3 s̄ > 2s   ：K−Ξ  = %+0.4f ± %+0.4f (%.2fσ)" % (dKX, sKX, dKX / sKX))
    say("    ⟹ 3 條方向性邊界，全部同向（遞增）；但只有 M1 過 3σ。")
    OUT["H3"] = {"a_uu_minus_a_us": str(d_), "z": str(d_ / sd),
                 "decaying_kernel_dead": True, "t_independent": True,
                 "M1": [str(dS - dN), str(sd), str((dS - dN) / sd)],
                 "M2": [str(d2), str(s2), str(d2 / s2)],
                 "M3": [str(dKX), str(sKX), str(dKX / sKX)]}


def main():
    h1()
    h2()
    h3()
    with open(os.path.join(ROOT, "data", "04_boundary_vertex.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    say("")
    say("[written] data/04_boundary_vertex.json , data/run_04_boundary_vertex.log")

    with open(os.path.join(ROOT, "data", "run_04_boundary_vertex.log"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()