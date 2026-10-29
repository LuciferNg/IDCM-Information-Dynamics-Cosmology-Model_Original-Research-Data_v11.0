#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2-B —— 可分性、解析度預算、與符號庫密度閘
2026-10-28

三問：
  G5 可分（factorised / product）多體形式：Σ² = N·Ξ 之「幾何中項」恆等式可否成立？
  G6 解析度：現有四錨能把 t 收到幾窄？窗口 0.126300 相對 σ 係幾多？
  G7 符號庫密度：窗口 [0.262900, 0.389200] 內有幾多「已註冊數字」？
     覆蓋率 ≥ 30% ⟹ 「落入窗口」不可作命中判據（symbol-library explosion）。

跑法：/usr/bin/python3 scripts/08_resolution_density.py
"""
from __future__ import annotations
import json
import os
from fractions import Fraction as Fr

from mpmath import mp, mpf, sqrt, pi, sin, cos, acos

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


DN, SN = mpf("0.1263"), mpf("0.1848")
DS, SS = mpf("0.9090"), mpf("0.1575")
DX, SX = mpf("1.1719"), mpf("0.1038")
DK, SK = mpf("1.4213"), mpf("0.1007")
ON, OS, OX, OK = DN + 1, DS + 1, DX + 1, DK + 1
T_LO, T_HI = mpf("0.262900"), mpf("0.389200")

PHI = (1 + mp.sqrt(5)) / 2


# ─────────── G5 可分性（幾何中項）───────────
def g5():
    say("=" * 80)
    say("G5  可分／乘積形式之恆等式檢定：Σ² = N·Ξ")
    say("")
    say("  若響應可分解為（輕環境）（重環境）之積，則")
    say("      N = U·D²   Σ = U·D·S   Ξ = U·S²   ⟹  Σ² = N·Ξ  （精確恆等式）")
    say("")
    pred = sqrt(ON * OX)
    say("  預言 Σ = √(N·Ξ) = %0.6f" % pred)
    say("  實測 Σ = %0.6f ± %0.6f   差 = %+0.6f" % (OS, SS, OS - pred))
    # 傳播：σ(√(NΞ))
    s_pred = pred * sqrt((SN / ON) ** 2 + (SX / OX) ** 2) / 2
    say("  σ(√(NΞ)) = %0.6f ⟹ z = %+0.3f"
        % (s_pred, (OS - pred) / sqrt(s_pred ** 2 + SS ** 2)))
    say("")
    say("  代數比 N·Ξ / Σ² = %0.6f（若可分應 = 1）" % (ON * OX / OS ** 2))
    say("  ⟹ 代數上差 %.1f%%，但 z 僅 ≲2σ ⟹ 按 V3 **不可裁決**，只作趨勢記錄。"
        % (abs(float(ON * OX / OS ** 2) - 1) * 100))
    OUT["G5"] = {"pred_sigma": str(pred), "obs_sigma": str(OS),
                 "ratio": str(ON * OX / OS ** 2),
                 "z": str((OS - pred) / sqrt(s_pred ** 2 + SS ** 2)),
                 "decisive": False}


# ─────────── G6 解析度預算 ───────────
def g6():
    say("")
    say("=" * 80)
    say("G6  解析度預算：四錨能收幾窄？")
    say("")
    w = T_HI - T_LO
    say("  正性段寬 = %0.6f = ⟨O⟩_N − 1 = %0.6f ± %0.6f ⟹ 段寬本身 %0.3fσ"
        % (w, DN, SN, w / SN))
    say("  ⟹ 兩條邊界（a_uu=0、a_ud=0）在統計上**近乎重合**：")
    say("     若 N 通道無響應（⟨O⟩_N = 1），則兩邊界完全重合 ⟹ 段退化為一點。")
    say("     實測 ⟨O⟩_N − 1 = %+0.6f ± %0.6f ⟹ 距「無響應」僅 %0.2fσ。"
        % (DN, SN, DN / SN))
    say("")
    say("  ⟹ 承重之錨 = n−p（唯一分開兩條正性邊界者），但其自身僅 %0.2fσ。" % (DN / SN))
    say("  ⟹ 去掉 n−p：3 錨 / 5 未知 ⟹ 2 維解集 ⟹ t 完全不受限。")
    say("  ⟹ 保留 n−p：t ∈ [%0.4f, %0.4f]，而 t 之 σ ≈ %0.3f ⟹ 未分辨。"
        % (T_LO, T_HI, SN))
    say("")
    say("  要把窗口收窄至 4× 以下（可測），需要 σ(t) ≲ %0.4f；" % (w / 4))
    say("  現時 σ 主項來自 n−p 之 em 分割讀出（σ = %0.4f，佔相對 146%%）。" % SN)
    OUT["G6"] = {"window": str(w), "window_in_sigma": str(w / SN),
                 "n_minus_p_z": str(DN / SN), "load_bearing": "n-p",
                 "sigma_t": str(SN), "unresolved": True}


# ─────────── G7 符號庫密度閘 ───────────
def catalogue():
    C = set()
    # (1) 有理數 p/q, q<=16
    for q in range(2, 17):
        for p in range(1, q):
            C.add(mpf(p) / q)
    # (2) φ 冪
    for n in range(1, 9):
        C.add(PHI ** (-n))
    # (3) Lucas / Fibonacci 比
    luc = [2, 1, 3, 4, 7, 11, 18, 29, 47]
    fib = [1, 1, 2, 3, 5, 8, 13, 21, 34]
    for a in luc:
        for b in fib:
            if b:
                v = mpf(a) / b
                if 0 < v < 1:
                    C.add(v)
    # (4) 譜 γ_n/4 = sin²(nπ/88)
    for n in range(1, 45):
        C.add(sin(n * pi / 88) ** 2)
    # (5) 計數比 k/88, k/90, k/45, k/60, k/41, k/19
    for den in (88, 90, 45, 60, 41, 19, 32, 44, 16):
        for k in range(1, den):
            C.add(mpf(k) / den)
    # (6) √5 型
    C.add(mp.sqrt(5) - 2)
    C.add((3 - mp.sqrt(5)) / 2)
    C.add((mp.sqrt(5) - 1) / 3)
    return C


def g7():
    say("")
    say("=" * 80)
    say("G7  符號庫密度閘（skill §15(f)：宣稱結構數前必報密度）")
    C = catalogue()
    say("  符號庫大小 = %d（有理 p/q q≤16 ／ φ⁻ⁿ ／ L/F 比 ／ sin²(nπ/88) ／ k/88,k/90,k/45,k/60,k/41,k/19 ／ √5 型）"
        % len(C))
    inside = sorted(v for v in C if T_LO <= v <= T_HI)
    say("  落於窗口 [%0.4f, %0.4f] 之成員 = %d" % (T_LO, T_HI, len(inside)))
    for v in inside:
        say("      %0.6f" % v)
    # 覆蓋率：窗口內任一點距最近成員 ≤ halfσ 之比例（以 0.02 網格數值估計）
    step = mpf("0.002")
    tol = SN / 2
    tot = hit = 0
    x = T_LO
    while x <= T_HI:
        tot += 1
        if any(abs(x - v) <= tol for v in C):
            hit += 1
        x += step
    cov = mpf(hit) / tot
    say("  ±σ/2 (=%0.4f) 覆蓋率 = %0.1f%%  ⟹ %s"
        % (tol, cov * 100,
           "🔴 密度過高，「落入窗口」不可作命中（symbol-library explosion）"
           if cov > mpf("0.3") else "🟡 密度可接受"))
    say("  以 ±0.02 嚴格容差計覆蓋率 = %0.1f%%"
        % (100 * mpf(sum(1 for i in range(int((T_HI - T_LO) / mpf("0.002")) + 1)
                         if any(abs(T_LO + i * mpf("0.002") - v) <= mpf("0.02") for v in C)))
           / (int((T_HI - T_LO) / mpf("0.002")) + 1)))
    OUT["G7"] = {"catalogue_size": len(C), "n_inside": len(inside),
                 "inside": [str(v) for v in inside],
                 "coverage_sigma_half": str(cov),
                 "density_gate_fails": bool(cov > mpf("0.3"))}


def main():
    g5()
    g6()
    g7()
    say("")
    say("=" * 80)
    say("R2-B 判詞")
    say("")
    say("  🟡 G5 可分形式：代數差 22%，但 z ≈ 2σ ⟹ 不可裁決（V3）")
    say("  🔴 G6 解析度：窗口寬僅 0.68σ ⟹ t 未分辨；承重者係最弱錨（n−p, 0.68σ）")
    say("  🔴 G7 密度：窗口內符號庫覆蓋率過高 ⟹ 「落入窗口」= 無效命中")
    say("")
    say("  ⟹ R2 不可靠『預測 t 落入窗口』收工。要麼(A)改讀數使 σ(t) 收窄 ≥4×，")
    say("     要麼(B)結構直接出精確值**同時**預言三條盲通道（V5）。")


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/08_resolution_density.json , data/run_08_resolution_density.log")
    with open(os.path.join(ROOT, "data", "08_resolution_density.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    with open(os.path.join(ROOT, "data", "run_08_resolution_density.log"), "w") as f:
        f.write("\n".join(L) + "\n")