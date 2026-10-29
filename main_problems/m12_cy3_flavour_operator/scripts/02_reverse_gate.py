#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2 方法閘（非 R2 推導本身）—— 反向必要性法（reverse-necessity）之可行性計量
2026-10-28

目的：使用者提出「由所有粒子所有味反向假設崩潰」取代剪枝式推導。本腳本只做
      計量判定：該法在 M12 現況下是否有足夠約束，以及典型崩潰判定是否 model-free。

⚠ 本腳本**不出任何 R2 預言**，只出「約束數 vs 結構位元數」與「崩潰判定之模型依賴度」。
⚠ D9 之 a_{type} 係 0 dof 後驗量；此處僅用作**敏感度對照**，禁作結構判據。

跑法：/usr/bin/python3 scripts/02_reverse_gate.py
"""
from __future__ import annotations
import json
import os
import sys
from mpmath import mp, mpf, sqrt

mp.dps = 30

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = {}
LINES = []


def say(s=""):
    LINES.append(s)
    print(s)


# ─────────────── 輸入：四錨（registered，contract 值）───────────────
# (name, ⟨O⟩, σ)  —— 一律「質量差通道 / 對應 em 分割」之讀數，見 docs/02 §4.5
# σ 取自 scripts/01 之 propag 值（0.1848 / 0.1575 / 0.1038 / 0.1007），
# 與 docs/02 所報 z = 0.68 / 5.77 / 11.29 / 14.11 一致。
ANCH = [
    ("n-p",      mpf("1.1263"), mpf("0.1848")),
    ("Sigma-0",  mpf("1.9090"), mpf("0.1575")),
    ("Xi-0",     mpf("2.1719"), mpf("0.1038")),
    ("K+0",      mpf("2.4213"), mpf("0.1007")),
]
NAME = {a[0]: a[1] for a in ANCH}
SIG = {a[0]: a[2] for a in ANCH}


def z(x, s):
    return x / s if s != 0 else mpf("inf")


def combine(*terms):
    """terms = [(value, sigma), ...] 線性相加，獨立誤差。"""
    v = sum(t[0] for t in terms)
    s = sqrt(sum(t[1] ** 2 for t in terms))
    return v, s


# ─────────── G1 兩個層次：model-free 讀數 vs D9 ansatz 分解 ───────────
def g1():
    say("=" * 80)
    say("G1  兩層量：model-free 讀數 vs D9 ansatz 分解（0 dof 後驗）")
    say("")
    say("  (1) model-free：四錨相對 1 之增量，唔作任何配對假設")
    raw = {}
    for n, v, s in ANCH:
        d, sd = v - 1, s
        raw[n] = (d, sd, z(d, sd))
        say("      %-8s Δ = %+0.6f ± %0.6f   z(Δ=0) = %+0.3f" % (n, d, sd, z(d, sd)))
    say("")
    say("  (2) D9 ansatz：N=uu+ud, Σ=ud+us, Ξ=2us, K=usbar（四未知四錨 ⟹ dof=0）")
    a = {}
    a["us"] = ((NAME["Xi-0"] - 1) / 2, SIG["Xi-0"] / 2)
    a["ud"] = combine(((NAME["Sigma-0"] - 1), SIG["Sigma-0"]), (-a["us"][0], a["us"][1]))
    a["uu"] = combine(((NAME["n-p"] - 1), SIG["n-p"]), (-a["ud"][0], a["ud"][1]))
    a["usbar"] = ((NAME["K+0"] - 1), SIG["K+0"])
    for k in ("uu", "ud", "us", "usbar"):
        v, s = a[k]
        say("      a_%-5s = %+0.6f ± %0.6f   z(=0) = %+0.3f  ⟹ %s"
            % (k, v, s, z(v, s),
               "穩健非零" if abs(z(v, s)) >= 3 else "⚠ 不可作約束"))
    r = a["usbar"][0] / a["us"][0]
    sr = r * sqrt((a["usbar"][1] / a["usbar"][0]) ** 2 + (a["us"][1] / a["us"][0]) ** 2)
    say("      a_usbar/a_us = %0.6f ± %0.6f   z(=1) = %+0.3f   z(=2) = %+0.3f"
        % (r, sr, z(r - 1, sr), z(r - 2, sr)))
    OUT["G1"] = {
        "raw": {k: [str(v[0]), str(v[1]), str(v[2])] for k, v in raw.items()},
        "ansatz": {k: [str(v[0]), str(v[1]), str(z(v[0], v[1]))] for k, v in a.items()},
        "ratio": [str(r), str(sr), str(z(r - 1, sr)), str(z(r - 2, sr))],
    }
    return a, raw


# ─────────── G2 崩潰測試：反夸克 = ±夸克 是否成立 ───────────
def g2(a, raw):
    say("")
    say("=" * 80)
    say("G2  崩潰測試（collapse test）：反夸克通道 = ∓夸克通道？")
    say("")
    say("  T-A  s̄ = −s（「反粒子取反號」之自然讀法）⟹ 預言 a_usbar = −a_us")
    num, sig = combine((a["usbar"][0], a["usbar"][1]), (a["us"][0], a["us"][1]))
    zA = z(num, sig)
    say("      偏離 = %+0.6f ± %0.6f ⟹ z = %+0.3f   %s"
        % (num, sig, zA, "🔴 崩潰" if abs(zA) >= 5 else "未崩潰"))
    say("  T-B  s̄ = +s（反夸克同類）⟹ 預言 a_usbar = +a_us")
    num, sig = combine((a["usbar"][0], a["usbar"][1]), (-a["us"][0], a["us"][1]))
    zB = z(num, sig)
    say("      偏離 = %+0.6f ± %0.6f ⟹ z = %+0.3f   %s"
        % (num, sig, zB, "🔴 崩潰" if abs(zB) >= 5 else "未崩潰"))
    say("  ⟹ %s ⟹ 「反夸克 = Z₂ 號複本」此一類（含 M5 docs/03 §六.3 之 s_A s_B B^-d）"
        % ("兩讀皆死" if (abs(zA) >= 5 and abs(zB) >= 5) else "未兩讀皆死"))
    say("     在 D9 ansatz 內**不可行**；須 ℓ_s̄ ≠ ℓ_s（層位移）或獨立幅度。")
    say("")
    say("  ⚠ 但同一結論之 model-free 版本：K 與 Ξ 各含恰好一個 s，原始 ⟨O⟩ 差 = ")
    num, sig = combine((NAME["K+0"], SIG["K+0"]), (-NAME["Xi-0"], SIG["Xi-0"]))
    zR = z(num, sig)
    say("      ⟨O⟩_K − ⟨O⟩_Ξ = %+0.6f ± %0.6f ⟹ z = %+0.3f  %s"
        % (num, sig, zR, "（<3σ，**不構成崩潰**）" if abs(zR) < 3 else ""))
    say("      ⟹ %0.1fσ／%0.1fσ 全部由『Ξ 算 2 個 s、K 算 1 個 s』之歸一化選擇產生；"
        % (float(abs(zA)), float(abs(zB))))
    say("         該歸一化係 ansatz 的一部分，**未經推導** ⟹ 崩潰判定 **非 model-free**。")
    OUT["G2"] = {"T_A_z": str(zA), "T_B_z": str(zB), "raw_K_minus_Xi": [str(num), str(sig), str(zR)]}


# ─────────── G3 位元計數閘：約束數 vs 結構位元數 ───────────
def g3():
    say("")
    say("=" * 80)
    say("G3  位元計數閘（RN4）：model-free 約束數 vs 結構位元數")
    say("")
    say("  model-free 約束（唔依賴 ansatz）：")
    say("    C1  ⟨O⟩_np − 1 > 0            (+0.68σ … 僅方向)")
    say("    C2  ⟨O⟩_Σ − 1 = +0.909±0.158  (5.75σ)")
    say("    C3  ⟨O⟩_Ξ − 1 = +1.172±0.104  (11.29σ)")
    say("    C4  ⟨O⟩_K − 1 = +1.421±0.101  (14.11σ)")
    say("    ⟹ 有效約束 = 4 個數（C1 方向僅 0.68σ，弱）")
    say("")
    say("  結構位元（D9 四味型）= 4 個幅度；若再加 Z₂ 號 = 4 bit")
    say("    ⟹ 4 約束 vs 4 幅度：**剛好飽和 → dof = 0**")
    say("    ⟹ 4 約束 vs 4 幅度 + 4 bit = 8 位元：**欠定 4 位元**")
    say("  ⟹ 反向法在現況**無餘裕**：任何崩潰判定的自由度都由 ansatz 補上 = 滑向擬合。")
    say("")
    say("  唯一不佔連續位元之額外約束 = 二值拓撲/對稱約束：")
    say("    ① G-parity（π⁰: ⟨O⟩ = 0，已作錨輸入）")
    say("    ② A₅ 圖論（見 G4）")
    say("    ③ 三條盲通道（BLIND，不可用）")
    OUT["G3"] = {"free_constraints": 4, "amplitudes": 4, "z2_bits": 4, "dof": 0, "shortfall_bits": 4}


# ─────────── G4 A₅ 圖論能否提供味↔頂點指派 ───────────
def g4():
    say("")
    say("=" * 80)
    say("G4  A₅/RH 圖論可否提供味↔頂點（或味→號）之典範指派？")
    say("")
    say("  已註冊事實（polyhedral-graph-analysis §4/§5c/§15(h)，非本輪新算）：")
    say("    V=60, E=120, 4-regular, 面 20△+30□+12⬠, Euler 60−120+62=2 ✅")
    say("    邊集可分割為 20 個邊不相交三角形 + 12 個邊不相交五邊形 ⟹ 3·20 + 5·12 = %d ✅" % (3 * 20 + 5 * 12))
    say("    ⟹ MaxCut ≤ E − #odd-faces = 120 − 32 = %d（上界定理）" % (120 - 32))
    say("    ⟹ MaxCut = 88 = 2N_d（O137 MILP 精確達界）⟹ Ising H = −(88−32) = −56 ✅")
    say("    飽和 2-色 = **30 員等結構族**（60 root BFS 去重 mod 全域翻轉）")
    say("")
    say("  ★ 決定性：RH 係 vertex-transitive（Archimedean solid，|Aut| = 120，")
    say("    對頂點傳遞）⟹ 對任意兩頂點存在自同構互換")
    say("    ⟹ 任一具體 2-色 c 與其像 g·c 皆為容許飽和色 ⟹ 同一頂點在族內可為 + 亦可為 −")
    say("    ⟹ **冇典範（canonical）味↔頂點 / 味→號 指派**")
    say("    ⟹ skill §15(d) 同一定理語言：「圖不能選；orbit-mates 之差異必來自非圖輸入」")
    say("")
    say("  再者：MaxCut = 88 唔唯一 —— 4-regular A₅ Cayley 圖共 5 個非同構，其中 **2 個**")
    say("    都得 MaxCut = 88（唯一判別者 = planarity）⟹ 88 本身唔能釘死 RH")
    say("    ⟹ 「88 = 2N_d 係結構」仍未證（docs/02 §九 已標 ⚠待核，本輪維持未證）")
    say("")
    say("  ⟹ 結論：A₅ 圖論可提供 **全域一個 bit**（互補色），但唔可以提供**每味一個 bit**")
    say("         除非另有味↔頂點規則 = R1 所稱「缺橋」。")
    OUT["G4"] = {"V": 60, "E": 120, "odd_faces": 32, "maxcut_bound": 88,
                 "maxcut_exact": 88, "ising_H": -56, "saturated_colourings": 30,
                 "vertex_transitive": True, "aut_order_geometric": 120,
                 "maxcut88_graphs_noniso": 2, "canonical_flavour_map": False}


# ─────────── G5 判詞彙總 ───────────
def g5():
    say("")
    say("=" * 80)
    say("G5  判詞（本腳本只判方法可行性，不出 R2 預言）")
    say("")
    say("  ✅ 反向必要性法**可行**，但四閘必須齊：")
    say("     RN1 只准用 model-free 量（禁 D9 之 a_type 作判據）")
    say("     RN2 每個假設須為真二分（X vs ¬X 窮盡），否則退化成剪枝")
    say("     RN3 崩潰判據須二值（存在/不存在），禁 χ² 級")
    say("     RN4 位元計數：約束數 ≥ 結構位元數，否則宣告數據側不能釘死")
    say("")
    say("  🔴 M12 現況在 RN4 失守：4 約束 vs 4 幅度（dof=0）或 8 位元（欠定 4）")
    say("  🔴 A₅ 圖論在 vertex-transitivity 下唔提供每味 bit（G4）")
    say("  ⟹ 反向法之合法產物 = 把問題**收窄為一**：缺嘅係「味↔結構」對應（橋），")
    say("     而唔係「核」；此與 R1 §4.6 二分之結論一致。")
    say("  ⟹ 反向法不能取代正向構造；能取代嘅係「候選空間之窮盡枚舉 + 不變量切割」。")
    OUT["G5"] = {"reverse_necessity_feasible": True, "RN4_fails_now": True,
                 "A5_gives_per_flavour_bit": False,
                 "verdict": "reverse narrows to the bridge, does not replace construction"}


def main():
    a, raw = g1()
    g2(a, raw)
    g3()
    g4()
    g5()
    os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)
    with open(os.path.join(ROOT, "data", "02_reverse_gate.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    say("")
    say("[written] data/02_reverse_gate.json , data/run_02_reverse_gate.log")

    with open(os.path.join(ROOT, "data", "run_02_reverse_gate.log"), "w") as f:
        f.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()