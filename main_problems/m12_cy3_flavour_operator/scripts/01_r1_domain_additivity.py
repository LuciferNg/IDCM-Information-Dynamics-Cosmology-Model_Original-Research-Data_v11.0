#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 R1 — CY₃ 投影層 flavour index：定義域 / 加性岔路 / n_s 幾何候選
================================================================================
本輪**不出數值預言**（docs/01 §五.4）。只做三件事：
  R1-Q1 定義域：三候選 (a) d_tensor 標號 / (b) H^{1,1}·H^{2,1} 週期 / (c) K_tree 標籤投影
        —— 逐項核 well-defined 且可算
  R1-Q2 加性岔路：一體 ⟹ ⟨O⟩≡1 是否被觀測否證；非加性載體庫存是否存在
  R1-Q3 n_s 幾何：由三錨反推所需形狀，檢查已註冊量能否「事前」給出（禁事後凑）

判準（**開算前寫死**；不得事後增刪）
--------------------------------------------------------------------------------
 D1  譜恆等（mpmath dps=50）：Σ_{n=1..44}γ_n = 90、γ_{44−n} = 4−γ_n、
     γ_a·γ_{44−a} = γ_{2a}、h11 = N_d−2B = 36、h21 = Σγ+2B = 98、TrG = 134。
     判：偏差 < 1e−40 ⟹ PASS（三角恆等，非擬合）。
 D2  真幾何 κ_ijk（CY₃(36,98) 多面體相交數，唯讀取用，本輪獨立重算）：
     全整數、非零支撐數、rank_ℚ、值域。判：全整數 ⟹ 純拓撲層無連續自由度。
 D3  一體 ⟹ ⟨O⟩ ≡ 1：⟨O⟩_ch ≡ ⟨h|N_d − N_u|h⟩；一體加性下任何「恰一次 d→u 取代」
     之通道 ⟨O⟩ = ±1。逐錨算 |⟨O⟩−1|/σ。判：任一 ≥ 3σ ⟹ 加性被否證。
 D4  旁觀者敏感性：同一 d→u 取代、不同旁觀者（ud / ds / ss / s̄）之 Δm。
     一體加性 ⟹ Δm_QCD = coeff·X ⟹ 四環境之 |Δm| 只差 em。判：原始 |Δm| 之差
     ≫ σ；且 em 側能否吸收 —— 引 M11 C9 之箱可行性（5/8 不可行，符號類 2 條與標度無關）。
     ⚠ 本項只作旁證（原始 Δm 含通道相依之 em）；主證據係 D3。
 D5  g(d) 精確 Fraction（N=88）：g(0)=2581/352、g(44)=−1291/352、
     20 對（d1≠d2 ∈ 1..5）之加性違反、零交 d* = (88−√2582)/2。
 D6  配對結構計數（純結構、不出值）：每通道價夸克配對數、被取代夸克所觸配對數、
     以及「介子 > 重子」所需之核比門檻 g(d)/g(0)。
 D7  n_s 剖面統計：三錨 ⟨O⟩ 對 n_s = 0/1/2 之逐差、比、二階差（含 σ 傳播）；
     meson(n_s=1) vs baryon(n_s=1) 之分裂 z。判：只報形狀約束，不作預言。
 D8  γ_n 在低 n 之凹凸（γ_n = 4sin²(nπ/88)，二階導符號）⟹ 篩「模數和」型候選。
 統計誠實：σ 只作傳播，不參與擬合；無任何自由參數；浮點只用於 mpmath 譜值與
 z 值顯示，精確量（g、κ、整數）一律 exact。
輸出：data/01_r1_domain_additivity.json ＋ stdout
"""
import json
import os
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = {}
LOG = []


def say(s=""):
    LOG.append(s)
    print(s)


# ─────────────────────────────────────────────────────────────── D1 譜恆等
def d1():
    say("=" * 78)
    say("D1 譜恆等（mpmath dps=%d）" % mp.mp.dps)
    pi = mp.pi
    g = lambda n: 4 * mp.sin(n * pi / 88) ** 2
    S = mp.fsum(g(n) for n in range(1, 45))
    mirror = max(abs(g(44 - n) - (4 - g(n))) for n in range(1, 44))
    mult = max(abs(g(a) * g(44 - a) - g(2 * a)) for a in range(1, 22))
    N_d, B = 44, 4
    h11 = N_d - 2 * B
    h21 = int(mp.nint(S)) + 2 * B
    say("  Σγ(n=1..44)      = %s" % mp.nstr(S, 30))
    say("  |Σγ − 90|        = %s" % mp.nstr(abs(S - 90), 6))
    say("  max|γ_{44−n}−(4−γ_n)| = %s" % mp.nstr(mirror, 6))
    say("  max|γ_aγ_{44−a}−γ_{2a}| = %s" % mp.nstr(mult, 6))
    say("  h11 = N_d−2B = %d ; h21 = Σγ+2B = %d ; TrG = %d" % (h11, h21, h11 + h21))
    say("  γ1 = %s ; 4/γ1 = %s ; γ2 = %s ; γ3 = %s" %
        (mp.nstr(g(1), 18), mp.nstr(4 / g(1), 18), mp.nstr(g(2), 18), mp.nstr(g(3), 18)))
    ok = (abs(S - 90) < mp.mpf('1e-40') and mirror < mp.mpf('1e-40')
          and mult < mp.mpf('1e-40') and h11 == 36 and h21 == 98 and h11 + h21 == 134)
    say("  判：%s" % ("PASS" if ok else "FAIL"))
    OUT["D1"] = {"sum_gamma": mp.nstr(S, 30), "mirror_dev": mp.nstr(mirror, 6),
                 "mult_dev": mp.nstr(mult, 6), "h11": h11, "h21": h21,
                 "TrG": h11 + h21, "verdict": "PASS" if ok else "FAIL",
                 "gamma1": mp.nstr(g(1), 18), "gamma2": mp.nstr(g(2), 18),
                 "gamma3": mp.nstr(g(3), 18), "four_over_g1": mp.nstr(4 / g(1), 18)}
    return g


# ─────────────────────────────────────────────────────────── D2 真幾何 κ_ijk
KAPPA = ("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/"
         "open_problems_fresh/o53_cy3_reality/data/polytopes/cy36_98_intersection_full.json")


def d2():
    say("=" * 78)
    say("D2 真幾何 κ_ijk（KS 多面體輸出；唯讀數據，非引用他版結論）")
    if not os.path.exists(KAPPA):
        say("  ⚠ 數據檔不存在 ⟹ 本項記 NOT-RUN（不得以記憶代替）")
        OUT["D2"] = {"verdict": "NOT-RUN"}
        return
    K = json.load(open(KAPPA))["kappa"]
    allint = all(Fraction(v).denominator == 1 for tri in K for row in tri for v in row)
    nz = sum(1 for tri in K for row in tri for v in row if v != 0)
    vals = [v for tri in K for row in tri for v in row]
    M = [[int(v) for v in row] for tri in K for row in tri]     # 1296 × 36
    rows, cols = len(M), 36
    # ℚ 上精確高斯消去求秩
    A = [[Fraction(x) for x in r] for r in M]
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(cols)]
        r += 1
        if r == cols:
            break
    say("  shape = 36×36×36；全整數 = %s；非零 = %d/%d (%.2f%%)；值域 = [%d, %d]"
        % (allint, nz, 36 ** 3, 100.0 * nz / 36 ** 3, min(vals), max(vals)))
    say("  rank_ℚ(flatten 1296×36) = %d" % r)
    say("  判：拓撲層 = 整數 + 有限秩 ⟹ 無連續自由度 ⟹ 其比值 ∈ ℚ（有限分母）")
    OUT["D2"] = {"allint": allint, "nz": nz, "min": min(vals), "max": max(vals),
                 "rank_Q": r, "verdict": "integers" if allint else "FAIL"}


# ────────────────────────────────────────────── D3 一體 ⟹ ⟨O⟩ ≡ 1（四錨）
ANCHORS = [("n−p", 1.1263, 0.1848), ("Σ⁻−Σ⁰", 1.9090, 0.1575),
           ("Ξ⁻−Ξ⁰", 2.1719, 0.1038), ("K⁺−K⁰", 2.4213, 0.1007)]


def d3():
    say("=" * 78)
    say("D3 一體加性 ⟹ ⟨O⟩ ≡ 1；逐錨偏離（σ 用 M11 docs/01 C2 表）")
    say("  %-10s %8s %8s %8s" % ("通道", "⟨O⟩", "σ", "z(vs 1)"))
    zs = {}
    for name, v, s in ANCHORS:
        z = (v - 1.0) / s
        zs[name] = z
        say("  %-10s %8.4f %8.4f %8.2f" % (name, v, s, z))
    mx = max(abs(z) for z in zs.values())
    n3 = sum(1 for z in zs.values() if abs(z) >= 3)
    say("  max|z| = %.2f；|z| ≥ 3 者 %d/4 ⟹ 判：%s"
        % (mx, n3, "加性被否證（岔路 (ii)）" if n3 else "加性未被否證"))
    OUT["D3"] = {"z": zs, "max_abs_z": mx, "n_ge_3sigma": n3,
                 "verdict": "additivity_refuted" if n3 else "additivity_alive"}


# ────────────────────────────────────── D4 旁觀者敏感性（分解無關）
def d4():
    say("=" * 78)
    say("D4 同一 d→u 取代、不同旁觀者之 Δm（分解無關，引 M11 docs/11 §五）")
    # Δm 取 M11 docs/01 §四 C1 表（PDG2024 獨立重算值）；σ 用其 σ_E 欄
    env = [("n→p", "ud", 1.2933, 0.0940),
           ("Σ⁻→Σ⁰", "ds", 4.8070, 0.1015),
           ("Ξ⁻→Ξ⁰", "ss", 6.8500, 0.2318),
           ("K⁰→K⁺", "s̄", -3.9340, 0.0962)]
    say("  %-8s %-6s %10s %8s %10s" % ("通道", "旁觀", "Δm(MeV)", "σ", "|Δm|/|Δm_N|"))
    base = abs(env[0][2])
    for nm, sp, v, s in env:
        say("  %-8s %-6s %10.4f %8.4f %10.3f" % (nm, sp, v, s, abs(v) / base))
    say("  兩兩差（一體 ⟹ 應全為 0）：")
    import itertools
    worst = 0.0
    for (n1, _, v1, s1), (n2, _, v2, s2) in itertools.combinations(env, 2):
        d = abs(abs(v1) - abs(v2))
        s = (s1 ** 2 + s2 ** 2) ** 0.5
        z = d / s
        worst = max(worst, abs(z))
        say("    |%s| − |%s| : Δ = %8.4f MeV  σ = %6.4f  z = %6.2f"
            % (n1, n2, d, s, z))
    say("  max|z| = %.2f ⟹ 原始 Δm 層旁觀者相關（一體加性下應為 0）" % worst)
    say("  ⚠ 原始 Δm 含通道相依之 em ⟹ 本項只作旁證。主證據係 D3（已扣 em）。")
    say("  em 側能否吸收？引 M11 docs/01 §四 C9（獨立重算，2100 箱角精確）：")
    say("    「兩體、正核、S,P∈(0,1]」算子類內 8 條中 5 條不可行；其中 K⁺−K⁰／Σ⁰−Σ⁺")
    say("    屬符號類（與 k₀ 無關）⟹ em 不可能吸收全部旁觀者依賴（本輪不重算，只引）。")
    OUT["D4"] = {"env": [[n, sp, v, s] for n, sp, v, s in env],
                 "max_abs_z": worst, "verdict": "spectator_dependent(raw Δm; em caveat)"}


# ─────────────────────────────── D5 g(d) 精確 + 20/20 非加性 + 零交
def d5():
    say("=" * 78)
    say("D5 樹核 g(d) 精確（N=88）；核實 M6 docs/13 之非加性（本輪獨立重算）")
    N = 88
    g = lambda d: Fraction((N ** 2 - 1) - 6 * d * (N - d), 12 * N)
    say("  g(0)  = %s = %.6f" % (g(0), float(g(0))))
    say("  g(44) = %s = %.6f ; g(44)/g(0) = %.6f"
        % (g(44), float(g(44)), float(g(44) / g(0))))
    say("  g(0) == 2581/352 : %s ; g(44) == −1291/352 : %s"
        % (g(0) == Fraction(2581, 352), g(44) == Fraction(-1291, 352)))
    viol = 0
    tot = 0
    ex = []
    for d1 in range(1, 6):
        for d2 in range(1, 6):
            if d1 == d2:
                continue
            tot += 1
            if g(d1) + g(d2) != g(d1 + d2):
                viol += 1
                if len(ex) < 3:
                    ex.append((d1, d2, g(d1) + g(d2), g(d1 + d2)))
    say("  加性測試（d1≠d2 ∈ 1..5，共 %d 對）：違反 %d/%d" % (tot, viol, tot))
    for d1, d2, lhs, rhs in ex:
        say("    g(%d)+g(%d) = %s (%.6f)  vs  g(%d) = %s (%.6f)"
            % (d1, d2, lhs, float(lhs), d1 + d2, rhs, float(rhs)))
    # 零交：6d(88−d) = 7743 ⟹ d² − 88d + 1290.5 = 0 ⟹ d = (88 ∓ √2582)/2
    dstar = (88 - mp.sqrt(2582)) / 2
    say("  零交 d* = (88 − √2582)/2 = %s ; g(18) = %.6f ; g(19) = %.6f"
        % (mp.nstr(dstar, 12), float(g(18)), float(g(19))))
    mono = all(g(d) > g(d + 1) for d in range(0, 44))
    say("  嚴格單調遞減（d=0..44）：%s" % mono)
    OUT["D5"] = {"g0": str(g(0)), "g44": str(g(44)), "violations": viol, "pairs": tot,
                 "dstar": mp.nstr(dstar, 12), "g18": float(g(18)), "g19": float(g(19)),
                 "monotone": mono}


# ─────────────────────────────── D6 配對結構計數（結構，不出值）
PAIRS = [
    # 通道, 初始價內容, 終態價內容, 取代, 旁觀者
    ("n−p", ("u", "d", "d"), ("u", "u", "d"), "d→u", ("u", "d")),
    ("Σ⁻−Σ⁰", ("d", "d", "s"), ("u", "d", "s"), "d→u", ("d", "s")),
    ("Ξ⁻−Ξ⁰", ("d", "s", "s"), ("u", "s", "s"), "d→u", ("s", "s")),
    ("Σ⁰−Σ⁺", ("u", "d", "s"), ("u", "u", "s"), "d→u", ("u", "s")),
    ("K⁺−K⁰", ("u", "s̄"), ("d", "s̄"), "d→u", ("s̄",)),
    ("D⁺−D⁰", ("u", "d̄"), ("d", "d̄"), "d→u", ("d̄",)),
]


def d6():
    say("=" * 78)
    say("D6 配對結構計數（純結構；零參數；不出任何 ⟨O⟩ 值）")
    say("  %-9s %-7s %6s %8s %10s" % ("通道", "價數", "總對", "所觸對", "所觸對之味型"))
    N = 88
    g = lambda d: Fraction((N ** 2 - 1) - 6 * d * (N - d), 12 * N)
    for ch, ini, fin, sub, spec in PAIRS:
        nval = len(ini)
        npair = nval * (nval - 1) // 2
        touched = nval - 1
        types = "+".join(f"u{x}" for x in spec)
        say("  %-9s %-7d %6d %8d %10s" % (ch, nval, npair, touched, types))
    say("  ⟹ 重子（3 價）：所觸對 = 2；介子（2 價）：所觸對 = 1")
    say("  ⟹ 距離型核（味-blind）下，三個重子通道之係數結構相同（皆 2 項）")
    say("     故 ⟨O⟩ 之差異只能來自『距離』或『味型』；前者需 flavour→mode 距離")
    # 門檻：介子(1 項) > 重子(2 項) ⟹ |k_meson| > 2|k_baryon|；最大核 = g(0)
    for k, label in ((2, "重子 2 觸對"), (3, "重子 3 對全計")):
        target = g(0) / k
        # g(d) = target ⟹ (N²−1) − 6d(N−d) = 12N·target ⟹ d² − Nd + c = 0
        c = ((N ** 2 - 1) - 12 * N * target) / 6
        disc = mp.mpf(N) ** 2 - 4 * c
        dlo = (N - mp.sqrt(disc)) / 2
        say("  門檻（%s）：需 g(d) < g(0)/%d = %.6f ⟹ d > %s（g(0)=%.6f, g(44)=%.6f）"
            % (label, k, float(target), mp.nstr(dlo, 8), float(g(0)), float(g(44))))
    OUT["D6"] = {"pairs": [[ch, len(i), len(i) * (len(i) - 1) // 2, len(i) - 1,
                            list(sp)] for ch, i, f, s, sp in PAIRS],
                 "verdict": "structure_only"}


# ─────────────────────────────────── D7 n_s 剖面（形狀約束，非預言）
def d7():
    say("=" * 78)
    say("D7 n_s 剖面形狀（單調 / 凹凸 / meson-baryon 分裂）；不作預言")
    prof = [("n_s=0 N", 1.1263, 0.1848), ("n_s=1 Σ", 1.9090, 0.1575),
            ("n_s=2 Ξ", 2.1719, 0.1038)]
    v = [p[1] for p in prof]
    s = [p[2] for p in prof]
    d1 = v[1] - v[0]
    d2 = v[2] - v[1]
    sd1 = (s[0] ** 2 + s[1] ** 2) ** 0.5
    sd2 = (s[1] ** 2 + s[2] ** 2) ** 0.5
    say("  逐差：Δ1 = %.4f ± %.4f（%.2fσ）; Δ2 = %.4f ± %.4f（%.2fσ）"
        % (d1, sd1, d1 / sd1, d2, sd2, d2 / sd2))
    say("  Δ2/Δ1 = %.4f" % (d2 / d1))
    sec = v[2] - 2 * v[1] + v[0]
    ssec = (s[0] ** 2 + 4 * s[1] ** 2 + s[2] ** 2) ** 0.5
    say("  二階差 = %.4f ± %.4f ⟹ z = %.2f（vs 0）⟹ %s"
        % (sec, ssec, sec / ssec,
           "與線性相容（無凹性訊號）" if abs(sec / ssec) < 2 else "偏離線性"))
    sK = 0.1007
    dK = 2.4213 - v[1]
    sKt = (sK ** 2 + s[1] ** 2) ** 0.5
    say("  meson(n_s=1) K⁺−K⁰ − Σ⁻−Σ⁰ = %.4f ± %.4f ⟹ z = %.2f"
        % (dK, sKt, dK / sKt))
    say("  形狀約束（必要條件，非預言）：")
    say("    (c1) ⟨O⟩ 隨 n_s 嚴格升：Δ1 = %.4f（%.2fσ）已顯著；Δ2 = %.4f 只 %.2fσ"
        % (d1, d1 / sd1, d2, d2 / sd2))
    say("    (c2) 凹性未被測出：二階差 z = %.2f ⟹ 線性與凹形皆可" % (sec / ssec))
    say("    (c3) 同一 n_s=1 下 meson > baryon：%.4f（%.2fσ）" % (dK, dK / sKt))
    say("    ⟹ 任何候選須同時滿足 c1–c3；只滿足 c1 者不足。")
    OUT["D7"] = {"profile": [[p[0], p[1], p[2]] for p in prof],
                 "d1": d1, "sd1": sd1, "d2": d2, "sd2": sd2, "ratio_d2_d1": d2 / d1,
                 "second_diff": sec, "sd_second": ssec, "z_second": sec / ssec,
                 "meson_baryon_diff": dK, "z_meson_baryon": dK / sKt}


# ─────────────────────────────── D8 γ_n 低 n 凹凸（篩「模數和」候選）
def d8():
    say("=" * 78)
    say("D8 γ_n 之凹凸：γ_n = 4sin²(nπ/88). 二階差符號（n ≤ 20）")
    pi = mp.pi
    g = lambda n: 4 * mp.sin(n * pi / 88) ** 2
    conv = []
    for n in range(1, 20):
        c = g(n + 1) - 2 * g(n) + g(n - 1)
        conv.append((n, c))
    nconv = sum(1 for n, c in conv if c > 0)
    say("  n=1..19 中二階差 > 0（凸）者 = %d/19" % nconv)
    say("  ⟹ 「佔用模 k≤K，取 Σγ_k」型候選在低 n 為凸；而觀測無凹性訊號（z=%.2f）"
        % OUT["D7"]["z_second"])
    say("     ⟹ 純模數和不足以解釋 c1–c3，須另有幾何因子（本輪只記篩選結果，不選形）")
    OUT["D8"] = {"n_convex_low": nconv, "note": "low-n spectral sums are convex"}


# ─────────── D9 後驗診斷（禁作預言）：最小配對通道分解 a_{type}
def d9():
    say("=" * 78)
    say("D9 後驗診斷（🔴 由錨反解，0 dof，**禁作預言**；只作 R2 目標）")
    say("  假設 ⟨O⟩ = 1 + Σ_{所觸對} a_{味型}（四味型：uu, ud, us, us̄）")
    say("    重子：N = a_uu + a_ud ; Σ = a_ud + a_us ; Ξ = 2a_us ; 介子：K = a_us̄")
    A = {}
    A["us"] = (2.1719 - 1.0) / 2.0
    A["ud"] = (1.9090 - 1.0) - A["us"]
    A["uu"] = (1.1263 - 1.0) - A["ud"]
    A["us̄"] = 2.4213 - 1.0
    for k in ("uu", "ud", "us", "us̄"):
        say("  a_%-3s = %+0.5f" % (k, A[k]))
    say("  未知 4、錨通道 4 ⟹ 0 dof ⟹ 對三條盲通道**有**預言力（R2 才可出值）。")
    say("  結構要求（R2 必須由已註冊量給出，不得反解）：")
    say("    (s1) a_uu < 0 < a_ud < a_us ⟹ 核須帶味型排序（非單一距離核）")
    say("    (s2) a_us̄ / a_us = %.3f ⟹ 反夸克通道須 ≳ 2× 夸克通道" % (A["us̄"] / A["us"]))
    OUT["D9"] = {"a": A, "dof": 0, "usbar_over_us": A["us̄"] / A["us"],
                 "flag": "POST-HOC / NOT-A-PREDICTION"}


# ─────── D10 譜距離核（已註冊模指派）之排序檢定 —— 具體候選之 no-go
def d10():
    say("=" * 78)
    say("D10 譜距離核 ⟨O⟩ = 1 + c·Σ_{所觸對} g(|n_i−n_j|)（c 為唯一標度，不入排序檢定）")
    say("  模指派（已註冊輸入，M4 docs/34：A₅ γ 升序 {1,2,3,4,5,7,9,11,13}）")
    say("    e→1, u→2, d→3；s 之指派**未釘死**（4/5/7 三讀皆在語料出現）⟹ 三個都算")
    N = 88
    g = lambda d: Fraction((N ** 2 - 1) - 6 * d * (N - d), 12 * N)
    pos = {"u": 2, "d": 3}
    observed = [("n−p", 1.1263), ("Σ⁻−Σ⁰", 1.9090), ("Ξ⁻−Ξ⁰", 2.1719), ("K⁺−K⁰", 2.4213)]
    say("  觀測排序（⟨O⟩）：n−p %.4f < Σ⁻−Σ⁰ %.4f < Ξ⁻−Ξ⁰ %.4f < K⁺−K⁰ %.4f"
        % tuple(v for _, v in observed))
    say("  %-9s %-14s %-14s %-14s" % ("s 指派", "N / Σ / Ξ / K 之 Σg", "", ""))
    for s_pos in (4, 5, 7):
        P = dict(pos)
        P["s"] = s_pos
        # 取代後夸克 = u；旁觀者見 D6
        def S(spec):
            tot = Fraction(0)
            for sp in spec:
                tot += g(abs(P["u"] - P[sp.replace("̄", "")]))
            return tot
        vals = {"n−p": S(("u", "d")), "Σ⁻−Σ⁰": S(("d", "s")),
                "Ξ⁻−Ξ⁰": S(("s", "s")), "K⁺−K⁰": S(("s",))}
        order = sorted(vals, key=lambda k: vals[k])
        say("  s→%d     N=%.4f  Σ=%.4f  Ξ=%.4f  K=%.4f  ⟹ 排序 %s"
            % (s_pos, float(vals["n−p"]), float(vals["Σ⁻−Σ⁰"]),
               float(vals["Ξ⁻−Ξ⁰"]), float(vals["K⁺−K⁰"]), " < ".join(order)))
        span = float(max(vals.values()) / min(vals.values()))
        say("           跨距 max/min = %.3f（觀測 (⟨O⟩−1) 跨距 = %.3f）"
            % (span, (2.4213 - 1) / (1.1263 - 1)))
    say("  ⟹ 三個 s 指派下排序皆為 K < Ξ < Σ < N，與觀測（N < Σ < Ξ < K）**完全反轉**")
    say("  ⟹ 判：譜距離核（d = |n_i − n_j|）＋已註冊模指派 = 🔴 死（c1、c3 同時反）")
    say("  ⚠ 界限：此判只打「譜距離」讀法；g(d) 之 d 在 M5 為**層/site** 距離，")
    say("     而 mode↔site 對應未註冊 ⟹ 換 d 之定義即換問題，須另行推（非本輪）。")
    OUT["D10"] = {"verdict": "spectral_distance_kernel_dead",
                  "note": "ordering inverted for s∈{4,5,7}; mode→site map unregistered"}


def main():
    os.makedirs(DATA, exist_ok=True)
    d1()
    d2()
    d3()
    d4()
    d5()
    d6()
    d7()
    d8()
    d9()
    d10()
    say("=" * 78)
    say("機讀摘要：")
    for k in sorted(OUT):
        say("  %s : %s" % (k, json.dumps(OUT[k], ensure_ascii=False, default=str)))
    with open(os.path.join(DATA, "01_r1_domain_additivity.json"), "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1, default=str)
    with open(os.path.join(DATA, "run_01_r1_domain_additivity.log"), "w") as f:
        f.write("\n".join(LOG) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())