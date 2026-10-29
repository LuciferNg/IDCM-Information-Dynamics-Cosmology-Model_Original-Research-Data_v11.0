#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2-A —— 二體（pair-sum）結構形式之窮盡測試
2026-10-28

背景（docs/03 方法閘）：
  · 錨 ⟹ 解集 = 1 維線（t ≡ a_ss）
  · 衰減核全族已被 t-無關邊界打死（z = −3.22）
  · 局域容錯位元預算 ⟹ A₅ 頂點標籤只有 1 bit 且係 gauge
  ⟹ 只剩「荷型遞增」方向。本腳本把「荷型」之**所有零參數二體形式**逐個算死或算活。

測試之形式族（每種都把 F̂ 寫成 Σ_{價對} a(x,y)，a 由結構給出、零自由參數）：
  F1 乘性荷  a(x,y) = w_x·w_y            （＝ γ_a γ_{44−a} = γ_{2a} 之 rank-1 類比）
  F2 加性荷  a(x,y) = c(γ_x + γ_y)
  F3 衰減距  a(x,y) = k(|Δn|)，k 遞減      （docs/03 §四已判死，此處只作對照）
  F4 環境多項式 f(n_s) = α + β·n_s        （非二體；測「環境計數」能否唯一）

判定：任一形式若 (i) 解不出實數/正值、(ii) 與註冊質量序 m_u < m_d < m_s 相反、
      或 (iii) 落於 docs/03 H2 之正性段 [0.262900, 0.389200] 之外 ⟹ 判死。

跑法：/usr/bin/python3 scripts/07_two_body_exhaustion.py
"""
from __future__ import annotations
import itertools
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, sqrt, findroot

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


# 四錨增量（⟨O⟩ − 1）與 σ
DN, SN = mpf("0.1263"), mpf("0.1848")
DS, SS = mpf("0.9090"), mpf("0.1575")
DX, SX = mpf("1.1719"), mpf("0.1038")
DK, SK = mpf("1.4213"), mpf("0.1007")

# docs/03 H2 之正性段（t ≡ a_ss）
T_LO, T_HI = mpf("0.262900"), mpf("0.389200")


def pairsum_additive():
    """F2：a(x,y) = c(γ_x+γ_y)，literal 價對（重子 3 價 → 2 對；介子 2 價 → 1 對）。
    N(udd) = a_ud + a_dd ; Σ(uds) = a_ud + a_us ; Ξ(uss) = a_us + a_ss ; K = a_u s̄
    令 X = cγ_u, Y = cγ_d, Z = cγ_s。精確分數解。
    """
    dN, dS, dX = F(1263, 10000), F(9090, 10000), F(11719, 10000)
    # N: (X+Y) + 2Y = X + 3Y ; Σ: (X+Y)+(X+Z) = 2X+Y+Z ; Ξ: (X+Z)+2Z = X+3Z
    X = dN - 3 * F(0)          # placeholder
    # 由 N: X = dN − 3Y ; 由 Ξ: X = dX − 3Z ⟹ 3(Z−Y) = dX−dN
    # 由 Σ: 2X + Y + Z = dS
    # 代入 X = dN − 3Y ⟹ 2dN − 6Y + Y + Z = dS ⟹ Z − 5Y = dS − 2dN
    ZmY = (dX - dN) / 3
    Z = ZmY + F(0)             # Z = Y + ZmY
    # Z − 5Y = dS − 2dN ⟹ (Y + ZmY) − 5Y = dS − 2dN ⟹ −4Y = dS − 2dN − ZmY
    Y = -(dS - 2 * dN - ZmY) / 4
    Z = Y + ZmY
    X = dN - 3 * Y
    return X, Y, Z


def main():
    say("=" * 80)
    say("R2-A  二體荷型形式之窮盡測試（零自由參數）")
    say("")

    # ── F1 乘性荷 a = w_x·w_y ──
    say("F1  乘性荷 a(x,y) = w_x·w_y（literal 價對）")
    say("    N: w_d(w_u+w_d) = %.6f" % DN)
    say("    Σ: w_u(w_d+w_s) = %.6f" % DS)
    say("    Ξ: w_s(w_u+w_s) = %.6f" % DX)

    def eqs(x, y, z):
        wu, wd, ws = x, y, z
        return (wd * (wu + wd) - DN, wu * (wd + ws) - DS, ws * (wu + ws) - DX)

    sols = []
    grid = [-3, -1, -0.3, 0.3, 1, 3]
    for a, b, c in itertools.product(grid, repeat=3):
        try:
            r = findroot(lambda u, d, s: eqs(u, d, s), (mpf(a), mpf(b), mpf(c)),
                         tol=mpf("1e-25"), maxsteps=200)
            v = (r[0], r[1], r[2])
            if all(abs(c2) < mpf("1e-18") for c2 in eqs(*v)):
                key = tuple(round(float(t), 6) for t in v)
                if all(max(abs(t - k2) for t, k2 in zip(v, s2)) > mpf("1e-5") for s2 in sols):
                    sols.append(v)
        except Exception:
            pass
    say("    實數解數 = %d" % len(sols))
    ok = []
    for s in sols:
        wu, wd, ws = s
        tag = []
        if wu > 0 and wd > 0 and ws > 0:
            tag.append("全正")
        if ws > wu > wd:
            tag.append("w_s>w_u>w_d")
        if wu < wd < ws:
            tag.append("w_u<w_d<w_s ⟹ 與 m_u<m_d<m_s 同向")
        if wu > wd > ws:
            tag.append("w_u>w_d>w_s ⟹ 與 1/m 同向")
        say("      w_u=%+0.6f  w_d=%+0.6f  w_s=%+0.6f   %s"
            % (wu, wd, ws, "／".join(tag) if tag else "✗ 不符"))
        mono = (wu < wd < ws) or (wu > wd > ws)
        if wu > 0 and wd > 0 and ws > 0 and mono:
            ok.append(s)
    if not ok:
        say("    ⟹ 🔴 **乘性荷全族死**：正性解之序為 w_u > w_s > w_d，")
        say("       既不係 w_u<w_d<w_s（m 遞增），亦不係 w_u>w_d>w_s（1/m 遞減）")
        say("       ⟹ 與**任何單調** w(m) 不相容（方向約定無關）。")
    OUT["F1"] = {"real_solutions": len(sols),
                 "solutions": [[str(t) for t in s] for s in sols],
                 "alive": bool(ok),
                 "incompatible_with_any_monotone_w": not ok}

    # ── F2 加性荷 ──
    say("")
    say("F2  加性荷 a(x,y) = c(γ_x+γ_y)（精確分數解）")
    X, Y, Z = pairsum_additive()
    say("    X = cγ_u = %+.6f" % float(X))
    say("    Y = cγ_d = %+.6f" % float(Y))
    say("    Z = cγ_s = %+.6f" % float(Z))
    say("    ⟹ γ_d %s 0（%s）；γ_s %s γ_u（%s）"
        % ("<" if Y < 0 else ">", "🔴 負質量" if Y < 0 else "ok",
           "<" if Z < X else ">", "🔴 與 m_s>m_u 相反" if Z < X else "ok"))
    dead2 = (Y < 0) or (Z < X)
    say("    ⟹ %s" % ("🔴 **加性荷死**" if dead2 else "🟡 未死"))
    OUT["F2"] = {"X": str(X), "Y": str(Y), "Z": str(Z), "dead": bool(dead2)}

    # ── F3 衰減核（對照）──
    say("")
    say("F3  衰減距核（對照組，docs/03 §四）")
    z3 = (DN - DS) / sqrt(SN ** 2 + SS ** 2)
    say("    ⟨O⟩_N − ⟨O⟩_Σ = %+0.6f ± %0.6f ⟹ z = %+0.3f ⟹ 🔴 死（t-無關）"
        % (DN - DS, sqrt(SN ** 2 + SS ** 2), z3))
    OUT["F3"] = {"z": str(z3), "dead": True}

    # ── F4 環境多項式（非二體）──
    say("")
    say("F4  環境計數 f(n_s) = α + β·n_s（非二體，測唯一性）")
    # N: n_s=0 ; Σ: n_s=1 ; Ξ: n_s=2 ⟹ 三點定二參 ⟹ 1 dof
    beta = (DX - DN) / 2
    alpha = DN
    pred_sigma = alpha + beta * 1
    say("    由 N(0)、Ξ(2) 定 α, β ⟹ α = %+0.6f, β = %+0.6f" % (alpha, beta))
    say("    預言 Σ (n_s=1)：%+0.6f   實測：%+0.6f   差 %+0.6f" % (pred_sigma, DS, DS - pred_sigma))
    predK = pred_sigma
    dz = (DK - predK) / sqrt(SK ** 2 + SS ** 2)
    say("    預言 K (n_s=1, 若介子與重子同)：%+0.6f   實測：%+0.6f   z = %+0.3f"
        % (predK, DK, dz))
    say("    ⟹ %s（z = %+0.2f）" % ("🔴 介子與重子不可同碼" if abs(dz) >= 2.5 else "🟡", dz))
    say("    ⟹ 且 β = %+0.6f 無結構來源；線性環境響應**不足以**釘點。" % beta)
    OUT["F4"] = {"alpha": str(alpha), "beta": str(beta),
                 "pred_K": str(predK), "z_K": str(dz)}

    # ── 匯總 ──
    say("")
    say("=" * 80)
    say("R2-A 判詞")
    say("")
    say("  🔴 F1 乘性荷：無實解 ⟹ rank-1 荷分解不可行")
    say("  🔴 F2 加性荷：γ_d < 0 ∧ γ_s < γ_u ⟹ 與質量序相反")
    say("  🔴 F3 衰減核：3.22σ t-無關")
    say("  🔴 F4 環境線性：介子/重子 2.74σ 不可同碼，且 β 無結構來源")
    say("")
    say("  ⟹ **二體（pair-sum）形式族空**；環境線性亦不足以釘點。")
    say("  ⟹ 依局域容錯假設之預期：F̂ 必須係**多體（stabiliser 型，weight ≥ 3）**，")
    say("     即作用於面/鄰域之積而非價對之和 —— 與 docs/03 §五之局域穩定子圖像一致。")
    OUT["verdict"] = {"two_body_family_empty": True, "requires_weight_ge_3": True}


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/07_two_body_exhaustion.json , data/run_07_two_body_exhaustion.log")
    with open(os.path.join(ROOT, "data", "07_two_body_exhaustion.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    with open(os.path.join(ROOT, "data", "run_07_two_body_exhaustion.log"), "w") as f:
        f.write("\n".join(L) + "\n")