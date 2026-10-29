#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2-E —— 目標重述：純質量分割比（σ 比 ⟨O⟩ 銳 8–10×）
2026-10-28

理由：⟨O⟩ = (Δm − em)/X 繼承 em 分解之 σ（16%/11%/8%/4%）。
      但 IDCM 從不分解 em/strong（M11 已判無通道分辨）⟹ 可測對象係**總分割之比值**，
      其 σ 由 PDG 質量精度支配（0.5%/0.8%/3.1%）⟹ 銳 8–10×。

本腳本：(A) 比值表與 σ；(B) 連分數收斂子掃描；(C) 密度計數；(D) 判詞。
不出任何預言。

跑法：/usr/bin/python3 scripts/11_splitting_ratios.py
"""
from __future__ import annotations
import json
import os
from fractions import Fraction

from mpmath import mp, mpf, sqrt, fabs

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


# PDG 2024 質量與 σ（來源：v11 M11 `scripts/01_r1_criteria.py` 表；K 之 σ 見 m7 `scripts/14`）
M = {
    "n": ("939.565421", "0.000000"),
    "p": ("938.272089", "0.000000"),
    "K+": ("493.677", "0.016"), "K0": ("497.611", "0.013"),
    "S-": ("1197.449", "0.030"), "S0": ("1192.642", "0.024"), "S+": ("1189.37", "0.07"),
    "X-": ("1321.71", "0.07"), "X0": ("1314.86", "0.20"),
    "D+": ("1869.66", "0.50"), "D0": ("1864.84", "0.50"),
}
MV = {k: mpf(v[0]) for k, v in M.items()}
MS = {k: mpf(v[1]) for k, v in M.items()}

CH = [
    ("n-p",   ("n", "p"),  +1),
    ("K+-K0", ("K+", "K0"), -1),   # 報告用 K+−K0（負值）
    ("S--S0", ("S-", "S0"), +1),
    ("X--X0", ("X-", "X0"), +1),
    ("S--S+", ("S-", "S+"), +1),
    ("D+-D0", ("D+", "D0"), +1),
]


def split(ch):
    a, b = ch[1]
    d = MV[a] - MV[b]
    s = sqrt(MS[a] ** 2 + MS[b] ** 2)
    return d, s


def cf_convergents(x, qmax=120):
    """回傳 (p, q) 收斂子列表，q <= qmax（標準遞推）。"""
    out = []
    a0 = x
    h1, h2 = 1, 0
    k1, k2 = 0, 1
    v = a0
    for _ in range(40):
        ai = int(mp.floor(v))
        h = ai * h1 + h2
        kk = ai * k1 + k2
        if kk > qmax:
            break
        out.append((h, kk))
        h2, k2 = h1, k1
        h1, k1 = h, kk
        fr = v - ai
        if fabs(fr) < mpf("1e-25"):
            break
        v = 1 / fr
    return out


def main():
    say("=" * 88)
    say("R2-E1  純分割比（目標重述）")
    say("")
    base, sbase = split(CH[0])
    say("  基準：Δm(n−p) = %.6f ± %.6f MeV" % (float(base), float(sbase)))
    say("")
    say("  %-8s %12s %10s %10s %10s %12s" % ("chan", "Δm (MeV)", "σ", "相對σ", "r = Δm/Δm_np", "σ_r"))
    rows = {}
    for name, ab, sg in CH:
        d, s = split((name, ab, sg))
        r = d / base
        sr = s / base
        rows[name] = (d, s, r, sr)
        say("  %-8s %12.6f %10.6f %9.3f%% %12.6f %12.6f  (%.2f%%)"
            % (name, float(d), float(s), float(s / fabs(d) * 100), float(r), float(sr),
               float(sr / fabs(r) * 100)))
    say("")
    say("  對照：⟨O⟩ 之相對 σ = 16.4% / 11.4% / 8.2% / 4.2%（docs/06 §二）")
    say("  ⟹ 純分割比之相對 σ = 0.00% / 0.52% / 0.80% / 3.09% ⟹ **銳 8–10×**（Ξ 只 1.4×）")
    OUT["E1"] = {k: {"dm": str(v[0]), "sig": str(v[1]), "r": str(v[2]), "sig_r": str(v[3])}
                 for k, v in rows.items()}

    say("")
    say("=" * 88)
    say("R2-E2  連分數收斂子掃描（q ≤ 120）—— 唯一可零參數検之形式")
    say("")
    ntest = 0
    hits = []
    for name in ["S--S0", "X--X0", "K+-K0", "S--S+"]:
        d, s, r, sr = rows[name]
        say("  【%s】r = %.6f ± %.6f" % (name, float(r), float(sr)))
        cands = cf_convergents(fabs(r), 120)
        ntest += len(cands)
        shown = 0
        for p, q in cands:
            val = mpf(p) / q
            z = (val - fabs(r)) / sr
            mark = ""
            if fabs(z) < 2:
                mark = "  ← |z|<2"
                hits.append((name, p, q, float(val), float(z)))
            if q >= 2 and (fabs(z) < 2 or shown < 3):
                say("     %4d/%-4d = %.8f   z = %+7.2f%s" % (p, q, float(val), float(z), mark))
                shown += 1
        say("     收斂子數 = %d" % len(cands))
    say("")
    say("  總檢定數 N_test = %d" % ntest)
    say("  |z|<2 之命中數 = %d" % len(hits))
    if hits:
        for h in hits:
            say("     %s：%d/%d = %.8f  z=%+.2f" % h)
    say("")
    # look-elsewhere：|z|<2 之單次機率 ≈ 0.9545 ⟹ 期望命中 = 0.9545*N
    exp_hits = mpf("0.9545") * ntest
    say("  Look-elsewhere：|z|<2 之單次機率 ≈ 0.9545 ⟹ 期望命中 = %.1f（N_test=%d）"
        % (float(exp_hits), ntest))
    say("  ⟹ 命中數 %d vs 期望 %.1f ⟹ 無訊號" % (len(hits), float(exp_hits)))
    say("")
    say("  🔴 **更重要：連分數掃描係恆真（tautology）** —— 任何實數之收斂子皆以 1/q² 收斂，")
    say("     故『某收斂子 |z|<0.1』必然成立（本例 197/53 得 z=+0.01、223/60 得 z=−0.00）。")
    say("     ⟹ 連分數命中**不可**作結構證據（同 M9 φ-power 教訓）。")
    say("     ⟹ 只可檢**事前指定**之代數恆等式（如 88/29、55/… 之類由結構給出者）。")
    OUT["E2"] = {"n_test": ntest, "n_hits": len(hits),
                 "expected_hits": str(exp_hits), "hits": hits, "cf_tautology": True}

    say("")
    say("=" * 88)
    say("R2-E3  密度對照（同一符號庫，窗口寬度比）")
    say("")
    say("  docs/05 G7：窗口寬 0.1263 內有 48 個庫成員 ⟹ 密度 380/單位")
    say("  本輪：最銳之 r(K) σ_r = %.6f ⟹ 窗口 2σ 寬 %.6f" % (float(rows["K+-K0"][3]), float(4 * rows["K+-K0"][3])))
    say("        以同一密度 380/單位外推 ⟹ 期望庫成員數 ≈ %.2f"
        % float(4 * rows["K+-K0"][3] * 380))
    say("  ⟹ 絕對窗口只收窄 ~%.1f× ⟹ **密度閘未回復**（仍 ~24 個成員在窗內）"
        % float(mpf("0.1263") / (4 * rows["K+-K0"][3])))
    say("  ⟹ 真正之增益係**相對精度**（0.52% vs 4.2%）：令『已推導候選』可被 0.5% 判死/判生，")
    say("     但**唔**令『掃描擬合』可行。窗口法繼續作廢。")
    OUT["E3"] = {"density_per_unit": 380,
                 "expected_members_in_2sigma_K": str(4 * rows["K+-K0"][3] * 380)}

    say("")
    say("=" * 88)
    say("R2-E 判詞")
    say("")
    say("  ✅ 目標重述：IDCM 唔分解 em/strong（M11 已判無通道分辨）⟹ 可測物 = 總分割比值")
    say("     相對 σ：K 0.52%｜Σ 0.80%｜Ξ 3.09%｜Σ⁻−Σ⁺ 0.94%｜n−p 0 ⟹ 比 ⟨O⟩ 銳 8–10×")
    say("  🔴 連分數掃描恆真 ⟹ 不可作證據（永久禁用為命中產生器）")
    say("  🔴 絕對窗口只收窄 ~2× ⟹ 密度閘仍失效（~24 成員）⟹ 窗口法不入 R3")
    say("  ✅ 唯一合法用法：**事前指定之代數式** → 以 0.5% 精度判死/判生")
    say("  ⟹ M12 之可交付物 = 尖銳目標表 ＋ 可辨識性地圖 ＋ 數據請求（非預言）。")


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/11_splitting_ratios.json , data/run_11_splitting_ratios.log")
    with open(os.path.join(ROOT, "data", "11_splitting_ratios.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    with open(os.path.join(ROOT, "data", "run_11_splitting_ratios.log"), "w") as f:
        f.write("\n".join(L) + "\n")