#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/16_trace_target_precision.py — 迹觀測量之裁決門檻（無 Σ 硬底）

由 `scripts/14/15`：迹組合 ρ_u+ρ_d = 2Δ_N/Δ_Ξ **完全不含 Δ_Σ** ⟹ `docs/08 §四` 之
Σ 硬底 (σ_ρd ≥ 0.01099) **對迹觀測量不適用**。
本腳本求：3σ 分辨 `3/8` vs 各競爭者所需之 `σ(Ξ⁰)`，並與 BESIII 可達精度對比。
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
PDG = ("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/main_problems/"
       "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt")
WANT = {"2112", "2212", "3312", "3322"}


def parse(path, wanted):
    out = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("*") or len(line) < 70:
            continue
        ids = [i.strip() for i in (line[0:8], line[8:16], line[16:24], line[24:32]) if i.strip()]
        for pid in wanted:
            if pid in ids:
                try:
                    m = F(line[33:51].strip().replace("E", "e"))
                    ep = F(line[52:60].strip().replace("E", "e"))
                    em = abs(F(line[61:69].strip().replace("E", "e")))
                except Exception:
                    continue
                out[pid] = (m, (ep+em)/2)
    return out


def mev(fr):
    return mp.mpf(fr.numerator)/fr.denominator*1000


pd = parse(PDG, WANT)
mN = mev(pd["2112"][0]) - mev(pd["2212"][0])
mX = mev(pd["3312"][0]) - mev(pd["3322"][0])
sXm = mev(pd["3312"][1])
sX0 = mev(pd["3322"][1])
sX = mp.sqrt(sXm**2 + sX0**2)
sSum = 2*mN/mX**2 * sX            # ∂/∂Δ_Ξ of 2Δ_N/Δ_Ξ

print("=== 0. 迹觀測量及其靈敏度 ===")
print(f"  Δ_N = {mp.nstr(mN,9)} MeV（σ ≈ 6e−7，可忽略）")
print(f"  Δ_Ξ = {mp.nstr(mX,9)} ± {mp.nstr(sX,5)} MeV")
print(f"  ∂(2Δ_N/Δ_Ξ)/∂Δ_Ξ = {mp.nstr(-2*mN/mX**2,6)} ⟹ σ_迹 = {mp.nstr(sSum,6)}")
print(f"  **Δ_Σ 完全唔入式** ⟹ `docs/08 §四` 之 Σ 硬底 0.01099 **不適用**")

print("\n=== 1. 3σ 分辨所需 σ(Ξ⁰) ===")
print(f"  {'對手':>8s} {'迹值':>10s} {'分離':>9s} {'現況σ':>8s} {'需σ_迹':>9s} {'需σ(Ξ⁰)':>9s} {'改善':>6s}")
for r_, lab in ((F(4), "3/8"), (F(15, 4), "15/4"), (F(11, 3), "11/3"), (F(7, 2), "7/2")):
    v = 3/(2*mp.mpf(r_.numerator)/r_.denominator)
    sep = mp.fabs(v-mp.mpf(3)/8)
    need_sum = sep/3
    # need_sum = (2Δ_N/Δ_Ξ²)·σ_ΔΞ  ⟹ σ_ΔΞ = need_sum·Δ_Ξ²/(2Δ_N)
    need_dX = need_sum*mX**2/(2*mN)
    need_X0 = mp.sqrt(max(need_dX**2 - sXm**2, mp.mpf(0)))
    impr = sX0/need_X0 if need_X0 > 0 else mp.inf
    lbl = "（基準）" if r_ == 4 else "🔴 不可達" if need_X0 < mp.mpf("0.013") else "✅ 可達"
    print(f"  {lab:>8s} {mp.nstr(v,7):>10s} {mp.nstr(sep,6):>9s} {mp.nstr(sep/sSum,5):>7s}σ "
          f"{mp.nstr(need_sum,7):>9s} {mp.nstr(need_X0,5):>9s} {mp.nstr(impr,4):>6s}×  {lbl}")

print("\n=== 2. 現況與 BESIII 可達性 ===")
print(f"  現 σ(Ξ⁰) = {mp.nstr(sX0,4)} MeV（NA48 1999，±0.06 stat ± 0.20 syst）")
print(f"  要 3σ 分辨 15/4（次佳對手）：σ(Ξ⁰) ≤ 0.134 ⟹ **只需 {mp.nstr(sX0/mp.mpf('0.134'),3)}×**")
print(f"  BESIII 系統標度 ~0.013 MeV ⟹ 遠低於門檻 ⟹ **實驗可裁決**")
print(f"  5σ 分辨 15/4：σ_迹 ≤ 0.005 ⟹ σ(Ξ⁰) ≤ 0.080（仍遠高於 0.013）✅")

print("\n=== 3. 反觀舊路線（ρ_d 單獨）===")
print("  舊：σ_ρd 硬底 0.01099（Δ_Σ 阻死）⟹ vs 15/4 需 0.00833 ⟹ **結構上不可能**")
print("  新：迹不含 Δ_Σ ⟹ 無硬底 ⟹ vs 15/4 只需 1.49× ⟹ **可裁決**")
print("  ⟹ 對稱（迹）觀測量唔止更精，仲**解除咗硬底**")

print("\n=== 4. 判 ===")
print("  ① `ρ_u + ρ_d = 3/8` 係零參數預言（由 α=(1,2,4) 推）；現況 z = +0.224σ")
print("  ② 迹上競爭者：15/4 −1.92σ、11/3 −2.69σ、7/2 −4.36σ")
print("  ③ 裁決門檻：σ(Ξ⁰) ≤ 0.134 MeV（1.5×，BESIII 可達）⟹ **本問題由 experiment-limited 變為 near-term decidable**")
print("  ④ 若 Ξ⁰ 更新後迹偏離 3/8 > 3σ ⟹ 候選死、α=(1,2,4) 死；反之升 🟢")
sys.exit(0)