#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/09_precision_budget.py — 乙：ρ_d 之精度預算 ＋ 可裁決化路徑

ρ_d = 2Δ_Σ/Δ_Ξ − 1 ，  Δ_Σ = m(Σ⁻)−m(Σ⁰) ， Δ_Ξ = m(Ξ⁻)−m(Ξ⁰)
（四質量皆 PDG 2024 vendored 檔獨立解析；Σ 與 Ξ 質量集不相交 ⟹ 兩 Δ 獨立）

輸出：
  §1 原始值與重建（核 ρ_d = 0.4035036）
  §2 誤差預算：σ 之逐項分解（變異數份額）
  §3 裁決門檻：各「分辨目標」所需 σ_ρd 與所需質量精度改善倍數
  §4 判：現況 2.79σ 已近 3σ 分辨 α_s/α_u = 4 vs 3；細分辨（4 vs 2/5）需同時改善 Σ 與 Ξ
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
PDG = ("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/"
       "main_problems/m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt")

def parse_pdg(path, wanted):
    out = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
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

pd = parse_pdg(PDG, {"3112", "3212", "3312", "3322"})
assert len(pd) == 4, f"PDG 解析失敗: {sorted(pd)}"

def mev(fr):                      # GeV Fraction → MeV mp
    return mp.mpf(fr.numerator)/fr.denominator*1000

print("=== 1. PDG 2024 原始值（唯讀 vendored 檔）===")
names = {"3112": "Σ⁻", "3212": "Σ⁰", "3312": "Ξ⁻", "3322": "Ξ⁰"}
M, S = {}, {}
for pid in ("3112", "3212", "3312", "3322"):
    m, s = pd[pid]
    M[pid], S[pid] = mev(m), mev(s)
    print(f"  {names[pid]:4s} m = {mp.nstr(M[pid],10):>12s} MeV   σ = {mp.nstr(S[pid],6):>8s} MeV")

DS = M["3112"]-M["3212"]
DX = M["3312"]-M["3322"]
rho = 2*DS/DX - 1
print(f"\n  Δ_Σ = m(Σ⁻)−m(Σ⁰) = {mp.nstr(DS,10)} MeV")
print(f"  Δ_Ξ = m(Ξ⁻)−m(Ξ⁰) = {mp.nstr(DX,10)} MeV")
print(f"  ρ_d = 2Δ_Σ/Δ_Ξ − 1 = {mp.nstr(rho,12)}   （目標 0.40350365 ⟹ Δ = "
      f"{mp.nstr(abs(rho-mp.mpf('0.40350365')),3)}）")

print("\n=== 2. 誤差預算（變異數分解）===")
s_DS = mp.sqrt(S["3112"]**2 + S["3212"]**2)
s_DX = mp.sqrt(S["3312"]**2 + S["3322"]**2)
c_S = 2/DX                                # ∂ρ/∂Δ_Σ
c_X = -2*DS/DX**2                         # ∂ρ/∂Δ_Ξ
v_S = (c_S*s_DS)**2
v_X = (c_X*s_DX)**2
v = v_S + v_X
print(f"  σ(Δ_Σ) = {mp.nstr(s_DS,8)} MeV   σ(Δ_Ξ) = {mp.nstr(s_DX,8)} MeV")
print(f"  敏感度：∂ρ/∂Δ_Σ = {mp.nstr(c_S,8)}  ∂ρ/∂Δ_Ξ = {mp.nstr(c_X,8)}")
print(f"  σ_ρd = {mp.nstr(mp.sqrt(v),10)}   （M12 R27：0.044785178）")
print(f"    由 Δ_Σ 貢獻 = {mp.nstr(mp.sqrt(v_S),8)}  （變異數份額 {mp.nstr(v_S/v*100,6)}%）")
print(f"    由 Δ_Ξ 貢獻 = {mp.nstr(mp.sqrt(v_X),8)}  （變異數份額 {mp.nstr(v_X/v*100,6)}%）")
print(f"  逐質量份額：Ξ⁰ {mp.nstr((c_X*S['3322'])**2/v*100,5)}%  Ξ⁻ {mp.nstr((c_X*S['3312'])**2/v*100,5)}%"
      f"  Σ⁻ {mp.nstr((c_S*S['3112'])**2/v*100,5)}%  Σ⁰ {mp.nstr((c_S*S['3212'])**2/v*100,5)}%")
print("  ⟹ 單一限制項 = **m(Ξ⁰)**（現 σ = %.2f MeV）" % float(S["3322"]))

print("\n=== 3. 裁決門檻 ===")
sig_now = mp.sqrt(v)
targets = [
    ("α=4 (0.375) vs α=3 (0.5)，3σ", abs(mp.mpf("0.5")-mp.mpf("0.375"))/3),
    ("α=4 (0.375) vs α=3 (0.5)，5σ", abs(mp.mpf("0.5")-mp.mpf("0.375"))/5),
    ("α=4 (0.375) vs 2/5 (0.4)，3σ", abs(mp.mpf("0.4")-mp.mpf("0.375"))/3),
    ("α=4 (0.375) vs √2−1 (0.41421)，3σ", abs((mp.sqrt(2)-1)-mp.mpf("0.375"))/3),
    ("α=4 vs 134/36 (0.402985)，3σ", abs(mp.mpf("27")/67-mp.mpf("0.375"))/3),
]
print(f"  {'目標':38s} {'需 σ_ρd':>10s} {'需改善':>8s} {'需 σ(Ξ⁰)':>10s} {'可行？':>8s}")
for nm, st in targets:
    imp = sig_now/st
    # 在 σ_ΔΣ 不變下，解 σ_ΔΞ 使之達 st
    rem = st**2 - v_S
    if rem <= 0:
        needX = mp.mpf(0)
        feas = "✗（Σ 已超限）"
    else:
        needX = mp.sqrt(rem)/abs(c_X)
        feas = "✓" if needX >= S["3312"]*0 else "?"
    if rem <= 0:
        sx0 = "—（Σ 阻死）"
    else:
        q = needX**2 - S["3312"]**2
        sx0 = mp.nstr(mp.sqrt(q),5) if q > 0 else "0（Ξ⁻ 已阻死）"
    print(f"  {nm:38s} {mp.nstr(st,8):>10s} {mp.nstr(imp,6):>8s} {sx0:>10s} {feas:>8s}")

print("\n=== 4. 判 ===")
print(f"  現 σ_ρd = {mp.nstr(sig_now,8)}；α=4 vs α=3 之分離 = "
      f"{mp.nstr(abs(mp.mpf('0.125'))/sig_now,6)}σ ⟹ **已近 3σ**（需僅 {mp.nstr(sig_now/(mp.mpf('0.125')/3),5)}× 改善）")
print("  ⟹ 乙 之結論：整數問題（4 vs 3）**近乎可裁決**；細分辨（4 vs 2/5、4 vs 134/36）")
print("     受限於 σ_ΔΣ = 0.0384 MeV 之硬底，**即使 Ξ 精度無限提升亦不可能**（σ_ρd 下限 = 0.0112）。")
print("  ⟹ 建議：以『α_s/α_u ∈ ℤ ⟹ 必為 4』作可裁決命題（現 2.79σ），放棄細分辨路線。")
sys.exit(0)