#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/14_symmetric_trace_emergence.py — ρ_u 與 ρ_d 之對稱演繹（u,d 無關公式）

調律者指示：「ρ_u 和 ρ_d 應該是同一套規則下的對稱性演繹；把 ρ_u 和 ρ_d 可以在和 u、d
完全無關的公式閉合或湧現」。

本腳本做三件事：
  ① 把 `(ρ_u, ρ_d)` 寫成 **A₂ Cartan 譜分解**下之**單一條對稱公式**（特徵值 1、3）；
  ② 由錨系統導出 **u,d 無關之精確恆等式** `ρ_u + ρ_d = 2Δ_N/Δ_Ξ`；
  ③ 檢查結構閉合 `Δ_N/Δ_Ξ = N_gen/B² = 3/16`（零參數）及 look-elsewhere。

錨系統（`docs/08 §〇`，`M12 docs/14 §三` 三錨唯一解、零自由度）：
  δ_s = Δ_Ξ/2          （Ξ⁻ = dss、Ξ⁰ = uss ⟹ 差 = δ_s ×2）
  δ_d = Δ_Σ − Δ_Ξ/2    （Σ⁻ = dds、Σ⁰ = uds）
  δ_u = Δ_N − δ_d      （n − p）
  Δ_N ≡ m_n − m_p、Δ_Σ ≡ m(Σ⁻) − m(Σ⁰)、Δ_Ξ ≡ m(Ξ⁻) − m(Ξ⁰)
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
PDG = ("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/main_problems/"
       "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt")
WANT = {"2112": "n", "2212": "p", "3112": "Σ⁻", "3212": "Σ⁰", "3312": "Ξ⁻", "3322": "Ξ⁰"}


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


print("=== 0. 輸入（PDG 2024，vendored；單位 MeV）===")
pd = parse(PDG, set(WANT))
for pid in ("2112", "2212", "3112", "3212", "3312", "3322"):
    m, s = pd[pid]
    print(f"  {WANT[pid]:>3s} ({pid}) = {mp.nstr(mev(m),10):>11s} ± {mp.nstr(mev(s),4)}")
mN  = mev(pd["2112"][0]) - mev(pd["2212"][0])
sN  = mp.sqrt(mev(pd["2112"][1])**2 + mev(pd["2212"][1])**2)
mS  = mev(pd["3112"][0]) - mev(pd["3212"][0])
sS  = mp.sqrt(mev(pd["3112"][1])**2 + mev(pd["3212"][1])**2)
mX  = mev(pd["3312"][0]) - mev(pd["3322"][0])
sX  = mp.sqrt(mev(pd["3312"][1])**2 + mev(pd["3322"][1])**2)
print(f"\n  Δ_N = m_n − m_p        = {mp.nstr(mN,10)} ± {mp.nstr(sN,4)}")
print(f"  Δ_Σ = m(Σ⁻) − m(Σ⁰)    = {mp.nstr(mS,10)} ± {mp.nstr(sS,4)}")
print(f"  Δ_Ξ = m(Ξ⁻) − m(Ξ⁰)    = {mp.nstr(mX,10)} ± {mp.nstr(sX,4)}")

print("\n=== 1. 三錨 ⟹ δ（零自由度）===")
dS = mX/2
dD = mS - mX/2
dU = mN - dD
print(f"  δ_s = Δ_Ξ/2            = {mp.nstr(dS,10)}")
print(f"  δ_d = Δ_Σ − Δ_Ξ/2      = {mp.nstr(dD,10)}")
print(f"  δ_u = Δ_N − δ_d        = {mp.nstr(dU,10)}")

print("\n=== 2. ρ_u、ρ_d 及其**對稱組合** ===")
rU, rD = dU/dS, dD/dS
print(f"  ρ_u = δ_u/δ_s = {mp.nstr(rU,12)}")
print(f"  ρ_d = δ_d/δ_s = {mp.nstr(rD,12)}")
sm, df = rU + rD, rD - rU
print(f"  ρ_u + ρ_d = {mp.nstr(sm,12)}")
print(f"  ρ_d − ρ_u = {mp.nstr(df,12)}")

print("\n=== 3. ★ 精確恆等式：和 = u,d 無關之單一比 ===")
print("  代數：ρ_u + ρ_d = (δ_u+δ_d)/δ_s = Δ_N/(Δ_Ξ/2) = **2Δ_N/Δ_Ξ**")
chk = 2*mN/mX
print(f"  驗算 2Δ_N/Δ_Ξ = {mp.nstr(chk,12)}  ⟹ 與 ρ_u+ρ_d 差 {mp.nstr(chk-sm,3)} （恆等）")
print("  而 ρ_d − ρ_u = (δ_d−δ_u)/δ_s = (2Δ_Σ−Δ_Ξ−Δ_N+Δ_Σ−Δ_Ξ/2)/(Δ_Ξ/2) …")
chk2 = (4*mS - 2*mN)/mX - 2
print(f"  驗算 = {mp.nstr(chk2,12)} ⟹ 與 ρ_d−ρ_u 差 {mp.nstr(chk2-df,3)} （恆等）")

print("\n=== 4. ★ 結構閉合：Δ_N/Δ_Ξ = N_gen/B² = 3/16 ===")
NGen, B = 3, 4
ratio = mN/mX
pred = mp.mpf(NGen)/B**2
sr = mp.sqrt((sN/mX)**2 + (mN*sX/mX**2)**2)
z = (ratio-pred)/sr
print(f"  觀測 Δ_N/Δ_Ξ = {mp.nstr(ratio,10)} ± {mp.nstr(sr,4)}")
print(f"  結構 N_gen/B² = {NGen}/{B**2} = {mp.nstr(pred,10)}")
print(f"  z = {mp.nstr(z,5)} σ  ⟹ **零參數閉合**")
print(f"  等價：ρ_u + ρ_d = 2·(3/16) = 3/8 = (B−1)/(2B)（因 N_gen = B−1 = 3）")

print("\n=== 5. A₂ Cartan 譜分解 ⟹ 單一對稱公式 ===")
print("  Cartan K = [[2,−1],[−1,2]]；特徵值 1（對稱方向）、3（反對稱方向）；tr = 2、det = 3")
s_, a_ = mp.mpf(1+2), mp.mpf(2-1)          # α=(1,2,4) ⟹ s=α_u+α_d=3、a=α_d−α_u=1
aS = mp.mpf(4)
print(f"  ρ_u = (s − 3a)/(4α_s) = {mp.nstr((s_-3*a_)/(4*aS),10)}")
print(f"  ρ_d = (s + 3a)/(4α_s) = {mp.nstr((s_+3*a_)/(4*aS),10)}")
print("  ⟹ 同一條公式、同一個 3（Cartan 特徵值）、只差 ± 號 = **對稱性演繹**")
print("  ⟹ 兩式相加 = 2s/(4α_s) = s/(2α_s) = 3/8（§4 之閉合）；相減 = 6a/(4α_s) = 3a/(2α_s)")

print("\n=== 6. 由結構反解 α（u,d 無關）===")
print("  ① A₂ 零方向：α_d = 2α_u（Cartan 元素 2，非參數）")
print("  ② 迹閉合：s/(2α_s) = (B−1)/(2B) ⟹ α_s/(α_u+α_d) = B/(B−1) = 4/3")
print("  ⟹ 令 α_u=1 ⟹ α = (1, 2, 4) **唯一**，無連續自由度")
print(f"  驗：α_s/(α_u+α_d) = 4/3；s/(2α_s) = 3/8；ρ_u = 0；ρ_d = 3/8")

print("\n=== 7. look-elsewhere（新觀測量 Δ_N/Δ_Ξ 之窗）===")
lo, hi = ratio-sr, ratio+sr
print(f"  1σ 窗 Δ_N/Δ_Ξ ∈ [{mp.nstr(lo,7)}, {mp.nstr(hi,7)}] （寬 {mp.nstr(100*(hi-lo)/ratio,3)}%）")
hits = []
for q in range(1, 21):
    for p in range(1, 21*q):
        r = F(p, q)
        if r.denominator != q:
            continue
        v = mp.mpf(p)/q
        if lo <= v <= hi:
            hits.append(r)
print(f"  分母 ≤ 20 之有理落窗：**{len(hits)} 個** → {', '.join(str(h) for h in hits[:8])}")
print(f"  其中 `N_gen/B² = 3/16` 係最簡結構式（其餘為未註冊之比）")

print("\n=== 8. 判 ===")
print("  ① 恆等式：ρ_u + ρ_d = 2Δ_N/Δ_Ξ（**完全 u,d 標籤無關**，純錨系統代數）")
print("  ② 閉合：Δ_N/Δ_Ξ = N_gen/B² = 3/16，z = +0.23σ（零參數）")
print("  ③ 對稱演繹：ρ_{u,d} = (s ∓ 3a)/(4α_s)，同一公式 ± 號")
print("  ④ 由 ①② 反解得 α = (1,2,4) 唯一 ⟹ ρ_d = 3/8 由結構湧現，非擬合")
print("  ⑤ 但：②③ 與 (ρ_u=0, ρ_d=3/8) 等價 ⟹ **資訊量不增**；增益在『形式』（u,d 無關）")
print("     與『可證偽對象』（單一比 Δ_N/Δ_Ξ），而 rival 15/4 對應 Δ_N/Δ_Ξ = 0.2（比率 5/25）")
sys.exit(0)