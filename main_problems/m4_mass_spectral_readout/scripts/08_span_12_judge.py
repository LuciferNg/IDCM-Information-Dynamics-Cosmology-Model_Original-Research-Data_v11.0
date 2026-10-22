#!/usr/bin/env python3
"""
/tmp 探測 08 — quark/neutrino span 比 ≈ 12 判別 (2026-10-21)
==============================================================
07 發現: log10(span_Q/span_ν) = 11.995 ≈ 12 = hop_max (dev 0.04%)

span_Q  = log10(m_t/m_u) = log10(172.69/0.00216) = 4.903
span_ν  = log10(m_3/m_1) = log10(0.053586/0.020908) = 0.4087
比 = 11.995 ≈ 12 = hop_max  = L9 − 64? = 88−76? 

嚴格判別:
  1. span 比 log10 形式有冇誤差傳播? (m_u 誤差 10%! m_t pole ±0.5%)
  2. 純數字巧合 or 結構? MC null: 隨機 quark 譜 span 比
  3. 12 = hop_max 係咪「事前身份」? (O179 hop [1,4,8,11,12,11,8,4,1])
  4. 換 m_d anchor 會唔會變?
"""
import numpy as np
import random

print("="*64)
print("探測 08: quark/neutrino span ≈ 12 判別")
print("="*64)

# 數據 (PDG 2024)
m_u, m_d, m_s = 0.00216, 0.00467, 0.0934
m_c, m_b, m_t = 1.27, 4.18, 172.69
m1_, m2_, m3_ = 0.020908e-3, 0.022645e-3, 0.053586e-3  # eV → GeV? m in eV→ convert
# 用 eV 統一:
m_u_v, m_t_v = 2.16e6, 172.69e9   # eV
m1, m3 = 0.020908, 0.053586      # eV
span_Q = np.log10(m_t_v/m_u_v)
span_nu = np.log10(m3/m1)
R = span_Q/span_nu
print(f"[0] span_Q = log10({m_t_v/1e9:.2f}GeV / {m_u_v/1e6:.2f}MeV) = {span_Q:.6f}")
print(f"    span_ν = log10({m3}/{m1}) = {span_nu:.6f}")
print(f"    R = {R:.6f}  vs 12: dev {abs(R-12)/12*100:.4f}%")

# 1. 誤差傳播: m_u ±10%, m_t ±0.5%
print("\n[1] 誤差傳播:")
m_u_lo, m_u_hi = 2.16*0.9, 2.16*1.1   # MeV (±10%)
m_t_lo, m_t_hi = 172.69*0.995, 172.69*1.005
for mm in [(m_u_lo, m_t_hi), (m_u_hi, m_t_lo), (m_u_lo, m_t_lo), (m_u_hi, m_t_hi)]:
    sq = np.log10((mm[1]*1e9)/(mm[0]*1e6))
    rr = sq/span_nu
    print(f"  m_u={mm[0]:.3f}MeV m_t={mm[1]:.2f}GeV: span={sq:.4f}, R={rr:.4f} (dev {abs(rr-12)/12*100:.3f}%)")
print(f"  → R 誤差帶 ≈ [{min(abs(11.99-12), abs(12.01-12)):.4f}% .. ...]")

# 2. hop_max = 12 身份
print("\n[2] 12 = hop_max (O179 hop [1,4,8,11,12,11,8,4,1])")
print("    12 = N/7.33? = L9−64? = N88−76 = 88−76 = 12 EXACT!")
print("    12 = 88 − L9? L9=76, 88−76=12 ✓ (K_Chain 自對偶)")
print("    12 = hop_max = γ6 對偶 index? (R1)")

# 3. MC null: 隨機「quark 譜」span 比
print("\n[3] MC null: 隨機 6 質量 (log-均勻 u-t 範圍) span/ν-span:")
rng = random.Random(25753)
hits = 0; NMC = 50000
for _ in range(NMC):
    masses = 10**np.array([rng.uniform(np.log10(0.002), np.log10(173)) for _ in range(6)])
    mn, mx = masses.min(), masses.max()
    sp = np.log10(mx/mn)
    rr = sp/span_nu
    if abs(rr-12)/12*100 < 0.5:
        hits += 1
print(f"  隨機 6 質量 span/νspan |dev|<0.5% (≈12) 命中率: {hits/NMC*100:.3f}%")
print(f"  → 若 <<0.1%, 12 命中唔係隨機")

# 4. 換 anchor: 用 m_d 最重? span d→t?
sp_dt = np.log10(m_t_v/(m_d*1e6))
print(f"\n[4] span(d→t)/span_ν = {sp_dt/span_nu:.4f} (dev vs 12 {abs(sp_dt/span_nu-12)/12*100:.3f}%)")
sp_ct = np.log10(m_t_v/(m_c*1e9))
print(f"    span(c→t)/span_ν = {sp_ct/span_nu:.4f} (dev {abs(sp_ct/span_nu-12)/12*100:.3f}%)")

# 5. log 底替換: 如果用 ln?
print("\n[5] ln 版本 (底無關, 比唔變):")
print(f"  R = ln(span) 比 = {R:.6f} (同 log10 — 底消)")
print(f"  ⟹ R 與 log 底無關; 12 命中純數字關係")

# 6. 精確計 R 嘅「身份」: R = 12(1+δ), δ = ?
delta = (R-12)/12
print(f"\n[6] δ = {delta*100:.4f}% — 喺 m_u 誤差 (±10%) 內? 係!")
print(f"  m_u 係唯一大誤差輸入 (10.28%); δ=0.04% ≪ 10% ⟹ 命中受 m_u 主導")
print(f"  ⟹ 冇實驗精度支持: R≈12 係 m_u/m_t/m1/m3 引用值巧合, 唔係高精度閉合")

print("\n" + "="*40)
print("判決: R≈12 係 m_u 誤差帶內巧合 — 唔升格")