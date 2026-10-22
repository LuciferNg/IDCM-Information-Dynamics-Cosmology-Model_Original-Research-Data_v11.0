#!/usr/bin/env python3
"""
/tmp 探測 07 — s/d=20 EXACT 追擊 + hop⇄質量 log 間距體系 (2026-10-21)
=======================================================================
06 爆出: s/d = 20.000000 (exact 數字) 命中 hop 結構池 20 = hop[2]+hop[4] = 8+12
         = hop 元素和; 另 t/c≈132=11×12 (3.0%), 但 132 唔喺池。

追擊:
  1. s/d = 20 係咪「巧合數字」(PDG m_s=93.4±2.5, m_d=4.67±0.48 — 誤差 2.7%/10%)
     → 20.0 唔係高精度量! 要睇誤差傳播
  2. hop 階梯同微子 log 間距: [0.0347, 0.3741] vs hop [0.0378, 0.1383]?
     0.0347 ≈ 0.0378 (hop 3→4 間距, dev 8%)? 中微子 m1→m2 係「近簡併」(8%),
     而 hop 9 個元素 1,4,8,11,12,11,8,4,1 — 喺 11→12 附近 log 間距最細!
  3. 駐波態交互: m2−m1 (近簡併 8%) ~ hop 11→12 階 (微間距) —「同步邊界」稠密區?
  4. 質量譜 log 間距和: lepton [2.315,1.226] sum=3.54; quark [0.335,1.301,1.134,0.517,1.616] sum=4.90
"""
import numpy as np

def pdg_range(v, rel):
    return v*(1-rel), v*(1+rel)

print("="*64)
print("探測 07: s/d=20 誤差判別 + hop⇄質量體系")
print("="*64)

# ── 1. s/d 誤差 ──
m_s, m_d = 93.4, 4.67  # MeV
rel_s, rel_d = 2.5/93.4, 0.48/4.67  # PDG 2024 relative (approx)
lo, hi = pdg_range(m_s/m_d, rel_s+rel_d)
print("[1] s/d = 20.0 誤差判別:")
print(f"  m_s = 93.4 ± {2.5} MeV ({rel_s*100:.2f}%), m_d = 4.67 ± {0.48} MeV ({rel_d*100:.2f}%)")
print(f"  s/d = {m_s/m_d:.4f}, 誤差帶 [{lo:.3f}, {hi:.3f}] ({(rel_s+rel_d)*100:.2f}%)")
print(f"  → hop 20 命中喺誤差帶內 ⇒ 唔係高精度閉合, 係「巧合於 PDG 引用值」")
print(f"  → 需要更準 m_s/m_d 先講 (lattice QCD: m_s/m_d = 27.2 ± 0.6 @500MeV?)")

# 2. 微子 log 間距 vs hop log 間距
print("\n[2] 微子 m1→m2→m3 log 間距 vs hop log 間距:")
m1, m2, m3 = 0.020908, 0.022645, 0.053586  # meV? no eV
l12, l23 = np.log10(m2/m1), np.log10(m3/m2)
hop = np.array([1,4,8,11,12,11,8,4,1], dtype=float)
hlog = np.abs(np.diff(np.log10(hop[hop>0])))
print(f"  微子: log(m2/m1) = {l12:.4f}, log(m3/m2) = {l23:.4f}")
print(f"  hop 細間距: {sorted(hlog)[:4]} → {np.round(sorted(hlog)[:4],4)}")
print(f"  l12 vs hop 細間距: dev {abs(l12-sorted(hlog)[0])/sorted(hlog)[0]*100:.1f}% ({sorted(hlog)[0]:.4f})")
print(f"  l23 vs hop 細間距: dev {abs(l23-sorted(hlog)[1])/sorted(hlog)[1]*100:.1f}% ({sorted(hlog)[1]:.4f})")

# 3. hop 得 9 元素: 1,4,8,11,12,11,8,4,1 — 11→12 (log 0.0380) vs 微子 m1→m2 (log 0.0347)
print("\n[3] hop 11→12 微間距 = 微子 m1→m2 近簡併?")
print(f"  hop 11→12: log10(12/11) = {np.log10(12/11):.5f}")
print(f"  微子 m1→m2: log10(1.0831) = {np.log10(1.083078):.5f}")
print(f"  dev = {abs(np.log10(12/11)-np.log10(1.083078))/np.log10(12/11)*100:.2f}%")
print(f"  (64/60 = 1.0667? hop peak 12 vs Σ=60: {np.log10(60/56.4):.5f}?)")

# 4. 質量譜 log 間距和
print("\n[4] log 間距和 (譜總跨度):")
print(f"  lepton (e→τ): {np.log10(1.77686/0.000510999):.4f}")
print(f"  quark (u→t): {np.log10(172.69/0.00216):.4f}")
print(f"  neutrino (m1→m3): {np.log10(0.053586/0.020908):.4f}")
print(f"  hop (1→12): {np.log10(12/1):.4f}")
print(f"  lepton/neutrino span 比: {np.log10(1.77686/0.000510999)/np.log10(0.053586/0.020908):.3f}")
print(f"  quark/neutrino span 比: {np.log10(172.69/0.00216)/np.log10(0.053586/0.020908):.3f}")