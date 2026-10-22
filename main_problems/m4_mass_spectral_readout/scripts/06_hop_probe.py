#!/usr/bin/env python3
"""
/tmp 探測 — hop 階梯 vs 費米子質量譜 + 雙 probe 讀數差 (2026-10-21)
=====================================================================
方向: 「駐波態之間交互嘅同步邊界」— K_Chain 有 hop 階梯結構 (O179):
  hop = [1, 4, 8, 11, 12, 11, 8, 4, 1], Σ=60=|A5|, <hop>=4=B
  hop=1 → cos(π/4), cos(π/5)

未試: hop 階梯元素/差值 對應質量譜 (比 γ 比池更結構化)
未試: 雙 probe 讀數差 (π em vs K strong 反相)

測試:
  [A] hop 元素/差/積 vs 費米子質量比 (log 空間)
  [B] 質量差絕對值 vs hop 階梯 × 標度
  [C] 雙 probe 讀數差結構: em/strong 分裂比
"""
import numpy as np

m_e = 0.510998950e-3; m_mu = 0.1056583755; m_tau = 1.77686  # GeV
m_u = 0.00216; m_d = 0.00467; m_s = 0.0934; m_c = 1.27; m_b = 4.18; m_t = 172.69
m1 = 0.020908e-3; m2 = 0.022645e-3; m3 = 0.053586e-3  # eV→GeV

hop = np.array([1, 4, 8, 11, 12, 11, 8, 4, 1], dtype=float)
print("="*64)
print("探測: hop 階梯 vs 質量譜")
print("="*64)
print(f"hop = {hop.tolist()}, Σ={hop.sum()}, <hop>={hop.mean()}")

# [A] hop 元素對 (和/差/積/比) 生成結構常數
print("\n[A] hop 結構池 (和/差/積/比):")
hop_struct = []
for i in range(len(hop)):
    for j in range(i, len(hop)):
        if hop[i]+hop[j] > 0: hop_struct.append(hop[i]+hop[j])
        if hop[i]-hop[j] > 0: hop_struct.append(hop[i]-hop[j])
        if hop[i]*hop[j] > 0: hop_struct.append(hop[i]*hop[j])
        if hop[j] != 0: hop_struct.append(hop[i]/hop[j])
hop_struct = sorted(set(round(h, 10) for h in hop_struct))
print(f"  {len(hop_struct)} 個結構常數: {[round(h,4) for h in hop_struct][:30]}...")

# 質量比對 log10: 邊啲比 ~ hop 結構?
print("\n  質量比 vs hop 結構池:")
def check(name, t):
    best = None
    for h in hop_struct:
        if h <= 0: continue
        # t ~ h or t ~ 10^h or log10(t) ~ h?
        d1 = abs(t/h - 1)*100
        d2 = abs(np.log10(t)/h - 1)*100 if h != 0 else 1e9
        d3 = abs(np.log10(t) - h)*100  # log ratio ~ hop count (差):
        if best is None or min(d1, d3) < best[0]:
            best = (min(d1, d3), d1, d3, h)
    print(f"  {name}: ratio={t:.6f} log10={np.log10(t):.4f}  best hop={best[3]} | 直接比 dev {best[1]:.3f}% | log差 dev {best[2]:.4f}%")

check("τ/μ", m_tau/m_mu)
check("μ/e", m_mu/m_e)
check("t/b", m_t/m_b)
check("b/c", m_b/m_c)
check("c/s", m_c/m_s)
check("s/d", m_s/m_d)
check("m3/m2", m3/m2)
check("m3/m1", m3/m1)
check("τ/e", m_tau/m_e)
check("t/c", m_t/m_c)

# [B] hop 階梯 log 間距 vs 質量 log 間距
print("\n[B] hop 階梯 log 間距 vs 質量譜 log 間距:")
hoplog = np.diff(np.log10(hop[hop > 0]))
print(f"  hop>0 log 間距: {np.round(hoplog, 4).tolist()}")
for name, mlist in [("leptons(e,μ,τ)", [m_e, m_mu, m_tau]),
                    ("quarks(u,d,s,c,b,t)", [m_u, m_d, m_s, m_c, m_b, m_t]),
                    ("neutrinos(m1,m2,m3)", [m1, m2, m3])]:
    lo = np.diff(np.log10(sorted(mlist)))
    print(f"  {name} log 間距: {np.round(lo, 4).tolist()}")

# [C] 雙 probe 讀數差: em/strong 分裂
print("\n[C] 雙 probe 讀數差 (同對象, 唔同同步邊界):")
print(f"  π⁺−π⁰ = +4.5936 MeV  (em probe 讀數差)")
print(f"  K⁺−K⁰ = −3.937 MeV   (strong probe 讀數差)")
print(f"  n−p = +1.2933 MeV    (strong+em 混合)")
print(f"  ratio |π/K| = {4.5936/3.937:.5f}")
print(f"  ratio |π/(n−p)| = {4.5936/1.2933:.5f}")
print(f"  ratio |K/(n−p)| = {3.937/1.2933:.5f}")
r1, r2, r3 = 4.5936/3.937, 4.5936/1.2933, 3.937/1.2933
for r, n in [(r1, "|π/K|"), (r2, "|π/np|"), (r3, "|K/np|")]:
    best = None
    for h in hop_struct:
        if h <= 0: continue
        d = abs(r/h - 1)*100
        if best is None or d < best[0]:
            best = (d, h, r/h)
    print(f"  {n} = {r:.5f}: 最優 hop 比 {best[1]} (目標 {r/best[1]:.4f}) dev {best[0]:.3f}%")

# hop 對: 60/13.09? 直接 scan hop_i/hop_j:
print("\n  hop 對比 (i/j) 對質量差比:")
for i in range(len(hop)):
    for j in range(len(hop)):
        if hop[j] == 0: continue
        hr = hop[i]/hop[j]
        for r, n in [(r1, "π/K"), (r2, "π/np"), (r3, "K/np")]:
            d = abs(hr/r - 1)*100
            if d < 1:
                print(f"    hop[{i}]/hop[{j}] = {hr:.4f} vs {n}={r:.5f}: dev {d:.3f}%")