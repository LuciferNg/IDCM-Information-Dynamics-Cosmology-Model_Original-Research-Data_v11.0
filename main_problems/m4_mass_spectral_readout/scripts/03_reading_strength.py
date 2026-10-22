#!/usr/bin/env python3
"""
/tmp 夸克差值閉合 03 — 讀數×譜強度形式測試 (m_q = γ_n²·A?)
========================================================================
R7 定案: m₃ = γ₆²·A (0.002%) — 讀數×譜強度係中微子成功形式。
夸克測試: m_q 有冇 m_q = γ²_{n_i}·A 形式?

核心問題:
  1. 6 個夸克質量 (跨越 10⁴ log) 對應 γ² ∈ [0,4] — 只可以 ×「唔同 A」
     → 若 A 逐代唔同就冇預言力 (fit 6 個 A = 0 自由度)
  2. 唯一有預言力: 質量差 m_q − m_q' = Δ(譜讀數)? 或質量列 = γ²譜 × 單一 A?
  3. 「同步邊界」語義測試: 用中微子 A=1.64044 eV 相同標度系, 夸克 m_t=172.69 GeV
     係咪 v²/A 或類似譜標度?

測試:
  [A] m_q = γ_n²·v (v=246 GeV 單一標度): 揀 n, 幾近?
  [B] m_q = γ_n·v  (線性): 6 夸克 → 揀 6 個 mode, 結果如何?
  [C] 中微子-夸克同譜: γ² 譜值 0.032666 (γ₆²) 對應 m₃=53.6 meV; 
      唔同粒子係咪「同一 γ² 值 + 唔同 A_sector」?
  [D] 量子數判準: 夸克分 two sectors (u,c,t / d,s,b) — sector 內比有冇譜身份?
"""
import numpy as np
N = 88
n = np.arange(1, 45, dtype=float)
gamma = 4*np.sin(n*np.pi/N)**2
g2 = gamma**2

v = 246.0
quarks = {"u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18, "t": 172.69}
print("="*64)
print("誇克 03: 讀數×譜強度形式")
print("="*64)

# [A] m = γ²·v
print("\n[A] m_q = γ_n²·v  (γ²∈[0,4], v=246GeV ⟹ m∈[0,984]GeV)")
print("   → 只有 t 喺範圍; u/d/s/c/b 全部 < 0.1 GeV 要 γ² tiny:")
for q, m in quarks.items():
    need = m/v  # γ² need
    if need > 4:
        print(f"  m_{q} = {m:8.4g} GeV: 需要 γ²={need:.4g} > 4 ✗ (唔可能)")
    else:
        # γ²<4 最低 γ1² = (4sin(π/88)²)² ≈ 0.00033
        # 所有輕夸克需要 γ² < γ1² 嘅話唔可能
        idx = np.argmin(np.abs(g2 - need))
        print(f"  m_{q} = {m:8.4g} GeV: 需要 γ²={need:.5g} → 最接近 γ{idx+1}²={g2[idx]:.5g} dev {abs(g2[idx]/need-1)*100:.2f}%")

# [B] m = γ·v (線性)
print("\n[B] m_q = γ_n·v:")
for q, m in quarks.items():
    need = m/v
    idx = np.argmin(np.abs(gamma - need))
    print(f"  m_{q}: 需要 γ={need:.5g} → γ{idx+1}={gamma[idx]:.5g} dev {abs(gamma[idx]/need-1)*100:.2f}%")

# [C] sector 比 (u,c,t / d,s,b)
print("\n[C] sector 內比 (up-type / down-type):")
print("  up: u=2.16MeV, c=1.27GeV, t=172.69GeV")
print("  dn: d=4.67MeV, s=93.4MeV, b=4.18GeV")
print(f"  t/c = {172.69/1.27:.5f}   c/u = {1.27/0.00216:.5f}")
print(f"  b/s = {4.18/0.0934:.5f}   s/d = {0.0934/0.00467:.5f}")
print(f"  t/b = {172.69/4.18:.5f}   c/s = {1.27/0.0934:.5f}")

# [D] 質量比 vs γ比 — 已 01 [B] 判: 30 比 14 個 <1%, 但 MC null 97.9% 撞池
# 呢度試「短階梯」: sector 內比都有 γ 身份?
print("\n[D] sector 內比 γ比池命中 (look-elsewhere):")
rat = gamma[:, None]/gamma[None, :]
for name, t in [("t/c", 172.69/1.27), ("c/u", 1.27/0.00216), ("b/s", 4.18/0.0934), ("s/d", 0.0934/0.00467)]:
    dev = np.abs(rat/t - 1)*100
    best = np.unravel_index(np.argmin(dev), dev.shape)
    print(f"  {name} = {t:10.4f}: γ{best[0]+1}/γ{best[1]+1} = {rat[best]:.4f} dev {dev[best]:.4f}%")

# [E] 中微子關聯: m₃ = γ₆²A = 0.053586 eV, A = 1.64044 eV
# 夸克: m_t = 172.69e9 eV。γ² 有邊個值 × 某標度?
print("\n[E] 中微子 m₃ = γ₆²·A (0.002%) — 夸克同構需要 A_q:")
# A_q = m_q/γ² — 對邊個 mode? 若 A 統一 = v²/M_comp·1e9 類似嘅「sector 標度」
# 反測試: m_q/γ² (揀最好 mode) 係咪 constant?
print("  6 夸克 m_q/γ² (每個 quark 揀 mode 令值最近某常數) — 冇預言力, 跳過")
print("  唯一可測: m_t/m_b vs γ-a/γ-b 得 1.1% (01) — 唔過")
print("\n結論: 誇克冇出現 R7 中微子嗰種「讀數×譜強度」單 hit 0.002% 級結構")
print("      γ比池 density 97.9% ⟹ 任何 0.04-1% 命中都係撞池後驗")