#!/usr/bin/env python3
"""
/tmp hadron 質量差測試 01 — 數據表 + 事前掃描 (2026-10-21)
=============================================================
調律者 framing: 質量 = 駐波態之間交互嘅同步邊界譜讀數差值
⟹ hadron 質量差 (π⁺−π⁰, n−p, K⁺−K⁰, D, B, m_d−m_u) 先係名副其實嘅「譜差讀數」。

物理層:
  π⁺−π⁰ = +4.5936 MeV — 主要電磁 self-energy (photon probe 讀數差!)
  n−p    = 1.2933 MeV — d−u 質量差 + 電磁 (strong sync boundary)
  K⁺−K⁰ = −3.937 MeV  — 電磁 vs u−d 反相!
  m_d−m_u ≈ 2.3 MeV   — isospin breaking (uncertainty ~10%)

協議 (事前定義, 唔准見 target 反揀):
  [A] 質量差比 Δm_i/Δm_j vs γ比池 (同 R7 A 族)
  [B] 絕對質量差 Δm vs γ_n × A_scale (A ∈ {m_e, f_π, Λ, v/1000...})
  [C] 質量差組合 (Σ/差) 對結構池
"""
import numpy as np

N = 88
nl = np.arange(1, 45, dtype=float)
gamma = 4*np.sin(nl*np.pi/N)**2
inv_g = 1/gamma
rat = gamma[:, None]/gamma[None, :]
pool = rat.flatten()

# PDG 2024 / CODATA 精確值
dm = {
    "π⁺−π⁰":  4.5936,      # MeV ±0.0005 (PDG: 4.59354)
    "n−p":     1.29333205,  # MeV ±5e-7 (CODATA)
    "K⁺−K⁰": -3.937,       # MeV ±0.003
    "m_d−m_u": 2.3,         # MeV MS̄(2GeV) ±0.6 (~10% unc)
}
print("="*64)
print("Hadron 質量差差值閉合 01: 事前掃描")
print("="*64)
for k, v in dm.items():
    print(f"  {k} = {v} MeV (± 誤差見上)")

# [A] 質量差比 (取正數對)
print("\n[A] 質量差比 vs γ比池:")
pairs = [("π−π0 vs n−p", abs(4.5936/1.29333205)),
         ("π−π0 vs d−u", abs(4.5936/2.3)),
         ("n−p vs d−u", abs(1.29333205/2.3)),
         ("K vs π", abs(3.937/4.5936)),
         ("K vs n−p", abs(3.937/1.29333205)),
         ("K vs d−u", abs(3.937/2.3))]
for name, t in pairs:
    dev = np.abs(pool/t - 1)*100
    best = np.argmin(dev)
    i, j = np.unravel_index(best, rat.shape)
    print(f"  {name} = {t:>10.6f}: γ{i+1}/γ{j+1} = {rat[i,j]:.6f} dev {dev[best]:.4f}%")

# MC null: 隨機質量差對 (0.5..6 MeV log-均勻)
import random
rng = random.Random(25741)
hits = 0; NMC = 3000
for _ in range(NMC):
    t = 10**rng.uniform(np.log10(0.1), np.log10(10))
    if np.min(np.abs(pool/t - 1)*100) < 0.1:
        hits += 1
print(f"\n  MC null: 隨機質量差比 (0.1-10) γ比 <0.1% 命中率: {hits/NMC*100:.1f}%")

# [B] 絕對質量差 = γ × A_scale
print("\n[B] Δm = γ_n × A (A = 標度):")
scales = {"m_e×1": 0.511, "m_e×2": 1.022, "m_e×10": 5.11, "f_π×1/30": 4.333,
          "f_π(130MeV)/28": 4.64, "Λ_QCD/100": 3.4, "v/1e5": 2.46}
for name, t in [("π−π0", 4.5936), ("n−p", 1.2933), ("K", 3.937), ("d−u", 2.3)]:
    row = []
    for sname, A in scales.items():
        need_g = t/A
        if need_g <= 4:
            idx = np.argmin(np.abs(gamma - need_g))
            row.append(f"{sname}:γ{idx+1}={gamma[idx]:.3f} dev {abs(gamma[idx]/need_g-1)*100:.2f}%")
        else:
            row.append(f"{sname}:γ>4 ✗")
    print(f"  Δm({name}) = {t} MeV: " + " | ".join(row))

# [C] 結構組合: Δm 對 γ_a+γ_b, γ_a−γ_b, Σγ/n... ×m_e?
print("\n[C] Δm = m_e × 譜函數 (m_e = 0.511 MeV):")
for name, t in [("π−π0", 4.5936), ("n−p", 1.2933), ("K", 3.937), ("d−u", 2.3)]:
    r = t/0.511  # Δm/m_e
    dev = np.abs(pool/r - 1)*100
    best = np.argmin(dev)
    i, j = np.unravel_index(best, rat.shape)
    print(f"  Δm({name})/m_e = {r:.5f}: γ{i+1}/γ{j+1} = {rat[i,j]:.6f} dev {dev[best]:.4f}%")

# [D] 電磁signature: π−π0 (+4.59) vs K⁺−K⁰ (−3.94) — 反相!
print("\n[D] 電磁質量差 signature (物理結構, 非掃描):")
print("  π⁺−π⁰ > 0 (π⁺ 重): photon self-energy +")
print("  K⁺−K⁰ < 0 (K⁺ 輕): u 輕於 d ⟹ K⁺(us) 受 u 質量拉輕 overcomes em")
print("  ⟹ 兩個效應 (em + strong isospin) 反相 — 呢個本身就係「雙 probe 讀數差」!")
print("  組合: π−π0 = em  |  K−K0 = u−d 質量差 − em")
print(f"  解: u−d = K 差 + em?  (3.937 = u−d_eff − em)")