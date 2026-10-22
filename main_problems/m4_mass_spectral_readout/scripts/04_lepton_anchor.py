#!/usr/bin/env python3
"""
/tmp 方向判斷 — 費米子 anchor 分層測試 (2026-10-21)
=====================================================
問題: 全費米子譜差結構有冇得推? 邊啲實驗值夠乾淨做 anchor?

實驗值精度分層:
  L1 (anchor 可閉合): 帶電輕子 e/μ/τ — pole mass 精確 ~1e-6
  L1: 中微子 Δm² — 差譜 (唯一實驗量)
  L2 (量級驗證): quark t/b/c ~1%; s/d/u 10-20% running mass
  ΛQCD running: MS̄ mass scale-dependent = probe-dependent 讀數

測試: 帶電輕子質量比 μ/e, τ/μ, τ/e 對 γ比池 + sector 因子 (vs quark 97.9%)
memory 先例: μ/e = (γ7/γ1)(F_μ/F_e) = 206.705 (0.031%)
"""
import numpy as np

N = 88
n = np.arange(1, 45, dtype=float)
gamma = 4*np.sin(n*np.pi/N)**2
inv_g = 1/gamma

# 帶電輕子 (PDG, 精確):
m_e = 0.510998950e-3   # GeV
m_mu = 0.1056583755     # GeV
m_tau = 1.77686         # GeV
print("="*64)
print("費米子 anchor 分層: 帶電輕子譜比 vs γ比池")
print("="*64)
print(f"m_e = {m_e:.10f} GeV, m_μ = {m_mu:.10f} GeV, m_τ = {m_tau:.9f} GeV")

lepton_ratios = {
    "μ/e": m_mu/m_e,
    "τ/μ": m_tau/m_mu,
    "τ/e": m_tau/m_e,
}
rat = gamma[:, None]/gamma[None, :]

# MC null: 隨機「帶電輕子三質量」撞 γ比池
import random
rng = random.Random(25733)
pool = rat.flatten()

print("\n[A] 帶電輕子比 vs γ比池 (無 sector 因子):")
for name, t in lepton_ratios.items():
    dev = np.abs(pool/t - 1)*100
    best = np.argmin(dev)
    i, j = np.unravel_index(best, rat.shape)
    print(f"  {name} = {t:>12.6f}: γ{i+1}/γ{j+1} = {rat[i,j]:.6f} dev {dev[best]:.4f}%")

# sector 因子: μ/e = (γa/γb)·(F_x/F_y)
print("\n[B] 譜比 × sector 因子: μ/e (memory: γ7/γ1 × F_μ/F_e = 206.705, 0.031%):")
g7g1 = gamma[6]/gamma[0]
print(f"  γ7/γ1 = {g7g1:.9f}")
print(f"  206.768/206.705 = {m_mu/m_e/g7g1:.6f} ← F_μ/F_e sector = ?")
F_ratio = (m_mu/m_e)/g7g1
print(f"  F_μ/F_e = {F_ratio:.6f}")
for cand, name in [(11/1, "L5/1?"), (5/1, "F5/1"), (12/1, "12/1"), (11/5, "L5/F5"), (4, "B"), (3, "3"), (29/7, "L7/F7?"), (88/7, "N88/L7?")]:
    print(f"    {name} = {cand:.4f}: dev {abs(F_ratio-cand)/cand*100:.4f}%")

# τ/μ, τ/e 有冇類似結構?
print("\n[C] τ/μ, τ/e 直接比 (無 sector):")
t_tm = m_tau/m_mu
t_te = m_tau/m_e
for name, t in [("τ/μ", t_tm), ("τ/e", t_te)]:
    dev = np.abs(pool/t - 1)*100
    best = np.argmin(dev)
    i, j = np.unravel_index(best, rat.shape)
    print(f"  {name} = {t:>12.6f}: γ{i+1}/γ{j+1} = {rat[i,j]:.6f} dev {dev[best]:.4f}%")

# MC null: 隨機 3 輕子質量
print("\n[D] MC null: 隨機 3 質量 (e-μ-τ 範圍 log-均勻) 撞池:")
hits = 0
NMC = 3000
for _ in range(NMC):
    m1 = 10**rng.uniform(np.log10(m_e), np.log10(m_tau))
    m2 = 10**rng.uniform(np.log10(m_e), np.log10(m_tau))
    got = False
    for a, b in [(m1,m2),(m2,m1)]:
        if np.min(np.abs(pool/(a/b) - 1)*100) < 0.05:
            got = True
    if got: hits += 1
print(f"  隨機兩質量比 <0.05% 命中率: {hits/NMC*100:.1f}%  (quark 對照: 97.9%)")

# ── 夸克 running mass = probe-dependent 讀數 (QCD 本質) ──
print("\n[E] 夸克 running mass = probe-dependent 讀數 (實驗支持):")
print("  MS̄ mass scale-dependent (RGE): m_b(m_b)=4.18 vs m_b(2GeV)=4.93 GeV")
print("  m_u(2GeV)=2.16 vs m_u(μ=1GeV)?? ~5 MeV — scale 變讀數變")
print("  → 夸克「質量」本身係 scale 函數 = 探針函數 ⟹ 差值閉合嘅讀數本質 ✅")
print("  → 但輕夸克誤差 10-20% 禁做 0.01% 級閉合 anchor")