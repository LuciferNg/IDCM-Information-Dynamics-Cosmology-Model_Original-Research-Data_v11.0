#!/usr/bin/env python3
"""
M4 R3b — m_d/m_u = 11/5 嘅第二實例搜索（2026-10-24）
======================================================
docs/04 T4 命中 m_d/m_u = 2.19929 ± 0.08605 vs L₅/F₅ = 11/5（dev 0.032%, 0.01σ），
但：① 後驗揀出 ② 密度 P1 1.85% ③ 單實例 ⟹ 只能 🟡。
本輪做「第二實例」判別（事前族，同 T4 同一整數比池）：
  族 K-ratio = { p/q : p,q ∈ {1,2,3,4,5,7,8,11,13,21,29,44,76,88}, 1 ≤ p/q ≤ 100 }
  對 FLAG24 Nf=2+1+1 全部夸克比（15 個）逐個掃 |dev|，確認 11/5 係唯一命中抑或成族。
判準：若 <1% 命中數 ≥ 2 且彼此獨立 ⟹ 族現象（可升 🟢 前提之一）；若只有 1 個 ⟹ 單例，維持 🟡。
"""
import numpy as np

K = np.array([1, 2, 3, 4, 5, 7, 8, 11, 13, 21, 29, 44, 76, 88], dtype=float)
fam = np.unique([p/q for p in K for q in K if 1.0 <= p/q <= 100.0])
# FLAG 2024 Nf=2+1+1 + PDG heavy
m = {"u": 2.14e-3, "d": 4.70e-3, "s": 93.46e-3, "c": 1.27, "b": 4.18, "t": 172.69}
sig = {"u": 0.08e-3, "d": 0.05e-3, "s": 0.58e-3, "c": 1.27*0.008, "b": 4.18*0.005, "t": 172.69*0.0030}

print("="*78)
print("M4 R3b: m_d/m_u = 11/5 第二實例搜索（同一整數比族）")
print("="*78)
print(f"族 K-ratio: {len(fam)} 值（域 1–100；K={K.astype(int).tolist()}）")
print(f"已知命中: m_d/m_u = {m['d']/m['u']:.5f} vs 2.2")

names = list(m)
rows = []
for i, a in enumerate(names):
    for b in names[i+1:]:
        hi, lo = (a, b) if m[a] >= m[b] else (b, a)
        t = m[hi]/m[lo]                      # 一律「重/輕」≥ 1
        rel = np.hypot(sig[hi]/m[hi], sig[lo]/m[lo])
        dv = np.abs(fam/t - 1)*100
        j = int(np.argmin(dv))
        nd = dv[j]
        rows.append((nd, f"m_{hi}/m_{lo}", t, fam[j], nd/(rel*100)))
rows.sort()
print("\n|dev| 排序（全部 15 比）:")
for nd, nm, t, fv, ns in rows:
    mark = " ⚠️ <1%" if nd < 1 else ""
    print(f"  {nm:<8} = {t:>12.6g}  最優 K-ratio {fv:>9.5g}  dev {nd:8.4f}%  = {ns:5.2f}σ{mark}")

n1 = sum(1 for r in rows if r[0] < 1.0)
print(f"\n<1% 命中數 = {n1} / 15")
hit1 = [r for r in rows if r[0] < 1.0]
for r in hit1:
    print(f"   {r[1]} = {r[2]:.6g} → {r[3]:.5g} (dev {r[0]:.4f}%, {r[4]:.2f}σ)")
print("\n判：")
if n1 >= 2:
    print(f"  族現象（{n1} 個獨立命中）⟹ 需檢查互為同一結構規則抑或族密度 artifact")
else:
    print("  只有 1 個命中 ⟹ 單例；唔構成第二實例 ⟹ m_d/m_u=11/5 維持 🟡（後驗、單例、scale 依賴未解）")

# ── 族密度：隨機目標喺同一域 [1,100] 落喺族 0.x% 內嘅機率 ⟹ 期望命中數 ──
rng = np.random.default_rng(20261029)
NMC = 1_000_000
tt = 10**rng.uniform(0, 2, NMC)          # log-均勻 [1,100]
logfam = np.log10(fam)
dmin = np.full(NMC, np.inf)
idx = np.searchsorted(logfam, np.log10(tt))
for k in range(-2, 3):
    ii = np.clip(idx+k, 0, fam.size-1)
    dmin = np.minimum(dmin, np.abs(10**(logfam[ii]-np.log10(tt)) - 1)*100)
print(f"\n族密度（域 log-均勻 [1,100]，族 {len(fam)} 值，NMC 1e6）:")
for thr in (0.1702, 0.749, 0.1, 1.0):
    p = float(np.mean(dmin < thr)*100)
    print(f"  thr {thr:>7.4f}%: P = {p:6.2f}%  ⟹ 15 比期望命中數 = {15*p/100:5.2f}")
print("  ⟹ 期望命中 ≈ 2 ⟹ 觀察到 2 個 <1% 命中 = 完全符合隨機族 ⟹ 唔構成「族現象」證據")

# 順帶：lepton 比是否落入同族（μ/e, τ/μ > 族上限 100 → 報告越界）
print("\nlepton 比檢查（族上限 100）:")
for nm, t in [("μ/e", 206.768283), ("τ/μ", 16.817029)]:
    dv = np.abs(fam/t - 1)*100
    print(f"  {nm} = {t:.6f}: 族內最優 {fam[int(np.argmin(dv))]:.5g} dev {dv.min():.3f}%")
