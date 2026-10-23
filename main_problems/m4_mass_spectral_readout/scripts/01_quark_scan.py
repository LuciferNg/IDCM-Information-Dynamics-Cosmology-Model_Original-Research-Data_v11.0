#!/usr/bin/env python3
"""
/tmp 夸克差值閉合測試 01 — 質量數據 + 譜比掃描 (2026-10-21)
==============================================================
相同方法 (R7 定案): 質量 = 譜讀數差值/比, 唔 fit 絕對標度。

夸克質量 (PDG 2024, MS̄):
  m_u = 2.16 MeV,  m_d = 4.67 MeV,  m_s = 93.4 MeV        @ 2 GeV
  m_c = 1.27 GeV                                          @ m_c
  m_b = 4.18 GeV                                          @ m_b
  m_t = 172.69 GeV                                        @ pole

問: 質量比 (m_i/m_j) 有冇 γ_n 譜結構? (同 R7 A/B/C/D 族)
    同步邊界讀數 = γ_n? 1/γ_n? γ_n²? Σ/γ?
"""
import numpy as np
import random

# ── 夸克質量 (GeV) ──
m_u = 0.00216; m_d = 0.00467; m_s = 0.0934; m_c = 1.27; m_b = 4.18; m_t = 172.69
quarks = {"u": m_u, "d": m_d, "s": m_s, "c": m_c, "b": m_b, "t": m_t}
print("="*70)
print("夸克差值閉合 01: 質量比 vs 譜結構 (事前掃描)")
print("="*70)

# ── γ 譜 + 擴展結構數池 ──
N = 88
n = np.arange(1, 45, dtype=float)
gamma = 4*np.sin(n*np.pi/N)**2
inv_g = 1/gamma

# 事前結構池 (純 K_Chain 量, 冇 M_comp)
S = {
    "γ6": gamma[5], "γ4": gamma[3], "γ5": gamma[4], "γ44": gamma[43],
    "Σγ": 88.0, "Πγ": 176.0, "B": 4.0, "Nd": 44.0, "N88": 88.0,
    "F5": 5.0, "F7": 13.0, "F8": 21.0, "L5": 11.0, "L7": 29.0, "L9": 76.0,
}
SNAME = {v: k for k, v in S.items()}

# 質量比 (全部 15 對)
names = list(quarks.keys())
print("\n[0] 夸克質量比矩陣 (log10):")
print("     " + "".join(f"{q:>10}" for q in names))
for i, qi in enumerate(names):
    row = f"  {qi:>2}: "
    for j, qj in enumerate(names):
        if i == j:
            row += f"{'-':>10}"
        else:
            row += f"{np.log10(quarks[qi]/quarks[qj]):>10.4f}"
    print(row)

# ── [A] 質量比 = γ_a/γ_b 掃描 (6 choose 2 = 15 比) ──
print("\n[A] 質量比 = γ_a/γ_b")
ga = gamma[:, None]; gb = gamma[None, :]
rat = ga/gb
pairs = [("t","b"), ("b","c"), ("c","s"), ("s","d"), ("d","u"), ("s","u"),
         ("t","c"), ("t","s"), ("c","d"), ("b","s")]
for qi, qj in pairs:
    tgt = quarks[qi]/quarks[qj]
    # 直接 γ 比 (範圍 0.008-322)
    dev = np.abs(rat/tgt - 1)*100
    best = np.unravel_index(np.argmin(dev), dev.shape)
    bd = dev[best]
    # 1/γ 差 (範圍 0-7.9)
    ia = inv_g[:, None] - inv_g[None, :]
    pos = ia[ia > 0]
    if pos.size:
        db = np.abs(pos/tgt - 1)*100
        bdb = db.min()
    else:
        bdb = 1e9
    # 結構常數比
    bestc = None
    for v1, k1 in SNAME.items():
        for v2, k2 in SNAME.items():
            if v2 != 0:
                dc = abs(v1/v2/tgt - 1)*100
                if bestc is None or dc < bestc[0]:
                    bestc = (dc, f"{k1}/{k2}")
    print(f"  m_{qi}/m_{qj} = {tgt:>12.6g}:  γ比最優 γ{best[0]+1}/γ{best[1]+1}={rat[best]:.6g} dev {bd:.4f}% | "
          f"1/γ差最優 dev {bdb:.4f}% | 常量最優 {bestc[1]} dev {bestc[0]:.4f}%")

# ── [B] 用「同步邊界」 = γ_n 喺 88 譜做前 10 命中 ──
print("\n[B] 跨比 look-elsewhere 密度 (全部 15 比, γ比池):")
pool = rat.flatten()
all_best_dev = []
for qi in names:
    for qj in names:
        if qi == qj: continue
        tgt = quarks[qi]/quarks[qj]
        dev = np.abs(pool/tgt - 1)*100
        all_best_dev.append((dev.min(), f"m_{qi}/m_{qj}"))
all_best_dev.sort()
print("  每比最優 dev (細 = 可疑):")
for d, name in all_best_dev:
    flag = " ⚠️" if d < 0.1 else ""
    print(f"    {name}: {d:.4f}%{flag}")
print(f"  |dev|<0.1% 嘅比數: {sum(1 for d,_ in all_best_dev if d<0.1)}/30")
print(f"  |dev|<1% 嘅比數: {sum(1 for d,_ in all_best_dev if d<1)}/30")

# ── [C] m_q = γ_n²·A_q 讀數×譜強度 (同 m₃=γ₆²A) ──
print("\n[C] m_q = γ_n²·A_q 讀數×譜強度形式:")
print("  ⚠️ 修正 (2026-10-24 R1 重建): /tmp 原版呢段用 A_q=m/γ² 嘅 44 值向量 同")
print("     inv_g / gamma 向量做「同 index 配對」再格式化向量 ⟹ TypeError 全支 crash。")
print("     改為 well-posed 自由 index 掃描: ∃(n,k) 使 γ_n²·c_k = m_q, c_k ∈ S ∪ {γ_k} ∪ {1/γ_k}")
g2 = gamma**2
cands = {}
for i, g in enumerate(gamma, 1):
    cands[f"γ{i}"] = g
    cands[f"1/γ{i}"] = 1.0/g
for k, val in S.items():
    cands[k] = val
cnames = list(cands)
cvals = np.array([cands[k] for k in cnames])
print(f"  掃描空間: 44 mode × {len(cnames)} 結構候選 = {len(cnames)*44} 組合/夸克")
for qname, m in quarks.items():
    dev = np.abs(g2[:, None]*cvals[None, :]/m - 1)*100
    bi = np.unravel_index(np.argmin(dev), dev.shape)
    Aq = m/g2
    print(f"  m_{qname} ({m:>8g} GeV): A_q=m/γ² ∈ [{Aq.min():.3g}, {Aq.max():.3g}] GeV | "
          f"最優 γ{bi[0]+1}²·{cnames[bi[1]]} = {g2[bi[0]]*cvals[bi[1]]:.6g} dev {dev[bi]:.4f}% | "
          f"<0.1% 組合 {int((dev<0.1).sum())} | <1% 組合 {int((dev<1).sum())}")
rng_c = random.Random(25711)
NMC_C = 200
h01 = h1 = 0
for _ in range(NMC_C):
    mm = 10**rng_c.uniform(np.log10(0.002), np.log10(173))
    dv = np.abs(g2[:, None]*cvals[None, :]/mm - 1)*100
    if dv.min() < 0.1: h01 += 1
    if dv.min() < 1.0: h1 += 1
print(f"  MC null ({NMC_C} 隨機夸克質量, log10∈[log10 2e-3, log10 173]): "
      f"best<0.1% {h01/NMC_C*100:.1f}% | best<1% {h1/NMC_C*100:.1f}%")
print("  ⟹ 密度閘: 池大到單點 <1% 命中 = 隨機級; A_q 逐夸克唔同值 ⟹ 零預言力")

print("\n" + "="*40)
print("Quark test 01 done — 見上")