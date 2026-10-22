#!/usr/bin/env python3
"""
/tmp 夸克差值閉合 02 — m_c/m_s = γ20/γ5 (0.042%) MC null 判別 + 盲測
=========================================================================
01 掃描唯一 0.1% 級命中: m_c/m_s = 13.5974 vs γ20/γ5 = 13.6032 (0.042%)

判別:
  1. m_c/m_s 係 15 個質量比入面唯一 <0.1% 命中 → 32/30 比 95% 唔夠 (Bonferroni)
  2. 但「穩健性」: γ20/γ5 嘅 mode 結構有冇身份?
     20/5 = 4 = B! γ20 = γ_{5B}, γ5 = γ_{F₅}?
     γ20: 20 = 2B·F₅ = 4×5 = B×F₅? 20 = B·F₅ = 4·5 ✓
     γ5:  5 = F₅ ✓
     → m_c/m_s ~ γ_{B·F₅}/γ_{F₅}?? 試其他 (B·k)/k 比:
        γ8/γ2 = ?, γ12/γ3 = ?, γ16/γ4 = ?, γ20/γ5 = ?, γ24/γ6 = ?, γ28/γ7, γ32/γ8, γ36/γ9, γ40/γ10, γ44/γ11
  3. MC null: 隨機質量比 (夸克質量範圍 log-均勻) 命中率
"""
import numpy as np
import random

N = 88
n = np.arange(1, 45, dtype=float)
gamma = 4*np.sin(n*np.pi/N)**2
inv_g = 1/gamma

m_c = 1.27; m_s = 0.0934
tgt = m_c/m_s
g20 = gamma[19]; g5 = gamma[4]
print("="*64)
print("誇克 02: m_c/m_s ~ γ20/γ5 判別")
print("="*64)
print(f"m_c/m_s = {tgt:.6f}")
print(f"γ20/γ5  = {g20/g5:.6f}  dev {abs(g20/g5/tgt-1)*100:.6f}%")
print(f"γ20 = {g20:.9f} (mode 20), γ5 = {g5:.9f} (mode 5)")
print(f"mode 20 = B·F₅ = 4×5 ✓; mode 5 = F₅ ✓")

# (B·k)/k 譜比族: k = 1..11
print("\n[A] (B·k)/k 比例族 (mode 4k vs k):")
for k in range(1, 12):
    a = 4*k
    if a > 44: break
    r = gamma[a-1]/gamma[k-1]
    print(f"  γ{4*k}/γ{k} = {r:.6f}  (m_c/m_s dev {abs(r/tgt-1)*100:.3f}%)")

# 全 (ab)/b 廣義: 32768 組合太大量, 試 (a·k)/(b·k)? 其實 γ_{ak}/γ_{bk}?
print("\n[B] 廣義保比族: γ_{ak}/γ_{bk} (a/b = 20/5 = 4):")
for k in range(1, 11):
    a, b = 4*k, 1*k
    if a > 44 or b > 44: break
    r = gamma[a-1]/gamma[b-1]
    print(f"  γ{4*k}/γ{k} = {r:.6f}  dev {abs(r/tgt-1)*100:.3f}%")

# m_c/m_s 幾何位置: log 空間 quarks 有冇「比例階梯」?
print("\n[C] 夸克 log 階梯:")
names = ["u","d","s","c","b","t"]
ms = [0.00216, 0.00467, 0.0934, 1.27, 4.18, 172.69]
logm = np.log10(ms)
for i in range(5):
    ratio = ms[i+1]/ms[i]
    print(f"  {names[i]}→{names[i+1]}: {names[i+1]}/{names[i]} = {ratio:9.4f}  log10 = {logm[i+1]-logm[i]:.4f}")

# 大問題: 呢 5 個比有冇 γ 譜身份?
print("\n[D] 5 個相鄰比 vs γ比池 top 命中:")
rat = gamma[:, None]/gamma[None, :]
pairs = [("u,d", 0.00467/0.00216), ("d,s", 0.0934/0.00467), ("s,c", 1.27/0.0934),
         ("c,b", 4.18/1.27), ("b,t", 172.69/4.18)]
for name, t in pairs:
    dev = np.abs(rat/t - 1)*100
    best = np.unravel_index(np.argmin(dev), dev.shape)
    print(f"  {name}: ratio={t:8.4f}  γ{best[0]+1}/γ{best[1]+1}={rat[best]:.4f} dev {dev[best]:.4f}%")

# MC null: 隨機夸克質量比 (log-均勻 6 質量) → γ比池 <0.1% 命中率
print("\n[E] MC null: 隨機 6 夸克質量 (log10 ∈ [log10(2e-3), log10(173)])")
rng = random.Random(25721)
hits = 0
NMC = 5000
pool = rat.flatten()
for _ in range(NMC):
    m6 = np.sort(10**np.array([rng.uniform(np.log10(0.002), np.log10(173)) for _ in range(6)]))
    # 15 個比
    got = False
    for i in range(6):
        for j in range(i+1, 6):
            r = m6[j]/m6[i]
            if np.min(np.abs(pool/r - 1)*100) < 0.1:
                got = True
    if got: hits += 1
print(f"  隨機夸克族 15 比中至少一個 <0.1% γ比命中: {hits/NMC*100:.1f}%")
print(f"  (若 >10%: m_c/m_s 命中唔特別)")