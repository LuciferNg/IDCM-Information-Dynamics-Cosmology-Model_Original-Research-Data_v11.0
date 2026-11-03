#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/29 — /88 之最細子鏈 /8（W1）＋ 二次域判準：ℚ(√5) 為何只屬 /90

────────────────────────────────────────────────────────────────────────
殘留：docs/27 §五 W1（/88 之 d=11 ⟹ /8 鏈 4 模之共振結構）
延伸：J1（為何恰 6 個 ℚ(√5) 角）之**域論根源**

判準（跑前凍結）
  S1 子鏈 /8：n = 11k，k = 1..4 ⟹ γ = 4sin²(kπ/8)；（2−√2, 2, 2+√2, 4）
  S2 共振：γ_a ± γ_b = γ_m，a,b,m ∈ {1,2,3,4} 互異（能量）；四波 = 兩不同價對同 S。
  S3 域判準（Kronecker–Weber）：ℚ(√d) ⊂ ℚ(ζ_n) ⟺ conductor(d) | n；
     conductor(√2) = 8、conductor(√5) = 5。
  S4 pslq 核最小多項式（2 次）：disc 8 ⟹ ℚ(√2)；disc 5 ⟹ ℚ(√5)。
性質：確定性。N_trials = 0。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools

import mpmath as mp

mp.mp.dps = 60
L = "=" * 78
MODES = [11, 22, 33, 44]          # /88 之 d=11 子鏈
KMODES = [1, 2, 3, 4]             # 對應 /8 鏈之 k


def gam(x, n):
    return 4 * mp.sin(n * mp.pi / x) ** 2


print(L)
print("M14 scripts/29 — /88 最細子鏈 /8（W1）＋ 二次域判準")
print(L)

# ============================================================================
print("\n【§一】/88 之 d=11 子鏈：n = 11k，γ = 4sin²(kπ/8)")
g = {k: gam(88, 11 * k) for k in KMODES}
for k in KMODES:
    print(f"    k={k}  n={11*k:>2}  γ = {mp.nstr(g[k], 16)}")
dev2 = max(abs(g[k] - 4 * mp.sin(k * mp.pi / 8) ** 2) for k in KMODES)
print(f"  γ_(11k) = 4sin²(kπ/8)：max|dev| = {mp.nstr(dev2, 3)}")
print(f"  精確值：2−√2 = {mp.nstr(2 - mp.sqrt(2), 16)}；2+√2 = {mp.nstr(2 + mp.sqrt(2), 16)}")
print(f"  核：|γ1 − (2−√2)| = {mp.nstr(abs(g[1] - (2 - mp.sqrt(2))), 3)}；"
      f"|γ3 − (2+√2)| = {mp.nstr(abs(g[3] - (2 + mp.sqrt(2))), 3)}")
print("  ⟹ 全 4 模皆在 **ℚ(√2)**")

# ============================================================================
print("\n【§二】/8 子鏈之共振（3 波、4 波）")
TOL = mp.mpf('1e-50')
three = []
for a, b, m in itertools.permutations(KMODES, 3):
    if abs(g[a] + g[b] - g[m]) < TOL:
        three.append(tuple(sorted((a, b, m))))
three = sorted(set(three))
print(f"  三波共振（γ_a+γ_b=γ_m，互異）= {len(three)} 個：{three}")
for t in three:
    a, b, m = t
    print(f"    驗算：γ{a}+γ{b} = {mp.nstr(g[a]+g[b], 16)}；γ{m} = {mp.nstr(g[m], 16)}")
print(f"    奇偶：{'偶' if sum(three[0]) % 2 == 0 else '奇'} ⟹ χ³ "
      f"{'暗' if sum(three[0]) % 2 == 0 else '亮'}" if three else "  （無）")
pairs = {}
for a, b in itertools.combinations(KMODES, 2):
    pairs[(a, b)] = g[a] + g[b]
vals = list(pairs.values())
dup = 0
for i in range(len(vals)):
    for j in range(i + 1, len(vals)):
        if abs(vals[i] - vals[j]) < TOL:
            dup += 1
print(f"  四波共振（兩不同價對同 S）= {dup} 個（6 個價對之 S 全互異 ⟹ 0）")
print(f"  價對 S：{[(k, mp.nstr(v, 8)) for k, v in pairs.items()]}")

# ============================================================================
print("\n【§三】★ 二次域判準（S3）：ℚ(√5) ⊄ ℚ(ζ_88)，但 ⊂ ℚ(ζ_90)")
print(L)
rows = [("√2", 8, 88, 90), ("√5", 5, 88, 90), ("√−11", 11, 88, 90), ("√−3", 3, 88, 90)]
print(f"  {'域':<8}{'conductor':>10}{'| 88 ?':>9}{'| 90 ?':>9}")
for name, cond, n1, n2 in rows:
    print(f"  {name:<8}{cond:>10}{'✓' if n1 % cond == 0 else '✗':>9}"
          f"{'✓' if n2 % cond == 0 else '✗':>9}")
print(f"\n  88 = 2³·11 ⟹ conductor 8 | 88 ✓（√2 在內）、11 | 88 ✓（√−11 在內）、**5 ∤ 88 ⟹ √5 在外**")
print(f"  90 = 2·3²·5 ⟹ conductor 5 | 90 ✓（√5 在內）、3 | 90 ✓（√−3 在內）")
print(f"  ⟹ /90 之 6 個 ℚ(√5) 角（18°,30°,36°,54°,60°,72° ⟹ γ ∈ ℚ(√5)）"
      f"**不可能**出現於 /88 之譜")
print(f"  ⟹ 此即 J1「為何恰 6 個」之**域論根源**：√5 只在 5 | n 之鏈上可表")

# ============================================================================
print("\n【§四】pslq 核最小多項式（S4）：兩個不同二次域")
print(L)


def minpoly(x):
    rel = mp.pslq([mp.mpf(1), x, x * x], maxcoeff=10 ** 12, maxsteps=10 ** 6)
    if rel is None:
        return None
    # 化為首一整數係數
    from math import gcd
    from functools import reduce
    a0, a1, a2 = [int(v) for v in rel]
    gg = reduce(gcd, [abs(a0), abs(a1), abs(a2)])
    a0, a1, a2 = a0 // gg, a1 // gg, a2 // gg
    if a2 < 0:                                  # 首一化（符號正規化）
        a0, a1, a2 = -a0, -a1, -a2
    disc = a1 * a1 - 4 * a0 * a2
    return (a0, a1, a2, disc)


for tag, x in (("/88 子鏈 γ1 = 2−√2", g[1]), ("/90 三元組 γ9 = (3−√5)/2", gam(90, 9)),
               ("/90 三元組 γ15 = 1", gam(90, 15)), ("/88 子鏈 γ3 = 2+√2", g[3])):
    r = minpoly(x)
    print(f"  {tag:<28} pslq 關係 {r[2]}x² + ({r[1]})x + ({r[0]})，disc = {r[3]}")
print(f"""  ⟹ disc 8 ⟹ 域 = ℚ(√2)（/88 側）
  ⟹ disc 5 ⟹ 域 = ℚ(√5)（/90 側）
  ⟹ γ15 = 1 為有理 ⟹ pslq 給退化關係 x²−x（disc 1，非平方自由）⟹ 屬 ℚ，非二次域
  ⟹ **兩鏈之譜落在不同二次域**，非同一數域之不同元素""")

# ============================================================================
print("\n【§五】判定（W1）")
print(L)
print(f"""  ✅ W1 得答：/88 最細子鏈 = /8（d=11，4 模：n = 11,22,33,44），γ ∈ {{2−√2, 2, 2+√2, 4}}
     —— 全部在 **ℚ(√2)**。
     三波共振 **{len(three)}** 個（{three}），四波共振 **{dup}** 個 ⟹ 4 模扇區近乎無共振。
  ✅ 域論結論（新）：ℚ(√5) 之 conductor 5 ∤ 88 ⟹ √5 ∉ ℚ(ζ_88)；
     而 5 | 90 ⟹ √5 ∈ ℚ(ζ_90)。pslq 核實 disc 8 vs disc 5。
  ⟹ **「ℚ(√5) 扇區」係 /90 專屬**；/88 最深只到 ℚ(√2)。
     此即 J1（恰 6 個 ℚ(√5) 角）之根源，並解釋 docs/15/dpcs17 之 φ-階梯為何只在 /90 出現。
  🔴 未做：/88 之 ℚ(√2) 扇區有無自身之「階梯」（類比 φ-階梯之 1:√2:2?）。
""")
print(L)
raise SystemExit(0)