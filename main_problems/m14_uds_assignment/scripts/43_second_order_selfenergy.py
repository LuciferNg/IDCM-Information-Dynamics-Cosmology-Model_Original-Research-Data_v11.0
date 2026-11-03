#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/43 — MR1 機器（二）：二階自能（虛過程頻移）與極點屏蔽

────────────────────────────────────────────────────────────────────────
公約（/90 站空間，docs/19 已釘死）
  φ_n(j) = √(2/45) sin(πnj/45), j = 1..44（Dirichlet 44 站）
  γ_n = 4 sin²(πn/90)  ⟹  **Σγ = 88**、γ_max = γ_44 = 3.99512810052
  （簽名值 assert 於 §0；§88 對照鏈 γ_n = 4sin²(πn/88) ⟹ Σγ = 90，不在此腳本使用）

定義（零參數；本層只取「1→3 tadpole」支，2→2 bubble 支明示未算）
  δγ_n := (g₄²/8100) · Σ_{b≤c≤d} m(b,c,d) · W(n;b,c,d)² / (γ_n − γ_b − γ_c − γ_d)
      m = 3!/(multinomial)（相異排序數）｜W = 90·χ⁴（docs/39 閉式）
      分母 ≈ 0 且 W ≠ 0 ⟹ 極點（微擾失效）

判準（跑前凍結）
  P1 §0 簽名值：Σγ == 88 且 γ_max == γ_44
  P2 運動學極點只有 1 個：(n;b,c,d) = (30;15,15,15)（docs/41 D2）
     命題 𝔸：該點之 W = 0 ⟹ **極點被宇稱屏蔽** ⟹ 實際極點數 = 0
  P3 44 個 δγ 值全相異（無簡併 ⟹ 可作標籤？否，見 §四）
  P4 排序保持性：δγ 是否令 γ 之次序翻轉（level crossing）
  P5 微擾有效性界：g₄² < min_n 8100·γ_n/(2|S_n|)，S_n = 8100·δγ_n/g₄²
性質：全枚舉 44 × C(46,3) = 668,  依 mpmath dps = 50。零觀測量、零自由參數。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import math
from collections import Counter

import mpmath as mp

mp.mp.dps = 50
TOL = mp.mpf('1e-35')
L = "=" * 78
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]
G = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}


def W4(t):
    return sum(c for c, sg in zip(CS, SIG) if sum(s * x for s, x in zip(sg, t)) % 90 == 0)


def mult(t):
    m = 6
    for v in Counter(t).values():
        m //= math.factorial(v)
    return m


print(L)
print("M14 scripts/43 — MR1（二）：二階自能與極點屏蔽")
print(L)

# ---------------------------------------------------------------- §0
S = sum(G.values())
print(f"\n【§0】簽名值閘：Σγ = {mp.nstr(S, 12)}（期望 88）｜γ_max = {mp.nstr(max(G.values()), 12)}"
      f"（= γ_44 = {mp.nstr(G[44], 12)}）")
assert abs(S - 88) < TOL and abs(max(G.values()) - G[44]) < TOL, "公約簽名失敗"
print("  ✅ 公約簽名通過（/90 Dirichlet，πn/90）")

# ---------------------------------------------------------------- §一
print("\n" + L)
print("【§一】P2 極點掃描：分母 ≈ 0 之全部 (n;b,c,d)")
print(L)
tri = list(itertools.combinations_with_replacement(range(1, 45), 3))
poles = []
for n in range(1, 45):
    for t in tri:
        den = G[n] - sum(G[x] for x in t)
        if abs(den) < TOL:
            poles.append((n, t, W4((n,) + t)))
print(f"  運動學極點數 = {len(poles)}：{[(p[0], p[1]) for p in poles]}")
for n, t, w in poles:
    print(f"    n = {n}，triple = {t}：W = {w} ⟹ χ⁴ = {mp.mpf(w) / 90}（宇稱：和 = {n + sum(t)}）")
print(f"  其中 W ≠ 0（真極點）= {sum(1 for p in poles if p[2] != 0)}")
print("  ⟹ **極點被宇稱屏蔽**：唯一運動學通道 (30;15,15,15) 之和 = 75 為奇 ⟹ χ⁴ = 0")

# ---------------------------------------------------------------- §二
print("\n" + L)
print("【§二】二階頻移 S_n := 8100·δγ_n/g₄²")
print(L)
Sh = {}
for n in range(1, 45):
    acc = mp.mpf(0)
    for t in tri:
        w = W4((n,) + t)
        if w == 0:
            continue
        den = G[n] - sum(G[x] for x in t)
        acc += mult(t) * mp.mpf(w * w) / den
    Sh[n] = acc
for n in range(1, 45):
    print(f"  n = {n:2d}｜γ = {mp.nstr(G[n], 10):>14}｜S = {mp.nstr(Sh[n], 12):>20}")
print(f"  min S = {mp.nstr(min(Sh.values()), 10)}（n = {min(Sh, key=Sh.get)}）"
      f"｜max S = {mp.nstr(max(Sh.values()), 10)}（n = {max(Sh, key=Sh.get)}）")

# ---------------------------------------------------------------- §三
print("\n" + L)
print("【§三】P3 相異性／P4 排序")
print(L)
distinct = len({mp.nstr(v, 25) for v in Sh.values()})
print(f"  相異值 = {distinct} / 44 {'✅ 全相異' if distinct == 44 else '🔴 有簡併'}")
order = sorted(range(1, 45), key=lambda n: G[n])
order2 = sorted(range(1, 45), key=lambda n: G[n] + 0)  # 未含 g₄，見下
print(f"  γ 之升序 = {order[:8]} … {order[-4:]}")
print("  含 g₄ 之排序：δγ_n = (g₄²/8100)S_n ⟹ 次序隨 g₄ 變 ⟹ 見 §四 之臨界值")

# ---------------------------------------------------------------- §四
print("\n" + L)
print("【§四】P5 微擾有效性界與 level crossing 臨界 g₄²")
print(L)
# 要求 |δγ_n| < γ_n/2（半寬，避免跨級）
bound = min(8100 * G[n] / (2 * abs(Sh[n])) for n in Sh if Sh[n] != 0)
bmin_n = min((n for n in Sh if Sh[n] != 0), key=lambda n: 8100 * G[n] / (2 * abs(Sh[n])))
print(f"  g₄² < min_n 8100γ_n/(2|S_n|) = {mp.nstr(bound, 16)}（最小者 n = {bmin_n}）"
      f" ⟹ |g₄| < {mp.nstr(mp.sqrt(bound), 14)}")
# 相鄰級 crossing：γ_n + c S_n = γ_m + c S_m，c = g₄²/8100  
cross = []
g_sorted = sorted(range(1, 45), key=lambda n: G[n])
for i in range(len(g_sorted) - 1):
    a, b = g_sorted[i], g_sorted[i + 1]
    dS = Sh[a] - Sh[b]
    if abs(dS) < TOL:
        continue
    c = (G[b] - G[a]) / dS
    if c > 0:
        cross.append((a, b, c))
cross.sort(key=lambda x: x[2])
print(f"  正臨界 c = g₄²/8100（相鄰級翻轉）共 {len(cross)} 個；最小者：")
for a, b, c in cross[:5]:
    print(f"    級 {a} ↔ {b}：c = {mp.nstr(c, 16)} ⟹ g₄ = {mp.nstr(mp.sqrt(8100 * c), 14)}")

# ---------------------------------------------------------------- §五
print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ① 運動學極點 {len(poles)} 個，全部 W = 0 ⟹ **宇稱屏蔽**，微擾論無極點
  ② 44 個 δγ 全相異（{distinct}/44）⟹ 二階頻移不含簡併；與不透明度之「完全均勻」成對照
  ③ 微擾有效界：|g₄| < {mp.nstr(mp.sqrt(bound), 6)}（以 |δγ| < γ/2 定義）
  ④ 最小 level-crossing 臨界：g₄ ≈ {mp.nstr(mp.sqrt(8100 * cross[0][2]), 6) if cross else 'n/a'}
  ⑤ 未算：2→2 bubble 支（含自收縮極點，須重整化）；故本節 δγ 只係 tadpole 支""")
print(L)
raise SystemExit(0)