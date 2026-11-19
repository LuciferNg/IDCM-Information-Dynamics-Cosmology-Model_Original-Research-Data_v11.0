#!/usr/bin/python3
"""
M20 — scripts/05_tradeoff_frontier.py
=====================================
把殘留 MC04（「壓低 ρ_d」之機制）**精確化為一條權衡前沿**，並判兩個類之聯合可達性。

類 A（difference-type，= docs/01 凍結族之推廣）：3606 成員
類 B（sum-type，= 加弦 ＋ 鏈上簽名之 gauge-frame 變體）：3376 有效成員
   —— 本類只作**放寬 A6 之後**之測試（見 docs/03 MC04），不作 M20 候選。

輸出：① 兩類之聯合最佳 ② 分支表（精確 Fraction）③ 權衡序列（ρ_d 下降時 |ρ_u| 如何變）
"""

import sys
from fractions import Fraction as F
import math

N, H = 44, 45
KD, KU, KS = 3, 2, 4

def g0(i, j):
    return F(min(i, j) * (H - max(i, j)), H)

def Gfun(i, j, lam, conv):
    if conv == "L":
        u = {k: g0(k, i) - g0(k, j) for k in range(1, N + 1)}
        uTu = g0(i, i) - 2 * g0(i, j) + g0(j, j)
        den = lam + uTu
        if den == 0:
            return None
        return lambda k, m: g0(k, m) - u[k] * u[m] / den
    a, b, c = g0(i, i), g0(i, j), g0(j, j)
    M2 = [[a, lam + b], [lam + b, c]]
    det = M2[0][0] * M2[1][1] - M2[0][1] * M2[1][0]
    if det == 0:
        return None
    Minv = [[M2[1][1] / det, -M2[0][1] / det], [-M2[1][0] / det, M2[0][0] / det]]
    ci = [g0(k, i) for k in range(1, N + 1)]
    cj = [g0(k, j) for k in range(1, N + 1)]
    def Gf(k, m):
        p = [ci[k - 1], cj[k - 1]]
        q = [ci[m - 1], cj[m - 1]]
        return g0(k, m) - sum(p[r] * Minv[r][s] * q[s] for r in range(2) for s in range(2))
    return Gf

TGT_D, TGT_U = F(40350365, 10**8), F(-2588837, 10**8)
SD, SU = F(351, 10000), F(3292, 100000)

clsA, clsB = [], []
for i in range(1, N + 1):
    for j in range(i + 2, N + 1):
        for lam in (1, -1):
            for conv in ("L", "M"):
                G = Gfun(i, j, lam, conv)
                if G is None:
                    continue
                A3, B3 = G(3, 3), G(2, 3)
                A4, B4 = G(3, 4), G(2, 4)
                A2, B2 = G(3, 2), G(2, 2)
                d = A4 - B4
                if d != 0:
                    clsA.append((conv, lam, i, j, (A3 - B3) / d, (A2 - B2) / d))
                s = A4 + B4
                if s != 0:
                    rd, ru = (A3 + B3) / s, -(A2 + B2) / s      # s3=s4, s2=−s3 ⟹ ρ_d=+rd, ρ_u=−ru
                    if rd > 0 and ru < 0:
                        clsB.append((conv, lam, i, j, rd, ru))

def z(rd, ru):
    return float((rd - TGT_D) / SD), float((ru - TGT_U) / SU)

print("=" * 88)
print("M20 / scripts/05 — 權衡前沿：ρ_d 與 ρ_u 之聯合可達性")
print("=" * 88)
for name, cls in (("A difference-type", clsA), ("B sum-type（放寬 A6）", clsB)):
    zs = [z(r[4], r[5]) for r in cls]
    best = min(range(len(cls)), key=lambda k: math.hypot(*zs[k]))
    nhit = sum(1 for a, b in zs if abs(a) < 1 and abs(b) < 1)
    print(f"\n{name}：成員 {len(cls)}")
    print(f"  聯合最佳：χ = {math.hypot(*zs[best]):.3f}  z_d = {zs[best][0]:+.2f}  z_u = {zs[best][1]:+.2f}"
          f"  ← {cls[best][0]} λ={cls[best][1]:+d} ({cls[best][2]},{cls[best][3]})")
    print(f"  ρ_d = {cls[best][4]} = {float(cls[best][4]):.8f}   ρ_u = {cls[best][5]} = {float(cls[best][5]):.8f}")
    print(f"  |z_d|<1 ∧ |z_u|<1 之成員 = {nhit} / {len(cls)}")

# ---------------------------------------------------------------- 分支表
print("\n" + "-" * 88)
print("分支表（精確）：ρ_d 與 ρ_u 之「單獨命中」與其代價")
print("-" * 88)
branches = [
    ("A 支線：ρ_u 命中", "M", +1, 1, 4, "A"),
    ("B 支線：ρ_d→16/41", "M", -1, 1, 4, "B"),
    ("B 支線：ρ_d→2/41", "M", -1, 1, 5, "B"),
    ("B 支線：(2,4) sum", "M", -1, 2, 4, "B"),
    ("A 支線：ρ_d→85/41", "L", +1, 2, 4, "A"),
]
for label, conv, lam, i, j, which in branches:
    src = clsA if which == "A" else clsB
    hit = [r for r in src if r[0] == conv and r[1] == lam and r[2] == i and r[3] == j]
    if not hit:
        continue
    rd, ru = hit[0][4], hit[0][5]
    zd, zu = z(rd, ru)
    print(f"  {label:<18} {conv} λ={lam:+d} ({i},{j})  ρ_d={str(rd):>10} (z {zd:+.2f})  ρ_u={str(ru):>10} (z {zu:+.2f})")

# ---------------------------------------------------------------- 權衡序列
print("\n" + "-" * 88)
print("權衡：ρ_d 落入目標鄰域（0.36–0.46）之成員，其 ρ_u 之分佈（精確值集）")
print("-" * 88)
band = [r for r in clsA + clsB if F(36, 100) <= r[4] <= F(46, 100) and r[5] < 0]
print(f"  成員數 = {len(band)}")
vu = sorted({r[5] for r in band})
print(f"  ρ_u 取值集（{len(vu)} 個）：")
for v in vu:
    print(f"    {str(v):>14} = {float(v):+.6f}   |ρ_u|/ρ_d(取 ρ_d=16/41) = {float(abs(v)*41/16):.4f}")
tgt_ratio = float(abs(TGT_U) / TGT_D)
print(f"\n  目標比值 |ρ_u|/ρ_d = {tgt_ratio:.6f}")
print(f"  但在 ρ_d ≈ 16/41 之成員，最小 |ρ_u| = {float(min(abs(v) for v in vu)):.6f}"
      f"（= 目標 |ρ_u| = {float(abs(TGT_U)):.6f} 之 {float(min(abs(v) for v in vu))/float(abs(TGT_U)):.1f}×）")

print("\n" + "=" * 88)
print("判：MC04 之精確形狀")
print("=" * 88)
print("  · ρ_d 單獨可達（16/41 z −0.38）｜ρ_u 單獨可達（−1/41 z +0.05）")
print("  · 聯合不可達（A 類 0/3606、B 類 0/3376；最佳 χ = 3.10 / 4.71）")
print("  · 機制：ρ_d 一壓低，|ρ_u| 急升 ⟹ 缺者 = **ρ_u 平坦之壓低機制**（非再多一個算子）")
print("=" * 88)
sys.exit(0)
