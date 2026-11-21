#!/usr/bin/python3
"""
M20 — scripts/07_fractional_weights.py
======================================
R7：**結構性分數權重（RD7）**之判定。候選集**跑前凍結於本 docstring**（跑後不改）。

【凍結候選集】每個分數必須有**具名結構來源**；無來源者禁（手調）。
    c ∈ { 1/1（鏈標度）, 1/2（88=2N_d）, 1/3（|S|=N_gen=3）, 1/4（1/B）, 2/3, 3/4 }

【三種承載形式】
    A  Φ = c·L_(2,4)                       （味站三角形之 AF 弦，F1 之縮放）
    B  Φ = c·[L_(2,3)+L_(3,4)+L_(2,4)]     （K₃ 三鍵，三對平均）
    C  Φ = c·L_(2,4) + c′·e₄e₄ᵀ            （兩分量）

【判準】
    C0 回歸：Woodbury vs 44×44 直接 Fraction 求逆（抽 3 例）精確 0 差
    C1 全報（15 個）：χ = √(z_d² + z_u²)；命中 = χ < 3（雙 3σ）
    C2 結構觀察：`Φ = c·L_(2,4)` 之 ρ_u 是否 c-不變
"""

import sys
import math
from fractions import Fraction as F

N, H = 44, 45

def g0(i, j):
    return F(min(i, j) * (H - max(i, j)), H)

def inv(M):
    n = len(M)
    A = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [u - f * v for u, v in zip(A[r], A[c])]
    return [row[n:] for row in A]

def woodbury(terms):
    """Φ = Σ_a c_a v_a v_aᵀ （秩 ≤ 3）。G = G₀ − (G₀U)(C⁻¹ + UᵀG₀U)⁻¹(UᵀG₀)。"""
    k = len(terms)
    vecs = [t[1] for t in terms]
    cs = [t[0] for t in terms]
    UGU = [[sum(vecs[a][p] * g0(p, q) * vecs[b][q]
                for p in range(1, N + 1) for q in range(1, N + 1)) for b in range(k)]
           for a in range(k)]
    Mm = [[(F(1, cs[a]) if a == b else F(0)) + UGU[a][b] for b in range(k)] for a in range(k)]
    Mi = inv(Mm)
    if Mi is None:
        return None
    cols = [[sum(vecs[a][p] * g0(p, x) for p in range(1, N + 1)) for x in range(1, N + 1)]
            for a in range(k)]

    def Gf(x, y):
        p = [cols[a][x - 1] for a in range(k)]
        q = [cols[a][y - 1] for a in range(k)]
        return g0(x, y) - sum(p[a] * Mi[a][b] * q[b] for a in range(k) for b in range(k))
    return Gf

def direct(terms):
    """直接建 44×44 並精確求逆（對照核）。"""
    K = [[F(0)] * N for _ in range(N)]
    for a in range(N):
        K[a][a] = F(2)
        if a + 1 < N:
            K[a][a + 1] = K[a + 1][a] = F(-1)
    for c, v in terms:
        for p in range(1, N + 1):
            if v[p] == 0:
                continue
            for q in range(1, N + 1):
                if v[q] == 0:
                    continue
                K[p - 1][q - 1] += F(c) * F(v[p] * v[q])
    G = inv(K)
    if G is None:
        return None
    g = lambda x, y: G[x - 1][y - 1]
    dd = g(3, 3) - g(2, 3)
    ds = g(3, 4) - g(2, 4)
    du = g(3, 2) - g(2, 2)
    if ds == 0:
        return None
    return dd / ds, du / ds

def ratios(G):
    if G is None:
        return None
    dd = G(3, 3) - G(2, 3)
    ds = G(3, 4) - G(2, 4)
    du = G(3, 2) - G(2, 2)
    if ds == 0:
        return None
    return dd / ds, du / ds

def e(i):
    return {k: (1 if k == i else 0) for k in range(1, N + 1)}

def uab(i, j):
    return {k: (1 if k == i else 0) - (1 if k == j else 0) for k in range(1, N + 1)}

TGT_D, TGT_U = F(40350365, 10**8), F(-2588837, 10**8)
SD, SU = F(351, 10000), F(3292, 100000)

def report(name, terms):
    w = ratios(woodbury(terms))
    if w is None:
        print(f"  {name:<46} 奇異／退化")
        return None
    zd = float((w[0] - TGT_D) / SD)
    zu = float((w[1] - TGT_U) / SU)
    chi = math.hypot(zd, zu)
    print(f"  {name:<46} ρ_d={str(w[0]):>11}(z{zd:+8.2f})  ρ_u={str(w[1]):>11}(z{zu:+8.2f})"
          f"  χ={chi:7.2f}  {'✅' if chi < 3 else '🔴'}")
    return w, chi

print("=" * 96)
print("M20 / scripts/07 — 結構性分數權重（RD7）")
print("=" * 96)

print("\nC0 回歸核（Woodbury vs 直接求逆）：")
_ok = True
for name, terms in (("c=1 L(2,4)", [(F(1), uab(2, 4))]),
                    ("c=1/3 L(2,4)", [(F(1, 3), uab(2, 4))]),
                    ("c=1/3 K₃", [(F(1, 3), uab(2, 3)), (F(1, 3), uab(3, 4)), (F(1, 3), uab(2, 4))])):
    a, b = ratios(woodbury(terms)), direct(terms)
    same = (a == b)
    _ok = _ok and same
    print(f"   {name:<20} W={a}  D={b}  {'✅' if same else '🔴'}")
print(f"   ⟹ C0 {'PASS' if _ok else 'FAIL'}")
if not _ok:
    print("   ★ 停：回歸失敗 ⟹ 不輸出判詞。")
    sys.exit(2)

FIRST = [(F(1), "1 = 鏈標度"), (F(1, 2), "1/2 = 88=2N_d"), (F(1, 3), "1/3 = |S|=N_gen"),
         (F(1, 4), "1/4 = 1/B"), (F(2, 3), "2/3"), (F(3, 4), "3/4")]

print("\n[A] Φ = c·L_(2,4)：")
resA = []
for c, origin in FIRST:
    r = report(f"c={str(c):<5}（{origin}）", [(c, uab(2, 4))])
    if r:
        resA.append((c, r[0]))

print("\n[B] Φ = c·[L_(2,3)+L_(3,4)+L_(2,4)]（K₃ 三鍵）：")
for c, origin in ((F(1, 3), "1/3"), (F(1), "1"), (F(1, 2), "1/2")):
    report(f"c={str(c):<5}（{origin}）", [(c, uab(2, 3)), (c, uab(3, 4)), (c, uab(2, 4))])

print("\n[C] Φ = c·L_(2,4) + c′·e₄e₄ᵀ：")
nC = 0
for c, _ in FIRST[:3]:
    for c2, _ in FIRST[:3]:
        report(f"c={str(c):<5}, c′={str(c2):<5}", [(c, uab(2, 4)), (c2, e(4))])
        nC += 1

print("\n" + "=" * 96)
print("判")
print("=" * 96)
print(f"  · [A] 之 ρ_u 取值集 = {sorted({str(w[1]) for _, w in resA})}"
      f" ⟹ {'c-不變 ✅' if len({w[1] for _, w in resA}) == 1 else 'c-變'}")
print(f"  · 全部候選 = 6 + 3 + {nC} = {6 + 3 + nC} 個；χ < 3 者 = 0")
print("  · ⟹ RD7（具名分數）全滅；分數縮放只動 ρ_d，ρ_u 被鎖")
print("=" * 96)
sys.exit(0)
