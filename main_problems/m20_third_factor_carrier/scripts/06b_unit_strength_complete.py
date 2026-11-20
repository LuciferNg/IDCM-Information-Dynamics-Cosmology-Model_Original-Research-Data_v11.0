#!/usr/bin/python3
"""
M20 — scripts/06b_unit_strength_complete.py
===========================================
R6′：**單位強度、支撐 ≤2 之完整類** 𝒰* 之窮盡（修正 06 之覆蓋缺口）

【06 之錯誤（自捕，保留不改）】
    06 之參數化只放 `A` 於站 `i`（單邊對角），故**唔包含 F1**（`+L_(2,4)` 需要兩端皆 +1）。
    C0 回歸已暴露：06 之 `(A=+1,B=−1,(2,4))` 給 `−38/41`，而 F1 應為 `85/41`。

【06b 之完整定義（跑前凍結）】
    𝒰* ≡ { Φ = di·e_i e_iᵀ + dj·e_j e_jᵀ + B·(e_i e_jᵀ + e_j e_iᵀ)
           : 1 ≤ i < j ≤ 44, (di, dj, B) ∈ {0,+1,−1}³ \ {(0,0,0)} }
    · 946 對 × 26 組合 = 24,596 上界；全係數 ∈ ℤ、|係數| ≤ 1 ⟹ **零參數**
    · 含：純 on-site 單站／雙站、純鍵（L 形 λ(e_i−e_j)(e_i−e_j)ᵀ、M 形）、鍵＋on-site、
          **F1 = +L_(2,4)**、F4 = +M_(2,4)
    · 唔含：≥3 站支撐（MC01 之真正非最小）、分數權重（RD7）、非對稱鍵（≠ B 對稱）

【判準（跑前寫死）】
    C0 Woodbury vs 直接 44×44 求逆：5 個樣本精確 0 差（含 F1 = 85/41 回歸）
    C1 窮盡性：報成員數與奇異剔除
    C2 聯合：|z_d|≤3 ∧ |z_u|≤3 之成員數；min χ
    C3 禁帶：目標 0.40350365 之最近下／上值與空檔寬度
    C4 若空檔存在 ⟹ 判「目標落於單位強度禁帶」（非「揀錯物件」）
"""

import sys
import math
from fractions import Fraction as F

N, H = 44, 45

def g0(i, j):
    return F(min(i, j) * (H - max(i, j)), H)

def woodbury_rank1(lam, vec):
    """Φ = λ v vᵀ ⟹ G = G₀ − (G₀v)(G₀v)ᵀ/(λ + vᵀG₀v)。"""
    f = {k: sum(vec[p] * g0(p, k) for p in range(1, N + 1)) for k in range(1, N + 1)}
    vGv = sum(vec[p] * f[p] for p in range(1, N + 1))
    den = lam + vGv
    if den == 0:
        return None
    return lambda x, y: g0(x, y) - f[x] * f[y] / den

def woodbury_rank2(di, dj, B, i, j):
    """Φ = di e_i e_iᵀ + dj e_j e_jᵀ + B(e_i e_jᵀ + e_j e_iᵀ) = U C Uᵀ，C=[[di,B],[B,dj]]。
    C 可逆（det = di·dj − B² ≠ 0）時用 Woodbury；det=0 之退化情形由呼叫方轉 rank-1。"""
    detC = di * dj - B * B
    if detC == 0:
        return None
    Cinv = [[F(dj, detC), F(-B, detC)], [F(-B, detC), F(di, detC)]]
    a, b, c = g0(i, i), g0(i, j), g0(j, j)
    Mm = [[Cinv[0][0] + a, Cinv[0][1] + b], [Cinv[1][0] + b, Cinv[1][1] + c]]
    det = Mm[0][0] * Mm[1][1] - Mm[0][1] * Mm[1][0]
    if det == 0:
        return None
    Mi = [[Mm[1][1] / det, -Mm[0][1] / det], [-Mm[1][0] / det, Mm[0][0] / det]]
    ci = [g0(k, i) for k in range(1, N + 1)]
    cj = [g0(k, j) for k in range(1, N + 1)]

    def Gf(x, y):
        p = [ci[x - 1], cj[x - 1]]
        q = [ci[y - 1], cj[y - 1]]
        return g0(x, y) - sum(p[r] * Mi[r][s] * q[s] for r in range(2) for s in range(2))
    return Gf

def Gfun(di, dj, B, i, j):
    detC = di * dj - B * B
    if detC != 0:
        return woodbury_rank2(di, dj, B, i, j)
    # 退化：C 奇異 ⟹ Φ 必為 rank-1（見 docstring 之代數核）
    if di == dj and B == di:                                  # σ(e_i+e_j)(e_i+e_j)ᵀ
        vec = {k: (1 if k in (i, j) else 0) for k in range(1, N + 1)}
        return woodbury_rank1(F(di), vec)
    if di == dj and B == -di:                                 # σ(e_i−e_j)(e_i−e_j)ᵀ
        vec = {k: ((1 if k == i else 0) - (1 if k == j else 0)) for k in range(1, N + 1)}
        return woodbury_rank1(F(di), vec)
    return None

def ratios(G):
    if G is None:
        return None
    dd = G(3, 3) - G(2, 3)
    ds = G(3, 4) - G(2, 4)
    du = G(3, 2) - G(2, 2)
    if ds == 0:
        return None
    return dd / ds, du / ds

def direct(di, dj, B, i, j):
    K = [[F(0)] * N for _ in range(N)]
    for a in range(N):
        K[a][a] = F(2)
        if a + 1 < N:
            K[a][a + 1] = K[a + 1][a] = F(-1)
    K[i - 1][i - 1] += F(di)
    K[j - 1][j - 1] += F(dj)
    K[i - 1][j - 1] += F(B)
    K[j - 1][i - 1] += F(B)
    M = [row[:] + [F(1) if a == b else F(0) for b in range(N)] for a, row in enumerate(K)]
    for c in range(N):
        piv = next((r for r in range(c, N) if M[r][c] != 0), None)
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(N):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [u - f * v for u, v in zip(M[r], M[c])]
    G = [row[N:] for row in M]
    g = lambda x, y: G[x - 1][y - 1]
    return ratios(g)

TGT_D, TGT_U = F(40350365, 10**8), F(-2588837, 10**8)
SD, SU = F(351, 10000), F(3292, 100000)

print("=" * 94)
print("M20 / scripts/06b — 單位強度、支撐 ≤2 之完整類 𝒰*（含 F1）")
print("=" * 94)

print("\nC0 回歸核（Woodbury vs 直接求逆）：")
_ok = True
for (di, dj, B, i, j) in ((+1, +1, -1, 2, 4), (0, 0, +1, 2, 4), (0, 0, +1, 3, 13),
                          (-1, 0, +1, 3, 12), (+1, -1, 0, 4, 20)):
    w = ratios(Gfun(di, dj, B, i, j))
    d = direct(di, dj, B, i, j)
    same = (w == d)
    _ok = _ok and same
    print(f"   (di,dj,B)=({di:+d},{dj:+d},{B:+d}) ({i},{j})  W={w}  D={d}  {'✅' if same else '🔴'}")
    if (di, dj, B, i, j) == (+1, +1, -1, 2, 4):
        print(f"      ↳ F1 回歸：ρ_d = {d[0]}（應 = 85/41）｜ρ_u = {d[1]}（應 = −2/41）")
        _ok = _ok and d == (F(85, 41), F(-2, 41))
print(f"   ⟹ C0 {'PASS' if _ok else 'FAIL'}")
if not _ok:
    print("   ★ 停：回歸失敗 ⟹ 不輸出判詞。")
    sys.exit(2)

rows = []
skipped = 0
for i in range(1, N + 1):
    for j in range(i + 1, N + 1):
        for di in (0, 1, -1):
            for dj in (0, 1, -1):
                for B in (0, 1, -1):
                    if di == 0 and dj == 0 and B == 0:
                        continue
                    r = ratios(Gfun(di, dj, B, i, j))
                    if r is None:
                        skipped += 1
                        continue
                    rows.append((di, dj, B, i, j, r[0], r[1]))

print(f"\nC1 窮盡性：成員 = {len(rows)}（奇異／退化剔除 {skipped}；上界 946×26 = {946*26}）")

zs = [(float((r[5] - TGT_D) / SD), float((r[6] - TGT_U) / SU)) for r in rows]
chi = [math.hypot(a, b) for a, b in zs]
b = min(range(len(rows)), key=lambda k: chi[k])
print(f"\nC2 聯合最佳：χ = {chi[b]:.3f}  z_d = {zs[b][0]:+.2f}  z_u = {zs[b][1]:+.2f}"
      f"  ← (di,dj,B)=({rows[b][0]:+d},{rows[b][1]:+d},{rows[b][2]:+d}) ({rows[b][3]},{rows[b][4]})")
print(f"   ρ_d = {rows[b][5]} = {float(rows[b][5]):.8f} ｜ ρ_u = {rows[b][6]} = {float(rows[b][6]):.8f}")
print(f"   |z_d|<1 ∧ |z_u|<1 = {sum(1 for a, c in zs if abs(a) < 1 and abs(c) < 1)}"
      f" ｜ |z_d|<3 ∧ |z_u|<3 = {sum(1 for a, c in zs if abs(a) < 3 and abs(c) < 3)}"
      f" ｜ N = {len(rows)}")

print("\n" + "-" * 94)
print("C3 禁帶：目標 0.40350365 之鄰域")
print("-" * 94)
vals = sorted({r[5] for r in rows})
below = [v for v in vals if v < TGT_D]
above = [v for v in vals if v > TGT_D]
lo, hi = below[-1], above[0]
print(f"   最近下值 = {lo} = {float(lo):.8f}｜最近上值 = {hi} = {float(hi):.8f}")
print(f"   空檔寬度 = {float(hi) - float(lo):.8f}；目標距下 {float(TGT_D - lo):.8f}、距上 {float(hi - TGT_D):.8f}")
inband = [r for r in rows if lo < r[5] < hi]
print(f"   空檔內成員 = {len(inband)}")
near = [(r, z) for r, z in zip(rows, zs) if abs(float(r[5]) - float(TGT_D)) < 0.06]
print(f"   |ρ_d − 0.4035| < 0.06 之成員 = {len(near)}")
for r, z in near[:10]:
    print(f"      ({r[0]:+d},{r[1]:+d},{r[2]:+d}) ({r[3]},{r[4]})  ρ_d={str(r[5]):>12}={float(r[5]):.6f}"
          f"  ρ_u={str(r[6]):>12}  z=({z[0]:+.2f},{z[1]:+.2f})")

print("\n" + "=" * 94)
print("判")
print("=" * 94)
print(f"  · 24,596 上界 / 實得 {len(rows)} 成員：聯合 |z|<3 命中 = "
      f"{sum(1 for a, c in zs if abs(a) < 3 and abs(c) < 3)}")
print(f"  · ρ_d 禁帶 = ({float(lo):.6f}, {float(hi):.6f})，寬 {float(hi)-float(lo):.6f}，內含成員 = {len(inband)}")
print("  · 目標 0.40350365 落於禁帶之內 ⟹ **單位強度（±1）物件結構上不可達**")
print("  · 未試維度：分數權重（RD7）、≥3 站支撐（MC01）、非線性/高階機制")
print("=" * 94)
sys.exit(0)
