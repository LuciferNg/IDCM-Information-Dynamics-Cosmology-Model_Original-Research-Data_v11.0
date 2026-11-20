#!/usr/bin/python3
"""
M20 — scripts/06_unit_strength_class.py

⚠️ **SUPERSEDED（2026-11-21）**：本檔有兩處缺陷，**判詞一律作廢**，保留作自捕記錄：
   ① `Gfun` 把 `(G₀U)` 誤寫成 `U`（unit vectors）⟹ Woodbury 修正項全錯（已由 C0 直接求逆暴露）；
   ② 參數化只放對角於站 `i` ⟹ **唔包含 F1**（`+L_(2,4)` 需兩端 +1）。
   正版 = `scripts/06b_unit_strength_complete.py`（C0 回歸含 F1=85/41；類上界 24,596）。

=======================================
R6：**單位強度局部擾動類**之窮盡（MC01 之結構理由測試 ＋ 可達 ρ_d 之禁帶）

【跑前凍結之類定義】（本 docstring 即凍結聲明；以下一字不改）

    𝒰 ≡ { Φ = A·e_i e_iᵀ + B·(e_i e_jᵀ + e_j e_iᵀ)
          : 1 ≤ i < j ≤ 44,  A ∈ {0, +1, −1},  B ∈ {0, +1, −1},  (A,B) ≠ (0,0) }

    · 係數全 ℤ、絕對值 ≤ 1 ⟹ **零參數**（強度 = 鏈自身標度，±1 = 唯一值）
    · 支撐 ≤ 2 站 ⟹ 含：純 on-site（A≠0,B=0）／純鍵（B≠0,A=0）／鍵＋on-site
    · 本類**包含** M20 docs/01 §七 之 F1（A=+1,B=−1,i=2,j=4）及各掃描子類
    · 本類**唔包含**：≥2 條鍵者（= MC01「非最小」）、分數權重者（記為新開放 RD7）

【判準（跑前寫死）】
    C1 窮盡性：成員數 = 946 × 9 − 946（(A,B)=(0,0) 排除）= 7528（另報奇異剔除數）
    C2 聯合：報 min χ；若有成員 |z_d|≤3 ∧ |z_u|≤3 ⟹ 命中
    C3 禁帶：把全部成員之 ρ_d 排序，找目標 0.40350365 所在之空檔（前後最近值）
    C4 對稱性：A=0 子類（純鍵）與 B=0 子類（純 on-site）各自報 min |z_d|

精確性：Fraction 全精確（秩 ≤2 Woodbury，G₀ 閉式）；浮點只用於排序／顯示與 z。
"""

import sys
import math
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
N, H = 44, 45

def g0(i, j):
    return F(min(i, j) * (H - max(i, j)), H)

def Gfun(A, B, i, j):
    """Φ = A e_i e_iᵀ + B (e_i e_jᵀ + e_j e_iᵀ)。
    若 A≠0 且 B≠0：秩 2，U=[e_i, e_j]，C=[[A,B],[B,0]]。
    若 A=0：僅 B 項 ⟹ U=[e_i,e_j]，C=[[0,B],[B,0]]。
    若 B=0：僅 A 項 ⟹ 秩 1。"""
    if A == 0:
        C = [[F(0), F(B)], [F(B), F(0)]]
    else:
        C = [[F(A), F(B)], [F(B), F(0)]]
    try:
        Cinv = [[C[1][1], -C[0][1]], [-C[1][0], C[0][0]]]
        detC = C[0][0] * C[1][1] - C[0][1] * C[1][0]
        if detC == 0:
            return None
        Cinv = [[x / detC for x in row] for row in Cinv]
    except Exception:
        return None
    # ★ 2026-11-21 自捕：修正前版把 (G₀U) 誤寫成 U（unit vectors）⟹ 全表無效。
    #   正解：G = G₀ − (G₀U)(C⁻¹ + UᵀG₀U)⁻¹(UᵀG₀)，故 p 須為 G₀ 之第 i／j 行。
    ci = [g0(k, i) for k in range(1, N + 1)]
    cj = [g0(k, j) for k in range(1, N + 1)]
    A2 = [[g0(i, i), g0(i, j)], [g0(j, i), g0(j, j)]]     # = UᵀG₀U
    Mm = [[Cinv[a][b] + A2[a][b] for b in range(2)] for a in range(2)]
    det = Mm[0][0] * Mm[1][1] - Mm[0][1] * Mm[1][0]
    if det == 0:
        return None
    Mi = [[Mm[1][1] / det, -Mm[0][1] / det], [-Mm[1][0] / det, Mm[0][0] / det]]

    def Gf(x, y):
        p = [ci[x - 1], cj[x - 1]]
        q = [ci[y - 1], cj[y - 1]]
        return g0(x, y) - sum(p[a] * Mi[a][b] * q[b] for a in range(2) for b in range(2))
    return Gf

TGT_D, TGT_U = F(40350365, 10**8), F(-2588837, 10**8)
SD, SU = F(351, 10000), F(3292, 100000)

print("=" * 92)
print("M20 / scripts/06 — 單位強度局部擾動類 𝒰 之窮盡")
print("=" * 92)

# ---------------------------------------------------------------- C0 回歸核（直接求逆對照）
def direct(A, B, i, j):
    """44×44 直接精確求逆（Fraction）—— 與 Woodbury 逐位對照。"""
    K = [[F(0)] * N for _ in range(N)]
    for a in range(N):
        K[a][a] = F(2)
        if a + 1 < N:
            K[a][a + 1] = K[a + 1][a] = F(-1)
    K[i - 1][i - 1] += F(A)
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
                M[r] = [a2 - f * b2 for a2, b2 in zip(M[r], M[c])]
    G = [row[N:] for row in M]
    g = lambda k, m: G[k - 1][m - 1]
    dd = g(3, 3) - g(2, 3)
    ds = g(3, 4) - g(2, 4)
    du = g(3, 2) - g(2, 2)
    if ds == 0:
        return None
    return dd / ds, du / ds

print("\nC0 回歸核（Woodbury vs 直接 44×44 求逆，精確 0 差）：")
_ok = True
for (A, B, i, j) in ((+1, -1, 2, 4), (0, +1, 3, 13), (0, -1, 3, 13), (-1, +1, 3, 12), (+1, +1, 4, 20)):
    w = Gfun(A, B, i, j)
    d = direct(A, B, i, j)
    if w is None or d is None:
        print(f"   A={A:+d} B={B:+d} ({i},{j}) 奇異／退化"); continue
    d_d = w(3, 3) - w(2, 3); d_s = w(3, 4) - w(2, 4); d_u = w(3, 2) - w(2, 2)
    ww = (d_d / d_s, d_u / d_s) if d_s != 0 else None
    same = (ww == d)
    _ok = _ok and same
    print(f"   A={A:+d} B={B:+d} ({i},{j})  Woodbury={ww}  direct={d}  {'✅' if same else '🔴'}")
    if (A, B, i, j) == (+1, -1, 2, 4):
        print(f"      ↳ F1 回歸：ρ_d = {d[0]}（應 = 85/41）、ρ_u = {d[1]}（應 = −2/41）")
print(f"   ⟹ C0 {'PASS' if _ok else 'FAIL'}")
if not _ok:
    print("   ★ 停止：Woodbury 與直接求逆不一致 ⟹ 不得輸出任何判詞。")
    sys.exit(2)

rows = []
skipped = 0
for i in range(1, N + 1):
    for j in range(i + 1, N + 1):
        for A in (0, 1, -1):
            for B in (0, 1, -1):
                if A == 0 and B == 0:
                    continue
                G = Gfun(A, B, i, j)
                if G is None:
                    skipped += 1
                    continue
                d_d = G(3, 3) - G(2, 3)
                d_s = G(3, 4) - G(2, 4)
                d_u = G(3, 2) - G(2, 2)
                if d_s == 0:
                    skipped += 1
                    continue
                rows.append((A, B, i, j, d_d / d_s, d_u / d_s))

print(f"\nC1 窮盡性：成員 = {len(rows)}（奇異／退化剔除 {skipped}；上界 946×9−946 = {946*9-946}）")

zs = [(float((r[4] - TGT_D) / SD), float((r[5] - TGT_U) / SU)) for r in rows]
chi = [math.hypot(a, b) for a, b in zs]
best = min(range(len(rows)), key=lambda k: chi[k])
print(f"\nC2 聯合最佳：χ = {chi[best]:.3f}  z_d = {zs[best][0]:+.2f}  z_u = {zs[best][1]:+.2f}"
      f"  ← A={rows[best][0]:+d} B={rows[best][1]:+d} (i,j)=({rows[best][2]},{rows[best][3]})")
print(f"   ρ_d = {rows[best][4]} = {float(rows[best][4]):.8f}   ρ_u = {rows[best][5]} = {float(rows[best][5]):.8f}")
n1 = sum(1 for a, b in zs if abs(a) < 1 and abs(b) < 1)
n3 = sum(1 for a, b in zs if abs(a) < 3 and abs(b) < 3)
print(f"   |z_d|<1 ∧ |z_u|<1 = {n1} ／ |z_d|<3 ∧ |z_u|<3 = {n3} ／ N = {len(rows)}")

# ---------------------------------------------------------------- C3 禁帶
print("\n" + "-" * 92)
print("C3 可達 ρ_d 譜：目標 0.40350365 所在之空檔")
print("-" * 92)
vals = sorted({r[4] for r in rows})
below = [v for v in vals if v < TGT_D]
above = [v for v in vals if v > TGT_D]
lo, hi = below[-1], above[0]
print(f"   低於目標之最近值 = {lo} = {float(lo):.8f}")
print(f"   高於目標之最近值 = {hi} = {float(hi):.8f}")
print(f"   空檔寬度 = {float(hi) - float(lo):.8f}（目標距下端 {float(TGT_D - lo):.8f}、"
      f"距上端 {float(hi - TGT_D):.8f}）")
band = [(r, z) for r, z in zip(rows, zs) if F(30, 100) <= r[4] <= F(55, 100)]
print(f"   ρ_d ∈ [0.30, 0.55] 之成員 = {len(band)} ／ {len(rows)}")
for r, z in band[:12]:
    print(f"      A={r[0]:+d} B={r[1]:+d} ({r[2]},{r[3]})  ρ_d={str(r[4]):>10}  ρ_u={str(r[5]):>12}  z=({z[0]:+.2f},{z[1]:+.2f})")

# ---------------------------------------------------------------- C4 子類
print("\n" + "-" * 92)
print("C4 子類：純鍵（A=0）／純 on-site（B=0）／混合")
print("-" * 92)
for name, sel in (("純鍵 A=0", lambda r: r[0] == 0),
                  ("純 on-site B=0", lambda r: r[1] == 0),
                  ("混合 A,B≠0", lambda r: r[0] != 0 and r[1] != 0)):
    sub = [(r, z) for r, z in zip(rows, zs) if sel(r)]
    if not sub:
        continue
    k = min(range(len(sub)), key=lambda t: abs(sub[t][1][0]))
    print(f"   {name:<16} n={len(sub):5d}  min|z_d| = {abs(sub[k][1][0]):.3f}"
          f" @ A={sub[k][0][0]:+d} B={sub[k][0][1]:+d} ({sub[k][0][2]},{sub[k][0][3]})"
          f"  ρ_d={sub[k][0][4]}  ρ_u={sub[k][0][5]}")

print("\n" + "=" * 92)
print("判")
print("=" * 92)
print(f"  · 單位強度、支撐 ≤2 之**全部** 7,528 成員：聯合 |z|<1 命中 = {n1}；|z|<3 = {n3}")
print(f"  · ρ_d 譜在目標附近之空檔 = [{float(lo):.6f}, {float(hi):.6f}] ⟹ 目標落於其中")
print("  · ⟹ 分數權重（≠±1）係唯一未試之維度 ⟹ 記為新開放 RD7（須結構性分數，非手調）")
print("=" * 92)
sys.exit(0)
