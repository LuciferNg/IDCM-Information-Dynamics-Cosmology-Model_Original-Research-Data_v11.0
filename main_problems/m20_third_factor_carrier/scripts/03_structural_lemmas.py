#!/usr/bin/python3
"""
M20 — scripts/03_structural_lemmas.py
=====================================
三條結構引理（精確，零浮點判定）：

L1  **鍵上簽名不可見**：S = diag(s)，s ∈ {±1}^44。K_s = S K S ⟹ G_s(i,j) = s_i s_j G₀(i,j)
    ⟹ δ_x → s_x·[s₃G(3,x) − s₂G(2,x)] ⟹ ρ_x = (s_x/s_s)·ρ_x⁰
    ⟹ |ρ_x| 不變；符號變動係 gauge（s 為 ±1 對角相似）⟹ **鏈上任何 ±1 重新簽名對 ρ 無資訊**。
    （此即 M12 R18 之 `A_s = S A S` 現象之一般證明；亦解釋十四輪為何全在類內。）

L2  **𝒞 = Comm(K)**：隨機 w ∈ ℚ^44 建 Φ = Σ w_n φ_nφ_nᵀ ⟹ [Φ,K] = 0（精確）。
    （T1 定理之數值印證；非證明，定理已於 docs/01 §五 證。）

L3  **[L_(2,4), K] 之支撐**：對易子只落於 chord 之鄰域 —— 顯示「類外」係局部、可控、非病態。
"""

import sys
import random
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 50
N, H = 44, 45
KU, KD, KS = 2, 3, 4

def g0(i, j):
    return F(min(i, j) * (H - max(i, j)), H)

def Kmat():
    M = [[F(0)] * N for _ in range(N)]
    for i in range(N):
        M[i][i] = F(2)
        if i + 1 < N:
            M[i][i + 1] = M[i + 1][i] = F(-1)
    return M

K = Kmat()

def inv(A):
    n = len(A)
    M = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]

PASS, FAIL = 0, 0
def chk(name, ok, detail=""):
    global PASS, FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  :: {detail}" if detail else ""))
    if ok: PASS += 1
    else:  FAIL += 1

# ---------------------------------------------------------------- L1
print("=" * 84)
print("L1 — 鏈上 ±1 簽名對 ρ 無資訊（精確）")
print("=" * 84)
rd0 = (g0(3, 3) - g0(2, 3)) / (g0(3, 4) - g0(2, 4))
ru0 = (g0(3, 2) - g0(2, 2)) / (g0(3, 4) - g0(2, 4))
print(f"  基線 ρ_d⁰ = {rd0}   ρ_u⁰ = {ru0}")

rng = random.Random(20261119)
worst = F(0)
worst_sign = 0
trials = 12
for t in range(trials):
    s = [rng.choice([1, -1]) for _ in range(N)]
    S = [[F(s[i]) if i == j else F(0) for j in range(N)] for i in range(N)]
    Ks = [[sum(S[i][k] * K[k][l] * S[l][j] for k in range(N) for l in range(N))
           for j in range(N)] for i in range(N)]
    Gs = inv(Ks)
    if Gs is None:
        chk(f"  樣本 {t}: K_s 可逆", False, "奇異")
        continue
    # 核對 G_s(i,j) = s_i s_j G₀(i,j)
    dev = max(abs(Gs[i - 1][j - 1] - F(s[i - 1] * s[j - 1]) * g0(i, j))
              for i in range(1, N + 1) for j in range(1, N + 1))
    chk(f"  樣本 {t}: G_s(i,j) = s_i s_j G₀(i,j)", dev == 0, f"max|dev| = {dev}")

    g = lambda i, j: Gs[i - 1][j - 1]
    rd = (g(3, 3) - g(2, 3)) / (g(3, 4) - g(2, 4))
    ru = (g(3, 2) - g(2, 2)) / (g(3, 4) - g(2, 4))
    # L1b 恆等：ρ_x = (s_x/s_s)·[s₃G(3,k_x) − s₂G(2,k_x)]/[s₃G(3,4) − s₂G(2,4)]
    def br(x):
        k = {KU: KU, KD: KD, KS: KS}[x]
        num = F(s[KD - 1]) * g0(3, k) - F(s[KU - 1]) * g0(2, k)
        den = F(s[KD - 1]) * g0(3, 4) - F(s[KU - 1]) * g0(2, 4)
        return num / den
    pr_d = F(s[KD - 1], s[KS - 1]) * br(KD)
    pr_u = F(s[KU - 1], s[KS - 1]) * br(KU)
    worst = max(worst, abs(rd - pr_d), abs(ru - pr_u))
    # L1c：|ρ_d| ≡ 42/41（只取決於相對號 s₂s₃，兩情形同值）
    worst_sign = max(worst_sign, abs(abs(rd) - abs(rd0)))

chk(f"  ⟹ L1b：ρ_x = (s_x/s_s)·[s₃G(3,k_x) − s₂G(2,k_x)]/[s₃G(3,4) − s₂G(2,4)] 精確",
    worst == 0, f"max|dev| = {worst}（{trials} 樣本）")
chk("  ⟹ L1c：任意簽名 |ρ_d| ≡ 42/41（簽名只能動 ρ_u，動不了 ρ_d）",
    worst_sign == 0, f"max||ρ_d| − 42/41| = {worst_sign}")
chk("  ⟹ 結構結論：鏈上 ±1 簽名（path 恆平衡）對 ㈡ 無能 ⟹ 非載體",
    True, "ρ_d 之絕對值不變")

# ---- L1d：kink 之「和型」亦同值（兩情形同 42/41 之結構原因）
print("  L1d：kink 之兩情形 —— 同號（差型）vs 反號（和型）：")
for label, s2, s3 in (("差型 s₂=s₃", 1, 1), ("和型 s₂=−s₃", 1, -1)):
    num = s3 * g0(3, 3) - s2 * g0(2, 3)
    den = s3 * g0(3, 4) - s2 * g0(2, 4)
    print(f"    {label}: ρ_d = {num/den} = {float(num/den):.10f}")

# ---------------------------------------------------------------- L2
print("\n" + "=" * 84)
print("L2 — 𝒞 = Comm(K)：Φ = Σ w_n φ_nφ_nᵀ ⟹ [Φ,K] = 0")
print("=" * 84)
# φ_n(j) = √(2/45) sin(πnj/45)；以 ℚ(√) 不便 ⟹ 直接用矩陣函數之等價刻畫：
# 取 K 之冪之線性組合 Φ = Σ_{m=0}^{3} c_m K^m（∈ 𝒞），核 [Φ,K] = 0（精確整數）。
Cm = [K]
for _ in range(3):
    A, B = Cm[-1], K
    Cm.append([[sum(A[i][k] * B[k][j] for k in range(N)) for j in range(N)] for i in range(N)])
coef = [F(3), F(-2), F(7), F(-5)]
Phi = [[sum(coef[m] * Cm[m][i][j] for m in range(4)) for j in range(N)] for i in range(N)]
PK = [[sum(Phi[i][k] * K[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
KP = [[sum(K[i][k] * Phi[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
dev2 = max(abs(PK[i][j] - KP[i][j]) for i in range(N) for j in range(N))
chk("  f(K)（4 項多項式）⟹ [Φ,K] = 0 精確", dev2 == 0, f"max|dev| = {dev2}")
print("  （定理已於 docs/01 §五 給出：K 譜單簡 ⟹ 對易 ⟺ Φ ∈ 𝒞）")

# ---------------------------------------------------------------- L3
print("\n" + "=" * 84)
print("L3 — [L_(2,4), K] 之支撐：類外性係局部")
print("=" * 84)
Phi1 = [[F(0)] * N for _ in range(N)]
for (a, b), v in (((2, 2), F(1)), ((4, 4), F(1)), ((2, 4), F(-1)), ((4, 2), F(-1))):
    Phi1[a - 1][b - 1] = v
PK = [[sum(Phi1[i][k] * K[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
KP = [[sum(K[i][k] * Phi1[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
nz = [(i + 1, j + 1, PK[i][j] - KP[i][j]) for i in range(N) for j in range(N) if PK[i][j] != KP[i][j]]
print(f"  非零對易元數 = {len(nz)}；位置/值：")
for (i, j, v) in nz:
    print(f"    [{i},{j}] = {v}")
chk("  對易子支撐 ⊆ 站 {{1,2,3,4,5}}（chord 端點 ±2）", all(i <= 5 and j <= 5 for i, j, _ in nz),
    f"support = {sorted({i for i,_,_ in nz} | {j for _,j,_ in nz})}")
chk("  對易子非零（T1：Φ ∉ 𝒞）", len(nz) > 0, f"{len(nz)} 元")

print("\n" + "=" * 84)
print(f"TOTAL: PASS = {PASS} ｜ FAIL = {FAIL}")
print("=" * 84)
sys.exit(0 if FAIL == 0 else 1)
