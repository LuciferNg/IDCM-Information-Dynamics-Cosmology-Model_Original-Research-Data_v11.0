#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/00 — 重現閘（SESSION_START §二 六項；跑 R1 之前必須 EXIT=0）

性質：零 import 他 M 腳本／零讀他 M 產物。全部自行由定義重建。
算術：整數 ＋ Fraction 精確；三角量用 mpmath dps=50（只作向量／Gram 幾何，不作判據主體）。

閘①  零 import 核 T1 定理之載體前提：K 譜單簡（γ_n = 4sin²(πn/90)，n=1..44 兩兩相異）
      ＋ 對易子判據 [Φ,K] ≠ 0 ⟺ Φ ∉ 𝒞（只核陳述與代數，不重跑 M20 scripts）
閘②  算術核 M20 F1 三數：δ_u = −2/131、δ_d = 85/131、δ_s = 41/131 ⟹ ρ_d = 85/41、ρ_u = −2/41
閘③  抽 3 個 ±1 簽名手核 L1c：|ρ_d| ≡ 42/41（差型與和型兩情形皆 42/41）
閘④  核 M12 R24 §五 之 Gram 幾何：cos(v_u,v_d) = −4/5、cos(v_d,v_s) = −1/√5、
      cos(v_u,v_s) = 1/(2√5)、‖v_u‖=‖v_d‖=0.23570226、‖v_s‖=0.21081851
閘⑤  核精度界限：σ_ρd ÷ (1/41) = 1.4391
閘⑥  唔准重跑 M20／M12／M13／M14 腳本（會寫入他 M 目錄）——本目錄 scripts/ 只允許本檔

判據全部跑前寫死：任何一項 FAIL ⟹ EXIT=1，不准開始 R1。
"""
from fractions import Fraction as F
import os
import sys

import mpmath as mp

mp.mp.dps = 50

N = 44                      # 站數（K_chain Dirichlet，P₄₄）
M = N + 1                   # 45（Dirichlet 邊界分母）
K_REF = 1                   # 44 站

TOL_MP = mp.mpf('1e-40')
TOL_EXACT = F(0)

results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))
    return bool(ok)


def chkf(name, value, target, tol, detail=""):
    """mp 值比較（顯示用 tol 由文檔有效位定）。"""
    ok = abs(value - target) <= tol
    return chk(name, ok, detail or "val=%s target=%s" % (mp.nstr(value, 20), mp.nstr(target, 20)))


# ─────────────────────────────────────────────────────────────────────
# 基礎物件：K = P₄₄ 之 Dirichlet Laplacian（整數矩陣）
# ─────────────────────────────────────────────────────────────────────
def build_K(n=N):
    A = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = F(2)
        if i + 1 < n:
            A[i][i + 1] = F(-1)
            A[i + 1][i] = F(-1)
    return A


def inv_exact(A):
    """Fraction Gauss–Jordan；回傳 (inverse, det)。n=44 規模足夠快。"""
    n = len(A)
    a = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(A)]
    det = F(1)
    for c in range(n):
        p = None
        for r in range(c, n):
            if a[r][c] != 0:
                p = r
                break
        if p is None:
            return None, F(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            det = -det
        piv = a[c][c]
        det *= piv
        inv_piv = F(1) / piv
        a[c] = [x * inv_piv for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    inv = [row[n:] for row in a]
    return inv, det


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def commutator(A, B):
    n = len(A)
    AB, BA = matmul(A, B), matmul(B, A)
    return [[AB[i][j] - BA[i][j] for j in range(n)] for i in range(n)]


def edge_lap(i, j, n=N):
    """L_(i,j) = (e_i − e_j)(e_i − e_j)ᵀ，站索引 1-based。"""
    A = [[F(0)] * n for _ in range(n)]
    for s in (i, j):
        A[s - 1][s - 1] += F(1)
    A[i - 1][j - 1] -= F(1)
    A[j - 1][i - 1] -= F(1)
    return A


K = build_K()
G0, detK = inv_exact(K)

print("=== M23 scripts/00 — 重現閘（SESSION_START §二） ===")
print("N = %d 站；det(K) = %s" % (N, detK))

# ── 基線：K G₀ = I 逐元精確 ─────────────────────────────────────────
KG = matmul(K, G0)
dev = max(abs(KG[i][j] - (F(1) if i == j else F(0))) for i in range(N) for j in range(N))
chk("G0.1   K·G₀ = I 逐元精確（1936 元）", dev == 0, "max|dev| = %s" % dev)

# 閉式 G₀(i,j) = min(i,j)(45−max(i,j))/45
dev_closed = F(0)
for i in range(1, N + 1):
    for j in range(1, N + 1):
        v = F(min(i, j) * (M - max(i, j)), M)
        dev_closed = max(dev_closed, abs(G0[i - 1][j - 1] - v))
chk("G0.2   閉式 G₀(i,j)=min(i,j)(45−max(i,j))/45 逐位", dev_closed == 0,
    "max|dev| = %s" % dev_closed)
print("       G₀(2,2)=%s  G₀(3,2)=%s  G₀(3,3)=%s  G₀(2,3)=%s  G₀(3,4)=%s  G₀(2,4)=%s"
      % (G0[1][1], G0[2][1], G0[2][2], G0[1][2], G0[2][3], G0[1][3]))

# ── 閘①：譜單簡 ＋ 對易子判據 ──────────────────────────────────────
# 特徵對 (γ_n, φ_n)：（K φ_n = γ_n φ_n）
def phi(n, j):
    return mp.sqrt(mp.mpf(2) / M) * mp.sin(mp.pi * n * j / M)


gam = [4 * mp.sin(mp.pi * n / (2 * M)) ** 2 for n in range(1, N + 1)]   # 4sin²(πn/90)
dev_eig = mp.mpf(0)
for n in range(1, N + 1):
    for j in range(1, N + 1):
        lhs = 2 * phi(n, j) - (phi(n, j - 1) if j > 1 else mp.mpf(0)) - (phi(n, j + 1) if j < N else mp.mpf(0))
        dev_eig = max(dev_eig, abs(lhs - gam[n - 1] * phi(n, j)))
chkf("G1.1   γ_n = 4sin²(πn/90) 滿足 K φ_n = γ_n φ_n（44×44）", dev_eig, mp.mpf(0), TOL_MP,
     "max|dev| = %s" % mp.nstr(dev_eig, 6))

gaps = [(abs(gam[n - 1] - gam[m - 1]), n, m) for n in range(1, N + 1) for m in range(n + 1, N + 1)]
mingap = min(gaps)
chk("G1.2   譜單簡：44 個 γ 兩兩相異", mingap[0] > mp.mpf('1e-30'),
    "min |γ_n−γ_m| = %s (n=%d,m=%d)" % (mp.nstr(mingap[0], 10), mingap[1], mingap[2]))
# 代數理由（唔靠浮點）：sin²a=sin²b ⟺ a=±b mod π；4sin²(πn/90)=4sin²(πm/90) ⟹ m=n 或 m=90−n
# 90−n ≥ 46 > 44 ⟹ 44 站內只可能 n=m。以整數模算術核：
chk("G1.3   代數唯一性：m ≡ ±n (mod 90) ∧ 1≤m,n≤44 ⟹ m=n",
    all((m % 90 in (n % 90, (-n) % 90)) == (m == n) for n in range(1, N + 1) for m in range(1, N + 1)
        if (m % 90 in (n % 90, (-n) % 90)) or m == n),
    "反例數 = 0（90−n ≥ 46 > 44）")

# [Φ,K] 判據：Φ = L_(2,4)（整數矩陣）⟹ 對易子整數，精確可算
L24 = edge_lap(2, 4)
C24 = commutator(L24, K)
nz = [(i + 1, j + 1, C24[i][j]) for i in range(N) for j in range(N) if C24[i][j] != 0]
chk("G1.4   [L_(2,4), K] ≠ 0（整數精確；非零元計數）", len(nz) > 0,
    "nonzero = %d，站支撐 = %s" % (len(nz), sorted({i for i, j, v in nz} | {j for i, j, v in nz})))

# [f(K), K] = 0：Φ 取整係數多項式 f(K) = K³ − 2K（整數矩陣，精確）
K2 = matmul(K, K)
K3 = matmul(K2, K)
fK = [[K3[i][j] - 2 * K[i][j] for j in range(N)] for i in range(N)]
C0 = commutator(fK, K)
chk("G1.5   [f(K), K] = 0（f 整係數；𝒞 ⊆ Comm(K) 之陳述核）",
    all(C0[i][j] == 0 for i in range(N) for j in range(N)),
    "max|dev| = 0；例 f(x)=x³−2x")

# ── 閘②：M20 F1 三數（Φ = +L_(2,4)）─────────────────────────────────
KPHI = [[K[i][j] + L24[i][j] for j in range(N)] for i in range(N)]
GPHI, detKPHI = inv_exact(KPHI)
chk("G2.1   K+Φ 可逆（det ≠ 0）", detKPHI != 0, "det(K+Φ) = %s" % detKPHI)


def dlt(G, kx):
    """δ_x = G(3,k_x) − G(2,k_x)（1-based 站索引）"""
    return G[2][kx - 1] - G[1][kx - 1]


du, dd, dsn = dlt(GPHI, 2), dlt(GPHI, 3), dlt(GPHI, 4)
print("       δ_u = %s = %s" % (du, float(du)))
print("       δ_d = %s = %s" % (dd, float(dd)))
print("       δ_s = %s = %s" % (dsn, float(dsn)))
chk("G2.2   δ_u = −2/131", du == F(-2, 131), "實值 %s" % du)
chk("G2.3   δ_d = 85/131", dd == F(85, 131), "實值 %s" % dd)
chk("G2.4   δ_s = 41/131", dsn == F(41, 131), "實值 %s" % dsn)
rho_d_f1 = dd / dsn
rho_u_f1 = du / dsn
chk("G2.5   ρ_d = 85/41 精確", rho_d_f1 == F(85, 41), "實值 %s = %s" % (rho_d_f1, float(rho_d_f1)))
chk("G2.6   ρ_u = −2/41 精確", rho_u_f1 == F(-2, 41), "實值 %s = %s" % (rho_u_f1, float(rho_u_f1)))

# ── 閘③：±1 簽名手核 L1c ────────────────────────────────────────────
# G_s(i,j) = s_i s_j G₀(i,j)（K_s = S K S ⟹ G_s = S G₀ S），全部 Fraction 精確。
def signed_rho(s):
    Gs = lambda i, j: s[i - 1] * s[j - 1] * G0[i - 1][j - 1]
    r = {}
    for x, kx in (("u", 2), ("d", 3), ("s", 4)):
        r[x] = Gs(3, kx) - Gs(2, kx)
    return r["d"] / r["s"], r["u"] / r["s"]


cases = {
    "C1 全 +1（差型 s₂=s₃=s₄=+1）": [F(1)] * N,
    "C2 s₃=s₂=+1, s₄=−1（差型）": [F(1)] * N,
    "C3 s₃=+1, s₂=−1, s₄=+1（和型）": [F(1)] * N,
}
cases["C2 s₃=s₂=+1, s₄=−1（差型）"][3] = F(-1)
cases["C3 s₃=+1, s₂=−1, s₄=+1（和型）"][1] = F(-1)
cases["C4 隨機 ±1（seed 固定，11 站翻轉）"] = [F(-1) if i in (5, 7, 9, 13, 21, 22, 30, 31, 38, 41, 44) else F(1)
                                              for i in range(1, N + 1)]

for name, s in cases.items():
    rd, ru = signed_rho(s)
    chk("G3  L1c |ρ_d| ≡ 42/41 — " + name, abs(rd) == F(42, 41),
        "ρ_d = %s = %s ｜ ρ_u = %s = %s" % (rd, float(rd), ru, float(ru)))

# 全 44 個單站簽名（手核擴充，非掃描：只作封閉性）
all_single = True
for j in range(1, N + 1):
    s = [F(1)] * N
    s[j - 1] = F(-1)
    all_single &= (abs(signed_rho(s)[0]) == F(42, 41))
chk("G3.5 全 44 個單站 ±1 簽名 |ρ_d| = 42/41（44/44）", all_single)

# ── 閘④：Gram 幾何（M12 R24 §五）────────────────────────────────────
# v_x(n) = φ_n(k_x)·Δ_n，Δ_n = φ_n(3) − φ_n(2)
D = [phi(n, 3) - phi(n, 2) for n in range(1, N + 1)]


def vx(kx):
    return [phi(n, kx) * D[n - 1] for n in range(1, N + 1)]


def ip(a, b):
    return sum(x * y for x, y in zip(a, b))


vu, vd, vs = vx(2), vx(3), vx(4)
nu, nd, ns = mp.sqrt(ip(vu, vu)), mp.sqrt(ip(vd, vd)), mp.sqrt(ip(vs, vs))
cos_ud, cos_ds, cos_us = ip(vu, vd) / (nu * nd), ip(vd, vs) / (nd * ns), ip(vu, vs) / (nu * ns)

chkf("G4.1  ‖v_u‖ = ‖v_d‖ = √2/6 = 0.23570226", nu, mp.sqrt(mp.mpf(2)) / 6, mp.mpf('1e-9'),
     "‖v_u‖=%s ‖v_d‖=%s" % (mp.nstr(nu, 12), mp.nstr(nd, 12)))
chkf("G4.2  ‖v_s‖ = (√2/6)(2/√5) = 0.21081851", ns,
     mp.sqrt(mp.mpf(2)) / 6 * 2 / mp.sqrt(5), mp.mpf('1e-9'), "‖v_s‖ = %s" % mp.nstr(ns, 12))
chkf("G4.3  cos(v_u,v_d) = −4/5", cos_ud, mp.mpf(-4) / 5, mp.mpf('1e-30'), "cos = %s" % mp.nstr(cos_ud, 25))
chkf("G4.4  cos(v_d,v_s) = −1/√5", cos_ds, -1 / mp.sqrt(5), mp.mpf('1e-30'), "cos = %s" % mp.nstr(cos_ds, 25))
chkf("G4.5  cos(v_u,v_s) = 1/(2√5)", cos_us, 1 / (2 * mp.sqrt(5)), mp.mpf('1e-30'),
     "cos = %s" % mp.nstr(cos_us, 25))

# Gram 譜
Gm = [[ip(vu, vu), ip(vu, vd), ip(vu, vs)],
      [ip(vd, vu), ip(vd, vd), ip(vd, vs)],
      [ip(vs, vu), ip(vs, vd), ip(vs, vs)]]
ev = sorted(mp.eigsy(mp.matrix(Gm), eigvals_only=True))
spec_ref = [mp.mpf('0.00900666'), mp.mpf('0.0378208'), mp.mpf('0.108728')]
chkf("G4.6  Gram 譜 = (0.00900666, 0.0378208, 0.108728)",
     max(abs(a - b) for a, b in zip(ev, spec_ref)), mp.mpf(0), mp.mpf('5e-7'),   # 文檔最少 6 位小數 ⟹ 半個末位 5e-7
     "譜 = %s" % [mp.nstr(x, 9) for x in ev])

# ── 閘⑤：精度界限（M12 R25）────────────────────────────────────────
sig_rd = mp.mpf('0.0351')
sig_ru = mp.mpf('0.03292')
grid = mp.mpf(1) / 41
chkf("G5.1  σ_ρd ÷ (1/41) = 1.4391", sig_rd / grid, mp.mpf('1.4391'), mp.mpf('5e-5'),
     "= %s" % mp.nstr(sig_rd / grid, 8))
chkf("G5.2  σ_ρu = 1.34972/41 = 0.03292", sig_ru, mp.mpf('1.34972') / 41, mp.mpf('1e-6'),
     "= %s" % mp.nstr(sig_ru, 8))
chk("G5.3  半格 0.0122 < σ_ρd ∧ < σ_ρu（落格點 = 零資訊）",
    grid / 2 < sig_rd and grid / 2 < sig_ru,
    "半格 = %s" % mp.nstr(grid / 2, 8))

# ── 閘⑥：唔准重跑他 M 腳本 ──────────────────────────────────────────
here = os.path.dirname(os.path.abspath(__file__))
pyfiles = sorted(f for f in os.listdir(here) if f.endswith(".py"))
chk("G6.1  本目錄 scripts/ 只有本檔（零他 M 腳本、零 subprocess、零 import 模組）",
    pyfiles == ["00_reproduction_gate.py"] or set(pyfiles) <= {"00_reproduction_gate.py"},
    "pyfiles = %s" % pyfiles)
chk("G6.2  本檔零 import 專案模組（只有 stdlib ＋ mpmath）",
    all(m in ("fractions", "os", "sys", "mpmath") for m in
        [x.split(" ")[1].split(".")[0] for x in open(os.path.abspath(__file__)).read().splitlines()
         if x.startswith("import ") or x.startswith("from ")]),
    "imports = %s" % [x for x in open(os.path.abspath(__file__)).read().splitlines()
                      if x.startswith("import ") or x.startswith("from ")])

# ── RESULT ─────────────────────────────────────────────────────────
npass = sum(1 for _, ok in results if ok)
ntot = len(results)
print("")
print("──────────────────────────────────────────")
print("結果: %d/%d PASS" % (npass, ntot))
if npass != ntot:
    print("FAIL 清單: %s" % [n for n, ok in results if not ok])
    print("──────────────────────────────────────────")
    sys.exit(1)
print("閘結論：EXIT=0 ⟹ 可開始 R1（讀出式語法凍結，零數值）")
print("──────────────────────────────────────────")
