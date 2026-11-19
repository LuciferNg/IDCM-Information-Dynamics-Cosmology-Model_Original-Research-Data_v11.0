#!/usr/bin/python3
"""
M20 — scripts/01_family_evaluation.py
=====================================
凍結族 F1–F5（見 docs/01 §七，計算前凍結）之精確評估。

零浮點判定：
  · 全部矩陣運算 = Fraction（ℚ 精確）
  · GΦ 以 Woodbury 秩-k 恆等求得，並以 `GΦ·(K+Φ) = I` 逐元精確回驗（偏差必須 = 0）
  · 阻挫指數 = 精確 max-cut（單環圖：環上 2^3 窮舉 ＋ 樹邊全割）
  · ρ/δ = Fraction；z 顯示才用 mpmath

判準見 docs/01 §八（跑前寫死）。
"""

import sys
from fractions import Fraction as F
from itertools import product
import mpmath as mp

mp.mp.dps = 40

N, H = 44, 45
KU, KD, KS = 2, 3, 4          # 已裁站指派（M14／O137→O233）

PASS, FAIL = 0, 0
def chk(name, ok, detail=""):
    global PASS, FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  :: {detail}" if detail else ""))
    if ok: PASS += 1
    else:  FAIL += 1
    return ok

# ------------------------------------------------------------------ 核心矩陣（整數）
def Kmat():
    """Dirichlet Laplacian on P₄₄：對角 2（端點 2），非對角 −1（i,i+1）。"""
    M = [[F(0)] * N for _ in range(N)]
    for i in range(N):
        M[i][i] = F(2)
        if i + 1 < N:
            M[i][i + 1] = M[i + 1][i] = F(-1)
    return M

K = Kmat()

G0 = {(i, j): F(min(i, j) * (H - max(i, j)), H)          # 1-based 索引
      for i in range(1, N + 1) for j in range(1, N + 1)}

# ------------------------------------------------------------------ 線性代數（Fraction）
def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]

def inv(A):
    n = len(A)
    M = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None                      # 奇異
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]

def Gfull(Phi):
    """G = (K+Φ)⁻¹，直接精確求逆（44×44，Fraction）。"""
    A = [[K[i][j] + Phi[i][j] for j in range(N)] for i in range(N)]
    return inv(A), A

def delta(G):
    """δ_x = G(3,k_x) − G(2,k_x)，1-based。"""
    g = lambda i, j: G[i - 1][j - 1]
    return {x: g(3, k) - g(2, k) for x, k in (('u', KU), ('d', KD), ('s', KS))}

# ------------------------------------------------------------------ Φ 族（docs/01 §七）
def Ledge(i, j):
    """邊 Laplacian：(e_i − e_j)(e_i − e_j)ᵀ（1-based）。"""
    M = [[F(0)] * N for _ in range(N)]
    for a, b in ((i, i), (j, j)):
        M[a - 1][b - 1] += F(1)
    for a, b in ((i, j), (j, i)):
        M[a - 1][b - 1] += F(-1)
    return M

def Medge(i, j):
    """純 hopping：e_i e_jᵀ + e_j e_iᵀ（1-based）。"""
    M = [[F(0)] * N for _ in range(N)]
    M[i - 1][j - 1] += F(1)
    M[j - 1][i - 1] += F(1)
    return M

def add(*Ms):
    out = [[F(0)] * N for _ in range(N)]
    for M in Ms:
        for i in range(N):
            for j in range(N):
                out[i][j] += M[i][j]
    return out

def scal(M, c):
    return [[M[i][j] * c for j in range(N)] for i in range(N)]

FAMILY = [
    ("F1", "+L(2,4)",  Ledge(2, 4),                    "加邊（graph-Laplacian）"),
    ("F2", "−L(2,4)",  scal(Ledge(2, 4), F(-1)),       "加邊但 FM 型"),
    ("F3", "−M(2,4)",  scal(Medge(2, 4), F(-1)),       "AF，tight-binding（只改非對角）"),
    ("F4", "+M(2,4)",  Medge(2, 4),                    "FM，tight-binding"),
    ("F5", "L23+L34+L24", add(Ledge(2, 3), Ledge(3, 4), Ledge(2, 4)), "全三角三鍵重加（違 A6，敏感度）"),
]

# ------------------------------------------------------------------ 圖結構
def support_edges(Phi):
    """supp(K+Φ) 之邊集：鏈邊 (i,i+1) 耦合 −1，Φ 之非對角併入同一鍵（權重相加）。
    回傳 dict{(i,j): 總耦合}，i<j。"""
    E = {(i, i + 1): F(-1) for i in range(1, N)}
    for i in range(1, N + 1):
        for j in range(i + 1, N + 1):
            w = Phi[i - 1][j - 1]
            if w != 0:
                E[(i, j)] = E.get((i, j), F(0)) + w
    return E

def is_bipartite(E, n=N):
    adj = {v: [] for v in range(1, n + 1)}
    for (i, j) in E:
        adj[i].append(j); adj[j].append(i)
    col = {}
    for s in range(1, n + 1):
        if s in col: continue
        col[s] = 0; st = [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if w not in col:
                    col[w] = 1 - col[v]; st.append(w)
                elif col[w] == col[v]:
                    return False, None
    return True, col

def max_satisfied_exact(E, n=N):
    """精確「最大滿足數」（符號圖，零近似）。
    符號語義：w<0 = AF（滿足 ⟺ 兩端異色）；w>0 = FM（滿足 ⟺ 兩端同色）。
    阻挫指數 = |E| − max_satisfied（= 最小未滿足鍵數；Harary 平衡 ⟹ 0）。
    方法：枚舉 chord 端點著色（≤2^6）× 沿徑 2-state DP（精確）。"""
    chords = [e for e in E if e[1] != e[0] + 1]
    ends = sorted({v for e in chords for v in e})
    best = -1
    for bits in product([0, 1], repeat=len(ends)):
        fix = dict(zip(ends, bits))
        ok = lambda w, a, b: (a != b) if w < 0 else (a == b)
        dp = {0: 0, 1: 0}
        for c in (0, 1):
            if 1 in fix and fix[1] != c:
                dp[c] = None
        for i in range(2, n + 1):
            nd = {}
            for c in (0, 1):
                if i in fix and fix[i] != c:
                    nd[c] = None; continue
                cand = []
                for p in (0, 1):
                    if dp[p] is None: continue
                    val = dp[p] + (1 if E.get((i - 1, i), F(0)) < 0 and p != c else 0)
                    for (a, b) in chords:
                        if b == i:
                            val += 1 if ok(E[(a, b)], fix[a], c) else 0
                    cand.append(val)
                nd[c] = max(cand) if cand else None
            dp = nd
        m = max(v for v in dp.values() if v is not None)
        best = max(best, m)
    return best

# ------------------------------------------------------------------ 對易子
def commutator(Phi):
    PK = matmul(Phi, K)
    KP = matmul(K, Phi)
    return max(abs(PK[i][j] - KP[i][j]) for i in range(N) for j in range(N))

# ------------------------------------------------------------------ 目標（判準用，非定義用）
T_D, S_D = F(40350365, 10**8), F(351, 10000)
T_U, S_U = F(-2588837, 10**8), F(3292, 100000)
RATIO_T = F(24782923, 10**7)

def zf(obs, tgt, sig):
    """z = (obs − tgt)/σ，以 mpmath 顯示。"""
    return (mp.mpf(obs.numerator) / obs.denominator - mp.mpf(tgt.numerator) / tgt.denominator) \
           / (mp.mpf(sig.numerator) / sig.denominator)

print("=" * 88)
print("M20 / scripts/01 — 凍結族 F1–F5 精確評估（docs/01 §七；判準 §八）")
print("=" * 88)

# ---------------------------------------------------------------- 基線（無 Φ）
G_base, _ = Gfull([[F(0)] * N for _ in range(N)])
db = delta(G_base)
rd_b, ru_b = db['d'] / db['s'], db['u'] / db['s']
print(f"\n基線（Φ=0, K⁻¹）：ρ_d = {rd_b}   ρ_u = {ru_b}   [1/ρ_d = {1/rd_b}]")
chk("基線 = 42/41 : 1 : −2/41（重現閘 §1b 之內嵌核對）",
    (rd_b, ru_b) == (F(42, 41), F(-2, 41)), "已核")

rows = []
for tag, name, Phi, note in FAMILY:
    print("\n" + "-" * 88)
    print(f"{tag}  Φ = {name}   （{note}）")
    print("-" * 88)

    # --- A1 離散
    a1 = all(isinstance(Phi[i][j], F) and Phi[i][j].denominator == 1
             for i in range(N) for j in range(N))
    # --- A2 非遞移
    cm = commutator(Phi)
    a2 = cm != 0
    # --- A3/A4 圖結構
    E = support_edges(Phi)
    bip, _ = is_bipartite(E)
    msat = max_satisfied_exact(E)
    frust = len(E) - msat
    a3 = not bip
    a4 = all(abs(w) == 1 for w in E.values()) and frust == 1
    # --- A5 帶味
    support_sites = sorted({i for i in range(1, N + 1) for j in range(1, N + 1)
                            if Phi[i - 1][j - 1] != 0})
    a5 = set(support_sites) <= {KU, KD, KS}
    # --- A6 附加性（不改既有鍵權重）
    a6 = all(Phi[i - 1][i] == 0 for i in range(1, N))       # 既有鍵 (i,i+1) 之非對角不得被改
    # --- A7 可傳播（正定）
    Gp, A = Gfull(Phi)
    a7 = Gp is not None
    if a7:
        prod = matmul(Gp, A)
        ident = all(prod[i][j] == (F(1) if i == j else F(0)) for i in range(N) for j in range(N))
        # 正定性：G 存在且 G = A⁻¹ ⟹ A 非奇異；再以 LDLᵀ 樞軸全正核 SPD
        def is_spd(M):
            n = len(M); Aa = [row[:] for row in M]
            for c in range(n):
                if Aa[c][c] <= 0: return False
                for r in range(c + 1, n):
                    f = Aa[r][c] / Aa[c][c]
                    if f != 0:
                        for k in range(c, n):
                            Aa[r][k] -= f * Aa[c][k]
            return True
        a7 = ident and is_spd(A)
    else:
        ident = False

    print(f"  A1 離散           : {'✅' if a1 else '🔴'}  係數全 ℤ")
    print(f"  A2 非遞移 [Φ,K]≠0 : {'✅' if a2 else '🔴'}  max|[Φ,K]| = {cm}")
    print(f"  A3 非二部（奇環） : {'✅' if a3 else '🔴'}  邊數 {len(E)}、最大滿足 {msat} ⟹ 阻挫 = {frust}")
    print(f"  A4 ±1 ＋ 阻挫=1   : {'✅' if a4 else '🔴'}")
    print(f"  A5 帶味 supp ⊆ S  : {'✅' if a5 else '🔴'}  supp = {support_sites}")
    print(f"  A6 只加新鍵       : {'✅' if a6 else '🔴'}")
    print(f"  A7 正定／G(K+Φ)=I : {'✅' if a7 else '🔴'}  回驗 = {ident}")
    admitted = all([a1, a2, a3, a4, a5, a6, a7])
    print(f"  ⟹ 准入 = {'✅ YES' if admitted else '🔴 NO'}")

    if Gp is None:
        print("  （K+Φ 奇異 ⟹ 無 Green 核 ⟹ 不計 ρ）")
        rows.append((tag, name, admitted, None, None, None, cm, frust))
        continue

    d = delta(Gp)
    rd, ru = d['d'] / d['s'], d['u'] / d['s']
    zd, zu = zf(rd, T_D, S_D), zf(ru, T_U, S_U)
    inv_rd = 1 / rd
    print(f"  δ_u = {d['u']} = {float(d['u']):+.8f}   （要求 < 0 ⟹ {'✅' if d['u'] < 0 else '🔴'}）")
    print(f"  δ_d = {d['d']} = {float(d['d']):+.8f}")
    print(f"  δ_s = {d['s']} = {float(d['s']):+.8f}")
    print(f"  ρ_d = {rd} = {float(rd):.8f}   z = {mp.nstr(zd,6)}")
    print(f"  ρ_u = {ru} = {float(ru):.8f}   z = {mp.nstr(zu,6)}")
    print(f"  1/ρ_d = {inv_rd} = {float(inv_rd):.8f}   （目標 2.4782923；相對差 "
          f"{mp.nstr(abs(mp.mpf(float(inv_rd))-mp.mpf('2.4782923'))/mp.mpf('2.4782923')*100,5)}%）")
    rows.append((tag, name, admitted, rd, ru, d, cm, frust))

# ---------------------------------------------------------------- 匯總
print("\n" + "=" * 88)
print("匯總表")
print("=" * 88)
print(f"{'#':<4}{'Φ':<16}{'准入':<7}{'ρ_d':>22}{'z_d':>9}{'ρ_u':>22}{'z_u':>9}{'1/ρ_d':>10}")
for tag, name, adm, rd, ru, d, cm, frust in rows:
    if rd is None:
        print(f"{tag:<4}{name:<16}{'NO':<7}{'— 奇異/不定 —':>22}")
        continue
    print(f"{tag:<4}{name:<16}{('YES' if adm else 'no'):<7}"
          f"{float(rd):>22.8f}{mp.nstr(zf(rd,T_D,S_D),4):>9}"
          f"{float(ru):>22.8f}{mp.nstr(zf(ru,T_U,S_U),4):>9}{float(1/rd):>10.5f}")

adm_rows = [r for r in rows if r[2]]
print(f"\n准入族大小 N_trials = {len(adm_rows)} ： {[r[0] for r in adm_rows]}")

# ---------------------------------------------------------------- 判準 C1–C6
print("\n" + "=" * 88)
print("判準（docs/01 §八）")
print("=" * 88)
if len(adm_rows) == 1:
    tag, name, adm, rd, ru, d, cm, frust = adm_rows[0]
    print(f"C1 准入：唯一 = {tag} {name}")
    chk("C2 阻挫指數 = 1", frust == 1, f"阻挫 = {frust}")
    chk("C3 ㈠ 非平移不變（[Φ,K]≠0）", cm != 0, f"max|[Φ,K]| = {cm}")
    rel = abs(mp.mpf(float(1 / rd)) - mp.mpf('2.4782923')) / mp.mpf('2.4782923')
    chk("C3 ㈡ 1/ρ_d 與 2.4782923 同向、相對差 ≤ 10%", rel <= mp.mpf('0.10'),
        f"1/ρ_d = {float(1/rd):.6f}；相對差 = {mp.nstr(rel*100,5)}%")
    chk("C4 T3 符號：δ_u < 0", d['u'] < 0, f"δ_u = {float(d['u']):+.8f}")
    chk("C4 T3 量：|ρ_d − 0.40350365| < 3σ_ρd", abs(zf(rd, T_D, S_D)) < 3,
        f"ρ_d = {float(rd):.8f}；z = {mp.nstr(zf(rd,T_D,S_D),5)}")
    print(f"     ρ_u = {float(ru):.8f}；z = {mp.nstr(zf(ru,T_U,S_U),5)}（報數，非判準）")
    chk("C6 精度：1/ρ_d 分母 ∉ 41ℤ ⟹ ℚ/41 格點界限不適用（預言可判別）",
        (1 / rd).denominator % 41 != 0,
        f"1/ρ_d = {1/rd}（分母 { (1/rd).denominator }；41·{(1/rd).denominator/41 if (1/rd).denominator%41==0 else 0}）")
else:
    print(f"C1 准入族非唯一（{len(adm_rows)}）⟹ 逐個報，不作單一預言")
    chk("C1 准入族唯一性", False, f"N_trials = {len(adm_rows)}")

print("\n" + "=" * 88)
print(f"TOTAL: PASS = {PASS} ｜ FAIL = {FAIL}")
print("=" * 88)
sys.exit(0 if FAIL == 0 else 1)
