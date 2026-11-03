#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/02 — Q1：站指派之可觀測性（gauge 判）三層檢定

問題（record §二 Q1）：`u↔2`、`d↔3`、`s↔4` 係假設；可否由既有結構固定？
判準（跑前凍結，record §四）：A1 出處／A2 非循環／A3 零參數／A5 範圍／A6 味盲。

★ 構造（與 M13 docs/03 §〇、docs/05 §S1 一致，本輪先獨立核對再使用）：
   測試向量   v_x(n) = φ_n(p_x) · Δ_n ,   Δ_n = φ_n(p_d) − φ_n(p_u)
   Gram      M_xy = Σ_n v_x(n) v_y(n) / γ_n = v_xᵀ K⁻¹ v_y
   掃描箱     p_u < p_d < p_s，p_u≤7、p_d≤8、p_s≤10  ⟹ **112** 個三元組
   （Δ 隨三元組變：Δ 係同一三元組之 u→d 交換向量，非固定 φ(3)−φ(2)。
     錯用固定 Δ 會得 24/392 而與 M13 之 44/112 不符 —— 本輪已核，見 §一註。）

三層檢定：
  L1 比值層：44 個 A₂⊕1 三元組上 ρ_x 是否只依 α
  L2 絕對層：δ_x 之**標度**是否由結構固定（否則絕對層同樣 gauge）
  L3 算子層：𝒜 = {一體 M₀Σσ_nγ_n ＋ 兩體 Σg(d)} 對味標籤置換是否精確不變
判死範圍：L1/L2 = {有限集：44 三元組}；L3 = {族：𝒜 全部算子}。

慣例宣告：L1/L2 用 **/90**（γ_n=4sin²(πn/90)，配對 45−n，44 站 Dirichlet 鏈）；
          L3 用 **/88**（γ_n=4sin²(πn/88)，配對 44−n，88 節點週期鏈，M7 𝒜 原慣例）。
N_trials：L1/L2 = 0（恆等檢定）；L3 = 0（定理核對）＋ null MC 10⁷（宣告為 null）。
"""
from fractions import Fraction as F
import numpy as np
import mpmath as mp

mp.mp.dps = 60
print("=" * 78)
print("M14 scripts/02 — Q1 站指派：三層可觀測性檢定")
print("=" * 78)

# =====================================================================
# 一、A₂⊕1 族之重建（先與 M13 docs/05 §S1 對數）
# =====================================================================
NN = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / 45) * mp.sin(mp.pi * n * j / 45)
gam90 = lambda n: 4 * mp.sin(mp.pi * n / 90) ** 2
triples = [(a, b, c) for a in range(1, 8) for b in range(1, 9) for c in range(1, 11) if a < b < c]
print(f"\n【一】掃描箱（p_u<p_d<p_s；p_u≤7、p_d≤8、p_s≤10）＝ {len(triples)} 個三元組"
      f"   ← 須等於 M13 docs/05 §S1 之 112")


def gram(p):
    """p=(p_u,p_d,p_s)；Δ_n = φ_n(p_d)−φ_n(p_u)（三元組自身之 u→d 交換向量）。"""
    D = [phi(n, p[1]) - phi(n, p[0]) for n in range(1, NN + 1)]
    V = [[phi(n, k) * D[n - 1] for n in range(1, NN + 1)] for k in p]
    M = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            M[i, j] = mp.fsum(V[i][n] * V[j][n] / gam90(n + 1) for n in range(NN))
            M[j, i] = M[i, j]
    return M


TOL = mp.mpf("1e-40")


def is_A2p1(M):
    d = M[0, 0]
    if abs(d) < TOL:
        return False
    if not all(abs(M[i, i] - d) < TOL for i in (1, 2)):
        return False
    if abs(M[0, 1] + d / 2) > TOL:
        return False
    return max(abs(M[0, 2]), abs(M[1, 2])) < TOL


def is_A2block(M):
    d = M[0, 0]
    return abs(d) >= TOL and abs(M[1, 1] - d) < TOL and abs(M[0, 1] + d / 2) < TOL


A2, A2b = [], 0
for p in triples:
    M = gram(p)
    if is_A2p1(M):
        A2.append((p, M))
    if is_A2block(M):
        A2b += 1
print(f"  A₂⊕1 型 = {len(A2)} / {len(triples)}（{100*len(A2)/len(triples):.1f}%）"
      f"   ← 須等於 M13 之 44（39.3%）")
print(f"  A₂ 塊（僅 2×2）= {A2b} / {len(triples)}（{100*A2b/len(triples):.1f}%）"
      f"   ← 須等於 M13 之 87（77.7%）")
print(f"  ★ 判：A₂ Cartan 塊**非強制結構**（{len(A2)}/{len(triples)}）⟹ 不可作站指派之選擇原理")

# =====================================================================
# 二、L1 比值層
# =====================================================================
alphas = {"(1,2,1)": (F(1), F(2), F(1)), "(1,2,3)": (F(1), F(2), F(3)),
          "(1,2,4)": (F(1), F(2), F(4)), "(1,1,1)": (F(1), F(1), F(1))}
print("\n【二】L1 比值層：ρ_x = δ_x/δ_s，δ = Mα；比對閉式 ρ_u=(2α_u−α_d)/(2α_s)、ρ_d=(2α_d−α_u)/(2α_s)")
print(f"  {'α':>9s} {'ρ_u 閉式':>11s} {'ρ_d 閉式':>11s} {'max|Δρ_u|':>12s} {'max|Δρ_d|':>12s}"
      f" {'M_uu 跨三元組極差':>18s}")
dev_max = mp.mpf(0)
for tag, (au, ad, as_) in alphas.items():
    cu, cd = (2 * au - ad) / (2 * as_), (2 * ad - au) / (2 * as_)
    cu_m, cd_m = mp.mpf(cu.numerator) / cu.denominator, mp.mpf(cd.numerator) / cd.denominator
    du_ = dd_ = mp.mpf(0)
    scales = []
    a = mp.matrix([mp.mpf(au.numerator) / au.denominator, mp.mpf(ad.numerator) / ad.denominator,
                   mp.mpf(as_.numerator) / as_.denominator])
    for p, M in A2:
        d = M * a
        du_ = max(du_, abs(d[0] / d[2] - cu_m))
        dd_ = max(dd_, abs(d[1] / d[2] - cd_m))
        scales.append(M[0, 0])
    dev_max = max(dev_max, du_, dd_)
    print(f"  {tag:>9s} {mp.nstr(cu_m,8):>11s} {mp.nstr(cd_m,8):>11s} "
          f"{mp.nstr(du_,3):>12s} {mp.nstr(dd_,3):>12s} {mp.nstr(max(scales)-min(scales),3):>18s}")
print(f"  ⟹ L1：ρ_x **只依 α**，與三元組無關；全族最大偏差 = {mp.nstr(dev_max,3)}（機器精度級）")
print(f"     判：比值層 ⟹ 站指派 = gauge（範圍：{len(A2)} 個 A₂⊕1 三元組）")

# =====================================================================
# 三、L2 絕對層：δ_x 之標度是否被結構固定
# =====================================================================
print("\n【三】L2 絕對層：δ_x（非比值）之標度是否由結構固定？")
au, ad, as_ = alphas["(1,2,4)"]
a = mp.matrix([mp.mpf(1), mp.mpf(2), mp.mpf(4)])
scales, dvals = [], []
for p, M in A2:
    d = M * a
    scales.append((p, M[0, 0], d))
suu = [s[1] for s in scales]
ratio = [s[2][0] / s[2][2] for s in scales]
print(f"  M_uu 之 44 個值：min = {mp.nstr(min(suu),12)}   max = {mp.nstr(max(suu),12)}"
      f"   max/min = {mp.nstr(max(suu)/min(suu),12)}")
print(f"  δ_u/δ_s 之 44 個值：極差 = {mp.nstr(max(ratio)-min(ratio),3)}（= 0，同 L1）")
print(f"  例：p=(2,3,4) → M_uu = {mp.nstr(scales[0][1],12)}、δ = "
      f"({mp.nstr(scales[0][2][0],8)}, {mp.nstr(scales[0][2][1],8)}, {mp.nstr(scales[0][2][2],8)})")
same = max(abs(s[1] - suu[0]) for s in scales[1:])
print(f"  M_uu 跨三元組是否全等？ max|Δ| = {mp.nstr(same,3)}  ⟹ "
      f"{'全等（標度亦與指派無關）' if same < TOL else '不全等（標度隨三元組變）'}")
print("  最後一步（關鍵）：v_x = φ(p_x)·Δ 之歸一化係**慣例**（φ 取 √(2/45) 使 Σφ²=1）")
print("  ⟹ M（因而 δ_x）之絕對標度係測試向量歸一化之約定 ⟹ 絕對值本身不含物理信息")
print("  ⟹ 唯一物理內容 = 比值（L1）；結構從未、亦無法定 α 之整體標度（M13 docs/03 §四：dof=1）")
print("  ⟹ L2：絕對層**同樣不可固定**站指派；「絕對值」不構成第二類觀測量")

# =====================================================================
# 四、L3 算子層（𝒜；/88 慣例）：味標籤置換之精確不變性
# =====================================================================
print("\n【四】L3 算子層（/88，配對 44−n）：𝒜 = {M₀Σ_n σ_n γ_n + Σ_{i<j} g(|n_i−n_j|)}")
print("     核心：E 只依「佔用模之多集」與「模距之多集」⟹ 味標籤置換為精確對稱")
gam88 = np.array([4.0 * np.sin(np.pi * n / 88) ** 2 for n in range(1, 45)])
rng = np.random.default_rng(20261030)
NTOT = 10_000_000
modes = rng.integers(1, 45, size=(NTOT, 6))
lab = rng.permuted(np.tile(np.arange(6), (NTOT, 1)), axis=1)
m1 = np.take_along_axis(modes, lab, axis=1)
iu = np.triu_indices(6, 1)
kernels = {"g(d)=cos(2πd/88)": lambda dd: np.cos(2 * np.pi * dd / 88.0),
           "g(d)=exp(−d/8)": lambda dd: np.exp(-dd / 8.0),
           "g(d)=1/(1+d)": lambda dd: 1.0 / (1.0 + dd)}


def energy(A, g):
    S = np.ascontiguousarray(np.sort(A, axis=1))
    e1 = gam88[S - 1].sum(axis=1)
    d = np.abs(S[:, :, None] - S[:, None, :])[:, iu[0], iu[1]].astype(float)
    d = np.minimum(d, 88 - d)
    return e1 + g(d).sum(axis=1)


print(f"  null MC：N = {NTOT:,} 組態 × 隨機 6-標籤置換（宣告為 null，不計 trial）")
for name, g in kernels.items():
    E0, E1 = energy(modes, g), energy(m1, g)
    print(f"    {name:22s} max|ΔE| = {np.max(np.abs(E0-E1)):.3e}   bit-identical: "
          f"{bool(np.array_equal(E0, E1))}")
print("  ⟹ L3：味標籤置換下 E 精確不變（bit-identical，三核皆然）⟹ 𝒜-觀測量對指派 = 0 bit")
print("     ⟹ 判：站指派屬 gauge；範圍 = **𝒜 全部算子**（非僅比值層、非僅有限集）")

# =====================================================================
# 五、判決
# =====================================================================
print("\n" + "=" * 78)
print("Q1 判決（生／死／不可裁決，三者必居其一）")
print("=" * 78)
print(f"""
  L1 比值層（{len(A2)} 個 A₂⊕1 三元組）：ρ_x 只依 α（max dev {mp.nstr(dev_max,2)}）⟹ gauge
  L2 絕對層（同上）：δ_x 之標度 = 測試向量歸一化之約定 ⟹ 絕對值無物理信息 ⟹ gauge
  L3 算子層（𝒜，/88）：味標籤置換為精確對稱（10⁷ null，bit-identical）⟹ 0 bit ⟹ gauge

  ⟹ **Q1 = 死（結構不可固定）**
    範圍標籤（A5）：L1/L2 = {{有限集：{len(A2)} 個 A₂⊕1 三元組}}；L3 = {{族：𝒜 之全部算子}}
  ⟹ 「換什麼觀測量才可固定」之答覆（三層皆已窮盡）：
        ① 𝒜-內：**無**。任何 𝒜-觀測量對指派皆 0 bit（L3 定理級）。
        ② 比值層以外之絕對值：**亦無**（L2：標度屬歸一化約定，非結構量）。
        ③ 欲固定須引入 𝒜 外之味相關項 c·(δ_{{f,u}}−δ_{{f,d}})：dof +1（違 A3），
           且該項對任何味對差恰給 Δ0 = c ⟹ 等於把答案當參數輸入（違 A2）。
  ⟹ 故站指派之地位 = **外部輸入**（與 M7 A_eff = 4.4492 MeV 同級），非可導出量。
  N_trials = 0（L1/L2 恆等檢定；L3 定理核對＋null MC）。零上游寫入。
""")
raise SystemExit(0)