#!/usr/bin/env python3
"""37 — R15c：entry-positive 之類別檢定（LP maximin，取代 35 之三支撐限制）

35 之限：Carathéodory 對「圓錐歸屬」完備，但 entry-positivity **唔**在削支撐下保持 ⟹
「387 個可行 triple 中 0 個 entry-positive」**不足以**判定類別不可行。本輪用 LP 正面判定。

預註冊
--------------------------------------------------
族       G(i,j) = Σ_n w_n φ_n(i)φ_n(j)（44 維 w）
約束     ①δ_s = 1、δ_d = ρ_d、δ_u = ρ_u（三等式，精確重現數據比值）
         ②G(i,j) ≥ t ∀ i,j（1936 不等式）
目標     max t
兩類     ㈠ w ≥ 0（PSD ∩ entry-positive，即 R8 之實際類別）
         ㈡ w 自由（entry-positive 但容許負譜權）
F1  t* > 0 ⟹ entry-positive 類**可行** ⟹ R8 之範疇判死對『比值約束』無效（須更正）
F2  t* ≤ 0 ⟹ entry-positive 係真類別約束 ⟹ R8 成立（首個類別級 theorem）
F3  以 mpmath dps=40 覆核 t* 之見證；並在 ρ_d、ρ_u 之 ±1σ 角點重測
"""
import mpmath as mp
import numpy as np
from scipy.optimize import linprog

mp.mp.dps = 40
NS = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
RU, RD, SDU, SDD = 0.02588837, 0.40350365, 0.03292, 0.03510

PH = np.array([[float(phi(n, j)) for j in range(1, NS + 1)] for n in range(1, NS + 1)])  # (n,j)
C = np.zeros((3, NS))
for i, j in enumerate(range(1, NS + 1)):
    C[0, i] = PH[i, 3] - PH[i, 2]      # δ_d（x=3，減 k_u=2）→ 用 φ_n(3)−φ_n(2)
C = np.zeros((3, NS))
for a, x in enumerate((4, 3, 2)):      # rows: δ_s, δ_d, δ_u
    for i in range(NS):
        C[a, i] = PH[i, x - 1] * (PH[i, 2] - PH[i, 1])   # φ_n(x)·(φ_n(3)−φ_n(2))
TGT = np.array([1.0, RD, -RU])

print("=" * 96)
print("37 — R15c：entry-positive 類別檢定（LP maximin min_ij G(i,j)）")
print("=" * 96)


def maximin(ru, rd, wfree, eps=1e-9):
    tgt = np.array([1.0, rd, -ru])
    # vars: w (44) , t (1)
    c = np.zeros(NS + 1); c[-1] = -1.0
    A_eq = np.hstack([C, np.zeros((3, 1))])
    b_eq = tgt
    # G(i,j) = Σ_n w_n φ_n(i)φ_n(j) ≥ t  ⟹  −Σ w φφ + t ≤ 0
    Gmat = np.array([np.outer(PH[:, i], PH[:, j]).diagonal().reshape(-1) for i, j in
                     [(i, j) for i in range(NS) for j in range(NS)]]) if False else None
    # 直接建 (1936, 44)：row (i,j) = [φ_n(i)φ_n(j)]_n
    cols = [(PH[:, i] * PH[:, j]) for i in range(NS) for j in range(NS)]
    A = np.array(cols)                                    # (1936, 44)
    A_ub = np.hstack([-A, np.ones((len(A), 1))])
    b_ub = np.zeros(len(A))
    bounds = [(0, None)] * NS + [(None, None)] if wfree else [(0, None)] * (NS + 1)
    if wfree:
        bounds = [(None, None)] * NS + [(None, None)]
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    return r


for tag, wfree in (("㈠ w ≥ 0（PSD ∩ entry-pos）", False), ("㈡ w 自由（entry-pos のみ）", True)):
    r = maximin(RU, RD, wfree)
    print("\n%s" % tag)
    if not r.success:
        print("   LP 不可行：%s ⟹ 類別空缺（F2 方向）" % r.message[:60])
        continue
    w = r.x[:NS]
    t = r.x[NS]
    G = PH.T @ w if False else (w @ PH) if False else np.einsum("n,ni,nj->ij", w, PH, PH)
    print("   t* = min_ij G(i,j) = %.9g ／ max G = %.9g" % (t, G.max()))
    sup = [n for n in range(1, NS + 1) if abs(w[n - 1]) > 1e-12]
    print("   w 非零個數 = %d ／ 44；w 之值域 [%.4g, %.4g]" % (len(sup), w.min(), w.max()))
    print("   ★ 支撐模式 = %s" % sup)
    print("     w 值 = %s" % ["%.6g" % w[n - 1] for n in sup])
    if t > 0:
        print("   ⟹ **F1：entry-positive 可行**（t* > 0）")
    else:
        print("   ⟹ **F2：entry-positive 不可行**（t* ≤ 0）")
    # mpmath 覆核
    wm = [mp.mpf(float(x)) for x in w]
    g = lambda i, j: sum(wm[n - 1] * phi(n, i) * phi(n, j) for n in range(1, NS + 1))
    ds = g(3, 4) - g(2, 4)
    print("   覆核（dps=40）：δ_u/δ_s = %s、δ_d/δ_s = %s、min G = %s"
          % (mp.nstr((g(3, 2) - g(2, 2)) / ds, 10), mp.nstr((g(3, 3) - g(2, 3)) / ds, 10),
             mp.nstr(min(g(i, j) for i in range(1, NS + 1) for j in range(1, NS + 1)), 8)))

print("\n[F3] 數據 ±1σ 角點（entry-positive ㈠）")
for dru, drd in ((+1, +1), (+1, -1), (-1, +1), (-1, -1)):
    ru = RU + dru * SDU
    rd = RD + drd * SDD
    r = maximin(ru, rd, False)
    print("   ρ_u = %+0.4f、ρ_d = %+0.4f ⟹ %s（t* = %s）"
          % (-ru, rd, "可行" if (r.success and r.x[NS] > 0) else "不可行",
             "%.4g" % r.x[NS] if r.success else "—"))

print("\n" + "=" * 96)
print("[註] 見證 w 係數據反解，不作預言；本輪只判定**類別**可行性。")
print("=" * 96)