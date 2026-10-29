#!/usr/bin/env python3
"""35 — R15b：entry-positive 子類是否仍可行？（R8「1936/1936 全正」之範圍檢定）

34 之果：PSD 類（w_n ≥ 0）於 (ρ_d, ρ_u) 上 **零約束**（cone = ℝ³、500/500 隨機方向命中）。
故 R7-A1／R8 之「範疇判死」只可能靠**額外**結構（entry-positive、單調、1/γ 型）成立。
本輪檢最強額外條件：**G(i,j) > 0 ∀ i,j（1936 全正 = R8 之實際判死依據）**。

預註冊
--------------------------------------------------
族      G = Σ_{n∈S} w_n φ_n φ_n，|S| ≤ 3（Carathéodory 完備），w ≥ 0，G(i,j) > 0 ∀ i,j
目標    t = (δ_u, δ_d, δ_s) ∝ (−0.02588837, 0.40350365, 1)
檢定    窮舉全部 triple ⟹ 篩 λ≥0 且殘差<1e-9 ⟹ 逐個建 G（44×44, numpy）⟹ 數 entry-positive 者
F1  存在 entry-positive 見證 ⟹ **R8 之範疇判死亦屬過寬**（entry-positive 對數據無約束）
F2  不存在 ⟹ **entry-positive 係真類別約束** ⟹ 首個成立之類別級 theorem
F3  見證只作結構診斷，不作預言
"""
import itertools

import mpmath as mp
import numpy as np

mp.mp.dps = 30
NS = 44
phi_m = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
D_m = lambda n: phi_m(n, 3) - phi_m(n, 2)
live = [n for n in range(1, NS + 1) if abs(D_m(n)) >= mp.mpf("1e-25")]

RU, RD = mp.mpf("-0.02588837"), mp.mpf("0.40350365")
tn = np.array([float(RU), float(RD), 1.0])

PHI = np.array([[float(phi_m(n, j)) for j in range(1, NS + 1)] for n in live])   # (Nn, 44)
vf = np.array([[float((1 if D_m(n) > 0 else -1) * phi_m(n, j)) for j in (2, 3, 4)] for n in live])
Dn = np.array([float(D_m(n)) for n in live])

print("=" * 96)
print("35 — R15b：entry-positive 子類可行性（R8 範圍檢定）")
print("=" * 96)
print("  有效模式 %d（δ-盲 n=9,27 已排除）；triple 總數 = %d"
      % (len(live), len(list(itertools.combinations(range(len(live)), 3)))))

TRIP = np.array(list(itertools.combinations(range(len(live)), 3)))
A = vf[TRIP]
det = np.linalg.det(A)
keep = np.abs(det) > 1e-12
TRIP, A = TRIP[keep], A[keep]
lam = np.linalg.solve(A, np.tile(tn.reshape(3, 1), (len(A), 1, 1)))[..., 0]
res = np.abs(np.einsum("tij,tj->ti", A, lam) - tn).max(axis=1)
feas = (lam.min(axis=1) >= -1e-9) & (res < 1e-9)
print("  可行 triple（λ≥0 且殘差<1e-9）= %d ／ %d" % (feas.sum(), len(TRIP)))

idx_feas = np.where(feas)[0]
pos_cnt = 0
best = None
for k in idx_feas:
    trip = TRIP[k]
    w = lam[k] / Dn[trip]                                   # w_n = λ/Δ_n ≥ 0 ?
    if w.min() < -1e-12:
        continue
    G = PHI[trip].T @ np.diag(w) @ PHI[trip]                # 44×44
    if G.min() > 0:
        pos_cnt += 1
        if best is None:
            best = (trip, w, G.min(), G.max())

print("\n[檢定] entry-positive（G 全正）之可行見證數 = %d" % pos_cnt)
if pos_cnt:
    trip, w, gmin, gmax = best
    print("  ⟹ **F1 成立**：entry-positive 子類**可行**，R8 之範疇判死亦屬過寬")
    print("     見證：n = %s" % [live[i] for i in trip])
    print("     w_n = %s" % ["%.6g" % x for x in w])
    print("     G 元素：min = %.6g、max = %.6g（全正）" % (gmin, gmax))
    print("     驗算 ρ_d、ρ_u（mpmath, dps=30）：")
    wm = {live[trip[c]]: mp.mpf(str(w[c])) for c in range(3)}
    gm = lambda i, j: sum(wm[n] * phi_m(n, i) * phi_m(n, j) for n in wm)
    ds = gm(3, 4) - gm(2, 4)
    print("       δ_u/δ_s = %s   δ_d/δ_s = %s"
          % (mp.nstr((gm(3, 2) - gm(2, 2)) / ds, 10), mp.nstr((gm(3, 3) - gm(2, 3)) / ds, 10)))
else:
    print("  ⟹ **F2 成立**：entry-positive 係真類別約束（首個成立之類別級 theorem）")

print("\n" + "=" * 96)
print("[註] 見證只作結構診斷。entry-positive 可行 ≠ 有結構規則揀到該 w。")
print("=" * 96)