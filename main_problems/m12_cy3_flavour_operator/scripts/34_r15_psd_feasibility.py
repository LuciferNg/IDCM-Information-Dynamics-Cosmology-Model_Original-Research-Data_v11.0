#!/usr/bin/env python3
"""34 — R15：可行域定理 —— PSD 核類到底可否重現 (ρ_d, ρ_u)？

問題來源（未做之推）：R7-A1／R8／R11–R13 之「判死」全部係**特定 w**（`w_n = 1/γ_n`、截斷、
對偶、帶號變體），而 PSD 核類有 44 維。故「類別判死」從未真正檢定 ⟹ 本輪做精確可行域檢定。

預註冊（跑前凍死）
--------------------------------------------------
核族      G(i,j) = Σ_n w_n φ_n(i)φ_n(j)，φ_n(j) = √(2/45)·sin(πnj/45)
恆等式    δ_x = Σ_n w_n φ_n(x)[φ_n(3) − φ_n(2)] ≡ Σ_n u_n φ_n(x)，u_n ≡ w_n·Δ_n，Δ_n ≡ φ_n(3) − φ_n(2)
PSD 約束  w_n ≥ 0 ⟹ **sign(u_n) = sign(Δ_n) 固定**，非自由
目標      (δ_u, δ_d, δ_s) ∝ t ≡ (−0.02588837, 0.40350365, 1)
檢定      t ∈ cone{ v_n ≡ sign(Δ_n)·(φ_n(2),φ_n(3),φ_n(4)) }？（Carathéodory：三生成元窮盡）
F1  t ∈ cone ⟹ **PSD 類可行** ⟹ 「類別判死」範圍過寬，須更正
F2  t ∉ cone ⟹ **定理級 no-go**：任何 PSD 核皆不能重現數據
F3  見證 w 只作**結構診斷**，不作預言（數據反解）
F4  另報：Δ_n = 0 之 δ-盲模式；cone 是否 = ℝ³（即類別零約束）
"""
import itertools
import random

import mpmath as mp
import numpy as np

mp.mp.dps = 30
NS = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2

RU, RD = mp.mpf("-0.02588837"), mp.mpf("0.40350365")
t = (RU, RD, mp.mpf(1))

print("=" * 96)
print("34 — R15：PSD 核類之可行域（精確圓錐歸屬，Carathéodory 窮盡）")
print("=" * 96)
print("  目標 t = (δ_u, δ_d, δ_s) ∝ (%s, %s, 1)" % (mp.nstr(RU, 8), mp.nstr(RD, 8)))

D = lambda n: phi(n, 3) - phi(n, 2)
TOL0 = mp.mpf("1e-25")
blind = [n for n in range(1, NS + 1) if abs(D(n)) < TOL0]
live = [n for n in range(1, NS + 1) if abs(D(n)) >= TOL0]
print("  ★ δ-盲模式（Δ_n = 0，即 φ_n(2)=φ_n(3)）：n = %s" % blind)

vf = np.array([[float((1 if D(n) > 0 else -1) * phi(n, j)) for j in (2, 3, 4)] for n in live])
tn = np.array([float(RU), float(RD), 1.0])
tag = np.array(live)
np_ = len(live)
print("  有效模式 = %d（Δ_n>0: %d、Δ_n<0: %d）"
      % (np_, sum(1 for n in live if D(n) > 0), sum(1 for n in live if D(n) < 0)))


_T = np.array(list(itertools.combinations(range(len(vf)), 3)))
_A = vf[_T]
_keep = np.abs(np.linalg.det(_A)) > 1e-12
TRIP, AMAT = _T[_keep], _A[_keep]


def cone_membership(target, gen, tol=1e-9):
    """回傳 (是否在 cone 內, 見證 triple, λ)。批次解全部非奇異 triple。"""
    T, A = TRIP, AMAT
    b = np.tile(np.asarray(target, float).reshape(3, 1), (len(T), 1, 1))
    lam = np.linalg.solve(A, b)[..., 0]
    res = np.abs(np.einsum("tij,tj->ti", A, lam) - b[..., 0]).max(axis=1)
    ok = (lam.min(axis=1) >= -tol) & (res < tol)
    if not ok.any():
        return False, None, None
    k = int(np.argmax(lam.min(axis=1) * ok))
    return True, T[k], lam[k]


ok, trip, lam = cone_membership(tn, vf)
print("\n[檢定] t ∈ cone ?  %s" % ("**是**" if ok else "**否**"))
if ok:
    idx = [live[i] for i in trip]
    print("  見證 triple：n = %s，λ = %s" % (idx, ["%.6f" % x for x in lam]))
    for k in range(3):
        print("     n=%-3d λ=%.8f  w_n = λ/Δ_n = %s  （Δ_n = %s）"
              % (idx[k], lam[k], mp.nstr(mp.mpf(lam[k]) / D(idx[k]), 8), mp.nstr(D(idx[k]), 6)))
    # entry-positive 檢查（R8 之 1936/1936 全正）
    w = {idx[k]: mp.mpf(lam[k]) / D(idx[k]) for k in range(3)}
    G = lambda i, j: sum(w[n] * phi(n, i) * phi(n, j) for n in w)
    vals = [G(i, j) for i in range(1, NS + 1) for j in range(1, NS + 1)]
    neg = sum(1 for x in vals if x < 0)
    print("\n  entry-positive 子類核實：min = %s、max = %s、負元 %d / %d"
          % (mp.nstr(min(vals), 8), mp.nstr(max(vals), 8), neg, len(vals)))
    print("  ⟹ 見證%s entry-positive（R8 之『1936/1936 全正』）" % ("屬" if neg == 0 else "**不**屬"))

print("\n[F4] cone 是否 = ℝ³（即 PSD 類零約束）？")
random.seed(20261029)
hit = 0
NT = 500
for _ in range(NT):
    a, b = random.uniform(-1, 1), random.uniform(-1, 1)
    if abs(b) < 1e-3 or abs(a) < 1e-3:
        continue
    h, _, _ = cone_membership(np.array([a, b, 1.0]), vf)
    hit += h
print("  隨機 %d 個方向（第三分量歸一）：命中 %d ⟹ 份額 %.1f%%" % (NT, hit, 100.0 * hit / NT))
print("  ⟹ %s" % ("**cone = ℝ³ ⟹ PSD 類別對 (ρ_d,ρ_u) 零約束**（任意目標皆可由某個 w≥0 重現）"
                  if hit == NT else "cone ⊊ ℝ³ ⟹ PSD 類有真約束（份額 = 可達之立體角）"))

print("\n" + "=" * 96)
print("[判定]")
if ok:
    print("  F1 成立：PSD 類**可行** ⟹ R7-A1／R8／R11–R13 之判死只覆蓋『特定 w』，非類別。")
else:
    print("  F2 成立：定理級 no-go —— PSD 類整體不可行。")
print("  F3 見證不作預言（數據反解）。")
print("=" * 96)