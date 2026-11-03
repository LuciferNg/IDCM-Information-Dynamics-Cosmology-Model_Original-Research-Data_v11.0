#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/03_null_capacity.py — 候選類 A（K ≥ 5 之離散 ±1 核）之**可裁決性上界**

目的：在耗費 34.75M 次全枚舉之前，先答一個更根本的問題 ——
      即使某一 ±1 核「命中」，該命中是否携帶資訊？

方法：以 M13 自身 pipeline（δ_x = Σ_n w_n φ_n(k_x)Δ_n）對隨機 ±1 核取樣，
      得 (z_d, z_u) 之經驗空分佈；再乘上類大小 N_class，得**預期命中數**。

判準（跑前凍結）：
  N_trials_frozen = 200000（MC 取樣；固定種子 20261029）
  K ∈ {2,3,4,5,6,8}
  類大小 N_class(K) = C(44,K)·2^K（精確整數，不取樣）
  命中窗 ε ∈ {0.10σ, 0.25σ, 0.50σ, 1.00σ}（雙邊，須同時 |z_d| 及 |z_u|）

結論形式：E[hits] = N_class × P̂(both within ε)。E ≫ 1 ⟹ 命中零資訊（不可正面裁決）；
          E ≪ 1 ⟹ 命中具資訊量（值得枚舉）。
本腳本不作候選判定，不宣稱任何命中。
"""
import mpmath as mp
import numpy as np
from math import comb

mp.mp.dps = 30
SEED = 20261029
N_TRIALS_FROZEN = 200_000
N = 44

phi = lambda n, j: float(mp.sqrt(mp.mpf(2)/45) * mp.sin(mp.pi*n*j/45))
D = np.array([phi(n, 3) - phi(n, 2) for n in range(1, N+1)])
V = np.array([[phi(n, k)*D[n-1] for n in range(1, N+1)] for k in (2, 3, 4)])  # u,d,s

RHO_D_OBS = 0.40350365
RHO_U_OBS = -0.02588837
SIG_RD = 0.044785178
SIG_RU = 0.033583899

rng = np.random.default_rng(SEED)

def zs(w):
    d = V @ w
    rd = d[1]/d[2]; ru = d[0]/d[2]
    return (rd-RHO_D_OBS)/SIG_RD, (ru-RHO_U_OBS)/SIG_RU

print("=== M13/03 null capacity: ±1 K-sparse kernels ===")
print(f"seed={SEED}  N_trials(MC)={N_TRIALS_FROZEN}  pipeline=M13/01 (validated to 1e-40)")

print("\n K   N_class=C(44,K)·2^K     MC |z_d| med   |z_u| med   frac|z_d|<1  corr(z_d,z_u)   frac(δ_s=0)")
table = {}
for K in (2, 3, 4, 5, 6, 8):
    W = np.zeros((N_TRIALS_FROZEN, N))
    for t in range(N_TRIALS_FROZEN):
        idx = rng.choice(N, K, replace=False)
        W[t, idx] = rng.choice([-1.0, 1.0], K)
    Dm = W @ V.T                       # (T,3)  → (δ_u, δ_d, δ_s)
    good = np.abs(Dm[:, 2]) > 1e-12    # δ_s = 0 之核：ρ 未定義（δ-盲支撐）
    frac_bad = 1.0 - good.mean()
    rd = np.where(good, Dm[:, 1]/np.where(good, Dm[:, 2], 1), np.nan)
    ru = np.where(good, Dm[:, 0]/np.where(good, Dm[:, 2], 1), np.nan)
    zd = (rd-RHO_D_OBS)/SIG_RD
    zu = (ru-RHO_U_OBS)/SIG_RU
    m = good & np.isfinite(zd) & np.isfinite(zu)
    zd, zu = zd[m], zu[m]
    ncls = comb(N, K)*2**K
    table[K] = (zd, zu, ncls, frac_bad)
    print(f" {K:2d}  {ncls:>20,d}   {np.median(abs(zd)):8.3f}   {np.median(abs(zu)):8.3f}   "
          f"{(abs(zd)<1).mean():7.4f}    {np.corrcoef(zd,zu)[0,1]:+.3f}   {frac_bad:.2e}")

print("\n=== 預期命中數 E[hits] = N_class × P̂(both < ε)（雙邊同取）===")
print(" K      ε=0.10σ        ε=0.25σ        ε=0.50σ        ε=1.00σ")
for K in (2, 3, 4, 5, 6, 8):
    zd, zu, ncls, _ = table[K]
    row = f" {K:2d} "
    for eps in (0.10, 0.25, 0.50, 1.00):
        p = float(np.mean((abs(zd) < eps) & (abs(zu) < eps)))
        row += f"  {ncls*p:12.4g}" + ("*" if p == 0 else " ")
    print(row)
print("  * = MC 200k 內 0 命中 ⟹ 僅得上界（95% CL：P < 1.5e-5 ⟹ E < N_class×1.5e-5）")

print("\n=== 資訊窗：E[hits] < 1 之最大 ε（由 P̂ 反演；受 MC 解析度限制）===")
print(" K      N_class      E<1 之最寬 ε    受 MC 限制？")
for K in (2, 3, 4, 5, 6, 8):
    zd, zu, ncls, _ = table[K]
    best = None
    for eps in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00):
        p = float(np.mean((abs(zd) < eps) & (abs(zu) < eps)))
        if p > 0 and ncls*p < 1:
            best = eps
    lim = ncls * 1.5e-5
    print(f" {K:2d}  {ncls:>12,d}   {('ε < ' + str(best)) if best else '無法判定（E 於 ε=0.05 已 ≳ 1）':>16s}"
          f"   {'上界 E<%.3g' % lim if lim > 1 else 'MC 足彀'}")

print("\n=== 判詞 ===")
zd5, zu5, nc5, fb5 = table[5]
p5_1 = float(np.mean((abs(zd5) < 1.0) & (abs(zu5) < 1.0)))
print(f"  K=5：N_class = {nc5:,d}；P̂(|z_d|<1σ ∧ |z_u|<1σ) = {p5_1:.4f}")
print(f"  ⟹ 預期「雙 1σ 命中」= {nc5*p5_1:,.0f} 個（≫ 1）")
print("  ⟹ 候選類 A 之任何 1σ 級命中皆**零資訊**；此類**不可正面裁決**，只可作排除（null）用途。")
print("  ⟹ 若日後枚舉 K=5 全族，唯一有意義之判據 = 「全族零命中」或「最小 z 顯著超越同族空分佈」，")
print("     報告前必須附 N_trials 與 MC null 分佈（M12 R28 P4 之方法）。")
print("\n  N_trials（本腳本）= MC 取樣 200,000；N_class 為精確組合數，不取樣。")
print("  註：δ_s = 0 之核（支撐落於 δ-盲模式 n=9,27）ρ 未定義，已剔除；其比例見 frac(δ_s=0) 欄。")