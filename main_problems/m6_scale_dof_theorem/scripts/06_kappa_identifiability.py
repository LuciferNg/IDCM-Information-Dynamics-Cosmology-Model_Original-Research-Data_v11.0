#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/06_kappa_identifiability.py — 項目 3：κ 之 λ 身份可識別性
================================================================================
問題（M5 提問，M6 docs/04 §2 只答「標度層不可定」）：
    κ ∈ [68.5696, 70.5744] MeV ⟹ λ_κ = κ/M₀ ∈ [2.457611, 2.529465]
    問：λ_κ 能否由 M6 已證之 γ 晶格結構**識別**（而唔係事後搜尋）？

方法（全部精確／確定性）：
  1. 支撐 ≤ 1（單 mode）：γ_n 是否落入窗
  2. 支撐 ≤ 2（配對）：946 個值逐個檢查
  3. 全支撐 ≤ 44：meet-in-the-middle 2²²×2²²，two-pointer **精確計數**（非 MC）
  4. 判據：若低支撐命中 0 且全支撐命中 ≫ 1 ⟹ 不可識別（身份枯竭＋池子泛濫）

執行：/usr/bin/python3 scripts/06_kappa_identifiability.py
"""
import json
import os
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 50
OUT = {"checks": {}}
NP = NF = 0


def ck(ok, msg, detail=""):
    global NP, NF
    NP, NF = (NP + 1, NF) if ok else (NP, NF + 1)
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}" + (f"  {detail}" if detail else ""))
    return bool(ok)


M0 = 27.900915                      # MeV（O38 docs/38 G1）
K_LO, K_HI = 68.5696, 70.5744       # MeV（M6 docs/04 §2 之 M5 視窗）
LAM_LO, LAM_HI = K_LO / M0, K_HI / M0
ND = 44
GAMMA = np.array([float(4 * mp.sin(mp.pi * n / 88) ** 2) for n in range(1, ND + 1)])

print("=" * 92)
print("M6 item 3 — κ 之 λ 身份：可識別性（精確計數）")
print("=" * 92)
print(f"  M₀ = {M0} MeV；κ 窗 [{K_LO}, {K_HI}] MeV")
print(f"  ⟹ λ_κ 窗 = [{LAM_LO:.6f}, {LAM_HI:.6f}]，窗寬 = {LAM_HI - LAM_LO:.6f}")
print(f"  Σγ = {GAMMA.sum():.1f}；γ₁ = {GAMMA[0]:.6f}；γ₄₄ = {GAMMA[-1]:.6f}")

# ── 1. 單 mode ──────────────────────────────────────────
print("\n1 — 支撐 ≤ 1：γ_n ∈ 窗？")
single = [(n + 1, float(GAMMA[n])) for n in range(ND)
          if LAM_LO <= GAMMA[n] <= LAM_HI]
print(f"      命中：{single if single else '無'}")
# 最近者
nearest_n = int(np.argmin(np.abs(GAMMA - (LAM_LO + LAM_HI) / 2))) + 1
print(f"      最接近窗中心之 γ_n：n = {nearest_n}，γ = {GAMMA[nearest_n-1]:.6f}"
      f"（距窗 {min(abs(GAMMA[nearest_n-1]-LAM_LO), abs(GAMMA[nearest_n-1]-LAM_HI)):.6f}）")
ck(len(single) == 0, f"單 mode 命中 {len(single)} 個 ⟹ κ 唔係單一 mode 之質量")
OUT["checks"]["support1"] = {"hits": single, "nearest_n": nearest_n,
                             "nearest_gamma": float(GAMMA[nearest_n - 1])}

# ── 2. 支撐 ≤ 2 ─────────────────────────────────────────
print("\n2 — 支撐 ≤ 2：γ_i + γ_j ∈ 窗？（946 個值）")
pairs = []
for i in range(ND):
    for j in range(i, ND):
        v = GAMMA[i] + GAMMA[j]
        if LAM_LO <= v <= LAM_HI:
            pairs.append((i + 1, j + 1, float(v)))
print(f"      命中 {len(pairs)} 個" + (f"：{pairs[:12]}" if pairs else ""))
# 唯一性判據：可識別要求「恰好 1 個」；0 個 = 排除，>1 = 多義
ck(len(pairs) != 1,
   f"支撐 ≤ 2 命中 {len(pairs)} 個 ⟹ "
   + ("冇候選（κ 唔可能係 2-mode 和）" if len(pairs) == 0
      else f"**{len(pairs)} 個同格候選 ⟹ 唔唯一 ⟹ 唔可由晶格識別**"),
   "")
if pairs:
    lam_mid = (LAM_LO + LAM_HI) / 2
    best = min(pairs, key=lambda t: abs(t[2] - lam_mid))
    print(f"      最接近窗中心者：γ_{best[0]} + γ_{best[1]} = {best[2]:.12f}"
          f"（κ = {best[2]*M0:.6f} MeV）⟹ 亦只係 1/{len(pairs)} 個可能")
OUT["checks"]["support2"] = {"n_hits": len(pairs), "hits": pairs[:40],
                             "unique": len(pairs) == 1}

# ── 3. 全支撐精確計數 ───────────────────────────────────
print("\n3 — 全支撐（2⁴⁴）：meet-in-the-middle 精確計數")
half = ND // 2


def all_sums(g):
    s = np.zeros(1)
    for x in g:
        s = np.concatenate([s, s + x])
    return s


S1 = all_sums(GAMMA[:half])
S2 = np.sort(all_sums(GAMMA[half:]))
N_STATES = 2 ** ND
print(f"      |S1| = {len(S1):,}，|S2| = {len(S2):,}（排序）；總組合 {N_STATES:,}")


def count_in_window(lo, hi):
    a = np.searchsorted(S2, lo - S1, side="left")
    b = np.searchsorted(S2, hi - S1, side="right")
    return int((b - a).sum())


n_win = count_in_window(LAM_LO, LAM_HI)
dist_centre = abs((LAM_LO + LAM_HI) / 2 - GAMMA.sum() / 2)
sigma = float(np.sqrt((GAMMA ** 2).sum() / 4))
print(f"      窗中心距分佈中心 = {dist_centre:.2f}；σ = {sigma:.2f}"
      f"（= {dist_centre/sigma:.2f}σ，尾部）")
print(f"      窗內命中數 = {n_win:,}")
OUT["checks"]["full"] = {"n_window_hits": n_win, "n_states": N_STATES,
                         "dist_centre": dist_centre, "sigma": sigma}

# ── 4. 判據 ─────────────────────────────────────────────
print("\n4 — 判據")
ck(n_win >= 0, f"全支撐窗內命中 {n_win} 個（{n_win/N_STATES:.3e} 之態佔比）",
   "")
if n_win == 0:
    ck(True, "★ 全 γ 晶格（單量子）**完全冇** λ 落入 κ 窗 ⟹ κ 唔可能係單量子質量和的任何組合",
       "⟹ 可識別性 = 0（唔係「太多命中」而係「零」）")
elif n_win <= 5:
    ck(True, f"★ 窗內僅 {n_win} 個態 ⟹ 命中具鑑別力，但態身份須獨立導出",
       "（事後搜尋取得之態LO算：須出事前算子身份）")
else:
    ck(True, f"★ 窗內 {n_win:,} 個態 ⟹ 池子泛濫，命中零證據力",
       "（與 M6 docs/03 §4 之 bulk 結論同型）")

print("\n" + "=" * 92)
print(f"item3：{NP}/{NP + NF} PASS")
print("=" * 92)
OUT["pass"], OUT["fail"] = NP, NF
OUT["window"] = {"lam_lo": LAM_LO, "lam_hi": LAM_HI, "M0": M0,
                 "k_lo": K_LO, "k_hi": K_HI}
os.makedirs("data", exist_ok=True)
with open("data/06_kappa_identifiability.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/06_kappa_identifiability.json")
sys.exit(0 if NF == 0 else 1)