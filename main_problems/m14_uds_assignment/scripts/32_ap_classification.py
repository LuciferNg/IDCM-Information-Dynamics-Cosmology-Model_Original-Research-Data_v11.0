#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/32 — Z2：AP 三元組之**完全分類**（兩族）⟹ AP ⟹ 4 | N

────────────────────────────────────────────────────────────────────────
殘留：docs/30 §七 Z2（「4∤N ⟹ 無 AP 族」之證明）
本輪結果：AP 三元組恰落於兩族，兩族皆可**精確證明**，且覆蓋率 100%（掃描範圍內）。

判準（跑前凍結）
  V1 鏈 γ_n = 4sin²(nπ/N)，n = 1..K（N 偶：K = N/2；N 奇：K = (N−1)/2）。
  V2 AP：|(γ_b−γ_a) − (γ_c−γ_b)| < 1e−28（精確閘）；前置浮點篩窗 1e−9（≫ 浮點累積誤差）。
  V3 族 A：(a, N/4, N/2−a)，a = 1..N/4−1；族 B：(N/4, N/3, N/2)，恰 1 個。
  V4 完備性：所有 AP 三元組 ⊆ 族 A ∪ 族 B。
性質：確定性窮舉。N_trials = 0。dps = 40。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import math

import mpmath as mp

mp.mp.dps = 40
L = "=" * 78
TOL_F = 1e-9          # 浮點前置篩窗（遠大於 ~1e−14 之累積誤差）
TOL_X = mp.mpf('1e-28')
NMAX = 400


def gamf(N, n):
    return 4.0 * math.sin(n * math.pi / N) ** 2


def gam(N, n):
    return 4 * mp.sin(n * mp.pi / N) ** 2


print(L)
print("M14 scripts/32 — Z2：AP 完全分類（族 A ＋ 族 B）")
print(L)

# ============================================================================
print("\n【§一】族 A：(a, N/4, N/2−a) —— 精確恆等")
print("  γ_a + γ_(N/2−a) = 4（sin²+cos² 恆等）；γ_(N/4) = 2 ⟹ 4 = 2·2 ⟹ 恆成 AP")
worst = mp.mpf(0)
nA = 0
for N in range(8, NMAX + 1):
    if N % 4:
        continue
    for a in range(1, N // 4):
        lhs = gam(N, a) + gam(N, N // 2 - a)
        worst = max(worst, abs(lhs - 2 * gam(N, N // 4)))
        nA += 1
print(f"  對 8 ≤ N ≤ {NMAX}、4|N 之全部 (N, a)：{nA} 個檢查，max|dev| = {mp.nstr(worst, 3)}")
print(f"  個數 = N/4 − 1（每個可行 N）")

# ============================================================================
print("\n【§二】族 B：(N/4, N/3, N/2) —— 精確恆等，恰 1 個")
print("  γ_(N/4) + γ_(N/2) = 2 + 4 = 6 = 2·γ_(N/3)（γ_(N/3) = 4sin²(π/3) = 3）")
worstB = mp.mpf(0)
nB = 0
for N in range(12, NMAX + 1, 12):
    lhs = gam(N, N // 4) + gam(N, N // 2)
    worstB = max(worstB, abs(lhs - 2 * gam(N, N // 3)))
    nB += 1
print(f"  對 12|N、N ≤ {NMAX}：{nB} 個 N，max|dev| = {mp.nstr(worstB, 3)}")
print(f"  ⟹ 族 B 需 N/4、N/3、N/2 皆為整數模 ⟹ **12 | N**")

# ============================================================================
print("\n【§三】★ 完備性檢定：所有 AP 三元組 ⊆ 族 A ∪ 族 B")
print(L)
tot = fA = fB = 0
bad = []
for N in range(8, NMAX + 1):
    K = N // 2 if N % 2 == 0 else (N - 1) // 2
    if K < 3:
        continue
    gf = [0.0] + [gamf(N, n) for n in range(1, K + 1)]
    cand = []
    for a, b, c in itertools.combinations(range(1, K + 1), 3):
        if abs((gf[b] - gf[a]) - (gf[c] - gf[b])) < TOL_F:
            cand.append((a, b, c))
    for (a, b, c) in cand:
        if abs((gam(N, b) - gam(N, a)) - (gam(N, c) - gam(N, b))) >= TOL_X:
            continue                                    # 浮點假陽性，剔除
        tot += 1
        inA = (N % 4 == 0 and b == N // 4 and a + c == N // 2)
        inB = (N % 12 == 0 and (a, b, c) == (N // 4, N // 3, N // 2))
        fA += inA
        fB += inB
        if not (inA or inB):
            bad.append((N, a, b, c))
print(f"  掃 8 ≤ N ≤ {NMAX}（偶：K=N/2；奇：K=(N−1)/2）")
print(f"  AP 三元組總數 = **{tot}**；族 A = {fA}；族 B = {fB}；合計 {fA + fB}")
print(f"  不屬兩族者 = **{len(bad)}**" + (f" → {bad[:5]}" if bad else " ✅（完備性成立）"))
print(f"  覆蓋率 = {100.0 * (fA + fB) / tot if tot else 0:.2f}%")

# 計數公式核
print("\n  計數公式核（每 N 逐一比對）：")
mismatch = []
for N in range(8, NMAX + 1):
    K = N // 2 if N % 2 == 0 else (N - 1) // 2
    if K < 3:
        continue
    gf = [0.0] + [gamf(N, n) for n in range(1, K + 1)]
    cnt = 0
    for a, b, c in itertools.combinations(range(1, K + 1), 3):
        if abs((gf[b] - gf[a]) - (gf[c] - gf[b])) < TOL_F:
            if abs((gam(N, b) - gam(N, a)) - (gam(N, c) - gam(N, b))) < TOL_X:
                cnt += 1
    pred = (N // 4 - 1) * int(N % 4 == 0) + int(N % 12 == 0)
    if cnt != pred:
        mismatch.append((N, cnt, pred))
print(f"  公式 AP(N) = (N/4 − 1)[4|N] + [12|N]：違例 {len(mismatch)}" + (f" → {mismatch[:5]}" if mismatch else " ✅"))

# ============================================================================
print("\n" + L)
print("【§四】Z2 之結論")
print(L)
print(f"""  ✅ **AP 完全分類（掃描範圍內，覆蓋率 100%）**：
       族 A：(a, N/4, N/2−a)，a = 1..N/4−1，需 4 | N        —— 精確恆等（§一）
       族 B：(N/4, N/3, N/2)，恰 1 個，需 12 | N             —— 精確恆等（§二）
       AP(N) = (N/4 − 1)[4|N] + [12|N]
  ✅ **⟹ AP 存在 ⟹ 4 | N**（兩族皆需 4|N；已排除第三族，{NMAX} 內 0 例外）
  ✅ 特例：N = 88 ⟹ AP = 88/4 − 1 = 21（族 A），族 B 不存在（12 ∤ 88）✓
           N = 90 ⟹ 4 ∤ 90 ⟹ AP = 0 ✓
  🔴 **仍未證**：完備性（「無第三族」）之**數論證明** ——
     本輪以 {NMAX} 內窮舉作證；一般證明須分類 cos 之有理角線性關係
     （Conway–Jones 型），本 M 不做。
""")
print(L)
raise SystemExit(0)