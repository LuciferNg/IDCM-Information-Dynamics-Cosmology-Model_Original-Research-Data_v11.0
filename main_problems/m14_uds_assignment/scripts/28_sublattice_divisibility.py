#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/28 — 子晶格存在性（V1）：mod-m 子晶格 ⟺ m | N

────────────────────────────────────────────────────────────────────────
殘留：docs/26 §五 V1（/88 之 mod-3 子晶格對應物）
背景：docs/18 §四 在 /90 上發現 γ_(3k) = 4sin²(kπ/30) ⟹ 14 模 /30 鏈。
問題：/88 有無對應物？一般判準為何？

判準（跑前凍結）
  R1 子晶格（步 d）：取 n = dk，k = 1..⌊44/d⌋；要求 γ_n = 4sin²(kπ/M)，M = N/d ∈ ℤ。
  R2 「成鏈」= K = ⌊44/d⌋ 等於 M/2（即取盡半譜）；否則為「部分子鏈」。
  R3 Σγ 以 dps=60 直接求和，並與 docs/21 之 2M∓2 律比對。
性質：確定性窮舉（d | N）。N_trials = 0。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import math

import mpmath as mp

mp.mp.dps = 60
L = "=" * 78
NMODE = 44


def gam(x, n):
    return 4 * mp.sin(n * mp.pi / x) ** 2


def divisors(n):
    return [d for d in range(2, n + 1) if n % d == 0]


print(L)
print("M14 scripts/28 — 子晶格存在性（V1）：mod-m 子晶格 ⟺ m | N")
print(L)

report = {}
for x, tag in ((88, "/88"), (90, "/90")):
    print(f"\n【§一{tag}】步 d | {x} 之子晶格")
    rows = []
    for d in divisors(x):
        M = x // d
        K = NMODE // d
        if K < 1:
            rows.append((d, M, K, None, "空"))
            continue
        # 核 γ_(dk) = 4sin²(kπ/M)
        dev = max(abs(gam(x, d * k) - 4 * mp.sin(k * mp.pi / M) ** 2) for k in range(1, K + 1))
        sig = mp.fsum(gam(x, d * k) for k in range(1, K + 1))
        full = (K == M // 2)   # 半譜模數 = ⌊M/2⌋（M 奇時 M/2 非模）
        rows.append((d, M, K, sig, f"max|dev|={mp.nstr(dev, 3)}", full))
    for r in rows:
        if r[3] is None:
            print(f"    d={r[0]:>2}  M={r[1]:>3}  K={r[2]}  空")
        else:
            _, M, K, sig, dvtxt, full = r
            print(f"    d={r[0]:>2}  M={M:>3}  K={K:>2}  Σγ={mp.nstr(sig, 10):>14s}  "
                  f"{'成鏈' if full else '部分'}  {dvtxt}")
    report[tag] = rows

print("\n" + L)
print("【§二】★ 判準：mod-m 子晶格存在 ⟺ m | N")
print(L)
print(f"  /88 之因子（>1）：{divisors(88)}")
print(f"  /90 之因子（>1）：{divisors(90)}")
print(f"  ⟹ 3 | 90 **但** 3 ∤ 88 ⟹ /88 **冇** mod-3 子晶格（γ_(3k) 之角 3kπ/88 不能寫成 kπ/M）")
print(f"  ⟹ /90 之「14 模 /30 鏈」係 3 | 90 之**專屬後果**，不可移植")
print(f"  ⟹ /88 之最近模數係 11（88 = 8×11）：d=11 ⟹ M=8, K=4")

print("\n" + L)
print("【§三】兩鏈子晶格數與成鏈數")
print(L)
for tag, rows in report.items():
    ok = [r for r in rows if r[3] is not None]
    full = [r for r in ok if r[5]]
    print(f"  {tag}：非空子晶格 {len(ok)} 個（d = {[r[0] for r in ok]}）；"
          f"成鏈 {len(full)} 個（d = {[r[0] for r in full]} ⟹ 鏈 {['/' + str(r[1]) for r in full]}）")

print("\n" + L)
print("【§四】/90 之 /30 鏈精確核（docs/18 §四 之獨立重算）")
print(L)
M, K = 30, 14
dev = max(abs(gam(90, 3 * k) - 4 * mp.sin(k * mp.pi / M) ** 2) for k in range(1, K + 1))
sig = mp.fsum(gam(90, 3 * k) for k in range(1, K + 1))
full_modes = M - 1                      # /30 Dirichlet 全譜：k = 1..29
sig_full = mp.fsum(4 * mp.sin(k * mp.pi / M) ** 2 for k in range(1, full_modes + 1))
sig_half = (sig_full - 4) / 2           # 半譜（k ≤ 14，排除 Nyquist k=15 之 γ=4）
print(f"  γ^(90)_(3k) = 4sin²(kπ/30)，k = 1..14：max|dev| = {mp.nstr(dev, 3)}")
print(f"  /30 全譜（k = 1..{full_modes}）Σγ = 2M = {mp.nstr(sig_full, 12)}"
      f"（Σ_{{k=1}}^{{M−1}} 4sin²(kπ/M) = 2M）；Nyquist k=15 之 γ = 4")
print(f"  /30 半譜（k = 1..14）Σγ = (2M − 4)/2 = {mp.nstr(sig_half, 12)}")
print(f"  本子晶格 Σγ = {mp.nstr(sig, 12)} ⟹ 與半譜完全一致 ⟹ **係完整半譜，非截斷**")
print(f"  下降至 /90 之條件：n = 3k ≤ 44 ⟹ k ≤ 14；Nyquist（k=15 ⟹ n=45）落於 /90 模域之外")
print(f"  ⟹ 與 docs/20 §三「和頻通道在 /90 關閉（落點 m=45 為邊界節點）」互為同一事實之兩面")
print(f"  註：§三「成鏈」判據為 K = ⌊M/2⌋（含 Nyquist）；d=3 之 K=14 ≠ 15 故標『部分』——"
      f"此與『除 Nyquist 外之半譜完整』並不矛盾（Nyquist 模 n=45 不在 /90 模域）")

print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ✅ V1 得答：mod-m 子晶格存在 ⟺ m | N（步 d ⟹ γ_n = 4sin²(kπ/(N/d))）
     /88：6 個非空（d = 2,4,8,11,22,44），**無 mod-3**；
     /90：9 個非空（d = 2,3,5,6,9,10,15,18,30），含 d=3 ⟹ /30 鏈。
  ✅ 14 模 /30 = /30 Dirichlet 之**完整半譜**（Σγ = 28 = (2M−4)/2 精確），非截斷；
     其 Nyquist 模落 n = 45（/90 模域外）⟹ 與 docs/20 §三 同源。來源為 **3 | 90 之算術專屬**，非兩鏈共有結構 ⟹ 不可移植至 /88。
  ✅ 與 docs/10 T1/T3（3 ∤ 88）、docs/24（i10 gcd 類 41）、docs/26（對偶）四者一致：
     「3」在 /88 側**完全缺席**，只以 3 | 90 之形式出現在 /90 側。
  🔴 未做：/88 之 d=11 子鏈（/8，4 模）之共振結構（是否為 /88 側「最細」之自洽扇區）。
""")
print(L)
raise SystemExit(0)