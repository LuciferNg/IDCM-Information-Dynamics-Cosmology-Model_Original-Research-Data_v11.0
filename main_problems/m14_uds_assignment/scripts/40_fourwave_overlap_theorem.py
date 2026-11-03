#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/40 — 四波（p = 4）重疊之**精確閉式**與非零性判準；兼 p 模宇稱定理

────────────────────────────────────────────────────────────────────────
待辦（`docs/18 §六 MR1` 之「四波之動能學」之前置）：
  三波重疊 χ³ 已有乘積閉式（docs/37）。四波 χ⁴ 一直只知「頂點常數全等」，不知選擇律。

本腳本推導（結構，非拟合）：
  ① 宇稱定理（任意 p）：χ⁽ᵖ⁾ = 0 除非 Σ n_i ≡ p (mod 2)
       反射 j ↦ 45−j：φ_n(45−j) = (−1)^{n+1} φ_n(j)
       ⟹ χ⁽ᵖ⁾ = (−1)^{Σn_i + p} χ⁽ᵖ⁾
  ② 鏈和（N = 45）：Σ_{j=1}^{44} cos(πkj/45) = 44（k ≡ 0 mod 90）｜−1（k 偶）｜0（k 奇）
  ③ 四正弦展開（係數已由數值擬合確認，見 §二）：
       sinA sinB sinC sinD = (1/8)[ +cos(A−B−C+D) +cos(A−B+C−D)
                                   −cos(A−B−C−D) −cos(A−B+C+D)
                                   −cos(A+B−C+D) −cos(A+B+C−D)
                                   +cos(A+B−C−D) +cos(A+B+C+D) ]
  ④ 代入 ②：Σ(even k) 得 −1，而 8 個係數和 = 0 ⟹ 只剩 k ≡ 0 (mod 90) 之項
       ⟹ **χ⁴ = W/90**，W = Σ_ε c_ε·[k_ε ≡ 0 (mod 90)] ∈ ℤ
  ⑤ 非零判準：χ⁴ ≠ 0 ⟺ W ≠ 0 ⟺ 某個帶號和 ∈ {0, ±90}

判準（跑前凍結）
  V1 鏈和恆等式（②）：k = 0..180 全檢
  V2 四正弦展開（③）：任意 200 組 (A,B,C,D) 數值對比
  V3 宇稱定理：p = 3,4,5 各抽 400 組
  V4 閉式 χ⁴ = W/90：n ≤ 12 全窮舉（20,736 組）＋ 隨機 20,000 組（含奇和）
  V5 W 之值域與分佈（有界性觀察）
  V6 χ⁴ ∈ (1/90)ℤ ⟹ 90·χ⁴ 必為整數（直接求和側獨立驗）
性質：確定性（隨機部分以固定 seed）。N_trials = 20000（隨機）＋全窮舉。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import random

import mpmath as mp

mp.mp.dps = 60
N = 45
C2 = mp.mpf(2) / 45
L = "=" * 78
R = [1, 1, -1, -1, -1, -1, 1, 1]
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]


def Scos(k):
    return mp.mpf(44) if k % 90 == 0 else (mp.mpf(-1) if k % 2 == 0 else mp.mpf(0))


def W4(a, b, c, d):
    t = (a, b, c, d)
    w = 0
    for cs, sg in zip(R, SIG):
        if sum(s * x for s, x in zip(sg, t)) % 90 == 0:
            w += cs
    return w


def chi4_direct(a, b, c, d):
    s = mp.mpf(0)
    for j in range(1, 45):
        s += (mp.sin(mp.pi * a * j / 45) * mp.sin(mp.pi * b * j / 45)
              * mp.sin(mp.pi * c * j / 45) * mp.sin(mp.pi * d * j / 45))
    return C2 ** 2 * s


def chi_p(ns):
    s = mp.mpf(0)
    for j in range(1, 45):
        p = mp.mpf(1)
        for n in ns:
            p *= mp.sin(mp.pi * n * j / 45)
        s += p
    return C2 ** (mp.mpf(len(ns)) / 2) * s


print(L)
print("M14 scripts/40 — 四波重疊閉式（χ⁴ = W/90）＋ p 模宇稱定理")
print(L)

# ---------------------------------------------------------------- §一
print("\n【§一】鏈和恆等式  Σ_{j=1}^{44} cos(πkj/45)")
bad = [k for k in range(0, 181)
       if abs(sum(mp.cos(mp.pi * k * j / 45) for j in range(1, 45)) - Scos(k)) > mp.mpf('1e-40')]
print(f"  k = 0..180 全檢：違例 = {len(bad)} {'✅' if not bad else bad[:5]}")
print("  值：44（k ≡ 0 mod 90）｜−1（k 偶，非 90 倍）｜0（k 奇）")

# ---------------------------------------------------------------- §二
print("\n" + L)
print("【§二】四正弦展開係數（① +cos +cos −cos −cos −cos −cos +cos +cos）")
random.seed(3)
bad2 = 0
for _ in range(200):
    A, B, Cc, D = [mp.mpf(random.randint(1, 44)) * mp.pi / 45 for _ in range(4)]
    lhs = mp.sin(A) * mp.sin(B) * mp.sin(Cc) * mp.sin(D)
    rhs = mp.mpf(1) / 8 * (mp.cos(A - B - Cc + D) + mp.cos(A - B + Cc - D)
                          - mp.cos(A - B - Cc - D) - mp.cos(A - B + Cc + D)
                          - mp.cos(A + B - Cc + D) - mp.cos(A + B + Cc - D)
                          + mp.cos(A + B - Cc - D) + mp.cos(A + B + Cc + D))
    if abs(lhs - rhs) > mp.mpf('1e-40'):
        bad2 += 1
print(f"  200 組隨機檢：違例 = {bad2} {'✅' if bad2 == 0 else '🔴'}")

# ---------------------------------------------------------------- §三
print("\n" + L)
print("【§三】宇稱定理（單向）：parity 不容（Σn ≢ p mod 2）⟹ χ⁽ᵖ⁾ = 0")
print("       注意：parity 容許**唔保證**非零 —— p = 4 另有 W = 0 之額外零（§四）")
for p in (3, 4, 5):
    badp = viol_nz = 0
    for _ in range(400):
        ns = tuple(random.randint(1, 44) for _ in range(p))
        forbidden = (sum(ns) - p) % 2 != 0
        is_zero = abs(chi_p(ns)) < mp.mpf('1e-40')
        if forbidden and not is_zero:
            badp += 1          # 真違例：parity 不容但非零
        if (not forbidden) and is_zero:
            viol_nz += 1       # 容許但為零（唔算違例，另行統計）
    print(f"  p = {p}：抽 400 組｜parity 不容而非零 = {badp} {'✅' if badp == 0 else '🔴'}"
          f"｜容許但為零 = {viol_nz}")

# ---------------------------------------------------------------- §四
print("\n" + L)
print("【§四】V4 閉式 χ⁴ = W/90")
print(L)
bad4 = 0
cnt = 0
for a, b, c, d in itertools.product(range(1, 13), repeat=4):
    v = chi4_direct(a, b, c, d)
    w = W4(a, b, c, d)
    cnt += 1
    exp = mp.mpf(0) if sum((a, b, c, d)) % 2 else mp.mpf(w) / 90
    if abs(v - exp) > mp.mpf('1e-40'):
        bad4 += 1
        if bad4 < 4:
            print(f"   反例 {(a, b, c, d)}  χ⁴={mp.nstr(v, 12)}  W={w}")
print(f"  全窮舉 n ≤ 12：{cnt} 組，違例 = {bad4} {'✅' if bad4 == 0 else '🔴'}")

random.seed(11)
bad5 = 0
for _ in range(20000):
    t = tuple(random.randint(1, 44) for _ in range(4))
    v = chi4_direct(*t)
    exp = mp.mpf(0) if sum(t) % 2 else mp.mpf(W4(*t)) / 90
    if abs(v - exp) > mp.mpf('1e-40'):
        bad5 += 1
        if bad5 < 4:
            print(f"   反例 {t}  χ⁴={mp.nstr(v, 12)}  W={W4(*t)}")
print(f"  隨機 20,000 組：違例 = {bad5} {'✅' if bad5 == 0 else '🔴'}")

# ---------------------------------------------------------------- §五
print("\n" + L)
print("【§五】V5 W 之值域（n ≤ 12 全窮舉）／V6 90χ⁴ ∈ ℤ")
print(L)
dist = {}
for a, b, c, d in itertools.product(range(1, 13), repeat=4):
    w = W4(a, b, c, d)
    dist[w] = dist.get(w, 0) + 1
print("  W 分佈（n ≤ 12）：", dict(sorted(dist.items())))
tot = sum(dist.values())
nz = tot - dist.get(0, 0)
print(f"  非零（四波耦合非零）比例 = {nz}/{tot} = {mp.nstr(mp.mpf(nz) / tot * 100, 4)}%")
bad6 = 0
for _ in range(2000):
    t = tuple(random.randint(1, 44) for _ in range(4))
    x = chi4_direct(*t) * 90
    if abs(x - mp.nint(x)) > mp.mpf('1e-40'):
        bad6 += 1
print(f"  90·χ⁴ 為整數（2,000 組）：違例 = {bad6} {'✅' if bad6 == 0 else '🔴'}")

# ---------------------------------------------------------------- §五b
print("\n" + L)
print("【§五b】例（數值落地，13 位）")
print(L)
for t in [(1, 2, 3, 4), (1, 1, 1, 1), (1, 2, 44, 43), (9, 15, 18, 20)]:
    v = chi4_direct(*t); w = W4(*t)
    print(f"  {str(t):>16}：W = {w:>2}  χ⁴ = {mp.nstr(v, 14)}")

# ---------------------------------------------------------------- §六
print("\n" + L)
print("【§六】判準與結論")
print(L)
print("""  ★ 定理（四波）
     χ⁴_{abcd} = 0                                   若 a+b+c+d 奇
     χ⁴_{abcd} = W/90                                若 a+b+c+d 偶
     其中 W = Σ_ε c_ε·[k_ε ≡ 0 (mod 90)] ∈ ℤ，
     k_ε 取 (±a±b±c±d) 之八個帶號和，c_ε = (+ + − − − − + +)

  ⟹ 非零判準：**存在一個帶號和 ∈ {0, ±90}**（|k| ≤ 176 ⟹ 只可能 0 或 ±90）
  ⟹ 四波耦合振幅量化於 1/90 網格（唔係連續可調）

  ★ 通用宇稱定理（任意 p）
     χ⁽ᵖ⁾ = 0 除非 Σ n_i ≡ p (mod 2)      （反射 j ↦ 45−j）
     p = 3 即 docs/37 之 K3 前半；p = 4 ⟹ 和偶才可能非零""")
print(L)
raise SystemExit(0)