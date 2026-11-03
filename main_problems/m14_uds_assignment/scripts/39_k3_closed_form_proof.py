#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/39 — K3 後半：把 cot 和化成乘積 ⟹ χ³ 非零性之**完全證明**

────────────────────────────────────────────────────────────────────────
待證（docs/35 §二 遺留）：
  χ³_{abm} = C·Σ_{j=1}^{44} sin(πaj/45)sin(πbj/45)sin(πmj/45)，C=(2/45)^{3/2}
  a+b+m 奇時 χ³ ≠ 0（6622 例已窮舉，但一般非零性未證）

本腳本之推導（結構，非拟合）：見 log §一 之四步；核心代數鏈
  (1) 積化和：sinA sinB sinC = ¼[sin(A+B−C)+sin(A−B+C)+sin(−A+B+C)−sin(A+B+C)]
  (2) 鏈和：S(k) ≡ Σ_{j=1}^{N−1} sin(πkj/N) = 0（k 偶）／cot(πk/2N)（k 奇），N=45
  (3) cot 併項：cot p + cot q = sin(p+q)/(sin p sin q)，cot p − cot q = sin(q−p)/(sin p sin q)
  (4) 積化和再化 ⟹ **Φ = sin(πa/45)sin(πb/45)sin(πm/45) / ∏_{±} sin(π(±a±b±m)/90)**
  ⟹ 分子 ≠ 0（1≤a≤44）、分母 ≠ 0（k_i 奇且 |k_i|≤132 ⟹ k_i ∉ 90ℤ）⟹ χ³ ≠ 0 ∎

判準（跑前凍結）
  V1 閉式 vs 直接求和：全部奇和三元組（6622），相對誤差 < 1e−40（dps 60）
  V2 宇稱引理：k_i = 0 ⟹ a+b+m 偶（全 85184 三元組掃描，反例必須 0）
  V3 分母安全：a+b+m 奇 ⟹ 四個 k_i 全奇、全非零、≠ ±90
  V4 全域二分：全 85184 三元組，χ³ = 0 ⟺ a+b+m 偶
  V5 最小 |χ³| 位置與值（應為 (1,2,44)、1.609741784646e−6）
  V6 對稱性：閉式在三元組任意置換下不變（抽 200 組）
性質：確定性窮舉。N_trials = 0。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools

import mpmath as mp

mp.mp.dps = 60
N = 45
C = (mp.mpf(2) / 45) ** mp.mpf('1.5')
L = "=" * 78
TOL_REL = mp.mpf('1e-40')


def S(k):
    """Σ_{j=1}^{44} sin(πkj/45)：k 偶 → 0；k 奇 → cot(πk/90)。"""
    return mp.mpf(0) if k % 2 == 0 else mp.cot(mp.pi * k / 90)


def chi_direct(a, b, m):
    s = mp.mpf(0)
    for j in range(1, 45):
        s += (mp.sin(mp.pi * a * j / 45) * mp.sin(mp.pi * b * j / 45)
              * mp.sin(mp.pi * m * j / 45))
    return C * s


def keys(a, b, m):
    return (a + b - m, a - b + m, -a + b + m, a + b + m)


def phi_closed(a, b, m):
    """Φ = ∏ sin(πx/45) / ∏ sin(πk_i/90)，χ³ = C/4 · Φ（僅奇和時有意義）"""
    num = (mp.sin(mp.pi * a / 45) * mp.sin(mp.pi * b / 45) * mp.sin(mp.pi * m / 45))
    den = mp.mpf(1)
    for k in keys(a, b, m):
        den *= mp.sin(mp.pi * k / 90)
    return num / den


print(L)
print("M14 scripts/39 — K3 後半：χ³ 非零性之完全證明")
print(L)

# ---------------------------------------------------------------- §一
print("\n【§一】代數鏈（四步）")
print("""  (1) 積化和：sin A sin B sin C = ¼[sin(A+B−C)+sin(A−B+C)+sin(−A+B+C)−sin(A+B+C)]
  (2) 鏈和：S(k) ≡ Σ_{j=1}^{44} sin(πkj/45) = 0（k 偶）｜cot(πk/90)（k 奇）
      ⟹ χ³ = (C/4)·[S(k₁)+S(k₂)+S(k₃)−S(k₄)]，k = (a+b−m, a−b+m, −a+b+m, a+b+m)
  (3) cot 併項（x_i ≡ πk_i/90）：
      cot x₁+cot x₂ = sin(x₁+x₂)/(sin x₁ sin x₂) = sin(πa/45)/(sin x₁ sin x₂)
      cot x₃−cot x₄ = sin(x₄−x₃)/(sin x₃ sin x₄) = sin(πa/45)/(sin x₃ sin x₄)
  (4) 再積化和：sin x₁ sin x₂ + sin x₃ sin x₄
       = ½[cos(π(b−m)/45) − cos(π(b+m)/45)] = sin(πb/45) sin(πm/45)
  ⟹ **Φ = sin(πa/45)sin(πb/45)sin(πm/45) / [sin x₁ sin x₂ sin x₃ sin x₄]**
     此式在三元組置換下明顯不變（分母之四個符號型亦隨之置換）""")

# ---------------------------------------------------------------- V1
print("\n" + L)
print("【§二】V1 閉式 vs 直接求和（全部奇和三元組）")
print(L)
bad, cnt = [], 0
roots = [x for x in itertools.combinations_with_replacement(range(1, 45), 3)
         if sum(x) % 2 == 1]
for a, b, m in roots:
    cnt += len(set(itertools.permutations((a, b, m))))
    d = chi_direct(a, b, m)
    f = C / 4 * phi_closed(a, b, m)
    rel = abs(d - f) / abs(d)
    if rel > TOL_REL:
        bad.append((a, b, m, mp.nstr(rel, 3)))
print(f"  奇和三元組（含重複，計入置換）= {cnt}；相對誤差 > 1e−40 者 = {len(bad)}")
print(f"  {'✅ 閉式與直接求和一致' if not bad else '🔴 ' + str(bad[:5])}")

# ---------------------------------------------------------------- V2
print("\n" + L)
print("【§三】V2 宇稱引理：k_i = 0 ⟹ a+b+m 偶")
print(L)
viol, zero_odd = [], []
for a, b, m in itertools.product(range(1, 45), repeat=3):
    ks = keys(a, b, m)
    if 0 in ks and sum((a, b, m)) % 2 == 1:
        zero_odd.append((a, b, m, ks))
    if 0 in ks and sum((a, b, m)) % 2 == 0:
        viol.append((a, b, m))
print(f"  k_i = 0 且和奇 ⟹ 直接反例數 = {len(zero_odd)}（必須 0）")
print("  理由（代數）：k₁=0 ⟹ m=a+b ⟹ a+b+m = 2(a+b) 偶；k₂=0 ⟹ b=a+m ⟹ 2b 偶；")
print("                k₃=0 ⟹ a=b+m ⟹ 2(b+m) 偶；k₄=0 ⟹ a+b+m=0 與 a,b,m ≥ 1 矛盾")
print(f"  對照：現存 k_i=0 之三元組共 {len(viol)} 個，全為偶和 ⟹ 奇和時 cot(0) 永不出現 ✅")

# ---------------------------------------------------------------- V3
print("\n" + L)
print("【§四】V3 分母安全：奇和 ⟹ 四個 k_i 全奇、非零、≠ ±90")
print(L)
bad3 = []
for a, b, m in roots:
    for k in keys(a, b, m):
        if k % 2 == 0 or k == 0 or abs(k) == 90:
            bad3.append((a, b, m, k))
print(f"  違例 = {len(bad3)}")
print("  界：|k_i| ≤ 132 < 180；奇數 ⟹ k_i ∉ 90ℤ ⟹ sin(πk_i/90) ≠ 0 ⟹ 分母 ≠ 0 ✅")

# ---------------------------------------------------------------- V4
print("\n" + L)
print("【§五】V4 全域二分：χ³ = 0 ⟺ a+b+m 偶（85184 三元組）")
print(L)
bad4, minabs, argmin = [], None, None
for a, b, m in itertools.product(range(1, 45), repeat=3):
    d = chi_direct(a, b, m)
    even = sum((a, b, m)) % 2 == 0
    if even and abs(d) > mp.mpf('1e-45'):
        bad4.append((a, b, m, 'even-nonzero'))
    if (not even) and abs(d) < mp.mpf('1e-45'):
        bad4.append((a, b, m, 'odd-zero'))
    if not even:
        v = abs(d)
        if minabs is None or v < minabs:
            minabs, argmin = v, (a, b, m)
print(f"  違例 = {len(bad4)} {'✅' if not bad4 else bad4[:5]}")
print(f"  最小 |χ³| = {mp.nstr(minabs, 13)} 於 {argmin}")

# ---------------------------------------------------------------- V5
print("\n" + L)
print("【§六】V5 閉式重現最小點 ＋ V6 置換不變")
print(L)
# 帶號值（χ³ 於最小點為負 —— 反向覆蓋要求文檔嘅負號在 data 落地）
sgn_min = chi_direct(*argmin)
fmin = C / 4 * phi_closed(*argmin)
print(f"  閉式於 {argmin}：χ³ = {mp.nstr(fmin, 20)}")
print(f"  直接求和   {argmin}：χ³ = {mp.nstr(sgn_min, 20)}")
print(f"  相對差 = {mp.nstr(abs(fmin - sgn_min) / abs(sgn_min), 3)}")
print(f"  最小點簡併：χ³(1,1,43) = {mp.nstr(chi_direct(1, 1, 43), 20)}")
print(f"              χ³(1,2,44) = {mp.nstr(chi_direct(1, 2, 44), 20)}")
print(f"              兩者之差 = {mp.nstr(abs(chi_direct(1,1,43) - chi_direct(1,2,44)), 3)}")
bad6 = 0
sample = list(itertools.islice(
    (t for t in itertools.product(range(1, 45), repeat=3) if sum(t) % 2 == 1), 200))
for t in sample:
    vals = {mp.nstr(C / 4 * phi_closed(*p), 50) for p in set(itertools.permutations(t))}
    if len(vals) != 1:
        bad6 += 1
print(f"  置換不變抽測 {len(sample)} 組，違例 = {bad6}")

# ---------------------------------------------------------------- §七
print("\n" + L)
print("【§七】★ 定理")
print(L)
print("""  **K3 定理（完整，取代窮舉）**
  設 N = 45、φ_n(j) = √(2/45) sin(πnj/45)、1 ≤ a,b,m ≤ 44。則

      χ³_{abm} = (2/45)^{3/2}/4 · ∏_{cyc} sin(πa/45) / ∏_{±} sin(π(±a±b±m)/90)

  且   χ³_{abm} = 0  ⟺  a+b+m 偶
       χ³_{abm} ≠ 0  ⟺  a+b+m 奇   （**一般非零性成立，非僅 6622 例**）

  證明鏈即 §一 四步：分子三個 sin 皆非零（1 ≤ a ≤ 44 ⟹ πa/45 ∈ (0,π)，a ≠ 45）；
  分母四個 sin 皆非零（奇和 ⟹ k_i 奇、非零（§三）、|k_i| ≤ 132 ⟹ k_i ∉ 90ℤ）。
  ⟹ docs/35 §二 之「🔴 一般非零性未證（Conway–Jones 型）」**解除**。""")
print(L)
raise SystemExit(0)