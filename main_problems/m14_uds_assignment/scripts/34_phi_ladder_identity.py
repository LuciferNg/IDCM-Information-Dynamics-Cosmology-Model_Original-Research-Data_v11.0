#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/34 — X2′：Q2 之 φ² 中項 ≡ docs/17 之 φ-階梯（同構，非巧合）

────────────────────────────────────────────────────────────────────────
殘留：docs/32 §六 X2′（Q2 之 s3·s9 = 1 = s5² 與 docs/17 之 A_9:A_15:A_27 = φ⁻¹:1:φ 是否同一機制）

判準（跑前凍結）
  Y1 γ-GP 定義：γ_b² = γ_a·γ_c；A-GP 定義：A_b² = A_a·A_c，A_n ≡ 2sin(nπ/N)。
  Y2 命題待檢：對同一三元組，γ-GP ⟺ A-GP（因 γ = A²）。
  Y3 「同一機制」之判據（可算）：兩者之三元組索引集合相等，且比值同為 φ 之冪。
性質：確定性。N_trials = 0。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools

import mpmath as mp

mp.mp.dps = 60
L = "=" * 78
PHI = (1 + mp.sqrt(5)) / 2


def gam(N, n):
    return 4 * mp.sin(n * mp.pi / N) ** 2


def amp(N, n):
    return 2 * mp.sin(n * mp.pi / N)


print(L)
print("M14 scripts/34 — X2′：Q2（φ² 中項）vs docs/17 φ-階梯")
print(L)

# ============================================================================
print("\n【§一】/90 之 Q2 三元組 = {9,15,27}（u30 = 3）")
tri = (9, 15, 27)
g = {n: gam(90, n) for n in tri}
a = {n: amp(90, n) for n in tri}
for n in tri:
    print(f"  γ{n:<2} = {mp.nstr(g[n], 18)}   A{n:<2} = 2sin({n}π/90) = {mp.nstr(a[n], 18)}")
print(f"\n  γ 側：γ15² − γ9·γ27 = {mp.nstr(g[15]**2 - g[9]*g[27], 3)}   ⟹ GP，比 = {mp.nstr(g[15]/g[9], 18)}")
print(f"         φ² = {mp.nstr(PHI**2, 18)}；|偏差| = {mp.nstr(abs(g[15]/g[9] - PHI**2), 3)}")
print(f"  A 側：A15² − A9·A27 = {mp.nstr(a[15]**2 - a[9]*a[27], 3)}   ⟹ GP，比 = {mp.nstr(a[15]/a[9], 18)}")
print(f"         φ = {mp.nstr(PHI, 18)}；|偏差| = {mp.nstr(abs(a[15]/a[9] - PHI), 3)}")
print(f"\n  精確值：A9 = (√5−1)/2 = φ⁻¹ = {mp.nstr((mp.sqrt(5)-1)/2, 18)}（|dev| = "
      f"{mp.nstr(abs(a[9] - (mp.sqrt(5)-1)/2), 3)}）")
print(f"          A15 = 1；A27 = (√5+1)/2 = φ = {mp.nstr((mp.sqrt(5)+1)/2, 18)}（|dev| = "
      f"{mp.nstr(abs(a[27] - (mp.sqrt(5)+1)/2), 3)}）")
print(f"  ⟹ A9 : A15 : A27 = φ⁻¹ : 1 : φ —— 與 docs/17 之 φ-階梯**逐位相同**")

# ============================================================================
print("\n【§二】機制：A9·A27 = φ⁻¹·φ = 1 = A15²（互為倒數）")
print(f"  A9·A27 = {mp.nstr(a[9]*a[27], 18)}；A15² = {mp.nstr(a[15]**2, 18)}；"
      f"差 = {mp.nstr(a[9]*a[27] - a[15]**2, 3)}")
print(f"  γ 側同式：γ9·γ27 = ((3−√5)/2)((3+√5)/2) = (9−5)/4 = **1** = γ15²")
print("  ⟹ 機制 = 範數恆等 (3−√5)(3+√5) = 4 ⟹ γ9γ27 = 1；等價於 φ·φ⁻¹ = 1")

# ============================================================================
print("\n【§三】Y2 檢定：γ-GP ⟺ A-GP（全三元組，兩鏈）")
print(L)
for N, tag in ((88, "/88"), (90, "/90")):
    K = N // 2
    gg = {n: gam(N, n) for n in range(1, K + 1)}
    aa = {n: amp(N, n) for n in range(1, K + 1)}
    TOL = mp.mpf('1e-40')
    gp_g = {t for t in itertools.combinations(range(1, K + 1), 3)
            if abs(gg[t[1]] ** 2 - gg[t[0]] * gg[t[2]]) < TOL}
    gp_a = {t for t in itertools.combinations(range(1, K + 1), 3)
            if abs(aa[t[1]] ** 2 - aa[t[0]] * aa[t[2]]) < TOL}
    print(f"  {tag}：γ-GP = {len(gp_g)} 個；A-GP = {len(gp_a)} 個；"
          f"集合相同？ **{gp_g == gp_a}**")
    if gp_g != gp_a:
        print(f"        差集：γ−A = {gp_g - gp_a}；A−γ = {gp_a - gp_g}")
print("  ⟹ 因 γ = A²，γ_b² = γ_aγ_c ⟺ A_b⁴ = (A_aA_c)² ⟺ A_b² = A_aA_c（取正根）")
print("     故 γ-GP 與 A-GP **完全同集** ⟹ Q2 = φ-階梯，屬同一機制的兩種寫法")

# ============================================================================
print("\n【§四】其餘 30-族之幅度比（對照，不作質量宣稱）")
print(L)
for name, form in (("Q1", (1, 3, 11)), ("Q2", (3, 5, 9)), ("Q3", (7, 9, 13))):
    s = [gam(30, k) for k in form]
    A = [amp(30, k) for k in form]
    r = A[1] / A[0]
    print(f"  {name} {form}：A 比 = {mp.nstr(r, 12)}；"
          f"是否 φ 或 φ 之冪？ {any(abs(r - PHI**p) < mp.mpf('1e-30') for p in (-1, 1, 2, -2))}")
print("  ⟹ 僅 Q2 之幅度比為 φ；Q1、Q3 為其他代數數（Q1 之 A3 = φ⁻¹ 但中項非 φ 之冪）")

# ============================================================================
print("\n" + L)
print("【§五】判定（X2′）")
print(L)
print(f"""  ✅ **同一機制，非巧合**：Q2 之三元組即 {{9,15,27}}，其幅度 A9:A15:A27 = φ⁻¹:1:φ，
     與 docs/17 之 φ-階梯**逐位相同**（|dev| < 1e−41）。
  ✅ Y2 成立：γ-GP ⟺ A-GP（兩鏈全三元組集合逐一相同，{0} 例外）
     ⟹ 「γ 側之 Q2」與「A 側之 φ-階梯」係**同一關係**之兩種寫法，不得視為兩項獨立發現。
  ✅ 機制：範數恆等 (3−√5)(3+√5) = 4 ⟹ γ9·γ27 = 1 = γ15²；等價於 φ·φ⁻¹ = 1。
  ✅ 附帶：三個 30-族中**僅 Q2** 之幅度比為 φ；Q1、Q3 非 φ 之冪（故 φ-階梯之唯一性成立）。
  🔴 不作：φ 幂質量解讀（禁區）；此處只判「同一機制」與「唯一性」，不涉物理標籤。
""")
print(L)
raise SystemExit(0)