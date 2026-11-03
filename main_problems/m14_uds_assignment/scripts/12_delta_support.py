#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/12 — E1：Δ 之有效支撐集（由已凍結之 P1/P2/P3/P5 約束給出）

────────────────────────────────────────────────────────────────────────
性質：對 946 個價對施加**已凍結**清單（docs/10 §一）之各組合約束，
      求「存活之 Δ 值集」及其基數。目的：檢查有無任何**已註冊**約束組合
      令 Δ 之支撐集恰為 3 個（即「3 個味」之最後可能來源）。
      零觀測量。N_trials = 7（7 個約束組合，全報）。
⚠ 若某組合得出 3，屬**後驗**（我已知「3」係目標）⟹ 只可作假設，不可作命中。
────────────────────────────────────────────────────────────────────────
"""
import math
from collections import Counter

N, ND, NSTAR = 88, 44, 27
L = "=" * 78


def gc(n):
    return math.gcd(n, N)


pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]
print(L)
print("M14 scripts/12 — E1：Δ 之有效支撐集（凍結約束組合）")
print(L)
print(f"  價對總數 = {len(pairs)}；Δ 範圍 = 1..{ND-1}")

FILTERS = [
    ("F0 無約束", lambda a, b: True),
    ("F1 P1 同 Galois 類（gcd 相等）", lambda a, b: gc(a) == gc(b)),
    ("F2 P2 同 mod-4 扇區（88 環）", lambda a, b: (a % 4) == (b % 4)),
    ("F3 P3 同奇偶", lambda a, b: (a - b) % 2 == 0),
    ("F4 P5 兩模皆 ≤ 27（T_pipe 投影段）", lambda a, b: b <= NSTAR),
    ("F5 F1 ∧ F4", lambda a, b: gc(a) == gc(b) and b <= NSTAR),
    ("F6 F2 ∧ F4", lambda a, b: (a % 4) == (b % 4) and b <= NSTAR),
]

rows = []
for name, f in FILTERS:
    live = [(a, b) for (a, b) in pairs if f(a, b)]
    ds = sorted({b - a for (a, b) in live})
    rows.append((name, len(live), ds))
    show = ds if len(ds) <= 30 else ds[:30] + ["…"]
    print(f"\n  {name}")
    print(f"    存活價對 = {len(live)}；Δ 集 = {len(ds)} 個：{show}")

print("\n" + L)
print("【判定】")
print(L)
print("  N_trials = 7（全報）。Δ 支撐集基數一覽：")
for name, cnt, ds in rows:
    flag = "★ 3" if len(ds) == 3 else ""
    print(f"    {name:<40s} |Δ| = {len(ds):>2d}  {flag}")
hit = [r for r in rows if len(r[2]) == 3]
if not hit:
    print("\n  ⟹ **無任何凍結約束組合令 Δ 支撐集 = 3**。")
    print("  ⟹ E2（「3」可否由 Δ 支撐給出）判 🔴：唔可以。")
else:
    print(f"\n  ⚠ 有 {len(hit)} 個組合得出 3 —— 屬後驗，只可作假設。")
print("\n  與 docs/10／docs/11 合併之總判定：")
print("    模層 1 類、站層 43 類、價對層 43 類（Δ 完備）、截斷段 1 類、")
print("    Δ 支撐集（7 個約束組合）皆 ≠ 3")
print("  ⟹ 「3」在 IDCM 載體上係**計數事實**（B−1 / 9 / 18 / 27），非分類結構。")
print(L)
raise SystemExit(0)