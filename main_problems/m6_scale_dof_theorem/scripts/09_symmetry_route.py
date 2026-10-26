#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/09_symmetry_route.py — 升級路線①：對稱能否強迫「單一標度」？
================================================================================
問題：T2 現時係定義性（scripts/08 B）。能否用**對稱**把它變成真定理？
     即：晶格嘅自同構群會唔會強迫所有 sector 用同一個標度？

方法（精確、群論）：
  1. 作用群 G = (ℤ/88)^* 作用於 mode 類 {n, 88−n}（因 γ_n = γ_{88−n}）
  2. 計 G 喺 44 個類上嘅軌道分解 ⟹ 軌道數 = 可由對稱強迫嘅「獨立 sector」上限
  3. 逐 n 核：軌道大小 == deg γ_n = φ(88/gcd(n,88))/2（Galois 一致性）
  4. 判決：若軌道數 = 1 ⟹ 對稱強迫單一標度（升級成功）
           若軌道數 > 1 ⟹ 對稱只給出「≤ 軌道數」個 sector 標度 ⟹ 20→1 需**額外公設**

執行：/usr/bin/python3 scripts/09_symmetry_route.py
"""
import json
import os
import sys
from math import gcd

import sympy as sp

OUT = {"checks": {}}
NP = NF = 0


def ck(ok, msg, detail=""):
    global NP, NF
    NP, NF = (NP + 1, NF) if ok else (NP, NF + 1)
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}" + (f"  {detail}" if detail else ""))
    return bool(ok)


N = 88
UNITS = [a for a in range(1, N) if gcd(a, N) == 1]
print("=" * 92)
print("升級路線① — 對稱（Galois 軌道）能否強迫單一標度？")
print("=" * 92)
print(f"  G = (ℤ/88)^*，|G| = φ(88) = {len(UNITS)}；作用於類 {{n,88−n}}（γ_n = γ_{{88−n}}）")

# ── 1. 軌道分解 ────────────────────────────────────────
cls = lambda n: min(n % N, (-n) % N)          # 類代表（1..44）
def act(a, n):
    return cls(a * n)
seen, orbits = set(), []
for n in range(1, 45):
    if n in seen:
        continue
    orb, frontier = set(), [n]
    while frontier:
        x = frontier.pop()
        if x in orb:
            continue
        orb.add(x)
        for a in UNITS:
            y = act(a, x)
            if y not in orb:
                frontier.append(y)
    seen |= orb
    orbits.append(sorted(orb))
orbits.sort(key=lambda o: (len(o), o[0]))
print(f"  軌道分解：{len(orbits)} 個軌道，大小分佈 = "
      f"{sorted([len(o) for o in orbits], reverse=True)}")
ck(sum(len(o) for o in orbits) == 44, f"44 個類完整分割（共 {sum(len(o) for o in orbits)}）")

# ── 2. Galois 一致性：軌道大小 == deg γ_n ──────────────
print("\n  Galois 一致性（軌道大小 vs deg γ_n = φ(88/gcd(n,88))/2）：")
bad = []
for orb in orbits:
    n0 = orb[0]
    m = N // gcd(n0, N)
    deg = 1 if m <= 2 else sp.totient(m) // 2
    if deg != len(orb):
        bad.append((n0, len(orb), deg))
print(f"    {len(orbits)} 個軌道：不一致 {len(bad)} 個"
      + (f" ⟹ {bad}" if bad else " ⟹ 軌道大小 ≡ 代數度 ✅"))
ck(len(bad) == 0, "軌道大小 = 代數度（Galois 一致性成立）⟹ 軌道即 Galois 共軛類")

# 樣本
print("    樣本：" + "；".join(f"n={o[0]} 大小 {len(o)}" for o in orbits[:6]) + " …")

# ── 3. 判決 ────────────────────────────────────────────
print("\n  判決：")
n_orb = len(orbits)
max_orb = max(len(o) for o in orbits)
ck(n_orb == 7, f"軌道數 = {n_orb}（實測；先前誤預期 20 — 已修正）")
ck(max_orb == 20,
   f"最大軌道大小 = {max_orb} = rank_ℚ{{γ_n}}（M6 docs/01 A4）⟹ 一般 mode 之代數度",
   "⟹ rank 係『最大軌道』而非『軌道數』：20 ≠ 7，兩者係唔同量")
if n_orb == 1:
    ck(True, "對稱強迫單一標度 ⟹ T2 可升級為真定理")
else:
    ck(True,
       f"對稱只強迫【每個 Galois 軌道內部】標度相等 ⟹ sector 標度自由度上界 = {n_orb}，**唔係 1**",
       "⟹ 7 → 1 嘅收縮需要**額外公設**，對稱路線單獨唔夠（但已由 44 收窄到 7）")

# 邊際：軌道內共軛 ⟹ 要 ℚ-定義就必須等權（呢個係真定理）
print("\n  副產品（真定理，可用）：")
print("    γ_n 與 γ_m Galois 共軛（同軌道）⟹ 任何 ℚ-定義（有理性）之泛函")
print("    必須對同軌道之 mode 給予**相同權重**。此為 M6 之 T1（單一 M₀）提供了")
print("    一個**部分**理由：至少在任何單一軌道內部，「一個標度」係被強迫嘅。")
OUT["checks"]["orbits"] = {"n_orbits": n_orb,
                           "sizes": sorted([len(o) for o in orbits], reverse=True),
                           "sample": [[o[0], len(o)] for o in orbits[:8]],
                           "galois_consistency": len(bad) == 0,
                           "verdict": ("symmetry forces single scale" if n_orb == 1
                                       else f"symmetry bounds sectors at {n_orb}, not 1")}

print("\n" + "=" * 92)
print(f"路線①：{NP}/{NP + NF} PASS")
print("=" * 92)
OUT["pass"], OUT["fail"] = NP, NF
os.makedirs("data", exist_ok=True)
with open("data/09_symmetry_route.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/09_symmetry_route.json")
sys.exit(0 if NF == 0 else 1)