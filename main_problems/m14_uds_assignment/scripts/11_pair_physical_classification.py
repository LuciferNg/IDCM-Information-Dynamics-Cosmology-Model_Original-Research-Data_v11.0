#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/11 — 價對層之**物理**對稱分類（D_88 作用）＋ 「3 | 27」缺口檢定

────────────────────────────────────────────────────────────────────────
性質：(a) 以物理對稱群 D_88（旋轉 ℤ_88 ＋ 反射）作用於價對空間 {(a,b)}，求軌道分解；
      (b) 檢定 T_pipe 截斷段 n ≤ 27 是否容許 canonical 3 等分（docs/10 §八.2 之唯一缺口）。
      **零觀測量**（無 ρ、無 δ、無質量表）。
N_trials：0（確定性；無掃描、無擬合）

⚠ 與 scripts/09 之分別：09 用 Galois（**唔保 γ**，算術共軛，非物理對稱）；
   本輪用 D_88（**保 γ**，真對稱）⟹ 分類具物理意義。
────────────────────────────────────────────────────────────────────────
"""
import math
from collections import Counter
import mpmath as mp

mp.mp.dps = 30
L = "=" * 78
N, ND = 88, 44


def fold(r):
    r %= N
    return 0 if r == 0 else (r if r <= ND else N - r)


def gam(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print(L)
print("M14 scripts/11 — 價對之物理對稱分類 ＋ 3|27 缺口")
print(L)

pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]

# ======================================================================
print("\n【§一】旋轉 ℤ_88 作用於價對：σ_s(a,b) = (a+s, b+s)")
seen, orbs = set(), []
for p in pairs:
    if p in seen:
        continue
    o = set()
    for s in range(N):
        a, b = (p[0] + s) % N, (p[1] + s) % N
        if a == 0 or b == 0:
            break                      # 0 ∉ 模集
        o.add((a, b) if a < b else (b, a))
    seen |= o
    orbs.append(sorted(o))
print(f"  價對 946；旋轉軌道數 = {len(orbs)}")
ds = {min(o)[1] - min(o)[0] for o in orbs}
print(f"  每軌之 Δ = b−a：不同值數 = {len(ds)}；範圍 {min(ds)}..{max(ds)}")
print(f"  Δ 係**完備不變量**？{'是（每軌唯一 Δ）' if len(ds) == len(orbs) else '否'}")
print(f"  軌道大小分佈 = {dict(sorted(Counter(len(o) for o in orbs).items()))}")

# ======================================================================
print("\n【§二】加上反射（D_88 全群）：τ(a,b) = (88−b, 88−a)")
seen2, orbs2 = set(), []
for p in pairs:
    if p in seen2:
        continue
    o = set()
    for s in range(N):
        for sgn in (1, -1):
            a = (sgn * p[0] + s) % N
            b = (sgn * p[1] + s) % N
            if a == 0 or b == 0 or a == b:
                continue
            o.add((a, b) if a < b else (b, a))
    if o:
        seen2 |= o
        orbs2.append(sorted(o))
print(f"  D_88 軌道數 = {len(orbs2)}")
ds2 = {min(o)[1] - min(o)[0] for o in orbs2}
print(f"  每軌之 Δ：不同值數 = {len(ds2)}；範圍 {min(ds2)}..{max(ds2)}")
print(f"  反射保持 Δ？{'是（Δ 不變）' if ds2 == ds else '否'}")
print(f"  ⟹ 價對之 D_88-分類 = **按 beat 指數 Δ**，共 {len(orbs2)} 類；Σ 係軌道內座標")

# ======================================================================
print("\n【§三】T_pipe 截斷段 n ≤ 27：有無 canonical 3 等分？")
nstar = 27
seg = list(range(1, nstar + 1))
g = [gam(n) for n in seg]
inc = all(g[i] < g[i + 1] for i in range(len(g) - 1))
print(f"  γ_n 在 1..27 嚴格遞增？{inc}（γ_1 = {mp.nstr(g[0],8)} … γ_27 = {mp.nstr(g[-1],8)}）")
print("  ⟹ 任何**保 γ** 之模置換必為恆等（γ 全相異）⟹ 段內無非平凡對稱群")
# 是否存在保 γ 之非平凡平移（段內，mod 27）？
nontriv = [s for s in range(1, nstar) if all(abs(gam(((n - 1 + s) % nstar) + 1) - gam(n)) < mp.mpf('1e-20') for n in seg)]
print(f"  段內 mod-27 平移保 γ 者（非恆等）：{nontriv if nontriv else '無'}")
print(f"  3 | 27 ？{'是（27 = 3³）' if nstar % 3 == 0 else '否'}——但整除只給**算術可能**，非結構分割")
# n mod 3 在段內之三類，其 γ 和
t = {0: mp.mpf(0), 1: mp.mpf(0), 2: mp.mpf(0)}
for n in seg:
    t[n % 3] += gam(n)
print(f"  n mod 3 之三類 γ 和：{ {k: mp.nstr(v, 8) for k, v in t.items()} }（不相等 ⟹ 非對稱扇區）")
print(f"  ⟹ 段內 3 等分**不 canonical**：mod-3 係外加著色，唔係載體之對稱")

# ======================================================================
print("\n【§四】判定")
print(L)
print(f"""  ① 價對之物理分類（D_88）：**{len(orbs2)} 類，由 beat Δ = b−a 唯一標定**
     ⟹ 「味對型」之自然座標 = Δ ∈ {{1,…,{len(orbs2)}}}；Σ 係軌道內自由度（相位/和頻）
  ② 呢個分類**唔**給 3 類（43 ≠ 3）⟹ 「3 個味」唔喺價對層
  ③ 「3 | 27」缺口：段內 γ 全相異 ⟹ 無非平凡對稱；mod-3 係外加著色（γ 和不等）
     ⟹ **T1 之唯一缺口亦封**（3 | 27 只係算術巧合，非結構）
  ④ 綜合 docs/10：3 類來源在「模/站/價對」三層皆無 canonical 結構支撐；
     唯餘「計數」（B−1／9／18／27）—— 即「3」係**計數事實**，非**分類結構**
""")
print(L)
raise SystemExit(0)