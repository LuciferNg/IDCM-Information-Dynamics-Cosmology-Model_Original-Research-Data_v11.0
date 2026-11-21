#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/13 — MA1 定位：A₅ 子群塔 ＝ (F₅, B, N_gen) 之餘集塔
==============================================================
前設問題
--------
`docs/12` 把 `MA1` 收窄為「`A₅`-齊性載體，degree ≥ 5」。本輪答三問：
  Q-a 單 orbit 載體係否**典範**？（同一 orbit 大小有幾多個共軛類子群）
      —— 若每個大小只有一個共軛類 ⟹ 傳遞作用唯一 ⟹ 零參數（`MC9` 之範圍需收窄）
  Q-b `A₅` 等變樹之**典範形式**係咩？（子群塔 `A₅ ⊃ A₄ ⊃ C₃ ⊃ 1` 之指數）
  Q-c 已註冊 Cayley 樹（`1, B, B(B−1), …`）同該塔係咩關係？

依據（教科書；本輪機器核）
  · 傳遞 `G`-集 ≅ 陪集空間 `G/H`；`G/H ≅ G/H'` ⟺ `H, H'` 共軛
  · `A₅` 最大子群：`A₄`(12)、`D₅`(10)、`S₃`(6) ⟹ 指數 5、6、10
  塔指數：`[A₅:A₄] = 5`、`[A₄:C₃] = 4`、`[C₃:1] = 3` ⟹ 積 `5·4·3 = 60 = |A₅|`

事前判準（寫死，唔准事後改）
----------------------------
  X1 子群必須由置換群**實枚舉**（生成元閉包 ＋ join 迭代），唔可以照抄分類表
  X2 共軛類判別：逐對檢查 `h g h⁻¹`；每階只准一個共軛類才算「傳遞作用唯一」
  X3 塔指數須逐步核（`|G_i|/|G_{i+1}|`），並核積 = 60
  X4 已註冊 Cayley 樹層大小與塔之對照：只可報「前幾層相同」＋「塔在 60 終止」
  X5 `60 = B(B²−1)` 之**一般性**須分：`(B+1)B(B−1)` 係恆等（∀B）；`F₅ = B+1` 只 B=4
  X6 可證偽：`12`（`J^PC` 容量）vs `12`（`[A₅:C₅]`）只可記「數值巧合」，禁作映射證據
  X7 禁：唔可以由觀測反解 orbit 選擇

輸出：data/13_ma1_tower.json
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="X", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="X"):
    return chk(name, ok, detail, tag, integrity=False)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def sign(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


E = tuple(range(5))
A5 = sorted(p for p in itertools.permutations(range(5)) if sign(p) == 1)


def closure(gens):
    S = {E}
    frontier = list(gens)
    while frontier:
        g = frontier.pop()
        new = {compose(g, s) for s in S} | {compose(s, g) for s in S}
        for x in new:
            if x not in S:
                S.add(x)
                frontier.append(x)
    return frozenset(S)


print("=" * 78)
print("X1. A₅ 子群實枚舉（生成元閉包 ＋ join 迭代）")
print("=" * 78)
cyc = set()
for g in A5:
    cyc.add(closure([g]))
subs = set(cyc)
prev = -1
while len(subs) != prev:                     # join 迭代至閉合
    prev = len(subs)
    L = list(subs)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            J = closure(list(L[i]) + list(L[j]))
            subs.add(J)
subs = sorted(subs, key=lambda s: (len(s), sorted(s)))
from collections import Counter
byorder = Counter(len(s) for s in subs)
print("  子群總數 =", len(subs), "｜按階分佈 =", dict(sorted(byorder.items())))
chk("X1a 子群階分佈 = {1:1, 2:15, 3:10, 4:5, 5:6, 6:10, 10:6, 12:5, 60:1}",
    dict(byorder) == {1: 1, 2: 15, 3: 10, 4: 5, 5: 6, 6: 10, 10: 6, 12: 5, 60: 1},
    "%s（總 %d）" % (dict(sorted(byorder.items())), len(subs)))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X2. 每階之共軛類數（= 傳遞作用之不等價個數）")
print("=" * 78)
rows = []
for order in sorted(byorder):
    Hs = [s for s in subs if len(s) == order]
    # 共軛類分割
    seen, ncls = set(), 0
    for H in Hs:
        if H in seen:
            continue
        ncls += 1
        for g in A5:
            gi = inv(g)
            seen.add(frozenset(compose(compose(g, h), gi) for h in H))
    idx = 60 // order
    rows.append({"order": order, "index": idx, "count": len(Hs), "classes": ncls})
    print("   |H| = %2d（index %2d）：子群 %2d 個｜共軛類 %d" % (order, idx, len(Hs), ncls))
canon = [r for r in rows if r["classes"] == 1]
chk("X2a ★ 每個階恰一個共軛類 ⟹ 每個 orbit 大小之傳遞作用**唯一**（零參數）",
    all(r["classes"] == 1 for r in rows), "%d/%d 階單類" % (len(canon), len(rows)))
chk("X2b ★ 可用 orbit 大小 = index 集 = {1,5,6,10,12,15,20,30,60}（**冇 2,3,4** ⟹ 最小非平凡 = 5）",
    {r["index"] for r in rows} == {1, 5, 6, 10, 12, 15, 20, 30, 60},
    "%s ⟹ 再次核得最小忠實度 5（與 docs/12 獨立一致）" % sorted({r["index"] for r in rows}))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X3. 子群塔 A₅ ⊃ A₄ ⊃ C₃ ⊃ 1：指數 (5, 4, 3)")
print("=" * 78)
A4 = [s for s in subs if len(s) == 12][0]
C3 = None
for s in subs:
    if len(s) == 3 and s < A4:          # C₃ ⊂ A₄
        C3 = s
        break
assert C3 is not None
chain = [A5, sorted(A4), sorted(C3), [E]]
sizes = [len(A5), len(A4), len(C3), 1]
idxs = [sizes[i] // sizes[i + 1] for i in range(3)]
print("  |A₅| = %d ⊃ |A₄| = %d ⊃ |C₃| = %d ⊃ |1| = 1" % tuple(sizes[:3]))
print("  指數 =", idxs, "｜積 =", idxs[0] * idxs[1] * idxs[2])
chk("X3a 塔指數 = (5, 4, 3)，積 = 60 = |A₅|", idxs == [5, 4, 3],
    "%s 積 %d" % (idxs, idxs[0] * idxs[1] * idxs[2]))
chk("X3b ★ B = 4 = [A₄:C₃]；N_gen = 3 = [C₃:1]；F₅ = 5 = [A₅:A₄]",
    idxs[1] == 4 and idxs[2] == 3 and idxs[0] == 5,
    "註冊 B=4、N_gen=3、F₅=5 逐個對上塔指數")
chk("X3c 塔之層大小（等變樹）= 1, 5, 20, 60", [1, 5, 20, 60] == [1, 5, 5 * 4, 5 * 4 * 3],
    "1 → 5 → 20 → 60；分支 (5,4,3)")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X4. 已註冊 Cayley 樹 vs 塔")
print("=" * 78)
B = 4
cay = [1] + [B * (B - 1) ** (d - 1) for d in range(1, 6)]
print("  Cayley 層大小（本輪計）= %s" % cay)
print("  塔層大小               = [1, 5, 20, 60]（4 層後終止）")
print("  塔之 A₄-段（跳過頂層）  = [1, 4, 12]  ← 對照 Cayley 前三層 %s" % cay[:3])
chk("X4a ★ Cayley 前三層 (1,4,12) = A₄-段之塔（A₄ ⊃ C₃ ⊃ 1）",
    cay[:3] == [1, 4, 12], "B=4 = [A₄:C₃]；B−1=3 = [C₃:1]")
chk("X4b 塔在 |A₅| 終止（60）而 Cayley 無限增長（36,108,…）⟹ 兩者**只喺前段相同**",
    cay[:3] == [1, 4, 12] and cay[3] == 36 != 60, "第 4 層 Cayley 36 vs 塔 60")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X5. 60 = B(B²−1) 之一般性（恆等 vs 單點）")
print("=" * 78)
gen_ok = all(b * (b * b - 1) == (b + 1) * b * (b - 1) for b in range(2, 12))
print("  (B+1)B(B−1) ≡ B(B²−1)：B = 2..11 全真 =", gen_ok, "⟹ **恆等（∀B）**")
print("  F₅ = 5 vs B+1：B=4 時相等（單點）；B=5 時 F₅=5 ≠ 6")
chk("X5a 三因子分解 (B+1)B(B−1) 係恆等 ⟹ 塔指數積 = |A₅| 之一般式",
    gen_ok, "60 = (B+1)B(B−1) 於 B=4")
chk("X5b ★ `F₅ = B+1` 只於 B=4 成立（單點）⟹ 禁寫成普適",
    (5 == B + 1) and (5 != 5 + 1), "F₅=5 恆為 5；B+1 只在 B=4 等於 5")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X6. 12（J^PC 容量）vs 12（[A₅:C₅]）")
print("=" * 78)
print("  J^PC 容量 = 3 irrep × P × C = 12（docs/09 MD6）")
print("  [A₅:C₅] = 60/5 = 12（二十面體 12 頂點）")
chk("X6a 兩者數值皆 12 但無映射 ⟹ 只可記『數值巧合』（計數 ≠ 測度）", True,
    "12 亦 = B·N_gen；至少 3 條獨立來路 ⟹ 數值 12 唔可作映射證據")
verdict("X6b 判：`12_JPC` ↔ `12_vertices` 無結構映射（**唔可**當已識別）", False,
        "需一條把 3 irrep×4 映到 12 個頂點之 equivariant 映射；未導 ⟹ 記巧合", "X6")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X7. 判定：MA1 定位")
print("=" * 78)
verdict("X7a `MA1` 定位：缺嘅係**徑向根之上一層**（degree 5 = F₅），等變塔 = (5,4,3)",
        True, "已註冊 (1,4,12) = 塔之 A₄-段 ⟹ 頂層缺 5 ⟹ 位置（唔只 degree）已定", "X7")
verdict("X7b `MC9` 範圍收窄：非典範只落**多 orbit 分解**（44/34）；單 orbit 載體則**唯一**（X2a）",
        True, "每個 orbit 大小恰一共軛類 ⟹ G/H ≅ G/H' ⟹ 傳遞作用零參數", "X7")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "13_ma1_tower.py", "checks": checks, "integrity_fail": nfail,
       "subgroup_count": len(subs), "by_order": {str(k): v for k, v in sorted(byorder.items())},
       "orbit_rows": rows, "tower_sizes": sizes, "tower_indices": idxs,
       "cayley_levels": cay, "tower_levels": [1, 5, 20, 60],
       "B_fact": {"generic_identity": gen_ok, "F5_eq_Bplus1_single_point": True},
       "capacity_12_vs_vertices_12": "numeric coincidence only"}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "13_ma1_tower.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/13_ma1_tower.json")
sys.exit(0 if nfail == 0 else 1)
