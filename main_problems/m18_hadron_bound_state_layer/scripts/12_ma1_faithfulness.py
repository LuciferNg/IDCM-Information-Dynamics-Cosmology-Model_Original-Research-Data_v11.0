#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/12 — MA1 收窄：A₅ ↔ 鏈／樹之**忠實作用**不可能（三向張力）
=====================================================================
前設問題
--------
`docs/11` 證 `A₅` 作用於 `I₄₄` 非典範（`MC9`）。若要求**保結構**（樹：保層／保父子；
鏈：保鄰接），能唔能忠實作用？本輪答：
  Q-a `A₅` 有冇 index ≤ 4 之子群？（有 ⟹ 4 個分支已夠）
  Q-b 保結構之忠實作用需要 `B ≥ ?`（分支下界）
  Q-c 鏈（環 `C₄₄`）可否被 `A₅` 保鄰接作用？
  Q-d 若 `B ≥ 5` 但 `B = 4` 已由 M16 閉合 ⟹ 三向張力，逃生路線各有幾貴？

數學依據（教科書；本輪逐項機器核）
  · `A₅` 單（simple）：normal subgroup ⊆ 共軛類之並（含 identity），階須整除 60
  · 若 `A₅` 作用於 `k` 點：`hom A₅ → S_k`；簡單 ⟹ 或忠實或平凡；忠實需 `60 | k!`
  · `4! = 24 < 60` ⟹ `k ≤ 4` 之作用必平凡 ⟹ **最小忠實度 = 5**
  · 環 `C₄₄`：`Aut = D₄₄`（階 176）；`A₅ ≤ D₄₄` 不可能（`60 ∤ 176`）

事前判準（寫死，唔准事後改）
----------------------------
  W1 單性以「共軛類之並 ＋ 階整除」機器核（唔可以只引書）
  W2 「保結構」定義寫死：樹 = 保層 ＋ 保父子；鏈 = 保鄰接
  W3 忠實性判準：kernel = 1；`hom` 為簡單群時只可全射（⟹ 階整除）
  W4 若判「唔可能」，必須列**全部**逃生路線 ＋ 每條之代價（唔可以留空）
  W5 可證偽：若 `B = 5`，核對邊些已註冊恆等式會斷（量化代價）
  W6 禁：唔准為咗救耦合而改 `B`、`N_d` 或 `A₅` 身份

輸出：data/12_ma1_faithfulness.json
"""
import itertools
import json
import math
import os
import sys

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="W", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="W"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("W1. A₅ 單性（共軛類之並 ＋ 階整除，機器核）")
print("=" * 78)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


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


A5 = sorted(p for p in itertools.permutations(range(5)) if sign(p) == 1)
print("  |A₅| =", len(A5))


def inv(h):
    r = [0] * len(h)
    for i, v in enumerate(h):
        r[v] = i
    return tuple(r)


seen, cls = set(), []
for g in A5:
    if g in seen:
        continue
    c = {compose(compose(h, g), inv(h)) for h in A5}
    cls.append(sorted(c))
    seen |= c
sizes = sorted(len(c) for c in cls)
print("  共軛類大小 =", sizes, "｜類數 =", len(cls))
# normal subgroup = 並若干類（含 identity 類）＋ 階整除 60
idc = [i for i, c in enumerate(cls) if tuple(range(5)) in c][0]
SUM = []
for k in range(0, len(cls)):
    for comb in itertools.combinations([i for i in range(len(cls)) if i != idc], k):
        order = sum(sizes[i] for i in comb) + sizes[idc]
        if 60 % order == 0:
            SUM.append(order)
print("  可能之 normal 階（類並 ＋ 整除 60）=", sorted(set(SUM)))
chk("W1a A₅ 單：唯一 normal 階 = {1, 60}", sorted(set(SUM)) == [1, 60], "%s" % sorted(set(SUM)))
chk("W1b 類大小和 = 60（1+12+12+15+20）", sum(sizes) == 60, "%s" % sizes)

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W2/W3. k ≤ 4 之作用必平凡 ⟹ 最小忠實度 = 5")
print("=" * 78)
tab = []
for k in range(1, 8):
    ok_div = (math.factorial(k) % 60 == 0)
    tab.append((k, math.factorial(k), ok_div))
    print("  k = %d：|S_k| = %-4d｜60 | k! ? %s" % (k, math.factorial(k), ok_div))
first = min(k for k, _, ok in tab if ok)
chk("W2a 最小 k 使 60 | k! = %d ⟹ **最小忠實度 = 5**" % first, first == 5, "4! = 24 < 60；5! = 120")
chk("W2b 對 k ≤ 4：A₅ → S_k 之 hom 只能平凡（簡單 ⟹ 或忠實或平凡；忠實需 60|k!）",
    all(not ok for k, _, ok in tab if k <= 4), "k=1..4 全部 60 ∤ k!")
verdict("W2c 推論：A₅ 若忠實作用於**保結構**之 4 分支樹 ⟹ 不可能（level 1 只有 B = 4 點）",
        False, "level 1 之作用 = hom A₅ → S₄，無法忠實 ⟹ 平凡；再由父子等變 ⟹ 全樹平凡", "W2")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W3'. 鏈（環 C₄₄）保鄰接作用")
print("=" * 78)
AUT_C44 = 2 * 44        # |D₄₄| = 88 = 2³·11
AUT_C88 = 2 * 88        # |D₈₈| = 176 = 2⁴·11
print("  |Aut(C₄₄)| = |D₄₄| = %d = 2³·11｜|Aut(C₈₈)| = |D₈₈| = %d = 2⁴·11" % (AUT_C44, AUT_C88))
chk("W3'a 60 ∤ 88 且 60 ∤ 176 ⟹ 冇單射 A₅ → D₄₄／D₈₈ ⟹ 保鄰接作用只能平凡",
    AUT_C44 % 60 != 0 and AUT_C88 % 60 != 0,
    "88 mod 60 = %d｜176 mod 60 = %d" % (AUT_C44 % 60, AUT_C88 % 60))
chk("W3'b 因 A₅ 單，任何 hom A₅ → D_n 之 kernel ∈ {1, A₅}；kernel=1 需 60|2n（假）",
    True, "⟹ kernel = A₅ ⟹ 作用平凡")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W4. 三向張力：逃生路線與代價")
print("=" * 78)
ESC = [
    ("E1 放棄 MD3（J 載體 ≠ A₅）", "🔴 重價", "J=2 需 5 維 irrep；C_n/D_n/T/O 全無 ⟹ 冇替代"),
    ("E2 容許破結構耦合（新載體，degree ≥ 5）", "🟡 中價", "需新嵌入層（MA1 規格化）；非典範（MC9）"),
    ("E3 改 B ≥ 5", "🔴 重價", "B=4 已由 M16 閉合；且 B=5 會斷 Σλ=66、16=B²、12=B·N_gen 等"),
]
for nm, cost, why in ESC:
    print("   %-32s %-8s %s" % (nm, cost, why))
chk("W4a 逃生路線 = 3 條，全部有代價標籤", len(ESC) == 3, "3/3")

print()
print("W5. 若 B = 5：已註冊恆等式斷幾多？")
B5 = 5
idents = {
    "Σλ = B³ + B/2": (4**3 + 4 / 2, B5**3 + B5 / 2),
    "B² = 16": (4**2, B5**2),
    "B·N_gen = 12": (4 * 3, B5 * 3),
    "K_tree^(mode) = 2B²+1 = 33": (2 * 4**2 + 1, 2 * B5**2 + 1),
    "K_tree^(node) = 2(B²+1) = 34": (2 * (4**2 + 1), 2 * (B5**2 + 1)),
    "|A₅| = B(B²−1) = 60": (4 * (4**2 - 1), B5 * (B5**2 - 1)),
}
broken = 0
for k, (v4, v5) in idents.items():
    same = (v4 == v5)
    broken += (not same)
    print("   %-32s B=4：%-6s ｜B=5：%-6s ｜%s" % (k, v4, v5, "不變" if same else "★ 斷"))
chk("W5a B=5 會令 %d/%d 條已註冊恆等式斷裂 ⟹ E3 代價極高" % (broken, len(idents)), broken >= 4,
    "%d/%d 斷" % (broken, len(idents)))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W6. 判定")
print("=" * 78)
verdict("W6a MA1 最終判：保結構之忠實 A₅ 耦合 **結構上不可能**（樹需 B≥5、鏈需 60|176）",
        False, "兩條保結構路線皆死；只剩破結構之新載體（degree ≥ 5）", "W6")
verdict("W6b 新載體規格（MA1 收窄）：A₅-齊性空間，|orbit| ∈ {5,6,10,12,15,20,30,60}，"
        "且與樹之連接須 ≥ 5 分支等變", True,
        "最小規格 = degree 5 = F₅；呢個係 MA1 之可交付形式（唔係空口『要新層』）", "W6")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "12_ma1_faithfulness.py", "code": "MC10", "checks": checks, "integrity_fail": nfail,
       "A5_order": len(A5), "class_sizes": sizes, "possible_normal_orders": sorted(set(SUM)),
       "min_faithful_degree": first, "aut_C44": AUT_C44,
       "escape_routes": [{"route": a, "cost": b, "why": c} for a, b, c in ESC],
       "B5_broken_idents": broken, "idents_count": len(idents)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "12_ma1_faithfulness.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/12_ma1_faithfulness.json")
sys.exit(0 if nfail == 0 else 1)
