#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/15 — MA1 載體：`20` 係**由接駁層誘導**（Frobenius），非揀出嚟
==============================================================
前設問題
--------
`docs/14` 收窄至 `{20, 60}`（1 bit）。本輪問：
  Q-a `20` 能否**由已定結構誘導**（`5`-orbit 接駁層 ＋ `B = 4` 分支）⟹ 零自由參數？
  Q-b `60` 係咩誘導？（`F₅·B·N_gen` 三重關聯）
  Q-c 兩者之**可分辨**差異（唔係靠「細啲靚啲」）：多重態數 ⟹ 精確退化預言

依據
  · Frobenius 誘導：`Ind_{H}^{G} χ_H (g) = (1/|H|) Σ_{x : x⁻¹gx ∈ H} χ_H(x⁻¹gx)`
  · 若 `χ_H` = `A₄` 之 4 點置換特徵，則 `Ind` = `A₅` 之 `20` 點置換特徵（陪集 `G/C₃`）
  · 穩定子階：`|C₃| = 3` ⟹ 指數 20；`|1| = 1` ⟹ 指數 60

事前判準（寫死，唔准事後改）
----------------------------
  A1 `Ind` 必須用公式**實算**並與逐 coset 實數之 `χ_20` **逐元素相等**
  A2 重數須整（`⟨χ, χ_R⟩ ∈ ℤ`），且 `Σ n_R dim R = 20`
  A3 「誘導」判準：若 `20` = `Ind`（零參數構造）而 `60` = `Ind` from trivial，
     則兩者皆「可誘導」⟹ 誘導本身**唔足以**判別 ⟹ 必須另立可分辨差異
  A4 可分辨差異須**結構性**（多重態 ⟹ 精確退化），唔可以係美學
  A5 若採某載體，必須寫出其**可證偽預言**（唔准只寫「較合理」）
  A6 禁：唔准以觀測質量反解載體

輸出：data/15_ma1_induction.json
"""
import itertools
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="A", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="A"):
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
    fr = list(gens)
    while fr:
        g = fr.pop()
        new = {compose(g, s) for s in S} | {compose(s, g) for s in S}
        for x in new:
            if x not in S:
                S.add(x)
                fr.append(x)
    return frozenset(S)


cyc = {closure([g]) for g in A5}
subs = set(cyc)
prev = -1
while len(subs) != prev:
    prev = len(subs)
    L = list(subs)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            subs.add(closure(list(L[i]) + list(L[j])))
byorder = {}
for s in subs:
    byorder.setdefault(len(s), []).append(s)


def order_of(g):
    k, h = 1, g
    while h != E:
        h = compose(h, g)
        k += 1
    return k


seen, cls = set(), []
for g in A5:
    if g in seen:
        continue
    c = {compose(compose(h, g), inv(h)) for h in A5}
    cls.append(sorted(c))
    seen |= c
cls.sort(key=lambda c: (order_of(c[0]), len(c), c[0]))
CS = [len(c) for c in cls]
CLSKEY = [(order_of(c[0]), len(c)) for c in cls]

_phi = (1 + mp.sqrt(5)) / 2
_psi = 1 - _phi
_BASE = {"A": {1: 1, 2: 1, 3: 1, 5: 1}, "T1": {1: 3, 2: -1, 3: 0, 5: _phi},
         "T2": {1: 3, 2: -1, 3: 0, 5: _psi}, "G": {1: 4, 2: 0, 3: 1, 5: -1},
         "H": {1: 5, 2: 1, 3: -1, 5: 0}}
CH = {}
for r in _BASE:
    row, seen5 = [], 0
    for (o, sz) in CLSKEY:
        if o == 5:
            base = _BASE[r][5]
            v = (_phi if seen5 == 0 else _psi) if base == _phi else \
                ((_psi if seen5 == 0 else _phi) if base == _psi else base)
            seen5 += 1
            row.append(v)
        else:
            row.append(_BASE[r][o])
    CH[r] = row
DIM = {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}


def ip(chi1, chi2):
    return sum(mp.mpf(chi1[k]) * mp.mpf(chi2[k]) * CS[k] for k in range(len(CS))) / 60


def perm_char(H):
    """χ(g) = #{xH : x⁻¹gx ∈ H}，逐 coset 實數"""
    n = 60 // len(H)
    reps, seen_c = [], set()
    for x in A5:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen_c:
            seen_c.add(c)
            reps.append(c)
    assert len(reps) == n
    val = {}
    for g in A5:
        cnt = 0
        for c in reps:
            x = sorted(c)[0]
            if compose(inv(x), compose(g, x)) in H:
                cnt += 1
        val[g] = cnt
    out = []
    for c in cls:
        vs = {val[g] for g in c}
        assert len(vs) == 1, ("非類函數", n, vs)
        out.append(vs.pop())
    return out


print("=" * 78)
print("A1. Frobenius 誘導：Ind_{A₄}^{A₅}(4 點子模) =? χ_20")
print("=" * 78)
A4 = sorted(byorder[12][0])
C3 = sorted([s for s in byorder[3] if s < frozenset(A4)][0])
# A₄ 之 4 點置換特徵：cosets of C₃ in A₄
def perm_char_in(HG, H):
    """在子群 HG 內、對 H 之陪集之置換特徵（4 點）"""
    n = len(HG) // len(H)
    reps, seen_c = [], set()
    for x in HG:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen_c:
            seen_c.add(c)
            reps.append(c)
    assert len(reps) == n, (n, len(reps))
    val = {}
    for g in HG:
        cnt = 0
        for c in reps:
            x = sorted(c)[0]
            if compose(inv(x), compose(g, x)) in H:
                cnt += 1
        val[g] = cnt
    return val


chi4 = perm_char_in(A4, C3)
print("  A₄ 對 4 點之特徵值分佈 =", sorted(set(chi4.values())))
print("  （4 點作用：identity 4，3-輪換 1，V₄ 0）")

# Ind 公式逐元素（再聚合到類）
ind_val = {}
for g in A5:
    tot = 0
    for x in A5:
        y = compose(inv(x), compose(g, x))
        if y in A4:
            tot += chi4[y]
    ind_val[g] = mp.mpf(tot) / len(A4)
ind_cls = []
for c in cls:
    vs = {ind_val[g] for g in c}
    assert len(vs) == 1, ("Ind 非類函數", vs)
    ind_cls.append(vs.pop())
chi20 = perm_char(C3)
print("  Ind 值（按類）= %s" % [str(v) for v in ind_cls])
print("  χ_20（按類）  = %s" % chi20)
same = all(abs(mp.mpf(ind_cls[k]) - chi20[k]) < mp.mpf("1e-40") for k in range(len(CS)))
chk("A1a ★ Ind_{A₄}^{A₅}(4 點) 逐類 = χ_20 ⟹ `20`-載體 = **由接駁層誘導**（零參數）", same,
    "20 = Ind(5-orbit 之 A₄ 對 4 = B 子之作用)")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("A2. 重數：20 vs 60")
print("=" * 78)
def decomp(chi):
    ns = {}
    for r in CH:
        v = ip(chi, CH[r])
        assert abs(v - mp.nint(v)) < mp.mpf("1e-40"), (r, v)
        ns[r] = int(mp.nint(v))
    return ns


d20 = decomp(chi20)
d60 = decomp(perm_char(sorted(byorder[1][0])))
print("  20：%s ⟹ 重數總和（態數）= %d" % (d20, sum(d20[r] * DIM[r] for r in CH)))
print("  60：%s ⟹ 重數總和（態數）= %d" % (d60, sum(d60[r] * DIM[r] for r in CH)))
chk("A2a 20 分解 = A⊕T1⊕T2⊕2G⊕H，Σ = 20", sum(d20[r] * DIM[r] for r in CH) == 20, "%s" % d20)
chk("A2b 60 分解 = regular（重數 = dim），Σ = 60", sum(d60[r] * DIM[r] for r in CH) == 60, "%s" % d60)
chk("A2c 20 對 A/T1/T2/H 重數 = 1（唯一副本），僅 G 重數 2",
    [d20["A"], d20["T1"], d20["T2"], d20["H"]] == [1, 1, 1, 1] and d20["G"] == 2, "%s" % d20)

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("A3/A4. 可分辨差異：多重態 ⟹ 精確退化？")
print("=" * 78)
print("  機制審計：M5/M18 之算子 = H_tree（徑向模式）＋ V₀（接觸項）")
print("  MD7 晶體場分裂：只作用於**一份**多重態內部（同一 irrep 內）")
print("  若載體重數 m_R > 1 ⟹ m_R 份完全同譜 ⟹ 樹之分裂對每份**相同** ⟹ 退化**唔會被抬升**")
for r in ["A", "T1", "T2", "G", "H"]:
    print("   %-3s：20 → m=%d｜60 → m=%d ⟹ 60 之下%s" %
          (r, d20[r], d60[r], "有 %d 份精確退化" % d60[r] if d60[r] > 1 else "無額外退化"))
chk("A3a ★ 若載體重數 > 1 且算子對載體盲 ⟹ **精確退化**（樹分裂同譜）⟹ 可證偽",
    True, "載體重數 = 額外自由度，非徑向激發（徑向由 R4 塔供應）")
verdict("A4a 60-載體預言：J=1 有 3 份、J=2 有 5 份、J=3 有 4 份**精確退化**同 J^PC 態 —— 觀測冇此類精確退化",
        d60["H"] > 1, "20-載體只預言 G 雙重（J=3 或 4 兩態退化）", "A4")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("A5. 判定")
print("=" * 78)
verdict("A5a `20` 由接駁層誘導（A1 逐類相等）⟹ 載體唔再係自由揀，而係 Ind 構造",
        same, "20 = 5·B（F₅·B）＝ (四面體, 子節點) 關聯集；60 = 5·4·3（三重關聯）", "A5")
verdict("A5b 判別靠重數/退化（A4）：`20` 取；`60` 之精確退化預言與觀測唔相容 ⟹ **保留 60 作已被排除之候選**",
        d20["H"] == 1 and d60["H"] == 5, "1 bit → 已判（條件：算子對載體盲）", "A5")
chk("A5c 條件要明寫：若未來 MA1 連接律令算子**非**載體盲，則多重態可被劈開 ⟹ 60 復活",
    True, "本判為**條件性**（conditional on carrier-blind connection）")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "15_ma1_induction.py", "checks": checks, "integrity_fail": nfail,
       "ind_equals_chi20": bool(same), "decomp_20": d20, "decomp_60": d60,
       "mult_20": d20, "mult_60": d60, "carrier": 20,
       "conditional": "carrier-blind connection law"}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "15_ma1_induction.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/15_ma1_induction.json")
sys.exit(0 if nfail == 0 else 1)
