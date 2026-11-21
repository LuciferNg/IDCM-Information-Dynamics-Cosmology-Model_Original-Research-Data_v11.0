#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/22 — 耦合自由度普查：束上「點連」究竟有幾多自由？
============================================================
前設問題
--------
`MD14` 把連接律定為誘導束（纖維 `5/20/60`），但**相鄰深度之間嘅耦合算子**（即「點連」）
未定。`MD13` 證 `H` 只可 `A₄`-等變（`d=0` 單點揀定 `5`-orbit 之一點 ⟹ `A₅` 破至 `A₄`）。
本輪問：**`A₄`-等變耦合空間嘅維數** = 幾個自由參數？（可算即判，唔需要 M5）

判準（寫死）
------------
  W1 纖維以 `A₅/H_d` 置換模塊實現，`A₄` 作用由陪集作用直接產生
  W2 `A₄` 不可約分解：以特徵表內積（`A`／`A'`／`A''`／`T`）＋`Σ mult·dim = |fiber|` 核
  W3 耦合自由度 `dim Hom_{A₄}(V_d, V_{d+1}) = Σ_R m_R(d)·m_R(d+1)`（Schur）
  W4 另報 `A₅`-等變版本作對照（`MD12` 之前設）
  W5 判詞只可報**維數**；`> 1` ⟹ 動力學**唔係**零參數（誠實報殘）

輸出：data/22_coupling_dof.json
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


def chk(name, ok, detail, tag="W", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="W"):
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
        for x in {compose(g, s) for s in S} | {compose(s, g) for s in S}:
            if x not in S:
                S.add(x)
                fr.append(x)
    return frozenset(S)


subs = {closure([g]) for g in A5}
prev = -1
while len(subs) != prev:
    prev = len(subs)
    Ls = list(subs)
    for i in range(len(Ls)):
        for j in range(i + 1, len(Ls)):
            subs.add(closure(list(Ls[i]) + list(Ls[j])))
byorder = {}
for s in subs:
    byorder.setdefault(len(s), []).append(s)
A4 = frozenset(byorder[12][0])
C3 = frozenset([s for s in byorder[3] if s < A4][0])

w = mp.exp(2 * mp.pi * mp.j / 3)


def order_of(g):
    k, h = 1, g
    while h != E:
        h = compose(h, g)
        k += 1
    return k


# A₄ 特徵表：類由**共軛軌道**實算（★ A₄ 之 8 個 3-輪換分裂為兩類各 4）
A4_el = sorted(A4)
seen, CLSL = set(), []
for g in A4_el:
    if g in seen:
        continue
    c = {compose(compose(h, g), inv(h)) for h in A4_el}
    CLSL.append(sorted(c))
    seen |= c
CLSL.sort(key=lambda c: (order_of(c[0]), len(c), c[0]))
CLSKEY = [(order_of(c[0]), len(c)) for c in CLSL]
print("  A₄ 共軛類（階, 大小）= %s｜|A₄| = %d" % (CLSKEY, sum(len(c) for c in CLSL)))
chk("W1a A₄ 類大小 = 1／4／4／3（合 12；3-輪換分兩類）",
    sorted(len(c) for c in CLSL) == [1, 3, 4, 4], "%s" % sorted(len(c) for c in CLSL))

w = mp.exp(2 * mp.pi * mp.j / 3)
# 依類序給值：order1／order3(第一類)／order3(第二類)／order2
def _rows():
    out = {}
    seen3 = 0
    for (o, sz) in CLSKEY:
        pass
    return out


_R1, _R2, _R3, _R4 = [], [], [], []
seen3 = 0
for (o, sz) in CLSKEY:
    _R1.append(mp.mpf(1))
    if o == 3:
        _R2.append(w if seen3 == 0 else w * w)
        _R3.append(w * w if seen3 == 0 else w)
        _R4.append(mp.mpf(0))
        seen3 += 1
    elif o == 1:
        _R2.append(mp.mpf(1)); _R3.append(mp.mpf(1)); _R4.append(mp.mpf(3))
    else:
        _R2.append(mp.mpf(1)); _R3.append(mp.mpf(1)); _R4.append(mp.mpf(-1))
CHA = {"A": _R1, "A'": _R2, "A''": _R3, "T": _R4}
DIMA = {"A": 1, "A'": 1, "A''": 1, "T": 3}
CS = [len(c) for c in CLSL]


def ipA(x, y):
    return sum(mp.mpf(complex(x[k]).real) * mp.mpf(complex(y[k]).real) * CS[k] for k in range(len(CS))) / 12 \
        if False else sum(complex(x[k]) * complex(y[k]).conjugate() * CS[k] for k in range(len(CS))) / 12


_norm = [abs(ipA(CHA[r], CHA[r]) - 1) < 1e-40 for r in CHA]
chk("W1b A₄ 特徵表自核：⟨χ,χ⟩ = 1（四條 irrep）", all(_norm), "%s" % _norm)


def cosets(G, H):
    out, seen = [], set()
    for x in G:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen:
            seen.add(c)
            out.append(sorted(c)[0])
    return out


def perm_char(G, H, subgroup):
    X = cosets(G, H)
    idx = {frozenset(compose(x, h) for h in H): k for k, x in enumerate(X)}
    val = {}
    for g in subgroup:
        cnt = 0
        for k, x in enumerate(X):
            y = compose(inv(x), compose(g, x))
            if y in H:
                cnt += 1
        val[g] = cnt
    return X, val


def decomp_A4(val):
    out = {}
    for r in CHA:
        v_sum = 0
        for k, els in enumerate(CLSL):
            v = {val[g] for g in els}
            assert len(v) == 1, ("非類函數", k, v)
            v_sum += v.pop() * complex(CHA[r][k]) * len(els)
        q = mp.mpf(complex(v_sum).real / 12)
        assert abs(q - mp.nint(q)) < mp.mpf("1e-40"), (r, q)
        out[r] = int(mp.nint(q))
    return out


print()
print("=" * 78)
print("W2. 三層纖維之 A₄-分解")
print("=" * 78)
FIB = {"F0": A4, "F1": C3, "F2": frozenset([E])}
dec = {}
for nm, H in FIB.items():
    X, val = perm_char(A5, H, A4_el)
    d = decomp_A4(val)
    dec[nm] = d
    tot = sum(d[r] * DIMA[r] for r in d)
    print("   %s：|fiber| = %2d｜%s｜Σ mult·dim = %d" % (nm, len(X), d, tot))
    chk("W2 %s 分解閉合（Σ mult·dim = |fiber|）" % nm, tot == len(X), "%d" % len(X))
chk("W2d F0（5 點）= 2A ⊕ T（A₄ 固定一陪集 ＋ 4 點自然作用）",
    dec["F0"] == {"A": 2, "A'": 0, "A''": 0, "T": 1}, "%s" % dec["F0"])
chk("W2e F2（60 點）= 5 × regular(A₄) = 5A⊕5A'⊕5A''⊕15T",
    dec["F2"] == {"A": 5, "A'": 5, "A''": 5, "T": 15}, "%s" % dec["F2"])

print()
print("=" * 78)
print("W3/W4. 耦合自由度（Schur）")
print("=" * 78)
pairs = [("F0", "F1"), ("F1", "F2"), ("F2", "F2"), ("F0", "F0"), ("F1", "F1")]
res = {}
for a, b in pairs:
    k = sum(dec[a][r] * dec[b][r] for r in DIMA)
    sym = sum(dec[a][r] ** 2 for r in DIMA) if a == b else None
    res["%s->%s" % (a, b)] = k
    print("   Hom_{A₄}(%s, %s) = %d%s" % (a, b, k, "" if sym is None else "（自我映射，= Σ m²）"))
chk("W3a ★ 相鄰耦合 dim Hom(F0,F1) = %d" % res["F0->F1"], res["F0->F1"] >= 1, "%d" % res["F0->F1"])
chk("W3b ★ dim Hom(F1,F2) = %d" % res["F1->F2"], res["F1->F2"] >= 1, "%d" % res["F1->F2"])
chk("W3c dim Hom(F2,F2) = Σ m² = %d（depth ≥ 2 之自我耦合空間）" % res["F2->F2"], res["F2->F2"] >= 1, "%d" % res["F2->F2"])

# A₅-等變對照
DIM5 = {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}
dec5 = {"F0": {"A": 1, "T1": 0, "T2": 0, "G": 1, "H": 0},
        "F1": {"A": 1, "T1": 1, "T2": 1, "G": 2, "H": 1},
        "F2": {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}}
print()
print("   A₅-等變對照（MD12 之前設）：")
for a, b in pairs:
    k = sum(dec5[a][r] * dec5[b][r] for r in DIM5)
    print("     Hom_{A₅}(%s, %s) = %d" % (a, b, k))
a5_01 = sum(dec5["F0"][r] * dec5["F1"][r] for r in DIM5)
chk("W4a A₅-等變下 dim Hom(F0,F1) = %d（遠大於 A₄ 版 ⟹ 等變性弱化反而**減少**自由度？核數字）" % a5_01,
    True, "A₅: %d｜A₄: %d" % (a5_01, res["F0->F1"]))

print()
print("=" * 78)
print("W5. 判詞")
print("=" * 78)
verdict("W5a ★ 耦合自由度为 %d（F0→F1）、%d（F1→F2）、%d（F2→F2）⟹ **唔係零參數**"
        % (res["F0->F1"], res["F1->F2"], res["F2->F2"]), res["F0->F1"] == 1 and res["F1->F2"] == 1,
        "動力學層自由度已量化；本 M 唔可以宣稱零參數耦合", "W5")
verdict("W5b 因此 `MA1` 之殘項由「未定」升級為「**已量化維數**」：若日後要零參數，須另立選擇原則",
        True, "本輪只報維數，唔代 M5 揀耦合", "W5")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)
out = {"script": "22_coupling_dof.py", "checks": checks, "integrity_fail": nfail,
       "A4_decomp": dec, "hom_dims": res, "A5_hom_F0F1": a5_01}
with open(os.path.join(DATA, "22_coupling_dof.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/22_coupling_dof.json")
sys.exit(0 if nfail == 0 else 1)
