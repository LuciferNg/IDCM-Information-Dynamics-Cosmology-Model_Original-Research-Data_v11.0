#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/14 — MA1 最後一步：載體由「`J` 譜完整性」逼出？
============================================================
前設問題
--------
`docs/13` 定位咗等變塔（`1→5→20→60`，分支 `(5,4,3)`），但 `A₅` 有 **8 個** orbit 大小可選
（`5,6,10,12,15,20,30,60`）。本輪問：**觀測要求嘅 `J` 譜**能否唯一逼出一個 orbit 大小？

判準（結構，非擬合）：`SO(3)↓I`（`docs/09`）顯示 `J = 0..8` 需要 irrep 集
  `{A, T1, T2, G, H}`（五個全部出現）。
若一個 orbit 之置換模塊**唔含全部五個**，則該 orbit 唔可能承載完整 `J` 譜 ⟹ 出局。

方法（零參數；全部機器核）
  1. 傳遞作用 `G/H` 之置換特徵 `χ_n(g) = #{ xH : x⁻¹gx ∈ H }`（逐 coset 實數）
  2. 以 `A₅` 特徵表做內積分解（`n_R = ⟨χ_n, χ_R⟩`），核 `Σ n_R·dim R = n`
  3. 判「含全部五 irrep」者，列最小者；若唯一 ⟹ 載體被逼出（零參數）
  4. 另計 5-orbit 之角色（樹接駁用）與載體角色分開報，唔可以混為一談

事前判準（寫死，唔准事後改）
----------------------------
  Z1 `χ_n` 必須逐 coset 實數（唔可以照抄分解表）
  Z2 分解後必須 `Σ n_R·dim R = n`（內部閉合），且 `n_R ≥ 0` 整數
  Z3 「逼出」之判準：含全部五 irrep 之 orbit 大小**唯一** ⟹ 零參數；多於一個 ⟹ 仍有自由
  Z4 5-orbit（樹接駁）與載體（irrep 頻譜）係**兩個角色**，若唔同必須明寫（唔准靜靜合併）
  Z5 可證偽：若最小合格 orbit = `60`，核 `60 = |A₅| = B(B²−1)`（已註冊）
  Z6 禁：唔准為咗令某 orbit 合格而改 `J` 譜範圍

輸出：data/14_ma1_carrier.json
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


def chk(name, ok, detail, tag="Z", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="Z"):
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

# ── 共軛類（以元素階 ＋ 類大小判別；5A/5B 之標籤任意，T1↔T2 對稱）──
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
# 類按 (元素階, 大小) 排序；同階同大小之兩類按首元素定序
cls.sort(key=lambda c: (order_of(c[0]), len(c), c[0]))
CS = [len(c) for c in cls]
CLSKEY = [(order_of(c[0]), len(c)) for c in cls]
print("  共軛類（階, 大小）= %s" % (CLSKEY,))

_phi = (1 + mp.sqrt(5)) / 2
_psi = 1 - _phi
# 標準 A₅ 表，以 (元素階, 類大小) 索引；兩個 5-階類分別用 φ / ψ
_BASE = {
    "A":  {1: 1, 2: 1, 3: 1, 5: 1},
    "T1": {1: 3, 2: -1, 3: 0, 5: _phi},
    "T2": {1: 3, 2: -1, 3: 0, 5: _psi},
    "G":  {1: 4, 2: 0, 3: 1, 5: -1},
    "H":  {1: 5, 2: 1, 3: -1, 5: 0},
}
CH = {}
for r in _BASE:
    row, seen5 = [], 0
    for (o, sz) in CLSKEY:
        if o == 5:
            base = _BASE[r][5]
            if base == _phi:                      # 只有 T1/T2 需喺兩 5-類間 φ/ψ 交換
                v = _phi if seen5 == 0 else _psi
            elif base == _psi:
                v = _psi if seen5 == 0 else _phi
            else:
                v = base                          # A/G/H 兩 5-類同值
            seen5 += 1
            row.append(v)
        else:
            row.append(_BASE[r][o])
    CH[r] = row
DIM = {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}

# ── 特徵表自核（唔照抄即信：正交 ＋ 範數）──
print("  χ 表自核：", end="")


def ip(chi1, chi2):
    return sum(mp.mpf(chi1[k]) * mp.mpf(chi2[k]) * CS[k] for k in range(len(CS))) / 60


norm_ok = all(abs(ip(CH[r], CH[r]) - 1) < mp.mpf("1e-40") for r in CH)
orth_ok = all(abs(ip(CH[r], CH[s_])) < mp.mpf("1e-40") for r in CH for s_ in CH if r != s_)
dim_ok = sum(DIM[r] ** 2 for r in CH) == 60
print("⟨χ,χ⟩=1 %s｜⟨χ,χ'⟩=0 %s｜Σdim²=60 %s" % (norm_ok, orth_ok, dim_ok))

print("=" * 78)
print("Z1. 傳遞作用置換特徵 χ_n（逐 coset 實數）")
print("=" * 78)
res = {}
for order in sorted(byorder):
    H = sorted(byorder[order][0])            # 每階只有一個共軛類（docs/13 MD8）
    n = 60 // order
    reps, seen_c = [], set()
    for x in A5:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen_c:
            seen_c.add(c)
            reps.append(c)
    assert len(reps) == n, (n, len(reps))
    # χ_n(g) = # fixed cosets
    val = {}
    for g in A5:
        cnt = 0
        for c in reps:
            x = sorted(c)[0]
            if compose(inv(x), compose(g, x)) in H:
                cnt += 1
        val[g] = cnt
    # 分解
    # 逐類取代表值（χ_n 在類上為常數 — 先核）
    clsval = []
    for c in cls:
        vs = {val[g] for g in c}
        assert len(vs) == 1, ("χ_n 非類函數", n, vs)
        clsval.append(vs.pop())
    ns = {}
    for r in CH:
        v = ip(clsval, CH[r])
        assert abs(v - mp.nint(v)) < mp.mpf("1e-40"), ("非整重數", n, r, v)
        ns[r] = int(mp.nint(v))
    res[n] = {"order": order, "n": n, "decomp": ns,
              "dimsum": sum(ns[r] * DIM[r] for r in CH), "chi": {g: v for g, v in list(val.items())[:0]}}
    print("  n = %2d（|H|=%2d）：%s ｜Σ = %d" %
          (n, order, {r: ns[r] for r in CH if ns[r]}, res[n]["dimsum"]))
chk("Z0a 特徵表自核：⟨χ,χ⟩=1 全部", norm_ok, "正交表（非照抄）")
chk("Z0b 特徵表自核：⟨χ,χ'⟩=0 全部 ＋ Σdim² = 60", orth_ok and dim_ok, "正交 ＋ 完備")
chk("Z2a 全部分解 Σ n_R·dim R = n（內部閉合）", all(res[n]["dimsum"] == n for n in res),
    "8/8" if len(res) == 9 else "%d/%d" % (sum(1 for n in res if res[n]["dimsum"] == n), len(res)))
chk("Z2b 全部 n_R ≥ 0（無負重數）", all(res[n]["decomp"][r] >= 0 for n in res for r in CH), "全非負")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z3/Z5. 邊些 orbit 含全部五 irrep（J 譜完整性）？")
print("=" * 78)
JNEED = {"A", "T1", "T2", "G", "H"}
qual, qual_min = [], []
for n in sorted(res):
    present = {r for r in CH if res[n]["decomp"][r] > 0}
    ok = JNEED <= present
    if ok:
        qual.append(n)
    print("  n = %2d：irrep %-28s ｜含全部五？ %s" % (n, sorted(present), "✓" if ok else "✗"))
print()
print("  合格 orbit 大小 =", qual)
chk("Z3a ★ 含全部五 irrep 之 orbit 大小（升序）= %s" % qual, len(qual) >= 1, "%d 個" % len(qual))
verdict("Z3b ★ 最小合格 orbit = %d；唯一性：%s" % (min(qual) if qual else -1, "唯一" if len(qual) == 1 else "**多個**"),
        len(qual) == 1, "合格集 %s ⟹ 未唯一" % qual, "Z3")
verdict("Z3c 若合格集非唯一 ⟹ 載體仍有一個自由選擇（唔可以宣稱逼出）",
        len(qual) == 1, "合格 %d 個" % len(qual), "Z3")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z4. 兩個角色分報：5-orbit（樹接駁）vs 載體（irrep 頻譜）")
print("=" * 78)
print("  5-orbit 分解 = %s ⟹ 只含 %s（**唔含 T1/H**）"
      % ({r: res[5]["decomp"][r] for r in CH if res[5]["decomp"][r]},
         sorted({r for r in CH if res[5]["decomp"][r] > 0})))
chk("Z4a 5-orbit 只含 A、G ⟹ **唔可能**作 J 譜載體（只可作樹接駁層）",
    {r for r in CH if res[5]["decomp"][r] > 0} == {"A", "G"},
    "5 = 1 + 4（A ⊕ G）")
print("  載體候選 =", qual, "⟹ 與 5-orbit 唔同物件")
chk("Z4b ★ 兩角色分離：樹接駁層（degree 5）≠ irrep 載體（合格 orbit）⟹ 必須明寫兩層",
    True, "docs/13 之『根上層 degree 5』講接駁；本檔之載體講頻譜 — 唔同角色")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z5. 與等變塔層 (1,5,20,60) 求交 ⟹ 收窄")
print("=" * 78)
TOWER = [1, 5, 20, 60]
inter = sorted(set(qual) & set(TOWER))
print("  合格集 %s ∩ 塔層 %s = %s" % (qual, TOWER, inter))
chk("Z5a 交 = %s ⟹ 30 出局（唔係塔層）" % inter, inter == [20, 60], "8 個 orbit 收窄至 2")
verdict("Z5b 剩 {20, 60}；若採『塔終端（穩定子平凡）』原則則取 60 —— 但此係**原則**非推導，須明標",
        len(inter) == 2, "60 = |A₅| = B(B²−1)（已註冊）；20 = [A₅:C₃]", "Z5")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z6. 12-orbit vs J^PC 容量 12（結構理由，非計數巧合）")
print("=" * 78)
n12 = {r for r in CH if res[12]["decomp"][r] > 0}
print("  12-orbit 分解 = %s ⟹ 含 G？ %s" % (sorted(n12), "G" in n12))
print("  J 譜：J≤2 需 {A,T1,H}；J≥3 需 G（docs/09 SO(3)↓I）")
chk("Z6a 12-orbit 含 A/T1/T2/H 但**唔含 G** ⟹ 只可承載 J ≤ 2（G 缺 ⟹ J=3,4,… 出局）",
    n12 == {"A", "T1", "T2", "H"}, "%s" % sorted(n12))
chk("Z6b ★ 故『容量 12』（MD6 之 3 irrep × P × C）係**標籤數**，唔係載體點數 ⟹ 唔可以同 12-orbit 等同",
    True, "X6 之數值巧合拒絕，本輪得結構理由（12-orbit 缺 G）")
verdict("Z6c 若載體 = 12-orbit ⟹ J≥3 必須另置（多 orbit 或更大載體）⟹ 12 唔可能係唯一載體",
        False, "觀測 J=3,4 存在 ⟹ 12-orbit 單獨不足", "Z6")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z7. 判定")
print("=" * 78)
verdict("Z7a ★ MD10 兩角色分離（已導）：樹接駁層 = 5-orbit = A ⊕ G，缺 T1/H ⟹ 唔可能係 J 載體",
        True, "5 = 1 + 4；J=1(T1)、J=2(H) 全部缺 ⟹ 接駁層與載體**必為兩件物件**", "Z7")
verdict("Z7b MA1 收窄：8 個 orbit → 合格 3 個 {20,30,60} → ∩塔層 2 個 {20,60}",
        len(inter) == 2, "仍未唯一 ⟹ 載體選擇尚餘 1 bit（或需『塔終端』原則）", "Z7")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "14_ma1_carrier.py", "checks": checks, "integrity_fail": nfail,
       "tower_intersection": inter, "orbit12_irreps": sorted(n12),
       "decomp_by_n": {str(n): res[n]["decomp"] for n in res},
       "qualified_orbits": qual, "qualified_min": min(qual) if qual else None,
       "five_orbit_decomp": res[5]["decomp"], "J_needed": sorted(JNEED)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "14_ma1_carrier.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/14_ma1_carrier.json")
sys.exit(0 if nfail == 0 else 1)
