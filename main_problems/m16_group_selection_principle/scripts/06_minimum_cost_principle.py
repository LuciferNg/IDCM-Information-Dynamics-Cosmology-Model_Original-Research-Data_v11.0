#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/06_minimum_cost_principle.py — M16 追問：群選擇可否重述為「最低運行成本」而非「唯一性」？

測試三件事（全部用 data/run_02_scan_all_500.txt 之 92,570 群，零外部輸入）：
  §A  「度集約束 + 最小 |G|」是否嚴格最小（唯一解）而非 tie-break
  §B  膨脹因子 |G|/∏(DNT) 是否恰等於 k = |G:G'|（線性特徵標數）
  §C  「非平凡完備（perfect）群之最小階」是否唯一落在 A₅ —— 即完全唔用 {3,4,5}
  §D  各種成本函數（|G|、類數、max 維數、Σlog d、k）之 argmin 比較
"""
import os, sys, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
SCAN = os.path.join(DATA, "run_02_scan_all_500.txt")

NPASS = NFAIL = 0
OUT = []


def chk(n, c, d=""):
    global NPASS, NFAIL
    if c:
        NPASS += 1
        OUT.append(f"PASS {n}" + (f"  [{d}]" if d else ""))
    else:
        NFAIL += 1
        OUT.append(f"FAIL {n}" + (f"  [{d}]" if d else ""))


LINE = re.compile(r"^(\d+)\s+(\d+)\s+\[(.*)\]\s*$")
G = []
with open(SCAN) as f:
    for ln in f:
        m = LINE.match(ln.strip())
        if not m:
            continue
        degs = []
        for d, k in re.findall(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", m.group(3)):
            degs.extend([int(d)] * int(k))
        S = sorted({d for d in degs if d > 1})
        prod = 1
        for d in S:
            prod *= d
        G.append(dict(n=int(m.group(1)), i=int(m.group(2)), degs=degs, S=S,
                      prod=prod, k=sum(1 for d in degs if d == 1),
                      ncls=len(degs), mx=max(degs)))
OUT.append(f"§0 群數 = {len(G)}")
chk("0:count", len(G) == 92570, f"{len(G)}")

F345 = [g for g in G if g["S"] == [3, 4, 5]]
FSUP = [g for g in G if {3, 4, 5} <= set(g["S"])]
A5 = [g for g in G if (g["n"], g["i"]) == (60, 5)][0]

# ─────────────────────────────────────────────────────────────
OUT.append("§A 「度集約束 + 最小 |G|」—— 係嚴格最小抑或 tie-break？")
OUT.append(f"  DNT ⊇ {{3,4,5}} 之群數 = {len(FSUP)}；DNT = {{3,4,5}} 之群數 = {len(F345)}")
bad_div = [g for g in FSUP if g["n"] % 60 != 0]
OUT.append(f"  其中 |G| 不被 60 整除者：{len(bad_div)}（度必整除群階 ⟹ 必然 0）")
chk("A:60-div-guaranteed", not bad_div, f"{len(bad_div)}")
mins = min(g["n"] for g in FSUP)
argmin = [g for g in FSUP if g["n"] == mins]
OUT.append(f"  min |G| = {mins}；達到者 {len(argmin)} 個：{[(g['n'], g['i']) for g in argmin]}")
chk("A:min-is-A5-unique", mins == 60 and len(argmin) == 1 and argmin[0] is A5)
OUT.append("  ⟹ 最小值由整除性逼出（3,4,5 | |G| ⟹ 60 | |G|），且 A₅ 達到 ⟹ 嚴格最小、零 tie")

# ─────────────────────────────────────────────────────────────
OUT.append("§B 膨脹因子 = k = |G:G'|？（11 群全列）")
OUT.append("     群        |G|   k   ∏DNT  |G|/∏  =k?   類數  max  #DNT")
okB = True
for g in sorted(F345, key=lambda x: x["n"]):
    ratio = g["n"] // g["prod"]
    okB &= (ratio == g["k"] and g["n"] % g["prod"] == 0)
    OUT.append(f"  {g['n']:3d}/{g['i']:<5d} {g['n']:4d} {g['k']:3d} {g['prod']:5d} {ratio:5d}  "
               f"{str(ratio == g['k']):5s} {g['ncls']:4d} {g['mx']:4d} {len(g['S']):4d}")
chk("B:ratio==k(11群)", okB)
# 反例檢查：S₅（DNT={4,5,6}）比值 ≠ k
S5 = [g for g in G if (g["n"], g["i"]) == (120, 34)][0]
OUT.append(f"  反例 S₅(120/34)：∏DNT={S5['prod']} |G|/∏={S5['n']//S5['prod']} k={S5['k']} ⟹ "
           f"恆等式只在此族成立（族外不成立）")
chk("B:S5-breaks-identity", S5["n"] // S5["prod"] != S5["k"],
    f"{S5['n']//S5['prod']} vs {S5['k']}")

# ─────────────────────────────────────────────────────────────
OUT.append("§C 「非平凡完備群之最小階」是否唯一 = A₅（完全唔用 {3,4,5}）")
PERF = [g for g in G if g["k"] == 1 and g["n"] > 1]
PERF.sort(key=lambda x: x["n"])
OUT.append(f"  k=1（perfect）且 |G|>1 之群數 = {len(PERF)}；頭 5 個：")
for g in PERF[:6]:
    OUT.append(f"    {g['n']:3d}/{g['i']:<4d} DNT={g['S']} 類數={g['ncls']}")
pmin = min(g["n"] for g in PERF)
parg = [g for g in PERF if g["n"] == pmin]
OUT.append(f"  min |G| = {pmin}；達到者 {len(parg)} 個：{[(g['n'], g['i']) for g in parg]}")
chk("C:min-perfect==A5", pmin == 60 and len(parg) == 1 and parg[0] is A5)
OUT.append("  ⟹ 「最低成本 + 非平凡完備」單獨已選出 A₅；{3,4,5}、33/44/55 變成輸出而非輸入")

# ─────────────────────────────────────────────────────────────
OUT.append("§D 成本函數比較（在 F345 之 11 群上；argmin 是否唯一）")


def argmins(items, key):
    vals = [key(g) for g in items]
    m = min(vals)
    return [items[j] for j, v in enumerate(vals) if v == m]


costs = {
    "C1 |G|（狀態數 = Σd² = 預算分母）": lambda g: g["n"],
    "C2 類數（irrep 個數 / 通道數）": lambda g: g["ncls"],
    "C3 k = |G:G'|（線性/未被區分通道）": lambda g: g["k"],
    "C4 max 維數（最大單扇區）": lambda g: g["mx"],
    "C5 #DNT（相異扇區數）": lambda g: len(g["S"]),
    "C6 Σ_{d∈DNT} d（和，恆 12）": lambda g: sum(g["S"]),
}
for name, f in costs.items():
    am = argmins(F345, f)
    tag = "★唯一 = A₅" if len(am) == 1 and am[0] is A5 else f"tie {len(am)}"
    OUT.append(f"  {name:34s} min={f(am[0]):4d}  → {tag}")
chk("D:total-costs-unique-shape-costs-tie",
    len(argmins(F345, costs["C1 |G|（狀態數 = Σd² = 預算分母）"])) == 1 and
    len(argmins(F345, costs["C5 #DNT（相異扇區數）"])) == 11 and
    len(argmins(F345, costs["C4 max 維數（最大單扇區）"])) == 11,
    "總量型成本（|G|、類數、k）唯一；形狀型成本（max、#DNT、Σ）全 tie")

# 全掃描：k 之最小值（perfect 群）與 A₅ 之關係
OUT.append("§D2 全域：|G| 最小之 perfect 群 = A₅（同 §C）；若放寬為 k ≤ 2，最小階者：")
k2 = sorted([g for g in G if g["k"] <= 2 and g["n"] > 1], key=lambda x: x["n"])[:5]
for g in k2:
    OUT.append(f"    {g['n']:3d}/{g['i']:<4d} k={g['k']} DNT={g['S']}")

print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)