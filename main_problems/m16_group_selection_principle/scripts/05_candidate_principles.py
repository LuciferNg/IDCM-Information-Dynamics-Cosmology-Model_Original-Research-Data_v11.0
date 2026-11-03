#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/05_candidate_principles.py — M16 R1：候選「群選擇原理」之逐一檢定
輸入：data/run_02_scan_all_500.txt（92,570 群）
每個原理皆「事前可述」（不引用 A₅ 之名），報告其命中集與基數（自由度）。
"""
import os, sys, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
SCAN = os.path.join(DATA, "run_02_scan_all_500.txt")

NPASS = 0
NFAIL = 0
OUT = []


def chk(name, cond, detail=""):
    global NPASS, NFAIL
    if cond:
        NPASS += 1
        OUT.append(f"PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        NFAIL += 1
        OUT.append(f"FAIL {name}" + (f"  [{detail}]" if detail else ""))


LINE = re.compile(r"^(\d+)\s+(\d+)\s+\[(.*)\]\s*$")
rows = []
with open(SCAN) as f:
    for ln in f:
        m = LINE.match(ln.strip())
        if not m:
            continue
        degs = []
        for d, k in re.findall(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", m.group(3)):
            degs.extend([int(d)] * int(k))
        rows.append((int(m.group(1)), int(m.group(2)), degs))
N = len(rows)
OUT.append(f"§0 群數 = {N}")
chk("0:count", N == 92570, f"{N}")

n_arr = np.array([r[0] for r in rows])
k_lin = np.zeros(N, dtype=np.int64)
prodDNT = np.ones(N, dtype=np.int64)
sumDNT = np.zeros(N, dtype=np.int64)
sqDNT = np.zeros(N, dtype=np.int64)
cntDNT = np.zeros(N, dtype=np.int64)
consec = np.zeros(N, dtype=bool)
sqALL = np.zeros(N, dtype=np.int64)
IS_A5 = np.zeros(N, dtype=bool)
CAP = 10 ** 7


def sat(vals):
    p = 1
    for v in vals:
        p *= v
        if p > CAP:
            return CAP
    return p


for j, (n, idx, degs) in enumerate(rows):
    if (n, idx) == (60, 5):
        IS_A5[j] = True
    k_lin[j] = sum(1 for d in degs if d == 1)
    sqALL[j] = sum(d * d for d in degs)
    S = sorted({d for d in degs if d > 1})
    prodDNT[j] = sat(S)
    sumDNT[j] = sum(S)
    sqDNT[j] = sum(d * d for d in S)
    cntDNT[j] = len(S)
    consec[j] = len(S) >= 3 and all(S[i + 1] - S[i] == 1 for i in range(len(S) - 1))

A5 = mask = int.from_bytes(np.packbits(IS_A5.astype(np.uint8), bitorder="little").tobytes(), "little")


def M(b):
    return int.from_bytes(np.packbits(np.asarray(b, dtype=bool).astype(np.uint8), bitorder="little").tobytes(), "little")


def rep(name, m):
    tag = "★A5 唯一" if m == A5 else ("⊇A5" if m & A5 else "∌A5")
    OUT.append(f"  {name:52s} |H|={m.bit_count():5d}  {tag}")
    return m


OUT.append("§1 單條件原理（每條皆事前可述）")
P = {}
P["P1 ∏(相異非平凡維數) = |G|"] = rep("P1 ∏(DNT) = |G|", M(prodDNT == n_arr))
P["P2 Σd² = |G|（Burnside）"] = rep("P2 Σd² = |G|（Burnside，全部 irrep）", M(sqALL == n_arr))
P["P3 |G| 之質因數分解可乘性"] = rep("P3 |G/G'| = 1（perfect）", M(k_lin == 1))
P["P4 DNT 為連續整數"] = rep("P4 DNT 連續", M(consec))
P["P5 DNT 連續 ∧ ∏=|G|"] = rep("P5 DNT 連續 ∧ ∏=|G|", M(consec & (prodDNT == n_arr)))
P["P6 ∏=|G| ∧ |G/G'|=1"] = rep("P6 ∏=|G| ∧ perfect", M((prodDNT == n_arr) & (k_lin == 1)))
P["P7 ∏=|G| ∧ ΣDNT=12（本 O 準則）"] = rep("P7 ∏=|G| ∧ ΣDNT=12【原準則】", M((prodDNT == n_arr) & (sumDNT == 12)))
P["P8 ΣDNT = 12"] = rep("P8 ΣDNT = 12", M(sumDNT == 12))
P["P9 ΣDNT = γ_max(B−1) = 12 ∧ ∏DNT=|G| ∧ 三扇區"] = rep("P9 ∏=|G| ∧ ΣDNT=12 ∧ #DNT=3", M((prodDNT == n_arr) & (sumDNT == 12) & (cntDNT == 3)))
P["P10 ∏=|G| ∧ #DNT=3"] = rep("P10 ∏=|G| ∧ #DNT=3", M((prodDNT == n_arr) & (cntDNT == 3)))
# P11 需真 min：
minDNT = np.zeros(N, dtype=np.int64)
for j, (n, idx, degs) in enumerate(rows):
    S = [d for d in set(degs) if d > 1]
    minDNT[j] = min(S) if S else 0
P["P11 ∏=|G| ∧ min DNT = B−1 = 3"] = rep("P11 ∏=|G| ∧ min(DNT)=3", M((prodDNT == n_arr) & (minDNT == 3)))

# P12 / P13：只講度集
IS345 = np.zeros(N, dtype=bool)
SUP345 = np.zeros(N, dtype=bool)
for j, (n, idx, degs) in enumerate(rows):
    S = {d for d in degs if d > 1}
    IS345[j] = (S == {3, 4, 5})
    SUP345[j] = {3, 4, 5} <= S
P["P12 DNT = {3,4,5} 恰好"] = rep("P12 DNT = {3,4,5} 恰好", M(IS345))
P["P13 {3,4,5} ⊆ DNT"] = rep("P13 {3,4,5} ⊆ DNT", M(SUP345))

OUT.append("§2 命中集內容（|H| ≤ 8 者列出；P1/P12 另列全部）")
for name, m in P.items():
    if 0 < m.bit_count() <= 8:
        idxs = [j for j in range(N) if (m >> j) & 1]
        desc = ", ".join(f"{rows[j][0]}/{rows[j][1]}" for j in idxs)
        OUT.append(f"  {name}: {desc}")

OUT.append("§2b P1（∏DNT=|G|，17 群）全列：")
idxs = [j for j in range(N) if (P["P1 ∏(相異非平凡維數) = |G|"] >> j) & 1]
for j in idxs:
    OUT.append(f"  {rows[j][0]:3d}/{rows[j][1]:<4d} DNT={sorted({d for d in rows[j][2] if d > 1})} "
               f"∏={prodDNT[j]} Σ={sumDNT[j]} k_lin={k_lin[j]}")
OUT.append("§2c P12（DNT 恰好 = {3,4,5}，11 群）全列：")
idxs = [j for j in range(N) if (P["P12 DNT = {3,4,5} 恰好"] >> j) & 1]
for j in idxs:
    OUT.append(f"  {rows[j][0]:3d}/{rows[j][1]:<4d} 維數={rows[j][2]} k_lin={k_lin[j]} "
               f"∏={prodDNT[j]} {'= |G|' if prodDNT[j] == rows[j][0] else '≠ |G|'}")

OUT.append("§3 判準（哪些原理對 A₅ 有 1 命中且無自由度）")
uniq = [name for name, m in P.items() if m == A5]
OUT.append(f"  唯一命中 A₅ 之原理：{uniq}")
chk("3:has-unique", len(uniq) >= 1, f"{uniq}")

# S₅ 近敵分析（S₅ = SmallGroup(120,34)，GAP 核實）
s5 = [j for j, r in enumerate(rows) if (r[0], r[1]) == (120, 34)]
OUT.append(f"§4 S₅ 在掃描中之索引：{s5}（SmallGroup(120,34) = S5，GAP 核實）")
if s5:
    j = s5[0]
    OUT.append(f"  S5: |G|={rows[j][0]} k_lin={k_lin[j]} DNT={sorted({d for d in rows[j][2] if d>1})} "
               f"∏={prodDNT[j]} Σ={sumDNT[j]} ∏=|G|? {prodDNT[j]==n_arr[j]} perfect={k_lin[j]==1}")
    chk("4:S5-vs-A5-perfect", k_lin[j] != 1 and k_lin[int(np.where(IS_A5)[0][0])] == 1,
        "S₅ 非 perfect（k=2），A₅ perfect（k=1）⟹ 「perfect」可分辨兩者")

print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)