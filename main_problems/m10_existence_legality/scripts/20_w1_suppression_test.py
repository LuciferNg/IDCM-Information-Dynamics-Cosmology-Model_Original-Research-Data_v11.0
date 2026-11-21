#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/20 — R21-B：W1 自能通道所需之 3.62× 抑制，有無**結構式**？

背景
  M5 自能通道（∝q²）之量級 = 4.5804 MeV（外部值），W1 需求 = 1.2667 MeV
  ⟹ 若要此通道成立，必須有一個 **0.27654** 之結構性抑制因子（零自由參數）。
  但「搵到一條式就當命中」係 φ-power 陷阱 ⟹ 必須做 **命中密度檢定**：
  問「目標值之命中數是否異常」，而非「有無命中」。

事前登記（寫死，先於計算；不得事後擴大）
  目標     t = 1.2667 / 4.5804
  式族 F（單一比，不含乘積）：
    F1 γ_m/γ_n（1≤m<n≤44，雙向）        γ_n = 2 − 2cos(2πn/88)
    F2 sin(πa/88)/sin(πb/88)（1≤a<b≤43，雙向）
    F3 k/88（k = 1…87）
    F4 1/k（k = 2…880）
    F5 φ^{±k}（k = 1…12）
    F6 L_k/L_j、F_k/F_j（k<j≤12）
  窗    ±1%（相對）與 ±0.1%
  零假設 500 個種子隨機目標（uniform 0.01–1，seed 20261031），同一式族、同一窗
  判準（事前）
    命中數 ≤ 零假設之中位數 ⟹ 無結構信號（巧合級別）
    命中數 ≥ 零假設 95 百分位 ⟹ 異常，值得查（但仍須獨立結構理由，否則不採）

執行：/usr/bin/python3 scripts/20_w1_suppression_test.py
"""
import json
import math
import os
import random
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


N = 88
ND = 44
PHI = (1 + 5 ** 0.5) / 2


def gamma(n):
    return 2.0 - 2.0 * math.cos(2.0 * math.pi * n / N)


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def lucas(k):
    a, b = 2, 1
    for _ in range(k):
        a, b = b, a + b
    return a


# ---------------------------------------------------------------- 式族 F
def build_family():
    fam = []
    # F1 γ 比
    g = {n: gamma(n) for n in range(1, ND + 1)}
    for m in range(1, ND + 1):
        for n in range(m + 1, ND + 1):
            fam.append(("γ%d/γ%d" % (m, n), g[m] / g[n]))
            fam.append(("γ%d/γ%d" % (n, m), g[n] / g[m]))
    # F2 sin 比
    s = {a: math.sin(math.pi * a / N) for a in range(1, 44)}
    for a in range(1, 44):
        for b in range(a + 1, 44):
            fam.append(("sin%d/sin%d" % (a, b), s[a] / s[b]))
            fam.append(("sin%d/sin%d" % (b, a), s[b] / s[a]))
    # F3 k/88
    for k in range(1, 88):
        fam.append(("%d/88" % k, k / N))
    # F4 1/k
    for k in range(2, 881):
        fam.append(("1/%d" % k, 1.0 / k))
    # F5 φ^±k
    for k in range(1, 13):
        fam.append(("φ^%d" % k, PHI ** k))
        fam.append(("φ^-%d" % k, PHI ** -k))
    # F6 Lucas / Fibonacci 比
    for k in range(1, 13):
        for j in range(k + 1, 13):
            fam.append(("L%d/L%d" % (k, j), lucas(k) / lucas(j)))
            fam.append(("L%d/L%d" % (j, k), lucas(j) / lucas(k)))
            fam.append(("F%d/F%d" % (k, j), fib(k) / fib(j)))
            fam.append(("F%d/F%d" % (j, k), fib(j) / fib(k)))
    # 去重（按值）
    seen, out = set(), []
    for nm, v in fam:
        key = round(v, 12)
        if key in seen or not (0.0 < v < 10.0):
            continue
        seen.add(key)
        out.append((nm, v))
    return out


FAM = build_family()
print("=" * 78)
print("M10 R21-B — W1 自能通道之抑制因子結構檢定")
print("=" * 78)
print(f"  式族 |F| = {len(FAM)}（單一比；不含乘積）")

TARGET = 1.2667 / 4.5804
print(f"  目標 t = 1.2667 / 4.5804 = {TARGET:.6f}")


def hits(x, rel):
    lo, hi = x * (1 - rel), x * (1 + rel)
    return [(nm, v) for nm, v in FAM if lo <= v <= hi]


for rel, lbl in ((0.01, "±1%"), (0.001, "±0.1%")):
    h = hits(TARGET, rel)
    print(f"\n[目標命中 {lbl}] {len(h)} 條")
    for nm, v in sorted(h, key=lambda z: abs(z[1] - TARGET))[:8]:
        print(f"      {nm:<14} = {v:.6f}   偏 {100*(v/TARGET-1):+.3f}%")
    ck(f"S.1 {lbl} 命中數已計", True, f"{len(h)} 條")

# ---------------------------------------------------------------- 零假設
print("\n[零假設] 500 個隨機目標（seed 20261031），同一式族同一窗")
random.seed(20261031)
for rel, lbl in ((0.01, "±1%"), (0.001, "±0.1%")):
    cnts = []
    for _ in range(500):
        x = random.uniform(0.01, 1.0)
        cnts.append(len(hits(x, rel)))
    cnts.sort()
    mean = sum(cnts) / len(cnts)
    med = cnts[len(cnts) // 2]
    p95 = cnts[int(0.95 * len(cnts))]
    actual = len(hits(TARGET, rel))
    pct = 100.0 * sum(1 for c in cnts if c <= actual) / len(cnts)
    print(f"      {lbl}: 實際 {actual} ｜ 隨機 mean {mean:.2f} / median {med} / p95 {p95} / max {cnts[-1]}"
          f" ｜ 實際所處百分位 {pct:.0f}%")
    ck(f"S.2 {lbl} 命中數 ≤ 零假設中位數（無結構信號）", actual <= med,
       f"實際 {actual} vs median {med} ⟹ " +
       ("巧合級別，**無結構信號**" if actual <= med else "異常，須獨立結構理由方可採"))

# ---------------------------------------------------------------- 鄰近對照
print("\n[對照] 已知結構常數落點（不參與判準，只顯示式族之覆蓋密度）")
for nm, v in (("4/15（A₅ 阻挫邊）", 4 / 15), ("γ₂/γ₁", gamma(2) / gamma(1)),
              ("1/4", 0.25), ("1/3", 1 / 3)):
    print(f"      {nm:<18} = {v:.6f}   命中 {len(hits(v, 0.01))} 條（±1%）")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "20_w1_suppression_test.py", "role": "R21-B W1 抑制因子結構檢定",
       "family_size": len(FAM), "target": TARGET,
       "hits_pct1": len(hits(TARGET, 0.01)), "hits_pct01": len(hits(TARGET, 0.001)),
       "hits_pct1_list": [[nm, v] for nm, v in sorted(hits(TARGET, 0.01), key=lambda z: abs(z[1] - TARGET))],
       "null": {}, "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
random.seed(20261031)
for rel, key in ((0.01, "pct1"), (0.001, "pct01")):
    cnts = sorted(len(hits(random.uniform(0.01, 1.0), rel)) for _ in range(500))
    out["null"][key] = {"mean": sum(cnts) / len(cnts), "median": cnts[len(cnts) // 2],
                        "p95": cnts[int(0.95 * len(cnts))], "max": cnts[-1],
                        "actual": len(hits(TARGET, rel))}
with open(os.path.join(DATA, "20_w1_suppression.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/20_w1_suppression.json 已寫")