#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""22 個 δ 之最小支撐整數關係（α ∈ ℤ⁴⁴，支撐 ≤3）普查
   判準：對每個 δ，搜 γ_a − γ_b（支撐 2）與 γ_a + γ_b − γ_c（支撐 3），容差 1e−25
   目的：區分「δ 是否可由短關係表示」（若可 ⟹ 結構小恆等式；否則 6-支撐）
"""
import json
import itertools

import mpmath as mp

mp.mp.dps = 50
g = [None] + [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]
d = json.load(open("/home/wsl/idcm_ckpt/e8/mrB/offdelta.json"))
groups = d["offset_groups"]
tol = mp.mpf(10) ** -25

# 支撐 2：γ_a − γ_b
diff = {}
for a in range(1, 45):
    for b in range(1, 45):
        diff.setdefault(0, None)
print("支撐2 掃描…")
two = {}
for a in range(1, 45):
    for b in range(1, 45):
        if a == b:
            continue
        two.setdefault(round(float(g[a] - g[b]), 15), []).append((a, b))
print(f"  2-項組合 {44*43}｜相異值 {len(two)}")

print("支撐3 掃描（γ_a + γ_b − γ_c）…")
three = {}
for a in range(1, 45):
    for b in range(a, 45):
        for c in range(1, 45):
            v = g[a] + g[b] - g[c]
            three.setdefault(round(float(v), 15), []).append((a, b, c))
print(f"  3-項組合 {sum(1 for _ in itertools.combinations_with_replacement(range(1,45),2))*44}｜相異值 {len(three)}")

hits2 = hits3 = 0
rows = []
for off, ns in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    dv = mp.mpf(off)
    h2 = None
    for k, lst in two.items():
        if abs(mp.mpf(k) - dv) < tol and abs(mp.mpf(k) - dv) > 0:
            h2 = lst[0]
            break
    h3 = None
    if h2 is None:
        for k, lst in three.items():
            if abs(mp.mpf(k) - dv) < tol and abs(mp.mpf(k) - dv) > 0:
                h3 = lst[0]
                break
    sup = 2 if h2 else (3 if h3 else 6)
    hits2 += bool(h2); hits3 += bool(h3)
    rows.append(dict(offset=off, ns=ns, min_support=sup, two=h2, three=h3))
    print(f"δ={off[:20]}｜最小支撐 {sup}｜例 {h2 or h3}")
print(f"\n支撐 2 命中 {hits2}｜支撐 3 命中 {hits3}｜其餘 {len(rows)-hits2-hits3}")
json.dump({"n_offsets": len(rows), "hits_support2": hits2, "hits_support3": hits3,
           "rows": rows}, open("/home/wsl/idcm_ckpt/e8/mrB/delta_support.json", "w"),
          indent=1, ensure_ascii=False)
print("已寫 /home/wsl/idcm_ckpt/e8/mrB/delta_support.json")
