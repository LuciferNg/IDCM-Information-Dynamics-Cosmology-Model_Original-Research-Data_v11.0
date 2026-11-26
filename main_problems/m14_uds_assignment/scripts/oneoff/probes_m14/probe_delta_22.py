#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""δ 枚舉規則之第三代假說：2↔2 等和替換（γ_a + γ_b = γ_c + γ_d）
   步驟：
   ① 建 1..44 之無序對和類（dps 60 精確聚類）
   ② 對每個 δ 組，檢其成員 5-集是否可由「單次 2↔2 等和替換」互相到達（圖連通性）
   ③ 對照：若全部組皆連通 ⟹ 規則 = 2↔2 替換；否則異質再細分
"""
import json
import itertools
from collections import defaultdict

import mpmath as mp

mp.mp.dps = 60
H = 45
g = lambda k: 4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2

# ① 對和類
pairs = list(itertools.combinations_with_replacement(range(1, 45), 2))
sums = [g(a) + g(b) for a, b in pairs]
classes = defaultdict(list)
tol = mp.mpf(10) ** -40
for (a, b), s in zip(pairs, sums):
    key = None
    for k in classes:
        if abs(k - s) < tol:
            key = k
            break
    if key is None:
        key = s
    classes[key].append((a, b))
multi = {k: v for k, v in classes.items() if len(v) > 1}
print(f"對數 {len(pairs)}｜不同對和 {len(classes)}｜多代表類（等和） {len(multi)}")
sizes = sorted((len(v) for v in multi.values()), reverse=True)
print("等和類之大小分佈（前 10）:", sizes[:10])
print("最大等和類之成員（前 12）:", list(multi[max(multi, key=lambda k: len(multi[k]))])[:12])

# ② 各組之 2↔2 連通性
d = json.load(open("/home/wsl/idcm_ckpt/e8/mrB/offdelta.json"))
rows = {r["n"]: r for r in d["rows"]}
eq_pairs = {}
for k, v in multi.items():
    for p in v:
        eq_pairs[p] = v


def neighbours(S):
    """S: tuple 5 元素（多重集）→ 所有經單次 2↔2 等和替換（含 γ_n 之位置不變）可得之 5-集集合"""
    out = set()
    idx = range(len(S))
    for i, j in itertools.combinations(idx, 2):
        a, b = S[i], S[j]
        for (c, dd) in eq_pairs.get((min(a, b), max(a, b)), [(a, b)]):
            if (c, dd) == (a, b):
                continue
            T = list(S)
            T[i], T[j] = c, dd
            out.add(tuple(sorted(T)))
    return out


print()
for off, ns in sorted(d["offset_groups"].items(), key=lambda kv: -len(kv[1])):
    sets = {n: tuple(rows[n]["set"]) for n in ns}
    # 只檢是否可由某成員經 2↔2 到達其他成員
    reached = set()
    stack = [sets[ns[0]]]
    seen = {sets[ns[0]]}
    while stack:
        cur = stack.pop()
        for nb in neighbours(cur):
            if nb in seen:
                continue
            seen.add(nb)
            stack.append(nb)
    hit = sum(1 for n in ns if sets[n] in seen)
    print(f"δ={off[:20]}｜成員 {len(ns)}｜2↔2 可達 {hit}/{len(ns)}｜連通 {'是' if hit==len(ns) else '否'}")
    if hit < len(ns):
        miss = [n for n in ns if sets[n] not in seen]
        print("   未達成員：", miss, [sets[n] for n in miss][:3])
