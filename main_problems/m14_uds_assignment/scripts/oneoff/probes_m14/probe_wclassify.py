#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""w ∈ {−1,+1,+2,3} 之簽名結構分類（承接 §二十四）"""
import itertools
from collections import Counter, defaultdict

N, MOD = 44, 90
EPS = list(itertools.product([1, -1], repeat=4))
SGN = [e[0] * e[1] * e[2] * e[3] for e in EPS]
ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))
H = MOD // 2  # 45

by_val = defaultdict(list)
for rem in range(1, N + 1):
    for a1, a2, a3 in ADD3:
        v = (rem, a1, a2, a3)
        acc = 0
        for e, s in zip(EPS, SGN):
            if (e[0] * v[0] + e[1] * v[1] + e[2] * v[2] + e[3] * v[3]) % MOD == 0:
                acc += s
        by_val[acc // 2].append(v)

for val in (-1, 1, 2, 3):
    vs = by_val[val]
    print(f"\n=== w = {val:+d} | 簽名數 {len(vs)}")
    mult = Counter()
    for v in vs:
        mult[tuple(sorted(Counter(v).values(), reverse=True))] += 1
    print("  重數型分佈:", dict(sorted(mult.items(), key=lambda kv: -kv[1])))
    nrefl = Counter()
    for v in vs:
        cnt = 0
        used = set()
        for i, x in enumerate(v):
            if i in used:
                continue
            for j in range(i + 1, len(v)):
                if j not in used and (v[i] + v[j]) % MOD == 0:
                    cnt += 1
                    used.add(i)
                    used.add(j)
                    break
        nrefl[cnt] += 1
    print("  模90反射對數分佈（x+y≡0）:", dict(nrefl))
    sums = Counter(sum(v) % MOD for v in vs)
    print("  Σ mod 90（前 6）:", sorted(sums.items(), key=lambda kv: -kv[1])[:6])
    print("  前 12 例:", vs[:12])
