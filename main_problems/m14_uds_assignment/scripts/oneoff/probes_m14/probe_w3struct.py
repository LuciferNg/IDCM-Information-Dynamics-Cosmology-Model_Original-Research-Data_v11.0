#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""w=+3 之 88 個簽名之結構（承接 §二十三 之疑）"""
import itertools
from collections import Counter

N, MOD = 44, 90
EPS = list(itertools.product([1, -1], repeat=4))
SGN = [e[0] * e[1] * e[2] * e[3] for e in EPS]
ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))

hits = []
for rem in range(1, N + 1):
    for a1, a2, a3 in ADD3:
        v = (rem, a1, a2, a3)
        acc = 0
        mask = []
        for e, s in zip(EPS, SGN):
            ok = (e[0] * v[0] + e[1] * v[1] + e[2] * v[2] + e[3] * v[3]) % MOD == 0
            mask.append(1 if ok else 0)
            if ok:
                acc += s
        if acc // 2 == 3:
            hits.append((rem, a1, a2, a3, tuple(mask)))
print("w=+3 之簽名數 =", len(hits))
print("前 20：")
for h in hits[:20]:
    print("  ", h[:4])
print("\nε 順序（{±1}⁴ 之 itertools 序）：", EPS)
# 檢查結構：是否 rem + a1 + a2 + a3 ≡ c (mod 90)
sums = Counter((r + a1 + a2 + a3) % MOD for r, a1, a2, a3, _ in hits)
print("\nΣ ≡ mod90 之分佈：", dict(sums))
# 檢查是否全部 rem ≡ a1+a2+a3 (mod 90) 或類似
c1 = sum(1 for r, a1, a2, a3, _ in hits if (a1 + a2 + a3 - r) % MOD == 0)
c2 = sum(1 for r, a1, a2, a3, _ in hits if (a1 + a2 + a3 + r) % MOD == 0)
c3 = sum(1 for r, a1, a2, a3, _ in hits if (a1 + a2 - a3 - r) % MOD == 0)
print("rem ≡ a1+a2+a3 :", c1, "｜rem ≡ −(a1+a2+a3):", c2, "｜rem ≡ a1+a2−a3 :", c3)
# 統計三 add 之和
ss = Counter(a1 + a2 + a3 for _, a1, a2, a3, _ in hits)
print("\nadds 和之範圍：", min(ss), "-", max(ss))
print("adds 和直方：", sorted(ss.items())[:20])
