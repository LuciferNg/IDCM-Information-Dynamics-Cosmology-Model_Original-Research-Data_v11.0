#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""核 CS 版 W_table 與 ∏ε 版 W_ident 在【簽名域】上是否恆等（§十二 只在 1..44 隨機 4-向量測）
   並取 W_cs 之值域／負份額，與 §十一 之 32.9260%／17 值 [−10,+30] 對照
"""
import itertools
from collections import Counter

N, MOD = 44, 90
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]
E4 = list(itertools.product([1, -1], repeat=4))
PE = [e[0] * e[1] * e[2] * e[3] for e in E4]
ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))


def w_cs(v):
    s = 0
    for c, sg in zip(CS, SIG):
        if sum(a * b for a, b in zip(v, sg)) % MOD == 0:
            s += c
    return s


def w_id(v):
    s = 0
    for e, p in zip(E4, PE):
        if sum(a * b for a, b in zip(v, e)) % MOD == 0:
            s += p
    return s // 2


cnt = Counter(); cnt2 = Counter(); mism = 0; tot = 0
for rem in range(1, N + 1):
    for a1, a2, a3 in ADD3:
        v = (rem, a1, a2, a3)
        a, b = w_cs(v), w_id(v)
        cnt[a] += 1; cnt2[b] += 1; tot += 1
        if a != b:
            mism += 1
print(f"簽名域 {tot}：W_cs 與 W_ident 不符 = {mism}")
print("W_cs 值集合（%d 個）: %s" % (len(cnt), sorted(cnt)))
print("  ", {k: cnt[k] for k in sorted(cnt)})
neg = sum(c for k, c in cnt.items() if k < 0)
print(f"W_cs < 0 之份額 = {neg}/{tot} = {neg/tot:.6f}")
print("W_ident 值集合:", sorted(cnt2), {k: cnt2[k] for k in sorted(cnt2)})
