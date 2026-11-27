#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R1 餘項：16 個實現遮罩之代數結構
   ① w(mask) 是否為單值函數（同一遮罩 ⟹ 同一 w）
   ② 是否對 ∪ / ∩ / 補 封閉（子格？布林子代數？）
   ③ 每遮罩之簽名計數（精確整數）
"""
import itertools
import json
from collections import Counter, defaultdict

N, MOD = 44, 90
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]
FULL = frozenset(range(8))

mk = defaultdict(Counter)
for rem in range(1, N + 1):
    for a in itertools.combinations_with_replacement(range(1, N + 1), 3):
        v = (rem,) + a
        m = frozenset(i for i, sg in enumerate(SIG)
                      if (v[0] * sg[0] + v[1] * sg[1] + v[2] * sg[2] + v[3] * sg[3]) % MOD == 0)
        w = sum(CS[i] for i in m)
        mk[m][w] += 1

m16 = sorted(mk)
print("實現遮罩數 =", len(m16))
for m in m16:
    ws = mk[m]
    print(f"  {sorted(m)}  w={dict(ws)}")

single = all(len(mk[m]) == 1 for m in m16)
wmap = {m: next(iter(mk[m])) for m in m16}
print("① w(mask) 單值:", single)

def closed(op, name):
    S = set(m16)
    bad = []
    for a in m16:
        for b in m16:
            r = frozenset(op(a, b))
            if r not in S:
                bad.append((sorted(a), sorted(b), sorted(r)))
    print(f"② {name} 封閉: {not bad}（違例 {len(bad)}）")
    return bad[:5]

bU = closed(lambda a, b: a | b, "∪")
bI = closed(lambda a, b: a & b, "∩")
comp_bad = [sorted(m) for m in m16 if frozenset(FULL - m) not in set(m16)]
print(f"② 補封閉: {not comp_bad}（違例 {len(comp_bad)}，例 {comp_bad[:3]}）")

# 三點結構：w 是否 = |mask| 之函數
len_w = defaultdict(set)
for m in m16:
    len_w[len(m)].add(wmap[m])
print("③ |mask| → w 集合:", {k: sorted(v) for k, v in sorted(len_w.items())})

out = {"n_masks": len(m16),
       "masks": {",".join(map(str, sorted(m))): {"w": wmap[m], "signatures": sum(mk[m].values())} for m in m16},
       "w_single_valued": single,
       "union_closed": not bU, "inter_closed": not bI, "complement_closed": not comp_bad,
       "union_violations_sample": bU, "inter_violations_sample": bI,
       "cardinality_to_w": {str(k): sorted(v) for k, v in sorted(len_w.items())}}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/r1_mask_algebra.json", "w"), indent=1, ensure_ascii=False)
print("→ 落 mrB/r1_mask_algebra.json")
