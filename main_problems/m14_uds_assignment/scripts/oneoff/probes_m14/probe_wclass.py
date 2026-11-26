#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
殘留「17 個 w 值之完整分類」：純整數精確枚舉
簽名 v = (rem, add1, add2, add3)，rem ∈ 1..44，adds ∈ 多重組合(1..44,3)
w(v) = ½ Σ_{ε∈{±1}⁴} (∏ε) [ε·v ≡ 0 mod 90]     （§十二 之恆等式，此處全整數實作）
輸出：每個 w 值之（簽名數、正負、以及 w 值集合是否封閉於某算術結構）
"""
import itertools
import json
from collections import Counter

N = 44
MOD = 90
EPS = list(itertools.product([1, -1], repeat=4))
SGN = [e[0] * e[1] * e[2] * e[3] for e in EPS]

ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))
print(f"簽名空間：44 × {len(ADD3)} = {44 * len(ADD3)}")

cnt = Counter()
pos = Counter()
neg = Counter()
zero = 0
for rem in range(1, N + 1):
    for a1, a2, a3 in ADD3:
        v = (rem, a1, a2, a3)
        acc = 0
        for e, s in zip(EPS, SGN):
            if (e[0] * v[0] + e[1] * v[1] + e[2] * v[2] + e[3] * v[3]) % MOD == 0:
                acc += s
        w = acc // 2
        cnt[w] += 1
        if w > 0:
            pos[w] += 1
        elif w < 0:
            neg[w] += 1
        else:
            zero += 1

vals = sorted(cnt)
print(f"\nw 值集合（{len(vals)} 個）：{vals}")
print(f"零值簽名 {zero}｜總 {44 * len(ADD3)}")
print("\n w ｜簽名數 ｜負 ｜正")
for w in vals:
    print(f"{w:3d}｜{cnt[w]:8d}｜{neg[w]:7d}｜{pos[w]:7d}")
tot = sum(cnt.values())
print(f"\n總計 {tot}")
print("w 之加權和 Σw·count =", sum(w * c for w, c in cnt.items()))
json.dump({"n_values": len(vals), "values": vals,
           "counts": {str(w): cnt[w] for w in vals},
           "neg": {str(w): neg[w] for w in vals},
           "pos": {str(w): pos[w] for w in vals},
           "zero_count": zero, "total": tot,
           "weighted_sum": sum(w * c for w, c in cnt.items())},
          open("/home/wsl/idcm_ckpt/e8/mrB/wclass.json", "w"), indent=1, ensure_ascii=False)
print("已寫 /home/wsl/idcm_ckpt/e8/mrB/wclass.json")
