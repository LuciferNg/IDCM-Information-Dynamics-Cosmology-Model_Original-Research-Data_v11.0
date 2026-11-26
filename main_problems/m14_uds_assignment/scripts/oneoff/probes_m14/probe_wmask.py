#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""四腿 W 之完整分類：w 值 ⟺ ε-條件遮罩類（8-bit）之精確計數"""
import itertools
import json
from collections import Counter, defaultdict

N, MOD = 44, 90
EPS = list(itertools.product([1, -1], repeat=4))
SGN = [e[0] * e[1] * e[2] * e[3] for e in EPS]
ADD3 = list(itertools.combinations_with_replacement(range(1, N + 1), 3))
H = MOD // 2

mask_of = defaultdict(Counter)   # w -> Counter(mask)
sample = defaultdict(dict)       # w -> mask -> 首例
for rem in range(1, N + 1):
    for a1, a2, a3 in ADD3:
        v = (rem, a1, a2, a3)
        m = 0
        acc = 0
        for i, (e, s) in enumerate(zip(EPS, SGN)):
            if (e[0] * v[0] + e[1] * v[1] + e[2] * v[2] + e[3] * v[3]) % MOD == 0:
                m |= 1 << i
                acc += s
        w = acc // 2
        mask_of[w][m] += 1
        sample[w].setdefault(m, v)

print("w ｜相異遮罩數 ｜總簽名")
out = {}
for w in sorted(mask_of):
    tot = sum(mask_of[w].values())
    print(f"{w:+d} ｜{len(mask_of[w]):10d} ｜{tot}")
    out[str(w)] = {"n_masks": len(mask_of[w]),
                   "masks": [{"mask_bits": format(m, "08b"), "count": c, "example": list(sample[w][m])}
                             for m, c in mask_of[w].most_common()],
                   "total": tot}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/wmask.json", "w"), indent=1, ensure_ascii=False)
print("\n已寫 /home/wsl/idcm_ckpt/e8/mrB/wmask.json")
for w in (2, 1, -1):
    print(f"\n=== w={w:+d} 之遮罩（前 5）")
    for m, c in mask_of[w].most_common(5):
        print("   mask {}｜count {:6d}｜例 {}".format(format(m, "08b"), c, sample[w][m]))
