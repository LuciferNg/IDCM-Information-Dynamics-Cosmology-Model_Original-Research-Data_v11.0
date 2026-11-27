#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R1：四腿頂點遮罩計數之獨立第二法（有序 44⁴ 全枚舉 → 多重集）＋ 計數之因數結構
   R5：m 之 2D 分位結構（gap* × 對比 之格）
   兩者共用一次載入，輸出 json
"""
import itertools
import json
import os
from collections import Counter, defaultdict

N, MOD = 44, 90
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]


def w_of(rem, a1, a2, a3):
    v = (rem, a1, a2, a3)
    return sum(c for c, sg in zip(CS, SIG)
               if (v[0] * sg[0] + v[1] * sg[1] + v[2] * sg[2] + v[3] * sg[3]) % MOD == 0)


def mask_of(rem, a1, a2, a3):
    v = (rem, a1, a2, a3)
    return tuple(i for i, sg in enumerate(SIG)
                 if (v[0] * sg[0] + v[1] * sg[1] + v[2] * sg[2] + v[3] * sg[3]) % MOD == 0)


# ---- R1：有序全枚舉 → 多重集去重 ----
seen = {}
order_cnt = Counter()      # 有序 (rem, a1,a2,a3) 之 W 分佈
orbit = Counter()          # 每個 4-向量之置換軌道大小
for rem in range(1, N + 1):
    for a1 in range(1, N + 1):
        for a2 in range(1, N + 1):
            for a3 in range(1, N + 1):
                key = (rem,) + tuple(sorted((a1, a2, a3)))
                w = w_of(rem, a1, a2, a3)
                order_cnt[w] += 1
                seen[key] = (w, mask_of(rem, a1, a2, a3))
                orbit[key] += 1
per = Counter(orbit.values())
ms_cnt = Counter(w for w, _ in seen.values())
mask_by_w = defaultdict(Counter)
for w, m in seen.values():
    mask_by_w[w][m] += 1

out = {
    "method": "有序 44^4 全枚舉 → 依 (rem, sorted adds) 去重",
    "n_ordered": N ** 4, "n_multiset": len(seen),
    "multiset_w_counts": {str(k): ms_cnt[k] for k in sorted(ms_cnt)},
    "ordered_w_counts": {str(k): order_cnt[k] for k in sorted(order_cnt)},
    "orbit_size_hist": {str(k): per[k] for k in sorted(per)},
    "n_masks": {str(w): len(mask_by_w[w]) for w in sorted(mask_by_w)},
    "masks_per_w": {str(w): {",".join(map(str, m)): c for m, c in mask_by_w[w].most_common()} for w in sorted(mask_by_w)},
}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/r1_masks_ordered.json", "w"), indent=1, ensure_ascii=False)
print("R1 多重集 W 計數:", dict(sorted(ms_cnt.items())))
print("R1 有序 W 計數:", dict(sorted(order_cnt.items())))
print("R1 軌道大小分佈:", dict(sorted(per.items())))
print("R1 各 W 之遮罩數:", dict(sorted(out["n_masks"].items())))

# ---- R5：m 之 2D 分位結構 ----
import numpy as np
d = json.load(open("/home/wsl/idcm_ckpt/e8/mrB/dominance.json"))
R = d["rows"]
x = np.log([r["cz_over_cu"] for r in R]); g = np.log([r["gap"] for r in R])
m = np.exp(np.log([r["E2"] for r in R]) - g - x)
q = lambda a, p: float(np.percentile(a, p))
grid = {}
for lo, hi, name in ((-9, 2, "low"), (2, 6, "mid"), (6, 9, "high")):
    sel = (x >= lo) & (x < hi)
    if sel.sum() >= 3:
        grid[name] = {"n": int(sel.sum()), "m_median": float(np.median(m[sel])),
                      "m_min": float(m[sel].min()), "m_max": float(m[sel].max())}
out5 = {"m_quantiles": {f"p{p}": q(m, p) for p in (5, 10, 25, 50, 75, 90, 95)},
        "contrast_bins": grid,
        "n_lt1": int((m < 1).sum()), "n_lt10": int((m < 10).sum()),
        "spearman_like": "見 docs/49 §二十八 corr(log gap, log|E2|) = −0.036817594735969815"}
json.dump(out5, open("/home/wsl/idcm_ckpt/e8/mrB/r5_distribution.json", "w"), indent=1, ensure_ascii=False)
print("R5 m 分位:", {k: round(v, 6) for k, v in out5["m_quantiles"].items()})
print("R5 對比分箱:", json.dumps(grid, ensure_ascii=False))
