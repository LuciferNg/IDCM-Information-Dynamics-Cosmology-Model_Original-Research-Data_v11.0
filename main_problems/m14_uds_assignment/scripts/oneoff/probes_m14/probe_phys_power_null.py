#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""單一冪假說之 null：隨機 9 質量（同跨度 log-uniform）下「可無重用注入」之期望值
   若觀測值 7 落於 null 分佈內側 ⟹ 該吻合純屬機會（判死）
"""
import json

import mpmath as mp
import numpy as np

mp.mp.dps = 30
Gmp = [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]
SB = [(a, b) for a in range(1, 45) for b in range(a + 1, 45)]
loggr = np.array([np.log(float(Gmp[b - 1] / Gmp[a - 1])) for a, b in SB])
grv = np.exp(loggr)
EPS = 0.001
LO, HI = np.log(0.00051099895), np.log(172.69)


def injectable(logm):
    pairs = [(i, j) for i in range(9) for j in range(i)]   # logm 升序
    logmr = np.array([logm[i] - logm[j] for i, j in pairs])
    S = logmr[:, None] / loggr[None, :]
    cand = np.unique(S.ravel())
    cand = cand[(cand > 0) & (cand <= 6)]
    cnt = np.zeros(cand.size, dtype=np.int32)
    for i in range(len(pairs)):
        row = np.sort(S[i])
        idx = np.searchsorted(row, cand)
        cnt += (np.abs(row[np.clip(idx - 1, 0, row.size - 1)] - cand) < EPS).astype(np.int32)
        cnt += (np.abs(row[np.clip(idx, 0, row.size - 1)] - cand) < EPS).astype(np.int32)
    k = int(cnt.argmax())
    pbest = float(cand[k])
    pair_ok = []
    for i, j in pairs:
        need = float(np.exp(logmr[i - j - 1] if False else (logm[i] - logm[j]) / pbest))
        ok = [(a, b) for (a, b), g in zip(SB, grv) if abs(g / need - 1) < EPS]
        pair_ok.append(ok)
    used, n = set(), 0
    for ok in pair_ok:
        pick = next(((a, b) for a, b in ok if a not in used and b not in used), None)
        if pick:
            used |= set(pick)
            n += 1
    return n, int(cnt[k])


REAL = np.log([0.00216, 0.00467, 0.0934, 1.27, 4.18, 172.69, 0.00051099895, 0.1056583755, 1.77686])
REAL.sort()
r_inj, r_cnt = injectable(REAL)

rng = np.random.default_rng(20261127)
N = 2000
inj = np.empty(N, dtype=np.int32)
cntm = np.empty(N, dtype=np.int32)
for t in range(N):
    lm = np.sort(rng.uniform(LO, HI, 9))
    inj[t], cntm[t] = injectable(lm)

out = {"real_injectable": int(r_inj), "real_max_matched_raw": int(r_cnt),
       "null_n": N, "null_inj_mean": float(inj.mean()), "null_inj_sd": float(inj.std(ddof=1)),
       "null_inj_p95": float(np.percentile(inj, 95)), "null_inj_max": int(inj.max()),
       "z": float((r_inj - inj.mean()) / inj.std(ddof=1)),
       "eps_log": EPS, "seed": 20261127}
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/phys_power_null.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(out, indent=1, ensure_ascii=False))
