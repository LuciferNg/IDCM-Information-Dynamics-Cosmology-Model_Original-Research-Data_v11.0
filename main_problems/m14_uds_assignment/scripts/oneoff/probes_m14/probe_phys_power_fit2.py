#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""m ∝ γ^p ＋ 標記：單一冪可能性（向量化精確搜尋）
   一致定義：|log(γ_b/γ_a) − p·log(m_hi/m_lo)| < eps（eps=0.001 ⟹ 相對 0.1%）
   候選 p = 所有 (質量對, 站對) 之 logmr/loggr
"""
import json

import mpmath as mp
import numpy as np

mp.mp.dps = 40
Gmp = [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]
SB = [(a, b) for a in range(1, 45) for b in range(a + 1, 45)]
loggr = np.array([np.log(float(Gmp[b - 1] / Gmp[a - 1])) for a, b in SB])
grv = np.exp(loggr)

SECT = {"lepton": {"e": "0.00051099895", "mu": "0.1056583755", "tau": "1.77686"},
        "up": {"u": "0.00216", "c": "1.27", "t": "172.69"},
        "down": {"d": "0.00467", "s": "0.0934", "b": "4.18"}}
SECT["all9"] = {**SECT["lepton"], **SECT["up"], **SECT["down"]}
EPS = 0.001


def analyse(masses):
    Mv = {k: mp.mpf(v) for k, v in masses.items()}
    names = sorted(Mv, key=lambda k: Mv[k])
    pairs = [(names[i], names[j]) for i in range(len(names)) for j in range(i)]
    logmr = np.array([float(mp.log(Mv[h] / Mv[l])) for h, l in pairs])
    S = logmr[:, None] / loggr[None, :]
    cand = np.unique(S.ravel())
    cand = cand[(cand > 0) & (cand <= 6)]
    cnt = np.zeros(cand.size, dtype=np.int32)
    for i in range(len(pairs)):
        row = np.sort(S[i])
        idx = np.searchsorted(row, cand)
        for off in (-1, 0):
            j = np.clip(idx + off, 0, row.size - 1)
            cnt += (np.abs(row[j] - cand) < EPS).astype(np.int32)
    k = int(cnt.argmax())
    pbest = float(cand[k])
    matched = []
    for h, l in pairs:
        need = float((Mv[h] / Mv[l]) ** (1 / mp.mpf(pbest)))
        ok = [(a, b) for (a, b), g in zip(SB, grv) if abs(g / need - 1) < EPS]
        if ok:
            matched.append((h, l, ok))
    used, assign = set(), []
    for h, l, ok in matched:
        pick = next(((a, b) for a, b in ok if a not in used and b not in used), None)
        if pick:
            used |= set(pick)
            assign.append((h, l, pick[0], pick[1]))
    return {"n_pairs": len(pairs), "max_matched": int(cnt[k]), "p_best": repr(pbest),
            "n_injectable_no_reuse": len(assign), "assignment": assign,
            "pairs": [f"{a}>{b}" for a, b in pairs], "n_candidates": int(cand.size)}


out = {}
for k, v in SECT.items():
    out[k] = analyse(v)
    r = out[k]
    print(f"{k}: 對 {r['n_pairs']}｜最大容納 {r['max_matched']}｜可無重用注入 {r['n_injectable_no_reuse']}｜p = {r['p_best']}")
    for a in r["assignment"]:
        print(f"     {a[0]:>4}/{a[1]:<4} → 站 {a[2]},{a[3]}")
out["eps_log"] = EPS
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/phys_power_fit.json", "w"), indent=1, ensure_ascii=False)
print("判詞:", "完整注入" if out["all9"]["n_injectable_no_reuse"] == out["all9"]["n_pairs"] else "不完整 ⟹ 單一冪不可能")
