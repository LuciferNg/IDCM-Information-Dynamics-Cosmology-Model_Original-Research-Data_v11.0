#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""物理應用：『m ∝ γ^p ＋ 標記指派』是否可能同時容納觀測質量？（修正版）
   每個質量對 (hi,lo) 只需【存在】一對站 (a,b) 使 (m_hi/m_lo)^{1/p} = γ_b/γ_a（容差 0.1%）
   計分：被容納之質量對數 / 總對數；同時要求站不重複（同一 p 下作 greedy 注入檢查）
"""
import itertools
import json

import mpmath as mp

mp.mp.dps = 40
G = [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]
GRlist = [((a, b), G[b - 1] / G[a - 1]) for a in range(1, 45) for b in range(a + 1, 45)]

SECT = {
    "lepton": {"e": "0.00051099895", "mu": "0.1056583755", "tau": "1.77686"},
    "up": {"u": "0.00216", "c": "1.27", "t": "172.69"},
    "down": {"d": "0.00467", "s": "0.0934", "b": "4.18"},
    "all9": None,
}
SECT["all9"] = {**SECT["lepton"], **SECT["up"], **SECT["down"]}
TOL = mp.mpf("0.001")


def analyse(masses):
    Mv = {k: mp.mpf(v) for k, v in masses.items()}
    names = sorted(Mv, key=lambda k: Mv[k])                 # 由輕到重
    pairs = [(names[i], names[j]) for i in range(len(names)) for j in range(i)]  # (hi, lo)，names 升序
    cand = []
    for hi, lo in pairs:
        mr = Mv[hi] / Mv[lo]
        for _, gr in GRlist:
            if gr > 1:
                cand.append(mp.log(mr) / mp.log(gr))
    best = (-1, None, [], 0)
    for p in cand:
        if p <= 0 or p > 4:
            continue
        matched = []
        for hi, lo in pairs:
            need = (Mv[hi] / Mv[lo]) ** (1 / p)
            ok = [(a, b) for (a, b), gr in GRlist if abs(gr / need - 1) < TOL]
            if ok:
                matched.append((hi, lo, ok[:2]))
        # 站重用檢查（貪心：每個質量對取第一個不衝突之站對）
        used, assign = set(), []
        for hi, lo, ok in matched:
            pick = next(((a, b) for a, b in ok if a not in used and b not in used), None)
            if pick:
                used |= set(pick)
                assign.append((hi, lo, pick[0], pick[1]))
        score = (len(matched), len(assign))
        if score > (best[0], best[3]):
            best = (len(matched), p, matched, len(assign), assign)
    return {"n_pairs": len(pairs), "n_matched": best[0], "best_p": mp.nstr(best[1], 12) if best[1] else None,
            "n_injectable_no_reuse": best[3], "assignment": best[4] if len(best) > 4 else [],
            "pairs": [f"{a}>{b}" for a, b in pairs]}


out = {}
for k, v in SECT.items():
    out[k] = analyse(v)
    r = out[k]
    print(f"{k}: 對 {r['n_pairs']}｜容納 {r['n_matched']}｜可無重用注入 {r['n_injectable_no_reuse']}｜p = {r['best_p']}")
    for a in r["assignment"]:
        print(f"    {a[0]:>4}/{a[1]:<4} → 站 {a[2]},{a[3]}")
out["tolerance"] = "0.001"
out["verdict"] = ("單一冪 p ＋ 標記：不可能同時容納（容納對數與可注入站數皆不足）"
                  if out["all9"]["n_injectable_no_reuse"] < out["all9"]["n_pairs"] else "全對可容納（需檢注入唯一性）")
json.dump(out, open("/home/wsl/idcm_ckpt/e8/mrB/phys_power_fit.json", "w"), indent=1, ensure_ascii=False)
print("判詞:", out["verdict"])
