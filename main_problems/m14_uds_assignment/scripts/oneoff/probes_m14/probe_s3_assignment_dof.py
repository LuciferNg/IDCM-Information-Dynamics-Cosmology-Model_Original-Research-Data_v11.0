#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""S3：「味站指派」良定義性之量化判定
   ① 指派自由度：把 k 個味注入 44 個站之單射數
   ② 預言散佈：對全部指派，最輕/最重比 (γ_b/γ_a) 之範圍（若跨數量級 ⟹ 無約束力）
   ③ Aut(H) = {e}（γ 嚴格遞增）⟹ 無對稱原則可選指派 ⟹ dof 無法由結構消去
"""
import itertools
import json

import mpmath as mp

mp.mp.dps = 30
G = [4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45)]


def n_inject(k):
    return mp.binomial(44, k) * mp.factorial(k)


out = {}
for k in (3, 6):
    out[f"n_assign_{k}"] = mp.nstr(n_inject(k), 20)
    print(f"k={k} 單射數 = {mp.nstr(n_inject(k), 20)}（= C(44,{k})·{k}!）")

# ② 三味（如 e,μ,τ）之預言散佈：取站 (a<b<c)，預言 m_c/m_a = sqrt(γ_c/γ_a)
lo, hi = None, None
cnt = 0
for a, b, c in itertools.combinations(range(1, 45), 3):
    r = mp.sqrt(G[c - 1] / G[a - 1])
    lo = r if lo is None else min(lo, r)
    hi = r if hi is None else max(hi, r)
    cnt += 1
out["triples"] = cnt
out["pred_ratio_min"] = mp.nstr(lo, 12)
out["pred_ratio_max"] = mp.nstr(hi, 12)
out["pred_span_decades"] = mp.nstr(mp.log10(hi / lo), 8)
print(f"三味組合 {cnt}：預言 m_c/m_a 之範圍 {mp.nstr(lo,8)} … {mp.nstr(hi,8)}（跨 {mp.nstr(mp.log10(hi/lo),5)} 個數量級）")

# ③ Aut(H)：γ 嚴格遞增 ⟹ 交換 H=diag(γ) 之置換只有恆等
strict = all(G[i] < G[i + 1] for i in range(43))
out["gamma_strictly_increasing"] = strict
out["aut_trivial"] = strict
print("γ 嚴格遞增:", strict, "⟹ Aut(H) = {e}")

out["verdict"] = ("ill-posed：指派自由度 79,464（k=3）／5,082,517,440（k=6）且 Aut(H)={e} 無對稱消去；"
                  "預言比散佈跨 " + mp.nstr(mp.log10(hi / lo), 5) + " 個數量級 ⟹ 無約束力")
print("判詞:", out["verdict"])
with open("/home/wsl/idcm_ckpt/e8/mrB/s3_assignment_dof.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("→ 落 mrB/s3_assignment_dof.json")
