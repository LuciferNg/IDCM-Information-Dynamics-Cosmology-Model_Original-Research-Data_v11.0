#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R2 —— 局域容錯架構：可共存之最大標籤位元數
2026-10-28

問題：若宇宙係局域容錯架構，則**任何標籤自由度**都必須不增加缺陷密度
      （缺陷集已飽和於 32 = #odd faces）。兩個 Z₂ 頂點着色 c₁, c₂ 之積為第三個着色；
      若 c₁ 已被 32-缺陷飽和佔用，則 c₂ 可共存 ⟺ c₁c₂ 亦飽和。

本腳本計算：族內兩兩乘積之飽和數 ⟹ 可共存之獨立 Z₂ 標籤位元數。

跑法：/usr/bin/python3 scripts/06_fault_tolerance_bits.py
"""
from __future__ import annotations
import importlib.util
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

spec = importlib.util.spec_from_file_location("ft", os.path.join(HERE, "05_fault_tolerance.py"))
ft = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ft)

L, OUT = [], {}


def say(s=""):
    L.append(s)
    print(s)


def main():
    P = ft.build_rh()
    adj, edges, ang = ft.graph(P)
    V = len(P)
    fam = {}
    for r in range(V):
        c = ft.bfs_colour(adj, r)
        fam.setdefault(ft.canon(c), c)
    keys = sorted(fam.keys())
    say("=" * 80)
    say("B1  族員 = %d（飽和於 32 缺陷）" % len(keys))

    def nfr(c):
        return sum(1 for e in edges if c[e[0]] == c[e[1]])

    assert all(nfr(fam[k]) == 32 for k in keys)
    say("    全族缺陷數 32/32 ✅")

    say("")
    say("=" * 80)
    say("B2  兩兩乘積 c_i·c_j 之飽和性（可共存判據）")
    hist = {}
    sat_pairs = []
    for i, ki in enumerate(keys):
        ci = fam[ki]
        for j, kj in enumerate(keys):
            cj = fam[kj]
            cp = [1 if ci[v] == cj[v] else 0 for v in range(V)]
            n = nfr(cp)
            hist[n] = hist.get(n, 0) + 1
            if n == 32:
                sat_pairs.append((i, j, n))
    say("    乘積缺陷數分佈 = %s" % dict(sorted(hist.items())))
    say("    飽和乘積數 = %d / %d" % (len(sat_pairs), len(keys) ** 2))
    # 判別：除了 i=j 與 j=complement(i) 之外的飽和乘積
    extra = [(i, j) for (i, j, n) in sat_pairs
             if i != j and ft.canon(tuple(1 - x for x in fam[keys[i]])) != keys[j]]
    say("    除 identity／全域翻轉之外之飽和乘積 = %d" % len(extra))
    say("")
    if len(extra) == 0:
        say("    ⟹ 族內**冇**第三個可共存之 Z₂ 標籤：")
        say("       c₁c₂ 飽和 ⟺ c₂ = c₁ 或 ¬c₁")
        say("    ⟹ 局域容錯要求下，A₅ 頂點只可供應 **1 個 bit**，且該 bit = 全域翻轉 = gauge")
        say("    ⟹ 味**唔可以**係第二個 Z₂ 頂點着色（會令缺陷 > 32 ⟹ 崩潰）")
        OUT["B2"] = {"product_hist": {str(k): v for k, v in hist.items()},
                     "saturated_products": len(sat_pairs),
                     "extra_saturated": 0,
                     "max_label_bits": 1, "gauge_only": True}
    else:
        say("    ⟹ 存在 %d 對額外可共存組合，需逐一檢查是否為 gauge。" % len(extra))
        OUT["B2"] = {"product_hist": {str(k): v for k, v in hist.items()},
                     "saturated_products": len(sat_pairs),
                     "extra_saturated": len(extra)}

    with open(os.path.join(ROOT, "data", "06_fault_tolerance_bits.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False)
    say("")
    say("[written] data/06_fault_tolerance_bits.json , data/run_06_fault_tolerance_bits.log")

    with open(os.path.join(ROOT, "data", "run_06_fault_tolerance_bits.log"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()