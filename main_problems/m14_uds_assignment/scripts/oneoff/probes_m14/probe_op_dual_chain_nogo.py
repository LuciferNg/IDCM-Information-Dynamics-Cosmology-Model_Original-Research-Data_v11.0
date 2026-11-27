#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""雙鏈對偶判定：四腿頂點（Dirichlet L=90）可否接 W-field 勢（環 C=88）？
   ① 譜重合檢定：{4sin²(kπ/88)} ∩ {4sin²(nπ/90)}（k,n = 1..44）是否有非平凡交點
   ② 和恆等式：環(M) Σγ = 2M+2；Dirichlet(L) Σγ = L−2（精確 Fraction）
   ③ ★四分之一週期模存在性：μ² = 2 需 γ_{C/4} = 4sin²(π/4) = 2 ⟹ 4 | C
        環 C=88 ⟹ 88/4 = 22 ∈ ℤ ✓；Dirichlet L=90 ⟹ 90/4 = 22.5 ∉ ℤ ✗
   ④ 判詞：不可譯 ⟹ 出口 ⑤ 關閉（no-go）
"""
import json
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 40

# ① 譜重合
R = [4 * mp.sin(mp.pi * mp.mpf(k) / 88) ** 2 for k in range(1, 45)]
D = [4 * mp.sin(mp.pi * mp.mpf(n) / 90) ** 2 for n in range(1, 45)]
tol = mp.mpf("1e-30")
overlaps = [(k + 1, n + 1) for k in range(44) for n in range(44) if abs(R[k] - D[n]) < tol]
print("① 譜重合（1e−30）：", overlaps if overlaps else "無非平凡交點")

# 有理比論證：γ 相等 ⟹ k/88 = n/90（或互補）⟹ 45k = 44n ⟹ k=44m, n=45m
print("   有理比論證：45k = 44n ⟹ (k,n) = (44m, 45m)；k≤44、n≤44 ⟹ 只有 m=0 ⟹ 無解")

# ② 和恆等式（精確 Fraction 版本：以 cos 之和證明）
# Σ_{n=1}^{M} cos(nπ/M) = −1 ⟹ Σγ = 2M − 2Σcos = 2M + 2（環）
# Dirichlet：Σ_{n=1}^{M-1} cos(nπ/M) = 0 ⟹ Σγ = 2(M−1)
def sum_cos(M, upto):
    return sum(mp.cos(mp.pi * mp.mpf(n) / M) for n in range(1, upto + 1))

# γ_k = 2 − 2cos(2θ_k)；環 θ_k = kπ/88 ⟹ 2θ = kπ/44；Dirichlet θ_n = nπ/90 ⟹ 2θ = nπ/45
sR = 2 * 44 - 2 * sum(mp.cos(mp.pi * mp.mpf(k) / 44) for k in range(1, 45))
sD = 2 * 44 - 2 * sum(mp.cos(mp.pi * mp.mpf(n) / 45) for n in range(1, 45))
print("② Σγ(環 C=88) =", mp.nstr(sR, 20), "（應 90；Σcos(kπ/44) = −1）")
print("   Σγ(Dir L=90) =", mp.nstr(sD, 20), "（應 88；Σcos(nπ/45) = 0）")

# ③ 四分之一週期模
print("③ 4 | 88 ?", 88 % 4 == 0, "（88/4 =", 88 // 4, "）｜γ_22(/88) =", mp.nstr(R[21], 20))
print("   4 | 90 ?", 90 % 4 == 0, "（90/4 =", mp.nstr(mp.mpf(90) / 4, 8), "）⟹ 無整數模 ⟹ 無 γ = 2")

out = {"spectral_overlap": overlaps, "overlap_nontrivial": bool(overlaps),
       "rational_ratio_obstruction": "45k = 44n ⟹ (k,n) = (44m, 45m) ⟹ m=0 only",
       "sum_gamma_ring_88": mp.nstr(sR, 20), "sum_gamma_dirichlet_90": mp.nstr(sD, 20),
       "div4_88": 88 % 4 == 0, "quarter_mode_ring": 88 // 4,
       "div4_90": 90 % 4 == 0, "quarter_mode_dirichlet": None,
       "verdict": ("不可譯（no-go）：兩譜無非平凡交點（45k=44n 無解）；μ²=2 需四分之一週期模 ⟹ 4|C，"
                   "而 4∤90 ⟹ Dirichlet 鏈無 γ=2 ⟹ λ 無對應物 ⟹ 四腿頂點（mod 90）不可接 W-field 勢"
                   if not overlaps and 88 % 4 == 0 and 90 % 4 != 0 else "需進一步分析")}
with open("/home/wsl/idcm_ckpt/e8/mrB/op_dual_chain_nogo.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("④ 判詞:", out["verdict"])
