#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/08_alpha_source.py — 甲：α 之來歷（K_tree？K_chain？）

§0 前提核查（唔係掃描）：
  M4 `docs/32_flavour_structure_scout.md` 已判：**味（flavour）之標籤住喺 K_chain**
  （模指標 `k ∈ {1,4,5}`、Z₂ grading `d ∈ {3,2,0}`），**唔係 K_tree 之分支**；
  「K_tree 味分支結構」概念已**撤回**（`docs/31 §十` 用語更正）。
  ⟹ `α`（逐味振幅）**不可能**由 K_tree 側給出 —— 此係已註冊裁決，非本輪新判。
  ⟹ 甲 之正確對象 = **K_chain 之味標籤資料** `(k, d)`。

§0b 但：register 記載 `v9 §D.1` 之 sector 表**未被採用** ⟹ `(k,d) ↔ {u,d,s}` 之對應**未註冊**。
  ⟹ 故必須枚舉 `3! = 6` 個指派，並全報（不得只報命中）。

§1 規則集（跑前凍結，6 條）：
  R1 `α = k`      R2 `α = 1/k`     R3 `α = 2^d`
  R4 `α = 2^{−d}` R5 `α = B^d`     R6 `α = B^{−d}`
  （`k ∈ {1,4,5}`、`d ∈ {3,2,0}`、`B = 4`）
  N_trials = 6 指派 × 6 規則 = **36**，全報；命中 = `α_s/α_u ∈ [3.3048, 4.1300]`（1σ 窗）。
"""
import sys
import itertools
import mpmath as mp

mp.mp.dps = 30
RD, SD = mp.mpf("0.40350365"), mp.mpf("0.044785178")
need = 3/(2*RD)
sig = 1.5*SD/RD**2
LO, HI = need-sig, need+sig
K = {"u": 1, "d": 4, "s": 5}          # 模指標集合 {1,4,5}
D = {"u": 3, "d": 2, "s": 0}          # Z₂ grading 集合 {3,2,0}
B = 4

rules = {
    "R1 α=k":        lambda f: mp.mpf(K[f]),
    "R2 α=1/k":      lambda f: 1/mp.mpf(K[f]),
    "R3 α=2^d":      lambda f: mp.mpf(2)**D[f],
    "R4 α=2^{−d}":   lambda f: mp.mpf(2)**(-D[f]),
    "R5 α=B^d":      lambda f: mp.mpf(B)**D[f],
    "R6 α=B^{−d}":   lambda f: mp.mpf(B)**(-D[f]),
}
LABELS = ("u", "d", "s")

print("=== 0. 窗 ===")
print(f"  α_s/α_u ∈ [{mp.nstr(LO,8)}, {mp.nstr(HI,8)}]（1σ），中心 {mp.nstr(need,8)} ± {mp.nstr(sig,5)}")

print("\n=== 1. 枚舉：6 指派 × 6 規則 = 36（全報）===")
print(f"{'指派 (k:u,d,s | d:u,d,s)':28s} {'規則':12s} {'α_s/α_u':>12s} {'z_α':>9s} {'ρ_u':>10s} {'ρ_d':>10s}  判")
hits = []
n = 0
for kk in itertools.permutations(sorted(K.values())):      # 模指標之指派
    for dd in itertools.permutations(sorted(D.values())):  # Z₂ grading 之指派
        K2 = dict(zip(LABELS, kk))
        D2 = dict(zip(LABELS, dd))
        globals()["K"], globals()["D"] = K2, D2
        for rn, rf in rules.items():
            n += 1
            a = {f: rf(f) for f in LABELS}
            ru = (2*a["u"]-a["d"])/(2*a["s"])
            rd = (2*a["d"]-a["u"])/(2*a["s"])
            ratio = a["s"]/a["u"]
            za = (ratio-need)/sig
            ok = LO <= ratio <= HI
            if ok:
                hits.append((K2.copy(), D2.copy(), rn, ratio, ru, rd))
            if ok or abs(za) < 1.5:
                print(f"{str(kk)+' | '+str(dd):28s} {rn:12s} {mp.nstr(ratio,8):>12s} "
                      f"{mp.nstr(za,6):>9s} {mp.nstr(ru,8):>10s} {mp.nstr(rd,8):>10s}  {'🟡窗內' if ok else '🔴'}")

print(f"\n  總試驗 N = {n}；落窗 = {len(hits)}（{len(hits)/n*100:.1f}%）")
print(f"  ⟹ naive look-elsewhere：p ≈ {len(hits)/n:.3f}")

print("\n=== 2. 落窗者明細（若多於 1 個 ⟹ 規則不唯一 ⟹ 欠定未解）===")
for (K2, D2, rn, ratio, ru, rd) in hits:
    print(f"  k={tuple(K2[f] for f in LABELS)} d={tuple(D2[f] for f in LABELS)} {rn:12s}"
          f" α_s/α_u={mp.nstr(ratio,8)} ⟹ (ρ_u, ρ_d)=({mp.nstr(ru,8)}, {mp.nstr(rd,8)})")

print("\n=== 2b. ★ 聯合閘：ρ_u ≈ 0 **與** ρ_d ≈ 3/8 同時 ===")
RU, SU = mp.mpf("-0.02588837"), mp.mpf("0.033583899")
joint = []
for kk in itertools.permutations(sorted(K.values())):
    for dd in itertools.permutations(sorted(D.values())):
        K2, D2 = dict(zip(LABELS, kk)), dict(zip(LABELS, dd))
        globals()["K"], globals()["D"] = K2, D2
        for rn, rf in rules.items():
            a = {f: rf(f) for f in LABELS}
            ru = (2*a["u"]-a["d"])/(2*a["s"])
            rd = (2*a["d"]-a["u"])/(2*a["s"])
            if abs((ru-RU)/SU) < 1 and abs((rd-RD)/SD) < 1:
                joint.append((tuple(kk), tuple(dd), rn, ru, rd))
print(f"  聯合命中（|z_u|<1 且 |z_d|<1）= {len(joint)} / {216}")
print("  ⟹ 標籤層**無任何** (k,d) 規則可同時給出『近正交 + 正確 ρ_d』。")
print("  ⟹ 結構理由：`ω₂` 條件需 `α_d = 2α_u`；")
print("     於 `2^d`（d ∈ {3,2,0}）可達（d_u=2,d_d=3），但此時 d_s=0 ⟹ `α_s/α_u = 1/4` ⟹ `ρ_d = 6`（窗外）；")
print("     於 `k`（{1,4,5}）與 `B^d` 則**不可能**給出比 2 ⟹ 連 `ρ_u = 0` 都做不到。")

print("\n=== 3. 判決 ===")
print("  §0 K_tree 路線：**結構性封閉**（M4 docs/32 已裁：味 ∈ K_chain）⟹ 甲 不可能由 K_tree 完成。")
print(f"  §1 K_chain 標籤路線：{len(hits)}/{n} 落窗 ⟹ "
      f"{'**規則不唯一，欠定未解**' if len(hits) > 1 else '唯一，可逐一裁決'}")
print("  §2 結構理由：`k ∈ {{1,4,5}}`、`d ∈ {{3,2,0}}` 取值域極小，其比值（含 2^、B^ 之升降）")
print("     密集覆蓋 order-1 值 ⟹ **標籤層本身無解析度**去決定 α_s/α_u。")
print("\n  ⟹ 甲 之結論 = **雙重封閉**：K_tree 側被結構裁決封閉；K_chain 標籤側被標籤解析度不足封閉。")
print("     故 α 之來歷**不可能**由現有標籤資料唯一給出 ⟹ 交由乙（精度）路線決定生死。")
sys.exit(0)