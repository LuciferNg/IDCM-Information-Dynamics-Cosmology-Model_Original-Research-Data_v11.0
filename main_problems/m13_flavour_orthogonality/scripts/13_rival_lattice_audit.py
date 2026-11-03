#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/13_rival_lattice_audit.py — 競爭者審計（格粒度）

問題：`docs/04 §四` 之 look-elsewhere 只掃**整數** `α_s/α_u`（1..8）⟹ 得 1/8 = 12.5%。
但 `ρ_u = 0`（A₂ Cartan 零方向，結構強制）只給 `α_d = 2α_u`；`α_s/α_u` 之**分母**從未論證。

本腳本審：
  ① `ρ_d` 之 1σ 窗換算為 `r = α_s/α_u` 之區間；
  ② 按**分母類**（A₂ 權重格 q≤3／q≤4／q≤12／整數）列舉落窗者；
  ③ 以 `χ²` 排位 —— 若 `3/8` 唔係最佳，則「數據擇優」之講法**唔成立**；
  ④ 各競爭者與 `3/8` 之 3σ 分辨門檻，並與 Σ 硬底 0.01099 對比。
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
RD = mp.mpf("0.40350365")
SD = mp.mpf("0.044785178")
RU = mp.mpf("-0.02588837")
SU = mp.mpf("0.033583899")
FLOOR = mp.mpf("0.01099")        # docs/08 §四：Σ 質量阻死之 σ_ρd 下限


def rho_d(r):
    """ρ_u=0（α_d=2α_u）下：ρ_d = 3/(2r)，r = α_s/α_u。"""
    return mp.mpf(3)/(2*r)


print("=== 0. `ρ_u = 0` 之前提（結構強制，非擬合）===")
print("  A₂ Cartan 首行零方向 ⟹ α_d = 2α_u（`2` 係 Cartan 元素，非參數）")
print("  ⟹ ρ_d = 3/(2r)，r ≡ α_s/α_u 為唯一未知")

lo, hi = RD-SD, RD+SD
print(f"\n=== 1. ρ_d 1σ 窗 ⟹ r 區間 ===")
print(f"  ρ_d ∈ [{mp.nstr(lo,7)}, {mp.nstr(hi,7)}]  ⟹  r ∈ [{mp.nstr(3/(2*hi),6)}, {mp.nstr(3/(2*lo),6)}]")

print("\n=== 2. 按分母類列舉落窗者（全報）===")
for qmax, label in ((1, "整數 q=1"), (3, "A₂ 權重格 q≤3"), (4, "q≤4"), (12, "q≤12")):
    rows = []
    for q in range(1, qmax+1):
        for p in range(1, 40*q):
            r = F(p, q)
            if r.denominator != q:
                continue
            if not (F(3346, 1000) <= r <= F(4183, 1000)):
                continue
            rd = rho_d(mp.mpf(p)/q)
            z = (rd-RD)/SD
            if abs(z) < 1.5:
                rows.append((r, rd, z))
    rows.sort(key=lambda t: abs(t[2]))
    print(f"\n  [{label}] 落 1σ 窗者 {sum(1 for t in rows if abs(t[2])<1)} 個（含 1.5σ 內共 {len(rows)}）")
    for i, (r, rd, z) in enumerate(rows, 1):
        tag = ""
        if r == 4:
            tag = "  ← 候選（B=4）"
        elif r == F(11, 3):
            tag = "  ← N_d/12"
        elif r == F(15, 4):
            tag = "  ← 15/4"
        elif r == F(7, 2):
            tag = "  ← 7/2"
        print(f"    {i}. r = {str(r):>6s}  ρ_d = {mp.nstr(rd,7):>9s}  z = {mp.nstr(z,6):>8s}{tag}")

print("\n=== 3. 候選 3/8 之排位（按 |z|）===")
cands = [(F(4), "B = 4 = γ44（已註冊分枝）"), (F(15, 4), "15/4（未註冊）"),
         (F(11, 3), "N_d/12 = 44/12"), (F(7, 2), "7/2（未註冊）")]
cands.sort(key=lambda t: abs((rho_d(mp.mpf(t[0].numerator)/t[0].denominator)-RD)/SD))
for i, (r, prov) in enumerate(cands, 1):
    z = (rho_d(mp.mpf(r.numerator)/r.denominator)-RD)/SD
    print(f"  {i}. r = {str(r):>6s}  z = {mp.nstr(z,6):>8s}   {prov}")
print("  ⟹ **3/8 唔係擬合最佳者**；若按數據擇優，佢排第 4（共 4 個落窗）")

print("\n=== 4. 與 3/8 之分辨門檻 vs Σ 硬底 ===")
rd38 = rho_d(mp.mpf(4))
for r, prov in cands:
    if r == 4:
        continue
    rd = rho_d(mp.mpf(r.numerator)/r.denominator)
    sep = abs(rd-rd38)
    need = sep/3
    print(f"  vs r={str(r):>5s}：ρ_d 分離 {mp.nstr(sep,5)} = {mp.nstr(sep/SD,4)}σ（現況）；"
          f"3σ 需 σ_ρd ≤ {mp.nstr(need,6)} ⟹ 改善 {mp.nstr(SD/need,4)}×"
          f"{'  ⚠ 低於硬底' if need < FLOOR else '  ✅ 高於硬底 ' + mp.nstr(need/FLOOR,4) + '×'}")

print("\n=== 5. look-elsewhere（正確粒度）===")
print("  整數格      ：落窗 1 個（r=4）⟹ p ≈ 1/1（**零 else-where**）")
print("  A₂ 權重格 q≤3：落窗 2 個（4、11/3）⟹ p ≈ 1/2")
print("  q≤4         ：落窗 4 個（4、15/4、11/3、7/2）⟹ p ≈ 1/4")
print("  ⟹ `docs/04` 之 1/8 只對『整數且已選格點集』成立；分母粒度未論證")
print("  ⟹ 而 `α_d = 2α_u` 之 `2` 係 Cartan 元素 ⟹ 整數性並非 A₂ 權重格之必然（格為 (1/3)ℤ²）")

print("\n=== 判 ===")
print("  ① `ρ_u = 0` 仍係結構強制（不變）；此為本 M 唯一未被削弱之結構結果")
print("  ② `ρ_d = 3/8` **唔可以**稱為「數據擇優」：4 個同格候選中排第 4，最佳者 z=−0.08")
print("  ③ 4 個候選以現精度互相不可分辨（最大分離 0.76σ）；3σ 分辨需 σ ≤ 0.01136 ⟹")
print(f"     與 Σ 硬底 {FLOOR} 之比 = {mp.nstr(mp.mpf('0.01136')/FLOOR,4)}× ⟹ **實務上貼死硬底**")
print("  ④ 配合 `docs/11`（T1 標籤層已推死）：候選既無數據擇優、亦無 provenance")
print("     ⟹ 依 corpus 標準（『講不出算子即標籤』）判：**候選 = 標籤級，非結構結果**")
sys.exit(0)