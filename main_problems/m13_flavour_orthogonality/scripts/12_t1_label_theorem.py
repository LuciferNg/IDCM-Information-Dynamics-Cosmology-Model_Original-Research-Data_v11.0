#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/12_t1_label_theorem.py — T1（結構出處）之標籤層完備性定理

背景：`docs/09 §五 R2` 曾把 T1 定性為「類型性不可推盡（provenance 無界）」。
本輪把「無界」拆成兩個域，並在**合法域**內推死：

  域 ①（非法）：任意表達式搜尋 —— 無界，但唔符合 corpus 標準（算術 ≠ 算子身份）
  域 ②（合法）：α 由**已註冊標籤**生成 —— 有限，且**可完全判定**

標籤語義（出處：M4 `docs/32 §一` 引 v9 `03_extended_derivation.md §D.1` + `paper/17` L217–231）：
  `k ∈ {1,4,5}` = **扇區**：k=1 → Up sector；k=4 → Down sector；k=5 → Lepton sector
  `d ∈ {3,2,0}` = **Z₂ grading = 世代**（三代）
⟹ 讀出三點落喺 3×3 標籤網格：
      u = (k=1, d=3)      ← 上型、第一代
      d = (k=4, d=3)      ← 下型、第一代（與 u **同代**）
      s = (k=4, d=2)      ← 下型、第二代（與 d **同扇區**）

定理（三個 case，全部精確）：
  A. α 只依 `k`：u 獨佔 k=1，d、s 同 k=4 ⟹ `α_d ≡ α_s`
     ρ_u = 0 ⟹ g₁ = g₄/2 ⟹ **ρ_d = 3/4 強制**（與 g₄ 無關）⟹ 對比觀測 z = +7.7σ ⟹ 死
  B. α 只依 `d`：u、d 同 d=3 ⟹ `α_u ≡ α_d`
     ρ_u = f₃/(2f₂) ⟹ ρ_u = 0 ⟺ f₃ = 0 ⟹ **ρ_d = 0 強制** ⟹ z = −9.0σ ⟹ 死
  C. α = h(k,d) 一般函數（h 於 9 格任意）
     三讀出佔三格、彼此無函數約束 ⟹ h 於該三格自由 ⟹ projective dim = 2 ⟹ **dof = 2** ⟹ 違 T2
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
RD = mp.mpf("0.40350365")
RU = mp.mpf("-0.02588837")
SD = mp.mpf("0.044785178")
SU = mp.mpf("0.033583899")


def rho(au, ad, as_):
    return (2*au-ad)/(2*as_), (2*ad-au)/(2*as_)


def z(ru, rd):
    zu = (mp.mpf(ru.numerator)/ru.denominator - RU)/SU
    zd = (mp.mpf(rd.numerator)/rd.denominator - RD)/SD
    return zu, zd


print("=== 0. 標籤網格與讀出三點 ===")
print("  k ∈ {1,4,5} = 扇區 (Up / Down / Lepton)")
print("  d ∈ {3,2,0} = Z₂ grading = 世代")
print("  u = (1,3)   d = (4,3)   s = (4,2)")
print("  ⟹ u 與 d **同代**（d 相同）；d 與 s **同扇區**（k 相同）")

print("\n=== A. α 只依 k：α_u = g(1), α_d = α_s = g(4) ===")
g1, g4 = mp.mpf("0.5"), mp.mpf("1")     # 任意取 g4=1, g1 由 ρ_u=0 反解（定理形式，非擬合）
au, ad, as_ = g1, g4, g4
ru, rd = rho(mp.mpf(au), mp.mpf(ad), mp.mpf(as_))
print(f"  ρ_u = {mp.nstr(ru,10)}  （已施 ρ_u = 0 ⟹ g₁ = g₄/2）")
print(f"  ρ_d = 3/4 = 0.75 **強制**，與 g₄ 之標度無關")
print(f"  z = (ρ_u, ρ_d) − 觀測 = ({mp.nstr((ru-RU)/SU,4)}, {mp.nstr((mp.mpf('0.75')-RD)/SD,4)}) σ")
print("  ⟹ 🔴 7.74σ 排除（且此值與 docs/04 §四 之 `α_s=2` 行一致，互相印證）")

print("\n=== B. α 只依 d：α_u = α_d = f(3), α_s = f(2) ===")
f3, f2 = mp.mpf("0"), mp.mpf("1")       # ρ_u = f3/(2f2) = 0 ⟹ f3 = 0
ru2, rd2 = rho(f3, f3, f2)
print(f"  ρ_u = f₃/(2f₂) ⟹ ρ_u = 0 ⟺ f₃ = 0；代回：ρ_d = {mp.nstr(rd2,10)}")
print(f"  z = ({mp.nstr((ru2-RU)/SU,4)}, {mp.nstr((rd2-RD)/SD,4)}) σ  ⟹ 🔴 −9.01σ 排除")
print("  （若 f₃ ≠ 0 則 ρ_u > 0，直接違反 ρ_u = 0 —— 兩難，無解）")

print("\n=== C. α = h(k,d) 一般函數：自由度計數 ===")
print("  h 於 9 格有 9 個值；三讀出佔 (1,3),(4,3),(4,2) 三格，彼此**無函數關係**")
print("  ⟹ 三格之值自由 ⟹ α 之 projective space（模標度）維度 = 3 − 1 = **2**")
print("  ⟹ dof = 2 ≠ 0 ⟹ 違反 T2（dof=0）⟹ 判死（無需枚舉）")

print("\n=== D. 候選 α = (1,2,4) 之來源定位 ===")
print("  (1,2,4) = 2^r，r = 讀出**次序**（0,1,2）—— 係**序數規則**，唔係標籤函數")
print("  依 docs/05 站無關性定理：換站／換味次序不影響 ρ ⟹ 序數係 **gauge**")
print("  ⟹ 候選之 α 無標籤層出處（T1① FAIL 屬結構性，非未搜到）")

print("\n=== E. 幾何性 ⟹ 零自由度（族內 ρ_u = 0 唯一釘死 x）===")
# 族 α=(1,x,x²)：ρ_u = (2−x)/(2x²) ⟹ ρ_u = 0 ⟺ x = 2 **精確**（無連續自由度）
ru_x = lambda x: (mp.mpf(2)-x)/(2*x*x)
print("  ρ_u(x) = (2−x)/(2x²) ⟹ ρ_u = 0 ⟺ **x = 2 精確**（族內唯一根，非擬合）")
for xv in ("1.5", "2", "2.5"):
    print(f"    x={xv}: ρ_u = {mp.nstr(ru_x(mp.mpf(xv)),6)}")
print("  ⟹ 幾何性 + ρ_u = 0 ⟹ α = (1,2,4) 唯一 ⟹ ρ_d = 3/8 **係預言，非擬合**")
print("  ⟹ look-elsewhere 只關於「族之選取」，**唔**關於 x（x 已被 ρ_u=0 釘死）")

print("\n=== E2. 反過來：若不假設幾何性，ρ_d 條件本身有幾鬆？ ===")
lo, hi = RD-SD, RD+SD
# ρ_d(x) = (2x−1)/(2x²) 在 [lo,hi] 內之 x 區間（解二次）
def xsolve(rho_val):
    a, b, c = 2*rho_val, mp.mpf(-2), mp.mpf(1)      # 2ρx² − 2x + 1 = 0
    disc = b*b - 4*a*c
    if disc <= 0:
        return None
    return ((-b - mp.sqrt(disc))/(2*a), (-b + mp.sqrt(disc))/(2*a))
print(f"  1σ 窗 ρ_d ∈ [{mp.nstr(lo,6)}, {mp.nstr(hi,6)}]")
for v in (lo, hi):
    r = xsolve(v)
    print(f"    ρ_d = {mp.nstr(v,6)} ⟹ x = {mp.nstr(r[0],5)} 或 {mp.nstr(r[1],5)}" if r else "    無實根")
cnt = 0
for den in range(1, 21):
    for num in range(1, 21*den):
        x = F(num, den)
        if x <= 0:
            continue
        rd = F(den*(2*num-den), 2*num*num)          # ρ_d = (2x−1)/(2x²)
        zz = (mp.mpf(rd.numerator)/rd.denominator - RD)/SD
        if abs(zz) < 1:
            cnt += 1
print(f"  分母 ≤ 20 之 x 中，**{cnt} 個**令 ρ_d 落 1σ 窗 ⟹ ρ_d 條件單獨幾乎無約束力")
print("  ⟹ 判別力全部來自 `ρ_u = 0`（族內唯一）；數值命中本身零資訊")

print("\n=== 判 ===")
print("  A 🔴：α 依 k ⟹ ρ_d 被強制為 3/4 ⟹ 7.74σ 排除")
print("  B 🔴：α 依 d ⟹ ρ_u = 0 逼令 ρ_d = 0 ⟹ 9.01σ 排除")
print("  C 🔴：α = h(k,d) ⟹ dof = 2 ⟹ 違 T2")
print("  ⟹ **標籤層（k,d）之任何函數皆不可能產生 (ρ_u, ρ_d) = (0, 3/8)**")
print("  ⟹ T1① 之不可滿足係**定理**，唔係「未搜夠」；範圍 = 標籤網格之全函數族（非 6 規則）")
print("  ⟹ 域外（K_tree 側新符號）已被 M4 docs/32 裁為不存在 ⟹ T1 於 v11 存貨內**空集**")
sys.exit(0)