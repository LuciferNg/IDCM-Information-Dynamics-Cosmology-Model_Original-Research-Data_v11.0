#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/11_sector_extension.py — 「同一套規則是否適用 d,s,c,t（及其他扇區）」

三問：
  Q-A 結構：M13 之 α 規則（α_ℓ ∝ (log₂B)^{ℓ−1}）作用於**幾個槽**？可否加第 4 槽？
  Q-B 觀測：重味（c、b、t）有無可算之 δ？現值為何？可否裁決？
  Q-C 判：同一套規則「適用」嘅三種意思分別判生／死。

數據（全部來自已註冊／已 vendored，零新量測）：
  M12 `docs/14` R5 二體荷（單位 MeV，唯讀引用）：
    δ_u = −0.0887 、 δ_d = +1.3820 、 δ_s = +3.4250（三錨唯一解，無 dof）
    δ_c = −0.2820 ± 0.2830（由 Σc⁰−Σc⁺ 解）
    Σb⁻−Σb⁺、Ξb⁻−Ξb⁰ 只有組合命中（C1 z = +0.27），δ_b **未隔離**
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
# 精確三錨（M13 scripts/01 之精確解）
DU = F(-4433383, 50000000)
DD = F(691, 500)
DS_ = F(137, 40)
DC, SC = mp.mpf("-0.2820"), mp.mpf("0.2830")
SD = mp.mpf("0.044785178")     # σ_ρd（R27）

def mpf(fr):
    return mp.mpf(fr.numerator)/fr.denominator

print("=== Q-A 結構：規則作用於幾個槽？ ===")
print("  K⁻¹-Gram = (1/90)(A₂ Cartan ⊕ 1) ⟹ 讀出空間 = **2 耦合槽 + 1 單態槽 = 3**")
print("  A₂ Cartan = [[2,−1],[−1,2]] ⟹ 只有 u,d 進入耦合塊；s 為單態（正交）")
print("  ⟹ 規則之「3」唔係巧合，係 rank(A₂)+1；第 4 個味只能二選一：")
print("     ① 進入耦合塊 ⟹ Gram 須升為 A₃ Cartan ⊕ 1（rank 3）⟹ **v11 未註冊，須新結構（dof > 0）**")
print("     ② 另立單態 ⟹ 與耦合塊正交 ⟹ **對 (ρ_u, ρ_d) 零約束**（同 M12 docs/14 §六 介子扇區之資訊論事實）")

print("\n=== Q-B 觀測：重味有無 δ？ ===")
print(f"  輕扇區（已裁決）： δ_u = {mp.nstr(mpf(DU),8):>10s}  δ_d = {mp.nstr(mpf(DD),8):>8s}  δ_s = {mp.nstr(mpf(DS_),8):>8s} MeV")
print(f"  帶粲：            δ_c = {mp.nstr(DC,6)} ± {mp.nstr(SC,6)} MeV")
zc = (DC-0)/SC
print(f"  z(δ_c vs 0) = {mp.nstr(zc,5)} ⟹ δ_c 與 **零相容（1.0σ）**")
print(f"  δ_c / δ_s = {mp.nstr(DC/mpf(DS_),5)}  （若 c 進入耦合塊，比值應為 O(0.1–1)）")
print(f"  |δ_s| / σ(δ_c) = {mp.nstr(mpf(DS_)/SC,5)} ⟹ δ_c 與 δ_s 之量級差 **{mp.nstr(mpf(DS_)/SC,4)}σ**")
print("  帶底：Σb⁻−Σb⁺、Ξb⁻−Ξb⁰ 只以組合 C1 命中（z=+0.27）⟹ **δ_b 未隔離**（1 方程 2 未知）")
print("  頂：t **無同味對**（無 isospin 伴侶）⟹ δ_t 在二體荷語言中**無定義**")

print("\n=== Q-C 判：『適用』之三義 ===")
print("  (1) 同一規則可唔可以作用於另一組 3 槽三元組？")
print("      ✅ 可以 —— 但係**空洞**：docs/05 站無關性定理證明 ρ **只依賴 α**、與站指派無關（偏差 ≤1e−34）")
print("         ⟹ 換一組站／換一組味，只要 α 相同，預言一字不變 ⟹ **唔係新檢定，係恆等**")
print("  (2) 同一規則可否**延伸至第 4 個味**（c、b）？")
print("      🔴 不可，且係結構性：")
print("         · 用 A₂ 原式外推 ⟹ 需把 c 放入耦合塊 ⟹ Gram 須變 A₃ ⟹ 新結構（dof>0），違 T2")
print("         · 留作單態 ⟹ 零約束（同介子扇區）")
print("         · 觀測已否決「c 進入耦合塊」：δ_c = −0.282 ± 0.283 與 0 相容，")
print(f"           而耦合塊成員之 |δ| 為 1.38（d）、3.43（s）⟹ 若 c 等量耦合，差 "
      f"{mp.nstr(mpf(DD)/SC,4)}σ（d）/{mp.nstr(mpf(DS_)/SC,4)}σ（s）")
print("  (3) 同一規則可否套到 d,s,c,t 之**其他組合**（如 (d,s,b) 或 (u,c,t)）？")
print("      🔴 不可：δ_x 之定義係 a(d,x)−a(u,x)（u 與 d 之電荷載體），")
print("         x 只可以係**旁觀者味**；(u,c,t) 之間無 isospin 關係 ⟹ 無 δ 可言。")

print("\n=== 附帶：命名陷阱 ===")
print("  M12 docs/14 之 ρ_s = 1.148613、ρ_c = 17.09 係 **介子/重子比值**（|δ^{qq̄}|/|δ^{qq}|），")
print("  與 M13 之 ρ_x = δ_x/δ_s（重子電荷比）**唔同物件** ⟹ 跨記錄引用必寫明定義，禁用同一字母。")

print("\n=== 判 ===")
print("  · 規則之內容 = 世代權重 α_ℓ ∝ (log₂B)^{ℓ−1}（味盲，合 M7 味盲定理 F̂ ∉ 𝒜）")
print("  · 槽數 = 3 係 A₂ rank+1 之結構結果，**唔可以外推第 4 槽**（否則 Gram 換代 ⟹ dof>0）")
print("  · 重味數據（δ_c ≈ 0，1.0σ）與「c 不入耦合塊」相容 ⟹ **M13 不必亦不可延伸**")
print("  · t 完全在語言之外（無同味對）")
sys.exit(0)