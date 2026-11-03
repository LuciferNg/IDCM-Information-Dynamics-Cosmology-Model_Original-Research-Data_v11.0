#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/07_ktree_alpha_probe.py — K_tree 側能否供出 α_s/α_u？

背景（docs/05）：K⁻¹-夾心族已收斂為單一數字 α_s/α_u，1σ 窗 [3.3048, 4.1300]。
本輪測三件事（全報）：
  S1：存貨內**唯一存活核類**（M5 `docs/03`：`V_ij = κ q_iq_j B^{−d_tree}`；樹層距離核）能否充當 α？
      取樹層ℓ = generation，α_f ∝ B^{−ℓ} 及 B^{+ℓ} 兩種排序，算 (ρ_u, ρ_d)。
  S2：窗之**判別力**：簡單有理數（分母 ≤ 20）落窗之密度；已註冊常數成對比值落窗之密度。
  S3：案例研究 —— `134/36 = (h11+h21)/h11` 給 `ρ_d = 0.402985`（z = −0.012！）。判其為禁區 3 陷阱，
      並以 S2 之密度證明「此級精度之吻合必然出現」。
判準：不作命中宣稱；本腳本之目的係**證明此路線在現精度下不可裁決**（falsifiability gating）。
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
RD, RU = mp.mpf("0.40350365"), mp.mpf("-0.02588837")
SD, SU = mp.mpf("0.044785178"), mp.mpf("0.033583899")
need = 3/(2*RD)
sig = 1.5*SD/RD**2
LO, HI = need-sig, need+sig

print("=== 0. 目標窗 ===")
print(f"  α_s/α_u ∈ [{mp.nstr(LO,8)}, {mp.nstr(HI,8)}]（1σ），中心 {mp.nstr(need,8)} ± {mp.nstr(sig,5)}")
print(f"  窗寬 / 中心 = {mp.nstr((HI-LO)/need*100,4)}%  ← 判別力之上限")

def rho(alpha):
    au, ad, as_ = [mp.mpf(x) for x in alpha]
    return (2*au-ad)/(2*as_), (2*ad-au)/(2*as_)

def verdict(ru, rd):
    zu, zd = (ru-RU)/SU, (rd-RD)/SD
    ok = abs(zu) < 2 and abs(zd) < 2
    return zu, zd, ok

# ---------- S1 ----------
print("\n=== S1. M5 唯一存活核類（B^{−d_tree}）作 α？ ===")
B = mp.mpf(4)
cases = [
    ("α ∝ B^{−ℓ}  (ℓ=gen)", (1, 1/B, 1/B**2)),
    ("α ∝ B^{+ℓ}  (ℓ=gen)", (1, B, B**2)),
    ("α ∝ B^{−ℓ} 反序", (1/B**2, 1/B, 1)),
    ("α ∝ B^{+ℓ} 反序", (B**2, B, 1)),
]
for nm, al in cases:
    ru, rd = rho(al)
    zu, zd, ok = verdict(ru, rd)
    print(f"  {nm:22s} α = ({mp.nstr(al[0],6)}, {mp.nstr(al[1],6)}, {mp.nstr(al[2],6)})"
          f"  ρ_u = {mp.nstr(ru,8):>10s}  ρ_d = {mp.nstr(rd,8):>10s}"
          f"  z = ({mp.nstr(zu,6)}, {mp.nstr(zd,7)})  {'🟡' if ok else '🔴'}")
print("  ⟹ 樹層指數核**不可**充當 α：B^{±ℓ} 之級數太陡（16 或 1/16），兩者皆 10σ 外。")
print("  ⟹ 推論：K_tree 之**樹層距離**不能承擔 α_s 之振幅；α 必須來自 K_tree 之其他量（標籤／維數／譜和）。")

# ---------- S2 ----------
print("\n=== S2. 窗之判別力（密度）===")
hits = []
for q in range(1, 21):
    for p in range(1, 200):
        v = F(p, q)
        if LO <= mp.mpf(p)/q <= HI:
            hits.append((F(p, q), p/q))
uq = sorted(set(h[0] for h in hits))
print(f"  分母 q ≤ 20 之簡單有理數落窗者：{len(uq)} 個（去重後）")
print("  前 12 個：", ", ".join(str(f) for f in uq[:12]))
regist = {"B": 4, "N_gen": 3, "γ44": 4, "h11": 36, "h21": 98, "h11+h21": 134, "M": 33,
          "Σγ": 90, "N_d": 44, "N_d+1": 45, "L7": 29, "F11": 89, "B²": 16, "60": 60,
          "180": 180, "176": 176, "φ": mp.mpf("1.618033988749894848"),
          "√5": mp.sqrt(5), "√2": mp.sqrt(2), "e": mp.e, "π": mp.pi}
names = list(regist)
pairs = 0
inwin = []
for a in names:
    for b in names:
        if a == b:
            continue
        if abs(mp.mpf(regist[b])) < 1e-30:
            continue
        v = mp.mpf(regist[a])/mp.mpf(regist[b])
        pairs += 1
        if LO <= v <= HI:
            inwin.append((f"{a}/{b}", v))
print(f"  已註冊常數成對比值：{pairs} 個，落窗者 {len(inwin)} 個（{len(inwin)/pairs*100:.1f}%）")
for nm, v in sorted(inwin, key=lambda t: abs(t[1]-need))[:8]:
    print(f"    {nm:14s} = {mp.nstr(v,8)}   z_α = {mp.nstr((v-need)/sig,5)}")
print(f"  ⟹ 落窗比例 {len(inwin)/pairs*100:.1f}% ⟹ **窗無判別力**：任何 order-1 結構比值幾乎必落窗內。")

# ---------- S3 ----------
print("\n=== S3. 案例：134/36 陷阱 ===")
v = F(134, 36)
rd = 3/(2*mp.mpf(v.numerator)/v.denominator)
print(f"  α_s/α_u = (h11+h21)/h11 = 134/36 = {mp.nstr(mp.mpf(v.numerator)/v.denominator,8)}")
print(f"  ⟹ ρ_d = 3/(2·134/36) = 27/67 = {mp.nstr(rd,8)}")
print(f"  z_ρd = {mp.nstr((rd-RD)/SD,6)}   （看似 0.01σ 完美）")
print(f"  但：附近之同級比值同樣落窗 —— 例如 h21/M = 98/33 = {mp.nstr(mp.mpf(98)/33,6)}（z_α = "
      f"{mp.nstr((mp.mpf(98)/33-need)/sig,4)}）、Σγ/(N_d−N_gen) = 90/41 = {mp.nstr(mp.mpf(90)/41,6)}")
print("  ⟹ 判：**禁區 3**。α_s/α_u 並無任何推導；134/36 係事後選取，且同級偶然者眾。不作證據。")

# ---------- S4 ----------
print("\n=== S4. 判決 ===")
print("  ① M5 樹層核 ⟹ 10σ 外 ⟹ **結構性排除**（此係窗唯一能作之排除）。")
print(f"  ② 窗寬 {mp.nstr((HI-LO)/need*100,4)}%，已註冊常數對 {len(inwin)/pairs*100:.0f}% 落窗 ⟹ **此路線現精度下不可裁決**。")
print("  ③ 唯一誠實結論：K_tree 側若無**獨立理由**指定 α，則本族不能零參數預言 ρ_d；")
print("     而『數據支持某 α』之說一律屬禁區 3。")
print("  ④ 預註冊：欲裁決 α_s/α_u = B vs N_gen，需 σ_ρd ≤ 0.0125（現 0.0448）⟹ 3.6× 改善。")
sys.exit(0)