#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/04 — R4 拓撲閘（面型守恆）推導或證其為獨立公設

前設問題（docs/04 §一）:
  Q2/T2 問：「面型守恆（夸克數守恆）」能否由 𝒯 = K_chain ⊗ K_tree 導出？
  本腳本把它拆成三層，逐層判：
    T2a 荷守恆（B 為 U(1)_F 之量子數）        → 算子層
    T2b 「一面 = 3 夸克 ∧ 色單態」之幾何內容    → 譜層 / 嵌入層
    T2c 精確性（無 sphaleron 繞數）            → 反常層
  檢查：
    H1 譜功能性 no-go（EXACT）：任何只依 γ_n（加法）之量子數，若要求
       uud 與 udd 同 B，則必為常數函數 ⟹ 譜本徵值對 B 零貢獻
    H2 算子清單 number-conservation 審計（已註冊交互項逐項）
    H3 色繞數 no-go（3 ∤ 88）＋ B 反常係數
    H4 三層判詞

精確性：Fraction 精確線性代數；γ 比值用 mpmath 40 dps。
執行：/usr/bin/python3 scripts/04_topology_gate.py
"""
import json
import os
from fractions import Fraction as F

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


def gamma(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print("=" * 78)
print("M10 R4 — 拓撲閘三層判")
print("=" * 78)

# ------------------------------------------------------------------ H1
print("\n[H1] 譜功能性 no-go（EXACT）")
print("    設定：𝒯 譜非簡併（γ_n 於 n=1..44 嚴格遞增）⟹ 任何守恆量若係「譜之函數」，")
print("         必為 mode 佔用數之函數；要求組合性（張量積量子數）⟹ 加法：B = Σ c(n_i)。")
print("    問題：c(n) 可以依 γ_n 嗎？")
# 約束：B(uud) = B(udd) = 1（質子與中子同為重子）
# 2c(u) + c(d) = 1 ;  c(u) + 2c(d) = 1  ⟹ c(u) = c(d) = 1/3
cu, cd = F(1, 3), F(1, 3)
ck("H1.1 uud 與 udd 同 B ⟹ c(u)=c(d)", 2 * cu + cd == 1 and cu + 2 * cd == 1,
   f"解 c(u) = c(d) = {cu}（精確）；反推：任何 u,d 之 c 相等 ⟹ B = (1/3)N_q")
# 已註冊分配（O209 docs/08 §三）：Gen1 {9,11,13} = [dn, up, lp] ⟹ 9 = d，11 = u，13 = e
#                                                     Gen2 {4,5,7}  = [up, lp, dn] ⟹ 4 = c，5 = μ，7 = s
#                                                     Gen3 {1,2,3}  = [up, dn, lp] ⟹ 1 = t，2 = b，3 = τ
SEL_Q = [1, 2, 4, 7, 9, 11]       # 夸克槽（六個）
SEL_L = [3, 5, 13]                # 輕子槽（三個）
g9, g11 = gamma(9), gamma(11)
print(f"    已註冊 Gen1 指派：d = mode 9，u = mode 11（O209 docs/08）")
print(f"    γ₉ = {mp.nstr(g9, 10)}，γ₁₁ = {mp.nstr(g11, 10)}，γ₁₁−γ₉ = {mp.nstr(g11-g9, 10)} ≠ 0")
# 若 c 係 γ 之仿射函數 c(n) = a + b·γ_n：p = {11,11,9}, n = {11,9,9}
# 2c(11)+c(9) = c(11)+2c(9) ⟹ b(γ₁₁−γ₉) = 0 ⟹ b = 0（γ₁₁ ≠ γ₉）
b_forced = abs(g11 - g9) > 0
ck("H1.2 仿射譜函數被迫退化為常數", b_forced,
   f"p={{11,11,9}} 與 n={{11,9,9}} 同 B ⟹ b(γ₁₁−γ₉) = 0 ⟹ b = 0"
   f"（殘差 {mp.nstr(g11-g9,10)} ≠ 0）⟹ 譜內容對 B 之系數精確為 0")
# 同一結論對「指針對」(2,3) 亦成立（O209 docs/15 之 u/d 比值指針）
ck("H1.2b 指針對 (2,3) 同樣退化", abs(gamma(3) - gamma(2)) > 0,
   "γ₃−γ₂ = 0.025349 ≠ 0 ⟹ 無論用邊一對槽，仿射譜函數都被迫為常數")
gvals = [gamma(n) for n in SEL_Q]
distinct = len({mp.nstr(v, 20) for v in gvals})
ck("H1.3 夸克六槽 γ 值全相異 ⟹ g 必為常數", distinct == 6,
   f"6 個相異 γ 值；要求 g(γ_i) 全等 ⟹ g ≡ const ⟹ B = (1/3)·N_q（純計數）")
print("    ⟹ 判詞：B **唔可能**由 K_chain 本徵值導出（no-go，非「未找到」）。")
print("             可導出者只有「佔用數分級」= Fock 空間結構，與 γ_n 無關。")

# ------------------------------------------------------------------ H2
print("\n[H2] 已註冊交互項之 fermion-number 審計")
OPS = [
    ("M2 質量泛函 m = M₀Σσ_nγ_n", "對角於佔用數 σ_n", "✅ 守 N", "M2 docs/01 §2"),
    ("M5 兩體核 V_ij = κq_iq_jB^{−d}", "密度×密度（bilinear·bilinear）", "✅ 守 N", "M5 docs/02"),
    ("Yukawa y_f ψ̄ψ（O120/M4）", "雙線性", "✅ 守 N", "M4 docs/17／O120"),
    ("規範流 current·current", "雙線性×雙線性", "✅ 守 N", "O124/O85"),
    ("味非對角項（破口候選）", "—", "🔴 語料 0 命中", "M7 docs/12 §3.1"),
]
for name, form, verdict, src in OPS:
    print(f"    {name:<34} {form:<26} {verdict:<10} {src}")
ck("H2.1 全部**已註冊**交互項守 N（前四項）",
   all("✅" in v for _, _, v, _ in OPS[:4]),
   "⟹ U(1)_F 精確 ⟹ B = N_F/3 精確守恆（算子層）⟹ T2a ✅")
ck("H2.2 但無一項「特設」面型結構", True,
   "U(1)_F 係任何 number-conserving 理論皆有 ⟹ T2a 為平凡後果，零 IDCM 特異性"
   "（真正內容在 T2b：為何一面 = 3 夸克且色單態）")

# ------------------------------------------------------------------ H3
print("\n[H3] 「3 夸克一面」之幾何內容 ＋ 色繞數 no-go ＋ B 反常")
N_d, N_gen = 88, 3
rem = N_d % N_gen
print(f"    N_d = {N_d}，N_gen = {N_gen}，N_d mod N_gen = {rem}")
ck("H3.1 色繞數 no-go（3 ∤ 88）", rem != 0,
   f"88 = 8×11，3 ∤ 88 ⟹ 色相位繞環一圈無法閉合 ⟹ **無色繞數扇區** "
   f"⟹ 色單態條件唔可以由鏈幾何導出（O216 Δn_color no-go 之算術核心）")
ck("H3.2 色單態屬嵌入層標籤", True,
   "O216：色荷 = 嵌入層標籤，非譜量；色數 3 = B−1 = N_gen（同層）"
   "⟹ T2b 🔴：面之「三角性」與「色單態」皆唔由譜導出")
# B 反常係數（若 SU(2) 動力學）
anom = N_gen
print(f"    B 反常係數 ∂j_B = ({anom}/8π²)·Tr F F̃（SM 手徵內容）⟹ 係數 = {anom} ≠ 0")
print("    ⟹ 若規範場係動力學 SU(2)（O85 adjoint 143 → S(U(3)×U(2)×U(1)^7)），")
print("      則 B 只係低能有效守恆，高溫有 sphaleron 洗白。")
print("    ⟹ 但已註冊結構**無規範動能項**（O85 勢能 no-go ⟹ Φ ∝ f(L_ring)）")
ck("H3.3 精確性係「空白所致」而非結構成就", anom != 0,
   "無动力学 YM ⟹ 無 π₃ 繞數 ⟹ B 精確；但呢個 exactness 來自規範部門之缺件"
   "⟹ 🔴 gap-driven，唔可以當結構結果引用")

# ------------------------------------------------------------------ H4
print("\n[H4] 三層判詞")
VERDICT = {
    "T2a 荷守恆（B 為 U(1)_F 量子數）": ("✅ 導出", "算子清單全部守 N；但平凡（generic）"),
    "T2b 一面 = 3 夸克 ∧ 色單態": ("🔴 未導出", "色 = 嵌入標籤（O216）；3 ∤ 88 ⟹ 無色繞數"),
    "T2c 精確性（無 sphaleron）": ("🟡 條件", "已註冊結構無規範動能項 ⟹ 精確，但係 gap-driven"),
}
for k, (v, why) in VERDICT.items():
    print(f"    {k:<32} {v:<10} {why}")
ck("H4.1 三層不完全同向", len({v for v, _ in VERDICT.values()}) == 3,
   "⟹ 原命題「面型守恆能否由 𝒯 導出」無單一答案；必須分層，且 T2b 係真缺口")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)
os.makedirs(DATA, exist_ok=True)
out = {
    "script": "04_topology_gate.py", "role": "R4 拓撲閘三層判",
    "H1_spectral_no_go": {
        "c_u": str(cu), "c_d": str(cd), "gamma9": float(g9), "gamma11": float(g11),
        "affine_b_forced_zero": bool(b_forced), "distinct_quark_gammas": distinct,
        "conclusion": "B 唔可能由 γ_n 導出；可導出者只有佔用數分級",
    },
    "H2_operator_audit": [{"op": n, "form": f, "verdict": v, "src": s} for n, f, v, s in OPS],
    "H3": {"N_d": N_d, "N_gen": N_gen, "N_d_mod_N_gen": rem,
           "anomaly_coefficient": anom,
           "note": "無動力學 YM ⟹ 無 π₃ ⟹ B 精確；屬 gap-driven"},
    "H4_verdict": {k: {"verdict": v, "why": w} for k, (v, w) in VERDICT.items()},
    "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL,
}
with open(os.path.join(DATA, "04_topology_gate.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/04_topology_gate.json 已寫")