"""M5 docs/22 配套 — D4 範圍審計：符號定理之**前提依賴**（A3 層指派）
================================================================================================
`docs/21` 之幾何符號定理以 A3 層指派為前提：`ℓ(u) = ℓ(c) = 0`、`ℓ(d) = ℓ(s) = 7`。
定理之核心不等式為：幾何抑制比 `S(K)/S(D) = B^{−(d_K − d_D)}`（`d ≡ d_pair`），
    ⟹ 比 < 1 **若且僅若** `d_K > d_D`。
A3 給 `d_K = 0`、`d_D = 7` ⟹ 比 = `B^7` ≫ 1（符號相反）。
但若 `ℓ(d) = 0` 而 `ℓ(s) = 7`，則 `d_K = 7`、`d_D = 0` ⟹ 比 = `B^{−7} = 1/16384`（符號**正確**）。

本腳本做四件事（判準寫死）：
  A1 枚舉 4 個與「up-type 在根」相容之指派（`ℓ(u) = ℓ(c) = 0`；`ℓ(d), ℓ(s) ∈ {0, 7}`），
     逐個計 `d_K`、`d_D`、幾何比 `S(K)/S(D)`（B^{−(d_K−d_D)}）與 `G1` 同點項。
  A2 逐個計「幾何 × 動力學」總比（動力學取 `docs/21` 之最佳點 κ = 0.2150：`P₀(K)/P₀(D) = 0.09173`），
     與數據 0.038298 比，報所需 `d_K − d_D`（非整數則記「無整數指派可中」）。
  A3 範圍陳述：定理只在 A3 指派下成立；不得外推為「態依賴不可能出自幾何」之普遍命題。
  A4 禁：不得藉此改 A3（A3 之 15 成員表與層釘已用 A3；改指派須另立 O）。

輸出：`data/24_assignment_scope.json`
"""
import json
import math
from fractions import Fraction as F

B, N = 4, 88
C_K, C_D = F(135, 10000), F(3525, 10000)
MCS, MS = F(127, 100), F(495, 1000)
DATA_RATIO = C_K / C_D                       # 0.038298
TARGET_TOTAL = DATA_RATIO / (MCS / MS)       # F_7/F_0 = 1/66.99
DYN = 0.09173                                # docs/21：κ = 0.2150 之 P₀(K)/P₀(D)（最佳）

print("=" * 104)
print("M5 docs/22 — D4 範圍審計：幾何符號定理之前提依賴（A3 層指派）")
print("=" * 104)
print(f"\n[D0] 數據 contact 比 = {DATA_RATIO} = {float(DATA_RATIO):.6f}")
print(f"     目標總態因子 F_7/F_0 = {TARGET_TOTAL} = 1/{float(1/TARGET_TOTAL):.4f}")
print(f"     動力學（最佳 κ = 0.2150）P₀(K)/P₀(D) = {DYN} = 1/{1/DYN:.2f}")
print(f" ⟹ 幾何需提供 {float(TARGET_TOTAL)/DYN:.6f} = 1/{DYN/float(TARGET_TOTAL):.3f}")
need_exp = math.log(DYN / float(TARGET_TOTAL), B)
print(f" ⟹ 所需 `d_K − d_D` = log_B(幾何需求之倒數) = {need_exp:.4f}（非整數）")

print(f"\n[A1–A2] 四種相容指派（ℓ(u) = ℓ(c) = 0）")
print(f"     {'ℓ(d)':>5s}{'ℓ(s)':>5s}{'d_K':>5s}{'d_D':>5s}{'幾何比':>14s}{'幾何×動力':>13s}{'vs 數據':>10s}{'判':>26s}")
rows = []
for ld in (7, 0):
    for ls in (7, 0):
        dK, dD = abs(ls - ld), abs(0 - ld)
        geo = F(B) ** (dD - dK)
        tot = float(geo) * DYN
        dev = tot / float(DATA_RATIO)
        if ld == 7 and ls == 7:
            tag = "A3（現行）⟹ 符號相反"
        elif abs(math.log(dev)) <= math.log(2):
            tag = "接近（2× 內）但非 A3"
        elif dev < 1:
            tag = "過強"
        else:
            tag = "過弱"
        rows.append({"ell_d": ld, "ell_s": ls, "d_K": dK, "d_D": dD, "geo": str(geo),
                     "geo_dyn": tot, "dev_vs_data": dev, "tag": tag})
        print(f"     {ld:5d}{ls:5d}{dK:5d}{dD:5d}{str(geo):>14s}{tot:13.6f}{dev:10.3f}{tag:>26s}")

print(f"\n[A3] 範圍陳述")
print(f"     幾何比 < 1 ⟺ d_K > d_D ⟹ 只有當『s 較 d 更深層』時，幾何才可能給出正確符號。")
print(f"     A3 給 d_K = 0 < d_D = 7 ⟹ 幾何比 = B^7 = {B**7}（`docs/21` 之定理）。")
print(f"     替代指派（ℓ(d) = 0、ℓ(s) = 7）⟹ 幾何比 = B^-7 = {B**-7:.6e}，連同動力學後 = {float(F(1,B**7))*DYN:.6e}"
      f"（數據 {float(DATA_RATIO):.6f} ⟹ 過強 {float(DATA_RATIO)/(float(F(1,B**7))*DYN):.1f}×）。")
print(f"     所需 d_K − d_D = {need_exp:.4f} 非整數 ⟹ **無任何整數指派可精確命中**（最近者 d = 1："
      f"B^-1×動力學 = {B**-1*DYN:.6f} = 數據之 {B**-1*DYN/float(DATA_RATIO)*100:.1f}%"
      f"（低 {100*(1-B**-1*DYN/float(DATA_RATIO)):.1f}% ⟹ 2× 內））。")
print(f"     ⟹ `docs/21` 之定理**不得外推**為「態依賴不可能出自幾何」之普遍命題；")
print(f"        正確表述：**在 A3 指派下**幾何符號相反；改指派須另立 O（A3 之 15 成員表與層釘已用 A3）。")

out = {"data_ratio": str(DATA_RATIO), "target_total": str(TARGET_TOTAL), "dyn": DYN,
       "needed_dK_minus_dD": need_exp, "rows": rows,
       "scope": "幾何符號定理僅在 A3 指派（ℓ(d)=ℓ(s)=7）下成立；替代指派可令幾何符號正確但量級過強",
       "nearest_integer": {"d": 1, "total": B ** -1 * DYN, "frac_of_data": B ** -1 * DYN / float(DATA_RATIO)},
       "realizable_example": {"ell_d": 3, "ell_s": 7, "d_K": 4, "d_D": 3, "note": "與 A3 相斥，須 A3/P10 側獨立推導"}}
with open("data/24_assignment_scope.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print(f"     ★ 可實現性：`d_K − d_D = 1` 可由 e.g. ℓ(d) = 3、ℓ(s) = 7 達成（d_K = 4、d_D = 3）；"
      f"但與 A3（ℓ(d) = ℓ(s) = 7）相斥 ⟹ 屬**指派層**問題，非幾何層。")
print(f"     ⟹ 分工：M5 只 register（殘餘 1/6.14 之來源 = A3 之 ℓ(d) = ℓ(s) 假設）；"
      f"改指派須 A3／P10 側獨立推導（A3 之 15 成員表與層釘全部建基於現指派）。")

print("\n-> data/24_assignment_scope.json")