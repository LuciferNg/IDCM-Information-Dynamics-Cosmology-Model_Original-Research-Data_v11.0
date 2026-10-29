#!/usr/bin/env python3
"""R4-P1 — 荷階判別（一體 vs 二體），零參數，無擬合。

預註冊判準（2026-10-29 跑前凍結，不得事後改）
--------------------------------------------------
D1 一體判別：若 Δm = Σ_i a(f_i)（每夸克荷），則恆等式
       r(Σ⁻−Σ⁺) = 2 · r(Ξ⁻−Ξ⁰)
   判死門檻：|z| ≥ 3

D2 二體判別：若 Δm = Σ_(價對) a(type_ij)（配對荷），則恆等式
       r(Σ⁻−Σ⁺) = r(n−p) + r(Ξ⁻−Ξ⁰)
   判死門檻：|z| ≥ 3

D3 樣本外（decuplet J=3/2⁺，未入 M12 任何擬合）：同 D1/D2 形式
       Σ*⁻=dds  Σ*⁰=uds  Ξ*⁻=dss  Ξ*⁰=uss
   解析度閘（跑前寫死）：若 σ(殘差) > 0.5 MeV ⟹ 標「解析度不足」，不作判據。

D4 介子（單一價對，無上面之恆等式）：只列值，供 R4-P2 交叉核。

紀律：禁止調參、禁止掃描、禁止改質量表、禁止逐通道揀方案。
      全部用 Fraction/mpmath 精確算術；浮點只用於最終顯示。
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30

# 質量表：與 M12 scripts/11 / l2 完全一致（協議：不新增/不改既有輸入）
MASS = {
    "n":   ("939.565421", "0.0000001"),
    "p":   ("938.272089", "0.0000001"),
    "S-":  ("1197.449",  "0.030"), "S0": ("1192.642", "0.024"), "S+": ("1189.37", "0.07"),
    "X-":  ("1321.71",   "0.07"),  "X0": ("1314.86",  "0.20"),
    "K+":  ("493.677",   "0.016"), "K0": ("497.611",  "0.013"),
    "D+":  ("1869.66",   "0.50"),  "D0": ("1864.84",  "0.50"),
    # --- 本輪新增（樣本外，PDG 2024 面值）---
    "S*-": ("1387.2",    "0.5"),   "S*0": ("1383.7",  "1.0"),  "S*+": ("1382.8", "0.4"),
    "X*-": ("1535.0",    "0.6"),   "X*0": ("1531.80", "0.32"),
    "K*+": ("891.67",    "0.26"),  "K*0": ("895.55",  "0.20"),
    "D*+": ("2010.26",   "0.05"),  "D*0": ("2006.85", "0.05"),
}


def dm(a, b):
    """Δm = m_a − m_b，及 σ（獨立誤差）"""
    ma, sa = mp.mpf(MASS[a][0]), mp.mpf(MASS[a][1])
    mb, sb = mp.mpf(MASS[b][0]), mp.mpf(MASS[b][1])
    return ma - mb, mp.sqrt(sa * sa + sb * sb)


def both(a, b, label):
    """回 (Δm, σ) 並印出"""
    d, s = dm(a, b)
    print("  %-24s Δm = %+12.5f ± %8.5f MeV" % (label, d, s))
    return d, s


print("=" * 78)
print("R4-P1 荷階判別 — 一體 vs 二體（零參數，無擬合）")
print("=" * 78)

print("\n[0] 錨通道（既有 M12 資料）")
d_np, s_np = both("n", "p", "n−p")
d_S, s_S = both("S-", "S0", "Σ⁻−Σ⁰")
d_X, s_X = both("X-", "X0", "Ξ⁻−Ξ⁰")
d_SSp, s_SSp = both("S-", "S+", "Σ⁻−Σ⁺")
d_K, s_K = both("K+", "K0", "K⁺−K⁰")
d_D, s_D = both("D+", "D0", "D⁺−D⁰")

print("\n[1] 樣本外通道（本輪新增，PDG 2024）")
d_S2, s_S2 = both("S*-", "S*0", "Σ*⁻−Σ*⁰")
d_X2, s_X2 = both("X*-", "X*0", "Ξ*⁻−Ξ*⁰")
d_S2Sp, s_S2Sp = both("S*-", "S*+", "Σ*⁻−Σ*⁺")
d_K2, s_K2 = both("K*+", "K*0", "K*⁺−K*⁰")
d_D2, s_D2 = both("D*+", "D*0", "D*⁺−D*⁰")

print("\n[2] 分割比 r ≡ Δm / Δm(n−p)")
rows = [
    ("r(n−p)", 1.0, 0.0),
    ("r(Σ⁻−Σ⁰)", d_S / d_np, s_S / d_np),
    ("r(Ξ⁻−Ξ⁰)", d_X / d_np, s_X / d_np),
    ("r(Σ⁻−Σ⁺)", d_SSp / d_np, s_SSp / d_np),
    ("r(K⁺−K⁰)", d_K / d_np, s_K / d_np),
    ("r(D⁺−D⁰)", d_D / d_np, s_D / d_np),
    ("r(Σ*⁻−Σ*⁰)", d_S2 / d_np, s_S2 / d_np),
    ("r(Ξ*⁻−Ξ*⁰)", d_X2 / d_np, s_X2 / d_np),
    ("r(Σ*⁻−Σ*⁺)", d_S2Sp / d_np, s_S2Sp / d_np),
    ("r(K*⁺−K*⁰)", d_K2 / d_np, s_K2 / d_np),
    ("r(D*⁺−D*⁰)", d_D2 / d_np, s_D2 / d_np),
]
for name, v, s in rows:
    print("  %-14s = %+10.6f ± %.6f" % (name, v, s))

# ---------------------------------------------------------------- D1 / D2
print("\n" + "=" * 78)
print("[D1] 一體判別：r(Σ⁻−Σ⁺) =? 2·r(Ξ⁻−Ξ⁰)")
pred = 2 * d_X
res = d_SSp - pred
sig = mp.sqrt(s_SSp ** 2 + (2 * s_X) ** 2)
z1 = res / sig
print("     預言 = %.5f MeV   實測 = %.5f MeV   殘差 = %+.5f ± %.5f" % (pred, d_SSp, res, sig))
print("     z = %+.2f   ⟹ %s" % (z1, "判死 🔴" if abs(z1) >= 3 else "未死 🟡"))

print("\n[D2] 二體判別：r(Σ⁻−Σ⁺) =? r(n−p) + r(Ξ⁻−Ξ⁰)")
pred = d_np + d_X
res = d_SSp - pred
sig = mp.sqrt(s_SSp ** 2 + s_X ** 2)
z2 = res / sig
print("     預言 = %.5f MeV   實測 = %.5f MeV   殘差 = %+.5f ± %.5f" % (pred, d_SSp, res, sig))
print("     z = %+.2f   ⟹ %s" % (z2, "判死 🔴" if abs(z2) >= 3 else "通過 🟢"))

# ---------------------------------------------------------------- D3 樣本外
print("\n" + "=" * 78)
print("[D3] 樣本外（decuplet J=3/2⁺）— 同兩式")
pred = 2 * d_X2
res = d_S2 - pred
sig = mp.sqrt(s_S2 ** 2 + (2 * s_X2) ** 2)
z3a = res / sig
print("     一體：預言 = %.5f   實測 = %.5f   殘差 = %+.5f ± %.5f   z = %+.2f"
      % (pred, d_S2, res, sig, z3a))
pred = d_np + d_X2
res = d_S2 - pred
sig = mp.sqrt(s_S2 ** 2 + s_X2 ** 2)
z3b = res / sig
print("     二體：預言 = %.5f   實測 = %.5f   殘差 = %+.5f ± %.5f   z = %+.2f"
      % (pred, d_S2, res, sig, z3b))
print("     解析度閘（σ > 0.5 MeV ⟹ 不作判據）：一體 σ=%.3f %s；二體 σ=%.3f %s"
      % (mp.sqrt(s_S2 ** 2 + (2 * s_X2) ** 2), "不足" if mp.sqrt(s_S2**2 + (2*s_X2)**2) > 0.5 else "足",
         sig, "不足" if sig > 0.5 else "足"))

# ---------------------------------------------------------------- D4 介子
print("\n" + "=" * 78)
print("[D4] 介子單對（無恆等式）— 只列值供 P2")
print("     r(K⁺−K⁰) = %+10.6f ± %.6f" % (d_K / d_np, s_K / d_np))
print("     r(K*⁺−K*⁰)= %+10.6f ± %.6f" % (d_K2 / d_np, s_K2 / d_np))
print("     r(D⁺−D⁰) = %+10.6f ± %.6f" % (d_D / d_np, s_D / d_np))
print("     r(D*⁺−D*⁰)= %+10.6f ± %.6f" % (d_D2 / d_np, s_D2 / d_np))

# ---------------------------------------------------------------- 判詞
print("\n" + "=" * 78)
print("[判詞]")
print("  D1 一體（每夸克荷）: z = %+.2f  ⟹ %s" % (z1, "排除 🔴" if abs(z1) >= 3 else "未排除"))
print("  D2 二體（配對荷）  : z = %+.2f  ⟹ %s" % (z2, "排除 🔴" if abs(z2) >= 3 else "相容 🟢"))
print("  D3 樣本外          : z = %+.2f（二體）  ⟹ %s"
      % (z3b, "解析度不足，不作判據 🟡" if sig > 0.5 else ("判死 🔴" if abs(z3b) >= 3 else "通過 🟢")))
print("  註：D2 為 ansatz 之必要條件（非充分），通過 ≠ 已識別形式。")
print("=" * 78)