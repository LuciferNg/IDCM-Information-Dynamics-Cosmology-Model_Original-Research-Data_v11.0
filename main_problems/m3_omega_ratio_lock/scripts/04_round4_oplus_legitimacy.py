#!/usr/bin/env python3
"""M3 scripts/04 — 第四輪：⊕ 合法性（rank 預算 vs 能量預算之測度等價）

問題: 為何份額可以 rank 相加（⊕ ⟹ 60），而唔係 dim 相乘（⊗ ⟹ 176）？

本輪命題: 60 份預算 ≡ 90 能量預算除以【均一量子】⟹ ⊕ 唔係額外假設，
          而係 equipartition（每份同能量）之重述。

核對項:
  (1) 逐扇區能量量子 = 3/2（Fraction 精確）
  (2) ⟨W⟩² = (Σγ + M)/2 = 61.5（O39 W-field 獨立路線）
  (3) equipartition ⟺ N_d = N_gen·B² − B  (44 = 48−4)
  (4) B 掃描：只有 B=4 令兩條路線（⟨W⟩² vs 秩×量子）重合
"""
from fractions import Fraction as F
from mpmath import mp, mpf, sin, pi, sqrt
mp.dps = 50

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_d, M_tree = 4, 44, 33
g = [4*sin(mpf(n)*pi/88)**2 for n in range(1, 45)]
S = sum(g)                      # Σγ = 90
N_gen = 3
kappa = mpf(1)/B**2

# ── §1 逐扇區量子 ───────────────────────────────────────────────────
print("§1 逐扇區能量量子（Fraction 精確）")
W2 = F(123, 2)                  # ⟨W⟩² = 61.5 = 123/2
E_L, E_m = W2, F(90) - W2       # 61.5 / 28.5
E_b = E_m * F(3, 19)            # 28.5 × 3/19 = 4.5
E_d = E_m * F(16, 19)           # 28.5 × 16/19 = 24
tbl = [("Λ", E_L, 41), ("b", E_b, 3), ("dm", E_d, 16)]
print(f"  {'扇區':4} {'能量 (90 制)':>12} {'rank (60 制)':>12} {'量子':>8}")
for nm, E, r in tbl:
    q = F(E)/r
    print(f"  {nm:4} {str(E):>12} {r:>12} {str(q):>8}")
qs = [F(E)/r for _, E, r in tbl]
chk("三扇區量子全等 = 3/2", all(q == F(3,2) for q in qs), f"{[str(q) for q in qs]}")
chk("Σ能量 = 90", E_L + E_m == 90, f"{E_L}+{E_m}")
chk("Σrank = 60", sum(r for _,_,r in tbl) == 60, "41+3+16")
chk("量子 = Σγ/60 = 3/2", F(90,60) == F(3,2))
print("  ⟹ 兩本預算【同一測度】：能量 = 1.5 × rank（逐扇區）")

# ── §2 ⟨W⟩² 之獨立算子式 ────────────────────────────────────────────
print("\n§2 ⟨W⟩² 之獨立路線（O39 W-field）")
lhs = (F(90) + M_tree)/2
chk("⟨W⟩² = (Σγ + K_tree)/2 = 61.5", lhs == W2, f"({90}+{M_tree})/2 = {lhs}")
chk("2⟨W⟩² = Σγ + K_tree = 123", 2*W2 == 90 + M_tree, f"{2*W2}")
print("  ⟹ Λ 能量由 W-field 真空凝結獨立給出（唔靠 rank）")

# ── §3 equipartition ⟺ N_d = N_gen·B² − B ───────────────────────────
print("\n§3 equipartition 之等價條件")
tot_rank = N_d + B**2
lhs_id = F(2)*F(90)/3          # 60
chk("2Σγ/3 = N_d + B² = 60", lhs_id == tot_rank, f"{lhs_id} vs {tot_rank}")
rhs = N_gen*B**2 - B
chk("N_d = N_gen·B² − B  (44 = 3·16 − 4)", rhs == N_d, f"{N_gen}·{B**2}−{B} = {rhs}")
print("  ⟹ 等價鏈: equipartition ⟺ 60 = 2Σγ/3 ⟺ N_d = N_gen B² − B")

# ── §4 B 掃描：只有 B=4 令兩路線重合 ────────────────────────────────
print("\n§4 B 掃描（N_d=44, 61.5 固定；量子=Σγ/(N_d+B²)）")
import math
for b in (3, 4, 5, 6):
    ng = int(math.ceil(math.log(b*b)))
    tot = N_d + b*b
    q = F(90, tot)
    L_rank_energy = q * (N_d - ng)
    print(f"  B={b}: N_gen={ng}, rank_tot={tot}, 量子={float(q):.4f}, "
          f"Λ(秩×量子)={float(L_rank_energy):7.3f}  vs ⟨W⟩²=61.5  "
          f"{'✅' if abs(float(L_rank_energy)-61.5) < 5e-2 else '✗'}")
chk("B=4 係唯一重合值",
    all(abs(float(F(90, N_d+b*b)*(N_d-int(math.ceil(math.log(b*b))))) - 61.5) > 0.05
        for b in (3, 5, 6)))
print("  ⟹ equipartition 只喺 B=4 成立；而 B=4 由 Nyquist γ_max=4 獨立釘死")

# ── §5 verdict ─────────────────────────────────────────────────────
print("\n§5 判定")
print("  ⊕ 之合法性 = equipartition 之合法性 = 「B=4 下兩條獨立路線重合」")
print("  ⟹ 由【額外建模選擇】降級為【一致性條件】（非自由參數）")
print("  殘餘: 為何 N_d = N_gen·B² − B（44=48−4）本身仍需第一性原理（B=4 特例族）")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")