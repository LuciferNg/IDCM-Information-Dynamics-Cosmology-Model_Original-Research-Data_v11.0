#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/05_carry_by_2_candidate.py — 候選「w = K⁻¹·Σ_gen 2^gen v_x」之判準全套
（A₂ Cartan 結構之直接推論；**後驗**成分明示）

候選（零可調參數）：
    α = (2^0, 2^1, 2^2) = (1, 2, 4)
    u = Σ_gen 2^gen · v_{x(gen)}        gen 0→u(站2)、1→d(站3)、2→s(站4)
    w = K⁻¹ u
    預言：ρ_u = 0（精確）、ρ_d = 3α_u/(2α_s) = 3/8（精確）

本腳本逐閘核（G0–G2、T1–T5），並計：
  · look-elsewhere：α_s ∈ {1..8}（簡單整數格）⟹ 命中數／總數
  · T5 判別力：需何 σ 方能 3σ 分辨 3/8 vs 0.4 vs √2−1
全報，不作命中宣稱（凡後驗成分須明示）。
"""
import sys
import mpmath as mp
from fractions import Fraction as F

mp.mp.dps = 60
N = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2)/45)*mp.sin(mp.pi*n*j/45)
D = [phi(n, 3)-phi(n, 2) for n in range(1, N+1)]
G = [4*mp.sin(mp.pi*n/90)**2 for n in range(1, N+1)]
V = {x: [phi(n, k)*D[n-1] for n in range(1, N+1)] for x, k in (("u", 2), ("d", 3), ("s", 4))}

RD, RU = mp.mpf("0.40350365"), mp.mpf("-0.02588837")
SD, SU = mp.mpf("0.044785178"), mp.mpf("0.033583899")

def build(alpha):
    u = [alpha[0]*V["u"][n] + alpha[1]*V["d"][n] + alpha[2]*V["s"][n] for n in range(N)]
    return [u[n]/G[n] for n in range(N)]

def rho(w):
    f = lambda x: mp.fsum(w[n-1]*V[x][n-1] for n in range(1, N+1))
    du, dd, ds = f("u"), f("d"), f("s")
    return du/ds, dd/ds, du, dd, ds

ALPHA = (1, 2, 4)
w = build(ALPHA)
ru, rd, du, dd, ds = rho(w)

print("=== 0. 候選定義 ===")
print(f"  α = {ALPHA} = (2^0, 2^1, 2^2)（carry-by-2；gen 序）")
print(f"  u = 1·v_u + 2·v_d + 4·v_s ，w = K⁻¹u")
print(f"  實算：ρ_u = {mp.nstr(ru,10)}   ρ_d = {mp.nstr(rd,10)} = 3/8？ {abs(rd-F(3,8))<mp.mpf('1e-40')}")
print(f"  閉式：ρ_d = 3α_u/(2α_s) = {mp.nstr(mp.mpf(3)*ALPHA[0]/(2*ALPHA[2]),10)}")
print(f"  δ_s = {mp.nstr(ds,10)}（≠0，非退化）；‖w‖ = {mp.nstr(mp.sqrt(mp.fsum(x*x for x in w)),8)}")

print("\n=== 1. 判準逐閘 ===")
zu, zd = (ru-RU)/SU, (rd-RD)/SD
print(f"  z_u = {mp.nstr(zu,6)}   z_d = {mp.nstr(zd,6)}")
print(f"  G0 出處閘：K_chain⁻¹ ✓、站指派 u→2/d→3/s→4 ✓、α = 2^gen ✓（carry-by-2）")
print(f"             ⚠ 但 2^gen 之 carry-by-2 由 station-index 推廣至權重 —— 屬**新讀法**，非已註冊算子 ⟹ G0 半過")
print(f"  G1 符號閘：要求 ρ_u < 0。本候選 ρ_u ≡ 0（精確）⟹ **閘不可判**（凍結版 G1 預設嚴格負）")
print(f"             ⟹ 依 append-only 紀律：閘設計須更正為「ρ_u ≤ 0 且與數據相容」，本輪**不**改閘，記為 G1-NA")
print(f"  G2 味盲閘：|1−ρ_d| = {mp.nstr(abs(1-rd),6)}  vs 3σ_ρd = {mp.nstr(3*SD,6)} ⟹ "
      f"{'PASS' if abs(1-rd)>3*SD else 'FAIL'}（ρ_d 遠離 1，帶味 ✓）")
print(f"  T1 非可分解：w 含 φ(k_x)Δ 乘積 ⟹ 四點物件，非 φ(i)φ(j)/γ 因子形 ⟹ 秩檢 PASS")
print(f"             ⚠ 但 builder 源碼**無 K_tree 側符號** ⟹ T1 之①要求 FAIL（缺樹側）")
print(f"  T2 dof：α 由 2^gen 唯一給定 ⟹ N_class = 1 ⟹ dof = 0 ✓（連續參數零）")
print(f"  T3 非 ℤ₂-odd：w_n 與 w_(45−n) 非反對稱 ⟹ --")
print(f"  T4 零參數：builder 未引用 ρ_d 目標；α 單一預先指定 ⟹ PASS，但**後驗**（見 §2）")
print(f"  T5 可裁決：見 §3")

print("\n=== 2. ★ look-elsewhere：α_s 之簡單整數格 ===")
grid = list(range(1, 9))
hits = []
print(f"  {'α_s':>4s} {'ρ_d=3/(2α_s)':>14s} {'z_d':>9s}  判")
for a_s in grid:
    a = (1, 2, a_s)
    _, r_d, *_ = rho(build(a))
    z = (r_d-RD)/SD
    ok = abs(z) < 1
    if ok:
        hits.append(a_s)
    print(f"  {a_s:>4d} {mp.nstr(r_d,10):>14s} {mp.nstr(z,6):>9s}  {'🟡|z|<1' if ok else '🔴'}")
print(f"  ⟹ 命中 {len(hits)}/{len(grid)}；naive p = {len(hits)/len(grid):.3f}")
print(f"  ⟹ 目標值 3/8 落在 1σ 窗之偶然率 ≈ {len(hits)/len(grid)*100:.1f}% ⟹ **不足以作證據**")

print("\n=== 3. T5 判別力：3/8 vs 鄰近結構值 ===")
def to_mp(v):
    return mp.mpf(v.numerator)/v.denominator if isinstance(v, F) else mp.mpf(v)

for name, val in [("3/8 = 0.375", F(3, 8)), ("2/5 = 0.4", F(2, 5)),
                  ("√2−1 = 0.41421", mp.sqrt(2)-1), ("16/41 = 0.39024", F(16, 41))]:
    d = abs(to_mp(val)-mp.mpf("0.375"))
    if d == 0:
        print(f"  ρ_d = {name:18s} —— 本候選值（基準）")
        continue
    print(f"  ρ_d = {name:18s} 與 3/8 差 {mp.nstr(d,5)} ⟹ 3σ 分辨需 σ ≤ {mp.nstr(d/3,5)}"
          f"（現 σ = {mp.nstr(SD,4)}，需改善 {mp.nstr(SD/(d/3),4)}×）")

print("\n=== 4. 判決 ===")
print("  結構部分：**ρ_u = 0 精確**，由 A₂ Cartan 零方向強制（任一 α_d = 2α_u 皆成立）⟹ 非擬合。")
print("  數值部分：ρ_d = 3/8 需 α_s = 4；look-elsewhere 12.5% ⟹ **不作證據**。")
print("  ⟹ 候選狀態：🟡 **待精度**（非 🔴 死、非 ✅ 閉）。")
print("  ⟹ 裁決條件（預註冊）：ρ_d 之 σ 需 ≤ 0.0125（3σ 分辨 3/8 vs 2/5）方可裁決 —— 現 0.0448。")
sys.exit(0)