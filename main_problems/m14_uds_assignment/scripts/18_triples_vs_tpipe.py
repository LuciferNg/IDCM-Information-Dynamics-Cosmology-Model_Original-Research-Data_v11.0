#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/18 — 三個共振三元組 vs T_pipe 接收窗（n* = 27）

────────────────────────────────────────────────────────────────────────
已註冊（skill idcm-88site-dual-chain）：
    n* = N_d − L₆ + 1 = 44 − 18 + 1 = 27
    接收端 σ_n = 1（n ≤ n*）／0（n ≥ n*+1）⟹ n ≥ 28 之模接收失效
    讀數 V_n = σ_n · A_n，A_n = 2sin(nπ/90)（/90 Dirichlet 色散 ω_n = 2sin(nπ/90)）
本腳本：把 docs/15 之 3 個共振三元組放進接收窗，逐模計 σ、A_n、V_n。
性質：確定性，N_trials = 0。零觀測量。dps = 40。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp

mp.mp.dps = 40
L = "=" * 78
ND, NSTAR = 44, 27

TRIP = [(9, 15, 18), (9, 27, 30), (15, 27, 36)]


def sigma(n):
    return 1 if n <= NSTAR else 0


def A(n):
    return 2 * mp.sin(n * mp.pi / 90)


print(L)
print("M14 scripts/18 — 共振三元組 vs T_pipe 接收窗")
print(L)
print(f"\n  n* = N_d − L₆ + 1 = {ND} − 18 + 1 = {NSTAR}；σ_n = 1 (n≤{NSTAR}) / 0 (n≥{NSTAR+1})")

print("\n【§一】逐模之接收狀態")
modes = sorted({m for t in TRIP for m in t})
print(f"  {'n':>3s} {'σ':>2s} {'γ_n':>12s} {'A_n':>12s} {'V_n':>12s}  角色")
role = {9: "共用（T1,T2）", 15: "共用（T1,T3）", 27: "共用（T2,T3）",
        18: "合成（T1）", 30: "合成（T2）", 36: "合成（T3）"}
for n in modes:
    print(f"  {n:>3d} {sigma(n):>2d} {mp.nstr(4*mp.sin(n*mp.pi/90)**2, 10):>12s}"
          f" {mp.nstr(A(n), 10):>12s} {mp.nstr(sigma(n)*A(n), 10):>12s}  {role[n]}")

print("\n【§二】逐三元組之存活")
for t in TRIP:
    n_in = [n for n in t if sigma(n)]
    comp = max(t)
    print(f"  {t}: 窗內 {len(n_in)}/3 {n_in}；合成模 {comp} "
          f"{'在窗內 ✅' if sigma(comp) else '在窗外 🔴（σ=0 ⟹ 接收失效）'}")

# ======================================================================
print("\n【§四】★ 共用模之振幅恰為 φ-階梯（精確）")
from fractions import Fraction as F


class Q5:
    def __init__(self, a, b=0):
        self.a, self.b = F(a), F(b)

    def __add__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    def __eq__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return self.a == o.a and self.b == o.b

    def __repr__(self):
        return f"({self.a}{'+' if self.b >= 0 else '−'}{abs(self.b)}√5)"


# A_n = 2 sin(nπ/90) = 2 sin(2n°)
# 2sin18° = (√5−1)/2 = φ⁻¹ ; 2sin30° = 1 = φ⁰ ; 2sin54° = 2cos36° = φ
exact = {9: Q5(F(-1, 2), F(1, 2)),      # (√5−1)/2 = φ⁻¹
         15: Q5(1, 0),                  # 1 = φ⁰
         27: Q5(F(1, 2), F(1, 2))}      # (1+√5)/2 = φ
print("  精確值（a + b√5）：")
names = {9: "φ⁻¹", 15: "φ⁰ = 1", 27: "φ"}
for n in (9, 15, 27):
    num = mp.nstr(A(n), 12)
    print(f"    A_{n:<2d} = 2sin(2·{n}°) = {num:>14s} ＝ {names[n]:<8s} 精確 {exact[n]}")
phi = (1 + mp.sqrt(5)) / 2
dev = max(abs(A(9) - 1 / phi), abs(A(15) - 1), abs(A(27) - phi))
print(f"  數值核：max|A − φ^(−1,0,1)| = {mp.nstr(dev, 3)}")
print("  ⟹ **三個共用模之振幅 = φ⁻¹ : 1 : φ（等比 φ）**，精確，非擬合")
print("\n  合成模之能量（γ）：")
for n in (18, 30, 36):
    print(f"    γ_{n:<2d} = {mp.nstr(4*mp.sin(n*mp.pi/90)**2, 12)}"
          f"   ＝ {['3−φ', '3', '2+φ'][(18,30,36).index(n)]}")
print(f"  γ_18 + γ_36 = {mp.nstr(4*mp.sin(18*mp.pi/90)**2 + 4*mp.sin(36*mp.pi/90)**2, 12)}（= 5，精確）")
print(f"  γ_18 + γ_30 + γ_36 = {mp.nstr(sum(4*mp.sin(n*mp.pi/90)**2 for n in (18,30,36)), 12)}（= 8 = 2B）")

print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ① 三個共振三元組只有 **1 個**（{{9,15,18}}）完整落於接收窗內；
     T2、T3 之合成模（30、36）σ = 0 ⟹ 其「和頻產物」在觀察者端**唔存在**。
  ② 合成模之 V 值：V_18 = {mp.nstr(A(18), 8)}；V_30 = 0；V_36 = 0
     ⟹ 若「三元組 ↔ 世代」，則 3 個世代之中 2 個無質量讀數 —— 與觀測（3 個有質量）不符
     ⟹ **R1／R2 讀法被 T_pipe 結構 disfavour**（非否證：需先有映射規則）。
  ③ 巧合值得記錄：n* = {NSTAR} 本身係 T2、T3 之成員（邊界模 = 27 = 3³）。
  ④ 剩餘讀法：三元組唔係世代；世代須由別處（質量譜 K_chain 之 M45 帶）給出。
""")
print(L)
raise SystemExit(0)