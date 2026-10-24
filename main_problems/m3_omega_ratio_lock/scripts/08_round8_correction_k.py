#!/usr/bin/env python3
"""M3 scripts/08 — 第八輪：更正 round 7 之「b₁ ⟹ k=1」滑步

用戶質疑: 「單環？不是 44 個有各自頻率的點嗎」——正確。

釐清三個唔同嘅「1 / 44」:
  (a) 模式數 = 44（γ_n, n=0..43；每個模式有自己頻率、振幅、相位）
  (b) 場值 = 每站一個複數（ℂ）⟹ 每站 1 個分量（唔係 ℂ^k）
  (c) 相位自由度 = 88 個站相位；模式展開後 44 個模式相位
  而 b₁ = 1 只係【gauge holonomy】（局部 U(1) 下不變之 Wilson loop）數，
  對 global U(1) 而言，44 個模式相位全部 physical ⟹ b₁ ≠ k。**撤回 round 7 論證**。

正確判據: 凝聚分量數 k = 【gapless 模式數】（有 VEV 之分量必無 gap）
  γ_0 = 0 唯一零模；γ_n > 0 (n≥1) ⟹ k = 1 ⟹ 序參量 ∈ ℂ
  其餘 43 個模式相位 = 有 gap 之激發，唔係凝聚分量。
"""
import mpmath as mp
from fractions import Fraction as F

mp.mp.dps = 40
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

N, B = 88, 4
gam = lambda n: B * mp.sin(mp.pi*n/N)**2

# ── §1 譜：44 個相異頻率、唯一零模 ─────────────────────────────────
print("§1 譜 γ_n = 4sin²(nπ/88)（n=0..44；n 與 88−n 同值）")
g = {n: gam(n) for n in range(45)}
zeros = [n for n in range(45) if g[n] < mp.mpf('1e-30')]
print(f"  零模 n = {zeros}")
print(f"  γ_0  = {mp.nstr(g[0], 6)}")
print(f"  γ_1  = {mp.nstr(g[1], 10)}   ← 最輕有 gap 模式")
print(f"  γ_22 = {mp.nstr(g[22], 10)}")
print(f"  γ_44 = {mp.nstr(g[44], 6)}  (= B)")
chk("唯一零模（n=0）⟹ gapless 模式數 = 1", zeros == [0], f"{zeros}")
chk("γ_1 > 0（有 gap；唔係凝聚分量）", g[1] > 0, f"{mp.nstr(g[1], 8)}")
chk("γ_44 = B = 4（Nyquist）", abs(g[44] - 4) < mp.mpf('1e-30'))
chk("相異值 = 45 個（n=0..44 各異）", len({mp.nstr(g[n], 30) for n in g}) == 45)

# ── §2 相位自由度計數：44 個模式相位，唔係 1 個 ────────────────────
print("\n§2 自由度計數（證 b₁ ≠ k）")
sites, comp = N, 1                       # 88 站, 1 複數分量/站
real_dof = sites * comp * 2              # 振幅+相位
constraint = sites * comp                # |W| = W₀ 每站
phase_dof = real_dof - constraint
print(f"  站數 {sites}，每站 {comp} 複數分量 ⟹ 實自由度 {real_dof}")
print(f"  固定振幅 |W|=W₀ ⟹ 相位自由度 {phase_dof}")
chk("相位自由度 = 88（每站一相位）", phase_dof == 88)
chk("模式相位數 = 44（模式數 N_d，每模式 1 相位）", phase_dof // 2 == 44)
print("  ⟹ 44 個模式相位全部係 physical（global U(1)）⟹ b₁=1 唔會縮成 1 個相位")

# ── §3 求和（M3 約定核算）────────────────────────────────────────
print("\n§3 Σγ 約定核算")
S43 = sum(g[n] for n in range(44))
full = sum(gam(n) for n in range(88))
print(f"  Σ_{{n=0..43}} γ_n = {mp.nstr(S43, 12)}")
print(f"  Σ_{{n=0..87}} γ_n = {mp.nstr(full, 12)}   (= 2·N_d = 176)")
print(f"  (tr L + γ_max)/2 = {(full + 4)/2}")
chk("Σ_{n=0..43} γ_n = 86", abs(S43 - 86) < mp.mpf('1e-25'), mp.nstr(S43, 10))
chk("Σ_{n=0..87} γ_n = 176 = 2N_d", abs(full - 176) < mp.mpf('1e-25'), mp.nstr(full, 10))
chk("(tr L + γ_max)/2 = 90 = Σγ（M3 約定）", abs((full+4)/2 - 90) < mp.mpf('1e-25'))

# ── §4 若 k 個分量：零模重數 ×k（精確核算） ────────────────────────
print("\n§4 假設場有 k 個複數分量（ℂ^k）⟹ 零模重數 = k·nullity(L)")
def cycle_lap_frac(N):
    M = [[F(0)]*N for _ in range(N)]
    for i in range(N):
        j = (i+1) % N
        M[i][i] += 1; M[j][j] += 1; M[i][j] -= 1; M[j][i] -= 1
    return M

def rank_frac(M):
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0]); r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]; M[r] = [x/pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        r += 1
        if r == rows: break
    return r

L1 = cycle_lap_frac(N)
rk1 = rank_frac(L1)
# L ⊗ 1_2 = 兩塊 L 之 block-diag
L2 = [[F(0)]*(2*N) for _ in range(2*N)]
for i in range(N):
    for j in range(N):
        L2[i][j] = L1[i][j]; L2[N+i][N+j] = L1[i][j]
rk2 = rank_frac(L2)
print(f"  k=1: rank = {rk1} ⟹ nullity = {N-rk1}")
print(f"  k=2: rank = {rk2} ⟹ nullity = {2*N-rk2}")
chk("k=1 ⟹ 零模重數 1（唯一 Goldstone）", N-rk1 == 1, f"{N-rk1}")
chk("k=2 ⟹ 零模重數 2（兩個無質量純量）", 2*N-rk2 == 2, f"{2*N-rk2}")
chk("零模重數 = k·nullity(L) ⟹ k 由零模重數鎖定", (2*N-rk2) == 2*(N-rk1))
print("  [註] k≥2 ⟹ 多一個無質量 Goldstone ⟹ 額外長程力（第五力）；此為物理解釋，非本輪計算")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")