#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/04_closure_scan.py — 閉包集 C1–C3 之完整枚舉（N_trials = 16，全報）
＋ ★ A₂ Cartan 結構（近正交之結構理由）＋ 剩餘 dof 定位

判準（跑前凍結）：scripts/02 之 G0–G2、T1–T5。
命中定義：通過 G1、G2 且 |z_d| < 3 且 |z_u| < 3。
N_trials = 16（全部列出，不得只報命中）。look-elsewhere：p ≈ N × P(hit)。

閉包（有限且完備，逐個附出處）：
  C1（7）w = K⁻¹·u，u ∈ {v_u, v_d, v_s, Δ, φ(2), φ(3), φ(4)}  —— 鏈傳播子作用於已註冊結構向量
  C2（6）w_n = φ_n(a)φ_n(b)/γ_n，2 ≤ a ≤ b ≤ 4               —— 三點譜核（傳播子夾兩個已註冊站）
  C3（3）w_n = Δ_n φ_n(a)/γ_n，a ∈ {2,3,4}                    —— Δ-加權三點核

結構副產物：
  · Dirichlet Green 閉式 G(i,j) = min(i,j)(45−max(i,j))/45 —— 精確 Fraction 核 L G = I（1936 格）
  · K⁻¹-Gram 矩陣 M_xy = v_xᵀK⁻¹v_y —— 60 dps ＋ PSLQ ＋ 120 dps 覆核
    ⟹ M = (1/90)·[[2,−1,0],[−1,2,0],[0,0,2]] = (1/90)·(A₂ Cartan ⊕ 1)
"""
import sys
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 60
N = 44
STATION = {"u": 2, "d": 3, "s": 4}

phi = lambda n, j: mp.sqrt(mp.mpf(2)/45) * mp.sin(mp.pi*n*j/45)
DELTA = [phi(n, 3) - phi(n, 2) for n in range(1, N+1)]
GAM = [4*mp.sin(mp.pi*n/90)**2 for n in range(1, N+1)]
V = {x: [phi(n, k)*DELTA[n-1] for n in range(1, N+1)] for x, k in STATION.items()}

RHO_D_OBS = mp.mpf("0.40350365")
RHO_U_OBS = mp.mpf("-0.02588837")
SIG_RD = mp.mpf("0.044785178")
SIG_RU = mp.mpf("0.033583899")

# ---------------- 0. Dirichlet Green 閉式（精確 Fraction）----------------
print("=== 0. Dirichlet Green 閉式：G(i,j) = min(i,j)(45−max(i,j))/45 ===")
def G_closed(i, j):
    lo, hi = min(i, j), max(i, j)
    return F(lo*(45-hi), 45)

maxdev, bad = 0, 0
for i in range(1, N+1):
    for j in range(1, N+1):
        gg = lambda a: (G_closed(a, j) if 1 <= a <= N else F(0))
        val = 2*G_closed(i, j) - gg(i-1) - gg(i+1)
        tgt = F(1) if i == j else F(0)
        d = abs(val - tgt)
        if d:
            bad += 1
        maxdev = max(maxdev, d)
print(f"  (LG = I) 最大偏差 = {maxdev}（Fraction 精確）；不符格數 = {bad}/{N*N}")
print(f"  G(2,2)={G_closed(2,2)}  G(2,3)={G_closed(2,3)}  G(3,3)={G_closed(3,3)}  "
      f"G(3,4)={G_closed(3,4)}  G(4,2)={G_closed(4,2)}")
sp = mp.fsum(phi(n, 2)*phi(n, 3)/GAM[n-1] for n in range(1, N+1))
print(f"  譜和覆核 Σ_n φ_n(2)φ_n(3)/γ_n = {mp.nstr(sp,20)}  vs  28/15 = {float(F(28,15)):.20f}")

# ---------------- 1. K⁻¹-Gram ----------------
print("\n=== 1. M_xy = v_xᵀ K⁻¹ v_y（60 dps ＋ PSLQ ＋ 120 dps 覆核）===")
def bilin(a, b, dps=60):
    old = mp.mp.dps
    mp.mp.dps = dps
    r = mp.fsum(V[a][n-1]*V[b][n-1]/GAM[n-1] for n in range(1, N+1))
    mp.mp.dps = old
    return r

Ms = {(a, b): bilin(a, b) for a in "uds" for b in "uds"}
for a in "uds":
    print("  " + a + ": " + "  ".join(f"{mp.nstr(Ms[(a,b)],20):>24s}" for b in "uds"))

exact = {}
for k, v in Ms.items():
    r = mp.pslq([v, mp.mpf(1)], maxcoeff=10**12, maxsteps=10**6)
    if r:
        a_, b_ = int(r[0]), int(r[1])           # a_·v + b_ = 0  ⟹  v = −b_/a_
        exact[k] = F(-b_, a_)
        ok = abs(mp.mpf(-b_)/a_ - v) < mp.mpf("1e-45")
        print(f"  M[{k[0]},{k[1]}] = {mp.nstr(v,20)} ⟹ {exact[k]}  {'✅' if ok else '⚠'}")
    else:
        exact[k] = F(0)
        print(f"  M[{k[0]},{k[1]}] = {mp.nstr(v,20)} ⟹ 精確 0（PSLQ 無有理解）")

print("\n  120 dps 覆核（精確零須隨 dps 下降，不停在地板）：")
for k in (("u", "s"), ("d", "s"), ("u", "u"), ("u", "d")):
    hi = bilin(k[0], k[1], dps=120)
    tgt = mp.mpf(exact[k].numerator)/exact[k].denominator
    print(f"    M[{k[0]},{k[1]}] = {mp.nstr(hi,10)}   |M−閉式| = {mp.nstr(abs(hi-tgt),4)}")

# ---------------- 2. 閉包枚舉 ----------------
def rho_of(w):
    f = lambda x: mp.fsum(w[n-1]*V[x][n-1] for n in range(1, N+1))
    du, dd, ds = f("u"), f("d"), f("s")
    if abs(ds) < mp.mpf("1e-30"):
        return None
    return (du/ds, dd/ds, du, dd, ds)

PHI = {j: [phi(n, j) for n in range(1, N+1)] for j in (2, 3, 4)}
base = {"v_u": V["u"], "v_d": V["d"], "v_s": V["s"], "Δ": DELTA,
        "φ(2)": PHI[2], "φ(3)": PHI[3], "φ(4)": PHI[4]}

cands = []
for name_, u in base.items():
    cands.append((f"C1 w=K⁻¹·{name_}", [u[n-1]/GAM[n-1] for n in range(1, N+1)]))
for a in (2, 3, 4):
    for b in (2, 3, 4):
        if b < a:
            continue
        cands.append((f"C2 w=φ({a})φ({b})/γ", [phi(n, a)*phi(n, b)/GAM[n-1] for n in range(1, N+1)]))
for a in (2, 3, 4):
    cands.append((f"C3 w=Δφ({a})/γ", [DELTA[n-1]*phi(n, a)/GAM[n-1] for n in range(1, N+1)]))

print(f"\n=== 2. 閉包枚舉：N_trials = {len(cands)}（全報）===")
print(f"{'候選':24s} {'ρ_u':>16s} {'ρ_d':>16s} {'z_u':>9s} {'z_d':>9s}  判")
hits = []
for nm, w in cands:
    r = rho_of(w)
    if r is None:
        print(f"{nm:24s} {'δ_s=0（退化）':>16s} {'—':>16s} {'—':>9s} {'—':>9s}  🔴退化")
        continue
    ru, rd, du, dd, ds = r
    zu = (ru - RHO_U_OBS)/SIG_RU
    zd = (rd - RHO_D_OBS)/SIG_RD
    ok = abs(zd) < 3 and abs(zu) < 3
    if ok:
        hits.append((nm, ru, rd, zu, zd))
    print(f"{nm:24s} {mp.nstr(ru,10):>16s} {mp.nstr(rd,10):>16s} "
          f"{mp.nstr(zu,6):>9s} {mp.nstr(zd,6):>9s}  {'🟡過閘' if ok else '🔴'}")
print(f"\n  命中（|z|<3 兩者）= {len(hits)} / {len(cands)}")

# ---------------- 3. A₂ Cartan 結構 ----------------
print("\n=== 3. ★★ A₂ Cartan 結構：近正交之結構理由 ＋ 剩餘 dof ===")
M = [[exact[("u","u")], exact[("u","d")], exact[("u","s")]],
     [exact[("d","u")], exact[("d","d")], exact[("d","s")]],
     [exact[("s","u")], exact[("s","d")], exact[("s","s")]]]
print("  M（精確）= 1/90 · [[2,−1,0],[−1,2,0],[0,0,2]]：")
for r in M:
    print("     ", r)
print("  ⟹ 於 K_chain 傳播子度量下：{v_u, v_d, v_s} = A₂ Cartan ⊕ 單態")
print("      · v_s ⊥_{K⁻¹} span{v_u,v_d}（兩交叉項精確 0）")
print("      · 2×2 塊 [[2,−1],[−1,2]]/90 = A₂ Cartan / 90")

print("\n  對 u = α_u v_u + α_d v_d + α_s v_s，核 w = K⁻¹u ⟹ δ = Mα：")
print("     δ_u = (2α_u − α_d)/90 ,  δ_d = (2α_d − α_u)/90 ,  δ_s = 2α_s/90")
print("     ρ_u = (2α_u − α_d)/(2α_s) ,  ρ_d = (2α_d − α_u)/(2α_s)")
print("\n  ★ 推論 1（近正交之結構理由）：ρ_u = 0 ⟺ α_d = 2α_u（Cartan 首行之零方向）")
print("     ⟹ α ∝ (1, 2, α_s)，**無需調參**，由 A₂ Cartan 強制。")
print("     （M12 docs/37：『真正缺者 = 使 ⟨w,v_u⟩ 近零之結構理由』⟹ 本輪給出：Cartan 零方向。）")
print("\n  ★ 推論 2（剩餘自由度）：代入 α_d = 2α_u ⟹ ρ_d = 3α_u/(2α_s)")
need = mp.mpf(3)/(2*RHO_D_OBS)
print(f"     需 α_s/α_u = 3/(2ρ_d) = {mp.nstr(need,10)}  ⟹ **恰 1 個連續自由度**")
print("     已註冊常數中無此值；鄰近之簡單有理/代數值全部落禁區 3：")
for nm, val in [("α_s=B=4 ⟹ ρ_d=3/8", mp.mpf(3)/8),
                ("α_s=15/4 ⟹ ρ_d=0.4", mp.mpf("0.4")),
                ("α_s=√5 ⟹ ρ_d", mp.mpf(3)/(2*mp.sqrt(5))),
                ("α_s=12/5 ⟹ ρ_d=5/8", mp.mpf("0.625"))]:
    print(f"       {nm:26s} ρ_d = {mp.nstr(val,8)}   z = {mp.nstr((val-RHO_D_OBS)/SIG_RD,4)}")
print(f"       （禁）r_Σ = 3.716756 與 needed {mp.nstr(need,8)} 差 0.02% ⟹ 禁區 3 後驗凑數，拒絕")
print("\n  ⟹ 判：**T2（dof = 0）FAIL** —— 缺一個固定 α_s 之結構量。")
print("     該量須來自 K_tree 側（A₂⊕1 之塊結構完全由 K_chain 給出，鏈內無此標度）。")
print("  ⟹ 缺口定位（本 M 最精確表述）：**A₂ 面內方向 vs 單態方向之相對歸一化**。")

# ---------------- 4. 判決 ----------------
print("\n=== 4. 判決 ===")
print(f"  C1–C3 全族 N_trials = {len(cands)}；命中 = {len(hits)}")
print(f"  look-elsewhere（若隨機族 P(雙1σ)≈3e-4）：族級 p ≈ {len(cands)*3e-4:.3g}")
print("  ⟹ 本輪**不作命中宣稱**：scoped negative（三個閉包）＋ 結構結果（A₂ Cartan ⊕ 1）")
print("     ＋ dof 定位（缺 α_s 之結構歸一化，須出自 K_tree 側）。")
print("\n  N_trials 聲明：C1 = 7、C2 = 6、C3 = 3，總 16，全部列出（無篩選）。")
sys.exit(0)