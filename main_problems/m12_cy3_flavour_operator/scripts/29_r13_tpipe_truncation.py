#!/usr/bin/env python3
"""R13 — S2 四攻：T_pipe 截斷核（結構性截斷，非自由形式；`N_trials = 2`）

為何係截斷核（結構依據，先於計算）
--------------------------------------------------
R12b 結論收窄為兩條**同時**成立：
  (i) 非平移不變（否則 `ρ_u = −ρ_d` 恆等，實測偏 15.58× ⟹ 全 TI 族死）
  (ii) 短程差分比達 `2.4783`（Green 族只得 0.967–0.984 ⟹ 差 2.5×，需 `n=1→2` 有「崖」）
載體中**唯一已註冊**之「崖」= `T_pipe` 之模式截斷：
  `n* = N_d − L_6 + 1 = 44 − 18 + 1 = 27`（`L_6 = φ⁶+φ⁻⁶ = 18`）
  `ker(T_pipe) = {n ≥ 28}` = **17 模式**（`F_proj = N_d − n* = 17`）
沿用 R11/R12 之判準結構（先驗形式、`N_trials=2`、含負控制），**唔掃 `n*`**。

預註冊協議（跑前凍結）
--------------------------------------------------
A1 世界側投影核：`G_proj(i,j) = Σ_{n=1}^{27} φ_n(i)φ_n(j)/γ_n`
   `δ_x = C[G_proj(3,k_x) − G_proj(2,k_x)]`，`k_u=2, k_d=3, k_s=4`
A2 互補（核）側：`G_ker(i,j) = Σ_{n=28}^{44} φ_n(i)φ_n(j)/γ_n`，同式
判準
  K1：A1 之 `|z(ρ_d)| < 3` 且 `|z(ρ_u)| < 3` ⟹ 命中
  K2：A2 全報（第二協議；`N_trials = 2`，若僅 A2 中須報選擇效應）
  K3：負控制 = 全 `G`（R11 已判死 `z=+17.96`）**必須仍判死**，否則 pipeline 失效
  K4：`C` 之結構性；及**勿**掃 `n*`（掃即擬合）；只可報「達標所需 `n*`」作敏感度，不作證據
禁：不得改 `n*`、不得改 BC、不得選 `k_s`、不得加任何權重函數。
"""
import mpmath as mp

mp.mp.dps = 40
N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
DS, DB, DC = mp.mpf("3.4250"), mp.mpf("1.8933"), mp.mpf("-0.2820")
KU, KD, KS = 2, 3, 4
NSTAR = 27          # = N_d − L_6 + 1 = 44 − 18 + 1

print("=" * 96)
print("R13 — S2 四攻：T_pipe 截斷核（N_trials = 2，含負控制；不掃 n*）")
print("=" * 96)
print("  n* = N_d − L_6 + 1 = 44 − 18 + 1 = %d ；ker(T_pipe)={n≥28} = %d 模式（F_proj=17）"
      % (NSTAR, N - NSTAR))
print("  模式能量分解核實：Σ_{n≤27}γ_n = %s ；E_proj=28.5、E_sync=61.5"
      % mp.nstr(sum(gam(n) for n in range(1, NSTAR + 1)), 8))


def make_G(lo, hi):
    cache = {}

    def G(i, j):
        if (i, j) not in cache:
            cache[(i, j)] = sum(phi(n, i) * phi(n, j) / gam(n) for n in range(lo, hi + 1))
        return cache[(i, j)]
    return G


def run(tag, G, ks=KS):
    du = G(KD, KU) - G(KU, KU)
    dd = G(KD, KD) - G(KU, KD)
    ds = G(KD, ks) - G(KU, ks)
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("\n  [%s] k_u=%d k_d=%d k_s=%d" % (tag, KU, KD, ks))
    print("     δ_u=%s  δ_d=%s  δ_s=%s" % (mp.nstr(du, 8), mp.nstr(dd, 8), mp.nstr(ds, 8)))
    print("     ρ_d = %s（z %+7.2f）   ρ_u = %s（z %+7.2f） ⟹ %s"
          % (mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "不中"))
    return du, dd, ds, pd, pu, zd, zu, ok


print("\n" + "=" * 96)
Gproj = make_G(1, NSTAR)
Gker = make_G(NSTAR + 1, N)
Gfull = make_G(1, N)
a1 = run("A1 世界側 G_proj（n≤27）", Gproj)
a2 = run("A2 核側 G_ker（n≥28，17 模式）", Gker)
FPROJ = (mp.mpf(57) / 2 - sum(gam(n) for n in range(1, 27))) / gam(27)
WT = {n: (mp.mpf(1) if n <= 26 else (FPROJ if n == 27 else mp.mpf(0))) for n in range(1, N + 1)}


def Gtp(i, j):
    return sum(WT[n] * phi(n, i) * phi(n, j) / gam(n) for n in range(1, N + 1))


print("\n  ★ 精確 T_pipe 權重核實：f_proj = (57/2 − Σ_{1..26}γ)/γ_27 = %s（註冊 0.771960）"
      % mp.nstr(FPROJ, 8))
a3 = run("A3 精確 T_pipe 權重（n≤26 全權 + 0.77196×mode27）", Gtp)
c0 = run("負控制 全 G（R11 已死）", Gfull)

print("\n" + "=" * 96)
print("[判定]")
print("   K1（A1）：%s" % ("命中 🟢" if a1[7] else "不中 🔴（z_d = %+0.2f）" % a1[5]))
print("   K2（A2）：%s" % ("命中（第二協議，須報選擇效應）" if a2[7] else "不中（z_d = %+0.2f）" % a2[5]))
print("   K2b（A3 精確 T_pipe）：%s" % ("命中" if a3[7] else "不中（z_d = %+0.2f、z_u = %+0.2f）" % (a3[5], a3[6])))
print("   K3（負控制）：%s" % ("判死 ✓ pipeline 有效" if ((abs(c0[5]) >= 3) or (abs(c0[6]) >= 3)) else "**反而通 ⟹ pipeline 失效**"))

print("\n" + "=" * 96)
print("[K4] 標度與敏感度（後驗，不作證據）")
C = DS / a1[2]
print("   C = δ_s/[G_proj(3,4) − G_proj(2,4)] = %s MeV" % mp.nstr(C, 8))
for nm, dv in (("δ_b", DB), ("δ_c", DC)):
    tgt = dv / C
    best = min(range(1, 45), key=lambda k: abs((Gproj(KD, k) - Gproj(KU, k)) - tgt) / abs(tgt))
    rel = (Gproj(KD, best) - Gproj(KU, best)) / tgt
    print("   由 %s 反解：最近整數 k = %d（值/目標 = %s）⟹ %s"
          % (nm, best, mp.nstr(rel, 6), "整數解存在" if abs(rel - 1) < mp.mpf("0.05") else "**無整數解**"))

print("\n   截斷點敏感度（**僅供參考，掃 n* 即擬合、不得作判據**）：")
for ns in (7, 13, 18, 22, 27, 33, 39, 44):
    G = make_G(1, ns)
    du = G(KD, KU) - G(KU, KU)
    dd = G(KD, KD) - G(KU, KD)
    ds = G(KD, KS) - G(KU, KS)
    pd, pu = dd / ds, du / ds
    print("     n*=%2d：ρ_d = %-12s（z %+7.2f）  ρ_u = %-12s（z %+7.2f）"
          % (ns, mp.nstr(pd, 8), (pd - RD) / RD_S, mp.nstr(pu, 8), (pu - RU) / RU_S))

print("\n" + "=" * 96)
print("[匯總 — 只報事實]")
print("   A1（n≤27）：ρ_d=%s（z %+0.2f）、ρ_u=%s（z %+0.2f）" % (mp.nstr(a1[3], 8), a1[5], mp.nstr(a1[4], 8), a1[6]))
print("   A2（n≥28）：ρ_d=%s（z %+0.2f）、ρ_u=%s（z %+0.2f）" % (mp.nstr(a2[3], 8), a2[5], mp.nstr(a2[4], 8), a2[6]))
print("   A3（精確 T_pipe 權重）：ρ_d=%s（z %+0.2f）、ρ_u=%s（z %+0.2f）" % (mp.nstr(a3[3], 8), a3[5], mp.nstr(a3[4], 8), a3[6]))
print("   ⟹ 判詞：%s" % ("A1 命中 ⟹ 零參數形狀成立（標度待解）" if a1[7] else
                        "A1 不中 ⟹ 「T_pipe 截斷核」判死；S2 再收窄"))
print("=" * 96)