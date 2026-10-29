#!/usr/bin/env python3
"""R4-P2 — 定向判別（對稱荷 vs 有向荷），零參數，無擬合。

預註冊判準（2026-10-29 跑前凍結，不得事後改）
--------------------------------------------------
O1 對稱二體荷判別：若配對荷對味對正反不敏感（a(q q̄) = a(q q)，等價於相關子
   G_ij = G_ji 實對稱），則由價對多重集
       Δm(Ξ⁻−Ξ⁰) = 2[G(d,s) − G(u,s)]
       Δm(K⁺−K⁰) = G(u,s̄) − G(d,s̄) = G(u,s) − G(d,s) = −½ Δm(Ξ⁻−Ξ⁰)
   恆等式： r(K⁺−K⁰) = −½ · r(Ξ⁻−Ξ⁰)
   判死門檻：|z| ≥ 3

O2 若 O1 死 ⟹ 「實對稱相關子」形式族全部出局（含一切 G(|i−j|) 距離核）。
   同時量出必須之破缺量 ρ ≡ |r_K| / |r_Ξ|/2，作未來有向結構之靶。

O3 可辨識性計數（唔准用擬合繞過）：未知量 = 13 個味對型荷（6 qq + 7 q q̄）
   之相對值（絕對標度錨定於 n−p）；方程 = 5 個獨立比值（n−p 已作標度）。
   欠定度 = 13 − 5。以此判「自由荷」是否可唯一化。

O4 光滑核檢：G_ij = Σ_(n=1..44) φ_n(i) φ_n(j) / γ_n，
       φ_n(j) = √(2/45)·sin(πnj/45)（44-site 鏈 Dirichlet 標準基），
       γ_n = 4sin²(nπ/88)（IDCM 44 模態譜）。
   模式指派（O209/O233 已裁）：u→2, d→3, s→4（敏感度：s∈{4,5,7} 只報）。
   判死門檻：任一比值 |dev| > 10%（相對目標）且無單一變號機制 ⟹ 光滑核死。

紀律：禁止調參、禁止掃描、禁止擬合；O4 只報不調。浮點僅作顯示，
      核心用 mpmath dps=40。
"""
import mpmath as mp

mp.mp.dps = 40

# 目標比值（R4-P1 已獨立算出；此處重算，不讀取 M12 產物）
MASS = {
    "n": ("939.565421", "0.0000001"), "p": ("938.272089", "0.0000001"),
    "S-": ("1197.449", "0.030"), "S0": ("1192.642", "0.024"), "S+": ("1189.37", "0.07"),
    "X-": ("1321.71", "0.07"), "X0": ("1314.86", "0.20"),
    "K+": ("493.677", "0.016"), "K0": ("497.611", "0.013"),
}


def dm(a, b):
    return mp.mpf(MASS[a][0]) - mp.mpf(MASS[b][0])


d_np = dm("n", "p")
TGT = {
    "r(Σ⁻−Σ⁰)": dm("S-", "S0") / d_np,
    "r(Ξ⁻−Ξ⁰)": dm("X-", "X0") / d_np,
    "r(Σ⁻−Σ⁺)": dm("S-", "S+") / d_np,
    "r(K⁺−K⁰)": dm("K+", "K0") / d_np,
}
SIG = {"r(Σ⁻−Σ⁰)": mp.mpf("0.029705"), "r(Ξ⁻−Ξ⁰)": mp.mpf("0.163837"),
       "r(Σ⁻−Σ⁺)": mp.mpf("0.058885"), "r(K⁺−K⁰)": mp.mpf("0.015940")}

print("=" * 78)
print("R4-P2 定向判別 — 對稱荷（G 對稱）vs 有向荷（q q̄ 有別）")
print("=" * 78)

print("\n[目標] 5 個獨立比值（本腳本獨立重算）")
for k in ["r(Σ⁻−Σ⁰)", "r(Ξ⁻−Ξ⁰)", "r(Σ⁻−Σ⁺)", "r(K⁺−K⁰)"]:
    print("  %-12s = %+10.6f ± %.6f" % (k, TGT[k], SIG[k]))

# ---------------------------------------------------------------- O1
print("\n" + "=" * 78)
print("[O1] 對稱二體荷恆等式： r(K⁺−K⁰) = −½·r(Ξ⁻−Ξ⁰)")
pred = -mp.mpf("0.5") * TGT["r(Ξ⁻−Ξ⁰)"]
res = TGT["r(K⁺−K⁰)"] - pred
sig = mp.sqrt(SIG["r(K⁺−K⁰)"] ** 2 + (mp.mpf("0.5") * SIG["r(Ξ⁻−Ξ⁰)"]) ** 2)
z1 = res / sig
print("  預言 = %+10.6f   實測 = %+10.6f   殘差 = %+10.6f ± %.6f" % (pred, TGT["r(K⁺−K⁰)"], res, sig))
print("  z = %+.2f   ⟹ %s" % (z1, "對稱荷族判死 🔴" if abs(z1) >= 3 else "未死 🟡"))

# O2 破缺量
rho = abs(TGT["r(K⁺−K⁰)"]) / (abs(TGT["r(Ξ⁻−Ξ⁰)"]) / 2)
print("\n[O2] 必要之正反破缺量 ρ ≡ |r_K| / (|r_Ξ|/2) = %.6f" % rho)
print("     即 a(q q̄) 之正反不對稱須達 %.2f%% 才可同時給 −3.0418 與 +5.2964。" % (abs(rho - 1) * 100))

# ---------------------------------------------------------------- O3
print("\n" + "=" * 78)
print("[O3] 可辨識性（自由味對荷）")
n_type = 6 + 7          # qq: dd,ds,du,ss,su,uu ; q q̄: Dc,Dd,Du,Sd,Su,Uc,Uu
n_eq = 4                # 4 個非標度比值（n−p 作標度）
print("     未知味對型 = %d（6 qq + 7 q q̄）" % n_type)
print("     方程       = %d 個獨立比值（n−p 已作絕對標度）" % n_eq)
print("     欠定度     = %d  ⟹ 自由荷不可唯一化（違 G4 零參數）" % (n_type - n_eq))

# ---------------------------------------------------------------- O4
print("\n" + "=" * 78)
print("[O4] 光滑核檢：44-site 鏈 Dirichlet Green 函數（IDCM γ 譜）")


def gamma(n):
    """IDCM 44 模態譜 γ_n = 4 sin²(nπ/88), n=1..44"""
    return 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2


def phi(n, j, N=44):
    """44-site 鏈 Dirichlet 標準正交基（DST-I）"""
    return mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * mp.mpf(j) / (N + 1))


def Gloc(i, j, N=44):
    """G_ij = Σ_n φ_n(i) φ_n(j) / γ_n"""
    s = mp.mpf(0)
    for n in range(1, N + 1):
        s += phi(n, i, N) * phi(n, j, N) / gamma(n)
    return s


ASSIGN = {"u": 2, "d": 3, "s": 4}


def ratios(asg):
    iu, idd, iss = asg["u"], asg["d"], asg["s"]
    X = Gloc(idd, idd) - Gloc(iu, iu)                       # n−p
    dS = Gloc(idd, idd) + Gloc(idd, iss) - Gloc(iu, idd) - Gloc(iu, iss)
    dX = 2 * (Gloc(idd, iss) - Gloc(iu, iss))               # 對稱 ⟹ Ξ 通道
    dSSp = X + dX                                           # 恆等式，必然成立
    dK = -(dX / 2)                                          # 對稱 ⟹ K 通道 = −Ξ/2
    return {"r(Σ⁻−Σ⁰)": dS / X, "r(Ξ⁻−Ξ⁰)": dX / X,
            "r(Σ⁻−Σ⁺)": dSSp / X, "r(K⁺−K⁰)": dK / X}


for label, asg in [("主方案 s=4", ASSIGN), ("敏感度 s=5", {"u": 2, "d": 3, "s": 5}),
                   ("敏感度 s=7", {"u": 2, "d": 3, "s": 7})]:
    R = ratios(asg)
    print("\n  ── %s ──" % label)
    worst = mp.mpf(0)
    for k in ["r(Σ⁻−Σ⁰)", "r(Ξ⁻−Ξ⁰)", "r(Σ⁻−Σ⁺)", "r(K⁺−K⁰)"]:
        dev = (R[k] - TGT[k]) / abs(TGT[k]) * 100
        worst = max(worst, abs(dev))
        print("     %-12s 預言 %+12.6f   目標 %+10.6f   dev = %+9.2f%%" % (k, R[k], TGT[k], dev))
    print("     max |dev| = %.2f%%  ⟹ %s" % (worst, "光滑核判死 🔴" if worst > 10 else "未死 🟡"))

# ---------------------------------------------------------------- 判詞
print("\n" + "=" * 78)
print("[判詞]")
print("  O1 對稱二體荷: z = %+.2f ⟹ %s" % (z1, "排除 🔴" if abs(z1) >= 3 else "未排除"))
print("  O2 正反破缺量須 ρ = %.4f（偏離 1 約 %.1f%%）" % (rho, abs(rho - 1) * 100))
print("  O3 自由荷欠定 %d ⟹ 不可作結構（違零參數）" % (n_type - n_eq))
print("  O4 光滑 Green 核：見上表（判死門檻 10%）")
print("=" * 78)