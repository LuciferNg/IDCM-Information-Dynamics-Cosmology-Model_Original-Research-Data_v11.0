#!/usr/bin/env python3
"""32 — G 路線：Gemini「干涉」之**可算**形式（振幅雙線性項）

說明：docs/24 已把 Gemini 之「非對角交叉干涉」五個**矩陣元**讀法判死（R0–R4）。
      本腳本測其唯一未被覆蓋之數學形式：**振幅之雙線性干涉項**（而非矩陣元之差、亦非 G 之線性組合）。

定義（先於計算）
--------------------------------------------------
振幅：  A(i,j) = Σ_n φ_n(i)φ_n(j)/√γ_n            （G = Σ φφ/γ ⟹ 干涉項對 G 係非線性）
G1：    δ_x = C·2·A(k_d,x)·A(k_u,x)                 （實振幅之干涉叉項）
G2：    Ã(i,j) = Σ_n φ_n(i)φ_n(j)·e^{iθ_n}/√γ_n，θ_n ≡ 2πn/88（環相位，結構性非自由）
        δ_x = C·2·Re[Ã(k_d,x)·conj(Ã(k_u,x))]      （帶相位之干涉）

預註冊判準（跑前凍死）
--------------------------------------------------
F1：|z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中
F2：`N_trials = 2`（G1／G2 全報）
F3：負控制 = R11（全 G 之頂點差分，已判死 +17.96）**必須仍判死**
F4：不得再加相位形式（θ_n 已固定為環相位；改 θ 即擬合）
"""
import mpmath as mp

mp.mp.dps = 40
N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

Amp = lambda i, j: sum(phi(n, i) * phi(n, j) / mp.sqrt(gam(n)) for n in range(1, N + 1))
Gfull = lambda i, j: sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, N + 1))
th = lambda n: 2 * mp.pi * n / 88


def Atil(i, j):
    return sum(phi(n, i) * phi(n, j) * mp.e ** (1j * th(n)) / mp.sqrt(gam(n)) for n in range(1, N + 1))


print("=" * 96)
print("32 — G 路線：振幅雙線性干涉（N_trials = 2，含負控制）")
print("=" * 96)
print("  振幅 A(3,2) = %s ； A(3,3) = %s ； A(3,4) = %s"
      % (mp.nstr(Amp(3, 2), 8), mp.nstr(Amp(3, 3), 8), mp.nstr(Amp(3, 4), 8)))


def run(tag, d_of, note=""):
    du, dd, ds = d_of(KU), d_of(KD), d_of(KS)
    if ds == 0:
        print("   %-38s δ_s = 0 ⟹ 未定義 ⟹ 判死 %s" % (tag, note))
        return None
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-38s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


print("\n[G1] δ_x = 2·A(3,x)·A(2,x)（實振幅干涉）")
g1 = run("G1 實振幅叉項", lambda kx: 2 * Amp(KD, kx) * Amp(KU, kx))

print("\n[G2] δ_x = 2·Re[Ã(3,x)·Ã(2,x)*]（環相位 θ_n = 2πn/88）")
g2 = run("G2 帶相位干涉", lambda kx: 2 * mp.re(Atil(KD, kx) * mp.conj(Atil(KU, kx))))

print("\n[負控制] R11：全 G 之頂點差分（已判死 +17.96）")
c = run("CTRL 全 G 頂點差分", lambda kx: Gfull(KD, kx) - Gfull(KU, kx))

print("\n" + "=" * 96)
print("[判定]")
for tag, r in (("G1", g1), ("G2", g2)):
    print("   %s：%s" % (tag, "— 未定義" if r is None else ("★ 命中" if r[4] else "判死（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]))))
print("   負控制：%s" % ("判死 ✓ pipeline 有效" if (c is None or not c[4]) else "**反而通 ⟹ pipeline 失效**"))
print("\n[結構備註]")
print("   G1／G2 對 G 係**非線性**（A ~ G^{1/2}），故不在 R10（譜函數）、R11-R13（G 之線性組合）之已殺族內。")
print("   G2 之相位 θ_n 係環之自然相位（非自由參數）；若 G2 命中，須另立選擇性檢定（相位族是否唯一）。")
print("=" * 96)