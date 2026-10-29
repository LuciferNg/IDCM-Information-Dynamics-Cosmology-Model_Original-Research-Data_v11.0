#!/usr/bin/env python3
"""46 — R25：半整數宇稱 ℤ₂（「−」可否係對稱本徵值？）

承 R24（scripts/45）：四個物件層讀法（規範／鏡像／映射／相位）皆判死或空。
本輪檢「−」之一個**精確**對稱讀法：
    模反射 R：n → N_d+1−n = 45−n（即半整數零點對偶 n ↔ 45−n，已註冊結構）
    站空間作用：φ_{45−n}(j) = −(−1)^j · φ_n(j)  （站反射 × staggered）
    ⟹ 核之宇稱 w_{45−n} = ±w_n 給出 ℤ₂ 本徵值 ±1 —— **離散、無參數**

已註冊恆等：γ_n + γ_{45−n} = 4（γ_n = 4sin²(nπ/90)）
⟹ 宇稱定核之自然基底：
    P+ ：w_n = 1/γ_n + 1/γ_{45−n} = 4/(γ_n γ_{45−n})          （偶）
    P− ：w_n = 1/γ_n − 1/γ_{45−n} = (γ_{45−n} − γ_n)/(γ_n γ_{45−n})（奇）

預註冊（跑前凍死；N_trials = 2 ＋ 1 控制）
判準：|z(ρ_d)| < 3 且 |z(ρ_u)| < 3 ⟹ 命中
附註：若 P− 命中，則「−」= ℤ₂ 宇稱本徵值；若全判死，則「−」非任何物件層對稱。
"""
import mpmath as mp

mp.mp.dps = 30
NS = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

print("=" * 96)
print("46 — R25：半整數宇稱 ℤ₂（模反射 n → 45−n）")
print("=" * 96)

print("\n[A] 已註冊恆等核實：γ_n + γ_{45−n} = 4")
bad = 0
for n in range(1, NS + 1):
    if abs(gam(n) + gam(45 - n) - 4) > mp.mpf("1e-25"):
        bad += 1
print("   44 個模式逐一核實：不符 %d 個 ⟹ %s" % (bad, "恆等成立 ✓" if bad == 0 else "**不成立**"))

print("\n[B] 站空間作用：φ_{45−n}(j) = −(−1)^j φ_n(j)")
dev = max(abs(phi(45 - n, j) - (-1) ** j * (-1) * phi(n, j)) for n in range(1, NS + 1) for j in (2, 3, 4))
print("   站 j=2,3,4 之最大偏差 = %s ⟹ %s" % (mp.nstr(dev, 6), "恆等成立 ✓" if dev < mp.mpf("1e-25") else "不成立"))
print("   ⟹ 模反射 = **站反射 × staggered**；且 Δ_n = φ_n(3)−φ_n(2) 反射後變 φ_n(3)+φ_n(2)（差型 → 和型）")


def run(tag, w, note=""):
    G = lambda i, j: mp.fsum(w(n) * phi(n, i) * phi(n, j) for n in range(1, NS + 1))
    ds = G(KD, KS) - G(KU, KS)
    if ds == 0:
        print("   %-22s δ_s = 0 ⟹ 未定義 ⟹ 判死" % tag)
        return None
    pd = (G(KD, KD) - G(KU, KD)) / ds
    pu = (G(KD, KU) - G(KU, KU)) / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = abs(zd) < 3 and abs(zu) < 3
    print("   %-22s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


print("\n[C] 兩個宇稱定核 ＋ 控制")
r_plus = run("P+ 1/γ_n+1/γ_45−n", lambda n: 1 / gam(n) + 1 / gam(45 - n), "（偶）")
r_minus = run("P− 1/γ_n−1/γ_45−n", lambda n: 1 / gam(n) - 1 / gam(45 - n), "（奇）")
r_ctrl = run("C 1/γ_n", lambda n: 1 / gam(n), "（= R17 f0，已知判死）")

print("\n[D] 結構附註")
Wp = [1 / gam(n) + 1 / gam(45 - n) for n in range(1, NS + 1)]
Wm = [1 / gam(n) - 1 / gam(45 - n) for n in range(1, NS + 1)]
print("   P+ 等價形 4/(γ_n γ_{45−n})：與上式最大偏差 = %s"
      % mp.nstr(max(abs(Wp[i] - 4 / (gam(i + 1) * gam(45 - i - 1))) for i in range(NS)), 6))
print("   P+ 之值域 [%s, %s]；P− 之值域 [%s, %s]"
      % (mp.nstr(min(Wp), 6), mp.nstr(max(Wp), 6), mp.nstr(min(Wm), 6), mp.nstr(max(Wm), 6)))
print("   P− 之符號：正 %d、負 %d ⟹ %s"
      % (sum(1 for x in Wm if x > 0), sum(1 for x in Wm if x < 0),
         "**本質帶號（奇宇稱）**" if any(x < 0 for x in Wm) else "全正"))

print("\n[E] 精度界限：載體自然格 ℚ/41 與量測解析度")
lat = mp.mpf(1) / 41
print("   本載體之比值恆落於 ℚ/41（實測：1/γ → 42/41、−2/41；P+ → −63/41、43/41；P− → −28/41、+28/41）")
print("   格距 1/41 = %s" % mp.nstr(lat, 10))
print("   σ_ρd = %s ⟹ σ_d / 格距 = %s" % (mp.nstr(RD_S, 6), mp.nstr(RD_S / lat, 6)))
print("   σ_ρu = %s ⟹ σ_u / 格距 = %s" % (mp.nstr(RU_S, 6), mp.nstr(RU_S / lat, 6)))
print("   1σ 窗寬 = %s 格 ⟹ 每 1σ 覆蓋約 %s 個格點"
      % (mp.nstr(2 * RD_S / lat, 4), mp.nstr(2 * RD_S / lat + 1, 4)))
print("   ρ_d 之最近格點：16/41 = %s（z %+.4f）、17/41 = %s（z %+.4f）"
      % (mp.nstr(mp.mpf(16) / 41, 8), (mp.mpf(16) / 41 - RD) / RD_S,
         mp.nstr(mp.mpf(17) / 41, 8), (mp.mpf(17) / 41 - RD) / RD_S))
print("   ρ_u 之最近格點：−1/41 = %s（z %+.4f）、−2/41 = %s（z %+.4f）"
      % (mp.nstr(-mp.mpf(1) / 41, 8), (-mp.mpf(1) / 41 - RU) / RU_S,
         mp.nstr(-mp.mpf(2) / 41, 8), (-mp.mpf(2) / 41 - RU) / RU_S))
print("   ★ 半格 = %s < σ_ρu = %s 且 < σ_ρd = %s" % (mp.nstr(lat / 2, 8), mp.nstr(RU_S, 6), mp.nstr(RD_S, 6)))
print("   ⟹ **任何目標值必落於某格點之半格內，而半格 < σ ⟹ 必在 1σ 內**")
print("   ⟹ 格點吻合之機率 = **100%（必然）** ⟹ 「ρ 等於某結構值」之吻合**零資訊**")

print("\n" + "=" * 96)
print("[判定]")
if r_minus and r_minus[4]:
    print("   ★ P−（奇宇稱）命中 ⟹ 「−」= 半整數位移之 ℤ₂ 宇稱本徵值（物件層對稱）")
else:
    print("   P−（奇宇稱）：判死 ⟹ 「−」**不是**半整數宇稱本徵值")
print("   判定：%s" % ("宇稱讀法成立" if (r_plus and r_plus[4]) or (r_minus and r_minus[4]) else "宇稱讀法判死；「−」非物件層對稱，續歸『配對之區分』（docs/37）"))
print("=" * 96)