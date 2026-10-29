#!/usr/bin/env python3
"""33 — E 路線：味居於**邊（割）**而非頂點

本體論依據（用戶）：Reality = Relations + Distinctions；**區分係割，唔係點**。
故 u／d 不應係鏈上模式（頂點 2、3），而應係**割（bond）**。本輪將 u／d／s 重定義為 bond。

預註冊定義（跑後不得改）
--------------------------------------------------
bond 態：  |b_p⟩ = (|p⟩ + |p+1⟩)/√2 ，p = 1..43
bond 格林：G_b(p,q) = ⟨b_p|K⁻¹|b_q⟩ = Σ_n (φ_n(p)+φ_n(p+1))(φ_n(q)+φ_n(q+1))/(2γ_n)
位址：     u → bond 2 （= 邊 2–3）；d → bond 3（= 邊 3–4）；s → bond 4（= 邊 4–5）
E1（割能量差）：δ_x = C[G_b(k_d,x) − G_b(k_u,x)]
E2（2×2 塊對角化／避免交叉）：Δ_x = G_b(k_d,x) − G_b(k_u,x)，t ≡ G_b(k_u,k_d)（固定，非自由）
              δ_x = sgn(Δ_x)·√(Δ_x² + 4t²)      ← 唯一分支規則：絕熱跟隨（不另試）
E3（plaquette ℤ₂，結構控制）：鏈係二部 ⟹ 無奇環 ⟹ 三點通量恆 +1 ⟹ δ ≡ 0

預註冊判準
--------------------------------------------------
F1：|z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中
F2：`N_trials = 2`（E1／E2 全報；E3 屬結構控制，不計入）
F3：E3 必須**恆等零**（否則鏈非二部 ⟹ 載體定義有誤）
F4：分支規則只得「絕熱跟隨」一種；失敗**不得**改規則或改 bond 位址
"""
import mpmath as mp

mp.mp.dps = 40
N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4          # 現為 bond 編號


def Gb(p, q):
    return sum((phi(n, p) + phi(n, p + 1)) * (phi(n, q) + phi(n, q + 1)) / (2 * gam(n))
               for n in range(1, N + 1))


T_UD = Gb(KU, KD)             # u–d 耦合（固定，非自由）

print("=" * 96)
print("33 — E 路線：味居於邊（bond）；N_trials = 2")
print("=" * 96)
print("  bond 態 |b_p⟩ = (|p⟩+|p+1⟩)/√2 ；u→bond2, d→bond3, s→bond4")
print("  G_b(2,2)=%s  G_b(3,3)=%s  G_b(2,3)=%s  G_b(2,4)=%s  G_b(3,4)=%s"
      % (mp.nstr(Gb(2, 2), 8), mp.nstr(Gb(3, 3), 8), mp.nstr(Gb(2, 3), 8),
         mp.nstr(Gb(2, 4), 8), mp.nstr(Gb(3, 4), 8)))
print("  ★ u–d 耦合 t = G_b(2,3) = %s" % mp.nstr(T_UD, 8))


def judge(tag, du, dd, ds, note=""):
    if ds == 0:
        print("   %-34s δ_s = 0 ⟹ 未定義 ⟹ 判死 %s" % (tag, note))
        return None
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-34s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


print("\n[E1] 割能量差：δ_x = G_b(3,x) − G_b(2,x)")
du1 = Gb(KD, KU) - Gb(KU, KU)
dd1 = Gb(KD, KD) - Gb(KU, KD)
ds1 = Gb(KD, KS) - Gb(KU, KS)
e1 = judge("E1 割能量差", du1, dd1, ds1)

print("\n[E2] 2×2 塊對角化（避免交叉）：δ_x = sgn(Δ_x)·√(Δ_x² + 4t²)")
dlt = lambda kx: Gb(KD, kx) - Gb(KU, kx)
ddof = lambda kx: (1 if dlt(kx) >= 0 else -1) * mp.sqrt(dlt(kx) ** 2 + 4 * T_UD ** 2)
du2, dd2, ds2 = ddof(KU), ddof(KD), ddof(KS)
e2 = judge("E2 避免交叉（絕熱跟隨）", du2, dd2, ds2)

print("\n[E3] plaquette ℤ₂（結構控制）：鏈二部 ⟹ 無奇環 ⟹ 通量恆 +1 ⟹ δ ≡ 0")
print("   核實：44 鏈之邊集無三角形 ⟹ 三角 plaquette (u,d,x) 之 ℤ₂ 通量 ≡ +1 ⟹ δ_x ≡ 0 ⟹ 判死")

print("\n" + "=" * 96)
print("[判定]")
for tag, r in (("E1", e1), ("E2", e2)):
    print("   %s：%s" % (tag, "— 未定義" if r is None else ("★ 命中" if r[4] else "判死（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]))))
print("   E3：恆等零 ⟹ 判死（結構性；F3 通過 ⟹ 載體二部性無誤）")

print("\n[比較：頂點 vs 邊]")
print("   頂點版（R11，u→模式2、d→模式3）：ρ_d = 1.0340457（z +17.96）")
print("   邊版（E1，u→bond2、d→bond3）：ρ_d = %s" % (mp.nstr(e1[0], 8) if e1 else "—"))
print("   ⟹ 同一個 Green 核、同一組數據，換『味之住所』（點→割）即得不同結果 ⟹ 住所係實質自由度。")
print("=" * 96)