#!/usr/bin/env python3
"""R12b — R12 之結構旁證：四條精確陳述（無新協議，不加 N_trials）

本輪**唔測新形式**。只把 R12 觀測到嘅「三個載體都給 ρ_d ≈ 1.03」提升為結構陳述，
並用精確代數核實。全部 ρ 值在本腳本內由載體譜重算（不引用外部硬編碼），
可與 R11／R12 之 log 逐位對照。判準只有「等式是否精確成立」，無門檻可調。

S1  TI 核 G(i,j)=f(|i−j|) ⟹ δ_u = −δ_d 恆等 ⟹ ρ_u = −ρ_d。
S2  三個已測載體之 ρ_d 集中度（本腳本重算）。
S3  單模 PSD 核 f(n)=cos(nθ) 之閉式 ρ_d = 1/(1+2cosθ)、ρ_u = −ρ_d。
S4  曲率條件：需 (f1−f2)/(f0−f1) = 1/ρ_d = 2.4783；實得 1/ρ_d(載體)。
"""
import mpmath as mp

mp.mp.dps = 40

N44, N88 = 44, 88
g44 = lambda n: 4 * mp.sin(n * mp.pi / 88) ** 2
g88 = lambda k: 4 * mp.sin(k * mp.pi / 176) ** 2
ph44 = lambda n, j: mp.sqrt(mp.mpf(2) / (N44 + 1)) * mp.sin(n * mp.pi * j / (N44 + 1))
ph88 = lambda k, j: mp.sqrt(mp.mpf(2) / (N88 + 1)) * mp.sin(k * mp.pi * j / (N88 + 1))

G44 = lambda i, j: sum(ph44(n, i) * ph44(n, j) / g44(n) for n in range(1, N44 + 1))
G88 = lambda i, j: sum(ph88(k, i) * ph88(k, j) / g88(k) for k in range(1, N88 + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
DU, DD = mp.mpf("-0.0887"), mp.mpf("1.3820")

print("=" * 96)
print("R12b — R12 結構旁證：四條精確陳述（無新協議）")
print("=" * 96)

# ---------------- S1
print("\n[S1] 平移不變核之恆等 δ_u = −δ_d（δ_x = f(|k_d−x|) − f(|k_u−x|)）")
for nm, f in (("f(n)=e^{-0.3n}", lambda n: mp.e ** (-mp.mpf("0.3") * n)),
              ("f(n)=1/(1+n)^2", lambda n: 1 / (1 + n) ** 2),
              ("f(n)=cos(nπ/3)", lambda n: mp.cos(n * mp.pi / 3))):
    du, dd = f(1) - f(0), f(0) - f(1)
    print("   %-18s δ_u = %-24s δ_d = %-24s δ_u+δ_d = %s ⟹ %s"
          % (nm, mp.nstr(du, 8), mp.nstr(dd, 8), mp.nstr(du + dd, 3),
             "恆等 ✓" if abs(du + dd) < mp.mpf("1e-30") else "✗"))
ratio_obs = abs(DU) / DD
print("   實測：|δ_u|/δ_d = %s（δ_u=%s、δ_d=%s）" % (mp.nstr(ratio_obs, 8), DU, DD))
print("   TI 預言此比值 = 1 ⟹ 實測偏 %s 倍 ⟹ **TI 全族判死（與核形狀、載體大小無關）**"
      % mp.nstr(1 / ratio_obs, 6))

# ---------------- S2
print("\n[S2] 三個載體之 ρ_d（本腳本由譜重算，與 R11／R12 log 對照）")
cfgs = [
    ("44-site Dirichlet（R11 P1）", G44, 2, 3, 4),
    ("88-site carry-by-2（R12 A1）", G88, 4, 6, 8),
    ("88-site direct（R12 A2）", G88, 2, 3, 4),
]
res = []
for nm, G, ku, kd, ks in cfgs:
    du = G(kd, ku) - G(ku, ku)
    dd = G(kd, kd) - G(ku, kd)
    ds = G(kd, ks) - G(ku, ks)
    pd, pu = dd / ds, du / ds
    res.append((nm, pd, pu))
    print("   %-30s ρ_d = %-12s（z %+6.2f）  ρ_u = %-12s（z %+6.2f）"
          % (nm, mp.nstr(pd, 8), (pd - RD) / RD_S, mp.nstr(pu, 8), (pu - RU) / RU_S))
pds = [r[1] for r in res]
spread = (max(pds) - min(pds)) / max(pds)
print("   三者極差 = %s%% ⟹ **載體無關**（44/88、兩種 carry 皆同）" % mp.nstr(spread * 100, 4))
print("   對目標之偏離 = %s–%s 倍 ⟹ 失效係**族性質**，非載體選擇問題"
      % (mp.nstr(min(pds) / RD, 5), mp.nstr(max(pds) / RD, 5)))

# ---------------- S3
print("\n[S3] 單模 PSD 核 f(n) = cos(nθ)（Toeplitz 正定，屬 TI 族）")
print("     核實恆等 ρ_d = 1/(1+2cosθ)：")
mx = mp.mpf(0)
for k in (1, 7, 22, 44, 88, 120, 170):
    th = mp.pi * k / 176
    f = lambda n, th=th: mp.cos(n * th)
    lhs = (f(0) - f(1)) / (f(1) - f(2))
    rhs = 1 / (1 + 2 * mp.cos(th))
    mx = max(mx, abs(lhs - rhs))
    print("     θ=π·%3d/176：ρ_d = %-26s 1/(1+2cosθ) = %-26s |Δ| = %s"
          % (k, mp.nstr(lhs, 10), mp.nstr(rhs, 10), mp.nstr(abs(lhs - rhs), 3)))
print("     最大偏差 = %s ⟹ 閉式%s" % (mp.nstr(mx, 3),
                                    "精確 ✓" if mx < mp.mpf("1e-30") else "**不成立**"))
print("     可達範圍：θ→0 ⟹ ρ_d→1/3；θ→2π/3⁻ ⟹ +∞ ⟹ ρ_d ∈ (1/3, ∞)")
th_star = mp.acos((1 / RD - 1) / 2)
print("     反解 ρ_d = %s ⟹ cosθ = %s ⟹ θ = %s°；同核同時給 ρ_u = −%s（仍被 S1 判死）"
      % (mp.nstr(RD, 8), mp.nstr(mp.cos(th_star), 8),
         mp.nstr(th_star * 180 / mp.pi, 8), mp.nstr(RD, 8)))

# ---------------- S4
print("\n[S4] 曲率條件：ρ_d ≈ [f(0)−f(1)]/[f(1)−f(2)] ⟹ 需 (f1−f2)/(f0−f1) = 1/ρ_d")
print("     目標 1/ρ_d = %s" % mp.nstr(1 / RD, 8))
for nm, pd, pu in res:
    print("     %-30s 1/ρ_d = %-10s ⟹ 距目標差 %s 倍"
          % (nm, mp.nstr(1 / pd, 6), mp.nstr((1 / RD) / (1 / pd), 5)))
print("     另兩核（TI 參考）：")
for nm, f in (("e^{-0.3n}", lambda n: mp.e ** (-mp.mpf("0.3") * n)),
              ("e^{-n²/4}", lambda n: mp.e ** (-mp.mpf(n) ** 2 / 4))):
    r = (f(1) - f(2)) / (f(0) - f(1))
    print("     %-30s 1/ρ_d = %-10s ⟹ ρ_d = %s" % (nm, mp.nstr(r, 6), mp.nstr(1 / r, 6)))
print("     三載體（Green 族）實得 0.967–0.984、指數 0.741、Gauss 1.858 ⟹ 皆 < 2.478")
print("     ⟹ 達標需 n=1→2 有「崖」；單模 PSD 可達（S3 反解 θ=42.34°）但被 S1 判死")
print("     ⟹ 存活須**非平移不變**且具 ρ_u/ρ_d ≈ 0.064 之 15.6× 不對稱（帶號性非由 S4 推出）")

print("\n" + "=" * 96)
print("[結論 — 四條皆精確陳述，無判準可調]")
print("   S1  TI ⟹ ρ_u = −ρ_d 恆等；實測比 0.0640 ⟹ TI 全族判死（15.6×）")
print("   S2  三載體 ρ_d 極差 %s%% ⟹ 載體無關" % mp.nstr(spread * 100, 4))
print("   S3  單模 PSD 閉式 ρ_d = 1/(1+2cosθ) 精確；ρ_u = −ρ_d 同被 S1 殺")
print("   S4  需 2.4783；Green 族實得 %s–%s（差 2.5×）、指數 0.741、Gauss 1.858 ⟹ Green 短程結構不足"
      % (mp.nstr(min(1 / p for p in pds), 5), mp.nstr(max(1 / p for p in pds), 5)))
print("   ⟹ 非由 S4 推出「帶號」；S1 已獨立判死全部 TI（含 PSD）⟹ 存活須非 TI ＋ 15.6× 不對稱")
print("=" * 96)