#!/usr/bin/env python3
"""41 — R20：球面三三角形提案之結構判定（Euler／符號／均勻性／p）

使用者提案（2026-10-29）：「網是三個三角形組成一個球面三維立體網，而且這三個三角形內角和
為了達成球面投影並不是 180 度，還帶有一個 p 的阻挫來提供味」。

本輪逐項判定（判準跑前凍死）
--------------------------------------------------
A Euler：F 個三角面閉合球面之可行性（V − E + F = 2、3F = 2E）
B 符號閘：球面 excess E_i ≥ 0（Gauss-Bonnet）；若 δ ∝ E ⟹ ρ_u = δ_u/δ_s ≥ 0 ⟹ 與數據比對
C 均勻性：若對稱群在面上遞移 ⟹ 所有 E_i 相等 ⟹ δ ≡ 0
D 最小非對稱（等差）模型之零參數關係：E ∝ (1−p, 1, 1+p) ⟹ ρ_d = (1+ρ_u)/2
E 已註冊三元組之置換窮舉（N_trials = 18）：(3,16,41)、(5,8,13)、(11,18,29)
F p 之反解與已註冊量之比對

判準：|z| < 3 為過；兩觀測量皆過才算命中。E 之 18 組合須全報（選擇效應）。
"""
import itertools

import mpmath as mp

mp.mp.dps = 30
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")

print("=" * 96)
print("41 — R20：球面三三角形提案（excess ＋ p 阻挫）之結構判定")
print("=" * 96)

# ---------- A：Euler ----------
print("\n[A] Euler 可行性（V − E + F = 2，三角面 ⟹ 3F = 2E ⟹ V = 2 + F/2）")
for F in (2, 3, 4, 6, 8, 20):
    E = mp.mpf(3 * F) / 2
    V = 2 + mp.mpf(F) / 2
    ok = (E == int(E)) and (V == int(V))
    print("   F = %-2d ⟹ E = %-5s、V = %-5s ⟹ %s"
          % (F, mp.nstr(E, 4), mp.nstr(V, 4), "可閉合" if ok else "**不可能（非整數）**"))
print("   ⟹ **三個三角形閉合球面不可能**（3F=2E ⟹ E=4.5）。")
print("      最小修補：① 四面體（F=4，每面 excess = 4π/4 = π）")
print("                ② 三個 lune（二角形；V=2、E=3、F=3 ⟹ Euler 2 ✓，每個 lune 角 = 2π/3）")

# ---------- B：符號閘 ----------
print("\n[B] 符號閘：球面 excess ≥ 0（Gauss-Bonnet：Σ E_i = 4π > 0）")
print("   若 δ ∝ E（每味一面），則 δ_x ≥ 0 ∀x ⟹ ρ_u = δ_u/δ_s ≥ 0")
print("   數據：ρ_u = %s ± %s ⟹ 距 0 = %s σ（弱，但號為負）"
      % (mp.nstr(RU, 8), mp.nstr(RU_S, 6), mp.nstr(abs(RU) / RU_S, 4)))
print("   反解所需 excess 比：E_d/E_s = %s、E_u/E_s = %s" % (mp.nstr(RD, 8), mp.nstr(RU, 8)))
print("   ⟹ **E_u/E_s < 0 ⟹ 非負 excess 不可能** ⟹ `δ ∝ excess` 之字面讀法**判死**")
print("   ⟹ 唯一存活讀法：δ_x = C(E_{x;d} − E_{x;u})（兩個 excess 之差，帶號）—— 但差值形式")
print("      無 Gauss-Bonnet 約束 ⟹ 六個 excess 只有一個和約束 ⟹ 欠定（同 docs/28 之結論）")

# ---------- C：均勻性 ----------
print("\n[C] 均勻性定理")
print("   若對稱群 G 在面上遞移 ⟹ 所有 E_i = 4π/F 相等 ⟹ 任何面內差或面間差 ≡ 0")
print("   四面體（T_d 遞移於 4 面）、二十面體（I 遞移於 20 面）、三個 lune（D₃ 遞移）皆遞移")
print("   ⟹ **對稱版本必死**（與 docs/29 之 uniform-flux 同一機制）")

# ---------- D：等差模型之零參數關係 ----------
print("\n[D] 最小非對稱（等差）模型：E ∝ (1−p, 1, 1+p)  [u, d, s]；Σ ∝ 3 = 常數 ✓")
print("   ⟹ ρ_u = (1−p)/(1+p)、ρ_d = 1/(1+p) ⟹ **消去 p 得零參數關係 ρ_d = (1+ρ_u)/2**")
pred = (1 + RU) / 2
z = (pred - RD) / RD_S
print("   預言 ρ_d = %s（由實測 ρ_u 反推） vs 實測 %s ⟹ z = %+.2f" % (mp.nstr(pred, 8), mp.nstr(RD, 8), z))
print("   ⟹ %s（門檻 |z| < 3）" % ("**過閘**（但邊緣）" if abs(z) < 3 else "判死"))

# ---------- E：已註冊三元組之置換 ----------
print("\n[E] 已註冊三元組之置換窮舉（N_trials = 18，全報）")
triples = {"(3,16,41) Ω 家族": (3, 16, 41), "(5,8,13) Fibonacci": (5, 8, 13), "(11,18,29) Lucas": (11, 18, 29)}
hits = []
for name, trip in triples.items():
    print("   %s：" % name)
    for perm in itertools.permutations(trip):
        au, ad, as_ = perm
        if as_ == 0:
            continue
        pd, pu = mp.mpf(ad) / as_, mp.mpf(au) / as_
        zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
        ok = abs(zd) < 3 and abs(zu) < 3
        if ok:
            hits.append((name, perm, pd, pu, zd, zu))
        print("      (u,d,s)=%-12s ρ_d=%8s（z %+7.2f）  ρ_u=%8s（z %+7.2f） %s"
              % (str(perm), mp.nstr(pd, 6), zd, mp.nstr(pu, 6), zu, "★" if ok else ""))
print("   命中（兩觀測量皆 |z|<3）= %d ／ 18" % len(hits))
for h in hits:
    print("     ★ %s (u,d,s)=%s ⟹ ρ_d = %s（z %+.2f）、ρ_u = %s（z %+.2f）"
          % (h[0], h[1], mp.nstr(h[2], 8), h[4], mp.nstr(h[3], 8), h[5]))

# ---------- F：p 之反解 ----------
print("\n[F] p 之反解與已註冊量比對（等差模型）")
p_star = (1 - RU) / (1 + RU)
print("   由 ρ_u 反解 p = (1−ρ_u)/(1+ρ_u) = %s" % mp.nstr(p_star, 10))
cands = {"32/120 = 4/15（RH 阻挫比）": mp.mpf(4) / 15, "1/5": mp.mpf(1) / 5, "3/19（Ω_b/Ω_m）": mp.mpf(3) / 19,
         "19/60（Ω_m）": mp.mpf(19) / 60, "41/60（Ω_Λ）": mp.mpf(41) / 60, "1/φ = 0.618": 1 / mp.phi}
for k, v in cands.items():
    print("   %-24s = %-10s   差 = %s" % (k, mp.nstr(v, 8), mp.nstr(abs(v - p_star), 6)))
print("   ⟹ 無任何已註冊量匹配（最近者差 %s）⟹ p 無結構來歷"
      % mp.nstr(min(abs(v - p_star) for v in cands.values()), 6))

print("\n" + "=" * 96)
print("[判定]")
print("   A 三三角形閉合：**Euler 不可能**（E=4.5）⟹ 須改四面體（F=4）或三 lune（F=3）")
print("   B δ ∝ excess：**判死**（數據要求 E_u < 0，而球面 excess ≥ 0）")
print("   C 對稱版本：**必死**（遞移 ⟹ 全等 ⟹ δ ≡ 0）")
print("   D 等差模型：零參數關係 ρ_d = (1+ρ_u)/2 ⟹ 預言 %s vs 實測 %s ⟹ z = %+.2f（過閘但邊緣）"
      % (mp.nstr(pred, 7), mp.nstr(RD, 7), z))
print("   E 已註冊三元組置換：命中 %d ／ 18 ⟹ %s" % (len(hits), "有命中，須報選擇效應" if hits else "零命中"))
print("   F p 反解 = %s ⟹ 與已註冊量皆不匹配 ⟹ 無結構來歷" % mp.nstr(p_star, 8))
print("   ⟹ 提案之唯一存活形式 = **兩個 excess 之差**（欠定），或**等差模型+邊緣命中**（z=+2.38σ，")
print("      1 dof，未達 5σ；且 p 無結構來歷）")
print("=" * 96)