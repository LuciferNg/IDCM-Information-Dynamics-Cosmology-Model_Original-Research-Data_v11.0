#!/usr/bin/env python3
"""44 — R22：Ω 帳本之真結構式（180 分母／能量帳本）＋ 對 R21 之更正

使用者指正：`41/60 = 61.5/90 = 123/180`，且 `123 = 2⟨W⟩² = E_J`、`Σγ = 90`、`E_sync = 61.5` 皆已註冊。
⟹ R21 之「41 = 殘餘、60 = |A₅| 巧合」**錨錯對象**（應錨 90 = Σγ 與 123 = E_J，非 60 = |A₅|）。

本輪：
 A 更正：核實 41/60 = 61.5/90 = 123/180 為**精確有理恆等**（非巧合），並指出真錨
 B 完整帳本：180 分母之三項分解；及帳本鏈
 C 味比值：以帳本條目之兩兩比值測 ρ_d（N = 20，全報）
 D ρ_u 之顆粒 no-go：帳本條目之差之最小非零顆粒 vs |ρ_u|
 E 有理數密度對照（說明為何任何有理命中共皆不足以作證據）

判準（跑前凍死）
--------------------------------------------------
F1 A/B 用 `Fraction` 精確；任何「相等」須為精確有理恆等，否則不宣稱
F2 C 須全報 20 個比值；命中 = |z| < 3
F3 D：若最小非零顆粒 > |ρ_u| ⟹ 該帳本形式**結構上不可表示** ρ_u（不依賴有理數密度）
F4 E：任何有理吻合須同時報鄰域密度；否則不作證據
"""
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 30
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")

print("=" * 96)
print("44 — R22：Ω 帳本之真結構式 ＋ 對 R21 之更正")
print("=" * 96)

# ---------- A：更正 ----------
print("\n[A] 更正：41/60 = 61.5/90 = 123/180 係**精確有理恆等**")
b = Fraction(41, 60)
print("   41/60 = %s ；61.5/90 = %s ；123/180 = %s ⟹ 全部相等 = %s"
      % (Fraction(41, 60), Fraction(615, 10) / 90, Fraction(123, 180), b))
print("   已註冊量：Σγ = 90、E_sync = 61.5 = 123/2（= 88 − 4L，L = 53/8）、123 = 2⟨W⟩² = E_J")
print("   ⟹ **真錨 = Σγ = 90 與 E_J = 123**，非 60 = |A₅|。")
print("   ⟹ R21 §A「41 為殘餘」與 §B「60 = |A₅| 單點相交 ⟹ 判巧合」**撤回**：")
print("      41/60 = E_sync/Σγ 係**精確結構式**，41 有身份（= 123·60/180 = E_J/3 + ... 見 B）。")

# ---------- B：完整帳本 ----------
print("\n[B] 完整帳本（180 分母 = 2Σγ）")
print("   %-8s %-14s %-14s %s" % ("密度", "180 分母", "Σγ 單位", "算子身份"))
tab = [("Ω_b", 9, Fraction(9, 2), "N_gen² / 2 = 4.5"),
       ("Ω_DM", 48, Fraction(48, 2), "2·B·N_gen = 24"),
       ("Ω_Λ", 123, Fraction(123, 2), "E_sync = 123/2 = 61.5"),
       ("Ω_m", 57, Fraction(57, 2), "E_proj = 57/2 = 28.5")]
for name, num, e, ident in tab:
    print("   %-8s %-14s %-14s %s" % (name, "%d/180" % num, "%s/90" % e, ident))
tot = sum(t[1] for t in tab[:3])
print("   核實：9 + 48 + 123 = %d = 180 = 2Σγ ✓（Ω_b + Ω_DM + Ω_Λ = 1）" % tot)
print("         Ω_b + Ω_DM = 9 + 48 = 57 = Ω_m 之分子 ✓ ⟹ **E_b + E_DM = E_proj**（4.5 + 24 = 28.5）")
print("\n   帳本鏈（精確）：")
print("      E_sync − E_proj = 61.5 − 28.5 = %s = K_tree = 33 ✓" % (Fraction(615, 10) - Fraction(285, 10)))
print("      E_proj − E_DM  = 28.5 − 24   = %s = E_b = 4.5 ✓" % (Fraction(285, 10) - 24))
print("      E_b + E_DM = 28.5 = E_proj ✓ ； E_proj + E_sync = 90 = Σγ ✓")
print("   ★ 新結構陳述：Ω 家族 = **N_gen² : 2BN_gen : 2⟨W⟩²** / 2Σγ，且帳本滿足")
print("      `E_sync − E_proj = K_tree`、`E_proj − E_DM = E_b`。")

# ---------- C：味比值 vs 帳本兩兩比值 ----------
print("\n[C] 帳本條目之兩兩比值 vs ρ_d（N = 20，全報）")
Q = {"E_b": Fraction(9, 2), "E_DM": Fraction(24), "E_proj": Fraction(57, 2),
     "K_tree": Fraction(33), "E_sync": Fraction(123, 2)}
print("   帳本條目：%s" % {k: str(v) for k, v in Q.items()})
hits = []
rows = []
for k1 in Q:
    for k2 in Q:
        if k1 == k2:
            continue
        r = Q[k1] / Q[k2]
        pd = mp.mpf(r.numerator) / r.denominator
        z = (pd - RD) / RD_S
        rows.append((k1, k2, r, pd, z, abs(z) < 3))
        if abs(z) < 3:
            hits.append((k1, k2, r, pd, z))
for k1, k2, r, pd, z, ok in sorted(rows, key=lambda x: abs(x[4]))[:8]:
    print("   %-7s / %-7s = %-10s = %-12s z = %+7.2f %s"
          % (k1, k2, str(r), mp.nstr(pd, 8), z, "★" if ok else ""))
print("   ...（其餘 %d 個 |z| ≥ 3）" % (len(rows) - len(hits)))
print("   命中 = %d ／ 20" % len(hits))
for h in hits:
    print("     ★ %s / %s = %s ⟹ ρ_d = %s（z %+.2f）" % (h[0], h[1], str(h[2]), mp.nstr(h[3], 8), h[4]))

# ---------- D：ρ_u 之顆粒 no-go ----------
print("\n[D] ρ_u 之顆粒 no-go（不依賴有理數密度）")
print("   數據 ρ_u = %s ± %s ⟹ |ρ_u| = %s" % (mp.nstr(RU, 8), mp.nstr(RU_S, 6), mp.nstr(abs(RU), 8)))
vals = sorted(Q.values())
diffs = sorted({abs(a - b) for i, a in enumerate(vals) for b in vals[i + 1:]})
print("   帳本條目之非零差（升序）：%s" % [str(d) for d in diffs])
dmin = diffs[0]
print("   最小非零差 = %s ；除 Σγ = 90 ⟹ 最小非零比值顆粒 = %s" % (dmin, dmin / 90))
print("   需求 |ρ_u| = %s ⟹ 顆粒 / |ρ_u| = %s"
      % (mp.nstr(abs(RU), 8), mp.nstr((mp.mpf(dmin.numerator) / dmin.denominator / 90) / abs(RU), 6)))
print("   ⟹ **最小顆粒（%s）比需求（%s）大 %s 倍** ⟹ 任何「帳本條目之差／和 ÷ Σγ」之形式"
      % (dmin / 90, mp.nstr(abs(RU), 6), mp.nstr((mp.mpf(dmin.numerator) / dmin.denominator / 90) / abs(RU), 4)))
print("      **結構上不可能**表示 ρ_u（除差為零）⟹ F3 成立")

# ---------- E：有理數密度對照 ----------
print("\n[E] 有理數密度對照（為何任何有理吻合皆不足）")
print("   ρ_d：分母 ≤ 180 之最簡分數中最接近者")
best = None
for c in range(2, 181):
    for a in range(0, c + 1):
        f = Fraction(a, c)
        if f.denominator != c:
            continue
        z = (mp.mpf(a) / c - RD) / RD_S
        if best is None or abs(z) < abs(best[2]):
            best = (f, mp.mpf(a) / c, z)
print("      %s = %s ⟹ z = %+.5f" % (best[0], mp.nstr(best[1], 10), best[2]))
cnt = sum(1 for c in range(2, 181) for a in range(0, c + 1)
          if Fraction(a, c).denominator == c and abs(mp.mpf(a) / c - RD) <= RD_S)
tot = sum(1 for c in range(2, 181) for a in range(0, c + 1) if Fraction(a, c).denominator == c)
print("      分母 ≤ 180 之最簡分數共 %d 個，±1σ 內 %d 個（%.1f%%）" % (tot, cnt, 100.0 * cnt / tot))
print("   ⟹ 任何單一有理吻合若無結構來源，皆屬此 %.1f%% 之必然" % (100.0 * cnt / tot))

# ---------- F：精確換算（供文檔溯源） ----------
print("\n[F] 精確換算（本節值供文檔引用；以小數打印）")
print("   1/180（鬆形式顆粒）= %s ⟹ 文檔作 0.00556" % mp.nstr(mp.mpf(1) / 180, 20))
print("   4.5/57（比值之差之最小值）= %s ⟹ 文檔作 0.0789" % mp.nstr(mp.mpf("4.5") / 57, 20))
print("   4.5/90（差形式顆粒）= %s ⟹ 文檔作 0.05" % mp.nstr(mp.mpf("4.5") / 90, 20))
print("   |ρ_u| 之精確值 = %s（文檔作 0.02588837）" % mp.nstr(abs(RU), 20))

print("\n" + "=" * 96)
print("[判定]")
print("   A 更正：41/60 = E_sync/Σγ 係**精確結構式** ⟹ R21 §A（41 殘餘）與 §B（60 巧合）撤回")
print("   B 新結構：Ω 家族 = N_gen² : 2BN_gen : 2⟨W⟩² / 2Σγ；且 E_sync−E_proj = K_tree、E_proj−E_DM = E_b")
print("   C 味比值：ρ_d 之帳本候選命中 %d ／ 20（最佳 %s/%s，z = %+.2f）"
      % (len(hits), hits[0][0], hits[0][1], hits[0][4]) if hits else "   C 味比值：命中 0 ／ 20")
print("   D 顆粒 no-go：最小非零顆粒 = %s > |ρ_u| = %s ⟹ 帳本形式不可表示 ρ_u" % (dmin / 90, mp.nstr(abs(RU), 6)))
print("   E 有理吻合之密度 = %.1f%% ⟹ 不作證據" % (100.0 * cnt / tot))
print("=" * 96)