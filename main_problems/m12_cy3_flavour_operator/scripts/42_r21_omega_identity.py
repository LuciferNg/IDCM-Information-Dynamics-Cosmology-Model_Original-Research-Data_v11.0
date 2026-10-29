#!/usr/bin/env python3
"""42 — R21：(3,16,41) 之算子身份、B-擾動測試、有理數密度 null test

背景：docs/33 §E 發現 (3,16,41) 配 (u,d,s)=(3,16,41) 時
      ρ_d = 16/41 = 0.390244（z −0.38，極準）、ρ_u = 3/41 = 0.0731707（z +3.01，號錯）。
      本輪推三件事：
      A 算子身份：3、16、41、60 各係咩？
      B B-擾動：60 = B(B²−1) 係結構定巧合？（skill pitfall 5/6 之正式方法）
      C look-elsewhere 修正：單調指派原則 ⟹ N_trials 由 18 降為 3
      D 有理數密度 null test：0.4σ 命中係否稀奇？

判準（跑前凍死）
--------------------------------------------------
F1 A 部分：只接受「已註冊量之直接算子身份」；殘餘者標 residual，不宣稱結構
F2 B 部分：若 B(B²−1) 對 B 無一般恆等（只在 B=4 與 |A₅| 相交）⟹ 判巧合
F3 C 部分：單調指派 = excess 隨味世代單調 ⟹ 唯一指派，N_trials = 3（三個已註冊三元組）
F4 D 部分：若「±1σ 內之小分母有理數」甚多 ⟹ 該命中無資訊（同 R10 之教訓：數值只有否證權）
"""
import itertools
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 30
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")

print("=" * 96)
print("42 — R21：(3,16,41) 算子身份 ＋ B-擾動 ＋ 有理數密度 null test")
print("=" * 96)

# ---------- A：算子身份 ----------
print("\n[A] 算子身份（B = 4 = 分支數；N_gen = 3 = B−1；N_c = 3）")
B = 4
items = [
    ("3", Fraction(3), "N_gen = N_c = B−1"),
    ("16", Fraction(16), "B² = 16"),
    ("60", Fraction(60), "B(B²−1) = 4·15 = 60 = |A₅| = V(RH)"),
    ("41", Fraction(41), "B(B²−1) − B² − (B−1) = 60 − 16 − 3（殘餘 residual）"),
    ("19", Fraction(19), "Ω_m 之分子 = 60 − 41 = B² + (B−1) = 16 + 3"),
]
for name, val, ident in items:
    print("   %-4s = %-6s  %s" % (name, val, ident))
print("   核實：3 + 16 + 41 = %d ✓ ；19 + 41 = %d = B(B²−1) ✓" % (3 + 16 + 41, 19 + 41))
print("   ⟹ 三個分子中 **3 與 16 有直接算子身份**；41 為殘餘（60 − 19），其身份依附於 60 之身份")
print("   ⟹ 60 之身份（|A₅|／V(RH)／B(B²−1)）係關鍵：若 60 巧合，則 41 亦無結構")

# ---------- B：B-擾動 ----------
print("\n[B] B-擾動測試：B(B²−1) 有無一般恆等？（B = 2..7）")
print("   %-3s %-10s %-12s %-14s %s" % ("B", "B(B²−1)", "(B−1):B²:殘餘", "Σ", "備註"))
for b in range(2, 8):
    tot = b * (b * b - 1)
    rem = tot - b * b - (b - 1)
    note = ""
    if b == 4:
        note = "|A₅| = 60（B=4 專有）"
    elif b == 3:
        note = "|S₄| = 24"
    elif b == 5:
        note = "|S₅| = 120"
    print("   %-3d %-10d %-12s %-14d %s" % (b, tot, "%d:%d:%d" % (b - 1, b * b, rem), b * b + (b - 1) + rem, note))
print("   ⟹ `B(B²−1)` 對所有 B 皆成立（代數恆等），**但** `= |A₅|` 只在 B=4；")
print("      `B(B²−1) = |S_{B+1}|` 只在 B=2、3、5（B=4 時為 |S₅|/2 = |A₅|）。")
print("   ⟹ 依 F2：`60 = |A₅|` 屬 **B=4 之單點相交** ⟹ 判巧合（不得作為結構依據）")

# ---------- C：單調指派 ----------
print("\n[C] 單調指派原則：excess 隨味世代（u<d<s）單調 ⟹ 唯一指派")
triples = {"(3,16,41)": (3, 16, 41), "(5,8,13)": (5, 8, 13), "(11,18,29)": (11, 18, 29)}
print("   %-12s %-16s %-22s %-22s %s" % ("三元組", "單調指派", "ρ_d（z）", "ρ_u（z）", "判"))
for name, trip in triples.items():
    au, ad, as_ = sorted(trip)                      # 單調：u 最小、s 最大
    pd, pu = mp.mpf(ad) / as_, mp.mpf(au) / as_
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = abs(zd) < 3 and abs(zu) < 3
    print("   %-12s %-16s %8s（%+7.2f）      %8s（%+7.2f）      %s"
          % (name, str((au, ad, as_)), mp.nstr(pd, 6), zd, mp.nstr(pu, 6), zu, "★" if ok else "判死"))
print("   ⟹ N_trials 由 18 降為 3（F3），命中 0/3；惟 (3,16,41) 之 ρ_d 達 z = −0.38")

# ---------- D：有理數密度 null test ----------
print("\n[D] 有理數密度 null test：ρ_d 之 0.4σ 命中係否稀奇？")
print("   目標 ρ_d = %s ± %s" % (mp.nstr(RD, 10), mp.nstr(RD_S, 6)))
for CMAX in (41, 60, 100):
    fracs = []
    for c in range(2, CMAX + 1):
        for a in range(1, c):
            if Fraction(a, c).denominator == c:
                fracs.append(Fraction(a, c))
    in1 = [f for f in fracs if abs(mp.mpf(f.numerator) / f.denominator - RD) <= RD_S]
    in_half = [f for f in fracs if abs(mp.mpf(f.numerator) / f.denominator - RD) <= RD_S / 2]
    best = min(fracs, key=lambda f: abs(mp.mpf(f.numerator) / f.denominator - RD))
    zb = (mp.mpf(best.numerator) / best.denominator - RD) / RD_S
    print("   分母 ≤ %-3d：最簡分數共 %5d 個；±1σ 內 %4d 個（%.1f%%）；±0.5σ 內 %4d 個"
          % (CMAX, len(fracs), len(in1), 100.0 * len(in1) / len(fracs), len(in_half)))
    print("                   最接近者 = %s = %s ⟹ z = %+.4f" % (best, mp.nstr(mp.mpf(best.numerator) / best.denominator, 10), zb))
print("   ⟹ 0.4035 之鄰域有大量小分母有理數 ⟹ **16/41 之 0.4σ 命中毫不稀奇**")

print("\n[D2] 聯合命中（ρ_d 與 ρ_u 同時）之機率上界")
print("   ρ_u 之目標 %s ± %s（負值）⟹ 任何正分子三元組給 ρ_u > 0 ⟹ **結構上不可能**" % (mp.nstr(RU, 8), mp.nstr(RU_S, 6)))
print("   ⟹ 分母為正整數之 a/c 族**全部**落在 ρ_u > 0 ⟹ 該族對 ρ_u 之命中率 = 0")

print("\n" + "=" * 96)
print("[判定]")
print("   A：3 = N_gen、16 = B² 有直接身份；41 = 殘餘（依附 60）")
print("   B：60 = B(B²−1) 係代數恆等，但 `= |A₅|` 於 B=4 單點相交 ⟹ 判巧合（F2）")
print("   C：單調指派下 N_trials = 3，命中 0/3；(3,16,41) 之 ρ_d z = −0.38 但 ρ_u 差 3.01σ")
print("   D：±1σ 內小分母有理數甚多 ⟹ 單一觀測量之 0.4σ 命中無資訊（F4）")
print("   D2：ρ_u > 0 對正整數三元組係結構必然 ⟹ (3,16,41) 族**不可能**同時命中")
print("   ⟹ (3,16,41) 之近乎命中**不成立為證據**；(3:16:41)/60 之算子結構維持『41 為殘餘』")
print("=" * 96)