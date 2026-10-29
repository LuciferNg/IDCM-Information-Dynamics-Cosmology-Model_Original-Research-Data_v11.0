#!/usr/bin/env python3
"""43 — R22：ℤ₄ 扇區核（已註冊之離散帶號結構）—— 「帶號譜量 ＋ 結構規則」之最後候選

四輪匯合（docs/28、29、31、34）得出唯一未封之形式：
    **離散、已註冊、不可分離、非遞移、本質帶號** 之結構。
v11 內唯一符合「離散＋已註冊＋帶號」者 = **ℤ₄ 橋之 mod-4 扇區**（O210 charter，`176 = B·N_d = 4·44`）：
    `Σγ^(c) = N_d + (μ², −sec(π/88), 0, +sec(π/88))`  ⟹ 扇區 **帶號**，且**離散**、**已註冊**
    （`docs/28` 之 R15 已測 `w ≥ 0` 之一族；`docs/25` 之 S3a 只測 `(−1)^n`；**mod-4 扇區從未測**）

預註冊（跑前凍死；N_trials = 6）
--------------------------------------------------
 S1  w = 1_{n≡0 mod 4}/γ_n          S2  w = 1_{n≡1}/γ_n
 S3  w = 1_{n≡2}/γ_n                S4  w = 1_{n≡3}/γ_n
 S5  w = (1_{0} − 1_{2})/γ_n        S6  w = (1_{1} − 1_{3})/γ_n
 （S5 = 自共軛扇區對之差；S6 = ±sec 對之差 —— 二者於 charter 內有明確符號結構）
 控制 C1 = (−1)^n/γ_n（= S3a 之整數週期版，已知判死 z +9.87／+15.97）
 控制 C2 = 1/γ_n（R17 之 f0，已知判死 z +17.69）

判準：|z(ρ_d)| < 3 且 |z(ρ_u)| < 3 ⟹ 命中；全報，不得只報命中者。
"""
import mpmath as mp

mp.mp.dps = 30
NS = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2   # 44 站 Dirichlet（/90，一致約定）
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

print("=" * 96)
print("43 — R22：ℤ₄ 扇區核（N_trials = 6 ＋ 2 控制）")
print("=" * 96)
S = {c: [n for n in range(1, NS + 1) if n % 4 == c] for c in range(4)}
print("   扇區：c=0 → %d 個、c=1 → %d、c=2 → %d、c=3 → %d（44 = 4×11）"
      % tuple(len(S[c]) for c in range(4)))
print("   Σγ 之扇區和（/90 約定，非 176-L 之精確值，僅作自洽核對）：")
for c in range(4):
    print("      c=%d：Σγ = %s" % (c, mp.nstr(sum(gam(n) for n in S[c]), 10)))
print("   （charter 之精確式 `N_d + (μ², −sec, 0, +sec)` 屬 176-L；本輪用 44 站 /90，結構同型不同值）")


def run(tag, w, note=""):
    G = lambda i, j: mp.fsum(w(n) * phi(n, i) * phi(n, j) for n in range(1, NS + 1))
    ds = G(KD, KS) - G(KU, KS)
    if ds == 0:
        print("   %-26s δ_s = 0 ⟹ 未定義 ⟹ 判死 %s" % (tag, note))
        return None
    pd, pu = (G(KD, KD) - G(KU, KD)) / ds, (G(KD, KU) - G(KU, KU)) / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-26s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


ind = lambda c: (lambda n: 1 if n % 4 == c else 0)
res = {}
print("\n[六個 ℤ₄ 扇區核]")
res["S1"] = run("S1 1_{0}/γ", lambda n: mp.mpf(ind(0)(n)) / gam(n), "（扇區 0）")
res["S2"] = run("S2 1_{1}/γ", lambda n: mp.mpf(ind(1)(n)) / gam(n), "（扇區 1）")
res["S3"] = run("S3 1_{2}/γ", lambda n: mp.mpf(ind(2)(n)) / gam(n), "（扇區 2，Σ 最小）")
res["S4"] = run("S4 1_{3}/γ", lambda n: mp.mpf(ind(3)(n)) / gam(n), "（扇區 3）")
res["S5"] = run("S5 (1_0 − 1_2)/γ", lambda n: (mp.mpf(ind(0)(n)) - mp.mpf(ind(2)(n))) / gam(n), "（自共軛對差）")
res["S6"] = run("S6 (1_1 − 1_3)/γ", lambda n: (mp.mpf(ind(1)(n)) - mp.mpf(ind(3)(n))) / gam(n), "（±sec 對差）")

print("\n[控制]")
c1 = run("C1 (−1)^n/γ", lambda n: (mp.mpf(-1) ** n) / gam(n), "（= S3a，已知判死）")
c2 = run("C2 1/γ", lambda n: 1 / gam(n), "（= R17 f0，已知判死）")

print("\n" + "=" * 96)
print("[判定]")
hit = [k for k in ("S1", "S2", "S3", "S4", "S5", "S6") if res[k] and res[k][4]]
for k in ("S1", "S2", "S3", "S4", "S5", "S6"):
    r = res[k]
    print("   %s：%s" % (k, "— 未定義" if r is None else ("★ 命中（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]) if r[4] else "判死（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]))))
print("   命中 = %d ／ 6" % len(hit))
print("   控制 C1：%s ；C2：%s"
      % ("判死 ✓" if not (c1 and c1[4]) else "★ 異常", "判死 ✓" if not (c2 and c2[4]) else "★ 異常"))
print("\n[註] 若六個扇區全判死 ⟹ v11 內『離散＋已註冊＋帶號』之結構已窮盡；")
print("     此時 M12 之唯一剩餘形式 = 需**新註冊**帶號算子（非本 O 可造）。")
print("=" * 96)