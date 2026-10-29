#!/usr/bin/env python3
"""R9-附 — 後驗核快篩（POST-HOC，無預註冊判準，**不作判據**）

用途：只為令 `final_result.md §殘餘` 內兩個「後驗快篩」數字可溯源。
本腳本**不**構成任何檢定：形式係在看到目標比之後才提出，依禁區規則（後驗指數／
φ-power 同族）不得作為證據。判準欄留空係刻意。

目標（由 R7 之錨解出，scripts/20）：
  ρ_d = δ_d/δ_s = +0.40350365 ；ρ_b = +0.55279820 ；ρ_u = −0.02588837

形式（皆由 `u→2, d→3` 已註冊指數直接使用；s 取三讀 4/5/7）：
  F1 平移不變核：a(k₁,k₂) = γ_{|k₁−k₂|}          ⟹ δ_x = γ_{|3−k_x|} − γ_{|2−k_x|}
  F2 二分奇核  ：a(k₁,k₂) = (−1)^{k₁+k₂} γ_{|k₁−k₂|}
  γ_n = 4 sin²(πn/88)
"""
import mpmath as mp

mp.mp.dps = 30
g = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
RD, RB, RU = mp.mpf("0.40350365"), mp.mpf("0.55279820"), mp.mpf("-0.02588837")
KU, KD, KB = 2, 3, None  # b 無註冊指數 ⟹ ρ_b 不可評

print("=" * 88)
print("R9-附 — 後驗核快篩（POST-HOC，無判準，不作判據）")
print("=" * 88)
print("  γ_0..γ_5 =", [mp.nstr(g(n), 8) for n in range(6)])
print("  目標：ρ_d = %s  ρ_b = %s  ρ_u = %s" % (mp.nstr(RD, 8), mp.nstr(RB, 8), mp.nstr(RU, 8)))
print("  ⚠ ρ_b 不可評：b 之模式指數在 v11 未註冊（u→2、d→3 已裁）")

def kern(k1, k2, form):
    if form == "F1":
        return g(abs(k1 - k2))
    return (-1) ** (k1 + k2) * g(abs(k1 - k2))

for form in ("F1", "F2"):
    print("\n── %s ──" % ("F1 平移不變核 γ_|Δk|" if form == "F1" else "F2 二分奇核 (−1)^{k1+k2}γ_|Δk|"))
    for ks in (4, 5, 7):
        du = kern(KD, KU, form) - kern(KU, KU, form)
        dd = kern(KD, KD, form) - kern(KU, KD, form)
        ds = kern(KD, ks, form) - kern(KU, ks, form)
        ru = du / ds
        rd = dd / ds
        print("  k_s=%d: δ_u=%+.8f δ_d=%+.8f δ_s=%+.8f" % (ks, du, dd, ds))
        print("         ρ_u = %+.5f（目標 %+.5f，值/目標 = %6.2f）   ρ_d = %+.5f（目標 %+.5f，值/目標 = %5.2f）"
              % (ru, RU, abs(ru / RU) if RU != 0 else mp.inf,
                 rd, RD, abs(rd / RD) if RD != 0 else mp.inf))

print("\n" + "=" * 88)
print("[結論 — 後驗，僅供記錄]")
print("  F1：ρ_u 值/目標 = 12.94／7.76／4.37（全部大 4–13×）；ρ_d 值/目標 = 0.83／0.50／0.28（全部細 1.2–3.6×）⟹ 不合。")
print("  F2：ρ_d 值/目標 = 0.50／0.19／0.06（細 1.9–16×）；ρ_u 只有 k_s=7 吻合（0.95），")
print("      但同一 k_s 之 ρ_d 僅為目標 6% ⟹ 單點吻合不構成證據，整體不合。")
print("  兩形式皆屬 sign-definite 或二分奇之簡單核；真正存活者係 S2（δ 之 isovector 算子身份）。")
print("=" * 88)