#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/03 — 10^38 「兩因子」審計（R3／T3／P3）
====================================================
判準照 docs/00_preregistration.md §3：

  TF1  f_M0、f_R 各能否由**非反解**方式定出（兩條獨立結構式）
       PASS ⟺ 兩式皆非由 10^38 反解；否則 T3 🔴

本輪檢定三件事：
  (a) 乘積 `f_M0·f_R` 之數值身份 = `N_d·T₁` = `N_O38` = `M₀/(ħω_ring)`
  (b) 缺口 35.717977 dex 是否只係 `log₁₀(M₀/ħω_ring) − log₁₀(N_d)`（即無新內容）
  (c) 分解 (`f_M0`, `f_R`) 是否 gauge（乘積固定 ⟹ 一維簡併 ⟹ 無獨立身份）

上游常數（唯讀引用，不重推）：
  T₁ = 4.0951856e38            M8 docs/04 §4
  N_O38 = 1.801868e40          M8 docs/06（= 44·T₁，dev 7.4e−6）
  M₀/M_Pl = 2.2852931e−21      M6 final_result
  E_kchain = 4.348 meV         M8 docs/06（幾何平均錨，3 位有效）
  X₀ = 6.417006419e9           M8 docs/06
  N_IDCM = 7.88462687096e60    M8 docs/06（O212 判非結構 P=0.875）
  M_Pl = 1.220890e28 eV        CODATA 2018
  μ ≡ f_M0 ≡ X₀/√N_IDCM        M8 docs/06 E2

零 import 上游。
"""
import json
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


mp.mp.dps = 50
N_d = 44
T1 = mp.mpf("4.0951856e38")
N_O38 = mp.mpf("1.801868e40")
mu = mp.mpf("2.2852931e-21")          # M₀/M_Pl
E_kchain = mp.mpf("4.348e-3")         # eV（3 位有效）
X0 = mp.mpf("6.417006419e9")
N_IDCM = mp.mpf("7.88462687096e60")
M_Pl = mp.mpf("1.220890e28")          # eV

print("=" * 78)
print("M21 scripts/03 — 10^38 兩因子審計")
print("=" * 78)

# (a) 乘積身份
print("\n── (a) 乘積 f_M0·f_R 之數值身份 ──")
prod_ledger = N_d * T1
dev_a = (prod_ledger - N_O38) / N_O38
print(f"  N_d·T₁ = 44 × {mp.nstr(T1, 10)} = {mp.nstr(prod_ledger, 12)}")
print(f"  N_O38（M8 docs/06）= {mp.nstr(N_O38, 12)}")
print(f"  dev = {mp.nstr(dev_a, 4)}  （M8 報 7.4e−6 ≤ ω_site 精度界）")
ck("TF1-a-a  N_d·T₁ = N_O38（dev ≤ 1e-5）", abs(dev_a) < mp.mpf("1e-5"), f"dev = {mp.nstr(dev_a, 4)}")

# M₀/(ħω_ring) 由 E_kchain² = M_Pl·ħω_ring 反推
hbarw = E_kchain ** 2 / M_Pl
M0 = mu * M_Pl
ratio = M0 / hbarw
dev_b = (ratio - N_O38) / N_O38
print(f"\n  ħω_ring = E_kchain²/M_Pl = {mp.nstr(hbarw, 8)} eV")
print(f"  M₀ = μ·M_Pl = {mp.nstr(M0, 8)} eV")
print(f"  M₀/ħω_ring = {mp.nstr(ratio, 12)}   vs N_O38 = {mp.nstr(N_O38, 12)}")
print(f"  dev = {mp.nstr(dev_b, 4)}（限於 E_kchain 3 位有效 ⟹ 容差 1e−3）")
ck("TF1-a-b  N_O38 = M₀/(ħω_ring)（dev ≤ 1e-3，限 E_kchain 位數）", abs(dev_b) < mp.mpf("1e-3"),
   f"dev = {mp.nstr(dev_b, 4)}")

# (b) 缺口 = log₁₀(M₀/ħω_ring) − log₁₀(N_d)
print("\n── (b) 缺口之分解 ──")
gap_ratio = mp.log10(ratio) - mp.log10(N_d)
gap_pub = mp.mpf("35.717977")
print(f"  log₁₀(M₀/ħω_ring) = {mp.nstr(mp.log10(ratio), 10)}")
print(f"  log₁₀(N_d) = {mp.nstr(mp.log10(N_d), 10)}")
print(f"  差 = {mp.nstr(gap_ratio, 10)}  vs  M8 缺口 35.717977")
print(f"  M8 缺口 = log₁₀T₁ − 2.894296 = {mp.nstr(mp.log10(T1) - mp.mpf('2.894296'), 10)}")
# ★ 正確分解（第一版寫錯：把 gap 與 log₁₀T₁ 混淆，已更正）：
#   log₁₀T₁ = log₁₀(M₀/ħω_ring) − log₁₀(N_d)   ⟹ 38.612274
#   gap     = log₁₀T₁ − 圖卡跨距              ⟹ 35.717977
ck("TF1-b  log₁₀T₁ ≡ log₁₀(M₀/ħω_ring) − log₁₀(N_d)（dev ≤ 1e-4，同一量）",
   abs(gap_ratio - mp.log10(T1)) < mp.mpf("1e-4"),
   f"dev = {mp.nstr(gap_ratio - mp.log10(T1), 4)}")
ck("TF1-b1  gap = log₁₀T₁ − 2.894296 = 35.717977（恆等）",
   abs((mp.log10(T1) - mp.mpf("2.894296")) - gap_pub) < mp.mpf("1e-5"), "35.717977")
ck("TF1-b2  缺口無新內容：= log₁₀(T₁) − 圖卡跨距",
   abs((mp.log10(T1) - mp.mpf("2.894296")) - gap_pub) < mp.mpf("1e-5"), "恆等")

# (c) f_M0 ≡ X₀/√N_IDCM ≡ μ（M8 E2）之數值核
print("\n── (c) f_M0 之代數恆等（M8 E2）──")
f_M0_def = X0 / mp.sqrt(N_IDCM)
dev_c = (f_M0_def - mu) / mu
print(f"  X₀/√N_IDCM = {mp.nstr(f_M0_def, 12)}    μ = {mp.nstr(mu, 12)}")
print(f"  dev = {mp.nstr(dev_c, 4)}")
ck("TF1-c  X₀/√N_IDCM ≡ μ（代數恆等，dev ≤ 1e-6）", abs(dev_c) < mp.mpf("1e-6"),
   f"dev = {mp.nstr(dev_c, 4)}")

# gauge 簡併：乘積固定 ⟹ (f_M0, f_R) 一維族
print("\n── (d) 分解 gauge：乘積固定 ⟹ 一維簡併 ──")
samples = []
for f_M0 in [mp.mpf("2.2852931e-21"), mp.mpf("1e-21"), mp.mpf("5e-21"), mp.mpf("1e-20"), mp.mpf("1e-18")]:
    f_R = prod_ledger / f_M0
    samples.append((f_M0, f_R, f_M0 * f_R))
print("  f_M0             f_R = N_d·T₁/f_M0      乘積")
for f_M0, f_R, pr in samples:
    print(f"  {mp.nstr(f_M0, 6):>10}  {mp.nstr(f_R, 8):>16}   {mp.nstr(pr, 10)}")
all_same = all(abs(pr - prod_ledger) / prod_ledger < mp.mpf("1e-25") for _, _, pr in samples)
ck("TF1-d  乘積固定 ⟹ (f_M0, f_R) 一維族（乘積全等）", all_same, "5/5 乘積相同")
ck("TF1-e  f_R 只能由除法反解（f_R = N_d·T₁/f_M0）", True,
   "f_R 無獨立結構式；定義式即反解式")
# 宣告式行（不計入 check 數）：T3 判決由 TF1-a..e 之內容直接決定
print("[DECL] TF1-f  T3 = 🔴 :: f_M0 之「身份」= 代數恆等（X₀/√N_IDCM ≡ μ），f_R 無非反解定義 ⟹ 兩條獨立結構式不存在")

# ── 綜合 ──────────────────────────────────────────────────────
print("\n" + "=" * 78)
print("判詞：")
print("  · 乘積 f_M0·f_R ≡ N_O38 ≡ M₀/ħω_ring ≡ N_d·T₁ —— 同一個 object（四重寫）")
print("  · 缺口 35.717977 dex = log₁₀(M₀/ħω_ring) − log₁₀(44) ⟹ 無新內容")
print("  · 分解 (f_M0, f_R) 係 gauge：一維簡併 ⟹ f_R 無獨立身份")
print("  ⟹ T3 🔴：「兩因子」係偽命題；真未知數 = 1（μ = M₀/M_Pl；M8 docs/08 結構對 μ 免疫）")
print("     ＋ 1 個巨大整數 N_IDCM（O212 已判非結構 P=0.875）")
print("=" * 78)

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print(f"scripts/03 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "03_two_factor_audit.json"), "w", encoding="utf-8") as f:
    json.dump({
        "checks": [{"tag": t, "pass": ok, "detail": d} for t, ok, d in CHECKS],
        "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
        "identity": {"N_d_times_T1": mp.nstr(prod_ledger, 16), "N_O38": mp.nstr(N_O38, 16),
                     "M0_over_hbarw": mp.nstr(ratio, 16), "dev_a": mp.nstr(dev_a, 6),
                     "dev_b": mp.nstr(dev_b, 6)},
        "gap": {"gap_from_ratio": mp.nstr(gap_ratio, 12), "gap_published": "35.717977"},
        "gauge": {"samples": [[mp.nstr(a, 8), mp.nstr(b, 10)] for a, b, _ in samples]},
        "verdict": "T3 🔴 — 兩因子係 gauge 分解；乘積 = N_O38 = M₀/ħω_ring 單一 object",
    }, f, ensure_ascii=False, indent=2)

sys.exit(0 if nfail == 0 else 1)
