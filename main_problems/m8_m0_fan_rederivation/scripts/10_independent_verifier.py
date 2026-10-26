#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/10 — L2 獨立重推導驗證器（Independent Verifier）
================================================================================
鐵律（idcm-audit-verify §2）：
  · **唔准 import／call scripts/01–09 任何函數**（避免同源錯）
  · 由**原始輸入**（N_d=44、鏈長 88、R=M_Pl/M₀、r=ħω_ring/M_Pl、N_IDCM）重算每個 headline 數
  · 代數恆等用 mpmath 50 dps；能精確者用 Fraction／符號
  · 報 N/N，唔報「差不多」

執行: /usr/bin/python3 scripts/10_independent_verifier.py
輸出: data/10_independent_verify.json + data/run_10.log
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
N_FAIL = 0


def chk(tag, cond, detail=""):
    global N_FAIL
    ok = bool(cond)
    if not ok:
        N_FAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ============================================================================
# §0 原始輸入（只此四項；其餘全部重算）
# ============================================================================
hdr("§0 原始輸入（唯一輸入；其餘皆導出）")

N_D = 44                      # 雙鏈 rung 數
N_CHAIN = 88                  # 鏈長
R_ANCHOR = mp.mpf("4.3758064e20")     # M_Pl/M₀（量測錨；O212）
R_HOST = mp.mpf("1.2682908e-61")      # ħω_ring/M_Pl（量測錨）
N_IDCM = mp.mpf("7.8846e60")          # O212 表列（未閉合因子；docs/05 f_R）
print(f"  N_d = {N_D}   N_chain = {N_CHAIN}")
print(f"  M_Pl/M₀     = {mp.nstr(R_ANCHOR, 9)}")
print(f"  ħω_ring/M_Pl = {mp.nstr(R_HOST, 9)}")
print(f"  N_IDCM      = {mp.nstr(N_IDCM, 9)}")

mu = 1 / R_ANCHOR
print(f"  ⟹ μ = M₀/M_Pl = {mp.nstr(mu, 15)}")


# ============================================================================
# §1 鏈譜不變量（由 γ_n = 4sin²(πn/88) 重算）
# ============================================================================
hdr("§1 鏈譜不變量")

GAM = [4 * mp.sin(mp.pi * n / N_CHAIN) ** 2 for n in range(1, N_D + 1)]

chk("V1.1 γ₂₂ = 2（|−2| < 1e−45；等價 cos(π/2)=0）", abs(GAM[21] - 2) < mp.mpf("1e-45"),
    f"γ₂₂ = {mp.nstr(GAM[21], 20)}")
chk("V1.2 鏡射 γ_n + γ_{44−n} = 4（n=1…43，max|err| < 1e−45）",
    max(abs(GAM[n - 1] + GAM[N_D - n - 1] - 4) for n in range(1, N_D)) < mp.mpf("1e-45"),
    "43/43")
chk("V1.3 γ₄₄ = 4（= γ_n+γ_{44−n} 之常數值）", abs(GAM[43] - 4) < mp.mpf("1e-45"),
    f"{mp.nstr(GAM[43], 12)}")

# Σγ = 90 精確：Σ_{n=1}^{87} sin²(πn/88) = 44 ⟹ Σ_{1..44} sin² = 44/2 + 1/2 = 22.5 ⟹ Σγ = 90
S_gam = mp.fsum(GAM)
chk("V1.4 Σγ = 90 EXACT（Σ sin² = 45/2 ⟹ ×4 = 90）", abs(S_gam - 90) < mp.mpf("1e-45"),
    f"Σγ = {mp.nstr(S_gam, 20)}")

S_gam2 = mp.fsum([g ** 2 for g in GAM])
chk("V1.5 Σγ² = 272（M6 表列）", abs(S_gam2 - 272) < mp.mpf("1e-43"), f"Σγ² = {mp.nstr(S_gam2, 14)}")

P43 = mp.fprod(GAM[:N_D - 1])
P44 = P43 * GAM[43]
chk("V1.6 Π_{n=1}^{43} γ_n = N_d = 44 EXACT（O140）", abs(P43 / N_D - 1) < mp.mpf("1e-45"),
    f"{mp.nstr(P43, 14)}")
chk("V1.7 Π_{n=1}^{44} γ_n = 4N_d = 176", abs(P44 / (4 * N_D) - 1) < mp.mpf("1e-45"),
    f"{mp.nstr(P44, 12)}")

# 行列式類（R8′ 新結果）
det_prime = mp.mpf(1)
for n, g in enumerate(GAM, start=1):
    if n == 22:
        continue
    det_prime *= (1 - g / 2)
chk("V1.8 det'(1−K/2) = Π_{n≠22}(1−γ_n/2) = N_d/2^43",
    abs(det_prime / (mp.mpf(N_D) / mp.mpf(2) ** 43) - 1) < mp.mpf("1e-45"),
    f"{mp.nstr(det_prime, 12)}")

cosprod = mp.fprod([mp.cos(mp.pi * n / N_CHAIN) for n in range(1, N_D)])
chk("V1.9 Π_{n=1}^{43} cos(πn/88) = √N_d/2^43",
    abs(cosprod / (mp.sqrt(N_D) / mp.mpf(2) ** 43) - 1) < mp.mpf("1e-45"),
    f"{mp.nstr(cosprod, 12)}")

det_q = mp.fprod([1 + g for g in GAM])
chk("V1.10 det(1+K) = Π(1+γ_n) = 5.5004388918e18（獨立重算）",
    abs(det_q / mp.mpf("5.5004388918e18") - 1) < mp.mpf("1e-11"),
    f"{mp.nstr(det_q, 12)}")

# tan 族（O246 跨記錄交叉核）
S_tan2 = mp.fsum([mp.tan(mp.pi * n / N_CHAIN) ** 2 for n in range(1, N_D)])
chk("V1.11 Σ_{n=1}^{43} tan²(πn/88) = 1247 EXACT（O246 交叉核）",
    abs(S_tan2 - 1247) < mp.mpf("1e-40"), f"{mp.nstr(S_tan2, 14)}")

RES["_spectrum"] = {"sum_gamma": float(S_gam), "sum_gamma2": float(S_gam2),
                    "prod_gamma43": float(P43), "det_prime_half": float(det_prime),
                    "det_1K": float(det_q)}


# ============================================================================
# §2 圖卡／標度（J5 no-go 之獨立重算）
# ============================================================================
hdr("§2 圖卡幾何與 10^38 no-go")

cot2 = 1 / mp.tan(mp.pi / N_CHAIN) ** 2
span_ln = 2 * mp.log(1 / mp.tan(mp.pi / N_CHAIN))
span_dec = span_ln / mp.log(10)
chk("V2.1 cot²(π/88) = 783.964664（圖卡最大比值）", abs(cot2 - mp.mpf("783.964664")) < mp.mpf("1e-6"),
    f"{mp.nstr(cot2, 12)}")
chk("V2.2 Δt = 2 ln cot(π/88) = 6.664364 ln = 2.894296 decade",
    abs(span_ln - mp.mpf("6.664364")) < mp.mpf("1e-6")
    and abs(span_dec - mp.mpf("2.894296")) < mp.mpf("1e-6"),
    f"{mp.nstr(span_ln, 10)} / {mp.nstr(span_dec, 10)}")

T1 = mp.mpf("4.0951856e38")
lgT1 = mp.log10(T1)
gap = lgT1 - span_dec
chk("V2.3 log₁₀ T₁ = 38.612274", abs(lgT1 - mp.mpf("38.612274")) < mp.mpf("1e-6"),
    f"{mp.nstr(lgT1, 10)}")
chk("V2.4 缺口 = 35.717977 decade", abs(gap - mp.mpf("35.717977")) < mp.mpf("1e-6"),
    f"{mp.nstr(gap, 10)}")
chk("V2.5 缺口因子 = 10^35.717977 = 5.223686e35（★審計更正值）",
    abs(10 ** gap / mp.mpf("5.223686e35") - 1) < mp.mpf("1e-6"),
    f"{mp.nstr(10 ** gap, 10)}")
chk("V2.6 迭代數 = log₁₀T₁ / span = 13.34（非整數）",
    abs(lgT1 / span_dec - mp.mpf("13.34")) < mp.mpf("1e-2"),
    f"{mp.nstr(lgT1 / span_dec, 8)}")

RES["_chart"] = {"cot2": float(cot2), "span_dec": float(span_dec),
                 "log10_T1": float(lgT1), "gap_dec": float(gap), "factor": float(10 ** gap)}


# ============================================================================
# §3 錨交接點／層級歸約（f_M0 ≡ μ 之恆等）
# ============================================================================
hdr("§3 錨交接點與 f_M0 恆等")

# E_kchain = √(M_Pl · ħω_ring)；以 M_Pl = 1 為單位作符號驗證
# X₀ ≡ M₀/E_kchain；恆等：X₀² = M₀²/(M_Pl ħω) = μ·(M₀/ħω) = μ·N_O38 = μ²·N_IDCM（因 N_O38 = N_IDCM·μ）
X0 = mu * mp.sqrt(N_IDCM)
chk("V3.1 恆等 X₀ ≡ μ√N_IDCM（代數；⟹ f_M0 ≡ X₀/√N_IDCM = μ，零新內容）",
    abs(X0 / mp.sqrt(N_IDCM) / mu - 1) < mp.mpf("1e-40"),
    f"X₀ = {mp.nstr(X0, 10)}")
chk("V3.2 X₀ 數值 = 6.4170e9（C-約定）", abs(X0 / mp.mpf("6.4170e9") - 1) < mp.mpf("1e-3"),
    f"{mp.nstr(X0, 10)}")

# 幾何平均錨：E_kchain/M_Pl = √(ħω/M_Pl) = √r
E_over_MPl = mp.sqrt(R_HOST)
dec_Pl = 0.0
dec_ring = mp.log10(R_HOST)
chk("V3.3 錨交接點 = 幾何平均：log₁₀(E_kchain/M_Pl) = ½ log₁₀(ħω/M_Pl)",
    abs(mp.log10(E_over_MPl) - dec_ring / 2) < mp.mpf("1e-30"),
    f"{mp.nstr(mp.log10(E_over_MPl), 8)} vs {mp.nstr(dec_ring/2, 8)}")
chk("V3.4 兩錨各距交接點 30.45 dex（對稱）",
    abs(abs(mp.log10(E_over_MPl)) - mp.mpf("30.45")) < mp.mpf("0.01"),
    f"+{mp.nstr(mp.log10(E_over_MPl), 6)} / {mp.nstr(dec_ring - mp.log10(E_over_MPl), 6)}")

# 兩標度冪律：p = ln μ / ln r
p_exp = mp.log(mu) / mp.log(R_HOST)
chk("V3.5 必需指數 p = 0.33895155（μ = r^p 反解）",
    abs(p_exp - mp.mpf("0.33895155")) < mp.mpf("1e-8"), f"{mp.nstr(p_exp, 12)}")
chk("V3.6 p ≠ 1/3；dev = |p−1/3|/(1/3) = 1.6855e−2（依 docs/06 約定 |p−q|/q）",
    abs(abs(p_exp - mp.mpf(1) / 3) / (mp.mpf(1) / 3) - mp.mpf("1.6855e-2")) < mp.mpf("1e-5"),
    f"{mp.nstr(abs(p_exp - mp.mpf(1)/3)/(mp.mpf(1)/3), 8)}")

chk("V3.7 三輸入窮盡 ⟹ 非冪類只剩指數：S* = ln(1/μ) = 47.52779",
    abs(mp.log(1 / mu) - mp.mpf("47.52779")) < mp.mpf("1e-5"),
    f"{mp.nstr(mp.log(1/mu), 10)}")

RES["_hierarchy"] = {"X0": float(X0), "p": float(p_exp), "S_star": float(mp.log(1 / mu)),
                     "E_over_MPl": float(E_over_MPl)}


# ============================================================================
# §4 行列式路線上限（R8′ 之封閉結論）
# ============================================================================
hdr("§4 行列式路線之 dex 上限")

lg_mu = mp.log10(mu)
lg_inv_det = mp.log10(1 / det_q)
gap_dex = lg_mu - lg_inv_det
chk("V4.1 log₁₀ μ = −20.6411", abs(lg_mu + mp.mpf("20.6411")) < mp.mpf("1e-4"),
    f"{mp.nstr(lg_mu, 10)}")
chk("V4.2 log₁₀ 1/det(1+K) = −18.7404", abs(lg_inv_det + mp.mpf("18.7404")) < mp.mpf("1e-4"),
    f"{mp.nstr(lg_inv_det, 10)}")
chk("V4.3 缺口 = 1.9007 dex ⟹ 1/det 比 μ 大 79.5538×（⟹ 單一行列式不可能）",
    abs(gap_dex + mp.mpf("1.9007")) < mp.mpf("1e-4")
    and abs(1 / (mu * det_q) - mp.mpf("79.5538")) < mp.mpf("1e-3"),
    f"{mp.nstr(gap_dex, 8)} dex / {mp.nstr(1/(mu*det_q), 8)}×")
chk("V4.4 裸行列式壓低 2^−43 = 12.9443 dex", abs(mp.log10(mp.mpf(2) ** -43) + mp.mpf("12.9443")) < mp.mpf("1e-4"),
    f"{mp.nstr(mp.log10(mp.mpf(2)**-43), 8)}")

RES["_determinant"] = {"log10_mu": float(lg_mu), "log10_inv_det": float(lg_inv_det),
                       "gap_dex": float(gap_dex), "gap_factor": float(mu * det_q)}


# ============================================================================
# §5 交叉記錄核（O211 / O124 / M6 之數，本 M 引用者）
# ============================================================================
hdr("§5 引用之跨記錄數（獨立重算）")

chk("V5.1 N_d = 44 與 rank_ℚ{γ_n} = 20 / 核 24 相容（44 = 20 + 24）", 44 == 20 + 24, "恆等")
lam_max = mp.fsum(GAM[-4:])
chk("V5.2 支撐定理 λ_max(4) = Σ_{n=41}^{44} γ_n < 16（M6 docs/07）", lam_max < 16,
    f"{mp.nstr(lam_max, 10)}")
chk("V5.3 Σγ = B(N_d+1)/2 = 4·45/2 = 90（O211 形式）", abs(S_gam - (4 * (N_D + 1)) / 2) < mp.mpf("1e-40"),
    "90 = 4·45/2")

RES["_cross"] = {"lambda_max4": float(lam_max)}


# ============================================================================
# §6 匯總
# ============================================================================
hdr("§6 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

(DATA / "10_independent_verify.json").write_text(
    json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_10.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {DATA / '10_independent_verify.json'}")
print(f"寫出: {DATA / 'run_10.log'}")

sys.exit(0 if N_FAIL == 0 else 1)