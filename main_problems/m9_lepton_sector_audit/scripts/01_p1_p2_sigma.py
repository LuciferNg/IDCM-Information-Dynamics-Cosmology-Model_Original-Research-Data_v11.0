#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/01 — P1 觀測精度權威核實 ＋ P2 三聲稱 σ 級距表

零 import 上游腳本；零網絡（所有常數同來源 sha256 喺 docs/01 §一 凍結）。
精確度：mpmath dps = 60。
輸出：data/01_p1_p2_sigma.json ＋ stdout 摘要（PASS/FAIL 計數）
"""
import json
import os
import re

from mpmath import mp, mpf, log, sqrt

mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

PHI = (1 + sqrt(5)) / 2
LN_PHI = log(PHI)

# ════════════════════════════════════════════════════════════════════
# §一 權威輸入 — CODATA 2022 (NIST ascii 表)
#   來源檔 sha256 = 77fb90e66c40db3e6eb16630bc9c88e4c7c8beddbe5e71be406f2f26e3f67e67
#   行號：107-110 m_e、215 直接比 μ/e、221-224 m_μ、336-342 m_τ
# ════════════════════════════════════════════════════════════════════
M_E = mpf("0.51099895069")          # MeV      (line 109)
U_M_E = mpf("0.00000000016")        # MeV
M_MU = mpf("105.6583755")           # MeV      (line 223)
U_M_MU = mpf("0.0000023")           # MeV
M_TAU_KG = mpf("3.16754e-27")       # kg       (line 338)
U_M_TAU_KG = mpf("0.00021e-27")     # kg
C_LIGHT = mpf("299792458")          # m/s exact
E_CHARGE = mpf("1.602176634e-13")   # J/MeV exact
M_TAU = M_TAU_KG * C_LIGHT ** 2 / E_CHARGE          # MeV
U_M_TAU = U_M_TAU_KG * C_LIGHT ** 2 / E_CHARGE

# 直接比（CODATA 自己嘅 ratio entry；已含相關性）
R_MUE = mpf("206.7682827");   U_R_MUE = mpf("0.0000046")     # line 215
R_TAUMU = mpf("16.8170");     U_R_TAUMU = mpf("0.0011")      # line 342
R_TAUE = mpf("3477.23");      U_R_TAUE = mpf("0.23")         # line 336

# ════════════════════════════════════════════════════════════════════
# §二 譜（擁有者 M4/O137；本 M 只引用定義，不重推）
# ════════════════════════════════════════════════════════════════════
def gamma(n: int) -> "mpf":
    """γ_n = 4 sin²(πn/88)，n=1..44（n=44 = 4 精確）"""
    return 4 * mp.sin(mp.pi * n / 88) ** 2


G1 = gamma(1)
G7 = gamma(7)

results = {}
checks = []


def chk(name, cond, detail=""):
    checks.append({"check": name, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


print("=" * 72)
print("M9 scripts/01 — P1 觀測精度 ＋ P2 σ 級距表 (mp.dps=60)")
print("=" * 72)

# ── §三 P1：觀測量與不確定度 ──────────────────────────────────────
r_mue_ratio_of_masses = M_MU / M_E
u_r_mue_prop = r_mue_ratio_of_masses * sqrt((U_M_MU / M_MU) ** 2 + (U_M_E / M_E) ** 2)
r_taumu_ratio_of_masses = M_TAU / M_MU
u_r_taumu_prop = r_taumu_ratio_of_masses * sqrt((U_M_TAU / M_TAU) ** 2 + (U_M_MU / M_MU) ** 2)
r_taue_ratio_of_masses = M_TAU / M_E
u_r_taue_prop = r_taue_ratio_of_masses * sqrt((U_M_TAU / M_TAU) ** 2 + (U_M_E / M_E) ** 2)

print("\n[P1] 權威質量與比值 (CODATA 2022)")
print(f"  m_e   = {mp.nstr(M_E, 12)} MeV  u = {mp.nstr(U_M_E, 3)}  (rel {mp.nstr(U_M_E / M_E, 4)})")
print(f"  m_mu  = {mp.nstr(M_MU, 12)} MeV  u = {mp.nstr(U_M_MU, 3)}  (rel {mp.nstr(U_M_MU / M_MU, 4)})")
print(f"  m_tau = {mp.nstr(M_TAU, 12)} MeV  u = {mp.nstr(U_M_TAU, 4)}  (rel {mp.nstr(U_M_TAU / M_TAU, 4)})")
print(f"  ratio-of-masses  mu/e {mp.nstr(r_mue_ratio_of_masses, 12)} +- {mp.nstr(u_r_mue_prop, 3)}"
      f"   (rel {mp.nstr(u_r_mue_prop / r_mue_ratio_of_masses, 4)})")
print(f"  direct-entry     mu/e {mp.nstr(R_MUE, 11)} +- {mp.nstr(U_R_MUE, 3)}"
      f"   (rel {mp.nstr(U_R_MUE / R_MUE, 4)})")
print(f"  ratio-of-masses  tau/mu {mp.nstr(r_taumu_ratio_of_masses, 10)} +- {mp.nstr(u_r_taumu_prop, 3)}"
      f"   (rel {mp.nstr(u_r_taumu_prop / r_taumu_ratio_of_masses, 3)})")
print(f"  direct-entry     tau/mu {mp.nstr(R_TAUMU, 7)} +- {mp.nstr(U_R_TAUMU, 3)}"
      f"   (rel {mp.nstr(U_R_TAUMU / R_TAUMU, 3)})")
print(f"  direct-entry     tau/e  {mp.nstr(R_TAUE, 8)} +- {mp.nstr(U_R_TAUE, 3)}"
      f"   (rel {mp.nstr(U_R_TAUE / R_TAUE, 3)})")

# 兩個獨立路徑（比值 vs 比值之比）必須一致喺 σ 內 —— 否則有相關性未計
d_consistency = abs(r_mue_ratio_of_masses - R_MUE) / sqrt(u_r_mue_prop ** 2 + U_R_MUE ** 2)
print(f"  consistency |mu/e(ratio-of-masses) - mu/e(direct)| = "
      f"{mp.nstr(abs(r_mue_ratio_of_masses - R_MUE), 4)} = {mp.nstr(d_consistency, 3)}σ")
chk("P1-a 兩路徑 μ/e 一致 (<1σ)", d_consistency < 1,
    f"{mp.nstr(d_consistency, 3)}σ")

# M9 record 用嘅值 206.768282988（= m_μ/m_e from MeV），需標明差異
rec_val = mpf("206.768282988")
print(f"  record §三 用值 {mp.nstr(rec_val, 12)} vs direct entry {mp.nstr(R_MUE, 11)}"
      f"  Δ = {mp.nstr(rec_val - R_MUE, 4)} = {mp.nstr(abs(rec_val - R_MUE) / U_R_MUE, 3)}σ_ratio")
chk("P1-b record 用值與 direct entry 相容 (<=1σ)",
    abs(rec_val - R_MUE) / U_R_MUE < 1,
    f"{mp.nstr(abs(rec_val - R_MUE) / U_R_MUE, 3)}σ")

# ★ 判準核心：σ 分母 = 直接比 entry 之不確定度
SIG_MUE = U_R_MUE                     # 絕對
REL_MUE = U_R_MUE / R_MUE             # 相對
print(f"\n  ★ σ 判準（μ/e）: 絕對 u = {mp.nstr(SIG_MUE, 3)}，相對 u = {mp.nstr(REL_MUE, 4)}")
print(f"  ★ 5σ 相容所需相對精度 = {mp.nstr(5 * REL_MUE, 4)}  (即 {mp.nstr(5 * REL_MUE * 100, 4)} %)")

# ── §四 P2：三聲稱 σ 級距 ────────────────────────────────────────
print("\n[P2] 三聲稱 σ 級距（分母 = CODATA 直接比 u，唔係 % 級標準）")

zeta = 11 + mpf(7) / 88
C1 = PHI ** zeta
C3 = mpf(207)
g_ratio = G7 / G1
f_mu_over_e_quoted = mpf("4.3069")     # memory 式引用值（4 位）
f_mu_over_e_exact = R_MUE / g_ratio    # 反解（= 定義）
C2_quoted = g_ratio * f_mu_over_e_quoted
C2_exact = g_ratio * f_mu_over_e_exact


def sigma_row(label, value, target, u):
    dev = value - target
    rel = dev / target
    return {
        "label": label,
        "value": mp.nstr(value, 12),
        "target": mp.nstr(target, 11),
        "dev_abs": mp.nstr(dev, 10),
        "dev_rel": mp.nstr(rel, 10),
        "dev_pct": mp.nstr(rel * 100, 10),
        "sigma": mp.nstr(abs(dev) / u, 10),
    }


rows = []
rows.append(sigma_row("C1 phi^(11+7/88)", C1, R_MUE, SIG_MUE))
rows.append(sigma_row("C2 (g7/g1)(F_mu/F_e) quoted 4.3069", C2_quoted, R_MUE, SIG_MUE))
rows.append(sigma_row("C2' (g7/g1)(F_mu/F_e) exact-def", C2_exact, R_MUE, SIG_MUE))
rows.append(sigma_row("C3 integer 207", C3, R_MUE, SIG_MUE))

for r in rows:
    print(f"  {r['label']:38s} val={r['value']:>15s}  dev%={r['dev_pct']:>12s}  σ={r['sigma']:>12s}")

# C2 tautology 判據：exact-def 版本 dev ≡ 0（因 F_mu/F_e 由目標反解）
chk("P2-a C2 exact-def 版本 dev = 0（= 定義，1 dof 消耗，非預言）",
    abs(C2_exact - R_MUE) < mpf("1e-40"),
    f"|Δ| = {mp.nstr(abs(C2_exact - R_MUE), 3)}")
chk("P2-b C1 排除 (>5σ)", abs(C1 - R_MUE) / SIG_MUE > 5,
    f"{mp.nstr(abs(C1 - R_MUE) / SIG_MUE, 6)}σ")
chk("P2-c C3 排除 (>5σ)", abs(C3 - R_MUE) / SIG_MUE > 5,
    f"{mp.nstr(abs(C3 - R_MUE) / SIG_MUE, 6)}σ")

# 「要幾鬆嘅門檻才會 🟢」——解釋三個 🟢 並存之機制
print("\n  門檻敏感度（每個聲稱第一次會顯示為 🟢 所需之 % 門檻）:")
thresh = []
for r in rows:
    if r["label"].startswith("C2'"):
        continue
    pct = abs(mpf(r["dev_pct"]))
    thresh.append({"claim": r["label"], "min_pct_threshold": mp.nstr(pct, 6),
                   "sigma": r["sigma"]})
    print(f"    {r['label']:32s} 需 ≥ {mp.nstr(pct, 6):>12s} %   但實際 {r['sigma']:>12s} σ")

# ── §五 重現閘（SESSION_START §四 逐項核）───────────────────────
print("\n[重現閘] SESSION_START §四 表逐項")
logphi_obs = log(R_MUE) / LN_PHI
dev_x = logphi_obs - zeta
rep = {
    "m_mu/m_e (CODATA direct)": (mp.nstr(R_MUE, 11), "206.7682827"),
    "m_tau/m_mu (CODATA direct)": (mp.nstr(R_TAUMU, 7), "16.8170"),
    "phi^(11+7/88)": (mp.nstr(C1, 12), "206.770253211"),
    "dev_rel": (mp.nstr(C1 / R_MUE - 1, 10), "9.53e-6"),
    "log_phi(obs)": (mp.nstr(logphi_obs, 10), "11.07952565"),
    "11+7/88": (mp.nstr(zeta, 10), "11.07954545"),
    "x residual": (mp.nstr(dev_x, 6), "1.98e-5"),
    "gamma7/gamma1": (mp.nstr(g_ratio, 9), "48.008835"),
}
for k, (got, doc) in rep.items():
    print(f"  {k:28s} got={got:>16s}   doc={doc}")
chk("重現-a φ^(11+7/88) = 206.770253211", abs(C1 - mpf("206.770253211")) < mpf("1e-9"))
chk("重現-b dev_rel = 9.53e-6", abs((C1 / R_MUE - 1) - mpf("9.53e-6")) < mpf("1e-7"),
    mp.nstr(C1 / R_MUE - 1, 10))
chk("重現-c log_φ(obs) = 11.07952565", abs(logphi_obs - mpf("11.07952565")) < mpf("1e-9"))
chk("重現-d γ7/γ1 = 48.008835", abs(g_ratio - mpf("48.008835")) < mpf("1e-6"),
    mp.nstr(g_ratio, 12))
chk("重現-e log_φ(m_τ/m_μ) = 5.865175951（用 record 輸入 16.8170293324）",
    abs(log(mpf("16.8170293324")) / LN_PHI - mpf("5.865175951")) < mpf("1e-9"),
    mp.nstr(log(mpf("16.8170293324")) / LN_PHI, 12))
_x_codata = log(R_TAUMU) / LN_PHI
_x_record = log(mpf("16.8170293324")) / LN_PHI
print(f"  ⚠ record 之 5.865175951 建基於 16.8170293324 (=1776.86/105.6583755)，"
      f"CODATA 直接比 16.8170 給 {mp.nstr(_x_codata, 12)}")
print(f"    Δ_x = {mp.nstr(abs(_x_codata - _x_record), 3)} = "
      f"{mp.nstr(abs(_x_codata - _x_record) * LN_PHI / (U_R_TAUMU / R_TAUMU), 3)}σ（τ/μ 精度 6.5e−5）"
      f" ⟹ record 數字尾部無物理意義（超權威精度）")
chk("重現-f record τ/μ 值與 CODATA 直接比在 1σ 內",
    abs(_x_codata - _x_record) * LN_PHI / (U_R_TAUMU / R_TAUMU) < 1,
    f"{mp.nstr(abs(_x_codata - _x_record) * LN_PHI / (U_R_TAUMU / R_TAUMU), 3)}σ")

# C2 出處核：206.705 可否由 (γ₇/γ₁)(F_μ/F_e = 4.306880) 重現？
_c2_doc = mpf("206.705")
F_from_doc = _c2_doc / g_ratio
print(f"  ⚠ C2 記載值 206.705 反推 F_μ/F_e = {mp.nstr(F_from_doc, 8)}"
      f"（≠ M4 docs/04 反解值 4.306880）")
print(f"    以 4.306880 計：(γ₇/γ₁)·4.306880 = {mp.nstr(C2_exact, 12)} = 目標（構造恆等）")
chk("重現-g C2 記載 206.705 不可由 (γ₇/γ₁)·4.306880 重現",
    abs(g_ratio * mpf("4.306880") - _c2_doc) > mpf("0.01"),
    f"{mp.nstr(g_ratio * mpf('4.306880'), 10)} vs 206.705")

# ── §六 匯出 ────────────────────────────────────────────────────
results = {
    "source": {
        "name": "CODATA 2022 (NIST allascii.txt)",
        "sha256": "77fb90e66c40db3e6eb16630bc9c88e4c7c8beddbe5e71be406f2f26e3f67e67",
        "lines": {"m_e": 109, "mu/e direct": 215, "m_mu": 223, "m_tau": 338,
                  "tau/mu": 342, "tau/e": 336},
    },
    "observables": {
        "m_e_MeV": [mp.nstr(M_E, 12), mp.nstr(U_M_E, 3)],
        "m_mu_MeV": [mp.nstr(M_MU, 12), mp.nstr(U_M_MU, 3)],
        "m_tau_MeV": [mp.nstr(M_TAU, 12), mp.nstr(U_M_TAU, 4)],
        "mu_e_direct": [mp.nstr(R_MUE, 11), mp.nstr(U_R_MUE, 3)],
        "mu_e_ratio_of_masses": [mp.nstr(r_mue_ratio_of_masses, 12), mp.nstr(u_r_mue_prop, 3)],
        "tau_mu_direct": [mp.nstr(R_TAUMU, 7), mp.nstr(U_R_TAUMU, 3)],
        "tau_e_direct": [mp.nstr(R_TAUE, 8), mp.nstr(U_R_TAUE, 3)],
    },
    "sigma_criterion": {"abs_u": mp.nstr(SIG_MUE, 6), "rel_u": mp.nstr(REL_MUE, 6),
                        "rel_u_5sigma_need": mp.nstr(5 * REL_MUE, 6)},
    "claims": rows,
    "threshold_sensitivity": thresh,
    "reproduction": rep,
    "checks": checks,
    "pass": sum(1 for c in checks if c["pass"]),
    "total": len(checks),
}

out = os.path.join(DATA, "01_p1_p2_sigma.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2)

print("\n" + "=" * 72)
print(f"PASS {results['pass']}/{results['total']}   -> {out}")
print("=" * 72)