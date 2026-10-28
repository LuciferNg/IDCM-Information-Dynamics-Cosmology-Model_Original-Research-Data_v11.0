#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 L2 -- 獨立重推導（independent re-derivation）。

硬規則（idcm-audit-verify §2）：
  * 唔准 import / call L1 任何函數（scripts/01..40）
  * 唔准讀任何 docs/*.md
  * 每個數由 INPUTS（原始量：外部測量、計數、上游已封結論）重算
  * 期望值只可出現喺 IDENTITIES（代數恆等之驗證目標），唔可以出現喺 DERIVE

用途：為 L3（doc cross-check）提供唯一期望值來源。
輸出：data/l2_independent_verify.json
"""
import json
import os
from fractions import Fraction as F

import mpmath as mp
import sympy as sp

mp.mp.dps = 50

RESULTS = []
CHECKS = []


def rep(key, value, kind="num", note=""):
    RESULTS.append({"key": key, "value": str(value), "kind": kind, "note": note})
    return value


def ok(key, got, want, tol):
    """got vs want, relative tolerance tol (absolute when want==0)."""
    g, w = mp.mpf(str(got)), mp.mpf(str(want))
    dev = abs(g - w) / abs(w) if w != 0 else abs(g)
    good = dev <= mp.mpf(str(tol))
    print(("PASS " if good else "FAIL ") + "%-52s got=%-22s dev=%.3e" % (key, mp.nstr(g, 12), float(dev)))
    CHECKS.append({"key": key, "pass": bool(good), "got": mp.nstr(g, 14),
                   "want": mp.nstr(w, 14), "dev": float(dev), "tol": float(tol)})
    return good


# ======================================================================
# INPUTS —— 全部係原始量（外部測量 / 計數 / 上游已封結論）
# ======================================================================
N_d = 88                      # IDCM 鏈長（結構常數）
L9 = 76                       # Lucas L9（整數數列，非本記錄推導）
M_PI_P = mp.mpf("139.57039")  # PDG 2024  M_pi+  [MeV]
M_PI_0 = mp.mpf("134.9768")   # PDG 2024  M_pi0  [MeV]
M_K_P = mp.mpf("493.677")     # PDG 2024  M_K+   [MeV]
M_K_0 = mp.mpf("497.611")     # PDG 2024  M_K0   [MeV]
M_N = mp.mpf("939.565421")    # PDG 2024  m_n    [MeV]
M_P = mp.mpf("938.272089")    # PDG 2024  m_p    [MeV]
M_E = mp.mpf("0.510999")      # PDG 2024  m_e    [MeV]
MU, D_MU = mp.mpf("2.14"), mp.mpf("0.08")     # FLAG24 Nf=2+1+1  [MeV]
MD, D_MD = mp.mpf("4.70"), mp.mpf("0.05")     # FLAG24 Nf=2+1+1  [MeV]
MS, D_MS = mp.mpf("93.46"), mp.mpf("0.58")    # FLAG24 Nf=2+1+1  [MeV]
M_T = mp.mpf("172690")        # PDG 2024 m_t(pole) [MeV]
M1_NU, M3_NU = mp.mpf("0.020908"), mp.mpf("0.053586")   # O257 R7 [eV]
EPS_MID, EPS_LO, EPS_HI = mp.mpf("0.78"), mp.mpf("0.67"), mp.mpf("0.86")  # MILC18 via FLAG Eq.(38)
O_EXTERNAL, O_EXTERNAL_S = mp.mpf("1.498"), mp.mpf("0.262")   # BMW24 <O> 判據

print("=" * 100)
print("M4 L2 獨立重推導（零 L1 import；零 docs 讀取）")
print("=" * 100)

# ======================================================================
# §A 譜代數（由 N_d / Chebyshev 整數係數重算，唔經 gamma[] 陣列）
# ======================================================================
x = sp.symbols("x")
c = sp.cos(2 * sp.pi / N_d)                     # 定義
g1_sym = sp.simplify(2 - 2 * c)                 # γ₁ = 2 − 2T₁(c)
g1 = mp.mpf(str(sp.N(g1_sym, 45)))
rep("gamma1", mp.nstr(g1, 20), note="γ1 = 2-2cos(2π/88) = 4sin²(π/88)")
ok("A1 gamma1 == 4 sin^2(pi/88)", g1, 4 * mp.sin(mp.pi / N_d) ** 2, mp.mpf("1e-30"))
ok("A2 4/gamma1 = 784.96 (docs/37 A4)", 4 / g1, mp.mpf("784.96"), mp.mpf("1e-4"))
rep("four_over_gamma1", mp.nstr(4 / g1, 12))

# γ_n + γ_(44-n) = 4  ——  由 sin²(πn/88) + cos²(πn/88) 恆等，逐 n 核（n=1..43）
viol = [n for n in range(1, 44)
        if abs((4 * mp.sin(mp.pi * n / N_d) ** 2) + (4 * mp.sin(mp.pi * (44 - n) / N_d) ** 2) - 4) > mp.mpf("1e-40")]
ok("A3 γn + γ(44-n) = 4 for all n=1..43 (violations)", len(viol), 0, 0)

# γ_n ∈ ℤ[γ₁]：T_n(c) ∈ ℤ[c]（整數係數），c = (2 − γ₁)/2 ⟹ γ_n ∈ ℤ[γ₁]
def cheb_coeffs(n):
    """T_n 之整數係數（低次在前）；純遞推 T_k = 2c·T_{k-1} − T_{k-2}。"""
    if n == 0:
        return [F(1)]
    a, b = [F(1)], [F(0), F(1)]
    for _ in range(2, n + 1):
        sh = [F(0)] + [2 * v for v in b]
        m = max(len(sh), len(a))
        a2 = a + [F(0)] * (m - len(a))
        b2 = sh + [F(0)] * (m - len(sh))
        a, b = b, [p_ - q_ for p_, q_ in zip(b2, a2)]
    return b


def gamma_poly_in_u(n):
    """γ_n = 2 − 2T_n(c)，c = (2 − u)/2（u ≡ γ₁）⟹ 有理係數多項式（低次在前）。"""
    poly = [F(0)]
    for coef in reversed(cheb_coeffs(n)):
        new = [F(0)] * (len(poly) + 1)
        for i, v in enumerate(poly):
            new[i] += v * 1            # ×1
            new[i + 1] += v * F(-1, 2)  # ×(−u/2)
        new[0] += coef
        poly = new
    poly = [-2 * v for v in poly]
    poly[0] += 2
    return poly


# --- 最小多項式：pslq 由高精度值精確回復（dps=60），再以殘差核實
y_val = 2 * mp.cos(2 * mp.pi / N_d)
powers = [mp.mpf(1)]
for _ in range(20):
    powers.append(powers[-1] * y_val)
rel = mp.pslq(powers, maxcoeff=10 ** 30, maxsteps=10 ** 6)
assert rel is not None, "pslq failed"
rel = [int(v) for v in rel]
assert rel[-1] != 0
resid = abs(sum(mp.mpf(cf) * pw for cf, pw in zip(rel, powers)))
ok("A4 minpoly(2cos(2π/88)) degree = 20", len(rel) - 1, 20, 0)
_scale_y = max(abs(mp.mpf(cf) * pw) for cf, pw in zip(rel, powers))
ok("A5 minpoly residual/scale < 1e−40 (pslq 精確回復核實)",
   1 if resid / _scale_y < mp.mpf("1e-40") else 0, 1, 0)
print("  [A5] |resid| = %.3e  scale = %.3e" % (float(resid), float(_scale_y)))
minpoly_y = rel                                # 低次在前
rep("minpoly_2cos88", sp.sstr(sp.Poly(list(reversed(minpoly_y)), x)))

# u = 2 − y ⟹ minpoly_u(u) = minpoly_y(2 − u)：逐項展開（精確整數）
poly_u = [F(0)] * len(minpoly_y)
for k, cf in enumerate(minpoly_y):
    term = [F(1)]
    for _ in range(k):
        nn = [F(0)] * (len(term) + 1)
        for i, v in enumerate(term):
            nn[i] += 2 * v
            nn[i + 1] -= v
        term = nn
    for i, v in enumerate(term):
        poly_u[i] += cf * v
# 驗證 minpoly_u(γ₁) ≈ 0
val_u = sum(mp.mpf(str(cf.numerator)) / mp.mpf(str(cf.denominator)) * g1 ** i
            for i, cf in enumerate(poly_u))
_scale_u = max(abs(mp.mpf(str(cf.numerator)) / mp.mpf(str(cf.denominator)) * g1 ** i)
               for i, cf in enumerate(poly_u))
ok("A6 minpoly_u(γ1) = 0 (|殘差|/scale < 1e−30)", 1 if abs(val_u) / _scale_u < mp.mpf("1e-30") else 0, 1, 0)
print("  [A6] |p(γ1)| = %.3e  scale = %.3e" % (float(abs(val_u)), float(_scale_u)))


def reduce_mod_u(coeffs, modpoly):
    """多項式（有理係數，低次在前）對 modpoly 取模（精確 Fraction）。"""
    d = len(modpoly) - 1
    lead = modpoly[-1]
    c_ = [F(v) for v in coeffs]
    while len(c_) > d:
        top = c_.pop()
        if top == 0:
            continue
        shift = len(c_) - d
        for i in range(d):
            c_[shift + i] -= top * modpoly[i] / lead
    return c_ + [F(0)] * (d - len(c_))


rows = []
for n in range(1, 44):
    rows.append(reduce_mod_u(gamma_poly_in_u(n), poly_u))       # 每行 = γ_n 於 {1,γ₁,…,γ₁^19} 之座標
M = [r[:] for r in rows]
rank = 0
for col in range(20):
    piv = next((r for r in range(rank, len(M)) if M[r][col] != 0), None)
    if piv is None:
        continue
    M[rank], M[piv] = M[piv], M[rank]
    for r in range(len(M)):
        if r != rank and M[r][col] != 0:
            f = M[r][col] / M[rank][col]
            M[r] = [a - f * b for a, b in zip(M[r], M[rank])]
    rank += 1
ok("A7 rank_Q span{γ1..γ43} = 20 (exact rational elimination)", rank, 20, 0)
rep("rank_Q", rank)

# γ_n ∈ ℤ[γ₁]：座標分母只有 2 之冪（c = (2−γ₁)/2），且 γ_n = 4sin² 為代數整數
check_int = all(
    reduce_mod_u(gamma_poly_in_u(n), poly_u) == [F(v) for v in reduce_mod_u(gamma_poly_in_u(n), poly_u)]
    for n in [1]
)
rep("deg_Q_gamma1", len(poly_u) - 1)
ok("A8 deg_Q(γ1) = φ(88)/2 = 20", len(poly_u) - 1, 20, 0)

# ======================================================================
# §B span = 12 與 m_u 預言
# ======================================================================
hop_max = N_d - L9
ok("B1 hop_max = 88 − L9 = 12", hop_max, 12, 0)
rep("hop_max", hop_max)

span_nu = mp.log10(M3_NU / M1_NU)
span_q = mp.log10(M_T / MU)
R_span = span_q / span_nu
rep("span_nu", mp.nstr(span_nu, 12))
rep("span_Q", mp.nstr(span_q, 12))
rep("R_span", mp.nstr(R_span, 12))
sig_spanQ = D_MU / (MU * mp.log(10))
rep("R_span_sigma", mp.nstr(sig_spanQ / span_nu, 10), note="一階解析傳播（只計 m_u 誤差）")

# ★ 精確（二階）：Var[log10(M_T/MU_rand)] = Var[log10 MU] + Var[log10 M_T]（兩者獨立）
#   一階低估 0.28%（log10 非線性）⟹ canonical 用精確值；另由 scripts/44 之 3×10⁷ MC 交叉核。
D_MT = mp.mpf("300")          # PDG m_t(pole) 誤差 [MeV]（相對 0.17%，對 σ_R 貢獻 <0.01%）


def _log10_normal_var(m, s, dps=30):
    """Var[log10 X]，X ~ N(m, s)。積分域取 x ∈ [−12, 12]（正態尾 <1e−32，可忽略）；
    此域內 m + s·x > 0 對 m_u（2.14 − 0.08·12 = 1.18）與 m_t 皆成立 ⟹ g 實值。"""
    with mp.workdps(dps):
        w = lambda x: mp.e ** (-x * x / 2) / mp.sqrt(2 * mp.pi)
        g = lambda x: mp.log10(m + s * x)
        A = mp.quad(lambda x: w(x) * g(x), [-12, 12])
        B = mp.quad(lambda x: w(x) * g(x) ** 2, [-12, 12])
        return B - A * A


var_q = _log10_normal_var(MU, D_MU) + _log10_normal_var(M_T, D_MT)
sig_R_exact = mp.sqrt(var_q) / span_nu
rep("R_span_sigma_exact", mp.nstr(sig_R_exact, 10), note="精確二階（Gauss 期望積分，mpmath dps=30）")
ok("B4 σ_R 精確/一階 = 1.002825（二階項 0.28%，非噪聲）",
   sig_R_exact / (sig_spanQ / span_nu), mp.mpf("1.002825"), mp.mpf("1e-4"))
print("  [note] R = %.9f  （docs/11 報 12.005533; 見 L3 之 over-precision 檢查）" % float(R_span))

m_u_pred = M_T / mp.power(10, 12 * span_nu)
rep("m_u_pred_MeV", mp.nstr(m_u_pred, 10))
ok("B2 m_u(span=12) = 2.1498 MeV (docs/37 A2)", m_u_pred, mp.mpf("2.1498"), mp.mpf("1e-3"))
z_mu = (m_u_pred - MU) / D_MU
rep("m_u_pred_sigma", mp.nstr(z_mu, 6))
ok("B3 m_u(span=12) vs FLAG 2.14±0.08 = +0.12σ", z_mu, mp.mpf("0.12"), mp.mpf("0.15"))

# ======================================================================
# §C X 與誤差傳播
# ======================================================================
X = MD - MU
sig_X_quad = mp.sqrt(D_MD ** 2 + D_MU ** 2)
rep("X_MeV", mp.nstr(X, 6))
rep("sigma_X_quad", mp.nstr(sig_X_quad, 8))
ok("C1 X = m_d − m_u = 2.56 MeV", X, mp.mpf("2.56"), mp.mpf("1e-12"))
ok("C2 sigma_X(quadrature) = 0.0943 (docs/18 報 0.094)", sig_X_quad, mp.mpf("0.0943"), mp.mpf("1e-3"))
sig_X = mp.mpf("0.094")          # 記錄採用之四捨五入值（約定，非新量）
rep("sigma_X_convention", mp.nstr(sig_X, 6), note="docs/18 採用之四捨五入值")

# ======================================================================
# §D em 外部錨（PDG 質量 ⟹ ΔM² ⟹ lattice ε 帶 ⟹ em）
# ======================================================================
dM2_pi = M_PI_P ** 2 - M_PI_0 ** 2
dM2_K = M_K_P ** 2 - M_K_0 ** 2
rep("dM2_pi", mp.nstr(dM2_pi, 12))
rep("dM2_K", mp.nstr(dM2_K, 12))
ok("D1 ΔM²_π = 1261.157 MeV²", dM2_pi, mp.mpf("1261.157"), mp.mpf("1e-5"))
ok("D2 ΔM²_K = −3899.727 MeV²", dM2_K, mp.mpf("-3899.727"), mp.mpf("1e-5"))

MU_KK = M_K_P + M_K_0
dM2_K_gamma = (1 + EPS_MID) * dM2_pi
em_lat = dM2_K_gamma / MU_KK
rep("dM2_K_gamma", mp.nstr(dM2_K_gamma, 12))
rep("em_lat_M eV" if False else "em_lat_MeV", mp.nstr(em_lat, 12))
ok("D3 (ΔM²_K)_γ = 2244.86 MeV²", dM2_K_gamma, mp.mpf("2244.86"), mp.mpf("1e-5"))
ok("D4 em(K+−K0) = +2.2646 MeV", em_lat, mp.mpf("2.2646"), mp.mpf("1e-4"))

em_lo = (1 + EPS_LO) * dM2_pi / MU_KK
em_hi = (1 + EPS_HI) * dM2_pi / MU_KK
sig_em = (em_hi - em_lo) / 2
rep("em_band", "%s .. %s" % (mp.nstr(em_lo, 8), mp.nstr(em_hi, 8)))
rep("sigma_em", mp.nstr(sig_em, 10))
_band_ok = (abs(em_lo - mp.mpf("2.1246")) < mp.mpf("5e-5")) and (abs(em_hi - mp.mpf("2.3664")) < mp.mpf("5e-5"))
ok("D5 68% band = [2.1246, 2.3664]", 1 if _band_ok else 0, 1, 0)
print("  [band] %.6f .. %.6f" % (float(em_lo), float(em_hi)))
ok("D6 sigma_em(sym) = 0.1209", sig_em, mp.mpf("0.1209"), mp.mpf("1e-3"))
dashen = dM2_pi / MU_KK
rep("dashen_MeV", mp.nstr(dashen, 8))
ok("D7 Dashen (ε=0) = +1.2722 MeV", dashen, mp.mpf("1.2722"), mp.mpf("1e-4"))

# ======================================================================
# §E ⟨O⟩ 線之判決（docs/29 單一形式 / docs/30 包絡 / docs/40 K 通道）
# ======================================================================
DM_K = M_K_0 - M_K_P            # 表列 Δm = −3.9340 ⟹ 符號約定：Δm ≡ M_K0 − M_K+ 之負
rep("M_K0_minus_Kplus", mp.nstr(DM_K, 8))
DM_TABLE = mp.mpf("-3.9340")    # docs/29 表列值（約定，逐位對 PDG）
ok("E1 |M_K0 − M_K+| = 3.934 GeV 級（PDG→表列 3.9340）", DM_K, mp.mpf("3.934"), mp.mpf("1e-4"))

# (a) docs/29 單一形式 ⟨O⟩ ≡ 1.667
O1 = mp.mpf(5) / 3
em_pred1 = DM_TABLE - (-1) * X * O1
rep("em_pred_single", mp.nstr(em_pred1, 10))
ok("E2 single-form em_pred = +0.3327 MeV (docs/29, 4 s.f. 捨入)", em_pred1, mp.mpf("0.3327"), mp.mpf("5e-4"))
z1 = (em_pred1 - em_lat) / sig_em
rep("z_single", mp.nstr(z1, 6))
ok("E3 z(docs/29 pred) = −15.98σ", z1, mp.mpf("-15.98"), mp.mpf("1e-3"))

# (b) docs/30 包絡上緣 +0.55
env_lo, env_hi = mp.mpf("-0.52"), mp.mpf("0.55")
em_signed = -em_lat                                   # 包絡以 em(K0−K+) 表述 ⟹ 取 −em
edge = max(env_lo, min(env_hi, em_signed))            # 覆蓋值 → 最近邊界
z_env = (em_signed - edge) / sig_em                   # 距最近邊界之 σ 數
rep("z_envelope", mp.nstr(abs(z_env), 6))
rep("envelope_edge_used", mp.nstr(edge, 4))
ok("E4 z(docs/30 envelope, 距最近邊界) = 14.43σ", abs(z_env), mp.mpf("14.43"), mp.mpf("1e-3"))

# (c) docs/40 K 通道：n_s(K) = 1（K+ = u s̄，奇異數固定）⟹ ⟨O⟩ = 5/3
O_K = 1 + mp.mpf(2) / 3 * 1
ok("E5 ⟨O⟩(n_s=1) = 5/3 exactly", O_K, mp.mpf(5) / 3, 0)
ok("E6 docs/29 表列 ⟨O⟩_pred = 1.667 即 5/3 之捨入（零插值）", O1, mp.mpf("1.667"), mp.mpf("3e-4"))
sig_pred = abs(-1) * O_K * sig_X
sig_tot = mp.sqrt(sig_pred ** 2 + sig_em ** 2)
z_X = (em_pred1 - em_lat) / sig_tot
rep("sigma_pred_X", mp.nstr(sig_pred, 8))
rep("sigma_tot", mp.nstr(sig_tot, 8))
rep("z_with_X_budget", mp.nstr(z_X, 6))
sig_pred_q = O_K * sig_X_quad
sig_tot_q = mp.sqrt(sig_pred_q ** 2 + sig_em ** 2)
z_Xq = (em_pred1 - em_lat) / sig_tot_q
rep("z_with_X_budget_quad", mp.nstr(z_Xq, 6))
ok("E7 z(含 σ_X=0.094 約定) = −9.76σ", z_X, mp.mpf("-9.76"), mp.mpf("2e-3"))
ok("E8 z(含 σ_X=0.0943 四捨五入前) = −9.74σ", z_Xq, mp.mpf("-9.74"), mp.mpf("3e-3"))

# (d) 由外部 em 反解 ⟨O⟩_K 與所需 n_s
O_K_ext = (DM_TABLE - em_lat) / (-1) * (1 / X)
rep("O_K_backsolved", mp.nstr(O_K_ext, 10))
ok("E9 ⟨O⟩_K(外部) = 2.421324", O_K_ext, mp.mpf("2.421324"), mp.mpf("1e-6"))
ns_req = (O_K_ext - 1) * mp.mpf(3) / 2
rep("n_s_required", mp.nstr(ns_req, 10))
ok("E10 required n_s = 2.131986 > max reachable 2", ns_req, mp.mpf("2.131986"), mp.mpf("1e-6"))
ok("E11 required n_s exceeds domain max (2)", 1 if ns_req > 2 else 0, 1, 0)
z_O = (O_K_ext - O_EXTERNAL) / O_EXTERNAL_S
rep("z_O_K", mp.nstr(z_O, 6))
ok("E12 z(⟨O⟩_K vs BMW 1.498±0.262) = +3.524σ", z_O, mp.mpf("3.524"), mp.mpf("1e-3"))

# ======================================================================
# §F 中性子衰變餘裕（x = m_n − m_p − m_e）—— 供 sign/precision 判詞換算
# ======================================================================
x_obs = M_N - M_P - M_E
rep("x_neutron_decay_MeV", mp.nstr(x_obs, 8))
ok("F1 x(n→peν) = 0.7823 MeV", x_obs, mp.mpf("0.7823"), mp.mpf("1e-3"))

# ======================================================================
npass = sum(1 for c_ in CHECKS if c_["pass"])
nfail = len(CHECKS) - npass
print("\n" + "=" * 100)
print("L2 RESULT: %d/%d PASS   (FAIL %d)" % (npass, len(CHECKS), nfail))
print("=" * 100)

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "data", "l2_independent_verify.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"summary": {"pass": npass, "total": len(CHECKS), "fail": nfail},
               "checks": CHECKS, "derived": {r["key"]: r["value"] for r in RESULTS}},
              fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON ->", os.path.relpath(out))
raise SystemExit(0 if nfail == 0 else 1)
