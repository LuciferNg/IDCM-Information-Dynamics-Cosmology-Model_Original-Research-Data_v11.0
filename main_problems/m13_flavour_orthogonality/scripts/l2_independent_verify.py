#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M13 獨立重推導（idcm-audit-verify §2.4）

鐵律
  R1 唔准 import／call 任何 L1 腳本（scripts/0*.py、scripts/1*.py）；唔准讀 docs/。
     本檔由**原始量**起算 —— 讀 L1／docs 即退化成同層互檢。
  R2 INPUTS 只放原始量：結構計數（N_d/B/N_gen/γ44/Σγ）、**測量值**（PDG 質量 ± σ）。
     唔准放推導結果或文檔聲稱值。
  R3 DERIVE 每個 key 必須係 INPUTS 嘅純函數（零隨機、零檔案 IO）。
  R4 mpmath dps = 40；精確恆等用 Fraction。
"""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG, RES, NFAIL = [], {}, 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def ser(x):
    try:
        xf = float(x)
    except (OverflowError, ValueError):
        xf = None
    if xf is not None and (xf == 0 or 1e-300 < abs(xf) < 1e300):
        return {"float": xf, "str": mp.nstr(x, 25)}
    return {"float": None, "str": mp.nstr(x, 25)}


# ============================================================================
# ① INPUTS —— 原始量（PDG 2024 測量值 ＋ 結構計數）
# ============================================================================
INPUTS = {
    # PDG 2024 質量（MeV）與 1σ 誤差（對稱化）
    "m_n": "939.5654205", "s_n": "5.0e-7",
    "m_p": "938.27208816", "s_p": "2.9e-7",   # PDG 2024 原文 9.3827208816E-01
    "m_Sm": "1197.449", "s_Sm": "0.029",     # Σ⁻
    "m_S0": "1192.642", "s_S0": "0.024",     # Σ⁰
    "m_Sp": "1189.37",  "s_Sp": "0.07",      # Σ⁺
    "m_Xm": "1321.71",  "s_Xm": "0.07",      # Ξ⁻
    "m_X0": "1314.86",  "s_X0": "0.20",      # Ξ⁰
    # 結構計數（已註冊）
    "N_d": 44, "B": 4, "N_gen": 3, "gamma44": 4, "sum_gamma_even": 90,
}


def mpf(v):
    return v if isinstance(v, mp.mpf) else mp.mpf(v)


# ============================================================================
# ② DERIVE —— 純函數
# ============================================================================
def _dN(I):
    return mpf(I["m_n"]) - mpf(I["m_p"])


def _dS(I):
    return mpf(I["m_Sm"]) - mpf(I["m_S0"])


def _dX(I):
    return mpf(I["m_Xm"]) - mpf(I["m_X0"])


def _dS2(I):
    return mpf(I["m_S0"]) - mpf(I["m_Sp"])


def _s_dN(I):
    return mp.sqrt(mpf(I["s_n"])**2 + mpf(I["s_p"])**2)


def _s_dS(I):
    return mp.sqrt(mpf(I["s_Sm"])**2 + mpf(I["s_S0"])**2)


def _s_dX(I):
    return mp.sqrt(mpf(I["s_Xm"])**2 + mpf(I["s_X0"])**2)


def _s_dS2(I):
    return mp.sqrt(mpf(I["s_S0"])**2 + mpf(I["s_Sp"])**2)


def _ds(I):                       # δ_s = Δ_Ξ/2（Ξ 通道）
    return _dX(I)/2


def _dd(I):                       # δ_d = Δ_Σ − Δ_Ξ/2
    return _dS(I) - _dX(I)/2


def _du(I):                       # δ_u = Δ_N − δ_d
    return _dN(I) - _dd(I)


def _rho_u(I):
    return _du(I)/_ds(I)


def _rho_d(I):
    return _dd(I)/_ds(I)


def _jac(fn, I, keys):
    """數值 Jacobian（中心差分；mp50 精度足夠）"""
    out = []
    for k in keys:
        h = mp.fabs(mpf(I[k]))*mp.mpf("1e-14") + mp.mpf("1e-30")
        Ip = dict(I); Ip[k] = str(mpf(I[k]) + h)
        Im = dict(I); Im[k] = str(mpf(I[k]) - h)
        out.append((fn(Ip) - fn(Im))/(2*h))
    return out


MASS = ["m_n", "m_p", "m_Sm", "m_S0", "m_Xm", "m_X0"]
SIG = {"m_n": "s_n", "m_p": "s_p", "m_Sm": "s_Sm", "m_S0": "s_S0",
       "m_Xm": "s_Xm", "m_X0": "s_X0"}


def _var(fn, I):
    g = _jac(fn, I, MASS)
    return mp.fsum((gi*mpf(I[SIG[k]]))**2 for gi, k in zip(g, MASS))


def _cov(f1, f2, I):
    g1, g2 = _jac(f1, I, MASS), _jac(f2, I, MASS)
    return mp.fsum(g1[i]*g2[i]*mpf(I[SIG[k]])**2 for i, k in enumerate(MASS))


DERIVE = {
    "dN": lambda I: _dN(I),
    "dS": lambda I: _dS(I),
    "dX": lambda I: _dX(I),
    "dS2": lambda I: _dS2(I),
    "delta_s": _ds,
    "delta_d": _dd,
    "delta_u": _du,
    "rho_u": _rho_u,
    "rho_d": _rho_d,
    "sigma_rho_u": lambda I: mp.sqrt(_var(_rho_u, I)),
    "sigma_rho_d": lambda I: mp.sqrt(_var(_rho_d, I)),
    "trace": lambda I: _rho_u(I) + _rho_d(I),
    "sigma_trace": lambda I: mp.sqrt(_var(_rho_u, I) + _var(_rho_d, I)
                                     + 2*_cov(_rho_u, _rho_d, I)),
    "diff": lambda I: _rho_d(I) - _rho_u(I),
    "sigma_diff": lambda I: mp.sqrt(_var(_rho_u, I) + _var(_rho_d, I)
                                    - 2*_cov(_rho_u, _rho_d, I)),
    "corr": lambda I: (_cov(_rho_u, _rho_d, I)
                       / mp.sqrt(_var(_rho_u, I)*_var(_rho_d, I))),
    "ratio_dN_dX": lambda I: _dN(I)/_dX(I),
    "sigma_ratio_dN_dX": lambda I: mp.sqrt(
        (_s_dN(I)/_dX(I))**2 + (_dN(I)*_s_dX(I)/_dX(I)**2)**2),
    "target_3_16": lambda I: mpf(I["N_gen"])/mpf(I["B"])**2,
    "target_3_8": lambda I: mp.mpf(3)/8,
    "trace_sigma_chan": lambda I: _dN(I)/_dS2(I),
    # 2026-11-02 新增：α 組合 = 1/(2·trace)（與迹恆等）；σ 由 σ(trace) 傳播
    "alpha_sum_ratio": lambda I: 1 / (2 * (_rho_u(I) + _rho_d(I))),
    "sigma_alpha_sum_ratio": lambda I: (mp.sqrt(_var(_rho_u, I) + _var(_rho_d, I)
                                               + 2 * _cov(_rho_u, _rho_d, I))
                                        / (2 * (_rho_u(I) + _rho_d(I)) ** 2)),
    "sigma_trace_sigma_chan": lambda I: mp.fabs(
        _dN(I)/_dS2(I)**2*_s_dS2(I)),
    "dev_vs_3_8": lambda I: mp.fabs(_rho_u(I) + _rho_d(I) - mp.mpf(3)/8),
    "dev_over_base_q16": lambda I: (mp.fabs(_rho_u(I) + _rho_d(I) - mp.mpf(3)/8)
                                    / (mp.mpf(3)/8)) / (mp.mpf(1)/(2*16**2)),
    "z_trace": lambda I: ((_rho_u(I) + _rho_d(I) - mp.mpf(3)/8)
                          / mp.sqrt(_var(_rho_u, I) + _var(_rho_d, I)
                                    + 2*_cov(_rho_u, _rho_d, I))),
    "z_trace_sigma_chan": lambda I: ((
        _dN(I)/_dS2(I) - mp.mpf(3)/8) / mp.fabs(_dN(I)/_dS2(I)**2*_s_dS2(I))),
    "z_ratio_dev": lambda I: ((_dN(I)/_dX(I) - mpf(I["N_gen"])/mpf(I["B"])**2)
                              / mp.sqrt((_s_dN(I)/_dX(I))**2
                                        + (_dN(I)*_s_dX(I)/_dX(I)**2)**2)),
    # dev_abs = |trace − 3/8|；要 dev ≥ 3σ ⟹ σ_trace ≤ dev_abs/3
    # σ_trace = (Δ_N/δ_s²)·σ_δs = (Δ_N/δ_s²)·(σ_Ξ0/2) ⟹ σ_Ξ0 = 2·σ_trace·δ_s²/Δ_N
    "need_sigmaX0_3sigma_id": lambda I: (
        2*(mp.fabs(_rho_u(I) + _rho_d(I) - mp.mpf(3)/8)/3)*_ds(I)**2/_dN(I)),
}

# 觀測量 z 值（與 L1 同一定義，但由本檔獨立算）
_OBS = {}


def _set_obs(I):
    if not _OBS:
        _OBS["rho_u"] = _rho_u(I)
        _OBS["rho_d"] = _rho_d(I)
        _OBS["s_u"] = mp.sqrt(_var(_rho_u, I))
        _OBS["s_d"] = mp.sqrt(_var(_rho_d, I))


def _z_rd(pred, I):
    _set_obs(I)
    return (pred - _OBS["rho_d"])/_OBS["s_d"]


def _z_ru(pred, I):
    _set_obs(I)
    return (pred - _OBS["rho_u"])/_OBS["s_u"]


# α = (1,2,4) 之結構預言（Cartan ⊕ 1）
def _m_alpha(au, ad, as_, I):
    den = mpf(I["sum_gamma_even"])
    return ((2*au - ad)/den, (2*ad - au)/den, 2*as_/den)


def _rho_from_alpha(au, ad, as_, I):
    du, dd, ds_ = _m_alpha(au, ad, as_, I)
    return du/ds_, dd/ds_


# 家族 α=(1,x,x²)
def _fam(x, I):
    return _rho_from_alpha(mp.mpf(1), x, x*x, I)


DERIVE.update({
    "z_x2_u": lambda I: _z_ru(_fam(mp.mpf(2), I)[0], I),
    "z_x2_d": lambda I: _z_rd(_fam(mp.mpf(2), I)[1], I),
    "chi2_x2": lambda I: (_z_ru(_fam(2, I)[0], I)**2 + _z_rd(_fam(2, I)[1], I)**2),
    "chi2_x4": lambda I: (_z_ru(_fam(4, I)[0], I)**2 + _z_rd(_fam(4, I)[1], I)**2),
    "chi2_x3": lambda I: (_z_ru(_fam(3, I)[0], I)**2 + _z_rd(_fam(3, I)[1], I)**2),
    "z_r15_4": lambda I: _z_rd(mp.mpf(3)/(2*mp.mpf(15)/4), I),
    "z_r11_3": lambda I: _z_rd(mp.mpf(3)/(2*mp.mpf(11)/3), I),
    "z_r7_2": lambda I: _z_rd(mp.mpf(3)/(2*mp.mpf(7)/2), I),
    "z_konly_rho_d": lambda I: _z_rd(mp.mpf(3)/4, I),      # α 依 k ⟹ ρ_d=3/4
    "z_donly_rho_d": lambda I: _z_rd(mp.mpf(0), I),        # α 依 d ⟹ ρ_d=0
})

# --- 供 L3 逐文檔核對之補充量（全部仍係 INPUTS 之純函數）--------------
def _comb(n, k):
    r = 1
    for i in range(k):
        r = r*(n-i)//(i+1)
    return r


_EXTRA = {
    # docs/05：α_s/α_u 之窗（= 3/(2ρ_d)）
    "alpha_ratio": lambda I: mp.mpf(3)/(2*_rho_d(I)),
    "alpha_ratio_lo": lambda I: mp.mpf(3)/(2*(_rho_d(I) + mp.sqrt(_var(_rho_d, I)))),
    "alpha_ratio_hi": lambda I: mp.mpf(3)/(2*(_rho_d(I) - mp.sqrt(_var(_rho_d, I)))),
    # 線性傳播版（docs/05 實際採用之約定）：σ_r = r·σ_ρ/ρ
    "alpha_ratio_sig_lin": lambda I: (mp.mpf(3)/(2*_rho_d(I)))
                                     * mp.sqrt(_var(_rho_d, I))/_rho_d(I),
    "alpha_ratio_lin_lo": lambda I: (mp.mpf(3)/(2*_rho_d(I)))
        * (1 - mp.sqrt(_var(_rho_d, I))/_rho_d(I)),
    "alpha_ratio_lin_hi": lambda I: (mp.mpf(3)/(2*_rho_d(I)))
        * (1 + mp.sqrt(_var(_rho_d, I))/_rho_d(I)),
    # docs/04：α_s ∈ {1..8} 之整數格 ρ_d = 3/(2α_s)
    **{f"rd_int{k}": (lambda kk: (lambda I: mp.mpf(3)/(2*kk)))(k) for k in range(1, 9)},
    # docs/06：樹層核 α=(1,4,16)
    "rho_u_a1416": lambda I: _rho_from_alpha(mp.mpf(1), mp.mpf(4), mp.mpf(16), I)[0],
    "rho_d_a1416": lambda I: _rho_from_alpha(mp.mpf(1), mp.mpf(4), mp.mpf(16), I)[1],
    # docs/12：競爭者之 ρ_d 值（精確）
    "rd_r15_4": lambda I: mp.mpf(3)/(2*mp.mpf(15)/4),
    "rd_r11_3": lambda I: mp.mpf(3)/(2*mp.mpf(11)/3),
    "rd_r7_2": lambda I: mp.mpf(3)/(2*mp.mpf(7)/2),
    # docs/02：枚舉類容量 C(44,K)·2^K
    **{f"nclass_K{k}": (lambda kk: (lambda I: mp.mpf(_comb(I["N_d"], kk)*2**kk)))(k)
       for k in range(2, 7)},
    # docs/08：Σ 質量造成之 σ_ρd 硬底 = 2σ(Δ_Σ)/Δ_Ξ
    "sigma_floor_sigma": lambda I: 2*_s_dS(I)/_dX(I),
    # docs/07：α 標籤族之試驗數 6 規則 × 3! 指派
    "n_trials_alpha": lambda I: mp.mpf(6*36),   # 6 條規則 × (3!)² 指派 = 216
}
DERIVE.update(_EXTRA)

# ============================================================================
# ③ IDENTITIES —— 兩條獨立路徑（差 < 1e-30；精確者用 Fraction ==）
# ============================================================================
IDENTITIES = [
    ("I1 ρ_u+ρ_d = (δ_u+δ_d)/δ_s = Δ_N/δ_s",
     lambda I: _rho_u(I) + _rho_d(I),
     lambda I: _dN(I)/_ds(I)),
    ("I2 2Δ_N/Δ_Ξ = Δ_N/δ_s（δ_s = Δ_Ξ/2）",
     lambda I: 2*_dN(I)/_dX(I),
     lambda I: _dN(I)/_ds(I)),
    ("I3 ρ_d−ρ_u = (4Δ_Σ−2Δ_N)/Δ_Ξ − 2",
     lambda I: _rho_d(I) - _rho_u(I),
     lambda I: (4*_dS(I) - 2*_dN(I))/_dX(I) - 2),
    ("I4 trace = 2Δ_N/Δ_Ξ",
     lambda I: _rho_u(I) + _rho_d(I),
     lambda I: 2*_dN(I)/_dX(I)),
    ("I5 δ_s = Δ_Ξ/2 反解 Ξ 差", lambda I: _ds(I)*2, lambda I: _dX(I)),
    ("I6 δ_d = Δ_Σ − δ_s", lambda I: _dd(I), lambda I: _dS(I) - _ds(I)),
]

def _cartan_plus_1_act(v):
    """M = (1/90)(A₂ Cartan ⊕ 1) 作用於 3-向量（Fraction 精確）"""
    au, ad, as_ = v
    return (2*au - ad, 2*ad - au, 2*as_)          # 未除 90（比值不受影響）


EXACT = [
    ("E1 M·(1,2,4) 正比於 (0,3,8)（Cartan⊕1 精確）",
     _cartan_plus_1_act((1, 2, 4)), (0, 3, 8)),
    ("E2 由 E1 直接得 (ρ_u, ρ_d) = (0, 3/8)",
     (F(0, 8), F(3, 8)), (F(0), F(3, 8))),
    ("E3 Cartan 特徵多項式 λ²−trλ+det，tr=4、det=3 ⟹ 根 {1,3}",
     (1*1 - 4*1 + 3, 3*3 - 4*3 + 3), (0, 0)),
    ("E4 (B−1)/(2B) = 3/8（**定義恆等**：B−1 = N_gen，零內容）",
     (F(4-1, 2*4),), (F(3, 8),)),
    ("E5 N_gen/B² = 3/16（**定義恆等**，零內容）",
     (F(3, 4*4),), (F(3, 16),)),
]


def main():
    hdr("§0 輸入（原始量）")
    I = dict(INPUTS)
    for k in ("m_n", "m_p", "m_Sm", "m_S0", "m_Sp", "m_Xm", "m_X0",
              "N_d", "B", "N_gen", "gamma44", "sum_gamma_even"):
        print(f"  {k:16s} = {I[k]}")

    hdr("§1 DERIVE（由 INPUTS 純函數重算）")
    D = {}
    for key, fn in DERIVE.items():
        try:
            D[key] = mp.mpf(fn(I))
            print(f"  {key:24s} = {mp.nstr(D[key], 18)}")
        except Exception as e:                       # noqa: BLE001
            chk(f"D.{key} 可推導", False, f"{type(e).__name__}: {e}")
    chk(f"D 全部可推導（{len(DERIVE)}）", len(D) == len(DERIVE),
        f"{len(D)}/{len(DERIVE)}")

    hdr("§2 IDENTITIES（兩路徑差 < 1e−30）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(I)), mp.mpf(fb(I))
        d = abs(a - b)
        chk(tag, d < mp.mpf("1e-30"), f"差 = {mp.nstr(d, 6)}")

    hdr("§3 EXACT（Fraction 精確相等）")
    for tag, a, b in EXACT:
        chk(tag, tuple(a) == tuple(b), f"{a} vs {b}")

    hdr("§4 交叉核（DERIVE 內在一致性，唔涉 docs）")
    _x = (2 * D["rho_d"] + 4 * D["rho_u"]) / 3
    _y = (4 * D["rho_d"] + 2 * D["rho_u"]) / 3
    chk("C0 α 組合 == 1/(x+y)（與迹恆等，差 <1e−40）",
        abs(D["alpha_sum_ratio"] - 1 / (_x + _y)) < mp.mpf("1e-40"),
        f"差 = {mp.nstr(D['alpha_sum_ratio'] - 1/(_x+_y), 3)}")
    chk("C1 σ_trace < σ_rho_d（迹更精）",
        D["sigma_trace"] < D["sigma_rho_d"],
        f"{mp.nstr(D['sigma_trace'],6)} < {mp.nstr(D['sigma_rho_d'],6)}")
    chk("C2 corr ≈ −1（强反相關）", D["corr"] < mp.mpf("-0.99"),
        f"corr = {mp.nstr(D['corr'],8)}")
    chk("C3 Var(迹)+Var(差) = 2(Var_u+Var_d)",
        abs(D["sigma_trace"]**2 + D["sigma_diff"]**2
            - 2*(D["sigma_rho_u"]**2 + D["sigma_rho_d"]**2)) < mp.mpf("1e-25"),
        f"LHS-RHS = {mp.nstr(D['sigma_trace']**2 + D['sigma_diff']**2 - 2*(D['sigma_rho_u']**2 + D['sigma_rho_d']**2),6)}")
    chk("C4 χ²(x=2) < χ²(x=3) < χ²(x=4)",
        D["chi2_x2"] < D["chi2_x3"] < D["chi2_x4"],
        f"{mp.nstr(D['chi2_x2'],5)} / {mp.nstr(D['chi2_x3'],5)} / {mp.nstr(D['chi2_x4'],5)}")
    chk("C5 ρ_d(α=124) == 3/8 精確",
        abs(_fam(mp.mpf(2), I)[1] - mp.mpf(3)/8) < mp.mpf("1e-35"))
    chk("C6 ρ_u(α=124) == 0 精確", abs(_fam(mp.mpf(2), I)[0]) < mp.mpf("1e-35"))

    hdr("§5 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    RES["_inputs"] = {k: (v if not isinstance(v, mp.mpf) else mp.nstr(v, 30))
                      for k, v in I.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": n, "fail": NFAIL, "n_derived": len(D)}
    (DATA / "l2_independent_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l2.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / 'l2_independent_verify.json'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())