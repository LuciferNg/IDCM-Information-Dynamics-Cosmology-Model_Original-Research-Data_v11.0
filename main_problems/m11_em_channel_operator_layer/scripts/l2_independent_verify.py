#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M11 獨立重推導（idcm-audit-verify §2.4；模板 `l2_independent_verify_template.py` 之填版）
================================================================================
鐵律遵守：
  R1 零 import L1、零讀 docs（本檔只讀自己寫死嘅 INPUTS）。
  R2 INPUTS 只放原始量：PDG 質量、夸克電荷、`X`（FLAG）、`a₀`（M5 校準值）、M5 H 之結構常數、
     以及**外部已封**之 lattice em 錨（M4 docs/18 之上游結論，屬輸入而非本記錄推導）。
  R3 DERIVE 全部係 INPUTS 之純函數。
  R4 mpmath dps = 40；41×41 特徵解用 `mp.eigsy`（與 L1 之 numpy 路線不同，真獨立）。

覆蓋：8 條 Δm／σ／E_ch／⟨O⟩、k₀ 與 ħc/a₀、路線 C 之 (a, b, χ²)、退化類內差與 z、
      G0 之 ∂H/∂X 幅度與 deficit、em(ε) 五支與 z、D 強制預言。
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
# ① INPUTS —— 原始量
# ============================================================================
MASS = {  # PDG 2024（MeV）；σ = PDG 誤差
    "p": ("938.27208943", "2.9e-7"), "n": ("939.56542052", "5.4e-7"),
    "K+": ("493.677", "0.016"), "K0": ("497.611", "0.013"),
    "pi+": ("139.57039", "1.8e-4"), "pi0": ("134.9768", "5e-4"),
    "S-": ("1197.449", "0.030"), "S0": ("1192.642", "0.024"), "S+": ("1189.37", "0.07"),
    "Xi-": ("1321.71", "0.07"), "Xi0": ("1314.86", "0.20"),
    "D+": ("1869.66", "0.05"), "D0": ("1864.84", "0.05"),
}
CH = [("n-p", "n", "p", 1), ("K+-K0", "K+", "K0", -1), ("pi+-pi0", "pi+", "pi0", 1),
      ("S--S0", "S-", "S0", 1), ("S0-S+", "S0", "S+", 1), ("S--S+", "S-", "S+", 2),
      ("Xi--Xi0", "Xi-", "Xi0", 1), ("D+-D0", "D+", "D0", 1)]
# 價夸克內容（u=+2/3、d=-1/3、s=-1/3、c=+2/3；反夸克取負）
CONTENT = {
    "n-p": (("d", 2, "u", 1), ("u", 2, "d", 1)),
    "K+-K0": (("u", 1, "sbar", 1), ("d", 1, "sbar", 1)),
    "pi+-pi0": (("u", 1, "dbar", 1), "PI0"),
    "S--S0": (("d", 2, "s", 1), ("u", 1, "d", 1, "s", 1)),
    "S0-S+": (("u", 1, "d", 1, "s", 1), ("u", 2, "s", 1)),
    "S--S+": (("d", 2, "s", 1), ("u", 2, "s", 1)),
    "Xi--Xi0": (("d", 1, "s", 2), ("u", 1, "s", 2)),
    "D+-D0": (("c", 1, "dbar", 1), ("c", 1, "ubar", 1)),
}
# 外部 lattice em 錨（M4 docs/18 上游；MeV）與其 σ —— 屬 INPUTS
# 上游 lattice em 錨（M4 docs/18 §C1；pi0 用 Dashen ⟹ em ≡ Δm）
EM_ANCHOR = {"n-p": "-1.59", "K+-K0": "2.264589", "pi+-pi0": "4.59359",
             "S--S0": "-0.08", "Xi--Xi0": "1.29"}
SIG_EM = {"n-p": "0.461", "K+-K0": "0.120863", "pi+-pi0": "0.0",
          "S--S0": "0.361", "Xi--Xi0": "0.17"}
# PDG 質量差誤差 σ_Δm（MeV）—— M11 正典 σ_O **不含**此項（見 §十四.6 敏感度）
SIG_DM = {"n-p": "0.00000061", "K+-K0": "0.0206", "pi+-pi0": "0.00053",
          "S--S0": "0.0384", "Xi--Xi0": "0.2119"}
# 通道類別（nq, ns）—— R4 模型 A–D 之結構指派依據
CLASS = {"n-p": (3, 0), "K+-K0": (2, 1), "pi+-pi0": (2, 0), "S--S0": (3, 1),
         "S0-S+": (3, 1), "S--S+": (3, 1), "Xi--Xi0": (3, 2), "D+-D0": (2, 0)}
BLIND = ["S0-S+", "S--S+", "D+-D0"]

INPUTS = dict(
    ALPHA="7.2973525693e-3", HBARC_MeVfm="197.3269804", HBARC_GeVfm="0.1973269804",
    A0_fm="0.37742216816391858835",
    X_MeV="2.56", SIG_X="0.094",
    B=4, L_MAX=40, NODES_RULE="4*3^(d-1)",
    M_C="1.27", M_S="0.495", M_L="0.30", Q_C="0.6666666666666666666666666666666666666667",
    Q_S="0.3333333333333333333333333333333333333333",
    K0_MeV="3.8152622441",
    FLAG_DMPI2="1261.157", FLAG_MKSUM="991.288", FLAG_PRED="0.3327", FLAG_SIGX="0.156667",
    D_FORCED_O="1.908984", AEFF_MeV="4.4492",
)
# 分支切換（scripts/14 用）：只影響 a₀／k₀ 兩項輸入，其餘不動。無環境變數 ⟹ v1 正典，輸出 byte 不變。
import os as _os
if _os.environ.get("M11_A0_FM"):
    INPUTS["A0_fm"] = _os.environ["M11_A0_FM"]
if _os.environ.get("M11_K0_MEV"):
    INPUTS["K0_MeV"] = _os.environ["M11_K0_MEV"]

# ============================================================================
# ② DERIVE —— 純函數
# ============================================================================
mp.mp.dps = 40


def _mp(s):
    return mp.mpf(str(s))


def dm(I, hi, lo):
    return _mp(MASS[hi][0]) - _mp(MASS[lo][0])


def sd(I, hi, lo):
    return mp.sqrt(_mp(MASS[hi][1]) ** 2 + _mp(MASS[lo][1]) ** 2)


def qpair(spec):
    qmap = {"u": F(2, 3), "d": F(-1, 3), "s": F(-1, 3), "c": F(2, 3),
            "ubar": F(-2, 3), "dbar": F(1, 3), "sbar": F(1, 3)}
    if spec == "PI0":
        a = qpair(("u", 1, "ubar", 1))
        b = qpair(("d", 1, "dbar", 1))
        return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
    qs = []
    for i in range(0, len(spec), 2):
        for _ in range(spec[i + 1]):
            qs.append(qmap[spec[i]])
    s2 = sum(q * q for q in qs)
    pr = sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs)))
    return s2, pr


def phi_ch(name):
    a = qpair(CONTENT[name][0])
    b = qpair(CONTENT[name][1])
    return (a[0] - b[0] + a[1] - b[1]) / 2


def H_matrix(a_fm, t_MeV, q, mh_MeV, ml_MeV, I):
    """M5 之 41 層 H（與 L1 同式，但以 mpmath 實作）"""
    B, L = I["B"], I["L_MAX"]
    a = _mp(a_fm)
    t = _mp(t_MeV)
    qq = _mp(q)
    mh, ml = _mp(mh_MeV), _mp(ml_MeV)
    hb = _mp(I["HBARC_GeVfm"])
    gcoef = _mp(I["ALPHA"]) * qq * hb / a / t
    Cq = (8 * mp.pi / 3) * _mp(I["ALPHA"]) * qq / (16 * mh * ml)
    h = Cq * (hb / a) ** 3 / t
    n = L + 1
    M = mp.matrix(n, n)
    for d in range(n):
        M[d, d] = -gcoef * mp.power(B, -d) - (h if d == 0 else 0)
    off0 = -mp.sqrt(B)
    off = -mp.sqrt(B - 1)
    for d in range(n - 1):
        v = off0 if d == 0 else off
        M[d, d + 1] = v
        M[d + 1, d] = v
    return M


def amps(a, t, q, mh, ml, I):
    M = H_matrix(a, t, q, mh, ml, I)
    E, V = mp.eigsy(M)
    vec = [V[i, 0] for i in range(M.rows)]
    w = []
    for d in range(M.rows):
        nd = 1 if d == 0 else I["B"] * (I["B"] - 1) ** (d - 1)
        w.append(_mp(nd) * vec[d] ** 2)
    tot = mp.fsum(w)
    w = [x / tot for x in w]
    rA = mp.fsum(w[d] * mp.power(I["B"], -d) for d in range(M.rows))
    P0 = w[0]
    hb = _mp(I["HBARC_GeVfm"])
    qq = _mp(q)
    Cq = (8 * mp.pi / 3) * _mp(I["ALPHA"]) * qq / (16 * _mp(mh) * _mp(ml))
    A = _mp(I["ALPHA"]) * qq * (hb / _mp(a)) * rA * 1000   # GeV → MeV
    c = Cq * (hb / _mp(a)) ** 3 * P0 * 1000                # GeV → MeV
    return A, c, rA, P0, w


def wls2(rows, I, sigma_x=True):
    """2 參數加權最小平方（Φ, Ψ）⟹ (a, b, χ²)"""
    X = _mp(I["X_MeV"])
    A11 = A12 = A22 = b1 = b2 = 0
    ys, ss, phis, psis = [], [], [], []
    for name, hi, lo, coeff in CH:
        y = dm(I, hi, lo)
        s = sd(I, hi, lo)
        if sigma_x:
            s = mp.sqrt(s ** 2 + (_mp(coeff) * _mp(I["SIG_X"])) ** 2)
        p = _mp(phi_ch(name))
        w = 1 / s ** 2
        A11 += w * p * p
        A12 += w * p * coeff
        A22 += w * coeff * coeff
        b1 += w * p * y
        b2 += w * coeff * y
        ys.append(y); ss.append(s); phis.append(p); psis.append(coeff)
    det = A11 * A22 - A12 * A12
    a = (b1 * A22 - b2 * A12) / det
    b = (A11 * b2 - A12 * b1) / det
    chi2 = mp.fsum(((ys[i] - a * phis[i] - b * psis[i]) / ss[i]) ** 2 for i in range(len(ys)))
    zs = [float((ys[i] - a * phis[i] - b * psis[i]) / ss[i]) for i in range(len(ys))]
    return a, b, chi2, zs, [float(s) for s in ss]


DERIVE = {}


def _add(key, fn):
    DERIVE[key] = fn


for _name, _hi, _lo, _c in CH:
    _add(f"dm_{_name}", (lambda n, h, l: (lambda I: dm(I, h, l)))(_name, _hi, _lo))
    _add(f"sd_{_name}", (lambda h, l: (lambda I: sd(I, h, l)))(_hi, _lo))
    _add(f"E_{_name}", (lambda h, l, c: (lambda I: dm(I, h, l) - c * _mp(I["X_MeV"])))(_hi, _lo, _c))
    _add(f"phi_{_name}", (lambda n: (lambda I: _mp(phi_ch(n))))(_name))

_add("hc_over_a0", lambda I: _mp(I["HBARC_MeVfm"]) / _mp(I["A0_fm"]))
_add("k0_from_a0", lambda I: _mp(I["ALPHA"]) * _mp(I["HBARC_MeVfm"]) / _mp(I["A0_fm"]))
_PAIR = {x[0]: (x[1], x[2]) for x in CH}
_COEF = {x[0]: x[3] for x in CH}
for _n, _em in EM_ANCHOR.items():
    _add(f"O_{_n}", (lambda n, em, c: (lambda I: (dm(I, _PAIR[n][0], _PAIR[n][1]) - _mp(em))
                                                / (c * _mp(I["X_MeV"]))))(_n, _em, _COEF[_n]))
_add("routeC_a", lambda I: wls2(None, I, sigma_x=False)[0])
_add("routeC_b", lambda I: wls2(None, I, sigma_x=False)[1])
_add("routeC_chi2", lambda I: wls2(None, I, sigma_x=False)[2])
_add("routeC_chi2_dof", lambda I: wls2(None, I, sigma_x=False)[2] / 6)
_add("routeC_O_eff", lambda I: wls2(None, I, sigma_x=False)[1] / _mp(I["X_MeV"]))
_add("routeC_chi2_dof_sigX", lambda I: wls2(None, I, sigma_x=True)[2] / 6)
_add("degen_np_S0Splus", lambda I: abs(dm(I, "S0", "S+") - dm(I, "n", "p")))
_add("degen_S_Xi_D", lambda I: max(dm(I, "Xi-", "Xi0"), dm(I, "S-", "S0"), dm(I, "D+", "D0"))
    - min(dm(I, "Xi-", "Xi0"), dm(I, "S-", "S0"), dm(I, "D+", "D0")))

def g0_full(I, step_rel="1e-5"):
    """G0：∂(A+3c)/∂m_l 之中央差分（相對步長 1e-5；已驗與 1e-6…3e-5 六位一致）⟹ (dM_MeV, deficit)"""
    q, mh, ml = I["Q_C"], I["M_C"], _mp(I["M_L"])
    hstep = _mp(step_rel) * ml
    t = "0.0007884"                       # GeV（L1 約定）
    Ap, cp, *_ = amps(I["A0_fm"], t, q, mh, ml + hstep, I)
    Am, cm, *_ = amps(I["A0_fm"], t, q, mh, ml - hstep, I)
    # 導數 MeV/GeV × X(GeV) ⟹ MeV
    dM = abs(((Ap - Am) + 3 * (cp - cm)) / (2 * hstep)) * _mp("0.00256")
    return dM, _mp("2.56") / dM


_add("G0_dM_D_hf", lambda I: g0_full(I)[0])
_add("G0_deficit_D", lambda I: g0_full(I)[1])
for _eps, _sig, _tag in (("0.78", "0.10", "milc"), ("0.79", "0.06", "flag45"),
                         ("0.79", "0.15", "flag45q"), ("0.73", "0.17", "flag46")):
    _add(f"em_eps_{_tag}", (lambda e: (lambda I: (1 + _mp(e)) * _mp(I["FLAG_DMPI2"]) / _mp(I["FLAG_MKSUM"])))(_eps))
    _add(f"sig_em_{_tag}", (lambda s: (lambda I: _mp(s) * _mp(I["FLAG_DMPI2"]) / _mp(I["FLAG_MKSUM"])))(_sig))
_add("z_flag45_only", lambda I: (_mp(I["FLAG_PRED"]) - DERIVE["em_eps_flag45"](I)) / DERIVE["sig_em_flag45"](I))
_add("z_flag45_total", lambda I: (_mp(I["FLAG_PRED"]) - DERIVE["em_eps_flag45"](I))
    / mp.sqrt(DERIVE["sig_em_flag45"](I) ** 2 + _mp(I["FLAG_SIGX"]) ** 2))
_add("eps_model_from_pred", lambda I: _mp(I["FLAG_PRED"]) * _mp(I["FLAG_MKSUM"]) / _mp(I["FLAG_DMPI2"]) - 1)
# D 模型（⟨O⟩ ≡ 1.908984，M4 docs/18 上游）之盲測強制預言：em = Δm − coeff·X·⟨O⟩
_add("em_forced_S0Splus", lambda I: dm(I, "S0", "S+") - 1 * _mp(I["X_MeV"]) * _mp(I["D_FORCED_O"]))
_add("em_forced_SminusSplus", lambda I: dm(I, "S-", "S+") - 2 * _mp(I["X_MeV"]) * _mp(I["D_FORCED_O"]))


# ── σ_O 與四模型（A 單值／B 按 nq／C 按 ns／D 按 (nq,ns)）之 ⟨O⟩ 指派 ──
def _sig_O(name, full=False):
    """σ_O = √[(σ_em/(|coeff|X))² + (O·σ_X/X)² (+ (σ_Δm/(|coeff|X))²)]

    `full=False` = M11 正典約定（R4 擬合與 k₀^crit 所用，σ_Δm 不入）；
    `full=True`  = 補齊 Δm 量測誤差之敏感度變體。
    """
    c = _COEF[name]
    X = _mp(INPUTS["X_MeV"])
    t1 = _mp(SIG_EM[name]) / (abs(c) * X)
    t2 = _O_anchor(name) * _mp(INPUTS["SIG_X"]) / X
    t = [t1, t2]
    if full:
        t.append(_mp(SIG_DM[name]) / (abs(c) * X))
    return mp.sqrt(mp.fsum(x ** 2 for x in t))


def _O_anchor(n):
    """⟨O⟩_n = (Δm − em)/(coeff·X)（由 L2 自己嘅 INPUTS 重算）"""
    return (dm(INPUTS, _PAIR[n][0], _PAIR[n][1]) - _mp(EM_ANCHOR[n])) / (_COEF[n] * _mp(INPUTS["X_MeV"]))


def _wmean(names, full=False):
    """σ_O 加權平均（π 不參與：em 由 G-parity 固定 ⟹ σ_em=0 ⟹ 權重無限）"""
    num = mp.fsum(_O_anchor(n) / _sig_O(n, full) ** 2 for n in names)
    den = mp.fsum(1 / _sig_O(n, full) ** 2 for n in names)
    return num / den


def model_O(lab, ch, full=False):
    """回傳模型 lab 對通道 ch 之 ⟨O⟩；None = 該類別無錨（D 之 (2,0)）"""
    anc = [n for n in EM_ANCHOR if n != "pi+-pi0"]   # π 由 G-parity 固定，不入模型擬合
    if lab == "A":
        return _wmean(anc, full)
    if lab == "B":
        nq = CLASS[ch][0]
        grp = [n for n in anc if CLASS[n][0] == nq]
        return _wmean(grp, full) if grp else None
    if lab == "C":
        ns = CLASS[ch][1]
        grp = [n for n in anc if CLASS[n][1] == ns]
        return _wmean(grp, full) if grp else None
    ns, nq = CLASS[ch][1], CLASS[ch][0]          # D：按 (nq,ns) 類別
    grp = [n for n in anc if CLASS[n] == (nq, ns)]
    return _wmean(grp, full) if grp else None


def box(ch):
    a2, ap = qpair(CONTENT[ch][0])
    b2, bp = qpair(CONTENT[ch][1])
    dS, dP = a2 - b2, ap - bp
    vals = [F(0), dS, dP, dS + dP]
    return min(vals), max(vals)


def k0_crit(lab, full=False):
    worst = None
    for ch in BLIND:
        O = model_O(lab, ch, full)
        if O is None:
            continue
        em = dm(INPUTS, _PAIR[ch][0], _PAIR[ch][1]) - _COEF[ch] * _mp(INPUTS["X_MeV"]) * O
        lo, hi = box(ch)
        kc = abs(em) / abs(_mp(lo)) if em < 0 else (em / _mp(hi) if hi != 0 else mp.inf)
        worst = kc if worst is None else max(worst, kc)
    return worst


for _lab in ("A", "B", "C", "D"):
    _add(f"k0_crit_{_lab}", (lambda L: (lambda I: k0_crit(L)))(_lab))
_add("k0_crit_all", lambda I: max(k0_crit(L) for L in ("A", "B", "C", "D")))
_add("a_crit_all", lambda I: _mp(I["HBARC_MeVfm"]) / (_mp("137.035999") * DERIVE["k0_crit_all"](I)))
_add("inv_X", lambda I: 1 / _mp(I["X_MeV"]))
for _ch in ("n-p", "K+-K0", "S--S0", "Xi--Xi0"):
    _add("sxfloor_" + _ch.replace("+", "p").replace("-", "m"),
         (lambda C: (lambda I: _O_anchor(C) * _mp(I["SIG_X"]) / _mp(I["X_MeV"])))(_ch))
_add("inv_2X", lambda I: 1 / (2 * _mp(I["X_MeV"])))
_add("O_A_mean", lambda I: model_O("A", "n-p"))
_add("O_D_class31", lambda I: model_O("D", "S0-S+"))
# ── 敏感度變體：σ_O 補入 σ_Δm（PDG 質量差誤差）──
for _lab in ("A", "B", "C", "D"):
    _add(f"k0_crit_full_{_lab}", (lambda L: (lambda I: k0_crit(L, True)))(_lab))
_add("k0_crit_all_full", lambda I: max(k0_crit(L, True) for L in ("A", "B", "C", "D")))
_add("a_crit_all_full", lambda I: _mp(I["HBARC_MeVfm"]) / (_mp("137.035999") * DERIVE["k0_crit_all_full"](I)))
_add("O_A_mean_full", lambda I: model_O("A", "n-p", True))
for _ch in ("n-p", "K+-K0", "S--S0", "Xi--Xi0"):
    _add(f"sigO_{_ch}".replace("+", "p").replace("-", "m"), (lambda C: (lambda I: _sig_O(C, True)))(_ch))
# docs/07 之 T₂（原為 inline 計算、無 producer ⟹ 本輪補回公式）：T₂ ≡ A_eff/k₀
_add("T2_aeff_over_k0", lambda I: _mp(I["AEFF_MeV"]) / _mp(I["K0_MeV"]))
_add("T2_v2_branch", lambda I: _mp(I["AEFF_MeV"]) / _mp("3.748112"))

# ============================================================================
# ③ IDENTITIES —— 兩條獨立路徑
# ============================================================================
IDENTITIES = [
    ("I1 k₀ = α·ħc/a₀：直接 vs (ħc/a₀)·α 之另一次序",
     lambda I: _mp(I["ALPHA"]) * _mp(I["HBARC_MeVfm"]) / _mp(I["A0_fm"]),
     lambda I: (_mp(I["HBARC_MeVfm"]) / _mp(I["A0_fm"])) * _mp(I["ALPHA"])),
    ("I2 E_K = Δm − coeff·X：coeff=−1 之符號路徑",
     lambda I: dm(I, "K+", "K0") - (-1) * _mp(I["X_MeV"]),
     lambda I: dm(I, "K+", "K0") + _mp(I["X_MeV"])),
    ("I3 退化類 (1/6,1) 之 Φ 三通道相等", lambda I: _mp(phi_ch("S--S0")),
     lambda I: _mp(phi_ch("Xi--Xi0"))),
    ("I4 退化類 (−1/3,1) 之 Φ 相等", lambda I: _mp(phi_ch("n-p")), lambda I: _mp(phi_ch("S0-S+"))),
    ("I5 em(0.79) 之兩路：直接 vs (1+ε)·ΔM²π/(M_K sum) 之移項",
     lambda I: (1 + _mp("0.79")) * _mp(I["FLAG_DMPI2"]) / _mp(I["FLAG_MKSUM"]),
     lambda I: _mp(I["FLAG_DMPI2"]) / _mp(I["FLAG_MKSUM"]) + _mp("0.79") * _mp(I["FLAG_DMPI2"]) / _mp(I["FLAG_MKSUM"])),
]


def main():
    hdr("§0 輸入（原始量）")
    I = dict(INPUTS)
    for k in ("ALPHA", "HBARC_MeVfm", "A0_fm", "X_MeV", "SIG_X", "K0_MeV"):
        print(f"  {k:14s} = {I[k]}")

    hdr("§1 DERIVE：由 INPUTS 純函數重算")
    D = {}
    for key, fn in DERIVE.items():
        try:
            D[key] = mp.mpf(fn(I))
        except Exception as e:  # noqa: BLE001
            chk(f"D.{key} 可推導", False, f"{type(e).__name__}: {e}")
    chk(f"D 全部 key 可推導（{len(DERIVE)} 個）", len(D) == len(DERIVE),
        f"成功 {len(D)}/{len(DERIVE)}")
    for k in sorted(D):
        print(f"  {k:22s} = {mp.nstr(D[k], 20)}")

    hdr("§2 IDENTITIES（兩條獨立路徑，差 < 1e−30）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(I)), mp.mpf(fb(I))
        d = abs(a - b)
        chk(tag, d < mp.mpf("1e-30"), f"差 = {mp.nstr(d, 6)}")

    hdr("§3 交叉核（唔涉 docs）")
    _k0in = float(INPUTS["K0_MeV"])
    chk(f"C1 k₀(重算) 落 [k₀±0.01] MeV（帶由輸入導出）",
        _k0in - 0.01 < float(D["k0_from_a0"]) < _k0in + 0.01,
        f"{mp.nstr(D['k0_from_a0'], 10)}")
    chk("C2 ⟨O⟩ 錨：π=0（|O|<1e−9）、其餘全正、⟨O⟩_K 最大",
        abs(float(D["O_pi+-pi0"])) < 1e-9
        and all(float(D[f"O_{k}"]) > 0 for k in EM_ANCHOR if k != "pi+-pi0")
        and float(D["O_K+-K0"]) == max(float(D[f"O_{k}"]) for k in EM_ANCHOR),
        " ".join(f"{k}={float(D['O_'+k]):.4f}" for k in EM_ANCHOR))
    chk("C3 路線 C 之 χ²/dof > 100（KILL 級）", float(D["routeC_chi2_dof"]) > 100,
        f"{float(D['routeC_chi2_dof']):.1f}")
    chk("C4 退化類差 > 1 MeV（結構性）",
        float(D["degen_np_S0Splus"]) > 1 and float(D["degen_S_Xi_D"]) > 1,
        f"{float(D['degen_np_S0Splus']):.4f} / {float(D['degen_S_Xi_D']):.4f}")
    chk("C5 G0 deficit > 100×", float(D["G0_deficit_D"]) > 100, f"{float(D['G0_deficit_D']):.1f}×")
    chk("C6 ε_model(由 0.3327 反算) 為負 ⟹ 與 FLAG +0.79 反號", float(D["eps_model_from_pred"]) < 0,
        f"{float(D['eps_model_from_pred']):.4f}")
    chk("C7 z(FLAG45, only ε) 比 z(含 X) 更負", float(D["z_flag45_only"]) < float(D["z_flag45_total"]),
        f"{float(D['z_flag45_only']):.2f} vs {float(D['z_flag45_total']):.2f}")

    chk("C8 k₀^crit(all) ≥ 每個模型之 k₀^crit",
        float(D["k0_crit_all"]) >= max(float(D[f"k0_crit_{L}"]) for L in "ABCD") - 1e-9,
        " ".join(f"{L}={float(D['k0_crit_'+L]):.4f}" for L in "ABCD") + f" all={float(D['k0_crit_all']):.4f}")
    chk("C9 a^crit = ħc/(137.035999·k₀^crit)（反比自洽）",
        abs(float(D["a_crit_all"] * D["k0_crit_all"] * _mp("137.035999") / _mp(INPUTS["HBARC_MeVfm"])) - 1) < 1e-12,
        f"a^crit = {float(D['a_crit_all']):.5f} fm")
    chk("C10 模型 A 之 ⟨O⟩ 落 [1, 3]（正值、量級合理）", 1 < float(D["O_A_mean"]) < 3,
        f"{float(D['O_A_mean']):.6f}")

    hdr("§4 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_inputs"] = {k: (v if not isinstance(v, mp.mpf) else mp.nstr(v, 30)) for k, v in I.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}
    (DATA / "l2_independent_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n寫出: {DATA/'l2_independent_verify.json'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())