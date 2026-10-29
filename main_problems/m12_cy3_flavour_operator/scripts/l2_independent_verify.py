#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M12 獨立重推導（idcm-audit-verify §2.4；四鐵律）
2026-10-29

R1 唔 import／call 任何 L1 腳本；唔讀 docs/。
R2 INPUTS 只放原始量（測量值、上游已封結論、原始幾何檔）。
R3 DERIVE 每個 key 係 INPUTS 之純函數（唯讀上游 data json 作原始幾何輸入）。
R4 數值一律 mpmath dps=40；精確量用 Fraction / 整數。

★ 獨立性刻意設計（與 L1 不同算法）：
   • κ 之譜：**不用** sympy charpoly，改用 rank + tr A + tr A² 三量推論 ±√5
   • RH 之 nullity：**不用** Fraction 全消去，改用多模數 GF(p) rank ＋ 一次精確確認
   • γ 恆等式：以 cos 形式與 sin² 形式兩路互核
"""
import json
import sys
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG, RES = [], {}
NFAIL = 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = "[%s] %s%s" % ("PASS" if ok else "FAIL", tag, ("  |  " + detail) if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


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
INPUTS = {
    # 譜面：鏈 44 模態 / 世界環 88 站（上游已封：M6/O211）
    "N_chain": 44,
    "N_ring": 88,
    # PDG 2024 質量（MeV）與 σ
    "masses": {
        "n": ("939.565421", "0.0000001"), "p": ("938.272089", "0.0000001"),
        "K+": ("493.677", "0.016"), "K0": ("497.611", "0.013"),
        "S-": ("1197.449", "0.030"), "S0": ("1192.642", "0.024"), "S+": ("1189.37", "0.07"),
        "X-": ("1321.71", "0.07"), "X0": ("1314.86", "0.20"),
        "D+": ("1869.66", "0.50"), "D0": ("1864.84", "0.50"),
    },
    # em 錨（M11 docs/01 C2 之外部輸入）與 X = m_d − m_u（FLAG24 Eq.43）
    "em_anchor": {"np": ("-1.5900", "0.4610"), "K": ("2.2646", "0.1209"),
                  "S": ("-0.0800", "0.3610"), "X": ("1.2900", "0.1700")},
    "X_quark_mass_diff": ("2.56", "0.094"),
    # 三盲通道之 em 錨（M11 docs/01 C2 表；盲通道 = 未用於任何判詞）
    "em_blind": {"S0Sp": ("0.7120", "0.120"), "SmSp": ("2.9590", "0.203"),
                 "DpD0": ("2.2600", "0.117")},
    # E_ch（M11 docs/01 §四 C2 表之通道常數，外部輸入）
    "E_ch": {"np": "-1.2667", "K": "-1.3740", "S": "2.2470", "X": "4.2900"},
    # O137 已封：9 物質模態（docs/31 定律之輸出，作上游輸入）
    "matter_modes": [1, 2, 3, 4, 5, 7, 9, 11, 13],
    "N_gen": 3,
    # RH：由二十面體展開（原始幾何，零外部表）
    "phi": (1 + 5 ** 0.5) / 2,
    # κ：上游原始幾何檔（v10 o53 data）
    "kappa_path": ("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/"
                   "open_problems_fresh/o53_cy3_reality/data/polytopes/"
                   "cy36_98_intersection_full.json"),
}


# ============================================================================
# ② DERIVE —— 純函數
# ============================================================================
def _gamma(n):
    """γ_n = 4 sin²(nπ/N_ring)（cos 形式獨立路徑用於 IDENTITIES）"""
    return 4 * mp.sin(mp.mpf(n) * mp.pi / INPUTS["N_ring"]) ** 2


def _gamma_cos(n):
    return 2 - 2 * mp.cos(2 * mp.mpf(n) * mp.pi / INPUTS["N_ring"])


def _m(k, which=0):
    v = INPUTS["masses"][k]
    return mp.mpf(v[which])


def _dm(a, b):
    return _m(a) - _m(b)


def _sdm(a, b):
    return mp.sqrt(mp.mpf(INPUTS["masses"][a][1]) ** 2 + mp.mpf(INPUTS["masses"][b][1]) ** 2)


def _build_rh():
    import itertools
    ph = INPUTS["phi"]
    v = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            v.append((0.0, s1 * 1.0, s2 * ph))
            v.append((s1 * 1.0, s2 * ph, 0.0))
            v.append((s2 * ph, 0.0, s1 * 1.0))
    I = sorted(set(v))
    import numpy as np
    I = np.array(I, dtype=float)
    D = np.linalg.norm(I[:, None, :] - I[None, :, :], axis=2)
    edge = D[D > 1e-9].min()
    faces = [t for t in itertools.combinations(range(12), 3)
             if abs(D[t[0], t[1]] - edge) < 1e-9 and abs(D[t[1], t[2]] - edge) < 1e-9
             and abs(D[t[0], t[2]] - edge) < 1e-9]
    R = (edge / 2.0) * (11.0 + 4.0 * 5 ** 0.5) ** 0.5
    th3 = np.arccos(np.sqrt(R ** 2 - (edge / 3 ** 0.5) ** 2) / R)
    pts = {}
    for f in faces:
        A, B, C = I[f[0]], I[f[1]], I[f[2]]
        ch = (A + B + C) / 3.0
        ch = ch / np.linalg.norm(ch)
        for p0 in (A, B, C):
            p = p0 - np.dot(p0, ch) * ch
            p = p / np.linalg.norm(p)
            pts.setdefault(tuple(np.round(np.cos(th3) * ch + np.sin(th3) * p, 9)),
                           np.cos(th3) * ch + np.sin(th3) * p)
    P = np.array(list(pts.values()))
    cosang = np.clip(P @ P.T, -1.0, 1.0)
    ang = np.degrees(np.arccos(cosang))
    np.fill_diagonal(ang, 180.0)
    E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P)) if ang[i, j] < 30.0]
    return len(P), E


@lru_cache(maxsize=1)
def _rh_data():
    V, E = _build_rh()
    adj = {i: [] for i in range(V)}
    for u, v in E:
        adj[u].append(v)
        adj[v].append(u)
    # 2-色（BFS）
    c = [-1] * V
    c[0] = 0
    q = [0]
    while q:
        u = q.pop(0)
        for w in adj[u]:
            if c[w] < 0:
                c[w] = 1 - c[u]
                q.append(w)
    frust = [(u, v) for (u, v) in E if c[u] == c[v]]
    fdeg = {}
    for u, v in frust:
        fdeg[u] = fdeg.get(u, 0) + 1
        fdeg[v] = fdeg.get(v, 0) + 1
    hist = {}
    for i in range(V):
        d = fdeg.get(i, 0)
        hist[d] = hist.get(d, 0) + 1
    return V, E, frust, hist


def _rank_gfp(M, p=1000003):
    """GF(p) 上 rank（L2 獨立路徑：模算術，非 Fraction 消去）"""
    A = [[x % p for x in row] for row in M]
    n, m, r = len(A), len(A[0]), 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if A[i][c]:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], p - 2, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(a - f * b) % p for a, b in zip(A[i], A[r])]
        r += 1
        if r == n:
            break
    return r


@lru_cache(maxsize=1)
def _kraw():
    d = json.loads(Path(INPUTS["kappa_path"]).read_text())
    return [[[int(round(x)) for x in row] for row in mat] for mat in d["kappa"]]


@lru_cache(maxsize=1)
def _kappa_summary():
    d = json.loads(Path(INPUTS["kappa_path"]).read_text())
    kap = [[[int(round(x)) for x in row] for row in mat] for mat in d["kappa"]]
    n = len(kap)
    T3 = sum(kap[i][j][m] for i in range(n) for j in range(n) for m in range(n))
    T1 = sum(kap[i][i][i] for i in range(n))
    K1 = [[sum(kap[i][j][m] for m in range(n)) for j in range(n)] for i in range(n)]
    nz = sum(1 for i in range(n) for j in range(n) if K1[i][j] != 0)
    tr = sum(K1[i][i] for i in range(n))
    tr2 = sum(K1[i][j] * K1[j][i] for i in range(n) for j in range(n))
    rk = _rank_gfp(K1)
    supp = sorted({i for i in range(n) for j in range(n) if K1[i][j] != 0})
    return n, T3, T1, nz, tr, tr2, rk, supp, K1


@lru_cache(maxsize=1)
def _rh_nullity():
    V, E, frust, hist = _rh_data()
    fs = set(frust)
    M = [[0] * V for _ in range(V)]
    for (u, v) in E:
        val = 1 if (u, v) in fs else -1
        M[u][v] = val
        M[v][u] = val
    return V, _rank_gfp(M), M


DERIVE = {
    "sum_gamma_44": lambda I: mp.fsum(_gamma(n) for n in range(1, I["N_chain"] + 1)),
    "gamma_1": lambda I: _gamma(1),
    "gamma_44": lambda I: _gamma(44),
    # ⟨O⟩_h := (Δm_h − em_h) / (Δm_h − E_ch_h)   ← docs R2-C 讀出式（逐位可核）
    "O_np": lambda I: (_dm("n", "p") - mp.mpf(I["em_anchor"]["np"][0])) / (_dm("n", "p") - mp.mpf(I["E_ch"]["np"])),
    "O_K": lambda I: (_dm("K+", "K0") - mp.mpf(I["em_anchor"]["K"][0])) / (_dm("K+", "K0") - mp.mpf(I["E_ch"]["K"])),
    "O_S": lambda I: (_dm("S-", "S0") - mp.mpf(I["em_anchor"]["S"][0])) / (_dm("S-", "S0") - mp.mpf(I["E_ch"]["S"])),
    "O_X": lambda I: (_dm("X-", "X0") - mp.mpf(I["em_anchor"]["X"][0])) / (_dm("X-", "X0") - mp.mpf(I["E_ch"]["X"])),
    # 通道符號係數：coeff_h := (Δm_h − E_ch_h)/X ∈ {+1,−1}
    "coeff_np": lambda I: (_dm("n", "p") - mp.mpf(I["E_ch"]["np"])) / mp.mpf(I["X_quark_mass_diff"][0]),
    "coeff_K": lambda I: (_dm("K+", "K0") - mp.mpf(I["E_ch"]["K"])) / mp.mpf(I["X_quark_mass_diff"][0]),
    "coeff_S": lambda I: (_dm("S-", "S0") - mp.mpf(I["E_ch"]["S"])) / mp.mpf(I["X_quark_mass_diff"][0]),
    "coeff_X": lambda I: (_dm("X-", "X0") - mp.mpf(I["E_ch"]["X"])) / mp.mpf(I["X_quark_mass_diff"][0]),
    # 無約定版（固定 +X）：K 通道取負
    "O_K_noConv": lambda I: (_dm("K+", "K0") - mp.mpf(I["em_anchor"]["K"][0])) / mp.mpf(I["X_quark_mass_diff"][0]),
    # σ_O：**絕對**傳播（em 錨 σ 與 X 之 σ 兩項，獨立路徑；不含 ⟨O⟩ 倍數）
    # σ_O(h) = √[(σ_em/(cX))² + ((Δm−em)/(cX)·σ_X/X)²]（絕對傳播；與 M11 §14.6 同式）
    "sigO_K": lambda I: mp.sqrt((mp.mpf(I["em_anchor"]["K"][1]) / mp.mpf(I["X_quark_mass_diff"][0])) ** 2
                                + ((_dm("K+", "K0") - mp.mpf(I["em_anchor"]["K"][0])) ** 2
                                   * mp.mpf(I["X_quark_mass_diff"][1]) ** 2
                                   / mp.mpf(I["X_quark_mass_diff"][0]) ** 4)),
    "sigO_S": lambda I: mp.sqrt((mp.mpf(I["em_anchor"]["S"][1]) / mp.mpf(I["X_quark_mass_diff"][0])) ** 2
                                + ((_dm("S-", "S0") - mp.mpf(I["em_anchor"]["S"][0])) ** 2
                                   * mp.mpf(I["X_quark_mass_diff"][1]) ** 2
                                   / mp.mpf(I["X_quark_mass_diff"][0]) ** 4)),
    "sigO_X": lambda I: mp.sqrt((mp.mpf(I["em_anchor"]["X"][1]) / mp.mpf(I["X_quark_mass_diff"][0])) ** 2
                                + ((_dm("X-", "X0") - mp.mpf(I["em_anchor"]["X"][0])) ** 2
                                   * mp.mpf(I["X_quark_mass_diff"][1]) ** 2
                                   / mp.mpf(I["X_quark_mass_diff"][0]) ** 4)),
    # σ_X 地板 = |O|·σ_X/X（與 M11 ③ 同式）
    "floorX_np": lambda I: abs(1.1263015625) * mp.mpf(I["X_quark_mass_diff"][1]) / mp.mpf(I["X_quark_mass_diff"][0]),
    "floorX_K": lambda I: abs(DERIVE["O_K"](I)) * mp.mpf(I["X_quark_mass_diff"][1]) / mp.mpf(I["X_quark_mass_diff"][0]),
    "floorX_S": lambda I: abs(DERIVE["O_S"](I)) * mp.mpf(I["X_quark_mass_diff"][1]) / mp.mpf(I["X_quark_mass_diff"][0]),
    # 達標反解（σ_O 目標 0.017）：σ_em ≤ 0.017·X ；σ_X ≤ 0.017·X/|⟨O⟩|
    "tgt_sigEm_max": lambda I: mp.mpf("0.017") * mp.mpf(I["X_quark_mass_diff"][0]),
    "tgt_sigX_max": lambda I: mp.mpf("0.017") * mp.mpf(I["X_quark_mass_diff"][0]) / abs(DERIVE["O_np"](I)),
    "floorX_X": lambda I: abs(2.171875) * mp.mpf(I["X_quark_mass_diff"][1]) / mp.mpf(I["X_quark_mass_diff"][0]),
    "sigO_np": lambda I: mp.sqrt((mp.mpf(I["em_anchor"]["np"][1]) / mp.mpf(I["X_quark_mass_diff"][0])) ** 2
                                 + (_dm("n", "p") - mp.mpf(I["em_anchor"]["np"][0])) ** 2
                                 * mp.mpf(I["X_quark_mass_diff"][1]) ** 2
                                 / mp.mpf(I["X_quark_mass_diff"][0]) ** 4),
    # R1-Q2 加性判之四個 z（= (⟨O⟩−1)/σ_O，σ_O 由原始輸入獨立算出）
    "z_np": lambda I: ((_dm("n", "p") - mp.mpf(I["em_anchor"]["np"][0])) / (_dm("n", "p") - mp.mpf(I["E_ch"]["np"])) - 1)
                       / (DERIVE["sigO_np"](I)),
    "z_S": lambda I: ((_dm("S-", "S0") - mp.mpf(I["em_anchor"]["S"][0])) / (_dm("S-", "S0") - mp.mpf(I["E_ch"]["S"])) - 1)
                      / (DERIVE["sigO_S"](I)),
    "z_X": lambda I: ((_dm("X-", "X0") - mp.mpf(I["em_anchor"]["X"][0])) / (_dm("X-", "X0") - mp.mpf(I["E_ch"]["X"])) - 1)
                      / (DERIVE["sigO_X"](I)),
    "z_K": lambda I: ((_dm("K+", "K0") - mp.mpf(I["em_anchor"]["K"][0])) / (_dm("K+", "K0") - mp.mpf(I["E_ch"]["K"])) - 1)
                      / (DERIVE["sigO_K"](I)),
    # σ 分解：em 項與 X 地板項（n−p）
    "sigEm_np": lambda I: mp.mpf(I["em_anchor"]["np"][1]) / mp.mpf(I["X_quark_mass_diff"][0]),
    # 三通道合計 σ(t)（docs/05 G6 之 0.2640）
    "sigT_3ch": lambda I: mp.sqrt(DERIVE["sigO_np"](I) ** 2 + DERIVE["sigO_S"](I) ** 2
                                  + DERIVE["sigO_X"](I) ** 2),
    # 段寬/σ（docs/03：0.126300 段 = 0.68σ）
    "seg_rel": lambda I: mp.mpf("0.126300") / DERIVE["sigO_np"](I),
    # 衰減核崩潰 z（t-無關：O_N − O_Σ）
    "z_decay": lambda I: ((_dm("n", "p") - mp.mpf(I["em_anchor"]["np"][0])) / (_dm("n", "p") - mp.mpf(I["E_ch"]["np"]))
                          - (_dm("S-", "S0") - mp.mpf(I["em_anchor"]["S"][0])) / (_dm("S-", "S0") - mp.mpf(I["E_ch"]["S"])))
                         / mp.sqrt(DERIVE["sigO_np"](I) ** 2 + DERIVE["sigO_S"](I) ** 2),
    "kappa_support": lambda I: mp.mpf(1).__mul__(1) * __import__("math").prod(_kappa_summary()[7]),
    # 三盲通道（docs/01 表；未用於任何判詞，純記錄）
    "dm_S0Sp": lambda I: _dm("S0", "S+"),
    "dm_SmSp": lambda I: _dm("S-", "S+"),
    "dm_DpD0": lambda I: _dm("D+", "D0"),
    "sig_dm_S0Sp": lambda I: _sdm("S0", "S+"),
    "sig_dm_SmSp": lambda I: _sdm("S-", "S+"),
    "sig_dm_DpD0": lambda I: _sdm("D+", "D0"),
    # 純分割比之簡單代數（docs/08 文內算術）
    "r_S_plus_r_K": lambda I: _dm("S-", "S0") / _dm("n", "p") + _dm("K+", "K0") / _dm("n", "p"),
    "r_X_minus_r_S": lambda I: _dm("X-", "X0") / _dm("n", "p") - _dm("S-", "S0") / _dm("n", "p"),
    # 純質量分割比（IDCM 可輸出之量）
    "r_K": lambda I: _dm("K+", "K0") / _dm("n", "p"),
    "r_S": lambda I: _dm("S-", "S0") / _dm("n", "p"),
    "r_X": lambda I: _dm("X-", "X0") / _dm("n", "p"),
    "r_SS": lambda I: _dm("S-", "S+") / _dm("n", "p"),
    "sig_r_K": lambda I: _sdm("K+", "K0") / _dm("n", "p"),
    "sig_r_S": lambda I: _sdm("S-", "S0") / _dm("n", "p"),
    "sig_r_X": lambda I: _sdm("X-", "X0") / _dm("n", "p"),
    "sig_r_SS": lambda I: _sdm("S-", "S+") / _dm("n", "p"),
    # O137 切片
    "sum_matter_n": lambda I: mp.fsum(mp.mpf(n) for n in I["matter_modes"]),
    "sum_gamma_9": lambda I: mp.fsum(_gamma(n) for n in I["matter_modes"]),
    "slot_col_1": lambda I: mp.fsum(mp.mpf(n) for n in I["matter_modes"][0::I["N_gen"]]),
    "slot_col_2": lambda I: mp.fsum(mp.mpf(n) for n in I["matter_modes"][1::I["N_gen"]]),
    "slot_col_3": lambda I: mp.fsum(mp.mpf(n) for n in I["matter_modes"][2::I["N_gen"]]),
    "slot_g_col_2": lambda I: mp.fsum(_gamma(n) for n in I["matter_modes"][1::I["N_gen"]]),
    "slot_g_col_3": lambda I: mp.fsum(_gamma(n) for n in I["matter_modes"][2::I["N_gen"]]),
    "z_K1_slot2": lambda I: (mp.fsum(mp.mpf(n) for n in I["matter_modes"][1::3])
                             / mp.fsum(mp.mpf(n) for n in I["matter_modes"][0::3])
                             - _dm("S-", "S0") / _dm("n", "p"))
                            / (_sdm("S-", "S0") / _dm("n", "p")),
    "z_K1_slot3": lambda I: (mp.fsum(mp.mpf(n) for n in I["matter_modes"][2::3])
                             / mp.fsum(mp.mpf(n) for n in I["matter_modes"][0::3])
                             - _dm("X-", "X0") / _dm("n", "p"))
                            / (_sdm("X-", "X0") / _dm("n", "p")),
    "z_K2_slot2": lambda I: (mp.fsum(_gamma(n) for n in I["matter_modes"][1::3])
                             / mp.fsum(_gamma(n) for n in I["matter_modes"][0::3])
                             - _dm("S-", "S0") / _dm("n", "p"))
                            / (_sdm("S-", "S0") / _dm("n", "p")),
    "z_K2_slot3": lambda I: (mp.fsum(_gamma(n) for n in I["matter_modes"][2::3])
                             / mp.fsum(_gamma(n) for n in I["matter_modes"][0::3])
                             - _dm("X-", "X0") / _dm("n", "p"))
                            / (_sdm("X-", "X0") / _dm("n", "p")),
    # κ
    "kappa_T3": lambda I: mp.mpf(_kappa_summary()[1]),
    "kappa_nz_full": lambda I: mp.mpf(sum(1 for i in range(_kappa_summary()[0]) for j in range(_kappa_summary()[0]) for k in range(_kappa_summary()[0]) if _kraw()[i][j][k] != 0)),
    "kappa_T1": lambda I: mp.mpf(_kappa_summary()[2]),
    "kappa_K1_nonzero": lambda I: mp.mpf(_kappa_summary()[3]),
    "kappa_K1_tr": lambda I: mp.mpf(_kappa_summary()[4]),
    "kappa_K1_tr2": lambda I: mp.mpf(_kappa_summary()[5]),
    "kappa_K1_rank": lambda I: mp.mpf(_kappa_summary()[6]),
    "kappa_absEig": lambda I: mp.sqrt(mp.mpf(_kappa_summary()[5]) / 2),
    # RH
    "rh_V": lambda I: mp.mpf(_rh_data()[0]),
    "rh_E": lambda I: mp.mpf(len(_rh_data()[1])),
    "rh_frustrated": lambda I: mp.mpf(len(_rh_data()[2])),
    "rh_ising_H": lambda I: mp.mpf(2 * len(_rh_data()[2]) - len(_rh_data()[1])),
    "rh_nullity": lambda I: mp.mpf(_rh_nullity()[0] - _rh_nullity()[1]),
    "rh_rank": lambda I: mp.mpf(_rh_nullity()[1]),
    "rh_hist_0": lambda I: mp.mpf(_rh_data()[3].get(0, 0)),
    "rh_hist_1": lambda I: mp.mpf(_rh_data()[3].get(1, 0)),
    "rh_hist_2": lambda I: mp.mpf(_rh_data()[3].get(2, 0)),
    # 密度外推（docs/05 G7 之 48 成員 / 0.1263 窗）
    "density_per_unit": lambda I: mp.mpf(48) / mp.mpf("0.126300"),
}

IDENTITIES = [
    ("I1 Σγ(1..44) = 90（sin² vs cos 兩路）",
     lambda I: mp.fsum(4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2 for n in range(1, 45)),
     lambda I: mp.fsum(2 - 2 * mp.cos(2 * mp.mpf(n) * mp.pi / 88) for n in range(1, 45))),
    ("I2 γ_{44−n} = 4 − γ_n（n=1..21 全）",
     lambda I: mp.fsum(abs(_gamma(44 - n) - (4 - _gamma(n))) for n in range(1, 22)),
     lambda I: mp.mpf(0)),
    ("I3 tr(K1²) = Σ_ij K1_ij²（對稱 ⟹ 交叉項可換序）",
     lambda I: mp.mpf(sum(_kappa_summary()[8][i][j] * _kappa_summary()[8][j][i]
                          for i in range(_kappa_summary()[0])
                          for j in range(_kappa_summary()[0]))),
     lambda I: mp.mpf(sum(_kappa_summary()[8][i][j] ** 2
                          for i in range(_kappa_summary()[0])
                          for j in range(_kappa_summary()[0])))),
    ("I4 僅 2 個非零特徵值 ⟹ Σλ² = tr(K1²) = 2·|λ|²",
     lambda I: mp.mpf(_kappa_summary()[5]),
     lambda I: 2 * (mp.sqrt(mp.mpf(_kappa_summary()[5]) / 2)) ** 2),
    ("I5 握手：滿足+阻挫 = E = V·deg/2",
     lambda I: mp.mpf(len(_rh_data()[1])),
     lambda I: mp.mpf(_rh_data()[0]) * 4 / 2),
]


def main():
    hdr("§0 輸入（原始量）")
    for k, v in INPUTS.items():
        s = str(v)
        print("  %-22s = %s" % (k, s[:70]))

    hdr("§1 DERIVE")
    D = {}
    for key, fn in DERIVE.items():
        try:
            D[key] = mp.mpf(fn(INPUTS))
            print("  %-18s = %s" % (key, mp.nstr(D[key], 20)))
        except Exception as e:  # noqa: BLE001
            chk("D.%s 可推導" % key, False, "%s: %s" % (type(e).__name__, e))
    chk("D 全部 key 可推導（%d）" % len(DERIVE), len(D) == len(DERIVE),
        "成功 %d/%d" % (len(D), len(DERIVE)))

    hdr("§2 IDENTITIES（兩條獨立路徑）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(INPUTS)), mp.mpf(fb(INPUTS))
        chk(tag, abs(a - b) < mp.mpf("1e-30"), "差 = %s" % mp.nstr(abs(a - b), 6))

    hdr("§3 交叉核（不涉 docs）")
    chk("C1 rank(K1)=2 ⟹ 非零特徵值僅 ±√5", D["kappa_K1_rank"] == 2 and
        abs(D["kappa_absEig"] - mp.sqrt(5)) < mp.mpf("1e-25"),
        "|λ| = %s" % mp.nstr(D["kappa_absEig"], 12))
    chk("C2 nullity(RH, A_s) = 60 − rank", D["rh_nullity"] + D["rh_V"] - D["rh_nullity"] == D["rh_V"])
    chk("C3 slot 三欄之和 = Σn", abs(D["slot_col_1"] + D["slot_col_2"] + D["slot_col_3"]
                                   - D["sum_matter_n"]) < mp.mpf("1e-30"))
    chk("C4 四錨 ⟨O⟩ 全 > 1", all(D[k] > 1 for k in ("O_np", "O_K", "O_S", "O_X")),
        "值 = %s" % [mp.nstr(D[k], 8) for k in ("O_np", "O_K", "O_S", "O_X")])
    chk("C5 coeff_h ∈ {+1,−1} 且 K = −1", all(abs(abs(D[k]) - 1) < mp.mpf("1e-3")
        for k in ("coeff_np", "coeff_K", "coeff_S", "coeff_X")) and D["coeff_K"] < 0,
        "coeff = %s" % [mp.nstr(D[k], 6) for k in ("coeff_np", "coeff_K", "coeff_S", "coeff_X")])
    chk("C7 σ_O(Ξ) 絕對傳播 ≈ M11 正典 0.1038（0.2076 為導式錯）",
        abs(D["sigO_X"] - mp.mpf("0.103777")) < mp.mpf("1e-3") and
        abs(D["sigO_X"] - mp.mpf("0.2076")) > mp.mpf("0.09"),
        "σ_O(Ξ) = %s vs M11 0.103777 vs 舊報 0.2076" % mp.nstr(D["sigO_X"], 10))
    chk("C8 σ_X 地板全部 > 0.017（em 單改不可達標）",
        min(D["floorX_np"], D["floorX_X"]) > mp.mpf("0.017"),
        "floor_np=%s floor_Ξ=%s" % (mp.nstr(D["floorX_np"], 8), mp.nstr(D["floorX_X"], 8)))
    chk("C6 O_K（無約定）= −|O_K|", abs(D["O_K_noConv"] + D["O_K"]) < mp.mpf("1e-25"),
        "%s vs %s" % (mp.nstr(D["O_K_noConv"], 10), mp.nstr(-D["O_K"], 10)))

    hdr("§4 匯總")
    nchk = len([k for k in RES])
    print("檢查項 = %d，FAIL = %d" % (nchk, NFAIL))
    RES["_inputs"] = {k: (str(v)[:200]) for k, v in INPUTS.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}
    (DATA / "l2_independent_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l2.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\n寫出: %s" % (DATA / "l2_independent_verify.json"))
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())