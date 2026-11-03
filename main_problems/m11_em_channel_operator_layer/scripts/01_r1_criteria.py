#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M11 scripts/01 — R1：判準化 ＋ 輸入重建 ＋ 目標釘死 ＋ 可辨識性（不生模型）

判準（事前寫死於 docs/01，本腳本只執行）：
  C1  輸入重建：PDG2024 質量獨立重算 8 條 Δm、E_ch、σ_E；與 M4 docs/18 比 |Δ|≤1e-3 MeV
  C2  目標釘死：Ô 讀 ⟨h|d̄d−ūu|h⟩ ≡ ⟨O⟩_ch（強 IB 矩陣元）；em 部分由 5 個外部錨固定
  C4  可辨識性：13×16 約束矩陣之精確秩／核維（ℚ 上 Gauss 消去）
  C9  H1 可行性（精確箱測試）：em_ch = k₀[ΔΣq²·S + ΔΣ_{i<j}q_iq_j·P]，S,P ∈ (0,1]
      —— k₀ = α ħc/a₀；符號類判詞與 k₀ 無關，幅度類判詞明列 k_req
  C11 尺度：k_req/k₀（每通道所需），a₀ 敏感度已標；k₀ 屬 M5 輸入（非本記錄零參數）
  C12 紀律：Fraction 精確優先；JSON 與 stdout 內容一致；/usr/bin/python3

禁：不得引入任何擬合係數；不得以 E_ch 作 Ô 之目標（C9 先判其可行性）。
"""

import json
import math
import os
from fractions import Fraction as F

OUT_JSON = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "data", "01_r1_criteria.json")
RESULT = {}

X = F("2.56")           # m_d − m_u（FLAG24 中央值，外部固定）
SX = 0.094              # σ_X（FLAG24）

# ── 輸入 1：PDG 2024 質量（MeV）。docs/18 宣告其來源為 PDG 2024；此處為 reconstruction 用。
MASS = {
    "p":  (F("938.27208943"), 0.00000029),
    "n":  (F("939.56542052"), 0.00000054),
    "K+": (F("493.677"), 0.016),
    "K0": (F("497.611"), 0.013),
    "pi+": (F("139.57039"), 0.00018),
    "pi0": (F("134.9768"), 0.0005),
    "S-": (F("1197.449"), 0.030),
    "S0": (F("1192.642"), 0.024),
    "S+": (F("1189.37"), 0.07),
    "Xi-": (F("1321.71"), 0.07),
    "Xi0": (F("1314.86"), 0.20),
    "D+": (F("1869.66"), 0.05),
    "D0": (F("1864.84"), 0.05),
}

# ── 輸入 2：M4 docs/18 八通道表（上游，非本記錄結論）
DOCS18 = {
    "n-p":     dict(dm="1.2933",  coeff=1,  E="-1.2667", sE=0.0940),
    "K+-K0":   dict(dm="-3.9340", coeff=-1, E="-1.3740", sE=0.0963),
    "pi+-pi0": dict(dm="4.5936",  coeff=1,  E="2.0336",  sE=0.0940),
    "S--S0":   dict(dm="4.8070",  coeff=1,  E="2.2470",  sE=0.1010),
    "S0-S+":   dict(dm="3.2720",  coeff=1,  E="0.7120",  sE=0.1196),
    "S--S+":   dict(dm="8.0790",  coeff=2,  E="2.9590",  sE=0.2028),
    "Xi--Xi0": dict(dm="6.8500",  coeff=1,  E="4.2900",  sE=0.2301),
    "D+-D0":   dict(dm="4.8200",  coeff=1,  E="2.2600",  sE=0.1172),
}

# ── 通道定義：Δm = m(final) − m(initial)，取向照 docs/18 命名
CHANNELS = [
    ("n-p",     [("d", -1), ("d", -1), ("u", +2)],            [("u", +2), ("u", +2), ("d", -1)]),
    ("K+-K0",   [("u", +2), ("sbar", +1)],                     [("d", -1), ("sbar", +1)]),
    ("pi+-pi0", [("u", +2), ("dbar", +1)],                     "PI0"),
    ("S--S0",   [("d", -1), ("d", -1), ("s", -1)],             [("u", +2), ("d", -1), ("s", -1)]),
    ("S0-S+",   [("u", +2), ("d", -1), ("s", -1)],             [("u", +2), ("u", +2), ("s", -1)]),
    ("S--S+",   [("d", -1), ("d", -1), ("s", -1)],             [("u", +2), ("u", +2), ("s", -1)]),
    ("Xi--Xi0", [("d", -1), ("s", -1), ("s", -1)],             [("u", +2), ("s", -1), ("s", -1)]),
    ("D+-D0",   [("c", +2), ("dbar", +1)],                     [("c", +2), ("ubar", -2)]),
]


def qsum2(quarks):
    """Σ_i q_i² 與 Σ_{i<j} q_i q_j（單位 e²，Fraction）。π0 由 |uū>,|dd̄> 疊加處理。"""
    if quarks == "PI0":
        # π0 = (|u ū> − |d d̄>)/√2 ⟹ 期望值 = 平均（交叉項為 0）
        s2_1, p_1 = qsum2([("u", +2), ("ubar", -2)])
        s2_2, p_2 = qsum2([("d", -1), ("dbar", +1)])
        return (s2_1 + s2_2) / 2, (p_1 + p_2) / 2
    qs = [F(q, 3) for _, q in quarks]
    s2 = sum(q * q for q in qs)
    pr = sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs)))
    return s2, pr


def main():
    log = []
    def emit(s=""):
        log.append(s)
        print(s)

    emit("M11 R1 — 判準化／輸入重建／目標釘死／可辨識性（不生模型）")
    emit("=" * 96)
    emit(f"X = m_d − m_u = {float(X):.2f} ± {SX:.3f} MeV（FLAG24，外部固定）")
    emit("")

    # ══════════ C1：輸入重建 ══════════
    emit("C1｜輸入重建（PDG2024 質量 → Δm, E_ch, σ_E）；判準：與 docs/18 |Δ| ≤ 1e-3 MeV")
    emit(f"{'channel':10s}{'Δm(recalc)':>13s}{'Δm(docs18)':>12s}{'|Δ|':>9s}"
         f"{'E_ch':>9s}{'E(docs18)':>10s}{'|ΔE|':>9s}{'σ_E':>9s}{'σ(docs18)':>10s}")
    pairs = {
        "n-p": ("n", "p"), "K+-K0": ("K+", "K0"), "pi+-pi0": ("pi+", "pi0"),
        "S--S0": ("S-", "S0"), "S0-S+": ("S0", "S+"), "S--S+": ("S-", "S+"),
        "Xi--Xi0": ("Xi-", "Xi0"), "D+-D0": ("D+", "D0"),
    }
    c1_rows = {}
    worst_dm = F(0)
    worst_E = F(0)
    for name in pairs:
        a, b = pairs[name]
        dm = MASS[a][0] - MASS[b][0]
        d_docs = F(DOCS18[name]["dm"])
        coeff = DOCS18[name]["coeff"]
        E = dm - coeff * X
        E_docs = F(DOCS18[name]["E"])
        sE = math.sqrt((coeff * SX) ** 2 + MASS[a][1] ** 2 + MASS[b][1] ** 2)
        d1, d2 = abs(dm - d_docs), abs(E - E_docs)
        worst_dm = max(worst_dm, d1)
        worst_E = max(worst_E, d2)
        emit(f"{name:10s}{float(dm):13.4f}{float(d_docs):12.4f}{float(d1):9.5f}"
             f"{float(E):9.4f}{float(E_docs):10.4f}{float(d2):9.5f}"
             f"{sE:9.4f}{DOCS18[name]['sE']:10.4f}")
        c1_rows[name] = dict(dm=float(dm), E=float(E), coeff=coeff, sE=sE,
                             dm_docs=float(d_docs), E_docs=float(E_docs),
                             abs_d_dm=float(d1), abs_d_E=float(d2))
    emit("")
    emit(f"  max|ΔΔm| = {float(worst_dm):.2e} MeV ｜ max|ΔE| = {float(worst_E):.2e} MeV ｜ 判準 1e-3")
    c1_verdict = "PASS" if (worst_dm <= F("0.001") and worst_E <= F("0.001")) else "FAIL"
    emit(f"  C1 = {c1_verdict}")
    RESULT["C1"] = dict(verdict=c1_verdict, max_abs_d_dm=float(worst_dm),
                        max_abs_d_E=float(worst_E), rows=c1_rows)

    # ══════════ C2：目標釘死 ⟨O⟩_ch ══════════
    emit("")
    emit("C2｜目標釘死：⟨O⟩_ch = (Δm_ch − em_ch)/(coeff·X)")
    emit("外部錨（全部為 em 側；取向 = docs/18 通道取向 final−initial）：")
    emit("  π : Dashen/G-parity ⟹ 強 IB 分量 → 0 ⟹ ⟨O⟩_π = 0（錨於 ⟨O⟩ 而非 em）")
    emit("  K : FLAG24 ε=0.78(1)(+8/−11) ⟹ em(K⁺−K⁰) = +2.264589，σ=0.120863（M4 docs/36）")
    emit("  N : BMW PRL111,252001 TAB.I，QED(higher-Q 取向 p−n) = +1.59(30)(35) ⟹ n−p 取向 −1.59")
    emit("  Σ : BMW QED(Σ0−Σ−) = +0.08(12)(34) ⟹ Σ−−Σ0 取向 −0.08")
    emit("  Ξ : BMW QED(Ξ0−Ξ−) = −1.29(15)(8) ⟹ Ξ−−Ξ0 取向 +1.29")
    emit("  D : 無 em 錨 ⟹ ⟨O⟩ 不可解（屬盲測目標）")
    emit("")

    ANCHOR_EM = {           # docs/18 取向
        "n-p": (-1.59, 0.4610),
        "K+-K0": (2.264589, 0.120863),
        "S--S0": (-0.08, 0.3610),
        "Xi--Xi0": (1.29, 0.1700),
    }
    o_rows = {}
    emit(f"{'channel':10s}{'coeff':>6s}{'Δm':>9s}{'em(anchor)':>12s}{'σ_em':>8s}"
         f"{'QCD=Δm−em':>12s}{'⟨O⟩':>9s}{'σ⟨O⟩':>8s}  status")
    for name in pairs:
        coeff = DOCS18[name]["coeff"]
        dm = c1_rows[name]["dm"]
        if name == "pi+-pi0":
            O, sO, qcd, em, sem, st = 0.0, 0.0, dm, dm, 0.0, "anchor(Dashen)"
        elif name in ANCHOR_EM:
            em, sem = ANCHOR_EM[name]
            qcd = dm - em
            denom = coeff * float(X)
            O = qcd / denom
            # σ 傳播：⟨O⟩ = qcd/(coeff·X)，σ² = (σ_em/(coeff X))² + (qcd·σ_X/(coeff X²))²
            sO = math.sqrt((sem / denom) ** 2 + (qcd * SX / (coeff * float(X) ** 2)) ** 2)
            st = "anchor"
        else:
            em = sem = qcd = O = sO = None
            st = "BLIND（em 未定）"
        o_rows[name] = dict(coeff=coeff, dm=dm, em=em, sigma_em=sem, qcd=qcd,
                            O=O, sigma_O=sO, status=st)
        if st.startswith("anchor"):
            emit(f"{name:10s}{coeff:6d}{dm:9.4f}{em:12.4f}{sem:8.4f}"
                 f"{qcd:12.4f}{O:9.4f}{sO:8.4f}  {st}")
        else:
            emit(f"{name:10s}{coeff:6d}{dm:9.4f}{'—':>12s}{'—':>8s}{'—':>12s}{'—':>9s}{'—':>8s}  {st}")
    RESULT["C2"] = dict(X=float(X), sigma_X=SX, rows=o_rows)

    # docs/28 式（以 BMW 自身 ΔQCD 歸一）並列，顯示兩法之差
    emit("")
    emit("  並列｜docs/28 式 |ΔQCD^BMW|/(coeff·X)（用 lattice total 而非實驗 total）：")
    bmw = {"n-p": (2.28, 1), "S--S0": (7.67, 2), "Xi--Xi0": (5.87, 1)}
    bmw_rows = {}
    devs = []
    for k, (dQCD, c) in bmw.items():
        v = dQCD / (c * float(X))
        bmw_rows[k] = v
        dv = (o_rows[k]["O"] - v) / v * 100
        devs.append((k, dv))
        emit(f"    {k:10s} |ΔQCD|={dQCD:5.2f}  coeff={c}  ⟹ ⟨O⟩={v:.4f}"
             f"   （本輪實驗 total 法：{o_rows[k]['O']:.4f}，差 {dv:+.1f}%）")
    emit("    ⟹ 兩法之差源於 BMW total 與實驗 Δm 不相容（N: −0.68(39)(36) vs −1.2933）")
    emit("      本記錄一律用「實驗 total − lattice QED」法（data 先行），docs/28 式只作對照")
    RESULT["C2"]["docs28_deviation_pct"] = {k: d for k, d in devs}
    RESULT["C2"]["bmw_normalisation"] = bmw_rows

    # ══════════ C4：可辨識性（精確秩） ══════════
    emit("")
    emit("C4｜可辨識性：未知 = 8 通道 × (⟨O⟩, em) = 16；約束 = 8 條 Δm ＋ 5 個錨")
    names = list(pairs)
    idx = {}
    for i, n in enumerate(names):
        idx[(n, "O")] = 2 * i
        idx[(n, "em")] = 2 * i + 1
    rows = []
    for n in names:
        c = DOCS18[n]["coeff"]
        r = [F(0)] * 16
        r[idx[(n, "O")]] = c * X
        r[idx[(n, "em")]] = F(1)
        rows.append(("data:" + n, r))
    r = [F(0)] * 16
    r[idx[("pi+-pi0", "O")]] = F(1)
    rows.append(("anchor:O_pi=0", r))
    for n, val in (("K+-K0", "2.264589"), ("n-p", "-1.59"),
                   ("S--S0", "-0.08"), ("Xi--Xi0", "1.29")):
        r = [F(0)] * 16
        r[idx[(n, "em")]] = F(1)
        # 以「em 固定值」表示為非齊次；秩只看係數矩陣
        rows.append(("anchor:em_" + n, r))

    def rank_frac(mat):
        m = [row[:] for row in mat]
        nr, nc = len(m), len(m[0])
        rk = 0
        for c in range(nc):
            piv = None
            for i in range(rk, nr):
                if m[i][c] != 0:
                    piv = i
                    break
            if piv is None:
                continue
            m[rk], m[piv] = m[piv], m[rk]
            pv = m[rk][c]
            m[rk] = [v / pv for v in m[rk]]
            for i in range(nr):
                if i != rk and m[i][c] != 0:
                    f = m[i][c]
                    m[i] = [a - f * b for a, b in zip(m[i], m[rk])]
            rk += 1
            if rk == nr:
                break
        return rk

    rk = rank_frac([r for _, r in rows])
    kernel = 16 - rk
    emit(f"  約束矩陣 {len(rows)}×16；rank_ℚ = {rk}；核維 = {kernel}")
    emit(f"  判準：核維 > 0 ⟹ 數據本身不可唯一決定 (⟨O⟩, em) 分解")
    free = [n for n in names if n not in ANCHOR_EM and n != "pi+-pi0"]
    emit(f"  自由方向（每條 1 維）：{', '.join(free)}")
    for n in free:
        emit(f"    {n:10s} 核基：δ em = +1 ⟹ δ⟨O⟩ = −1/(coeff·X) = "
             f"{-1.0 / (DOCS18[n]['coeff'] * float(X)):+.4f}")
    emit(f"  ⟹ 任何 Ô 若只被要求「重現 8 條 Δm」，必然可被平移吸收 ⟹ 存在性問題空泛")
    emit(f"  ⟹ 可證偽結構：Ô 必須同時定出 5 個錨值 ⟨O⟩ 及 3 個自由通道之 (⟨O⟩,em) 分割")
    RESULT["C4"] = dict(n_unknowns=16, n_constraints=len(rows), rank=rk, kernel_dim=kernel,
                        free_channels=free)

    # ══════════ C9：H1 精確箱測試 ══════════
    emit("")
    emit("C9｜H1 可行性（E_ch = 純 em，正核兩體）")
    emit("   形式：em_ch = k₀[ΔΣq²·S + ΔΣ_{i<j}q_iq_j·P]，S,P ∈ (0,1]（M4 docs/26 H2/H3）")
    HC_MEV_FM = 197.3269804
    # a₀ 採 M5 docs/13「重校」值（現行正典；M6 docs/13／M7 重現閘同用）：
    #   v1 重校  a₀ = 0.37742216816391858835 fm ⟹ ħc/a₀ = 522.83 MeV  ← 本輪採用
    #   v2 (docs/28，已被重校取代) a₀ = 0.384184 fm ⟹ ħc/a₀ = 0.5136262 GeV
    # 舊版本腳本誤用 v2 之四捨五入值 514.0 MeV（k₀ 低 1.7%）；見 docs/01 §十 更正。
    A0_FM = 0.37742216816391858835
    ALPHA_INV = 137.035999
    k0 = (HC_MEV_FM / A0_FM) / ALPHA_INV
    emit(f"   k₀ = α ħc/a₀ ；a₀ = {A0_FM:.8f} fm（M5 docs/13 重校）⟹ ħc/a₀ = {HC_MEV_FM/A0_FM:.4f} MeV，"
         f"α⁻¹ = {ALPHA_INV} ⟹ k₀ = {k0:.6f} MeV")
    emit(f"   （對照：v2 docs/28 之 0.513626 GeV ⟹ k₀ = {(HC_MEV_FM/0.384184)/ALPHA_INV:.6f} MeV，已取代）")
    emit("")
    emit(f"{'channel':10s}{'ΔΣq²':>9s}{'ΔΣqq':>9s}{'E_ch':>9s}{'max(箱)':>9s}"
         f"{'t=E/k₀':>9s}{'k_req/k₀':>10s}  verdict")
    c9 = {}
    for name, quarks_fin, quarks_ini in CHANNELS:
        s2f, pf = qsum2(quarks_fin)
        s2i, pi_ = qsum2(quarks_ini)
        dS, dP = s2f - s2i, pf - pi_
        E = c1_rows[name]["E"]
        # 箱上最大／最小（S,P ∈ [0,1]）：線性 ⟹ 四角
        corners = [float(dS) * s + float(dP) * p for s in (0, 1) for p in (0, 1)]
        hi, lo = max(corners), min(corners)
        t = E / k0
        # 可行性：存在 (S,P) ∈ [0,1]² 使 dS·S + dP·P = t
        feasible = (lo - 1e-12) <= t <= (hi + 1e-12)
        # 所需尺度：k_req = E / h*，h* = 最接近目標之箱內值
        h_star = hi if t > 0 else lo
        k_req = abs(E / h_star) if h_star != 0 else float("inf")
        # 符號類：h* 與 t 反號 ⟹ 任何 k>0 都救不到
        sign_blocked = (t > 0 and hi <= 0) or (t < 0 and lo >= 0)
        if sign_blocked:
            verdict = "INFEASIBLE（符號類，k₀ 無關）"
        elif feasible:
            verdict = "feasible"
        else:
            verdict = f"INFEASIBLE（幅度，需 k ≥ {k_req:.3f} = {k_req / k0:.2f}×k₀）"
        emit(f"{name:10s}{float(dS):9.4f}{float(dP):9.4f}{E:9.4f}{hi:9.4f}"
             f"{t:9.4f}{k_req / k0:10.3f}  {verdict}")
        c9[name] = dict(dS=float(dS), dP=float(dP), E=E, box_max=hi, box_min=lo,
                        t=t, feasible=bool(feasible), sign_blocked=bool(sign_blocked),
                        k_req=k_req, k_req_over_k0=k_req / k0, verdict=verdict)
    infeas = [n for n, v in c9.items() if not v["feasible"]]
    sign_b = [n for n, v in c9.items() if v["sign_blocked"]]
    emit("")
    emit(f"  可行性總結：feasible {8 - len(infeas)}/8 ⟹ infeasible = {infeas}")
    emit(f"  其中符號類（與 k₀ 無關、任何正核皆不可能）= {sign_b}")
    emit(f"  ⟹ H1（E_ch = 純 em）在此算子類內被判死；E_ch 不可作 Ô 之目標")
    RESULT["C9"] = dict(k0_MeV=k0, a0_fm=A0_FM, hc_over_a0_MeV=HC_MEV_FM / A0_FM,
                        alpha_inv=ALPHA_INV,
                        a0_source="M5 docs/13 重校（現行正典）；v2 docs/28 之 0.513626 GeV 已取代",
                        a0_alt=dict(v2_a_fm=0.384184, v2_hc_over_a_MeV=HC_MEV_FM / 0.384184,
                                    v2_k0_MeV=(HC_MEV_FM / 0.384184) / ALPHA_INV),
                        channels=c9,
                        infeasible=infeas, sign_blocked=sign_b)

    # ══════════ C11：尺度診斷（a₀ 敏感度） ══════════
    emit("")
    emit("C11｜尺度診斷（條件式，a₀ 屬 M5 輸入）")
    emit("   k_req/k₀ 之意義：兩體正核要解釋該通道，需把 Coulomb 標度放大多少倍")
    ks = [(n, v["k_req_over_k0"]) for n, v in c9.items() if not v["sign_blocked"] and not v["feasible"]]
    if ks:
        emit(f"   幅度不足者：{', '.join(f'{n} {k:.2f}×' for n, k in ks)}")
        emit(f"   最惡 = {max(ks, key=lambda z: z[1])[0]} {max(k for _, k in ks):.2f}×k₀")
    emit("   註：π 之超界比例最小；a₀ 換支（v2 0.384184 ⟺ +1.8%）可令其變動 ⟹ 只 Ξ／Σ−−Σ+ 為穩健超界")
    RESULT["C11"] = dict(required=list(ks),
                         note="a0 = M5 docs/13 重校 0.377422 fm (hc/a0 = 522.83 MeV); conditional")

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(RESULT, f, ensure_ascii=False, indent=2)
    emit("")
    emit(f"[JSON] {OUT_JSON}")
    emit(f"[LINES] {len(log)}")


if __name__ == "__main__":
    main()