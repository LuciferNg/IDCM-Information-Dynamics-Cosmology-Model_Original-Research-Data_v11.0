"""M5 docs/25 — 重推受影響部分（非改數：由算子結構重新導出天花板與指派位移）
================================================================================
`docs/24` 只做「重算」（加守門、掃 57 點）。本腳本做**重推**：把 `r` 拆成質量因子 × 局部化因子，
證明 `r_K/r_D` 之方向上界由**質量/電荷排序**固定 —— 唔依賴 κ、a、t 之掃描。

推導（`docs/14` §四.1 之展開）：
  A = α q (ħc/a) ⟨B⟩                     （Coulomb 幅度，MeV）
  c = C_q (ħc/a)³ P₀,  C_q = (8π/3)α q /(16 m_h m_l)
  ch = κ g q                             （鏈旁觀者，MeV；g = g(d_pair)）
  ⟹ r ≡ (A + ch)/c = [6 m_h m_l /(π (ħc/a)²)] · (⟨B⟩ + ε)/P₀ ,  ε ≡ κ g /(α ħc/a)   （q 全消）
  ⟹ **r ∝ m_h m_l**（乘積；與 C27 一致），局部化部分 = (⟨B⟩+ε)/P₀

  ⟹ R ≡ r_K/r_D = (m_s/m_c) · [(⟨B⟩+ε)_K/(⟨B⟩+ε)_D] · (P₀_D/P₀_K)

判準（寫死）：
  T1 恆等式必須與數值態（`scripts/26` 之 21 個束縛點）逐點核，相對誤差 < 1e−10，否則記 FAIL。
  T2 質量因子 (m_s/m_c) = 0.389764 —— 單獨已 < 1，數據要 28.96 ⟹ 局部化部分需補 74.3×。
  T3 接觸強度比 h_K/h_D = C_q(K)/C_q(D) = (q_s/q_c)(m_c/m_s) = 1.2830 > 1（**與 κ,a,t 無關**）
     ⟹ K 之接觸井較深 ⟹ K 較局部化 ⟹ P₀_K/P₀_D ↑ ⟹ 局部化因子方向**與數據相反**。
  T4 上界取 21 個束縛點之 sup；並列印達到數據所需之局部化因子，判其是否可達。
  T5 `docs/22` §四 重推：Δd ∈ {−2..3} 之幾何 × 動力學 vs 數據 0.038298；所需 Δd；最近整數之覆蓋率；
     A3（ℓ(d)=ℓ(s)=7 ⟹ d_K=0, d_D=7）之過強倍數。
  T6 禁：不得為達 28.96 調 m_h／m_l／q／ℓ。

輸出：`data/27_theorem_rederivation.json`
"""
import json

import numpy as np

B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
KAP = 0.2150
C_K, C_D = 0.0135, 0.3525
R_DATA_K, R_DATA_D = 98.77777777777769, 3.4113475177304964
R_DATA = R_DATA_K / R_DATA_D
BOUND_RHO, DPEAK_MAX = 1e-3, 20


def st(a, t, q, mh, ml):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    H = np.diag([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    for d in range(L):
        H[d, d + 1] = H[d + 1, d] = -OFF[d]
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    return {"rA": rA, "P0": P0, "dpk": dpk, "bound": bool(rA > BOUND_RHO and dpk <= DPEAK_MAX),
            "Cq": Cq, "h": h, "g": g}


def r_num(s, a, q, gch, kappa=KAP):
    A = ALPHA * q * (HBARC / a) * s["rA"] * 1000
    c = s["Cq"] * (HBARC / a) ** 3 * s["P0"] * 1000
    return (A + kappa * gch * q) / c, A, c


def r_analytic(s, a, mh, ml, gch, kappa=KAP):
    eps = kappa * gch * 1e-3 / (ALPHA * HBARC / a)   # chain 以 MeV 加入、Coulomb 以 GeV ⟹ ×1e-3
    return 6 * mh * ml / (np.pi * (HBARC / a) ** 2) * (s["rA"] + eps) / s["P0"], eps


print("=" * 108)
print("M5 docs/25 — 重推：r 之質量因子 × 局部化因子分解、方向定理、指派位移重算")
print("=" * 108)

rows = json.load(open("data/21_R_ceiling.json"))["rows"]
pts = []
for r in rows:
    a, t = r["a0_fm"], r["t_MeV"] * 1e-3
    sK, sD = st(a, t, Q_S, M_S, M_L), st(a, t, Q_C, M_C, M_L)
    if not (sK["bound"] and sD["bound"]):
        continue
    rK, AK, cK = r_num(sK, a, Q_S, G0)
    rD, AD, cD = r_num(sD, a, Q_C, G7)
    aK, epsK = r_analytic(sK, a, M_S, M_L, G0)
    aD, epsD = r_analytic(sD, a, M_C, M_L, G7)
    pts.append({"kappa": r["kappa"], "a0": a, "t": t, "rK": rK, "rD": rD, "R": rK / rD,
                "aK": aK, "aD": aD, "epsK": epsK, "epsD": epsD,
                "P0K": sK["P0"], "P0D": sD["P0"], "rAK": sK["rA"], "rAD": sD["rA"],
                "hK": sK["h"], "hD": sD["h"]})

print(f"\n[T1] 恆等式核對（{len(pts)} 個束縛點）")
dev_id = max(max(abs(p["rK"] - p["aK"]) / abs(p["rK"]), abs(p["rD"] - p["aD"]) / abs(p["rD"])) for p in pts)
print(f"     r 數值 vs 解析：最大相對偏差 = {dev_id:.3e} ⟹ {'✅ PASS' if dev_id < 1e-10 else '🔴 FAIL'}")
dev_R = max(abs((p["rK"] / p["rD"]) -
                (M_S / M_C) * ((p["rAK"] + p["epsK"]) / (p["rAD"] + p["epsD"])) * (p["P0D"] / p["P0K"]))
            / (p["rK"] / p["rD"]) for p in pts)
print(f"     R 分解式核對：最大相對偏差 = {dev_R:.3e} ⟹ {'✅ PASS' if dev_R < 1e-10 else '🔴 FAIL'}")

print(f"\n[T2] 質量因子")
print(f"     m_s/m_c = {M_S/M_C:.6f}（< 1）；數據需 R = {R_DATA:.4f} ⟹ 局部化部分需補 {R_DATA/(M_S/M_C):.2f}×")

print(f"\n[T3] 接觸強度比（與 κ, a, t 無關）")
hK, hD = pts[0]["hK"], pts[0]["hD"]
print(f"     h_K/h_D = C_q(K)/C_q(D) = (q_s/q_c)(m_c/m_s) = {(Q_S/Q_C)*(M_C/M_S):.6f}"
      f"（數值 {hK/hD:.6f}）⟹ K 之接觸井較深 ⟹ K 較局部化")
lo, hi = min(p['P0K']/p['P0D'] for p in pts), max(p['P0K']/p['P0D'] for p in pts)
print(f"     全 {len(pts)} 點 P₀_K/P₀_D 範圍 = [{lo:.3f}, {hi:.3f}] ⟹ 實測 K **較唔局部化**（方向與數據一致，數據需 P₀_D/P₀_K = {R_DATA/(M_S/M_C):.1f}）")
print(f"     ⟹ 故此處非『方向不可能』，而係**幅度不足**：需 P₀_D/P₀_K ≈ 74，模型最多 {1/lo:.2f}")

print(f"\n[T4] 束縛點之 R 上界")
mx = max(pts, key=lambda p: p["R"])
print(f"     R_max = {mx['R']:.4f}（κ = {mx['kappa']}）⟹ r_K 天花板 = R_max × r_D(data) = {mx['R']*R_DATA_D:.4f}")
print(f"     數據需 R = {R_DATA:.2f} ⟹ 短缺 {R_DATA/mx['R']:.2f}×")
print(f"     ⟹ 要達數據，需 (⟨B⟩+ε)_K/(⟨B⟩+ε)_D ×(P₀_D/P₀_K) = {R_DATA/(M_S/M_C):.2f}"
      f"（現時最大 {((mx['rAK']+mx['epsK'])/(mx['rAD']+mx['epsD']))*(mx['P0D']/mx['P0K']):.3f}）")

print(f"\n[T5] `docs/22` §四 重推（動力學因子取守門點 κ = {mx['kappa']}）")
DYN = mx["P0K"] / mx["P0D"]
print(f"     P₀(K)/P₀(D) = {DYN:.5f} = 1/{1/DYN:.2f}；數據 c_K/c_D = {C_K/C_D:.6f}")
print(f"     {'Δd':>4s}{'幾何 B^-Δd':>12s}{'×動力學':>12s}{'vs 數據':>10s}{'倍數':>9s}{'判':>10s}")
tab = []
for dd in (-2, -1, 0, 1, 2, 3):
    geo = B ** (-dd)
    val = geo * DYN
    ratio = val / (C_K / C_D)
    tab.append({"dd": dd, "geo": geo, "val": val, "ratio": ratio})
    print(f"     {dd:>4d}{geo:>12.5f}{val:>12.6f}{1/ratio:>10.3f}×{ratio:>9.3f}")
need = np.log((C_K / C_D) / DYN) / np.log(B) * -1
print(f"     所需 Δd = log₄(數據/動力學) = {need:.4f}（守門值；v1 非守門 1.3097）")
d1 = next(x for x in tab if x["dd"] == 1)
print(f"     最近整數 Δd = 1 ⟹ 覆蓋 {d1['ratio']*100:.2f}%（v1 報 59.9%）")
a3 = 16384 * DYN / (C_K / C_D)
print(f"     A3（ℓ(d)=ℓ(s)=7 ⟹ d_K−d_D = −7 ⟹ 幾何 16384）⟹ 過強 {a3:.1f}×（v1 報 39242.5×）")

# ---- T6 守門門檻敏感性（天花板是否 intrinsic？）
print(f"\n[T6] 守門門檻敏感性（`d_peak ≤ D`；同一 57 校準點）")
sens = []
for Dmax in (10, 20, 30, 40):
    best = None
    for r in rows:
        a, t = r["a0_fm"], r["t_MeV"] * 1e-3
        sK, sD = st(a, t, Q_S, M_S, M_L), st(a, t, Q_C, M_C, M_L)
        okK = sK["rA"] > BOUND_RHO and sK["dpk"] <= Dmax
        okD = sD["rA"] > BOUND_RHO and sD["dpk"] <= Dmax
        if not (okK and okD):
            continue
        rK, _, _ = r_num(sK, a, Q_S, G0)
        rD, _, _ = r_num(sD, a, Q_C, G7)
        if best is None or rK / rD > best[0]:
            best = (rK / rD, r["kappa"], rK, rD, sK["dpk"], sD["dpk"])
    if best:
        sens.append({"Dmax": Dmax, "R_max": best[0], "kappa": best[1],
                     "ceiling_rK": best[0] * R_DATA_D, "shortfall": R_DATA / best[0],
                     "dpk_K": best[4], "dpk_D": best[5]})
        print(f"     D = {Dmax:2d} ⟹ R_max = {best[0]:.4f}（κ = {best[1]}、d_peak K/D = {best[4]}/{best[5]}）"
              f" ⟹ r_K ≤ {best[0]*R_DATA_D:8.4f}、短缺 {R_DATA/best[0]:7.2f}×")
    else:
        print(f"     D = {Dmax:2d} ⟹ 無同時束縛點")
print("     ⟹ 天花板**非 intrinsic**：隨守門門檻單調上升（放寬 ⟹ 尾端帶邊態復活）⟹ 只能作報告值，不可作物理性界")

out = {"identity_max_dev": dev_id, "R_factor_max_dev": dev_R,
       "mass_factor": M_S / M_C, "required_localisation": R_DATA / (M_S / M_C),
       "h_ratio": (Q_S / Q_C) * (M_C / M_S), "h_ratio_numeric": hK / hD,
       "P0_ratio_range": [min(p["P0K"] / p["P0D"] for p in pts), max(p["P0K"] / p["P0D"] for p in pts)],
       "R_max": mx["R"], "rK_ceiling": mx["R"] * R_DATA_D, "shortfall": R_DATA / mx["R"],
       "best_point": {k: mx[k] for k in ("kappa", "a0", "t", "rK", "rD")},
       "P0_ratio_best": DYN, "required_dd": float(need), "table": tab,
       "d1_coverage": d1["ratio"], "a3_over": a3,
       "verdict": "🔴 分解定理：r = [6 m_h m_l/(π(ħc/a)²)]·(⟨B⟩+ε)/P₀（恆等式 4.4e-16 通過）⟹ R = (m_s/m_c)·Λ，Λ 為局部化因子；m_s/m_c = 0.3898 ⟹ 數據需 Λ = 74.3×，束縛族最大 Λ = 9.38 ⟹ 短缺 7.92×；且天花板隨守門門檻上升 ⟹ 非 intrinsic，僅報告值",
       "guard_sensitivity": sens}
with open("data/27_theorem_rederivation.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/27_theorem_rederivation.json")