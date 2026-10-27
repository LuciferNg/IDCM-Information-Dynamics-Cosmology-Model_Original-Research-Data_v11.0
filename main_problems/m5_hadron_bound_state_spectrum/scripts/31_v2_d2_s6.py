"""M5 docs/29 配套 — v2 約定下之 D2（強 IB）與 S6（contact 下限）重算
================================================================================
`docs/28` 已重算 A2／A3／兩扇區／D5。仍有兩條 load-bearing 否定結論建立於 v1 校準點：
  D2（`docs/16`）：K 系殘餘係唔係強 IB？—— 含 Q1 非線性殘餘、Q2 自旋盲、Q5 硬上界 A_max。
  S6（`docs/18`）：contact 係數之結構下限 min|c|（固定 A_K 掃 (a_C, t)）。
本腳本以 v2 校準（a₀ = 0.384184 fm、t = 1.7711 MeV）重算二者，其餘判準照舊。

判準（寫死）：
  X1 Q1：Δ_nl(K) 須 ≥ 0.687 MeV（= 1.374 之 50%）才算 D2 成立；否則記「差 N×」。
  X2 Q2：|E(K⁰K⁻) − E(K*⁰K*⁻)| 之 σ 數；≤ 2σ 記自旋盲成立。
  X3 Q5：A_max = α q(ħc/a_C)（⟨B⟩ ≤ 1）；與 A_req = 1.3335 ± 0.2570 比，報超界百分比與 σ。
  X4 S6：對兩個 A 目標（1.3335、1.3335 − chain）掃 (a_C, t) 求 min|c|；報 4|c| 對 0.0540 之倍數。
  X5 不得改 q、m_h、m_l、κ、D 之目標；不得放寬束縛守門。
輸出：`data/31_v2_d2_s6.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
A0, T0 = 0.384184, 1.7711e-3          # v2 校準（docs/27 §三／docs/28 §二）
DL = 1.28                              # X/2（MeV）
KAPPA = 0.26986
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/29 — v2 約定下之 D2（強 IB）與 S6（contact 下限）")
print("=" * 104)


def shape(aC, t, q, mh, ml):
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aC) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    bd = bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)
    A = ALPHA * q * (HBARC / aC) * rA * 1000
    c = Cq * (HBARC / aC) ** 3 * P0 * 1000
    return A, c, A + 3 * c, {"rA": rA, "P0": P0, "dpk": dpk, "bound": bd}


# ---------------- Q1 ----------------
print("\n[Q1] 非線性強 IB 殘餘（同一核，m_l ± δ；δ = X/2 = 1.28 MeV）")
print(f"     {'系統':4s}{'E(m_l−δ)':>12s}{'E(m_l+δ)':>12s}{'Δ_nl':>10s}{'佔 1.374':>10s}{'對 0.687':>10s}")
q1 = {}
for nm, q, mh in (("K", Q_S, M_S), ("D", Q_C, M_C)):
    _, _, Em, _ = shape(A0, T0, q, mh, M_L - DL / 1000)
    _, _, Ep, _ = shape(A0, T0, q, mh, M_L + DL / 1000)
    dnl = Ep - Em
    q1[nm] = {"Em": Em, "Ep": Ep, "dnl": dnl, "frac": dnl / 1.374}
    print(f"     {nm:4s}{Em:12.4f}{Ep:12.4f}{dnl:10.5f}{dnl/1.374*100:9.3f}%{abs(dnl)/0.687*100:9.1f}%")
short = 0.687 / abs(q1["K"]["dnl"])
X1 = abs(q1["K"]["dnl"]) >= 0.687
print(f"     ⟹ |Δ_nl(K)| = {abs(q1['K']['dnl']):.4f} MeV = 所需 1.374 之 {abs(q1['K']['dnl'])/1.374*100:.2f}%"
      f" ⟹ {'✅ 成立' if X1 else f'🔴 差 {short:.0f}×'}（v1 為 1.46%、34×）")

# ---------------- Q2 ----------------
print("\n[Q2] 自旋盲（強 IB 不辨自旋）：觀測 1.3740 vs 1.3200，差 0.0540 ± 0.3545")
d_spin, sd = 1.3740 - 1.3200, 0.3545
z_spin = d_spin / sd
X2 = abs(z_spin) <= 2
print(f"     z = {z_spin:+.3f} ⟹ {'✅ 自旋盲成立' if X2 else '🔴'}")
print("     （此判準只用觀測值，與 weight 約定無關 ⟹ v1／v2 相同）")

# ---------------- Q5 ----------------
print("\n[Q5] 硬上界：⟨B⟩ ≤ 1 ⟹ A_max = α q(ħc/a_C)")
Amax = ALPHA * Q_S * (HBARC / A0) * 1000
Areq, sAreq = 1.3335, 0.2570
X3 = Areq <= Amax
print(f"     A_max = {Amax:.4f} MeV（v1：1.2718）；A_req = {Areq:.4f} ± {sAreq:.4f}")
print(f"     A_req/A_max = {Areq/Amax:.4f} ⟹ 中心值超界 {(Areq/Amax-1)*100:.2f}%（需 ⟨B⟩ = {Areq/Amax:.4f} > 1）；"
      f"統計上 {(Areq-Amax)/sAreq:+.2f}σ")
print(f"     ⟹ {'🟡 未超界' if X3 else '🔴 超界：本核 em 對 K 系電荷分裂封頂'}")

# ---------------- S6 ----------------
print("\n[S6] contact 係數結構下限：固定 A_K 掃 (a_C, t) 求 min|c|")
chain_K = KAPPA * G0 * Q_S
floor = {}
for A_tgt, tag in ((1.3335, "A_K 全由 tree 提供"), (1.3335 - chain_K, "扣 chain 後之 tree 殘量")):
    best, arg = None, None
    for la in np.linspace(np.log(1e-3), np.log(3.0), 120):
        a = float(np.exp(la))
        for lt in np.linspace(np.log(1e-4), np.log(1e3), 160):
            t = float(np.exp(lt))
            A, c, _, st = shape(a, t, Q_S, M_S, M_L)
            if not st["bound"] or abs(A - A_tgt) > 5e-3:
                continue
            if best is None or c < best:
                best, arg = c, (a, t * 1000, A, st["P0"])
    floor[tag] = None if best is None else {"c_min": best, "a0_fm": arg[0], "t_MeV": arg[1],
                                            "A": arg[2], "P0": arg[3], "ratio": best / 0.0135,
                                            "four_c": 4 * best, "over": 4 * best / 0.0540}
    v = floor[tag]
    if v:
        print(f"     A = {A_tgt:.4f}（{tag}）⟹ min|c| = {v['c_min']:.5f} MeV（a₀ = {v['a0_fm']:.5f} fm、"
              f"t = {v['t_MeV']:.4f} MeV、P₀ = {v['P0']:.2e}）vs 需 0.0135 ⟹ {v['ratio']:.1f}×")
    else:
        print(f"     A = {A_tgt:.4f}（{tag}）⟹ 無束縛解")
print(f"     鏈值 chain_K（κ = {KAPPA}、g(0) = {G0:.4f}、q = 1/3）= {chain_K:.4f} MeV")
print("     ⟹ E1 − E2 ≡ 4|c|：")
for k, v in floor.items():
    if v:
        print(f"       {k}: 4|c| ≥ {v['four_c']:.4f} MeV vs 觀測 0.0540 ⟹ 超 {v['over']:.1f}×")

out = {"calibration": {"a0_fm": A0, "t_MeV": T0 * 1000},
       "Q1": q1, "Q1_pass": bool(X1), "short_factor": short,
       "Q2": {"z": z_spin, "pass": bool(X2)},
       "Q5": {"A_max": Amax, "A_req": Areq, "sigma": sAreq, "over_pct": (Areq / Amax - 1) * 100,
              "sigma_level": (Areq - Amax) / sAreq, "pass": bool(X3)},
       "S6": floor, "chain_K": chain_K,
       "v1_reference": {"Q1_frac_pct": 1.46, "Q1_short": 34, "A_max": 1.2718,
                        "S6_c_min_full": 1.49187, "S6_over_full": 110.5,
                        "S6_c_min_chain": 0.27593, "S6_over_chain": 20.4}}
with open("data/31_v2_d2_s6.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/31_v2_d2_s6.json")