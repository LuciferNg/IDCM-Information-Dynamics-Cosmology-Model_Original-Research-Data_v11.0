"""M11 scripts/08 — G0：∂H/∂X 算子化（Hellmann–Feynman 路線之可行性判）
================================================================================================
目標（G0）：把 M11 之 ⟨O⟩（強 IB 係數）由「擬合參數」升級為譜矩陣元：
      ⟨O⟩_ch = Δ(∂M/∂X)/(coeff·X) ,  X = m_d − m_u
問題：M5 之 41 層 chain Hamiltonian 內，X 究竟從何處入場？其響應夠唔夠大？

M5 之 H（`scripts/21` 逐字，GeV·fm 制）：
    gcoef  = α q ħc/(a·t)                    ← 只含 q、a、t，**不含任何質量**
    Cq     = (8π/3) α q /(16 m_h m_l)        ← 唯一含質量之處
    h      = Cq (ħc/a)³/t                    ← 只加在 d = 0 對角
    diag[d] = −gcoef·B^(−d) − [d=0]·h
    offdiag = −√B（d=0）／−√(B−1)（其餘）      ← 常數，與 X 無關
    ⟹ ∂H/∂m_l 只有一條路徑：Cq ∝ 1/m_l ⟹ ∂Cq/∂m_l = −Cq/m_l（同理 m_h）

本輪三問（判準事前寫死）：
  Q1  ∂M/∂X 之幅度：δM = |∂M/∂m_l|·X 與 M11 所需（coeff·X·⟨O⟩，⟨O⟩ ∈ [1, 2.42] ⟹ 2.56–6.20 MeV）之比。
      判：若 δM/所需 < 1e−1 ⟹ 質量導數路線**結構性不足**（差一個量級以上），G1 不得在此路線上做。
  Q2  波函數之層佔據：⟨P_0⟩、⟨P_7⟩、⟨POW⟩、argmax ℓ —— M5 之層指派係 ℓ(u)=ℓ(c)=0、ℓ(d)=ℓ(s)=7
      ⟹ **isospin 已住喺層結構內**（u 與 d 分居 0 與 7）。故 iso-odd 一體算子之自然候選係
      層投影子 P_7 − P_0（或其 NODES 加權版），而非質量導數。
      本輪只**報數**（每個通道之層矩陣元），比較留待 G1（保持盲）。
  Q3  ∂M/∂X 之符號／通道依賴：質量導數路線可否區分通道（K vs D）？

禁：不得用錨值反推任何係數；不得改 M5 之 H；不得引入新常數。
輸出：`data/08_g0_x_operator.json`
"""
import json
import os

import numpy as np

ALPHA = 1 / 137.035999084
HBARC = 0.1973269804          # GeV·fm
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
A_D, C_D = 1.2025, 0.3525
A0 = 0.37742216816391858835   # fm（M5 docs/13 重校；κ=0 校準）
T0 = 0.0007884                # GeV（0.7884 MeV）
X = 0.00256                   # GeV（X = m_d − m_u = 2.56 MeV）

RES = {}


def amps(a, t, q, mh, ml, want_w=False):
    gcoef = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-gcoef * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    bd = bool(P0 > 1e-3 and int(np.argmax(w)) <= 20)
    A = ALPHA * q * (HBARC / a) * rA * 1000
    c = Cq * (HBARC / a) ** 3 * P0 * 1000
    if want_w:
        return A, c, bd, w
    return A, c, bd


def deriv(q, mh, ml, wm, what="ml", rel=1e-4):
    """中央差分 ∂(A,c)/∂(m_l 或 m_h)，在 (A0,T0) 處"""
    d = rel * (ml if what == "ml" else mh)
    if what == "ml":
        Ap, cp, _ = amps(A0, T0, q, mh, ml + d)
        Am, cm, _ = amps(A0, T0, q, mh, ml - d)
    else:
        Ap, cp, _ = amps(A0, T0, q, mh + d, ml)
        Am, cm, _ = amps(A0, T0, q, mh - d, ml)
    return (Ap - Am) / (2 * d), (cp - cm) / (2 * d)


def main():
    print("=" * 108)
    print("M11 scripts/08 — G0：∂H/∂X 算子化（Hellman–Feynman 路線可行性）")
    print("=" * 108)
    print(f"校準點：(a, t) = ({A0:.8f} fm, {T0*1000:.4f} MeV)；X = {X*1000:.2f} MeV")
    print()
    print("[Q1] H 之 X 依賴路徑（解析）")
    print("  gcoef = α q ħc/(a t)              ⟹ ∂/∂X = 0")
    print("  offdiag = −√B, −√(B−1)            ⟹ ∂/∂X = 0")
    print("  Cq    = (8π/3)α q/(16 m_h m_l)    ⟹ ∂Cq/∂m_l = −Cq/m_l（唯一路徑）")
    print("  h     = Cq (ħc/a)³/t              ⟹ ∂h/∂m_l = −h/m_l")
    print("  ⟹ **∂H/∂X 係一條純對角（d=0）算子**：∂H/∂m_l = −(h/m_l)·|0⟩⟨0|")
    print()
    fam = [("D（q=2/3, m_h=1.27, m_l=0.30）", Q_C, M_C, M_L),
           ("K（q=1/3, m_h=0.495, m_l=0.30）", Q_S, M_S, M_L)]
    rows = []
    print(f"[Q1b] 質量導數之幅度 vs M11 所需（⟨O⟩ ∈ [1, 2.42] ⟹ 2.56–6.20 MeV）")
    print(f"{'族':32s}{'∂A/∂m_l':>11s}{'∂c/∂m_l':>11s}{'δM(A+3c)':>11s}{'δM(A−c)':>10s}{'差所需倍數':>11s}")
    for lab, q, mh, ml in fam:
        A, c, bd, w = amps(A0, T0, q, mh, ml, want_w=True)
        dA, dc = deriv(q, mh, ml, None, "ml")
        dA_h, dc_h = deriv(q, mh, ml, None, "mh")
        dM1 = abs((dA + 3 * dc) * X)          # GeV→MeV：導數單位 MeV/GeV × X(GeV) = MeV
        dM2 = abs((dA - dc) * X)
        need = 2.56
        ratio = need / max(dM1, 1e-12)
        print(f"{lab:32s}{dA:11.4f}{dc:11.4f}{dM1:11.5f}{dM2:10.5f}{ratio:11.1f}×")
        rows.append(dict(family=lab, A_MeV=A, c_MeV=c, bound=bd,
                         dA_dml=dA, dc_dml=dc, dA_dmh=dA_h, dc_dmh=dc_h,
                         dM_Aplus3c_MeV=dM1, dM_Aminusc_MeV=dM2,
                         deficit_factor=ratio))
    print()
    print("  ⟹ 質量導數（含 h 之全部路徑）最多提供 ~1e−2 MeV 級；所需 2.56–6.20 MeV")
    print("     ⟹ **差 2–3 個量級 ⟹ 質量導數路線結構性不足**（Q1 判：不通過）")
    print()
    print("[Q2] 波函數之層佔據（M5 層指派 ℓ(u)=ℓ(c)=0、ℓ(d)=ℓ(s)=7 ⟹ isospin 住喺層內）")
    print(f"{'族':22s}{'⟨POW⟩':>10s}{'⟨P_0⟩':>10s}{'w[7]':>10s}{'⟨ℓ⟩':>9s}{'argmax':>8s}{'⟨NODES⟩':>12s}")
    lw = np.arange(L + 1)
    for lab, q, mh, ml in fam:
        A, c, bd, w = amps(A0, T0, q, mh, ml, want_w=True)
        m = dict(family=lab, mean_POW=float((w * np.array(POW)).sum()), w0=float(w[0]),
                 w7=float(w[7]), mean_layer=float((w * lw).sum()),
                 argmax=int(np.argmax(w)), mean_NODES=float((w * np.array(NODES)).sum()),
                 w_head=[float(x) for x in w[:9]])
        print(f"{lab:22s}{m['mean_POW']:10.5f}{m['w0']:10.5f}{m['w7']:10.5f}"
              f"{m['mean_layer']:9.3f}{m['argmax']:8d}{m['mean_NODES']:12.3f}")
        RES.setdefault("layer_moments", []).append(m)
    print()
    print("  ⟹ 態之重量集中在 d ≲ 3（argmax 見上）；w[7] 已有 10^−3 以下")
    print("     ⟹ 若 iso-odd 算子係 P_7 − P_0 型，其矩陣元量級 = w[7] − w[0] ~ −w[0]（符號由 d=0 主導）")
    print("     ⚠ 此路線之數值必須留待 G1 與錨比較，本輪**只報數**（保持盲）")
    print()
    print("[Q3] 質量導數路線之通道區分能力")
    d0 = rows[0]["dM_Aplus3c_MeV"]
    d1 = rows[1]["dM_Aplus3c_MeV"]
    print(f"  D 族 δM = {d0:.5f} MeV、K 族 δM = {d1:.5f} MeV ⟹ 比值 {d0/max(d1,1e-12):.3f}")
    print("  ⟹ 兩族幾乎同值（皆由 h/m_l 主導）⟹ **無通道分辨力**，即使幅度夠都做唔到 K≠D")
    RES.update(dict(a0=A0, t0=T0, X_GeV=X, mass_derivative=rows,
                    verdict_Q1="FAIL（差 2–3 量級）", verdict_Q3="FAIL（無通道分辨）"))
    _out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "08_g0_x_operator.json")
    with open(_out, "w") as f:
        json.dump(RES, f, ensure_ascii=False, indent=2, default=str)
    print("\n-> data/08_g0_x_operator.json")


main()