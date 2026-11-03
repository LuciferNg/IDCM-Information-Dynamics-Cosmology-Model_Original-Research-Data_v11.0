"""M11 scripts/09 — 路線 C：放棄 em／strong 分解，只測 8 通道**總量** Δm
================================================================================================
背景（`docs/10 §五 C`）：M11 之全部分析依賴 `Δm = strong + em` 之分解，而該分解係 LO χPT 約定
（FLAG 2024 自認），且 ε 有 QED quenching ~0.15（`docs/09`）。路線 C 直接測總量，避開約定。

模型（**結構固定，只有標度自由**）：
    Δm_ch = a·Φ_ch + b·Ψ_ch
    Φ_ch ≡ ½[ΔΣq² + ΔΣ_{i<j}q_iq_j]    ← 電荷結構（由夸克電荷精確算，無自由形狀參數）
    Ψ_ch ≡ coeff_ch                     ← 價夸克結構（整數，u↔d 取代數）
    a = k₀ 之標度、b = 強 IB 標度（結構期望 b = X = 2.56 MeV，即 ⟨O⟩ ≡ 1）

**判準（事前寫死，見下；跑之前已定，不得事後改）**：

  P1 主要（剛性，2 參數／8 數據 ⟹ 6 dof）：
       PASS ⟺ χ²/dof ≤ 2 且 8 條全部 |z| ≤ 2 且 a > 0 且 ⟨O⟩_eff ≡ b/X ∈ [0.5, 4]
       KILL ⟺ 任一 |z| ≥ 3，或 a ≤ 0，或 ⟨O⟩_eff ∉ [0.3, 5]
  P2 敏感度（形狀自由，4 參數 a,b,S,P ∈ (0,1]）：只報，不作判詞
  P3 對照（b 強制 = X，即 ⟨O⟩≡1，只擬 a）：報 χ²，用來判斷 P1 之改善係真結構抑或只係多一個標度

σ 用 **PDG 質量誤差**（純實驗）；**不用** σ_X ——因為路線 C 不把 X 當輸入，若混入 σ_X 會雙重計數。

禁：不得用錨值反推；不得事後加通道；不得改 σ；不得改 Φ/Ψ 之結構定義。
輸出：`data/09_route_c_total_dm.json`
"""
import json
import os

import numpy as np

X = 2.56
K0 = 3.8152622441

# PDG 質量（MeV）與誤差 —— 逐字照 scripts/01
M = {"p": (938.27208943, 2.9e-7), "n": (939.56542052, 5.4e-7),
     "K+": (493.677, 0.016), "K0": (497.611, 0.013),
     "pi+": (139.57039, 1.8e-4), "pi0": (134.9768, 5e-4),
     "S-": (1197.449, 0.030), "S0": (1192.642, 0.024), "S+": (1189.37, 0.07),
     "Xi-": (1321.71, 0.07), "Xi0": (1314.86, 0.20),
     "D+": (1869.66, 0.05), "D0": (1864.84, 0.05)}

CH = [("n-p", "n", "p", 1), ("K+-K0", "K+", "K0", -1), ("pi+-pi0", "pi+", "pi0", 1),
      ("S--S0", "S-", "S0", 1), ("S0-S+", "S0", "S+", 1), ("S--S+", "S-", "S+", 2),
      ("Xi--Xi0", "Xi-", "Xi0", 1), ("D+-D0", "D+", "D0", 1)]

QB = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3}
CONTENT = {
    "n-p": ([("d", 2), ("u", 1)], [("u", 2), ("d", 1)]),
    "K+-K0": ([("u", 1), ("sbar", 1)], [("d", 1), ("sbar", 1)]),
    "pi+-pi0": ([("u", 1), ("dbar", 1)], "PI0"),
    "S--S0": ([("d", 2), ("s", 1)], [("u", 1), ("d", 1), ("s", 1)]),
    "S0-S+": ([("u", 1), ("d", 1), ("s", 1)], [("u", 2), ("s", 1)]),
    "S--S+": ([("d", 2), ("s", 1)], [("u", 2), ("s", 1)]),
    "Xi--Xi0": ([("d", 1), ("s", 2)], [("u", 1), ("s", 2)]),
    "D+-D0": ([("c", 1), ("dbar", 1)], [("c", 1), ("ubar", 1)]),
}


def q2(spec):
    if spec == "PI0":
        a = q2([("u", 1), ("ubar", 1)])
        b = q2([("d", 1), ("dbar", 1)])
        return (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
    qs = [QB[k.rstrip("bar")] * (-1 if k.endswith("bar") else 1) for k, n in spec for _ in range(n)]
    s2 = sum(q * q for q in qs)
    pr = sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs)))
    return s2, pr


def wfit(A, y, s, extra=()):
    """加權最小平方；extra = [(col, fixed_value)] 之強制列"""
    n = A.shape[0]
    W = 1 / s
    Aw, yw = A * W[:, None], y * W
    if extra:
        rows = np.array([[v for _, v in extra]])
        Aw = np.vstack([Aw, rows * 1e8])
        yw = np.concatenate([yw, np.zeros(len(rows))])
    p, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    r = (y - A @ p) / s
    return p, float((r ** 2).sum()), r


def main():
    rows, P, Y, S, PSI = [], [], [], [], []
    print("=" * 104)
    print("M11 scripts/09 — 路線 C：8 通道總量 Δm 之 2 標度檢定（脫離 em/strong 分解）")
    print("=" * 104)
    print(f"{'通道':9s}{'Δm(MeV)':>10s}{'σ_PDG':>9s}{'Φ0':>9s}{'Ψ=coeff':>9s}{'strong=X·Ψ':>12s}")
    for name, hi, lo, coeff in CH:
        dm = M[hi][0] - M[lo][0]
        sd = (M[hi][1] ** 2 + M[lo][1] ** 2) ** 0.5
        a2, ap = q2(CONTENT[name][0])
        b2, bp = q2(CONTENT[name][1])
        phi = ((a2 - b2) + (ap - bp)) / 2
        rows.append(dict(ch=name, dm=dm, sd=sd, phi=phi, coeff=coeff))
        P.append(phi); Y.append(dm); S.append(sd); PSI.append(coeff)
        print(f"{name:9s}{dm:10.4f}{sd:9.4f}{phi:9.4f}{coeff:9d}{X*coeff:12.4f}")
    A = np.array([P, PSI]).T
    Y, S = np.array(Y), np.array(S)

    print()
    print("[P1] 剛性 2 參數擬合（a ≡ k₀ 標度、b ≡ 強 IB 標度）—— 6 dof")
    p, chi2, r = wfit(A, Y, S)
    a, b = p
    dof = len(Y) - 2
    print(f"  a = {a:.4f} MeV、b = {b:.4f} MeV ⟹ ⟨O⟩_eff = b/X = {b/X:.4f}")
    print(f"  χ² = {chi2:.2f}、dof = {dof}、χ²/dof = {chi2/dof:.2f}")
    print(f"  {'通道':9s}{'Δm_obs':>10s}{'Δm_pred':>10s}{'殘差':>9s}{'z':>8s}")
    ok = True
    for i, rw in enumerate(rows):
        pred = a * rw["phi"] + b * rw["coeff"]
        print(f"  {rw['ch']:9s}{rw['dm']:10.4f}{pred:10.4f}{rw['dm']-pred:9.4f}{r[i]:8.2f}")
        rw.update(pred=pred, z=float(r[i]))
        if abs(r[i]) >= 2:
            ok = False
    P1 = dict(a=a, b=b, O_eff=b / X, chi2=chi2, dof=dof, chi2_dof=chi2 / dof,
              all_z_lt2=ok,
              PASS=bool(ok and chi2 / dof <= 2 and a > 0 and 0.5 <= b / X <= 4),
              KILL=bool((not ok) or a <= 0 or not (0.3 <= b / X <= 5)))
    print(f"  ⟹ 判：**{'PASS' if P1['PASS'] else 'KILL'}**"
          f"（|z|≤2 全過={ok}、χ²/dof≤2={chi2/dof<=2}、a>0={a>0}、⟨O⟩_eff∈[0.5,4]={0.5<=b/X<=4}）")

    print()
    print("[P3] 對照：b 強制 = X（即 ⟨O⟩≡1，只擬 a）—— 呢個就係 C9 之總量版")
    pa, chi2a, ra = wfit(A[:, :1], Y - X * np.array(PSI), S)
    print(f"  a = {pa[0]:.4f} MeV（正典 k₀ = {K0:.4f}）⟹ χ² = {chi2a:.2f}、χ²/dof = {chi2a/(len(Y)-1):.2f}")
    print(f"  z：{' '.join(f'{x:+.1f}' for x in ra)}")
    P3 = dict(a=pa[0], chi2=chi2a, chi2_dof=chi2a / (len(Y) - 1), z=[float(x) for x in ra])

    print()
    print("[P2] 敏感度：形狀自由 (a,b,S,P)，S,P ∈ (0,1] —— 4 參數、4 dof（只報，不作判詞）")
    best = None
    for Si in [0.25, 0.5, 0.75, 1.0]:
        for Pi in [0.25, 0.5, 0.75, 1.0]:
            phiS = [((q2(CONTENT[r["ch"]][0])[0] - q2(CONTENT[r["ch"]][1])[0]) * Si
                     + (q2(CONTENT[r["ch"]][0])[1] - q2(CONTENT[r["ch"]][1])[1]) * Pi) for r in rows]
            A2 = np.array([phiS, PSI]).T
            p2, c2, r2 = wfit(A2, Y, S)
            if best is None or c2 < best[0]:
                best = (c2, Si, Pi, p2, r2)
    c2, Si, Pi, p2, r2 = best
    print(f"  最佳 S = {Si:.2f}、P = {Pi:.2f} ⟹ a = {p2[0]:.4f}、b = {p2[1]:.4f}（⟨O⟩_eff = {p2[1]/X:.3f}）")
    print(f"  χ² = {c2:.2f}、dof = 4、χ²/dof = {c2/4:.2f}、max|z| = {max(abs(r2)):.2f}")
    P2 = dict(S=Si, P=Pi, a=float(p2[0]), b=float(p2[1]), O_eff=float(p2[1] / X),
              chi2=c2, chi2_dof=c2 / 4, max_z=float(max(abs(r2))))

    _out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "09_route_c_total_dm.json")
    with open(_out, "w") as f:
        json.dump(dict(rows=rows, P1=P1, P2=P2, P3=P3, X=X, K0=K0), f, ensure_ascii=False, indent=2)
    print("\n-> data/09_route_c_total_dm.json")


main()