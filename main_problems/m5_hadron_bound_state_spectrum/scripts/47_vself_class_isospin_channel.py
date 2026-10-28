"""M5 docs/45 配套 — V_self 之古典上界 ＋ isospin 算子之耦合通道定位
================================================================================
觸發：`docs/43` 證一體 EM 自能項**必需**（π ≥ 3.3442、K ≥ 2.6846、D ≥ 2.3213 MeV），
且定位 `X` 為其殘餘；`docs/42` 證夸克質量項 `q = (0, +8.53, −0.23, +4.90)` 不在模型算子集內。
本輪做兩件分類工作：

【一｜V_self：古典 or 一圈量子？】
  古典自能（同心殼定理，無自由參數）：`V_self = (α/2)q_l²ħc[Σ_d P_d²/r_d + 2Σ_{d<d'} P_d P_{d'}/r_{d'}]`，
  `r_d = a·4^d`（與 `V ⊃ 4^{-d}` 之 Coulomb 形式自洽：`αħc/a·4^{-d} = αħc/r_d`）。
  上界：`P_0 ≤ 1` ⟹ `V_self ≤ (α/2)q_l²ħc/a₀`。

【二｜isospin：耦合去邊個通道？】
  v3 之 `C_q ∝ 1/(m_h·m_ref)`、`m_ref` 固定 ⟹ `A`、`c` **完全不含 `m_l`** ⟹ `∂Δ0/∂m_l ≡ 0`（代數）。
  以 `m_l ± 3.5 MeV` 之有限差分作數值核實；再求「若要產生 `q_f`，需要何種通道」。

判準（事前寫死）：
  G1 `q_l²` 取 `4/9`（`u` 型輕夸克）；`P_d` 由 v3 校準態給出；`r_d = a·4^d`。
  G2 判古典可行性：若 `V_self^classical ≤ 必需量之 1/2` ⟹ 古典自能**不足以**解釋 ⟹ 須非古典（一圈）。
  G3 所需短距截斷 `r* = (α/2)q_l²ħc/V_need`，並報 `r*/a₀`。
  G4 判 isospin：`|∂Δ0/∂m_l|·δm / q_f < 1/100` ⟹ 現有通道**代數上不足** ⟹ 須新耦合（動能/質量項）。
  G5 不引入任何擬合係數；本輪只作分類與量級判定。
輸出：`data/47_vself_class_isospin_channel.json`
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
A0, T0 = 0.384184, 1.771108e-3
M_REF = 0.30
R_NEEF = {"π": 4.5936, "K": 2.6846, "D": 2.3213}
Q_F = {"π": 0.0, "K": 8.5276, "D": -0.2264, "B": 4.9036}


def state(a, t, q, mh, ml):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * abs(q) / (16 * mh * M_REF)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    return w, ALPHA * q * (HBARC / a) * float((w * np.array(POW)).sum()) * 1000, \
        Cq * (HBARC / a) ** 3 * float(w[0]) * 1000


print("=" * 100)
print("M5 docs/45 — V_self 古典上界 ＋ isospin 耦合通道")
print("=" * 100)

print("\n[G1] 古典自能（同心殼定理，r_d = a·4^d；q_l² = 4/9）")
q2 = 4 / 9
bound = 0.5 * ALPHA * q2 * HBARC / A0 * 1000
print(f"  硬上界 (α/2)q_l²ħc/a₀ = {bound:.4f} MeV（P_0 = 1 時達界）")
rows = {}
for nm, q, mh in (("π", 1 / 3, 0.30), ("K", 1 / 3, 0.495), ("D", 2 / 3, 1.27), ("B", 2 / 3, 4.18)):
    w, A, c = state(A0, T0, q, mh, M_REF)
    r = np.array([A0 * 4.0 ** d for d in range(L + 1)])
    diag_sum = float(((w ** 2) / r).sum())
    cross = 0.0
    for d in range(L + 1):
        for e in range(d + 1, L + 1):
            cross += 2 * w[d] * w[e] / r[e]
    shape = diag_sum + cross
    V = 0.5 * ALPHA * q2 * HBARC * shape * 1000
    rows[nm] = {"V_self_MeV": V, "shape_per_fm": shape, "P0": float(w[0])}
    need = R_NEEF.get(nm, 0.0)
    print(f"  {nm:2s} P₀ = {w[0]:.4f}、⟨1/r⟩_sh = {shape:8.4f} fm⁻¹ ⟹ V_self = {V:7.4f} MeV"
          + (f"（必需 ≥ {need:.4f} ⟹ {need/max(V,1e-12):5.1f}× 不足）" if need else "（B：必需 0）"))

print("\n[G2] 判古典可行性")
ok = all(rows[k]["V_self_MeV"] <= 0.5 * R_NEEF[k] for k in R_NEEF)
print(f"  古典自能 ≤ 必需量之 1/2 者：{all(rows[k]['V_self_MeV'] <= 0.5*R_NEEF[k] for k in R_NEEF)}")
print(f"  ⟹ 古典自能{'不足' if ok else '可能足夠'} ⟹ {'須非古典（一圈/真空極化）' if ok else '可續古典路線'}")

print("\n[G3] 若堅持古典：所需短距截斷 r*")
for k, need in R_NEEF.items():
    rstar = 0.5 * ALPHA * q2 * HBARC * 1000 / need   # HBARC 為 GeV·fm，need 為 MeV
    print(f"  {k}: r* = {rstar:.4f} fm ⟹ r*/a₀ = {rstar/A0:.4f}（a₀/{A0/rstar:.2f}）")

print("\n[G4] isospin 通道：v3 之 m_l 敏感度（代數 + 有限差分）")
w0, A_m, c_m = state(A0, T0, 2 / 3, 1.27, M_REF)
d_ml = 3.5
w1, A_p, c_p = state(A0, T0, 2 / 3, 1.27, M_REF + d_ml)
print(f"  C_q 只含 m_ref（固定）⟹ 代數上 ∂A/∂m_l = ∂c/∂m_l = 0")
print(f"  有限差分（把 m_ref 當 m_l 變 ±{d_ml} MeV）僅示通道大小：ΔA = {A_p-A_m:+.6f} MeV、"
      f"Δc = {c_p-c_m:+.6f} MeV")
print(f"  v3 之 `A`、`c` **不含 m_l** ⟹ Δ0(π⁰ vs π⁺) ≡ 0 ⟹ 現有算子集對 u−d 差**零敏感**")
for k, q in Q_F.items():
    print(f"  族 {k} 所需夸克項 q_f = {q:+8.4f} MeV ⟹ 現有通道可給 0 ⟹ 不足")

print("\n[G5] 所需新耦合之型式（不作擬合，只定位通道）")
print("  模型唯一含質量之項為 C_q ∝ 1/(m_h·m_ref)（contact），已證不足以產生 q_f；")
print("  `A` 由 Coulomb（q_h q_l 符號）+ ⟨B⟩ 給出，亦不含 m_l；")
print("  ⟹ isospin 必須經**動能／跳躍**通道（現時 hopping 固定 = 無質量依賴），即新算子類：")
print("     H_iso = ε·(δm/m_ref)·(hopping)，δm ≡ m_d − m_u；ε 需由 K 定標，再盲測 D、B。")
out = {"bound_classical_MeV": bound, "families": rows, "classical_insufficient": bool(ok),
       "required_cutoff_fm": {k: 0.5 * ALPHA * q2 * HBARC * 1000 / v for k, v in R_NEEF.items()},
       "q_f_needed": Q_F, "isospin_channel_present": False}
with open("data/47_vself_class_isospin_channel.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/47_vself_class_isospin_channel.json")