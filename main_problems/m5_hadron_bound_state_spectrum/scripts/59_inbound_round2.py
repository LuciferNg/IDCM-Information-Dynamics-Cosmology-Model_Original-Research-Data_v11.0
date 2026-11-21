"""M5 docs/57 配套 — 第二輪入站裁定：M10 符號約定自糾、W1 σ 口徑、M18 P1–P3 核
================================================================================
觸發（register 未回覆入站）：
  M10（10-29 01:36／01:47／11-04 00:05）：指出 M5 `docs/54` 之「自能反號」係**標籤約定錯** ——
      M10 之 `em ≡ X − (m_n − m_p) = em(p) − em(n)`（電磁使**質子較重**），而 M5 寫成 `em(n) − em(p)`；
      並指需求 `1.2667 ± 0.094` ⟹ `A_max(v2) = 1.2494` 只差 **0.184σ** ⟹ **唔構成排除**。
  M18（11-21 01:53／15:20／23:17／23:42）：要求 M5 以既有 `|ψ_d|²` 核 P1–P3（`J = 1,2` 於 `d = 0` 振幅須 = 0）
      並提供根權重 `cap(J)` 對照（`J = 0` 100%；`J = 1,2,5` 0%；`J = 3/4/6/7/8` 57.14/44.44/38.46/26.67/23.53%）。

判準（事前寫死）：
  A1 態內電荷平方和 `Σq²` 逐態精確算（`π⁰ = (uū − dd̄)/√2` ⟹ `5/9`；`π⁺ = ud̄` ⟹ `5/9`）。
  A2 自能預言：`Δ0^SE ≡ κ[Σq²(中性) − Σq²(帶電)]`，`κ = (α/2)(ħc/a₀)g(0)`（v2），`κ > 0`。
     逐族判號：`Δ0^SE` 與實測 `Δ0` 同號 ⟹ 號對；異號 ⟹ 號錯；恰零 ⟹ 該族無貢獻。
  A3 重子：M10 約定 `em(p) − em(n) = κ[Σq²(p) − Σq²(n)]`，需求 `+1.2667`；判號與量級倍數。
  A4 W1 σ 口徑：需求 `1.2667 ± 0.094` ⟹ 缺口以 σ 計；`< 1σ` ⟹ **不得**稱排除。
  B1 M18 核法：取 M5 模型最低 6 個本徵態之 `P_d = |v_d|²`（v2 約定），報 `P₀`（`d = 0` 權重）與 `d_peak`。
  B2 判：`cap(J)` 係**上限** ⟹ 只有 `cap = 0` 之 `J = 1,2,5` 可證偽；若所有束縛態 `P₀ > 0` ⟹
       **M5 譜中不存在 `J = 1,2,5` 之載體** ⟹ P1 無反例（排除式驗證），但 M5 亦不提供該載體。
  B3 不引入新參數、不改 M18 檔案；本檔只交數字與判定。
輸出：`data/59_inbound_round2.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

ALPHA = 1 / 137.035999084
HBARC = 0.1973269804
A0, T0 = 0.384184, 1.7711e-3       # v2 正典：a₀ 用 fm、t 用 **GeV**（記錄值 1.7711 MeV）
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
G0 = 2581 / 352
KAPPA = 0.5 * ALPHA * (HBARC / A0) * G0 * 1000     # MeV：(α/2)(ħc/a₀)g(0)
M_C, M_S, M_L, M_B = 1.27, 0.495, 0.30, 4.18
Q = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3, "b": -1 / 3}

print("=" * 104)
print("M5 docs/57 — 第二輪入站裁定（M10 符號自糾 ／ W1 σ 口徑 ／ M18 P1–P3 核）")
print("=" * 104)

print("\n[A1] 態內 Σq²（含反夸克；π⁰ = (uū − dd̄)/√2）")
states = {
    "π⁰": (5 / 9, "½(q_u²+q_ū²+q_d²+q_d̄²) = 5/9"),
    "π⁺": (Q["u"] ** 2 + Q["d"] ** 2, "ud̄"),
    "K⁰": (Q["d"] ** 2 + Q["s"] ** 2, "d s̄"),
    "K⁺": (Q["u"] ** 2 + Q["s"] ** 2, "u s̄"),
    "D⁰": (Q["c"] ** 2 + Q["u"] ** 2, "c ū"),
    "D⁺": (Q["c"] ** 2 + Q["d"] ** 2, "c d̄"),
    "B⁰": (Q["d"] ** 2 + Q["b"] ** 2, "d b̄"),
    "B⁺": (Q["u"] ** 2 + Q["b"] ** 2, "u b̄"),
    "p": (2 * Q["u"] ** 2 + Q["d"] ** 2, "uud"),
    "n": (Q["u"] ** 2 + 2 * Q["d"] ** 2, "udd"),
}
for k, (v, c) in states.items():
    print(f"  Σq²({k}) = {v:.6f}  （{c}）")

print(f"\n[A2/A3] 自能預言（κ = (α/2)(ħc/a₀)g(0) = {KAPPA:.4f} MeV，v2）")
# Δ0 ≡ m(中性) − m(帶電)（PDG：π⁰ 134.9768/π± 139.5704；K⁰ 497.611/K± 493.677；
#      D⁰ 1864.84/D± 1869.66；B⁰ 5279.66/B± 5279.34）
DATA = {"π": (states["π⁰"][0], states["π⁺"][0], -4.5936, "π⁰ 較**輕**"),
        "K": (states["K⁰"][0], states["K⁺"][0], +3.9340, "K⁰ 較**重**"),
        "D": (states["D⁰"][0], states["D⁺"][0], -4.8200, "D⁰ 較**輕**"),
        "B": (states["B⁰"][0], states["B⁺"][0], +0.3100, "B⁰ 較**重**")}
rows = []
for fam, (sqN, sqC, obs_dN_minus_dC, note) in DATA.items():
    d_sq = sqN - sqC
    pred = KAPPA * d_sq
    if abs(d_sq) < 1e-12:
        verdict = "恰零（該族不可能由自能產生分裂）"
    elif (pred > 0) == (obs_dN_minus_dC > 0):
        verdict = f"號對；量級 {abs(pred):.4f} vs 觀測 {abs(obs_dN_minus_dC):.4f}"
    else:
        verdict = "🔴 號錯"
    rows.append({"fam": fam, "dq2": d_sq, "pred": pred, "obs": obs_dN_minus_dC, "verdict": verdict})
    print(f"  {fam}: Δ(Σq²) = 中性−帶電 = {d_sq:+.6f} ⟹ Δ0^SE = {pred:+.4f} MeV；觀測 {obs_dN_minus_dC:+.4f}"
          f"（{note}）⟹ {verdict}")
sqp, sqn = states["p"][0], states["n"][0]
pred_bary = KAPPA * (sqp - sqn)
print(f"  重子（M10 約定 em(p)−em(n)）：κ[Σq²(p)−Σq²(n)] = {pred_bary:+.4f} MeV vs 需求 +1.2667 ⟹ "
      f"號{'對' if pred_bary > 0 else '錯'}、超 {abs(pred_bary/1.2667):.2f}×")

print("\n[A4] W1 σ 口徑（M10 docs/19）")
NEED, SIG, AMAX_V2, AMAX_V1 = 1.2667, 0.094, 1.249371, 1.271754
gap = NEED - AMAX_V2
print(f"  需求 {NEED} ± {SIG}；A_max(v2) = {AMAX_V2:.6f} ⟹ 缺口 {gap:.4f} MeV = **{gap/SIG:.3f}σ**")
print(f"  A_max(v1) = {AMAX_V1:.6f} ⟹ 在上限內（餘裕 {100*(AMAX_V1/NEED-1):+.2f}%）")
print(f"  ⟹ [A4] 缺口 < 1σ ⟹ **不得稱排除、不得稱『超界 1.37%』**；W1 只可記"
      f"『結構封（IDCM 內無 em 算子）＋ 約定閘 ＋ 自能同號但 3.62× 數值否證』")

print("\n[B1] M5 模型最低 6 態之 P₀（d = 0 權重；v2 約定）")


def spectrum(mh, ml, q, aC=A0, aK=A0, t=T0, nstates=6):
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    g = ALPHA * q * HBARC / aC / t
    h = Cq * (HBARC / aK) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    ev, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, nstates - 1))
    out = []
    for i in range(nstates):
        v = evec[:, i] if evec.ndim == 2 else evec
        w = v ** 2
        w = w / w.sum()
        out.append({"P0": float(w[0]), "dpk": int(np.argmax(w)),
                    "rA": float((w * np.array(POW)).sum()), "E": float(ev[i])})
    return out


FAMS = [("π (0.30,0.30)", M_L, M_L, Q["u"]), ("K (0.495,0.30)", M_S, M_L, Q["s"]),
        ("D (1.27,0.30)", M_C, M_L, Q["c"]), ("B (4.18,0.30)", M_B, M_L, Q["b"])]
CAP = {0: 1.00, 1: 0.00, 2: 0.00, 5: 0.00, 3: 0.5714, 4: 0.4444, 6: 0.3846, 7: 0.2667, 8: 0.2353}
print(f"  cap(J)（M18 docs/23）：{CAP}")
allP0 = []
for name, mh, ml, q in FAMS:
    sp = spectrum(mh, ml, q)
    print(f"  {name}:")
    for i, s in enumerate(sp):
        allP0.append(s["P0"])
        print(f"    態{i}: P₀ = {s['P0']:.6f}  d_peak = {s['dpk']:>2}  ⟨B⟩ = {s['rA']:.4f}  E = {s['E']:+.6f}")
print(f"\n[B2] 判：全部 {len(allP0)} 態之 P₀ > 0 = {all(p > 0 for p in allP0)}；最小 P₀ = {min(allP0):.6e}")
print("  ⟹ cap(J=1,2,5) = 0 ⟹ 任何 P₀ > 0 之態**不可能**係 J = 1,2,5")
print("  ⟹ M5 譜中無 J=1,2,5 之載體 ⟹ P1（J=1,2 於 d=0 振幅 = 0）由 M5 側**排除式驗證**（無反例），")
print("     但 M5 亦唔提供 J=1,2 載體 ⟹ M18 之向量介子需求須另立載體（非本 M 譜）")
print("  其餘 cap（J=3/4/6/7/8 之 23.5–57.1%）皆係上限 ⟹ 本譜全部 P₀ ≤ 0.39 自動滿足 ⟹ 無鑑別力")

out = {"kappa": KAPPA, "meson_rows": rows, "baryon_pred": pred_bary,
       "w1_gap_sigma": gap / SIG, "A_max_v1": AMAX_V1, "A_max_v2": AMAX_V2,
       "P0_min": min(allP0), "all_P0_positive": all(p > 0 for p in allP0), "cap": CAP}
with open("data/59_inbound_round2.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/59_inbound_round2.json")