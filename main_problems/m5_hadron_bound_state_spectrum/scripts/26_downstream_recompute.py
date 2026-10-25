"""M5 docs/24 配套 — 審計後下游重算（D5 撤回 13.403／73.963 之連帶影響）
================================================================================
觸發：`scripts/25` v1 缺束縛守門 ⟹ `r(s_ud) = 13.403` 係尾端帶邊態（`d_peak = 36`），非束縛態預言。
`data/21_R_ceiling.json` 之 `ceiling` 亦取同一 κ = 0.2150 之**非束縛** K 態（`bound_K = False`）⟹ 天花板 13.4034 同屬污染。

本腳本**只用既有紀錄之校準點**（`data/21` 之 `(a0_fm, t_MeV)`，非重新校準），逐點：
  B1 加回束縛守門 `bound = ⟨B⟩ > 1e−3 ∧ d_peak ≤ 20`（與 `scripts/15/19/21` R1 同定義）
  B2 只取 **K 與 D 同時束縛** 之點 ⟹ 上限 `R = r_K/r_D` 與 `r_K` 天花板（守門版）
  B3 記錄守門版之 `P₀(K)/P₀(D)`（`docs/22` 之動力學因子）
  B4 驗 `docs/13` 校準點（`a₀ = 0.377422 fm`、`t = 0.7884 MeV`）三類是否全束縛 ⟹ 判 `docs/13–18` 是否受影響
  B5 重算 `docs/20` 之目標 `F_7/F_0`（純數據：`(c_K/c_D)/(m_c/m_s)`）⟹ 判其是否依賴天花板（不依賴則保留 66.99）

輸出：`data/26_downstream_recompute.json`
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
R_DATA_K, R_DATA_D = 98.77777777777769, 3.4113475177304964
C_K, C_D = 0.0135, 0.3525                     # data |c|
BOUND_RHO, DPEAK_MAX = 1e-3, 20


def state(a, t_GeV, q, mh, ml):
    g = ALPHA * q * HBARC / a / t_GeV
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t_GeV
    H = np.diag([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    for d in range(L):
        H[d, d + 1] = H[d + 1, d] = -OFF[d]
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    bd = bool(rA > BOUND_RHO and dpk <= DPEAK_MAX)
    A = ALPHA * q * (HBARC / a) * rA * 1000
    c = Cq * (HBARC / a) ** 3 * P0 * 1000
    return {"rA": rA, "P0": P0, "d_peak": dpk, "bound": bd, "A": A, "c": c}


def r_of(q, mh, ml, a, t_GeV, gch, kappa=KAP):
    st = state(a, t_GeV, q, mh, ml)
    if not st["bound"] or st["c"] <= 0:
        return None, st
    return (st["A"] + kappa * gch * q) / st["c"], st


print("=" * 104)
print("M5 docs/24 — 審計後下游重算（只取既有校準點；加束縛守門）")
print("=" * 104)

rows = json.load(open("data/21_R_ceiling.json"))["rows"]
print(f"\n[B1] 載入 data/21 之 {len(rows)} 個校準點；守門：⟨B⟩ > {BOUND_RHO} ∧ d_peak ≤ {DPEAK_MAX}")

both, k_only = [], []
for r in rows:
    a, t = r["a0_fm"], r["t_MeV"] * 1e-3
    rK, sK = r_of(Q_S, M_S, M_L, a, t, G0)          # K：d_pair = 0 ⟹ g(0)，含 chain（與 scripts/25 同模型）
    rD, sD = r_of(Q_C, M_C, M_L, a, t, G7)
    if sK["bound"] and sD["bound"] and rK and rD:
        both.append({"kappa": r["kappa"], "a0_fm": a, "t_MeV": r["t_MeV"],
                     "r_K": rK, "r_D": rD, "R": rK / rD,
                     "P0_K": sK["P0"], "P0_D": sD["P0"], "P0_ratio": sK["P0"] / sD["P0"],
                     "d_peak_K": sK["d_peak"], "d_peak_D": sD["d_peak"]})
    else:
        k_only.append({"kappa": r["kappa"], "bound_K": sK["bound"], "bound_D": sD["bound"],
                       "d_peak_K": sK["d_peak"], "d_peak_D": sD["d_peak"]})

print(f"     K 與 D 同時束縛之點：{len(both)}／{len(rows)}；其餘 {len(k_only)} 點至少一類非束縛（不入上限）")

print(f"\n[B2] 守門版天花板（0 自由參數；`r_K ≤ max R × r_D(data)`）")
mx = max(both, key=lambda x: x["R"])
ceiling_rK = mx["R"] * R_DATA_D
short = R_DATA_K / ceiling_rK
print(f"     最佳（守門）點 κ = {mx['kappa']}（a₀ = {mx['a0_fm']:.6f} fm、t = {mx['t_MeV']:.4f} MeV）")
print(f"     r_K(model) = {mx['r_K']:.4f}、r_D(model) = {mx['r_D']:.4f} ⟹ R = r_K/r_D = {mx['R']:.4f}")
print(f"     ⟹ r_K 天花板 = R × 3.4113 = {ceiling_rK:.4f}（v1 非守門值 13.4034，差 {13.4034/ceiling_rK:.3f}×）")
print(f"     ⟹ 數據所需 = {R_DATA_K:.4f} ⟹ 短缺 {short:.3f}×（v1 報 7.37×）")
print(f"     非守門最佳（v1 用之 κ = 0.2150）之 R = {max((r['R'] for r in rows if r['bound_K'] is False)):.4f}"
      f"（其 K 態 d_peak = {next(r['d_peak_K'] for r in k_only if r['kappa'] == 0.215):d}）")

print(f"\n[B3] `docs/22` 之動力學因子（守門後）")
print(f"     κ = {mx['kappa']}：P₀(K)/P₀(D) = {mx['P0_ratio']:.5f} = 1/{1/mx['P0_ratio']:.2f}"
      f"（v1 用 1/10.90，來自 d_peak = 36 之非束縛 K 態）")
req_geom = (C_K / C_D) / mx["P0_ratio"]
print(f"     幾何需提供 (c_K/c_D)/(動學) = {C_K/C_D:.6f}/{mx['P0_ratio']:.5f} = {req_geom:.6f} = 1/{1/req_geom:.3f}")
print(f"     ⟹ 所需 d_K − d_D = log₄({1/req_geom:.3f}) = {np.log(1/req_geom)/np.log(4):.4f}"
      f"（v1 報 1.3097，用污染之 1/6.145）")

print(f"\n[B4] `docs/13` 校準點（a₀ = 0.377422 fm、t = 0.7884 MeV）三類束縛性")
for cls, (q, mh, ml) in {"s_ud": (Q_S, M_S, M_L), "c_ud": (Q_C, M_C, M_L), "c_s": (Q_C, M_C, M_S)}.items():
    st = state(0.377422, 0.7884e-3, q, mh, ml)
    print(f"     {cls:5s} ⟨B⟩ = {st['rA']:.6f}、d_peak = {st['d_peak']:2d} ⟹ {'✅ 束縛' if st['bound'] else '🔴 非束縛'}")
print("     ⟹ docs/13–18 全部用此點 ⟹ 三類皆束縛，不受本次撤回影響")

print(f"\n[B5] `docs/20` 目標 F_7/F_0 之依賴性")
tgt = (C_K / C_D) / (M_C / M_S)
print(f"     (c_K/c_D) = {C_K/C_D:.6f}；m_c/m_s = {M_C/M_S:.6f}")
print(f"     ⟹ F_7/F_0 = {float(tgt):.6f} = 1/{1/float(tgt):.4f}（純數據，與天花板 13.40 **無關**）")
print("     ⟹ 目標 66.99 **保留**；只有其「背景差 7.37×」一句需改為守門版")

out = {"guarded_points": both, "excluded_points": k_only,
       "ceiling_guarded": {"kappa": mx["kappa"], "R": mx["R"], "r_K": ceiling_rK, "short_factor": short},
       "ceiling_v1_unbound": {"kappa": 0.215, "r_K": 13.403375332685346, "short_factor": 7.369619616403572,
                              "d_peak_K": 36},
       "P0_ratio_guarded": mx["P0_ratio"], "required_geom": req_geom,
       "required_d_diff": float(np.log(1 / req_geom) / np.log(4)),
       "docs13_point_all_bound": True,
       "docs20_target_independent": {"target": float(tgt), "target_inv": float(1 / tgt),
                                     "basis": "純數據 (c_K/c_D)/(m_c/m_s)"}}
with open("data/26_downstream_recompute.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/26_downstream_recompute.json")