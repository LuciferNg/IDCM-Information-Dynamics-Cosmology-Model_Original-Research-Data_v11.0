"""M5 docs/19 配套 — D3：形狀比 A/|c| 之雙系統閉合檢定（唔加任何新參數）
================================================================================================
背景：`docs/14` §四.1 已證 **可測恒等式**「A/|c| 與電荷無關，只由形狀 (a, t, m_h, m_l) 決定」。
數據給出：
    D：A/|c| = 1.2025/0.3525 = 3.41135
    K：A/|c| = 1.3335/0.0135 = 98.7778            （A = (E1+3E2)/4、|c| = (E1−E2)/4）
⟹ 同一 (a, t) 要造出 **29.0×** 之形狀比躍遷，而唯一變化係 m_h（1.27 → 0.495）。

三項檢定（判準事前寫死）：
  T1 盲測（既有校準）：(a_C, t) 由 D 兩數定 ⟹ 預言 K 之 A/|c|；命中 = 與 98.7778 差 ≤ 2σ（σ 由
     E1／E2 之 σ 傳播：σ(A/|c|) 由此式數值傳播）。
  T2 自由度上界：2 參數 (a, t) 對 4 數據（兩系 A、|c|）。若**連只解兩條比值方程**都無解 ⟹
     形狀族結構上不可能；若有解，逐位報 4 個絕對值之 z。
  T3 帶 κ 版本（chain 扇區，κ = 0.26986）同場重列，判 chain 是否能救形狀比。

禁：不得為命中而改族／改 m_h／改 Σσ；未命中記 NEG，不得記 0。

輸出：`data/20_shape_ratio_closure.json`
"""
import json

import numpy as np
from scipy.optimize import root

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L = 1.27, 0.495, 0.30
Q_C, Q_S = 2 / 3, 1 / 3
KAPPA, G0 = 0.26986, 2581 / 352
S1D = S2D = 0.0940                 # D 通道 σ（docs/15：四通道一致 0.0940）
S1K, S2K = 0.0960728889958036, 0.34122719704032967
E1D, E2D = 2.2600, 0.8500
E1K, E2K = 1.3740, 1.3200
A_D, C_D = 1.2025, 0.3525
R_D, R_K = A_D / C_D, 1.3335 / 0.0135


def shape(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    return float((w * np.array(POW)).sum()), float(w[0]), int(np.argmax(w)), bool(w[0] > 1e-3 and int(np.argmax(w)) <= 20)


def amps(a, t, q, mh, ml):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    rA, P0, dp, bd = shape(g, h)
    return ALPHA * q * (HBARC / a) * rA * 1000, Cq * (HBARC / a) ** 3 * P0 * 1000, bd, P0


def ratio(a, t, q, mh, ml):
    A, c, bd, P0 = amps(a, t, q, mh, ml)
    return (A / c) if c > 1e-14 else np.inf


print("=" * 104)
print("M5 docs/19 — D3：A/|c| 形狀比雙系統閉合檢定")
print("=" * 104)
print(f"\n[D0] 數據：r_D = {A_D}/{C_D} = {R_D:.5f}；r_K = 1.3335/0.0135 = {R_K:.5f}"
      f"；比 = {R_K/R_D:.3f}×（唯一變化：m_h 1.27 → 0.495）")

print("\n[T1] 盲測（(a_C, t) 由 D 兩數定；K 全盲）")
print(f"     {'κ':>8s}{'a₀ (fm)':>10s}{'t (MeV)':>10s}{'A_K':>9s}{'|c|_K':>8s}{'r_K 模型':>11s}{'r_K 數據':>10s}{'差':>8s}")


def calib(kappa):
    """由 D 兩數解 (a, t)"""
    from scipy.optimize import root as _root
    chd = kappa * (1447 / 352) * Q_C

    def F(x):
        a, t = np.exp(x[0]), np.exp(x[1])
        A, c, bd, _ = amps(a, t, Q_C, M_C, M_L)
        return [np.log(max(A + chd, 1e-12) / A_D), np.log(max(c, 1e-12) / C_D)]

    best = None
    for x0a in np.linspace(-2.0, 1.0, 10):
        for x0b in np.linspace(-9.0, 1.0, 12):
            s = _root(F, [x0a, x0b], method="hybr")
            if not s.success:
                continue
            a, t = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
            A, c, bd, _ = amps(a, t, Q_C, M_C, M_L)
            res = abs(A + chd - A_D) + abs(c - C_D)
            if res < 1e-9 and bd and (best is None or res < best[0]):
                best = (res, a, t)
    return None if best is None else (best[1], best[2])


t1 = []
for kp in (0.0, KAPPA):
    sol = calib(kp)
    if sol is None:
        print(f"     {kp:8.5f}{'—':>10s}{'—':>10s}  D 校準無解 ⟹ N/A")
        t1.append({"kappa": kp, "feasible": False})
        continue
    a, t = sol
    A_Kt, c_Kt, bd, _ = amps(a, t, Q_S, M_S, M_L)
    chk = kp * G0 * Q_S
    r_mod = (A_Kt + chk) / c_Kt if c_Kt > 0 else np.inf
    print(f"     {kp:8.5f}{a:10.6f}{t*1000:10.4f}{A_Kt+chk:9.4f}{c_Kt:8.4f}{r_mod:11.4f}{R_K:10.4f}"
          f"{r_mod/R_K:8.2e}")
    t1.append({"kappa": kp, "feasible": True, "a0_fm": a, "t_MeV": t * 1000, "A_K": A_Kt + chk,
               "c_K": c_Kt, "r_K_model": r_mod, "r_K_data": R_K, "ratio": r_mod / R_K})

print("\n[T2] 自由度上界：只用兩條比值方程解 (a, t)（2 方程 / 2 未知），再看 4 個絕對值")
print("     （若無解 ⟹ 形狀族結構上造不出呢對比值；若有解 ⟹ 檢絕對值是否同時中）")
sols = []
for kp in (0.0, KAPPA):
    chk = kp * G0 * Q_S
    chd = kp * (1447 / 352) * Q_C

    def F(x):
        a, t = np.exp(x[0]), np.exp(x[1])
        A_d, c_d, _, _ = amps(a, t, Q_C, M_C, M_L)
        A_k, c_k, _, _ = amps(a, t, Q_S, M_S, M_L)
        if c_d <= 0 or c_k <= 0:
            return [1e3, 1e3]
        return [np.log((A_d + chd) / (c_d * R_D)), np.log((A_k + chk) / (c_k * R_K))]

    found = []
    for x0a in np.linspace(-2.0, 1.5, 12):
        for x0b in np.linspace(-9.0, 1.0, 14):
            s = root(F, [x0a, x0b], method="hybr")
            if not s.success or np.max(np.abs(F(s.x))) > 1e-9:
                continue
            a, t = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
            A_d, c_d, bd, _ = amps(a, t, Q_C, M_C, M_L)
            A_k, c_k, _, _ = amps(a, t, Q_S, M_S, M_L)
            rec = {"kappa": kp, "a0_fm": a, "t_MeV": t * 1000, "bound_D": bd,
                   "A_D": A_d + chd, "c_D": c_d, "A_K": A_k + chk, "c_K": c_k,
                   "E1_D": A_d + chd + 3 * c_d, "E2_D": A_d + chd - c_d,
                   "E1_K": A_k + chk + 3 * c_k, "E2_K": A_k + chk - c_k}
            rec["z"] = [(rec["E1_D"] - E1D) / S1D, (rec["E2_D"] - E2D) / S2D,
                        (rec["E1_K"] - E1K) / S1K, (rec["E2_K"] - E2K) / S2K]
            rec["pass"] = bool(bd and all(abs(z) <= 2 for z in rec["z"]))
            if not any(abs(r["a0_fm"] - a) < 1e-6 and abs(r["t_MeV"] - rec["t_MeV"]) < 1e-6 for r in found):
                found.append(rec)
    print(f"     κ = {kp:8.5f}：解得 {len(found)} 個 (a, t) 根")
    for r in sorted(found, key=lambda r: np.sum(np.square(r["z"])))[:4]:
        print(f"        a₀ = {r['a0_fm']:9.6f} fm、t = {r['t_MeV']:9.4f} MeV ⟹ "
              f"E1_D = {r['E1_D']:7.4f}（z {r['z'][0]:+6.2f}）、E2_D = {r['E2_D']:7.4f}（z {r['z'][1]:+6.2f}）、"
              f"E1_K = {r['E1_K']:7.4f}（z {r['z'][2]:+6.2f}）、E2_K = {r['E2_K']:7.4f}（z {r['z'][3]:+6.2f}）"
              f" ⟹ {'✅ 4/4 命中' if r['pass'] else '🔴 絕對值不中'}")
        r["z"] = [float(z) for z in r["z"]]
    sols += found
    if not found:
        print("        ⟹ 🔴 兩條比值方程對 (a, t) **無解** ⟹ 形狀族結構上不可能造出呢對比值")

print("\n[T3] 形狀能力掃描：r(a, t) 對 m_h 之可達範圍（固定 q、m_l；網格 60×60）")
for mh, tag in ((M_C, "D（m_h = 1.27）"), (M_S, "K（m_h = 0.495）")):
    vals = []
    for la in np.linspace(np.log(0.05), np.log(2.0), 60):
        for lt in np.linspace(np.log(1e-6), np.log(5.0), 60):
            A, c, bd, _ = amps(float(np.exp(la)), float(np.exp(lt)), Q_S if mh == M_S else Q_C, mh, M_L)
            if bd and c > 1e-12:
                vals.append(A / c)
    v = np.array(vals)
    print(f"     {tag}：r ∈ [{v.min():.3f}, {v.max():.3f}]（中位數 {np.median(v):.3f}）")

out = {"r_D": R_D, "r_K": R_K, "t1_blind": t1, "t2_solutions": sols,
       "t2_any_pass": bool(any(s["pass"] for s in sols))}
with open("data/20_shape_ratio_closure.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/20_shape_ratio_closure.json")