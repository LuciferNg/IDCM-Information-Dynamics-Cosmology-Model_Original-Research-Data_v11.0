"""M5 docs/33 配套 — 第四前線：兩個長度（Coulomb 射程 a_C vs contact 尺度 a_X）
================================================================================
`docs/32` 定理：單尺度模型不可能產生 `A_K`（需 `a ≤ 0.35997 fm`，但 D 校準要 `a = 0.384184 fm`）。
必要延伸 = 分離兩個長度：Coulomb 用 `a_C`、contact 用 `a_X`，`r ≡ a_C/a_X`，定理給硬約束 `r < 0.9370`。

本輪問兩件事（分開報，唔混）：
  Q-A（**量測**，後驗）：對 `r ∈ [0.30, 1.00]` 逐點以 D 兩數校準 `(a_C, t)`，求 K 盲測 z 與 `R = r_K/r_D`，
       並定出「K 兩通道皆 |z| ≤ 2」所需之 `r*` 區間。此為**目標值**，唔係結果。
  Q-B（**對照**，明確標後驗）：把 `r*` 同結構候選（`B^−1/2, B^−1/4, (b−1)^−1/2, (b−1)^−1, …`）比對；
       若命中，只可記「疑似」，並須以 `docs/31` 之 MC 零假設 harness 判顯著性（本輪不判死，留待獨立輪）。

判準（事前寫死）：
  V1 校準一律用 D 兩數（A_D = 1.2025、|c|_D = 0.3525），K／Ξc 零新參數、零調整。
  V2 `r` 之掃描範圍寫死 [0.30, 1.00]，步長 0.005；不得因結果收窄或加點。
  V3 通過定義：`E(K⁰K⁻)` 與 `E(K*⁰K*⁺)` 皆 |z| ≤ 2。
  V4 `r*` 只作「數據要求」，唔可以宣稱結構；結構候選比對須另列並標「後驗」。
  V5 `R` 須同時報；若 `r*` 存在但 `R` 仍 ≪ 28.96，則記「A／c 分裂仍不足」。
輸出：`data/35_two_lengths.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
A_D, C_D = 1.2025, 0.3525
OBS = {"E1": 1.3740, "E2": 1.3200, "s": 0.3545}
R_DATA = 98.77777777777769 / (A_D / C_D)
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/33 — 第四前線：兩個長度 r = a_C/a_X（定理硬約束 r < 0.9370）")
print("=" * 104)
print(f"\n數據要求 R = r_K/r_D = {R_DATA:.4f}；觀測 E(K⁰K⁻) = 1.3740 ± 0.3545、E(K*⁰K*⁺) = 1.3200 ± 0.3545")


def amps(aC, aX, t, q, mh, ml):
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aX) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    return (ALPHA * q * (HBARC / aC) * rA * 1000,
            Cq * (HBARC / aX) ** 3 * P0 * 1000,
            {"rA": rA, "P0": P0, "dpk": dpk, "bound": bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)})


def bisect(f, lo, hi, n=200):
    flo = f(lo)
    if f(hi) * flo > 0:
        return None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-13:
            break
    return 0.5 * (lo + hi)


def calib(r, ngrid=360):
    """固定 r ⟹ 解 (a_C, t) 使 D 之 (A, |c|) = (A_D, C_D)。"""
    def t_solve(aC):
        x = bisect(lambda xx: amps(aC, aC / r, float(np.exp(xx)), Q_C, M_C, M_L)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        aC = float(np.exp(x))
        t = t_solve(aC)
        if t is None:
            prev = None
            continue
        k = amps(aC, aC / r, t, Q_C, M_C, M_L)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            rr = bisect(lambda xx: amps(float(np.exp(xx)), float(np.exp(xx)) / r,
                                        t_solve(float(np.exp(xx))), Q_C, M_C, M_L)[1] - C_D,
                        prev[0], x)
            if rr is not None:
                aa = float(np.exp(rr))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


rows = []
for r in np.arange(0.30, 1.0001, 0.005):
    sol = calib(float(r))
    if sol is None:
        rows.append({"r": float(r), "calibrated": False})
        continue
    aC, t = sol
    AD, cD, stD = amps(aC, aC / r, t, Q_C, M_C, M_L)
    AK, cK, stK = amps(aC, aC / r, t, Q_S, M_S, M_L)
    E1, E2 = AK + 3 * cK, AK - cK
    z1, z2 = (E1 - OBS["E1"]) / OBS["s"], (E2 - OBS["E2"]) / OBS["s"]
    rows.append({"r": float(r), "calibrated": True, "aC_fm": aC, "aX_fm": aC / r, "t_MeV": t * 1000,
                 "A_K": AK, "c_K": cK, "rK": AK / cK if cK else np.nan, "rD": AD / cD,
                 "R": (AK / cK) / (AD / cD) if cK else np.nan,
                 "z1": z1, "z2": z2, "pass": bool(abs(z1) <= 2 and abs(z2) <= 2),
                 "K": stK})
print(f"\n[Q-A] 掃描 {len(rows)} 點（r ∈ [0.30, 1.00]，步長 0.005，V2 寫死）")
hdr = f"  {'r':>6s}{'a_C(fm)':>10s}{'a_X(fm)':>10s}{'t(MeV)':>9s}{'A_K':>8s}{'c_K':>8s}{'r_K':>8s}{'R':>8s}{'z1':>7s}{'z2':>7s}  通"
print(hdr)
for row in rows:
    if not row["calibrated"] or not np.isfinite(row.get("R", np.nan)):
        continue
    mark = "✅" if row["pass"] else ""
    if abs(row["r"] * 100 % 5) < 1e-6 or row["pass"]:
        print(f"  {row['r']:6.3f}{row['aC_fm']:10.5f}{row['aX_fm']:10.5f}{row['t_MeV']:9.4f}"
              f"{row['A_K']:8.4f}{row['c_K']:8.5f}{row['rK']:8.3f}{row['R']:8.3f}"
              f"{row['z1']:+7.2f}{row['z2']:+7.2f}  {mark}")

good = [x for x in rows if x.get("calibrated") and x.get("pass")]
print("\n[Q-A 結果]")
if good:
    lo, hi = min(x["r"] for x in good), max(x["r"] for x in good)
    print(f"  通過區間 r* ∈ [{lo:.3f}, {hi:.3f}]（{len(good)} 點）")
    b = min(good, key=lambda x: abs(x["z1"]) + abs(x["z2"]))
    print(f"  最佳點 r = {b['r']:.3f}：A_K = {b['A_K']:.4f}、c_K = {b['c_K']:.5f}、r_K = {b['rK']:.3f}"
          f"（需 98.778）、R = {b['R']:.3f}（需 {R_DATA:.2f}）⟹ 仍缺 {R_DATA/b['R']:.1f}×")
else:
    lo = hi = None
    print("  🔴 無 r 令 K 兩通道同時 |z| ≤ 2")
    b = min((x for x in rows if x.get("calibrated")),
            key=lambda x: max(abs(x["z1"]), abs(x["z2"])))
    print(f"  最接近 r = {b['r']:.3f}：z1 = {b['z1']:+.2f}、z2 = {b['z2']:+.2f}、R = {b['R']:.3f}")

print("\n[Q-B] 結構候選比對（**後驗**，只記疑似，不判死）")
cands = {"B^−1/2 = 1/2": 0.5, "B^−1/4 = 1/√2": 2 ** -0.5, "(b−1)^−1/2 = 1/√3": 3 ** -0.5,
         "(b−1)^−1 = 1/3": 1 / 3, "3/4": 0.75, "d_H = ln3/ln4": np.log(3) / np.log(4)}
tgt = good[0]["r"] if good else (b["r"] if b else None)
for k, v in cands.items():
    print(f"  {k:22s} = {v:.6f}   偏差相對 r* = {tgt:.3f}：{(v/tgt-1)*100:+.2f}%" if tgt else k)

print("\n[判決]")
if good and b["R"] >= R_DATA / 3:
    print("  🟢 兩長度模型可同時通過 K 盲測並顯著收窄 ⟹ 須全鏈重算")
elif good:
    print("  🟡 K 盲測可通過，但 r_K 仍不足（A／c 分裂未夠）⟹ 兩長度必要但非充分")
else:
    print("  🔴 兩長度（單一 r、單一 t）仍不足 ⟹ 缺口更深（需 mode 佔用／非成對項）")

with open("data/35_two_lengths.json", "w") as f:
    json.dump({"R_data": R_DATA, "rows": rows,
               "pass_interval": [lo, hi] if good else None,
               "best": b, "candidates_posthoc": cands, "r_star": tgt},
              f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/35_two_lengths.json")