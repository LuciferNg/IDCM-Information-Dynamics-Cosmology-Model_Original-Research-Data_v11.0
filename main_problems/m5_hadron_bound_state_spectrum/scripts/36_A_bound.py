"""M5 docs/34 配套 — 第五前線：A 側先驗上界定理（證明類別封閉）＋ 逃生門量測
================================================================================
`docs/33` 掃描：`r = a_C/a_X` 由 0.30 到 1.00，`R` 只在 0.429–0.434（±0.6%），`A_K` 全部 ≪ 1.3335。
呢個暗示一條與所有幾何參數無關嘅上界。本輪證它，並量度唯一逃生門。

【定理 P】（純代數，兩扇區模型）
  `A = α q (ħc/a)⟨B⟩ + κ g(0) q`（tree ＋ chain 旁觀者），`⟨B⟩ ≤ 1`，`D` 之 `A_D` 已由校準固定。
  ⟹ `A_K/A_D = (q_K/q_D)·(1+Z)/(⟨B⟩_D+Z)`，其中 `Z ≡ κ g(0) a/(ħc) ≥ 0`。
  因 `⟨B⟩_D ≤ 1`，`(1+Z)/(⟨B⟩_D+Z)` 對 `Z` 單調遞減 ⟹ 上界在 `Z = 0` 取到：

      A_K ≤ (q_K/q_D)·A_D/⟨B⟩_D = (1/2)·1.2025/0.4812 = 1.2499 MeV

  與 `a_C`、`a_X`、`t`、ρ、`p`、支撐、`κ`、權重約定**全部無關**（只因 `q_K/q_D = 1/2` 與 `⟨B⟩ ≤ 1`）。
  觀測要求 `A_K = (E1 + 3E2)/4 = 1.3335` ⟹ 差 6.7% ⟹ **本類別在 A 側封閉**。

【逃生門】上界 ∝ 1/⟨B⟩_D ⟹ 若 D 之狀態更離域（⟨B⟩_D 更細）則上界升高。臨界值：

      ⟨B⟩_D ≤ (q_K/q_D)·A_D/A_K = 0.45089

本輪量測：(r, κ) 二維族（每點仍以 D 兩數校準）之 min⟨B⟩_D。

判準（事前寫死）：
  P1 定理只准用恆等式（`⟨B⟩ ≤ 1`、`A` 之定義），不得引入近似。
  P2 量測網格寫死：r ∈ [0.30,1.00] 步 0.02（36 點）、κ ∈ {0, 0.05, 0.10, 0.15, 0.20, 0.26986, 0.35}（7 點）。
  P3 逃生門成立 ⟺ min⟨B⟩_D < 0.45089（嚴格），且該點須束縛。
  P4 只報事實，不得為救模型改 q 比值（`q_K/q_D = 1/2` 係電荷指派）。
輸出：`data/36_A_bound.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0 = 2581 / 352
A_D, C_D = 1.2025, 0.3525
OBS = {"E1": 1.3740, "E2": 1.3200, "s": 0.3545}
A_K_REQ = (OBS["E1"] + 3 * OBS["E2"]) / 4
C_K_REQ = (OBS["E1"] - OBS["E2"]) / 4
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/34 — 定理 P：A 側先驗上界；逃生門 = ⟨B⟩_D ≤ 0.45089")
print("=" * 104)
print(f"\n觀測要求：A_K = (E1+3E2)/4 = {A_K_REQ:.5f}、c_K = (E1−E2)/4 = {C_K_REQ:.5f}"
      f"（r_K = {A_K_REQ/C_K_REQ:.4f}）")


def amps(aC, aX, t, kappa, q, mh, ml):
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aX) ** 3 / t
    chain = kappa * G0 * q / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) - chain for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    return (ALPHA * q * (HBARC / aC) * rA * 1000 + kappa * G0 * q * 1000,
            Cq * (HBARC / aX) ** 3 * P0 * 1000, rA, P0, dpk,
            bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX))


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


def calib(r, kappa, ngrid=280):
    def t_solve(aC):
        x = bisect(lambda xx: amps(aC, aC / r, float(np.exp(xx)), kappa, Q_C, M_C, M_L)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        aC = float(np.exp(x))
        t = t_solve(aC)
        if t is None:
            prev = None
            continue
        k = amps(aC, aC / r, t, kappa, Q_C, M_C, M_L)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            rr = bisect(lambda xx: amps(float(np.exp(xx)), float(np.exp(xx)) / r,
                                        t_solve(float(np.exp(xx))), kappa, Q_C, M_C, M_L)[1] - C_D,
                        prev[0], x)
            if rr is not None:
                aa = float(np.exp(rr))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


print("\n[定理 P] 上界（Z = 0 取到，與所有幾何參數無關）")
for kappa in (0.0, 0.26986):
    pass
# 用實際 ⟨B⟩_D 計算（下面掃描會填充）；先報解析式
print(f"     A_K ≤ (q_K/q_D)·A_D/⟨B⟩_D = (1/2)·{A_D}/⟨B⟩_D")
print(f"     要求 A_K = {A_K_REQ:.5f} ⟹ 需 ⟨B⟩_D ≤ {0.5*A_D/A_K_REQ:.5f}")
print(f"     現行 v2 校準 ⟨B⟩_D = 0.48120 ⟹ 上界 = {0.5*A_D/0.48120:.5f} MeV（差 {(A_K_REQ/(0.5*A_D/0.4812)-1)*100:.2f}%）")

rs = np.arange(0.30, 1.0001, 0.02)
kappas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.26986, 0.35]
print(f"\n[P2] 掃描 {len(rs)}×{len(kappas)} = {len(rs)*len(kappas)} 點（每點以 D 兩數校準）")
best = None
rows = []
for kappa in kappas:
    for r in rs:
        sol = calib(float(r), kappa)
        if sol is None:
            continue
        aC, t = sol
        AD, cD, rAD, pD, dpkD, bdD = amps(aC, aC / r, t, kappa, Q_C, M_C, M_L)
        A_K, c_K, rAK, pK, dpkK, bdK = amps(aC, aC / r, t, kappa, Q_S, M_S, M_L)
        bound = (Q_S / Q_C) * A_D / rAD
        row = {"r": float(r), "kappa": kappa, "aC_fm": aC, "aX_fm": aC / r, "t_MeV": t * 1000,
               "B_D": rAD, "P0_D": pD, "bound_A_K": bound, "A_K": A_K, "c_K": c_K,
               "rK": A_K / c_K if c_K else np.nan, "R": (A_K / c_K) / (AD / cD) if c_K else np.nan,
               "pass_A": bool(A_K >= A_K_REQ), "bound_ok": bool(bound >= A_K_REQ), "bdD": bdD}
        rows.append(row)
        if best is None or rAD < best["B_D"]:
            best = row
print(f"\n[P3] min⟨B⟩_D = {best['B_D']:.5f} @ (r = {best['r']:.2f}, κ = {best['kappa']})"
      f"；用此 ⟹ A_K 上界 = {0.5*A_D/best['B_D']:.5f} MeV")
print(f"     逃生門門檻 ⟨B⟩_D < {0.5*A_D/A_K_REQ:.5f} ⟹ "
      f"{'🟢 逃生門存在' if best['B_D'] < 0.5*A_D/A_K_REQ else '🔴 逃生門封閉'}")
print(f"     全掃描中 A_K ≥ A_K_REQ 之點數 = {sum(1 for x in rows if x['pass_A'])} / {len(rows)}")
print(f"     全掃描中 R ≥ 9.65 之點數 = {sum(1 for x in rows if np.isfinite(x['R']) and x['R'] >= 9.65)} / {len(rows)}")
print(f"     R 範圍（有 c_K>0 者）：[{min(x['R'] for x in rows if np.isfinite(x['R'])):.3f}, "
      f"{max(x['R'] for x in rows if np.isfinite(x['R'])):.3f}]；需 {A_K_REQ/C_K_REQ/(A_D/C_D):.2f}")

with open("data/36_A_bound.json", "w") as f:
    json.dump({"A_K_req": A_K_REQ, "c_K_req": C_K_REQ, "bound_formula": "A_K <= (q_K/q_D)*A_D/<B>_D",
               "B_D_critical": 0.5 * A_D / A_K_REQ, "min_B_D": best, "rows": rows},
              f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/36_A_bound.json")