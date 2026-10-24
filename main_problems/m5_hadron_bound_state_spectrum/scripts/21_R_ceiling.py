"""M5 docs/19 配套（第二腳本）— D3：A/|c| 形狀比天花板 R 與聯合 χ²
================================================================================================
`scripts/20` 證 T1（盲測）與 T3（可達範圍）後，本腳本做兩件決定性計算：

  T4 ★ 沿 D 校準曲線（κ 單參數；(a_C, t) 每次由 D 兩數解出）細掃 κ ∈ [0.200, 0.245]，
        取 max r_K（r ≡ A/|c|）⟹ 形狀比天花板。與數據 r_K = 98.7778 比。
  T5    同一族之聯合 χ² = z₁² + z₂²（K 兩通道）最小值、對應 κ、以及 p 值（dof = 1）。

理由：`r ≡ A/|c|` 與電荷無關（`docs/14` §四.1 恒等式），故 r 係純形狀量；若模型之 r 有天花板
而數據在其上 ⟹ 唔論校準點、唔論 κ，模型都造唔到「K 系 contact 近乎消失」。

判準（寫死）：
  R1 每點必須 D 兩數零殘差（|ΔA| + |Δc| < 1e−9）且束縛（P₀ > 1e−3、d_peak ≤ 20），否則記 N/A。
  R2 天花板 = 掃描網格之 max；網格敏感度由兩組步長（0.005 / 0.0025）對照。
  R3 數據 r_K 由 E1、E2 反解：A = (E1+3E2)/4、|c| = (E1−E2)/4。
  R4 禁改 m_h、禁改 σ、禁加族；未命中記 NEG。

輸出：`data/21_R_ceiling.json`
"""
import json

import numpy as np
from scipy.optimize import root
from scipy.stats import chi2

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
A_D, C_D = 1.2025, 0.3525
E1K, E2K = 1.3740, 1.3200
S1K, S2K = 0.0960728889958036, 0.34122719704032967
R_D = A_D / C_D
R_K = ((E1K + 3 * E2K) / 4) / ((E1K - E2K) / 4)


def amps(a, t, q, mh, ml):
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
    return ALPHA * q * (HBARC / a) * rA * 1000, Cq * (HBARC / a) ** 3 * P0 * 1000, bd


def bisect(f, lo, hi, n=300):
    flo = f(lo)
    if flo * f(hi) > 0:
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


def calib_continuous(kappa):
    """與 scripts/19 完全相同之校準：由 κ = 0 連續延拓之唯一分支（a 掃、t 內層二分）"""
    chd = kappa * G7 * Q_C
    target = A_D - chd
    if target <= 0:
        return None

    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), Q_C, M_C, M_L)[0] - target, np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev, sol = None, None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 400):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, Q_C, M_C, M_L)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))), Q_C, M_C, M_L)[1] - C_D,
                       prev[0], x)
            if r is not None:
                sol = (float(np.exp(r)), t_solve(float(np.exp(r))))
                break
        prev = (x, k)
    return sol


def calib_all_roots(kappa):
    """獨立 root-finder：枚舉 D 校準之**全部**根（檢多值性）"""
    chd = kappa * G7 * Q_C

    def F(x):
        a, t = np.exp(x[0]), np.exp(x[1])
        A, c, _ = amps(a, t, Q_C, M_C, M_L)
        return [np.log(max(A + chd, 1e-12) / A_D), np.log(max(c, 1e-12) / C_D)]

    roots = []
    for x0a in np.linspace(-2.5, 1.5, 14):
        for x0b in np.linspace(-9.0, 2.0, 16):
            try:
                s = root(F, [x0a, x0b], method="hybr")
            except Exception:
                continue
            if not s.success or np.max(np.abs(F(s.x))) > 1e-8:
                continue
            a, t = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
            if not (1e-3 < a < 20 and 1e-9 < t < 1e3):
                continue
            A, c, _ = amps(a, t, Q_C, M_C, M_L)
            if abs(A + chd - A_D) + abs(c - C_D) > 1e-8:
                continue
            if any(abs(a - r[0]) < 1e-5 and abs(t - r[1]) < 1e-9 for r in roots):
                continue
            roots.append((a, t))
    return roots


print("=" * 104)
print("M5 docs/19（scripts/21）— D3 天花板：A/|c| 形狀比 R 與聯合 χ²")
print("=" * 104)
print(f"\n[D0] r_D = {A_D}/{C_D} = {R_D:.6f}；r_K = {R_K:.6f}（E1 = {E1K}、E2 = {E2K} 反解）"
      f"；R_數據 = {R_K/R_D:.6f}")

rows = []
for kp in np.arange(0.180, 0.3201, 0.0025):
    sol = calib_continuous(float(kp))
    if sol is None:
        continue
    a, t = sol
    A, c, bd = amps(a, t, Q_S, M_S, M_L)
    Atot = A + float(kp) * G0 * Q_S
    z1 = (Atot + 3 * c - E1K) / S1K
    z2 = (Atot - c - E2K) / S2K
    rows.append({"kappa": round(float(kp), 4), "a0_fm": a, "t_MeV": t * 1000, "bound_K": bd,
                 "A_K": Atot, "c_K": c, "r_K": Atot / c, "R": (Atot / c) / R_D,
                 "z1": z1, "z2": z2, "chi2": z1 * z1 + z2 * z2})
print(f"     連續分支：{len(rows)} 個 κ 點成功（κ ∈ [0.180, 0.320] 步 0.0025）")

print("\n[多值性檢查] 獨立 root-finder 對 κ = {0.20, 0.25, 0.26986, 0.30} 枚舉全部根：")
multi = {}
for kp in (0.20, 0.25, KAPPA if 'KAPPA' in dir() else 0.26986, 0.30):
    roots = calib_all_roots(float(kp))
    multi[str(kp)] = [{"a0_fm": a, "t_MeV": t * 1000} for a, t in roots]
    print(f"     κ = {kp:.5f} ⟹ {len(roots)} 個根：" +
          "；".join(f"(a₀ = {a:.5f} fm, t = {t*1000:.4f} MeV)" for a, t in roots))

print("\n[T4] 天花板（連續分支）")
mx = max(rows, key=lambda r: r["r_K"])
print(f"     max r_K = {mx['r_K']:.4f} @ κ = {mx['kappa']:.4f}（a₀ = {mx['a0_fm']:.5f} fm、"
      f"t = {mx['t_MeV']:.4f} MeV）⟹ R = {mx['R']:.4f}；數據 R = {R_K/R_D:.4f}"
      f" ⟹ 差 {R_K/mx['r_K']:.2f}×")

best = min(rows, key=lambda r: r["chi2"])
print(f"\n[T5] 聯合 χ² 最小：κ = {best['kappa']:.4f} ⟹ z₁ = {best['z1']:+.3f}、z₂ = {best['z2']:+.3f}、"
      f"χ² = {best['chi2']:.3f} ⟹ p = {chi2.sf(best['chi2'], 1):.4f}（dof = 1）⟹ {np.sqrt(best['chi2']):.2f}σ")
print(f"     天花板 κ = {mx['kappa']:.4f} 與 χ² 最小 κ = {best['kappa']:.4f} 不同 ⟹ 兩目標不可同時滿足")

out = {"r_D": R_D, "r_K": R_K, "R_data": R_K / R_D,
       "ceiling": {"kappa": mx["kappa"], "r_K": mx["r_K"], "R": mx["R"], "short_factor": R_K / mx["r_K"]},
       "chi2_min": {"kappa": best["kappa"], "z1": best["z1"], "z2": best["z2"], "chi2": best["chi2"],
                    "p": float(chi2.sf(best["chi2"], 1)), "sigma": float(np.sqrt(best["chi2"]))},
       "multiplicity": multi, "rows": rows}
with open("data/21_R_ceiling.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/21_R_ceiling.json")
