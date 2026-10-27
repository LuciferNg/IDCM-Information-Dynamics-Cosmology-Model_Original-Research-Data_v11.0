"""M5 docs/32 配套 — 第三前線：配對勢之層衰減律 ρ（計數 vs 度量）
================================================================================
現行模型（`docs/09` 起）：`g(d) = −α q (ħc/a) B^{−d}`，`B = 4`。
但樹嘅幾何有兩個唔同嘅數：
  · 計數：`N_d = B (B−1)^{d−1}` ⟹ `B = 4` 係**鏈配位數**（每站 4 鄰）
  · 度量：分支因子 `b−1 = 3` ⟹ 層距自然標度 `r_d = a·3^d` ⟹ ρ = 1/3
  · `K_tree` 帶中 resolvent：Bethe 格衰減 `(b−1)^{−d/2} = 3^{−d/2}` ⟹ ρ = 1/√3
用 `4^d` 做距離而唔係 `3^d`，等於**用計數當度量**。呢個係結構問題，唔係擬合。

判準（事前寫死）：
  T1 ρ 只可取三個結構值：`1/4`（現行，鏈配位）、`1/3`（樹分支＝度量）、`3^{−1/2}`（K_tree 帶中）。不得取擬合值。
  T2 每變體以同一 D 兩數（A_D = 1.2025、|c|_D = 0.3525）校準 (a, t)；K／Ξc 零新參數。
  T3 通過 = K 盲測 `E(K⁰K⁻)` 與 `E(K*⁰K*⁺)` 皆 |z| ≤ 2（另報 `ΔE_K`）。
  T4 「顯著收窄」門檻（預註冊）：`R = r_K/r_D` 須 ≥ R_data/3 = 9.65。
  T5 只改 ρ 一個元素；跳躍結構（`−√b`、`−√(b−1)`）、contact 指數 3、質量、q 全部不動。
輸出：`data/34_metric_decay.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
A_D, C_D = 1.2025, 0.3525
R_DATA = (98.77777777777769 / (A_D / C_D))
OBS = {"E1": 1.3740, "E2": 1.3200, "dE": 0.0540, "s": 0.3545}
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/32 — 第三前線：配對勢之層衰減律 ρ（計數 4 vs 度量 3 vs 帶中 3^−1/2）")
print("=" * 104)
print(f"\n數據要求 R = r_K/r_D = {R_DATA:.4f}；顯著收窄門檻 R ≥ {R_DATA/3:.2f}")


def make(rho):
    POW = [rho ** d for d in range(L + 1)]

    def amps(a, t, q, mh, ml):
        g = ALPHA * q * HBARC / a / t
        Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
        h = Cq * (HBARC / a) ** 3 / t
        diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
        _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
        v = evec[:, 0] if evec.ndim == 2 else evec
        w = v ** 2
        w /= w.sum()
        rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
        dpk = int(np.argmax(w))
        bd = bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)
        return (ALPHA * q * (HBARC / a) * rA * 1000,
                Cq * (HBARC / a) ** 3 * P0 * 1000,
                {"rA": rA, "P0": P0, "dpk": dpk, "bound": bd})
    return amps


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


def calib(amps, ngrid=400):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), Q_C, M_C, M_L)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
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
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


variants = [("現行 ρ = 1/4（鏈配位＝計數）", 0.25),
            ("度量 ρ = 1/3（樹分支 b−1）", 1 / 3),
            ("帶中 ρ = 3^−1/2（K_tree resolvent）", 3 ** -0.5)]
out = {}
for name, rho in variants:
    print(f"\n[變體] {name}  （ρ = {rho:.9f}）")
    amps = make(rho)
    sol = calib(amps)
    if sol is None:
        print("     🔴 無法以 D 兩數校準")
        out[name] = {"rho": rho, "calibrated": False}
        continue
    a, t = sol
    AD, cD, stD = amps(a, t, Q_C, M_C, M_L)
    AK, cK, stK = amps(a, t, Q_S, M_S, M_L)
    E1, E2, dE = AK + 3 * cK, AK - cK, 4 * cK
    z = {k: (v - OBS[k]) / OBS["s"] for k, v in (("E1", E1), ("E2", E2), ("dE", dE))}
    R = (AK / cK) / (AD / cD)
    ok = abs(z["E1"]) <= 2 and abs(z["E2"]) <= 2
    print(f"     校準：a = {a:.6f} fm、t = {t*1000:.4f} MeV（D 零殘差）⟨B⟩_D = {stD['rA']:.4f}、P₀_D = {stD['P0']:.4f}")
    print(f"     K：⟨B⟩ = {stK['rA']:.4f}、P₀ = {stK['P0']:.4f}、d_peak = {stK['dpk']}、束縛 {stK['bound']}"
          f"  ⟹ r_K = {AK/cK:.4f}（需 98.7778）、r_D = {AD/cD:.4f}（需 3.4113）")
    print(f"     盲測：E(K⁰K⁻) = {E1:.4f} (z {z['E1']:+.2f})、E(K*⁰K*⁺) = {E2:.4f} (z {z['E2']:+.2f})、"
          f"ΔE_K = {dE:.4f} (z {z['dE']:+.2f})")
    print(f"     R = r_K/r_D = {R:.4f}（需 {R_DATA:.2f}）⟹ 缺 {R_DATA/R:.2f}×；T3 {'✅' if ok else '🔴'}；"
          f"T4 {'✅' if R >= R_DATA/3 else '🔴'}")
    out[name] = {"rho": rho, "calibrated": True, "a0_fm": a, "t_MeV": t * 1000,
                 "rD": AD / cD, "rK": AK / cK, "R": R, "z": z, "T3": bool(ok),
                 "T4": bool(R >= R_DATA / 3), "K": stK, "D": stD}

cal = {k: v for k, v in out.items() if v.get("calibrated")}
best = max(cal, key=lambda k: cal[k]["R"]) if cal else None
print("\n[判決]")
if best:
    print(f"     最佳 = {best}（ρ = {cal[best]['rho']:.6f}）⟹ R = {cal[best]['R']:.4f}、缺 {R_DATA/cal[best]['R']:.2f}×")
    impr = cal[best]["R"] / cal["現行 ρ = 1/4（鏈配位＝計數）"]["R"]
    print(f"     相對現行改善 = {impr:.3f}×；T3 通過者 = {[k for k,v in cal.items() if v['T3']] or '無'}")
    if cal[best]["T4"]:
        print("     🟢 達預註冊收窄門檻 ⟹ 須全鏈重算（A3／D5／天花板）")
    else:
        print("     🔴 未達預註冊收窄門檻 ⟹ 衰減律唔係缺口所在（本輪只改 ρ 一個元素）")

with open("data/34_metric_decay.json", "w") as f:
    json.dump({"R_data": R_DATA, "variants": out, "best": best}, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/34_metric_decay.json")