"""M5 docs/30 配套 — 攻天花板：contact 算子之 scaling dimension 與支撐（K_tree 幾何）
================================================================================
病灶（`docs/28` §五／`docs/29` §五）：模型要 `|c|_K/|c|_D = 0.0383`，實得 0.2465 ⟹ 天花板 17.58×。
`|c| ∝ C_q (ħc/a)³ P₀`。本輪攻**指數 3**：3 係 3D 空間維度，但 contact 項住喺 **K_tree（超度量）** 上，
其 Hausdorff 維度唔係 3。若支撐空間維度唔同，`δ`-函數之標度律 `a^{−d_H}` 亦唔同 ⟹ 呢個係**結構問題**，唔係擬合。

  Z1 由樹自己算維度：層 d 之站數 N_d = b(b−1)^{d−1}、層距 r_d = B^{−d}
     ⟹ N(r) ~ r^{−d_H}，d_H = ln(b−1)/ln B = ln3/ln4（精確，非擬合）。
  Z2 變體 A：contact 指數 p = d_H（樹維度）取代 3；其餘不變、零新參數。
  Z3 變體 B：contact 支撐 = 球（d ≤ 1，權重 1/2 各層）而非單站 d = 0。
  Z4 變體 C（對照）：p = 3（現行）。
  全部以同一 D 兩數（A_D = 1.2025、|c|_D = 0.3525）校準 (a, t)，K 通道零新參數盲測。

判準（事前寫死，不得事後放寬）：
  Y1 d_H 必須由 N_d 序列導出（不得取任意實數）。
  Y2 每變體只准改**一個**結構元素（指數 或 支撐），不得同時改兩個。
  Y3 通過門檻：K 盲測 `|ΔE_K|` 與 `E(K*⁰−K*⁺)` 皆 |z| ≤ 2。
  Y4 「收窄缺口」門檻（預註冊）：`R_max` 須 ≥ R_data/3 = 9.61 才算顯著；否則記「不足」。
  Y5 不得引入新參數、不得改 q／m／κ／D 目標；非束縛態一律 N/A。
輸出：`data/32_dimension_attack.json`
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
OBS = {"E1": (1.3740, 0.3545), "E2": (1.3200, 0.3545), "dE": (0.0540, 0.3545)}
R_DATA_K, R_DATA_D = 98.77777777777769, A_D / C_D
R_DATA = R_DATA_K / R_DATA_D
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/30 — 攻天花板：contact 算子之 scaling dimension 與支撐")
print("=" * 104)

# ---------------- Z1 ----------------
N_d = [1] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
dH = np.log(B - 1) / np.log(B)
print("\n[Z1] 由樹算 Hausdorff 維度（層站數 vs 層距）")
print(f"     N_d：d=1..5 ⟹ {[N_d[d] for d in range(1, 6)]}")
slopes = [np.log(N_d[d + 1] / N_d[d]) / np.log(B) for d in range(1, L)]
print(f"     ln(N_(d+1)/N_d)/ln B = {slopes[0]:.12f}（全層同值：{max(slopes)-min(slopes):.1e} 波動）")
print(f"     ⟹ d_H = ln(b−1)/ln B = ln3/ln4 = {dH:.12f}")
print(f"     ⟹ 現行模型用 (ħc/a)³（3D）；樹維度 {dH:.4f} ⟹ **維度不自洽**（呢個係結構缺陷，非調參）")


def amps(a, t, q, mh, ml, p=3.0, support="site"):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** p / t
    diag = np.empty(L + 1)
    for d in range(L + 1):
        diag[d] = -g * POW[d]
    if support == "site":
        diag[0] -= h
    elif support == "ball1":
        diag[0] -= 0.5 * h
        diag[1] -= 0.5 * h
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    bd = bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)
    A = ALPHA * q * (HBARC / a) * rA * 1000
    c = Cq * (HBARC / a) ** p * P0 * 1000
    return A, c, {"rA": rA, "P0": P0, "dpk": dpk, "bound": bd}


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


def calib(p, support, ngrid=400):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), Q_C, M_C, M_L, p, support)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, Q_C, M_C, M_L, p, support)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))), Q_C, M_C, M_L, p, support)[1] - C_D,
                       prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


# ---------------- Z2-Z4 ----------------
variants = [("C 現行 p=3、單站", 3.0, "site"), ("A 樹維度 p=d_H、單站", dH, "site"),
            ("B p=3、球支撐 d≤1", 3.0, "ball1")]
res = {}
for name, p, sup in variants:
    print(f"\n[變體] {name}")
    sol = calib(p, sup)
    if sol is None:
        print("     🔴 D 兩數無解（該變體不可校準）")
        res[name] = {"calibrated": False}
        continue
    a, t = sol
    AD, cD, stD = amps(a, t, Q_C, M_C, M_L, p, sup)
    AK, cK, stK = amps(a, t, Q_S, M_S, M_L, p, sup)
    E1, E2, dE = AK + 3 * cK, AK - cK, 4 * cK
    z1 = (E1 - OBS["E1"][0]) / OBS["E1"][1]
    z2 = (E2 - OBS["E2"][0]) / OBS["E2"][1]
    zd = (dE - OBS["dE"][0]) / OBS["dE"][1]
    rD, rK = AD / cD, AK / cK
    print(f"     校準：a = {a:.6f} fm、t = {t*1000:.4f} MeV（D 零殘差）、⟨B⟩_D = {stD['rA']:.4f}、P₀_D = {stD['P0']:.4f}")
    print(f"     K 態：⟨B⟩ = {stK['rA']:.4f}、P₀ = {stK['P0']:.4f}、d_peak = {stK['dpk']}、束縛 {stK['bound']}")
    print(f"     K 盲測：ΔE_K = {dE:.4f}（z = {zd:+.2f}）、E(K⁰K⁻) = {E1:.4f}（z = {z1:+.2f}）、"
          f"E(K*⁰K*⁻) = {E2:.4f}（z = {z2:+.2f}）")
    print(f"     形狀比：R = r_K/r_D = {rK/rD:.4f}（數據需 {R_DATA:.2f}）、r_K = {rK:.4f} vs 98.7778")
    res[name] = {"calibrated": True, "p": p, "support": sup, "a0_fm": a, "t_MeV": t * 1000,
                 "rD": rD, "rK": rK, "R": rK / rD, "z_dE": zd, "z_E1": z1, "z_E2": z2,
                 "K_state": stK, "D_state": stD,
                 "Y3_pass": bool(abs(zd) <= 2 and abs(z2) <= 2),
                 "Y4_pass": bool(rK / rD >= R_DATA / 3)}

print("\n[Z5] 逐變體判（Y3 盲測門檻、Y4 收窄門檻 R_max ≥ R_data/3 = %.2f）" % (R_DATA / 3))
for name, v in res.items():
    if not v.get("calibrated"):
        print(f"     {name:22s} 不可校準")
        continue
    print(f"     {name:22s} R = {v['R']:6.3f}（需 {R_DATA:.2f}）⟹ 缺 {R_DATA/v['R']:6.2f}×；"
          f"max|z| = {max(abs(v['z_dE']), abs(v['z_E1']), abs(v['z_E2'])):.2f}；"
          f"Y3 {'✅' if v['Y3_pass'] else '🔴'}；Y4 {'✅' if v['Y4_pass'] else '🔴'}")

best = max((v for v in res.values() if v.get("calibrated")), key=lambda v: v["R"])
print(f"\n[判決] 最佳變體 R = {best['R']:.4f}（p = {best['p']:.6f}、支撐 {best['support']}）"
      f" ⟹ 缺口由 {R_DATA / res['C 現行 p=3、單站']['R']:.2f}× 變 {R_DATA / best['R']:.2f}×")
verdict = ("🟢 有結構改善但仍未閉合" if best["R"] > res["C 現行 p=3、單站"]["R"] * 1.2
           else "🔴 兩種結構修正皆無效 ⟹ 缺口不在 contact 標度維度／支撐")
print("     " + verdict)

out = {"d_H": dH, "d_H_formula": "ln(b-1)/ln B = ln3/ln4",
       "variants": res, "R_data": R_DATA, "R_data_over_3": R_DATA / 3,
       "best": best.get("p"), "verdict": verdict}
with open("data/32_dimension_attack.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/32_dimension_attack.json")