"""M5 docs/17 配套 — 核更正第一驗：用 H1 精確核 g(d) 取代指數階梯 B^{-d}
================================================================================================
觸發：`docs/16` Q5 之「硬上界」`A_max = α q(ħc/a_C)·1` 其實只對 **指數核** `B^{-d}` 成立；
而 M5 自己嘅 H1（`docs/01`）已給出 **精確** 的 ⟨1/r⟩ 鏈上和：

    g(d) = [(N²−1) − 6d(N−d)]/(12N)   （Fraction 精確；N = 88）

`g(0) = 2581/352 = 7.332386…`、`g(7) = 1447/352 = 4.110795…`、`g(20) = −139/352`、`g(44) = −1291/352`。
指數核在 d = 7 為 `B^{−7} = 6.10e−5` ⟹ 與 `g(7)` 相差 **6.7×10⁴** ⟹ 舊振幅計算把跨層 Coulomb **壓死**。

判準（寫死）：
  R1 只用 D 通道（A_D = 1.2025、|c|_D = 0.3525）校準 (a_C, t)；K 系全盲。
  R2 命中 = K 兩通道 |z| ≤ 2（σ = 0.0961／0.3412）。
  R3 束縛判別沿用（⟨g⟩ > 1e−3 且 d_peak ≤ 20）。
  R4 若 D 通道校準無解（A/K 目標不可達）⟹ 標 N/A，不得報 0。
  R5 對照：指數核之舊結果（A_K = 0.2903、|c|_K = 0.1880、z = −5.41／−3.57）同場重列。

輸出：`data/18_g_kernel_first_pass.json`
"""
import json
from fractions import Fraction

import numpy as np

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
N = 88


def gexact(d):
    return (Fraction(N * N - 1) - 6 * d * (N - d)) / (12 * N)


GK = [float(gexact(d)) for d in range(L + 1)]
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L = 1.27, 0.495, 0.30
Q_C, Q_S = 2 / 3, 1 / 3


def shape(gcoef, h, kern):
    diag = np.array([-gcoef * kern[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    return float((w * np.array(kern)).sum()), float(w[0]), int(np.argmax(w)), bool(w[0] > 1e-3 and int(np.argmax(w)) <= 20)


def amps(aC, aK, t, q, mh, ml, kern):
    gcoef = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aK) ** 3 / t
    rg, P0, dp, bd = shape(gcoef, h, kern)
    return ALPHA * q * (HBARC / aC) * rg * 1000, Cq * (HBARC / aK) ** 3 * P0 * 1000, (rg, P0, dp, bd)


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


def calibrate(kern):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, a, float(np.exp(xx)), Q_C, M_C, M_L, kern)[0] - 1.2025,
                   np.log(1e-9), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev, sol = None, None
    for x in np.linspace(np.log(1e-3), np.log(5.0), 400):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, a, t, Q_C, M_C, M_L, kern)[1] - 0.3525
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), float(np.exp(xx)), t_solve(float(np.exp(xx))),
                                       Q_C, M_C, M_L, kern)[1] - 0.3525, prev[0], x)
            if r is not None:
                sol = (float(np.exp(r)), t_solve(float(np.exp(r))))
                break
        prev = (x, k)
    return sol


print("=" * 104)
print("M5 docs/17 — 核更正第一驗：g(d) 精確核 vs 指數核 B^{-d}")
print("=" * 104)
print("\n[R0] 兩核對照（Fraction 精確）")
print(f"     g(0) = {gexact(0)} = {GK[0]:.6f}；g(7) = {gexact(7)} = {GK[7]:.6f}；"
      f"g(20) = {gexact(20)} = {float(gexact(20)):.6f}；g(40) = {gexact(40)} = {GK[40]:.6f}；"
      f"g(44)〔鏈對極〕= {gexact(44)} = {float(gexact(44)):.6f}")
print(f"     指數核：B^0 = 1、B^-7 = {POW[7]:.3e} ⟹ g(7)/B^-7 = {GK[7]/POW[7]:.3e}")
print(f"     ⟹ docs/16 Q5 之「上限 A_max = α q(ħc/a_C)」只對指數核成立；g 核之上限 ×{GK[0]:.4f}")

for name, kern in (("g 核（H1 精確）", GK), ("指數核（舊）", POW)):
    print(f"\n[{name}] 校準（只用 D 通道）")
    sol = calibrate(kern)
    if sol is None:
        print("     🔴 D 通道校準無解（A/|c| 目標不可達）⟹ N/A，不報 0")
        continue
    a0, t0 = sol
    A_D, c_D, st_D = amps(a0, a0, t0, Q_C, M_C, M_L, kern)
    print(f"     a₀ = {a0:.6f} fm、t = {t0*1000:.6f} MeV ⟹ A_D = {A_D:.6f}、|c|_D = {c_D:.6f}"
          f"、⟨g or B⟩ = {st_D[0]:.6f}、P₀ = {st_D[1]:.6f}、束縛 = {st_D[3]}")
    A_K, c_K, st_K = amps(a0, a0, t0, Q_S, M_S, M_L, kern)
    E1, E2 = A_K + 3 * c_K, A_K - c_K
    z1 = (E1 - 1.3740) / 0.0960728889958036
    z2 = (E2 - 1.3200) / 0.34122719704032967
    print(f"     K 系盲測：A_K = {A_K:.4f}、|c|_K = {c_K:.4f}、⟨g⟩ = {st_K[0]:.4f}、P₀ = {st_K[1]:.4f}"
          f"、束縛 = {st_K[3]}")
    print(f"     E(K⁰−K⁺) = {E1:.4f} vs 1.3740（z = {z1:+.2f}）；"
          f"E(K*⁰−K*⁺) = {E2:.4f} vs 1.3200（z = {z2:+.2f}）⟹ "
          f"{'✅ 命中' if st_K[3] and abs(z1) <= 2 and abs(z2) <= 2 else '🔴 不中'}")
    if name.startswith("g"):
        out = {"a0_fm": a0, "t_MeV": t0 * 1000, "A_D": A_D, "c_D": c_D, "rg_D": st_D[0], "P0_D": st_D[1],
               "A_K": A_K, "c_K": c_K, "rg_K": st_K[0], "P0_K": st_K[1], "bound_K": st_K[3],
               "E1": E1, "z1": z1, "E2": E2, "z2": z2,
               "pass": bool(st_K[3] and abs(z1) <= 2 and abs(z2) <= 2)}

out["g_table"] = {f"g({d})": [str(gexact(d)), float(gexact(d))] for d in (0, 1, 7, 20, 40, 44)}
out["cap_ratio"] = GK[0]
with open("data/18_g_kernel_first_pass.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/18_g_kernel_first_pass.json")