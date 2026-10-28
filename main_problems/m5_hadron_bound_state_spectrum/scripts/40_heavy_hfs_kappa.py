"""M5 docs/38 配套 — 重區擴展（ΔI=0 超精細）＋ κ 兩約定統一
================================================================================
路線 E5／E6。本腳本做兩件互不相關但同屬「收尾」嘅事：

【一｜ΔI=0 超精細通道（新盲測）】`docs/35` 只測 ΔI=1 電荷通道（`c = (Δ1−Δ0)/4`）。
本輪加第二通道：`E_hfs ≡ M(J=1) − M(J=0)`（同電荷，故純自旋依賴）⟹ `E_hfs = 4|c|`。
模型零新參數（同一 `(a,t)`），比較 `R_s ≡ c(D_s)/c(D)` 與 `c(B_s)/c(B)`。

【二｜κ 約定統一】`docs/28`（`scripts/29`）把 `κ` 只加在**回報之 `r` 分子**；`scripts/36` 放入 Hamiltonian。
本輪定**正典約定**＝放入 Hamiltonian（因 `κ` 係勢之常數位移，物理上屬 H），並給該約定下：
  · 波函數不變（常數位移）⟹ 只有 `A`、`c` 平移
  · 校準硬界 `κ < A_D/(g(0)q_D)`（`⟨B⟩_D ≤ 1`）
  · `κ ∈ [0, 界)` 全掃之 `R_max`（與 `docs/28` 之 1.647 對照）

判準（事前寫死）：
  K1 `E_hfs` 一律用同電荷對（`D*⁰−D⁰`、`D_s*±−D_s±`、`B*⁰−B⁰`、`B_s*⁰−B_s⁰`），PDG 2024 直讀。
  K2 模型側 `c = C_q(ħc/a)³P₀`（`docs/09`）＋同一 D 校準 `(a₀,t)`；零新參數。
  K3 比值判：`|ln(R_model/R_data)| ≤ ln(1.5)` 記一致；否則報偏離倍數。
  K4 診斷：對 `R_data` 反解「所需 `m_l` 指數」`α`，使 `[P₀(D_s)/P₀(D)]·(m_l,D/m_l,Ds)^α = R_data`。
  K5 `κ` 正典約定下掃 25 點求 `R_max`；並列 `docs/28` 之回報式值作對照（不得混算）。
  K6 不得改 `q`、`m`、`X`、校準目標；非束縛態 N/A。
輸出：`data/40_heavy_hfs_kappa.json`
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
A0, T0 = 0.384184, 1.7711e-3
A_D, C_D = 1.2025, 0.3525
G0 = 2581 / 352
BOUND_RHO, DPEAK_MAX = 1e-3, 20

MASS, ERR = {}, {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    try:
        ev = abs(float(line[52:60]))
    except ValueError:
        ev = 0.0
    for k in ids:
        if k:
            MASS[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0


def amps(a, t, q, mh, ml, kappa=0.0):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    sh = kappa * G0 * q / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) - sh for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    # 正典約定：κ 已在 diag 內（均勻位移），A 由本徵向量直接給出，**不得再加一次**
    A = ALPHA * q * (HBARC / a) * rA * 1000
    c = Cq * (HBARC / a) ** 3 * P0 * 1000
    return A, c, P0, bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)


print("=" * 100)
print("M5 docs/38 — 重區 ΔI=0 超精細 ＋ κ 兩約定統一")
print("=" * 100)

print("\n[一｜K1] 同電荷對之超精細（PDG 2024）")
# 同電荷對：de = M(star) − M(ground)（K1）
PAIRS = [("D", 421, 423), ("D_s", 431, 433), ("B", 511, 513), ("B_s", 531, 533)]
data = {}
for nm, i0, i1 in PAIRS:
    de = MASS[i1] - MASS[i0]   # star − ground
    s = math.hypot(ERR[i0], ERR[i1])
    data[nm] = {"Ehfs": de, "sigma": s, "c": de / 4, "sigma_c": s / 4}
    print(f"  {nm:4s} E_hfs = {de:9.4f} ± {s:7.4f} MeV ⟹ c = {de/4:9.4f} ± {s/4:.4f}")

print("\n[K2] 模型（同一 D 校準、零新參數）")
MODEL = [("D", 1.27, 0.30), ("D_s", 1.27, 0.495), ("B", 4.18, 0.30), ("B_s", 4.18, 0.495)]
mod = {}
for nm, mh, ml in MODEL:
    A, c, P0, bd = amps(A0, T0, 2 / 3, mh, ml)
    mod[nm] = {"c": c, "P0": P0, "bound": bd}
    print(f"  {nm:4s} P₀ = {P0:.6f}、c_model = {c:9.5f} MeV（束縛 {bd}）")

print("\n[K3] 比值判（模型 / 數據；門檻 1.5×）")
res = {}
for pair, hi, lo in (("D_s/D", "D_s", "D"), ("B_s/B", "B_s", "B")):
    rm = (mod[hi]["c"] / mod[lo]["c"]) * (0.30 / 0.495)   # m_l 因子（C_q ∝ 1/(m_h m_l)）
    rmod = mod[hi]["c"] / mod[lo]["c"]
    rd = data[hi]["c"] / data[lo]["c"]
    sd = rd * math.hypot(data[hi]["sigma_c"] / data[hi]["c"], data[lo]["sigma_c"] / data[lo]["c"])
    off = max(rmod / rd, rd / rmod)
    res[pair] = {"R_model": rmod, "R_data": rd, "sigma_R": sd, "off_factor": off,
                 "z": abs(rmod - rd) / sd, "pass": off <= 1.5}
    print(f"  {pair:7s} 模型 {rmod:.4f} vs 數據 {rd:.4f} ± {sd:.4f}"
          f"（{abs(rmod-rd)/sd:.1f}σ）⟹ 偏離 {off:.2f}× {'✅' if off <= 1.5 else '🔴'}")

print("\n[K4] 診斷：反解所需 m_l 指數 α（0 = 與 m_l 無關；−1 = 模型現行）")
for pair, hi, lo, ml_hi, ml_lo in (("D_s/D", "D_s", "D", 0.495, 0.30), ("B_s/B", "B_s", "B", 0.495, 0.30)):
    p0r = mod[hi]["P0"] / mod[lo]["P0"]
    rd = data[hi]["c"] / data[lo]["c"]
    alpha = math.log(rd / p0r) / math.log(ml_lo / ml_hi)
    res[pair]["alpha_needed"] = alpha
    print(f"  {pair:7s} P₀ 比 = {p0r:.4f}、數據比 = {rd:.4f} ⟹ 需 α = {alpha:+.3f}"
          f"（模型現行 α = −1；heavy-quark 極限 α = 0）")

print("\n[二｜K5] κ 正典約定（放入 Hamiltonian）")


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


def calib_kappa(kappa, ngrid=320):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), 2 / 3, 1.27, 0.30, kappa)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, 2 / 3, 1.27, 0.30, kappa)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))),
                                       2 / 3, 1.27, 0.30, kappa)[1] - C_D, prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


kappa_max = A_D / (G0 * (2 / 3))
print(f"  硬界：κ < A_D/(g(0)q_D) = {A_D}/({G0:.4f}×2/3) = {kappa_max:.4f}")
kscan = np.linspace(0.0, kappa_max * 0.999, 25)
best = None
rows = []
for kp in kscan:
    sol = calib_kappa(float(kp))
    if sol is None:
        continue
    a, t = sol
    AD, cD, pD, _ = amps(a, t, 2 / 3, 1.27, 0.30, kp)
    AK, cK, pK, _ = amps(a, t, 1 / 3, 0.495, 0.30, kp)
    R = (AK / cK) / (AD / cD) if cK else float("nan")
    rows.append({"kappa": float(kp), "a0_fm": a, "t_MeV": t * 1000, "R": float(R),
                 "A_D": AD, "c_D": cD, "bound_D": _})
    if best is None or (np.isfinite(R) and R > best["R"]):
        best = {"kappa": float(kp), "R": float(R), "a0_fm": a, "t_MeV": t * 1000}
print(f"  可校準點 {len(rows)}／{len(kscan)}；校準範圍 κ ∈ [{rows[0]['kappa']:.4f}, {rows[-1]['kappa']:.4f}]")
print(f"  R 範圍 [{min(r['R'] for r in rows if np.isfinite(r['R'])):.4f}, "
      f"{max(r['R'] for r in rows if np.isfinite(r['R'])):.4f}]")
print(f"  R_max（正典約定）= {best['R']:.4f} @ κ = {best['kappa']:.4f}")
print(f"  對照：docs/28 回報式之 R_max = 1.6470（κ = 0.2250，不同約定，不可混算）")
res_ok = all(v["pass"] for v in res.values())
out = {"data": data, "model": mod, "ratios": res, "all_pass_K3": bool(res_ok),
       "kappa_max": kappa_max, "kappa_scan": rows, "kappa_best": best,
       "docs28_reported_R_max": 1.6470}
with open("data/40_heavy_hfs_kappa.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/40_heavy_hfs_kappa.json")