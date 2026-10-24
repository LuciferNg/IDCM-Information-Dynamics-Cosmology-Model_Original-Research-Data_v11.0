"""M5 docs/18 配套 — 兩扇區第一驗：tree 1/r（B^{-d} 階梯）＋ chain 旁觀者（κ g(d)）
================================================================================================
觸發（`docs/17` 之自我更正）：`g(d)` **唔可以**當 Coulomb 勢用 —— `docs/01` §五.2 已證
`g` 係緊緻環上嘅 Poisson 核（`L·G = I − J/N`，行和 0、`g(0)/g(44) = −2` 精確、**完全冇 1/r**），
即 chain 單獨不可能承載 QED 自能差。1/r 只住在 **tree 徑向**（`docs/02`：`1/r ⟺ w = b = 4`）。
故 `docs/17` §五.3 之「g 作 Toeplitz 勢」出口**關閉**（屬已否證路線，不得重推）。

但 `docs/01` §四 C/D 另有可用之物：**chain 嘅非加性旁觀者項**
`E_chain = Σ_{i<j} κ g(d_ij) q_i q_j`，`κ = 0.26986 MeV`（由 O216 之 pairwise 跨度 5.664 MeV 定）。
本輪把它與 tree 扇區並置，構成**兩扇區**模型：

    A_total(通道) = A_tree + A_chain ,   A_chain = κ·g(d_pair)·q_h·|Δq_l|
    |c|_total = |c|_tree（chain 自旋盲 ⟹ 不貢獻 ΔE_J）

層指派（P10 `docs/00` §6）：ℓ(c) = ℓ(u) = 0、ℓ(d) = ℓ(s) = 7 ⟹
    d_pair(D) = |0 − 7| = 7 ⟹ g(7) = 1447/352；d_pair(K) = |7 − 7| = 0 ⟹ g(0) = 2581/352

判準（寫死）：
  S1 零新參數：κ 取 `docs/01` §D 之 0.26986 MeV（外部既定），(a_C, t) 由 D 通道兩數定。
  S2 命中 = K 兩通道 |z| ≤ 2，且 D 通道零殘差。
  S3 κ 掃描（敏感度）：κ ∈ [0, 0.35] 逐點重校，報 K 之 z1／z2 曲線與可行窗。
  S4 若 D 校準無解 ⟹ N/A（不得報 0）。
  S5 界限：chain 項只取單對（pair）配置；族／三體配置未推。

輸出：`data/19_two_sector.json`
"""
import json
from fractions import Fraction

import numpy as np

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
N = 88
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
KAPPA = 0.26986          # docs/01 §D（O216 pairwise 跨度 5.664 MeV ⟹ κ = 5.664/20.989）
LAY = {"c": 0, "u": 0, "d": 7, "s": 7}


def gx(d):
    return (Fraction(N * N - 1) - 6 * d * (N - d)) / (12 * N)


def shape(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    return float((w * np.array(POW)).sum()), float(w[0]), int(np.argmax(w)), bool(w[0] > 1e-3 and int(np.argmax(w)) <= 20)


def tree(aC, aK, t, q, mh, ml):
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aK) ** 3 / t
    rA, P0, dp, bd = shape(g, h)
    A = ALPHA * q * (HBARC / aC) * rA * 1000
    c = Cq * (HBARC / aK) ** 3 * P0 * 1000
    return A, c, {"rA": rA, "P0": P0, "dpeak": dp, "bound": bd}


def chain(kappa, q, d_pair):
    return kappa * float(gx(d_pair)) * q


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


A_D_T, C_D_T = 1.2025, 0.3525
CH_D = chain(KAPPA, Q_C, abs(LAY["c"] - LAY["d"]))
CH_K = chain(KAPPA, Q_S, abs(LAY["s"] - LAY["d"]))
print("=" * 104)
print("M5 docs/18 — 兩扇區第一驗（tree 1/r ＋ chain 旁觀者 κ g(d)）")
print("=" * 104)
print(f"\n[S0] chain 扇區（κ = {KAPPA} MeV，docs/01 §D）")
print(f"     g(0) = {gx(0)}、g(7) = {gx(7)}（Fraction 精確）")
print(f"     D：d_pair = |ℓ_c − ℓ_d| = {abs(LAY['c']-LAY['d'])} ⟹ A_chain = κ·g(7)·q_c = {CH_D:.6f} MeV")
print(f"     K：d_pair = |ℓ_s − ℓ_d| = {abs(LAY['s']-LAY['d'])} ⟹ A_chain = κ·g(0)·q_s = {CH_K:.6f} MeV")
print(f"     ⟹ D 需要 tree 補 {A_D_T - CH_D:+.6f} MeV（正值 = 可行；負值 = κ 過大 ⟹ 排除）")


def solve(kappa):
    """以 D 通道兩數（A_D 總量、|c|_D）解 (a_C, t)"""
    chd = chain(kappa, Q_C, 7)
    target = A_D_T - chd
    if target <= 0:
        return None

    def t_solve(a):
        x = bisect(lambda xx: tree(a, a, float(np.exp(xx)), Q_C, M_C, M_L)[0] - target, np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev, sol = None, None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 400):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = tree(a, a, t, Q_C, M_C, M_L)[1] - C_D_T
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: tree(float(np.exp(xx)), float(np.exp(xx)), t_solve(float(np.exp(xx))),
                                       Q_C, M_C, M_L)[1] - C_D_T, prev[0], x)
            if r is not None:
                sol = (float(np.exp(r)), t_solve(float(np.exp(r))))
                break
        prev = (x, k)
    return sol


print("\n[S1] κ 取 docs/01 值（零新參數）")
sol = solve(KAPPA)
if sol is None:
    print(f"     🔴 D 通道校準無解：tree 需補 {A_D_T - CH_D:+.6f} ≤ 0（chain 已超出 A_D）⟹ S4 記 N/A")
    out = {"kappa": KAPPA, "chain_D": CH_D, "chain_K": CH_K, "feasible": False,
           "verdict": "🔴 κ(docs/01) = 0.26986 過大：D 系 chain 一項已給 %.4f > A_D = 1.2025 ⟹ 兩扇區模型在零新參數下不可行" % CH_D}
else:
    a0, t0 = sol
    A_Dt, c_Dt, stD = tree(a0, a0, t0, Q_C, M_C, M_L)
    A_Kt, c_Kt, stK = tree(a0, a0, t0, Q_S, M_S, M_L)
    A_Ktot = A_Kt + CH_K
    E1, E2 = A_Ktot + 3 * c_Kt, A_Ktot - c_Kt
    z1 = (E1 - 1.3740) / 0.0960728889958036
    z2 = (E2 - 1.3200) / 0.34122719704032967
    print(f"     a₀ = {a0:.6f} fm、t = {t0*1000:.6f} MeV；tree A_D = {A_Dt:.6f}（總 {A_Dt+CH_D:.6f}）"
          f"、|c|_D = {c_Dt:.6f}")
    print(f"     K：tree A = {A_Kt:.4f} + chain {CH_K:.4f} = {A_Ktot:.4f}；|c|_K = {c_Kt:.4f}")
    print(f"        E(K⁰−K⁺) = {E1:.4f} vs 1.3740（z = {z1:+.2f}）；E(K*⁰−K*⁺) = {E2:.4f} vs 1.3200（z = {z2:+.2f}）"
          f" ⟹ {'✅ 命中' if abs(z1) <= 2 and abs(z2) <= 2 else '🔴 不中'}")
    out = {"kappa": KAPPA, "chain_D": CH_D, "chain_K": CH_K, "feasible": True, "a0_fm": a0, "t_MeV": t0 * 1000,
           "A_D_tree": A_Dt, "c_D": c_Dt, "A_K_tree": A_Kt, "A_K_total": A_Ktot, "c_K": c_Kt,
           "E1": E1, "z1": z1, "E2": E2, "z2": z2}

print("\n[S3] κ 掃描（每點重校 D；報 K 之 z）：")
print(f"     {'κ (MeV)':>9s}{'a₀ (fm)':>10s}{'t (MeV)':>10s}{'A_K 總':>9s}{'|c|_K':>8s}{'E1':>8s}{'z1':>7s}{'E2':>8s}{'z2':>7s}{'判':>4s}")
scan = []
for kp in np.linspace(0.0, 0.35, 36):
    s2 = solve(float(kp))
    if s2 is None:
        scan.append({"kappa": float(kp), "feasible": False})
        print(f"     {kp:9.4f}{'—':>10s}{'—':>10s}{'—':>9s}{'—':>8s}{'—':>8s}{'—':>7s}{'—':>8s}{'—':>7s}{'N/A':>4s}")
        continue
    a, t = s2
    A_Kt, c_Kt, _ = tree(a, a, t, Q_S, M_S, M_L)
    Atot = A_Kt + chain(float(kp), Q_S, 0)
    E1, E2 = Atot + 3 * c_Kt, Atot - c_Kt
    z1 = (E1 - 1.3740) / 0.0960728889958036
    z2 = (E2 - 1.3200) / 0.34122719704032967
    ok = abs(z1) <= 2 and abs(z2) <= 2
    scan.append({"kappa": float(kp), "feasible": True, "a0_fm": a, "t_MeV": t * 1000, "A_K_total": Atot,
                 "c_K": c_Kt, "E1": E1, "z1": z1, "E2": E2, "z2": z2, "pass": bool(ok)})
    print(f"     {kp:9.4f}{a:10.5f}{t*1000:10.4f}{Atot:9.4f}{c_Kt:8.4f}{E1:8.4f}{z1:+7.2f}{E2:8.4f}{z2:+7.2f}"
          f"{'✅' if ok else '🔴':>4s}")
pass_k = [s["kappa"] for s in scan if s.get("pass")]
print(f"     ⟹ 命中窗 κ ∈ {['%.4f' % min(pass_k), '%.4f' % max(pass_k)] if pass_k else '空'}；"
      f"docs/01 之 κ = {KAPPA} {'在窗內' if pass_k and min(pass_k) <= KAPPA <= max(pass_k) else '唔在窗內'}")


print("\n[S6] 結構下限：contact 係數之不可避免下限（與 Coulomb 扇區無關）")
print("     對 K 系固定 A_K（tree 部分）掃 (a_C, t)，求 min |c|：")
floor = {}
for A_tgt, tag in ((1.3335, "A_K 全部由 tree 提供"), (1.3335 - CH_K, "扣除 chain 0.6596 後之 tree 殘量")):
    best, arg = None, None
    for la in np.linspace(np.log(1e-3), np.log(3.0), 120):
        a = float(np.exp(la))
        for lt in np.linspace(np.log(1e-4), np.log(1e3), 160):
            t = float(np.exp(lt))
            A, c, st = tree(a, a, t, Q_S, M_S, M_L)
            if not st["bound"] or abs(A - A_tgt) > 5e-3:
                continue
            if best is None or c < best:
                best, arg = c, (a, t * 1000, A, st["P0"])
    floor[tag] = None if best is None else {"c_min": best, "a0_fm": arg[0], "t_MeV": arg[1], "A": arg[2], "P0": arg[3]}
    print(f"       目標 A = {A_tgt:.4f}（{tag}）⟹ min |c| = {best:.5f} MeV"
          f"（a₀ = {arg[0]:.5f} fm、t = {arg[1]:.4f} MeV、P₀ = {arg[3]:.2e}）"
          f" vs 觀測需要 |c|_K = 0.013500 ⟹ 差 {best/0.0135:.1f}×")
print("     ⟹ E1 − E2 ≡ 4|c| ⟹ 模型下限 4×min|c| 對比觀測 0.0540 MeV")
for k, v in floor.items():
    if v:
        print(f"       {k}: 4|c| ≥ {4*v['c_min']:.4f} MeV vs 觀測 0.0540 ⟹ 超 {4*v['c_min']/0.0540:.1f}×")
out["contact_floor"] = floor

out["scan"] = scan
out["pass_window"] = [min(pass_k), max(pass_k)] if pass_k else None
with open("data/19_two_sector.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/19_two_sector.json")