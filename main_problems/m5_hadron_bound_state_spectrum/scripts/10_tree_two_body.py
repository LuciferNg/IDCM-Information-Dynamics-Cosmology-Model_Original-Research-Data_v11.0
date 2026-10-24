"""P16 docs/10 配套（A2a）— 樹層兩體問題：(a, t) 兩未知、兩方程（A_D, K_D），再盲測 K 通道
================================================================================
⚠️ 版本史（誠實記錄，append-only）：
  v1（本日首跑，已封存 `scripts/deprecated/10_v1_minE_bug.py`，log = `data/run_10_tree_two_body.v1_buggy.log`）
     用 `min|E|` 揀本徵態 ⟹ 吸引勢下基態能量大負 ⟹ 誤揀**激發態** ⟹ 假「簡併射線」(a→0, t→∞)。
     該 v1 結論**作廢**。
  v2（本檔）改用**基態**（升序索引 0）＋「外層 a 二分、內層 t 二分」穩健求解；
     掃描用 float64（41×41 三對角，精度 ~1e-14，遠優於報告位數），最終解與 P₀ 用 mpmath 30 dps 覆核。

模型：Cayley 樹（B = 4，深度 L = 40）徑向方程
      E ψ_d = −t[ψ_{d−1} + (B−1)ψ_{d+1}] + V_d ψ_d，V_d = −(α q_h Δq_l/a)B^{−d} − C δ_{d,0}/a³
      C = (8π/3)α q_h Δq_l/(16 m_h m_l)
觀測：A = α q_h Δq_l (ħc/a)⟨B^{−d}⟩、K = C (ħc/a)³ P₀（P₀ = |ψ_0|²/Σ N_d|ψ_d|²）
目標：A_D = 1.2025 MeV、K_D = 0.5287 MeV（docs/09 之 D/D* 分解）
輸出：`data/10_tree_two_body.json`
"""
import json
import numpy as np
import mpmath as mp

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
N = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
A_D, K_D = 1.2025, 0.5287
qh, mh, ml = 2 / 3, 1.27, 0.30
C_CONTACT = (8 * np.pi / 3) * ALPHA * qh / (16 * mh * ml)     # GeV^-2
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))


def shape_np(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]                                   # 基態（升序第一）
    w = np.array(N) * v ** 2
    nrm = w.sum()
    return float((w * np.array(POW)).sum() / nrm), float(w[0] / nrm)


def obs(a_fm, t_GeV):
    inv_a = HBARC / a_fm
    rA, rP = shape_np(ALPHA * qh * inv_a / t_GeV, C_CONTACT * inv_a ** 3 / t_GeV)
    return ALPHA * qh * inv_a * rA * 1000, C_CONTACT * inv_a ** 3 * rP * 1000, rA, rP


def solve_t(a_fm, target=A_D):
    lo, hi = np.log(1e-10), np.log(1e6)
    f = lambda x: obs(a_fm, np.exp(x))[0] - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-14:
            break
    return float(np.exp(0.5 * (lo + hi)))


print("=" * 96)
print("P16 docs/10 — A2a: 樹層兩體問題求解 (a, t) ⟹ 盲測 K 通道")
print("=" * 96)

scan = []
for x in np.linspace(np.log(1e-4), np.log(30.0), 121):
    a = float(np.exp(x))
    t = solve_t(a)
    if t is None:
        continue
    A, K, rA, rP = obs(a, t)
    scan.append((a, t, A, K, rA, rP))
print(f"可解 a ∈ [{scan[0][0]:.4g}, {scan[-1][0]:.4g}] fm（{len(scan)} 點）")
print(f"  K 範圍 {min(r[3] for r in scan):.4f} … {max(r[3] for r in scan):.4f} MeV（目標 {K_D}）")

sol = None
for r1, r2 in zip(scan, scan[1:]):
    if (r1[3] - K_D) * (r2[3] - K_D) <= 0:
        lo, hi = np.log(r1[0]), np.log(r2[0])
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            t = solve_t(float(np.exp(mid)))
            k = obs(float(np.exp(mid)), t)[1] - K_D
            if (obs(float(np.exp(lo)), solve_t(float(np.exp(lo))))[1] - K_D) * k <= 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-14:
                break
        a_sol = float(np.exp(0.5 * (lo + hi)))
        sol = (a_sol, solve_t(a_sol))
        break

if sol is None:
    print("\n🔴 (A_D, K_D) 無解 ⟹ 兩體樹層模型不足以同時承載 Coulomb 與接觸")
    out = {"status": "no_solution",
           "K_range": [min(r[3] for r in scan), max(r[3] for r in scan)]}
else:
    a_sol, t_sol = sol
    A, K, rA, rP = obs(a_sol, t_sol)
    print(f"\n[D 通道解] a = {a_sol:.6f} fm   t = {t_sol:.6g} GeV")
    print(f"  A = {A:.6f} MeV（目標 {A_D}）   K = {K:.6f} MeV（目標 {K_D}）")
    print(f"  rA = ⟨B^−d⟩ = {rA:.6f}   P₀ = {rP:.6f}")

    mp.mp.dps = 30
    Nm = [mp.mpf(1)] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
    inv_a = HBARC / a_sol
    g = ALPHA * qh * inv_a / t_sol
    h = C_CONTACT * inv_a ** 3 / t_sol
    Hm = mp.zeros(L + 1, L + 1)
    for d in range(L + 1):
        Hm[d, d] = -mp.mpf(g) * mp.power(B, -d) - (mp.mpf(h) if d == 0 else 0)
    for d in range(L):
        Hm[d, d + 1] = -mp.sqrt(B) if d == 0 else -mp.sqrt(B - 1)
        Hm[d + 1, d] = Hm[d, d + 1]
    ev, evec = mp.eigsy(Hm)
    v = [evec[d, 0] for d in range(L + 1)]
    nrm = sum(Nm[d] * v[d] ** 2 for d in range(L + 1))
    rA_mp = sum(Nm[d] * v[d] ** 2 * mp.power(B, -d) for d in range(L + 1)) / nrm
    rP_mp = v[0] ** 2 / nrm
    print(f"  mpmath 30dps 覆核：rA = {mp.nstr(rA_mp, 12)}  P₀ = {mp.nstr(rP_mp, 12)}"
          f"  (|ΔP₀| = {mp.nstr(abs(rP_mp - mp.mpf(rP)), 3)})")

    ms = 0.495
    inv_aK = HBARC / a_sol
    rAK, rPK = shape_np(ALPHA * (1 / 3) * inv_aK / t_sol,
                        (8 * np.pi / 3) * ALPHA * (1 / 3) / (16 * ms * ml) * inv_aK ** 3 / t_sol)
    A_K = ALPHA * (1 / 3) * inv_aK * rAK * 1000
    K_K = (8 * np.pi / 3) * ALPHA * (1 / 3) / (16 * ms * ml) * inv_aK ** 3 * rPK * 1000
    dE_K = 4 * K_K
    z = (dE_K - 0.0540) / 0.3545
    print(f"\n[盲測 K 通道（同一 a、t，零新參數）]")
    print(f"  A_K = {A_K:.4f} MeV   K_K = {K_K:.4f} MeV ⟹ |ΔE_K| = {dE_K:.4f} MeV")
    print(f"  實測 ΔE_K = 0.0540 ± 0.3545 ⟹ z = {z:+.2f}")
    print(f"  預言 E(K⁰−K⁺) = A_K + 3K_K = {A_K + 3 * K_K:.4f}（實測 1.3740 ± 0.3545）")
    print(f"  預言 E(K*⁰−K*⁺) = A_K − K_K = {A_K - K_K:.4f}（實測 1.3200 ± 0.3545）")
    out = {"status": "solved", "a_fm": a_sol, "t_GeV": t_sol, "A_MeV": A, "K_MeV": K,
           "rA": rA, "P0": rP, "P0_mpmath": mp.nstr(rP_mp, 12),
           "K_channel": {"A": A_K, "K": K_K, "dE_pred": dE_K, "dE_obs": 0.0540,
                         "sigma_obs": 0.3545, "z": z}}

with open("data/10_tree_two_body.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/10_tree_two_body.json")