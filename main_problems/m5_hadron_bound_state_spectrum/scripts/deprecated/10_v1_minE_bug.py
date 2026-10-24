"""P16 docs/10 配套 — A2：樹層兩體本徵問題（精確對角化）⟹ 由 D/D* 解 (a, t)，再盲測 K 通道
================================================================================
模型（全部由 docs/02 核 ＋ docs/09 磁接觸構成，零外加參數）：
  徑向方程（Cayley 樹，分支 B = 4）：
      E ψ_d = −t[ψ_{d−1} + (B−1)ψ_{d+1}] + V_d ψ_d,   d ≥ 1
      E ψ_0 = −t·B·ψ_1 + V_0 ψ_0
  勢：V_d = −(α q_h Δq_l / a)·B^{−d}          （Coulomb，全部 d）
      V_0 另加接觸項 −C·δ_{d,0}/a³，C = (8π/3)α q_h Δq_l/(16 m_h m_l)
  對稱化：N_0 = 1, N_d = B(B−1)^{d−1}；φ_d = √N_d ψ_d ⟹ 對稱三對角
  觀測：rA = ⟨B^{−d}⟩，rP = P₀ = |ψ_0|²/Σ N_d|ψ_d|²
兩個未知 (a[fm], t[GeV])，兩條方程（D 通道）：
  A_D = α q_h Δq_l (ħc/a) rA       （Coulomb 部分，實測 1.2025 MeV）
  K_D = C (ħc/a)³ rP               （接觸幅度，實測 0.5287 MeV）
再用同一 (a, t) **盲測** K 通道（m_s = 0.495）⟹ ΔE_K 預言 vs 0.0540 ± 0.3545
輸出：`data/10_tree_two_body.json`
"""
import json
import mpmath as mp

mp.mp.dps = 30
HBARC = mp.mpf("0.1973269804")        # GeV·fm
ALPHA = 1 / mp.mpf("137.035999084")
B = 4
L = 40


def tree_ground(a_fm, t_GeV, qh, mh, ml, dq=1):
    """回傳 (rA, rP)：⟨B^{−d}⟩、P₀（基態）"""
    inv_a = HBARC / a_fm                       # GeV
    g = ALPHA * qh * dq * inv_a / t_GeV        # Coulomb（單位 t）
    h = (8 * mp.pi / 3) * ALPHA * qh * dq / (16 * mh * ml) * inv_a ** 3 / t_GeV
    N = [mp.mpf(1)] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
    H = mp.zeros(L + 1, L + 1)
    for d in range(L + 1):
        Vd = -g * mp.power(B, -d)
        if d == 0:
            Vd -= h
        H[d, d] = Vd
    off1 = -mp.sqrt(B)
    off = -mp.sqrt(B - 1)
    for d in range(L):
        H[d, d + 1] = off1 if d == 0 else off
        H[d + 1, d] = H[d, d + 1]
    ev, evec = mp.eigsy(H)
    # ⚠ BUG FIX (2026-10-24): 原用 min|E| 會揀到**激發態**（吸引勢下基態能量大負）
    # eigsy 回升序 ⟹ 基態 = 索引 0
    idx = 0
    v = [evec[d, idx] for d in range(L + 1)]
    norm = sum(N[d] * v[d] ** 2 for d in range(L + 1))
    rP = v[0] ** 2 / norm
    rA = sum(N[d] * v[d] ** 2 * mp.power(B, -d) for d in range(L + 1)) / norm
    return rA, rP


def channel(a_fm, t_GeV, qh, mh, ml, dq=1):
    rA, rP = tree_ground(a_fm, t_GeV, qh, mh, ml, dq)
    inv_a = HBARC / a_fm
    A = ALPHA * qh * dq * inv_a * rA * 1000                      # MeV
    K = (8 * mp.pi / 3) * ALPHA * qh * dq / (16 * mh * ml) * inv_a ** 3 * rP * 1000  # MeV
    return A, K, rA, rP


# ---- 目標（docs/08/09 實測，D 通道）----
A_D, K_D = mp.mpf("1.2025"), mp.mpf("0.5287")
qh_D, mh_D, ml_D = mp.mpf(2) / 3, mp.mpf("1.27"), mp.mpf("0.30")

print("=" * 96)
print("P16 docs/10 — A2: 樹層兩體本徵問題 ⟹ 解 (a, t) ⟹ 盲測 K 通道")
print("=" * 96)


def resid(x, y):
    a_fm, t_GeV = mp.e ** x, mp.e ** y          # 正參數化
    A, K, _, _ = channel(a_fm, t_GeV, qh_D, mh_D, ml_D)
    return [mp.log(A / A_D), mp.log(K / K_D)]


sol = mp.findroot(resid, [mp.log("0.3"), mp.log("0.5")], tol=mp.mpf("1e-20"), maxsteps=200)
a_fm, t_GeV = mp.e ** sol[0], mp.e ** sol[1]
print(f"\n[D 通道解]  a = {mp.nstr(a_fm, 8)} fm   t = {mp.nstr(t_GeV, 8)} GeV")
A, K, rA, rP = channel(a_fm, t_GeV, qh_D, mh_D, ml_D)
print(f"  重現：A = {mp.nstr(A, 8)} MeV (目標 {A_D})   K = {mp.nstr(K, 8)} MeV (目標 {K_D})")
print(f"  ⟹ rA = {mp.nstr(rA, 8)}   P₀ = {mp.nstr(rP, 8)}   P₀/rA = {mp.nstr(rP/rA, 8)}")

# ---- 盲測：K 通道，同一 (a, t) ----
print("\n[盲測 K 通道（同一 a、t，零新參數）]")
AK, KK, rAK, rPK = channel(a_fm, t_GeV, mp.mpf(1) / 3, mp.mpf("0.495"), ml_D)
dE_K = 4 * KK                                    # |ΔE| = 4|K|（⟨σ·σ⟩: −3 → +1）
print(f"  A_K = {mp.nstr(AK, 8)} MeV   K_K = {mp.nstr(KK, 8)} MeV ⟹ |ΔE_K| 預言 = {mp.nstr(dE_K, 8)} MeV")
print(f"  實測 ΔE_K = 0.0540 ± 0.3545 MeV ⟹ z = {mp.nstr((dE_K - mp.mpf('0.0540')) / mp.mpf('0.3545'), 6)}")

# ---- 一致性：預言 E(K)、E(K*) ----
print(f"  ⟹ 預言 E(K⁰−K⁺) = A_K + 3K_K = {mp.nstr(AK + 3 * KK, 6)} MeV（實測 1.374 ± 0.35）")
print(f"  ⟹ 預言 E(K*⁰−K*⁺) = A_K − K_K = {mp.nstr(AK - KK, 6)} MeV（實測 1.320 ± 0.35）")

out = {"a_fm": mp.nstr(a_fm, 10), "t_GeV": mp.nstr(t_GeV, 10),
       "D_channel": {"A": mp.nstr(A, 10), "K": mp.nstr(K, 10),
                     "rA": mp.nstr(rA, 10), "P0": mp.nstr(rP, 10)},
       "K_blind": {"A": mp.nstr(AK, 10), "K": mp.nstr(KK, 10),
                   "dE_pred": mp.nstr(dE_K, 10), "dE_obs": 0.0540, "sigma_obs": 0.3545,
                   "z": mp.nstr((dE_K - mp.mpf('0.0540')) / mp.mpf('0.3545'), 6)},
       "note": "兩體樹層問題，B=4，L=40，mpmath 30 dps；未知 (a,t) 兩條方程，K 通道為真盲測"}
with open("data/10_tree_two_body.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/10_tree_two_body.json")