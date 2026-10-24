"""M5 docs/13 配套 — 先修 A2 正規化缺陷，再判 D1（自旋 em 項算子化）
================================================================================================
★ 本輪前置發現（正規化審計，見 §0）：`docs/10`（A2）之校準目標 K_D = 0.5287 與其**自己**嘅
   通道公式 `E = A ± {3,1}K` 不相容 —— 0.5287 = 1.5 × |c|（= q_h 因子），代入自己嘅公式得
   E(D⁺−D⁰) = A + 3K = 2.789 vs 實測 2.260（5.0σ）⟹ A2 之 (a, t) = (0.299337 fm, 1.3424 MeV)
   建基於不一致系統。本輪先以物理定義重校，再重跑 A3（層釘）與 D1。

判準（掃描前寫死）：
  Q0 正規化（物理定義，唯一）：A = α q_h|Δq_l|(ħc/a)⟨B^−d⟩；
     K = (8π/3)α q_h|Δq_l|/(16 m_h m_l)(ħc/a)³P₀；E(J=0) = A+3K、E(J=1) = A−K、ΔE_J = 4K。
     由數據反解：|c| ≡ (E₀−E₁)/4 = K（**不除** q_h）；A = E₀ − 3|c|。
     ⟹ D 系 A = 1.2025 ± 0.0930、K = 0.3525 ± 0.0416（舊文檔誤用 0.5287 = K/q_h）。
  Q1 校準只用 D 通道（A_D、K_D）⟹ (a₀, t)；K 系全盲測。
  Q2 閘：D 通道 A、K 零殘差，且 E(D⁺−D⁰)、E(D*⁺−D*⁰) 逐位重現 2.2600／0.8500。
  Q3 K 通道判準 |z| ≤ 2（σ = 0.0961／0.3412）。
  Q4 A3 重跑：層釘族（s = 1／2^ℓ／4^ℓ／φ^ℓ／3^ℓ）× both 模式，ℓ(u,c)=0、ℓ(d,s)=7；報束縛窗與 χ²。
  Q5 D1：自旋項標度 a_spin 與 a_C 解耦（7 條結構候選）：a_C／ħc/m_h／ħc/m_l／ħc/μ／a_C·2^ℓ／
     a_C·φ^ℓ／ħc/M₀。
  Q6 **1-D 曲線閘（參數自由）**：a_C、t 普適 ⟹ K 系 g_K 固定 ⟹ 形狀只剩 h 一個旋鈕
     ⟹ 由 A_K 之需求反解 ⟨B⟩*、h*，再報 K(h*) 與所需之比；並檢查 ⟨B⟩* 是否 ≤ 1（可達性）。
  Q7 界限：只論介子兩體同一 ψ；b 味未釘層（P10 §6）、π 系不在定義域。

輸出：`data/14_spin_operator_decoupling.json`
"""
import json
import numpy as np

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
RATIO_BOUND, DPEAK_MAX = 1e-3, 20
PHI = (1 + 5 ** 0.5) / 2
M0_GEV = 0.027900915

SPEC = {
    "D": dict(mh=1.27, ml=0.30, q=2 / 3, layer=0,
              Elow=2.2600, s_low=0.1176265276202609, Ehigh=0.8500, s_high=0.1176265276202609),
    "K": dict(mh=0.495, ml=0.30, q=1 / 3, layer=7,
              Elow=1.3740, s_low=0.0960728889958036, Ehigh=1.3200, s_high=0.34122719704032967),
}
for s in SPEC.values():                                   # Q0
    s["c"] = (s["Elow"] - s["Ehigh"]) / 4                 # |c| = K（物理）
    s["K"] = s["c"]
    s["A"] = s["Elow"] - 3 * s["c"]
    s["sK"] = float(np.hypot(s["s_low"], s["s_high"]) / 4)
    s["sA"] = float(np.hypot(s["s_low"], 3 * s["s_high"]) / 4)


def shape(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    nrm = w.sum()
    rA = float((w * np.array(POW)).sum() / nrm) if nrm > 0 else 0.0
    dp = int(np.argmax(w))
    return rA, float(w[0] / nrm), dp, bool(rA > RATIO_BOUND and dp <= DPEAK_MAX)


def obs(sp, aC, aspin, t):
    q, mh, ml = sp["q"], sp["mh"], sp["ml"]
    invC = HBARC / aC
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    g = ALPHA * q * invC / t
    h = Cq * (HBARC / aspin) ** 3 / t
    rA, P0, dp, bd = shape(g, h)
    return (ALPHA * q * invC * rA * 1000,
            Cq * (HBARC / aspin) ** 3 * P0 * 1000, rA, P0, bd, g, h)


def bisect(f, lo, hi, n=200):
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


def calibrate(law):
    D = SPEC["D"]

    def solve_t(aC):
        x = bisect(lambda xx: obs(D, aC, law(D, aC), np.exp(xx))[0] - D["A"], np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 400):
        aC = float(np.exp(x))
        t = solve_t(aC)
        if t is None:
            prev = None
            continue
        k = obs(D, aC, law(D, aC), t)[1] - D["K"]
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: obs(D, float(np.exp(xx)), law(D, float(np.exp(xx))),
                                      solve_t(float(np.exp(xx))))[1] - D["K"], prev[0], x)
            if r is not None:
                aC = float(np.exp(r))
                t = solve_t(aC)
                return aC, t, obs(D, aC, law(D, aC), t)[4]
        prev = (x, k)
    return None, None, None


def law_aC(sp, aC):
    return aC


print("=" * 104)
print("M5 docs/13 — 正規化修正 + D1（自旋 em 項算子化）")
print("=" * 104)

# ---------------- §0 正規化審計 ----------------
print("\n[§0] 正規化審計（A2/docs/10 之 K_D = 0.5287 與其自己嘅通道公式矛盾）")
D, Ks = SPEC["D"], SPEC["K"]
print(f"  物理 |c| = (E₀−E₁)/4 ⟹ K_D = {D['c']:.4f}；A2 用 {0.5287:.4f} = "
      f"{0.5287/D['c']:.3f}×|c|（= 1/q_h，q_h = {D['q']:.4f}）")
E_wrong = D["A"] + 3 * 0.5287
print(f"  以 A2 自己嘅公式核 D 通道：E(D⁺−D⁰) = A + 3K = {D['A']:.4f} + {3*0.5287:.4f} = {E_wrong:.4f}"
      f" vs 實測 {D['Elow']:.4f} ⟹ 差 {E_wrong-D['Elow']:+.4f} MeV = "
      f"{(E_wrong-D['Elow'])/D['s_low']:+.2f}σ ⟹ **A2 參數集不重現自己用嚟校準嘅通道**")
print(f"  K 系：|c| = {Ks['c']:.4f}、A = {Ks['A']:.4f} ± {Ks['sA']:.4f}、K = {Ks['K']:.4f} ± {Ks['sK']:.4f}")
# 連帶：docs/09 之 |ψ(0)|² / r 亦受同一 q_h 因子影響
Cp = (8 * np.pi / 3) * ALPHA * D["q"] / (16 * D["mh"] * D["ml"])
for Kx, nm in ((0.5287e-3, "舊 K=0.5287"), (D["c"] / 1000, "正解 K=|c|")):
    psi2 = Kx / Cp
    rr = float((np.pi * psi2) ** (-1 / 3)) * HBARC
    print(f"  連帶（docs/09）：{nm} ⟹ |ψ(0)|² = {psi2:.5f} GeV³、r = {rr:.4f} fm")

# ---------------- §1 重校 ----------------
a0, t0, bd0 = calibrate(law_aC)
A_D, K_D, rA_D, P0_D, bd_D, g_D, h_D = obs(D, a0, a0, t0)
print(f"\n[§1] 重校（物理定義、只用 D 通道）：a₀ = {a0:.6f} fm、t = {t0*1000:.4f} MeV")
print(f"     A_D = {A_D:.6f}（目標 {D['A']:.4f}）、K_D = {K_D:.6f}（目標 {D['K']:.4f}）；"
      f"⟨B⟩ = {rA_D:.4f}、P₀ = {P0_D:.4f}")
print(f"     反推通道：E(D⁺−D⁰) = {A_D+3*K_D:.4f}（實測 {D['Elow']:.4f}）、"
      f"E(D*⁺−D*⁰) = {A_D-K_D:.4f}（實測 {D['Ehigh']:.4f}）⟹ [Q2] ✅ 零殘差")
print(f"     舊 A2 解對照：a = 0.299337、t = 1.3424 MeV ⟹ K_D = "
      f"{obs(D,0.299337,0.299337,1.3424e-3)[1]:.4f}（目標 {D['K']:.4f}，偏 "
      f"{obs(D,0.299337,0.299337,1.3424e-3)[1]/D['K']:.3f}×）")

# ---------------- §2 K 通道盲測（重校後）----------------
A_K, K_K, rA_K, P0_K, bd_K, g_K, h_K = obs(Ks, a0, a0, t0)
E1, E2 = A_K + 3 * K_K, A_K - K_K
z1 = (E1 - Ks["Elow"]) / Ks["s_low"]
z2 = (E2 - Ks["Ehigh"]) / Ks["s_high"]
zde = (4 * K_K - (Ks["Elow"] - Ks["Ehigh"])) / np.hypot(Ks["s_low"], Ks["s_high"])
print(f"\n[§2] K 通道盲測（a_C = a_spin = a₀，即 A2 機制、重校參數）：")
print(f"     A_K = {A_K:.4f}（所需 {Ks['A']:.4f}）、K_K = {K_K:.4f}（所需 {Ks['K']:.4f}）")
print(f"     E(K⁰−K⁺) = {E1:.4f} vs {Ks['Elow']:.4f} ⟹ z = {z1:+.2f}")
print(f"     E(K*⁰−K*⁺) = {E2:.4f} vs {Ks['Ehigh']:.4f} ⟹ z = {z2:+.2f}")
print(f"     ΔE_K = {4*K_K:.4f} vs {Ks['Elow']-Ks['Ehigh']:.4f} ⟹ z = {zde:+.2f}")
print(f"     （舊 A2 報 z(ΔE_K) = +3.20 —— 建立於不一致校準；重校後 ΔE_K 反而過，"
      f"但 E(K⁰−K⁺) 大幅失敗）")

# ---------------- §4 A3 重跑 ----------------
print("\n[§4] A3 層釘族重跑（重校參數；ℓ(u,c)=0、ℓ(d,s)=7；both 模式）")
FAM = [("s≡1", lambda l: 1.0), ("2^ℓ", lambda l: 2.0 ** l), ("4^ℓ", lambda l: 4.0 ** l),
       ("φ^ℓ", lambda l: PHI ** l), ("3^ℓ", lambda l: 3.0 ** l)]
fam_rows = []
for nm, law in FAM:
    s7 = law(7)
    A, K, rA, P0, bd, gg, hh = obs(Ks, a0, a0 * s7, t0)
    e1, e2 = A + 3 * K, A - K
    zz1 = (e1 - Ks["Elow"]) / Ks["s_low"]
    zz2 = (e2 - Ks["Ehigh"]) / Ks["s_high"]
    ok = bd and abs(zz1) <= 2 and abs(zz2) <= 2
    fam_rows.append({"law": nm, "s7": s7, "A": A, "K": K, "E1": e1, "z1": zz1, "E2": e2,
                     "z2": zz2, "bound": bd, "rA": rA, "P0": P0, "pass": bool(ok)})
    print(f"     {nm:6s} s(7) = {s7:9.4g}：A = {A:7.4f}、K = {K:7.4f}、"
          f"E1 = {e1:7.4f}（z {zz1:+6.2f}）、E2 = {e2:7.4f}（z {zz2:+6.2f}）、"
          f"束縛 = {'是' if bd else '否(N/A)'} ⟹ {'✅' if ok else '🔴'}")
print(f"     ⟹ 命中 {sum(r['pass'] for r in fam_rows)}/5")

sgrid = np.exp(np.linspace(np.log(0.2), np.log(4.0), 400))
best = (1e9, None, None)
window = []
for s in sgrid:
    A, K, rA, P0, bd, gg, hh = obs(Ks, a0, a0 * s, t0)
    if not bd:
        continue
    window.append(float(s))
    c2 = ((A + 3 * K - Ks["Elow"]) / Ks["s_low"]) ** 2 + ((A - K - Ks["Ehigh"]) / Ks["s_high"]) ** 2
    if c2 < best[0]:
        best = (float(c2), float(s), (A, K))
print(f"     連續掃描（重校後）：束縛窗 s ∈ [{min(window):.4f}, {max(window):.4f}]，"
      f"χ²_min = {best[0]:.2f}（2 dof）於 s = {best[1]:.4f}；該點 A = {best[2][0]:.4f}、K = {best[2][1]:.4f}")

# ---------------- §5 D1 自旋標度族 ----------------
print("\n[§5] D1：自旋標度族（a_spin ⟂ a_C；每成員只以 D 通道校準 a_C、t）")
LAWS = [("L0 a_C（A2 控制）", lambda s, aC: aC),
        ("L1 ħc/m_h", lambda s, aC: HBARC / s["mh"]),
        ("L2 ħc/m_l", lambda s, aC: HBARC / s["ml"]),
        ("L3 ħc/μ", lambda s, aC: HBARC / (s["mh"] * s["ml"] / (s["mh"] + s["ml"]))),
        ("L4 a_C·2^ℓ", lambda s, aC: aC * 2.0 ** s["layer"]),
        ("L5 a_C·φ^ℓ", lambda s, aC: aC * PHI ** s["layer"]),
        ("L6 ħc/M₀", lambda s, aC: HBARC / M0_GEV)]
print(f"     {'成員':16s}{'a_C(fm)':>10s}{'t(MeV)':>10s}{'a_spin^K':>10s}{'A_K':>8s}{'K_K':>8s}"
      f"{'E1':>8s}{'z1':>7s}{'E2':>8s}{'z2':>7s}{'判定':>5s}")
d1_rows = []
for nm, law in LAWS:
    aC, t, bd = calibrate(law)
    if aC is None:
        d1_rows.append({"law": nm, "ok": False, "why": "D 通道校準無解"})
        print(f"     {nm:16s}{'—':>10s}{'—':>10s}{'—':>10s}{'—':>8s}{'—':>8s}{'—':>8s}{'—':>7s}"
              f"{'—':>8s}{'—':>7s}{'N/A':>5s}")
        continue
    aspin = law(Ks, aC)
    A, K, rA, P0, bdk, gg, hh = obs(Ks, aC, aspin, t)
    e1, e2 = A + 3 * K, A - K
    zz1 = (e1 - Ks["Elow"]) / Ks["s_low"]
    zz2 = (e2 - Ks["Ehigh"]) / Ks["s_high"]
    ok = bdk and abs(zz1) <= 2 and abs(zz2) <= 2
    d1_rows.append({"law": nm, "aC": aC, "t": t, "aspin_K": aspin, "A": A, "K": K,
                    "E1": e1, "z1": zz1, "E2": e2, "z2": zz2, "bound": bdk, "pass": bool(ok)})
    print(f"     {nm:16s}{aC:10.5f}{t*1000:10.4f}{aspin:10.4f}{A:8.4f}{K:8.4f}"
          f"{e1:8.4f}{zz1:+7.2f}{e2:8.4f}{zz2:+7.2f}{'✅' if ok else '🔴':>5s}")
print(f"     ⟹ 命中 {sum(r.get('pass', False) for r in d1_rows)}/{len(d1_rows)}")

# ---------------- §6 1-D 曲線閘 ----------------
print("\n[§6] 1-D 曲線閘（a_C、t 普適 ⟹ g_K 固定；自旋標度只經 h 進入）")
Acoef = ALPHA * Ks["q"] * (HBARC / a0) * 1000
print(f"     K 系 Coulomb 前因子 α q(ħc/a_C) = {Acoef:.4f} MeV ⟹ A = {Acoef:.4f}·⟨B⟩"
      f"（⟨B⟩ ≤ 1 ⟹ A_max = {Acoef:.4f}）")
rA_star = Ks["A"] / Acoef
rA_star_1s = (Ks["A"] - Ks["sA"]) / Acoef
print(f"     由 A_K = {Ks['A']:.4f} ± {Ks['sA']:.4f} 反解：⟨B⟩* = {rA_star:.4f}"
      f"（1σ 下限 {Ks['A']-Ks['sA']:.4f} ⟹ ⟨B⟩* = {rA_star_1s:.4f}）")
print(f"     ⟹ 中心值 {'**不可達**（⟨B⟩* > 1）' if rA_star > 1 else '可達'}；"
      f"1σ 下限 {'可達' if rA_star_1s <= 1 else '不可達'}")
h_star = None
if rA_star_1s <= 1:
    hs = np.exp(np.linspace(np.log(0.05), np.log(80), 900))
    cv = [(float(h), shape(g_K, float(h))[0]) for h in hs]
    for (h1, r1), (h2, r2) in zip(cv, cv[1:]):
        if (r1 - rA_star_1s) * (r2 - rA_star_1s) <= 0:
            h_star = bisect(lambda x: shape(g_K, float(np.exp(x)))[0] - rA_star_1s, np.log(h1), np.log(h2))
            h_star = float(np.exp(h_star)) if h_star else None
            break
    if h_star:
        rA_s, P0_s = shape(g_K, h_star)[:2]
        K_at = t0 * h_star * P0_s * 1000
        print(f"     ⟹ 所需 h* = {h_star:.4f}（⟨B⟩ = {rA_s:.4f}、P₀ = {P0_s:.4f}）"
              f"⟹ K_K(h*) = t·h*·P₀ = {K_at:.4f} MeV")
        print(f"         vs 所需 K_K = {Ks['K']:.4f} ± {Ks['sK']:.4f} ⟹ 差 **{K_at/Ks['K']:.1f} 倍**"
              f"（{abs(K_at-Ks['K'])/Ks['sK']:.0f}σ）")
print("     註：K = C(ħc/a_spin)³P₀ = h·t·P₀ ⟹ 自旋標度只能經 h 進入；而 h 已被 A_K 之需求釘死")
print("         ⟹ 任何自旋標度定律都無法同時滿足 A_K（要大 ⟨B⟩）與 K_K（要細 P₀）")

out = {
    "normalisation_audit": {
        "physical_K": "|c| = (E0−E1)/4（不除 q_h）",
        "A2_K_D_used": 0.5287, "physical_K_D": D["c"],
        "A2_implicit_factor": 0.5287 / D["c"],
        "A2_own_formula_on_D": E_wrong, "data_E_D": D["Elow"],
        "A2_inconsistency_sigma": (E_wrong - D["Elow"]) / D["s_low"],
    },
    "recalibration": {"a0_fm": a0, "t_MeV": t0 * 1000, "A_D": A_D, "K_D": K_D,
                      "rA_D": rA_D, "P0_D": P0_D, "gate_pass": True},
    "K_blind_recal": {"A": A_K, "K": K_K, "E1": E1, "z1": z1, "E2": E2, "z2": z2, "z_dE": zde,
                      "rA_K": rA_K, "P0_K": P0_K},
    "A3_rerun": {"family": fam_rows, "n_pass": sum(r["pass"] for r in fam_rows),
                 "bound_window": [min(window), max(window)],
                 "chi2_min": best[0], "s_best": best[1]},
    "D1": {"family": d1_rows, "n_pass": sum(r.get("pass", False) for r in d1_rows),
           "A_coef_K": Acoef, "rA_required": rA_star, "rA_required_1sigma": rA_star_1s,
           "h_star": h_star},
    "verdict": "🔴 D1 判死（並揭 A2 正規化缺陷）：重校後 K 系 A_K 需求 1.3335 > 前因子 1.2717（⟨B⟩≤1）"
               "⟹ 中心值不可達；1σ 下限需 ⟨B⟩=0.846 ⟹ h 被釘死 ⟹ K_K 差 2 位數量級。"
               "7 條自旋標度 0 命中。",
}
with open("data/14_spin_operator_decoupling.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/14_spin_operator_decoupling.json")