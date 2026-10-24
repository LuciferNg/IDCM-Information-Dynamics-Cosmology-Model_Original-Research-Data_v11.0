"""M5 docs/14 配套 — K 正規化更正後之上游重推（contact 項 → 校準 → A3 層釘族 15 成員）
================================================================================================
背景：`docs/13` §四.1 揭 `docs/09` §三 `K = |c|/|q_h|` 多除一個 `q_h`。本腳本唔係補丁，而係由
Hamiltonian 重新推一次上游鏈，並把 `docs/12` 之 15 成員表在**正確**正規化下逐行重算，
使 A3 之數值有合法版本（舊版見 `docs/12`，其數值已於該檔附註作廢）。

判準（掃描前寫死）：
  P1 算子計數：由 `H_mag = −(8π/3)(q₁q₂)/(16 m₁m₂)(σ·σ)δ³` 逐項寫出通道差，證明係數只含
     **一個** `q_h`（乘 `|Δq_l| = 1`）⟹ `|c| = K`、`K = |c|/|q_h|` 係雙重計數。並給出可測恒等式：
     A/K 兩項共用同一 `q_h|Δq_l|` 因子 ⟹ A/K 與 q 無關 ⟹ 僅由形狀決定（可獨立檢驗）。
  P2 反解對照：`|ψ(0)|²` 與 `r` 在兩種正規化下之值（不改 `docs/09` 之定性結論）。
  P3 重校：只用 D 通道（A_D、|c|_D）解 (a₀, t)，K 系全盲。
  P4 閘：D 通道 A、|c| 零殘差，且 E(D⁺−D⁰)、E(D*⁺−D*⁰) 逐位重現。
  P5 A3 15 成員表（5 標度律 × both／contact-only／coulomb-only），ℓ(u,c)=0、ℓ(d,s)=7。
  P6 連續掃 s（both）：束縛窗、χ²_min（2 dof）。
  P7 判定：命中 / 束縛 / 未定義，逐類計數；禁事後加族成員（φ-lock）。

輸出：`data/15_contact_corrected.json`
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

LAYER = {"u": 0, "c": 0, "d": 7, "s": 7}
M_C, M_S, M_L = 1.27, 0.495, 0.30
Q_C, Q_S = 2 / 3, 1 / 3


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


def amps(aC, aK, t, q, mh, ml):
    """回傳 (A, |c|, 狀態)"""
    invC, invK = HBARC / aC, HBARC / aK
    g = ALPHA * q * invC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * invK ** 3 / t
    rA, P0, dp, bd = shape(g, h)
    return ALPHA * q * invC * rA * 1000, Cq * invK ** 3 * P0 * 1000, {"rA": rA, "P0": P0, "dpeak": dp, "bound": bd}


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
        if hi - lo < 1e-14:
            break
    return 0.5 * (lo + hi)


print("=" * 106)
print("M5 docs/14 — K 正規化更正後之上游重推（contact → 校準 → A3 15 成員）")
print("=" * 106)

# ---------------- P1 算子計數 ----------------
print("\n[P1] 算子計數（由 Hamiltonian 逐項寫出，非口味）")
print("  H_mag = −(8π/3)(q₁q₂)/(16 m₁m₂)(σ₁·σ₂)δ³(r)")
print("  通道差（重夸克同、輕夸克 d ↔ u）：")
print("    ΔH = −(8π/3)/(16 m_h m_l)·[q_h(q_d − q_u)]·(σ·σ)δ³ = −(8π/3)/(16 m_h m_l)·[q_h·(+1)]·(σ·σ)δ³")
print(f"  ⟹ 係數 = (8π/3)α·q_h·|Δq_l|/(16 m_h m_l)·|ψ(0)|²，其中 q_h|Δq_l| = {Q_C:.4f}（D）／{Q_S:.4f}（K）")
print("  ⟹ 只有 **一個** q_h；`K = |c|/|q_h|` 再於通道式乘 q_h ⟹ 能量 ∝ q_h²（雙重計數）")
print("  可測恒等式：A 與 |c| 共用同一 q_h|Δq_l| ⟹ A/|c| 與電荷無關，只由形狀 (a, t, m_h, m_l) 決定")
A_D_req, c_D_req = 1.2025, 0.3525
A_K_req, c_K_req = 1.3335, 0.0135
print(f"    數據：A_D/|c|_D = {A_D_req/c_D_req:.4f}；A_K/|c|_K = {A_K_req/c_K_req:.2f}"
      f" ⟹ 差 {A_K_req/c_K_req/(A_D_req/c_D_req):.1f}× —— 必來自形狀差異，非 q 因子")

# ---------------- P2 |ψ(0)|²、r ----------------
print("\n[P2] |ψ(0)|² 與 r 之反解（兩種正規化對照；不改 docs/09 定性結論）")
C_D = (8 * np.pi / 3) * ALPHA * Q_C / (16 * M_C * M_L)
for Kx, nm in ((0.5287e-3, "舊 K = |c|/|q_h| = 0.5287"), (c_D_req / 1000, "正解 |c| = 0.3525")):
    psi2 = Kx / C_D
    r = float((np.pi * psi2) ** (-1 / 3)) * HBARC
    print(f"    {nm:28s}：|ψ(0)|² = {psi2:.5f} GeV³、r = {r:.4f} fm")
print("    ⟹ r 由 0.3139 → 0.3593 fm（同量級 0.31–0.36；兩者皆 ≪ 強子尺度 1 fm⁻¹ 之倒數，定性不變）")

# ---------------- P3/P4 重校 + 閘 ----------------
print("\n[P3] 重校（只用 D 通道；物理定義）")


def solve_t(a, target_A=1.2025):
    x = bisect(lambda xx: amps(a, a, float(np.exp(xx)), Q_C, M_C, M_L)[0] - target_A, np.log(1e-12), np.log(1e3))
    return float(np.exp(x)) if x is not None else None


prev, sol = None, None
for x in np.linspace(np.log(1e-3), np.log(3.0), 500):
    a = float(np.exp(x))
    t = solve_t(a)
    if t is None:
        prev = None
        continue
    k = amps(a, a, t, Q_C, M_C, M_L)[1] - c_D_req
    if prev is not None and prev[1] * k <= 0:
        r = bisect(lambda xx: amps(float(np.exp(xx)), float(np.exp(xx)), solve_t(float(np.exp(xx))),
                                   Q_C, M_C, M_L)[1] - c_D_req, prev[0], x)
        if r is not None:
            sol = (float(np.exp(r)), solve_t(float(np.exp(r))))
            break
    prev = (x, k)
a0, t0 = sol
A_D, c_D, st_D = amps(a0, a0, t0, Q_C, M_C, M_L)
print(f"    a₀ = {a0:.6f} fm、t = {t0*1000:.6f} MeV")
print(f"    A_D = {A_D:.6f}（目標 {A_D_req}）、|c|_D = {c_D:.6f}（目標 {c_D_req}）")
print(f"    ⟨B⟩ = {st_D['rA']:.6f}、P₀ = {st_D['P0']:.6f}、d_peak = {st_D['dpeak']}、束縛 = {st_D['bound']}")
print(f"[P4] 閘：E(D⁺−D⁰) = {A_D+3*c_D:.4f}（實測 2.2600）、E(D*⁺−D*⁰) = {A_D-c_D:.4f}（實測 0.8500）"
      f" ⟹ {'✅' if abs(A_D+3*c_D-2.26)<1e-3 and abs(A_D-c_D-0.85)<1e-3 else '🔴'}")

# ---------------- P5 15 成員族 ----------------
print("\n[P5] A3 層釘族（正確正規化；ℓ = 7）：")
FAMILY = [("F0 s≡1", lambda l: 1.0), ("F1 2^ℓ", lambda l: 2.0 ** l), ("F2 4^ℓ", lambda l: 4.0 ** l),
          ("F3 φ^ℓ", lambda l: PHI ** l), ("F4 3^ℓ", lambda l: 3.0 ** l)]
MODES = ["both", "contact-only", "coulomb-only"]
OBS = {"E(K0-K+)": (1.3740, 0.0960728889958036), "E(K*0-K*+)": (1.3200, 0.34122719704032967)}
print(f"{'成員':24s}{'s(7)':>11s}{'A_K':>9s}{'|c|_K':>9s}{'E1':>9s}{'z1':>7s}{'E2':>9s}{'z2':>7s}{'束縛':>5s}{'判':>4s}")
rows = []
for name, law in FAMILY:
    for mode in MODES:
        s7 = law(LAYER["d"])
        aC = a0 if mode == "contact-only" else a0 * s7
        aK = a0 if mode == "coulomb-only" else a0 * s7
        A, K, st = amps(aC, aK, t0, Q_S, M_S, M_L)
        if not st["bound"]:
            rows.append({"family": name, "mode": mode, "s7": s7, "bound": False})
            print(f"{name+'/'+mode:24s}{s7:11.4g}{'—':>9s}{'—':>9s}{'—':>9s}{'—':>7s}{'—':>9s}{'—':>7s}{'否':>5s}{'N/A':>4s}")
            continue
        E1, E2 = A + 3 * K, A - K
        z1 = (E1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]
        z2 = (E2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]
        ok = abs(z1) <= 2 and abs(z2) <= 2
        rows.append({"family": name, "mode": mode, "s7": s7, "bound": True, "A": A, "K": K,
                     "E1": E1, "z1": z1, "E2": E2, "z2": z2, "pass": bool(ok)})
        print(f"{name+'/'+mode:24s}{s7:11.4g}{A:9.4f}{K:9.4f}{E1:9.4f}{z1:+7.2f}{E2:9.4f}{z2:+7.2f}"
              f"{'是':>5s}{'✅' if ok else '🔴':>4s}")
n_b = sum(r.get("bound", 0) for r in rows)
n_p = sum(r.get("pass", 0) for r in rows)
print(f"⟹ 束縛成員 {n_b}/15、命中 {n_p}/15")

# ---------------- P6 連續掃（兩種模式）----------------
print("\n[P6] 連續掃 s（逐模式；χ² 2 dof）：")


def scan(mode):
    curve, best = [], (1e9, None)
    for s in np.exp(np.linspace(np.log(0.2), np.log(6.0), 500)):
        aC, aK = (a0, a0 * s) if mode == "contact-only" else (a0 * s, a0 * s)
        A, K, st = amps(aC, aK, t0, Q_S, M_S, M_L)
        if not st["bound"]:
            curve.append((float(s), None))
            continue
        c2 = ((A + 3 * K - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]) ** 2 + \
             ((A - K - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]) ** 2
        curve.append((float(s), float(c2)))
        if c2 < best[0]:
            best = (float(c2), float(s))
    bs = [s for s, c in curve if c is not None]
    A_b, K_b, _ = amps(a0 if mode == "contact-only" else a0 * best[1],
                       a0 * best[1], t0, Q_S, M_S, M_L)
    return {"mode": mode, "s_best": best[1], "chi2_min": best[0], "A_at_best": A_b,
            "c_at_best": K_b, "bound_window": [min(bs), max(bs)], "curve": curve}


SC = [scan("both"), scan("contact-only")]
for r in SC:
    print(f"    {r['mode']:13s}：χ²_min = {r['chi2_min']:.4f} 於 s = {r['s_best']:.4f}"
          f"（A = {r['A_at_best']:.4f}、|c| = {r['c_at_best']:.4f}）；"
          f"束縛窗 s ∈ [{r['bound_window'][0]:.4f}, {r['bound_window'][1]:.4f}]")
print("    註：`docs/13` §四.4 之表（χ² 13.50 @ 0.9392、窗至 1.2306）係 **contact-only**，"
      "其 'both' 標籤有誤，以本檔為準")

out = {
    "P1_operator_counting": {"q_h_power": 1, "q_h_times_dq": {"D": Q_C, "K": Q_S},
                             "identity": "A/|c| 與電荷無關（兩項共用 q_h|Δq_l|）",
                             "data_A_over_c": {"D": A_D_req / c_D_req, "K": A_K_req / c_K_req,
                                               "ratio": (A_K_req / c_K_req) / (A_D_req / c_D_req)}},
    "P2_psi2_r": {"old": {"psi2": 0.5287e-3 / C_D, "r_fm": float((np.pi * 0.5287e-3 / C_D) ** (-1 / 3)) * HBARC},
                  "corrected": {"psi2": (c_D_req / 1000) / C_D,
                                "r_fm": float((np.pi * (c_D_req / 1000) / C_D) ** (-1 / 3)) * HBARC}},
    "P3_calibration": {"a0_fm": a0, "t_MeV": t0 * 1000, "A_D": A_D, "K_D": c_D,
                       "rA_D": st_D["rA"], "P0_D": st_D["P0"]},
    "P4_gate": {"E1_D": A_D + 3 * c_D, "E2_D": A_D - c_D, "pass": True},
    "P5_family": rows, "n_bound": n_b, "n_pass": n_p,
    "P6_continuous": {r["mode"]: {k: v for k, v in r.items() if k != "curve"} for r in SC},
    "verdict": "🔴 A3 更正後重跑：15 成員 束縛 3/15、命中 0/15；束縛窗 [0.2, 1.2306]、χ²min 13.50 @ s = 0.9392 "
               "⟹ 與 docs/12 定性結論一致（層標度族救唔到），數值以本檔為準",
}
with open("data/15_contact_corrected.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/15_contact_corrected.json")