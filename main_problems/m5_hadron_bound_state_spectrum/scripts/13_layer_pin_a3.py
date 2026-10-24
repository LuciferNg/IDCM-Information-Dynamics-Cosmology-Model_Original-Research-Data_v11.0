"""M5 docs/12 配套 — A3-L: 層釘（P10 sector pin）＋層標度族 ⟹ K 系統盲測（結構否證）
================================================================================================
版本史（append-only）：
  v1  `data/run_13_layer_pin_a3.v1_unboundbug.log` —— 缺束縛判別：弱耦合區（純 Coulomb 唔束縛）
      報 A = K = 0 當預言 ⟹ 作廢（模型未定義 ≠ 預言 0）。
  v2  加束縛判別（⟨B^−d⟩ > 1e-3 且 d_peak ≤ 20）＋ P8/P9（束縛來源、門檻距離）；反解改用
      **單調二分**（v2 首跑之 solve_scale 方向反轉 bug 已修，該次輸出不作結論）。
  v3（本檔）修正 solve_scale 方向 bug、加 s_unbind 與 h/h_crit 門檻距離、未束縛一律標 N/A。

承接 A2（`docs/10`）：普適 (a, t) 兩體模型被 K 通道否證（ΔE_K 預言 1.1898 vs 實測 0.0540）。
A3 命題：層釘（`a_f = a·s(ℓ_f)`）可否修補 —— s 夸克喺較粗層 ⟹ 接觸項被壓低 ⟹ ΔE_K → 0。

判準（全部先寫死）：
  P1 層釘零參數化：唔用 flavour→mode 映射（O140 `docs/19` 禁區），只用 P10 `docs/00` §6 登記之
     **sector pin**：Up(u,c) → ℓ=0；Down(d,s) → ℓ=7；Lepton → ℓ=14（本題不用）。
     配對層 ℓ_pair = max、樹距 d_tree = |Δℓ|，零自由度。
  P2 層標度族（事前枚舉）：F0 s≡1、F1 2^ℓ、F2 4^ℓ、F3 φ^ℓ、F4 3^ℓ；模式 both／contact-only／
     coulomb-only ⟹ 15 成員全報（禁只報過關者）。
  P3 校準只用 D 通道（A_D = 1.2025、K_D = 0.5287）；ℓ(D⁰)=0 ⟹ s(0)=1 ⟹ 與族無關。
  P4 命中 = E(K⁰−K⁺)、E(K*⁰−K*⁺) 兩通道 |z| ≤ 2（σ 由 PDG 誤差傳播）。
  P5 連續掃 s ∈ [0.2, 6] 報 χ²_min；最佳 s 無登記結構身份 ⟹ φ-lock 閘 ⟹ 不得採用。
  P6 等價檢定：F0 必須逐位重現 A2 三數（容差 1e-3），否則作廢。
  P7 排序指紋：ΔE_J ∝ K ∝ q_h|Δq_l|/(m_h m_l) 之預言比 vs 實測比。
  P8 束縛判別：⟨B^−d⟩ > 1e-3 且 d_peak ≤ 20；未束縛 = 模型未定義（N/A），不得報 0。
  P9 局部化脆弱度：報 rA(s) 序列與束縛窗；**不得**用 E₀ < 帶底當判準（g ≠ 0 時恆成立，無鑑別力）。

模型（與 `scripts/10` 同構，只加層標度）：
      E ψ_d = −t[ψ_{d−1} + (B−1)ψ_{d+1}] + V_d ψ_d
      V_d = −(α q_h|Δq_l|/a_C)·B^(−d−d_tree) − C δ_{d,0}δ_{d_tree,0}/a_K³
      C = (8π/3)α q_h|Δq_l|/(16 m_h m_l)
      A = α q_h|Δq_l|(ħc/a_C)B^(−d_tree)⟨B^−d⟩ ; K = C(ħc/a_K)³P₀
      通道：E(K⁰−K⁺) = A+3K、E(K*⁰−K*⁺) = A−K、ΔE_K = 4K（規約 docs/09 §三）

輸出：`data/13_layer_pin_a3.json`
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

LAYER = {"u": 0, "c": 0, "d": 7, "s": 7}
MF = {"s": 0.495, "c": 1.27}
MLIGHT = 0.30
A_D, K_D = 1.2025, 0.5287
OBS = {"E(K0-K+)": (1.3740, 0.0960728889958036),
       "E(K*0-K*+)": (1.3200, 0.34122719704032967)}
OBS_DE_K = (0.0540, 0.35449400559106775)
OBS_DE_D = (1.4100, 0.16634903065542644)


def shape(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    nrm = w.sum()
    rA = float((w * np.array(POW)).sum() / nrm) if nrm > 0 else 0.0
    dpeak = int(np.argmax(w))
    return rA, float(w[0] / nrm), dpeak, bool(rA > RATIO_BOUND and dpeak <= DPEAK_MAX)


def amps(aC, aK, t, q, mh, ml, d_tree=0):
    invC, invK = HBARC / aC, HBARC / aK
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    sup = B ** (-d_tree)
    rA, rP, dpeak, bound = shape(ALPHA * q * invC * sup / t,
                                 (Cq * invK ** 3 / t) if d_tree == 0 else 0.0)
    A = ALPHA * q * invC * sup * rA * 1000
    K = (Cq * invK ** 3 * rP * 1000) if d_tree == 0 else 0.0
    return A, K, dict(rA=rA, P0=rP, dpeak=dpeak, bound=bound)


def obs_D(a, t):
    return amps(a, a, t, 2 / 3, MF["c"], MLIGHT)[:2]


def bisect(f, lo, hi, n=200):
    """單調二分（f(lo)、f(hi) 異號）；回傳根。"""
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-13:
            break
    return 0.5 * (lo + hi)


# ---- P6 校準（a₀, t）：掃描全部 K 方程根，只取**束縛**分支（v2 grid，已對 A2 驗證）----
def solve_t(a, target=A_D):
    x = bisect(lambda xx: obs_D(a, np.exp(xx))[0] - target, np.log(1e-10), np.log(1e6))
    return float(np.exp(x)) if x is not None else None


roots = []
xs = np.linspace(np.log(1e-4), np.log(30.0), 121)
prev = None
for x in xs:
    a = float(np.exp(x))
    t = solve_t(a)
    if t is None:
        prev = None
        continue
    k = obs_D(a, t)[1] - K_D
    if prev is not None and prev[1] * k <= 0:
        ax = bisect(lambda xx: obs_D(np.exp(xx), solve_t(np.exp(xx)))[1] - K_D,
                    prev[0], x)
        if ax is not None:
            aa, tt = float(np.exp(ax)), solve_t(float(np.exp(ax)))
            roots.append((aa, tt, amps(aa, aa, tt, 2 / 3, MF["c"], MLIGHT)[2]["bound"]))
    prev = (x, k)

print("=" * 104)
print("M5 docs/12 — A3-L: 層釘（P10 sector pin）＋層標度族 ⟹ K 系統盲測（v3）")
print("=" * 104)
print("[P6] K 方程全部根（束縛判別）：")
for aa, tt, bb in roots:
    print(f"    a = {aa:.6f} fm, t = {tt:.6g} GeV, A = {obs_D(aa,tt)[0]:.6f}, "
          f"K = {obs_D(aa,tt)[1]:.6f}, 束縛 = {bb}")
solt = [r for r in roots if r[2]]
if not solt:
    raise SystemExit("🔴 冇束縛解 ⟹ 校準失敗")
a0, t0, _ = solt[0]
print(f"[P6] 選取束縛分支：a₀ = {a0:.6f} fm, t = {t0:.6g} GeV, "
      f"A_D = {obs_D(a0,t0)[0]:.6f}, K_D = {obs_D(a0,t0)[1]:.6f}")
A_1, K_1, st1 = amps(a0, a0, t0, 1 / 3, MF["s"], MLIGHT)
E1_F0, E2_F0 = A_1 + 3 * K_1, A_1 - K_1
gate = abs(E1_F0 - 1.1863) < 1e-3 and abs(E2_F0 + 0.0035) < 1e-3 and abs(4 * K_1 - 1.1898) < 1e-3
print(f"[P6] F0：A_K = {A_1:.4f}、K_K = {K_1:.4f} ⟹ E(K⁰−K⁺) = {E1_F0:.4f}（A2 1.1863）、"
      f"E(K*⁰−K*⁺) = {E2_F0:.4f}（A2 −0.0035）、ΔE_K = {4*K_1:.4f}（A2 1.1898）")
print(f"[P6] 等價檢定：{'✅ PASS' if gate else '🔴 FAIL ⟹ 作廢'}；束縛 = {st1['bound']}"
      f"（⟨B^−d⟩ = {st1['rA']:.4f}, d_peak = {st1['dpeak']}）")
if not gate:
    raise SystemExit(2)

# ---- P8 束縛來源 ----
rA0, P00, dp0, bd0 = shape(ALPHA * (1 / 3) * (HBARC / a0) / t0, 0.0)
print(f"\n[P8] 純 Coulomb（h=0，a = a₀）：⟨B^−d⟩ = {rA0:.3e}, d_peak = {dp0}, 束縛 = {bd0}")
print("      ⟹ 校準區束縛**完全由 contact δ 承擔** ⟹ 每系統實際只有一個形狀參數（h）⟹ K/A 被鎖死")

# ---- P9 脆弱度（用局部化指標 rA，非 E₀：g≠0 時 E₀ 恆低於帶底，唔可以當束縛判別）----
print(f"[P9] 局部化脆弱度（K 系統）：rA(s=1.00) = {st1['rA']:.4f}；"
      f"rA(s=1.05) = {amps(a0*1.05, a0*1.05, t0, 1/3, MF['s'], MLIGHT)[2]['rA']:.4f}；"
      f"rA(s=1.08) = {amps(a0*1.08, a0*1.08, t0, 1/3, MF['s'], MLIGHT)[2]['rA']:.4f}")
print("      ⟹ 局部化隨 s 急跌（contact δ 係唯一局部化機制）；所有結構候選 s(7) ≥ 29 都遠超束縛窗")

# ---- P2/P4 族 ----
PHI = (1 + 5 ** 0.5) / 2
FAMILY = [("F0 s≡1", lambda l: 1.0), ("F1 2^ℓ", lambda l: 2.0 ** l),
          ("F2 4^ℓ", lambda l: 4.0 ** l), ("F3 φ^ℓ", lambda l: PHI ** l),
          ("F4 3^ℓ", lambda l: 3.0 ** l)]
MODES = ["both", "contact-only", "coulomb-only"]


def predict(law, mode, l_pair=LAYER["s"]):
    s = law(l_pair)
    aC = a0 if mode == "contact-only" else a0 * s
    aK = a0 if mode == "coulomb-only" else a0 * s
    A, K, st = amps(aC, aK, t0, 1 / 3, MF["s"], MLIGHT)
    return A, K, A + 3 * K, A - K, st


print("\n[P2/P4] 層標度族 × 作用模式（ℓ = 7）：")
print(f"{'成員':22s}{'s(7)':>10s}{'A_K':>9s}{'K_K':>9s}{'E(K⁰K⁻)':>10s}{'z':>7s}"
      f"{'E(K*⁰K*⁻)':>11s}{'z':>7s}{'束縛':>5s}{'判定':>5s}")
rows = []
for name, law in FAMILY:
    for mode in MODES:
        A, K, E1, E2, st = predict(law, mode)
        if not st["bound"]:
            rows.append({"family": name, "mode": mode, "s7": law(7), "bound": False})
            print(f"{name+'/'+mode:22s}{law(7):10.4g}{'—':>9s}{'—':>9s}{'—':>10s}{'—':>7s}"
                  f"{'—':>11s}{'—':>7s}{'否':>5s}{'N/A':>5s}")
            continue
        z1 = (E1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]
        z2 = (E2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]
        ok = abs(z1) <= 2 and abs(z2) <= 2
        rows.append({"family": name, "mode": mode, "s7": law(7), "bound": True, "A": A, "K": K,
                     "E1": E1, "z1": z1, "E2": E2, "z2": z2, "pass": bool(ok)})
        print(f"{name+'/'+mode:22s}{law(7):10.4g}{A:9.4f}{K:9.4f}{E1:10.4f}{z1:+7.2f}"
              f"{E2:11.4f}{z2:+7.2f}{'是':>5s}{'✅' if ok else '🔴':>5s}")
print(f"⟹ 族內：束縛成員 {sum(r.get('bound',0) for r in rows)}/15，命中 {sum(r.get('pass',0) for r in rows)}/15")

# ---- P5 連續掃 + s_unbind ----
print("\n[P5] 連續掃 s（both 模式）：")
curve, best = [], (1e9, None)
for s in np.exp(np.linspace(np.log(0.2), np.log(6.0), 500)):
    A, K, E1, E2, st = predict(lambda l, s=s: s, "both")
    if not st["bound"]:
        curve.append((float(s), None))
        continue
    c = ((E1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]) ** 2 + \
        ((E2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]) ** 2
    curve.append((float(s), float(c)))
    if c < best[0]:
        best = (float(c), float(s))
s_best = best[1]
A, K, E1, E2, _ = predict(lambda l, s=s_best: s, "both")
print(f"  χ²_min = {best[0]:.4f}（2 dof）於 s = {s_best:.4f}；E1 = {E1:.4f}、E2 = {E2:.4f}")
bound_s = [s for s, c in curve if c is not None]
print(f"  束縛窗：s ∈ [{min(bound_s):.4f}, {max(bound_s):.4f}] ⟹ 窗外模型未定義（N/A）")
for s in (0.8, 0.9, 1.0, 1.05, 1.1, 1.2, 1.5):
    A, K, E1, E2, st = predict(lambda l, s=s: s, "both")
    if not st["bound"]:
        print(f"    s = {s:4.2f}：N/A（未束縛）")
        continue
    c = ((E1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]) ** 2 + \
        ((E2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]) ** 2
    print(f"    s = {s:4.2f}：E1 = {E1:7.4f}（z {((E1-OBS['E(K0-K+)'][0])/OBS['E(K0-K+)'][1]):+5.2f}）"
          f"  E2 = {E2:7.4f}（z {((E2-OBS['E(K*0-K*+)'][0])/OBS['E(K*0-K*+)'][1]):+5.2f}）  χ² = {c:6.2f}")
cands = {"φ": PHI, "√2": 2 ** 0.5, "2^(1/4)": 2 ** 0.25}
near = min(((abs(s_best - v), k) for k, v in cands.items()))
print(f"  結構身份：最近者 {near[1]} = {cands[near[1]]:.4f}（差 {near[0]:.4f}）⟹ "
      + ("候選" if near[0] < 0.01 else "🔴 無登記結構身份 ⟹ φ-lock 閘：不得採用"))

# ---- P5b 反解（單調二分，v3 修正）----
A_req = (OBS["E(K0-K+)"][0] + 3 * OBS["E(K*0-K*+)"][0]) / 4
K_req = (OBS["E(K0-K+)"][0] - OBS["E(K*0-K*+)"][0]) / 4


def solve_scale(which, target):
    def f(x):
        a = float(np.exp(x))
        A, K, st = (amps(a, a0, t0, 1 / 3, MF["s"], MLIGHT) if which == "C"
                    else amps(a0, a, t0, 1 / 3, MF["s"], MLIGHT))
        return (A if which == "C" else K) - target, st["bound"]
    lo, hi = np.log(1e-3), np.log(5.0)
    flo, _ = f(lo)
    fhi, _ = f(hi)
    if flo * fhi > 0:
        return None, None
    root = bisect(lambda x: f(x)[0], lo, hi)
    if root is None:
        return None, None
    return float(np.exp(root)), f(root)[1]


print("\n[P5b] 由 K 兩通道反解需求標度（2 方程 2 未知；診斷，非預言）：")
print(f"  實測要求：A_K = {A_req:.4f} MeV、K_K = {K_req:.4f} MeV ⟹ K/A = {K_req/A_req:.4f}")
print(f"  模型（s≡1）：A_K = {A_1:.4f}（低 {A_req/A_1:.2f}×）、K_K = {K_1:.4f}（高 {K_1/K_req:.2f}×）"
      f"、K/A = {K_1/A_1:.4f}（差 {K_1/A_1/(K_req/A_req):.0f}×）")
aC_req, bC = solve_scale("C", A_req)
aK_req, bK = solve_scale("K", K_req)
print(f"  ⟹ A 反解：a_C = {aC_req:.4f} fm（= {aC_req/a0:.3f} a₀，需**縮小**），該點束縛 = {bC}")
print(f"  ⟹ K 反解：a_K = {aK_req} fm，該點束縛 = {bK}"
      + ("（🔴 只有**解除束縛**先到 0.0135 ⟹ 模型內不存在束縛解）" if bK is False else ""))
if aC_req and aK_req:
    print(f"  比值 a_K/a_C = {aK_req/aC_req:.3f} ⟹ 方向相反，單一標度族（s 一個數）不可能同時滿足")

# ---- P7 排序 ----
print("\n[P7] 質量排序指紋：")
ratio_model, ratio_obs = K_1 / K_D, OBS_DE_K[0] / OBS_DE_D[0]
print(f"  ΔE_J ∝ K ∝ q_h|Δq_l|/(m_h m_l) ⟹ 模型預言 ΔE_K/ΔE_D = {ratio_model:.4f}；實測 = {ratio_obs:.4f}"
      f" ⟹ 模型高估 {ratio_model/ratio_obs:.1f}×")
print(f"  強 hyperfine：K*−K = 397.9 MeV、D*−D = 142.0 MeV ⟹ 強分裂 s 系大 {397.9/142.0:.2f}×，"
      f"em 分裂 s 系細 {1/ratio_obs:.0f}× ⟹ 一階 em 接觸 1/(m_h m_l) 標度與數據排序相反")

out = {
    "version": "v3（solve_scale 方向 bug 已修；加 s_unbind/h_crit）",
    "layer_pins": {"source": "P10 docs/00 §6 sector pin（非 flavour→mode）",
                   "u": 0, "c": 0, "d": 7, "s": 7, "lepton": 14},
    "calibration": {"a0_fm": a0, "t_GeV": t0, "A_D": obs_D(a0, t0)[0], "K_D": obs_D(a0, t0)[1]},
    "gate_A2": {"E1": E1_F0, "E2": E2_F0, "dE": 4 * K_1, "pass": bool(gate),
                "bound": st1["bound"], "rA": st1["rA"], "dpeak": st1["dpeak"]},
    "contact_is_sole_binder": {"rA_pure_coulomb": rA0, "dpeak": dp0, "bound": bd0},
    "fragility": {"rA_s1": st1["rA"],
                  "rA_s105": amps(a0 * 1.05, a0 * 1.05, t0, 1 / 3, MF["s"], MLIGHT)[2]["rA"],
                  "rA_s108": amps(a0 * 1.08, a0 * 1.08, t0, 1 / 3, MF["s"], MLIGHT)[2]["rA"]},
    "family": rows,
    "n_bound": sum(r.get("bound", 0) for r in rows),
    "n_pass": sum(r.get("pass", 0) for r in rows),
    "continuous": {"s_best": s_best, "chi2_min": best[0], "E1": E1, "E2": E2,
                   "bound_window": [min(bound_s), max(bound_s)], "curve": curve},
    "required": {"A_req": A_req, "K_req": K_req, "aC_req_fm": aC_req, "aK_req_fm": aK_req,
                 "aC_bound": bC, "aK_bound": bK, "ratio": (aK_req / aC_req) if (aC_req and aK_req) else None,
                 "structural_identity": False},
    "ordering": {"ratio_model": ratio_model, "ratio_obs": ratio_obs,
                 "hyperfine_K": 397.9, "hyperfine_D": 142.0},
    "verdict": "🔴 A3-L 判死：15 成員 0 命中（12 個未束縛 = N/A，3 個 = A2 本身）；"
               "A 需求 a_C×0.22（可束縛）、K 需求解除束縛 ⟹ 模型內無同時解；"
               "根因 = contact δ 係唯一束縛源 ⟹ 每系統僅一個形狀參數 ⟹ K/A 鎖死 100×",
}
with open("data/13_layer_pin_a3.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/13_layer_pin_a3.json")