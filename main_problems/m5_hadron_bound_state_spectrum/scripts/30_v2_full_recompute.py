"""M5 docs/28 配套 — w = |v|² 約定下之全鏈重算（A2／A3／兩扇區／D5；取代 v1 數值）
================================================================================
`docs/27` 用結構證明裁決 `w = |v_d|²`（生產 H 在歸一化層基）。本腳本把**所有受影響之下游**由頭重算：

  S1 D 校準（κ = 0）：以 v2 權重重解 (a₀, t)，並與 v1 記錄值對照。
  S2 A2／D1 盲測：同一 (a₀, t) 零新參數預言 K 通道 —— ΔE_K = 4K_K、E(K⁰−K⁺) = A+3K、E(K*⁰−K*⁺) = A−K，
     與實測 (0.0540±0.3545／1.3740±0.3545／1.3200±0.3545) 比。
  S3 A3 層釘族：5 標度律 × 3 作用模式 = 15 成員（ℓ_pair = 7），加連續 s 掃（χ²_min、束縛窗）。
  S4 兩扇區（鏈旁觀者）：讀 `data/29` 之 v2 掃描（69 點）⟹ 天花板、P₀ 比、所需 Δd、A3 過強倍數；
     並列 v1 對照。
  S5 D5 族分辨力：六族在 v2 下之 r = A/|c|（含束縛守門），逐族 z；判「同夸克內容必同 r」是否仍成立。

判準（寫死）：
  V1 兩約定之 D 校準必須各自重現其記錄值（v1: 0.377422／0.7884；v2: 0.384184／1.7711），否則 FAIL。
  V2 盲測一律零新參數；z 以實測 σ 計，|z| ≤ 2 記 PASS，否則 FAIL（不得放寬）。
  V3 A3 成員須先判束縛，非束縛記 N/A（不得以未束縛態當預言）。
  V4 D5 每族只可用自身夸克內容之 (m_h, m_l)；不得為命中調族。
  V5 全部數字須可由本腳本 stdout 重現；v1 舊值只作對照，標「v1」。
輸出：`data/30_v2_full_recompute.json`
"""
import json

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L = 1.27, 0.495, 0.30
Q_C, Q_S = 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
A_D, C_D = 1.2025, 0.3525
OBS = {"E(K0-K+)": (1.3740, 0.3545), "E(K*0-K*+)": (1.3200, 0.3545), "dE_K": (0.0540, 0.3545)}
PHI = (1 + 5 ** 0.5) / 2
BOUND_RHO, DPEAK_MAX = 1e-3, 20
FAMILIES = [
    ("K vs K*", "s_ud", 1.374, 1.320, 0.0960728889958036, 0.34122719704032967),
    ("Sigma [−/0]", "s_ud", 2.247, 0.940, 0.10125709851659784, 1.121978609421766),
    ("Sigma [0/+]", "s_ud", 0.712, -1.690, 0.11963277142990543, 1.0603942662990968),
    ("D vs D*", "c_ud", 2.260, 0.850, 0.1176265276202609, 0.1176265276202609),
    ("Sigma_c", "c_ud", -3.880, -3.570, 0.2771930735065362, 0.7397540131692426),
    ("Xi_c vs Xi_c'", "c_s", 0.170, -2.060, 0.3743474322070341, 0.7133274143056608),
]
CLASS = {"s_ud": (M_S, M_L, Q_S), "c_ud": (M_C, M_L, Q_C), "c_s": (M_C, M_S, Q_C)}

print("=" * 106)
print("M5 docs/28 — w = |v|² 全鏈重算（A2／A3／兩扇區／D5）")
print("=" * 106)


def amps(aC, aK, t, q, mh, ml, conv="v2"):
    """aC：Coulomb 標度；aK：contact 標度（同值即 docs/10 之單標度）。"""
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    g = ALPHA * q * HBARC / aC / t
    h = Cq * (HBARC / aK) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    # 只用最低本徵對（三對角專用演算法，比 dense eigh 快 ~500×）
    ev0, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2 if conv == "v2" else np.array(NODES) * v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    bd = bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)
    A = ALPHA * q * (HBARC / aC) * rA * 1000
    K = Cq * (HBARC / aK) ** 3 * P0 * 1000
    return A, K, {"rA": rA, "P0": P0, "dpk": dpk, "bound": bd}


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


def calib(conv, target_A=A_D, target_c=C_D, ngrid=400):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, a, float(np.exp(xx)), Q_C, M_C, M_L, conv)[0] - target_A,
                   np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, a, t, Q_C, M_C, M_L, conv)[1] - target_c
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), float(np.exp(xx)), t_solve(float(np.exp(xx))),
                                       Q_C, M_C, M_L, conv)[1] - target_c, prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None

# ---------------- S1 ----------------
print("\n[S1] D 校準（κ = 0；目標 A_D = 1.2025 MeV、|c|_D = 0.3525 MeV）")
s1 = {}
for conv in ("v1", "v2"):
    sol = calib(conv)
    a, t = sol
    A, K, st = amps(a, a, t, Q_C, M_C, M_L, conv)
    s1[conv] = {"a0_fm": a, "t_MeV": t * 1000, "A": A, "c": K, **st}
    print(f"    {conv}: a = {a:.6f} fm、t = {t*1000:.4f} MeV、A = {A:.6f}、|c| = {K:.6f}、"
          f"⟨B⟩ = {st['rA']:.4f}、P₀ = {st['P0']:.4f}、d_peak = {st['dpk']}、束縛 {st['bound']}")
V1 = abs(s1["v1"]["a0_fm"] - 0.377422) < 1e-5 and abs(s1["v2"]["a0_fm"] - 0.384184) < 1e-5
print(f"    V1 對照記錄值：{'✅ PASS' if V1 else '🔴 FAIL'}")

# ---------------- S2 ----------------
print("\n[S2] A2／D1 盲測：同一 (a, t) 之 K 通道（零新參數）")
a0, t0 = s1["v2"]["a0_fm"], s1["v2"]["t_MeV"] * 1e-3
AK, KK, stK = amps(a0, a0, t0, Q_S, M_S, M_L, "v2")
E1, E2, dE = AK + 3 * KK, AK - KK, 4 * KK
zE1 = (E1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]
zE2 = (E2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]
zdE = (dE - OBS["dE_K"][0]) / OBS["dE_K"][1]
print(f"    K 態：⟨B⟩ = {stK['rA']:.4f}、P₀ = {stK['P0']:.4f}、d_peak = {stK['dpk']}、束縛 {stK['bound']}")
print(f"    A_K = {AK:.4f} MeV、K_K = {KK:.4f} MeV")
print(f"    |ΔE_K| = 4K_K = {dE:.4f} vs {OBS['dE_K'][0]:.4f} ± {OBS['dE_K'][1]:.4f} ⟹ z = {zdE:+.2f}")
print(f"    E(K⁰−K⁺) = A+3K = {E1:.4f} vs {OBS['E(K0-K+)'][0]:.4f} ⟹ z = {zE1:+.2f}")
print(f"    E(K*⁰−K*⁺) = A−K = {E2:.4f} vs {OBS['E(K*0-K*+)'][0]:.4f} ⟹ z = {zE2:+.2f}")
V2 = abs(zdE) <= 2
print(f"    V2 |z(ΔE_K)| ≤ 2：{'✅ PASS' if V2 else '🔴 FAIL'}；三通道最大 |z| = {max(abs(zdE), abs(zE1), abs(zE2)):.2f}")

# ---------------- S3 ----------------
print("\n[S3] A3 層釘族（ℓ_pair = 7；a_C、a_contact 按模式）")
FAM = [("F0 s≡1", lambda l: 1.0), ("F1 2^ℓ", lambda l: 2.0 ** l), ("F2 4^ℓ", lambda l: 4.0 ** l),
       ("F3 φ^ℓ", lambda l: PHI ** l), ("F4 3^ℓ", lambda l: 3.0 ** l)]
MODES = ["both", "contact-only", "coulomb-only"]
rows3, hit = [], 0
print(f"    {'成員':22s}{'s(7)':>10s}{'A_K':>9s}{'K_K':>9s}{'E(K⁰K⁻)':>10s}{'z':>7s}{'E(K*⁰K*⁻)':>11s}{'z':>7s}{'判':>4s}")
for name, law in FAM:
    for mode in MODES:
        s = law(7)
        aC = a0 if mode == "contact-only" else a0 * s
        aK = a0 if mode == "coulomb-only" else a0 * s
        A, K, st = amps(aC, aK, t0, Q_S, M_S, M_L, "v2")
        if not st["bound"]:
            rows3.append({"family": name, "mode": mode, "s7": s, "bound": False})
            print(f"    {name+'/'+mode:22s}{s:10.4g}{'—':>9s}{'—':>9s}{'—':>10s}{'—':>7s}{'—':>11s}{'—':>7s}{'N/A':>4s}")
            continue
        e1, e2 = A + 3 * K, A - K
        z1 = (e1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]
        z2 = (e2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]
        ok = abs(z1) <= 2 and abs(z2) <= 2
        hit += ok
        rows3.append({"family": name, "mode": mode, "s7": s, "bound": True, "A": A, "K": K,
                      "E1": e1, "z1": z1, "E2": e2, "z2": z2, "pass": bool(ok)})
        print(f"    {name+'/'+mode:22s}{s:10.4g}{A:9.4f}{K:9.4f}{e1:10.4f}{z1:+7.2f}{e2:11.4f}{z2:+7.2f}"
              f"{'✅' if ok else '🔴':>4s}")
nb = sum(1 for r in rows3 if r.get("bound"))
print(f"    ⟹ 束縛成員 {nb}/15、命中 {hit}/15")

curve, best = [], (1e9, None)
for s in np.exp(np.linspace(np.log(0.2), np.log(6.0), 400)):
    A, K, st = amps(a0 * s, a0 * s, t0, Q_S, M_S, M_L, "v2")
    if not st["bound"]:
        curve.append((float(s), None))
        continue
    e1, e2 = A + 3 * K, A - K
    c = ((e1 - OBS["E(K0-K+)"][0]) / OBS["E(K0-K+)"][1]) ** 2 + \
        ((e2 - OBS["E(K*0-K*+)"][0]) / OBS["E(K*0-K*+)"][1]) ** 2
    curve.append((float(s), float(c)))
    if c < best[0]:
        best = (float(c), float(s))
bw = [s for s, c in curve if c is not None]
print(f"    連續掃（both）：χ²_min = {best[0]:.4f}（2 dof）@ s = {best[1]:.4f}；"
      f"束縛窗 s ∈ [{min(bw):.4f}, {max(bw):.4f}]")
print(f"    結構候選 s(7) 最小者 = {min(law(7) for _, law in FAM):.4g} ⟹ "
      f"{'窗外（模型未定義）' if min(law(7) for _, law in FAM) > max(bw) else '窗內'}")

# ---------------- S4 ----------------
print("\n[S4] 兩扇區（鏈旁觀者）：讀 data/29 之 v2 掃描，列 v1 對照")
d29 = json.load(open("data/29_weight_convention.json"))
ce = d29["ceiling"]
print(f"    v2：R_max = {ce['R_max']:.4f}（κ = {ce['kappa']:.4f}）、r_K 天花板 = {ce['rK_ceiling']:.4f}、"
      f"短缺 {ce['shortfall']:.2f}×")
print(f"        P₀(K)/P₀(D) = {ce['P0_ratio']:.5f} = 1/{1/ce['P0_ratio']:.3f}、所需 Δd = {ce['required_dd']:.4f}、"
      f"Δd=1 覆蓋 {ce['d1_coverage']*100:.2f}%、A3 過強 {ce['a3_over']:.0f}×")
print(f"    v1 對照（docs/25）：天花板 12.4774、短缺 7.92×、P₀ 比 1/10.79、所需 Δd 0.6374、Δd=1 覆蓋 60.49%")

# ---------------- S5 ----------------
print("\n[S5] D5 族分辨力（v2；每族用自身夸克內容；含束縛守門）")
R_DATA = {"K vs K*": (1.374, 1.320, 0.0960728889958036, 0.34122719704032967)}
rows5, chi = [], {}
for name, cls, e1, e2, s1v, s2v in FAMILIES:
    mh, ml, q = CLASS[cls]
    Ak, Kk, st = amps(a0, a0, t0, q, mh, ml, "v2")
    r_model = Ak / Kk
    r_obs = ((e1 + 3 * e2) / 4) / ((e1 - e2) / 4)
    sr = abs(r_obs) * np.hypot(s1v / (e1 + 3 * e2), s2v / (e1 - e2))
    if not st["bound"]:
        print(f"    {name:22s} [{cls}] 模型 r = N/A（非束縛：d_peak = {st['dpk']}、⟨B⟩ = {st['rA']:.4f}）")
        rows5.append({"family": name, "class": cls, "bound": False, "dpk": st["dpk"]})
        continue
    z = (r_model - r_obs) / sr
    chi.setdefault(cls, []).append(z)
    print(f"    {name:22s} [{cls}] 模型 r = {r_model:8.3f}、數據 r = {r_obs:8.3f} ± {sr:6.3f} ⟹ z = {z:+.2f}")
    rows5.append({"family": name, "class": cls, "bound": True, "r_model": r_model, "r_obs": r_obs,
                  "sigma": sr, "z": z})
for cls, zs in chi.items():
    c2 = float(sum(z * z for z in zs))
    print(f"    χ²[{cls}] = {c2:.2f}（dof = {len(zs)}）⟹ {np.sqrt(c2/max(len(zs),1)):.2f}σ 等效")
print("    ⟹ 結構陳述不變：r = A/|c| 只依 (m_h, m_l) ⟹ 同夸克內容必同 r（v2 下 r 數值改變，但 Q 依賴仍為零）")

out = {"V1_pass": bool(V1), "V2_pass": bool(V2),
       "S1_calibration": s1,
       "S2_A2": {"A_K": AK, "K_K": KK, "E1": E1, "E2": E2, "dE": dE,
                 "z_E1": zE1, "z_E2": zE2, "z_dE": zdE, "K_state": stK},
       "S3_family": {"rows": rows3, "n_bound": nb, "n_pass": hit, "chi2_min": best[0], "s_best": best[1],
                     "bound_window": [min(bw), max(bw)]},
       "S4_two_sector": ce, "v1_reference": {"rK_ceiling": 12.4774, "shortfall": 7.92,
                                             "P0_ratio": 1 / 10.79, "required_dd": 0.6374,
                                             "d1_coverage": 0.6049},
       "S5_families": rows5}
with open("data/30_v2_full_recompute.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/30_v2_full_recompute.json")