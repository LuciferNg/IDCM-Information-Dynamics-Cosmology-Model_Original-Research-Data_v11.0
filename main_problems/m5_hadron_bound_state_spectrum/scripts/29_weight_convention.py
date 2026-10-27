"""M5 docs/27 配套 — weight 約定裁決：由 H 之層基結構證明 w = v²（非 N_d·v²），並重算下游
================================================================================
觸發：register handoff M7→M5（2026-10-26）——「w_d = N_d·v_d² 係雙重計數；正確 w_d = v_d²」。
本腳本唔靠裁決，而係由**算子結構**判定：

  A. 由 K_tree 顯式建 b-ary 樹（b = 4，深度 D），H = −t·A + diag(V_d)；
     對稱扇區基 |d⟩ = N_d^(−1/2) Σ_{i∈layer d} |i⟩ ⟹ 矩陣元
       ⟨d|H|d⟩ = V_d,
       ⟨0|H|1⟩ = −t·√N_1 = −t·√b,
       ⟨d|H|d+1⟩ = −t·√(b−1)  (d ≥ 1)
     ⟹ 與生產 H 之 OFF = [−√B, −√(B−1), …] **逐項相同** ⟹ 生產 H 已在**歸一化層基**中
     ⟹ 其本徵向量分量即 Φ_d = √N_d·ψ_d ⟹ **層機率 = |v_d|²**（w = N_d·v² 係把 Φ_d 當成 per-site 振幅 ⟹ 雙重計數）
  B. 以 w = v² 重校 D 兩數，與 M7 提供之重校值對照（a = 0.384184 fm、t = 1.7711 MeV、⟨B⟩ = 0.4812、P₀ = 0.3891）。
  C. 以 w = v² 重算下游：束縛守門、天花板 R、r_K 上限、P₀ 比、所需 Δd、A3 過強倍數；並列 v1（w = N_d v²）對照。

判準（寫死）：
  W1 樹投影與路徑 H_eff 之本徵值須一致（相對差 < 1e−10），否則結構論證不成立 ⟹ 記 FAIL。
  W2 生產 H 之 OFF 須逐項等於 [−√B, −√(B−1)×] ⟹ 否則無法聲稱層基。
  W3 重校須 D 兩數零殘差（|ΔA| + |Δc| < 1e−9）且束縛；與 M7 值相對差 < 1e−3 視為獨立確認（M7 用不同求解器）。
  W4 下游只報數，不作新擬合；v1 數字保留作對照並標作廢。

輸出：`data/29_weight_convention.json`
"""
import json

import numpy as np

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
A_D, C_D = 1.2025, 0.3525
E1K, E2K = 1.3740, 1.3200
R_DATA_D = A_D / C_D
R_DATA_K = ((E1K + 3 * E2K) / 4) / ((E1K - E2K) / 4)
KAPPA = 0.2150
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 106)
print("M5 docs/27 — weight 約定裁決（結構證明）＋下游重算")
print("=" * 106)

# ---------------- A 結構證明 ----------------
print("\n[A] 顯式 b-ary 樹投影 vs 路徑 H_eff（b = 4，深度 3）")
b, Dp = 4, 3
layers, idx = [], 0
for d in range(Dp + 1):
    n = 1 if d == 0 else b * (b - 1) ** (d - 1)
    layers.append(list(range(idx, idx + n)))
    idx += n
Ntot = idx
A = np.zeros((Ntot, Ntot))
for d in range(1, Dp + 1):               # 逐層接父：層 1 每父 b 子，之後每父 b−1 子
    kids = b if d == 1 else b - 1
    for k, j in enumerate(layers[d]):
        i = layers[d - 1][k // kids]
        A[i, j] = A[j, i] = 1.0
V = np.array([-1.0 * POW[d] - (3.0 if d == 0 else 0.0) for d in range(Dp + 1)])
Hfull = -1.0 * A + np.diag([V[d] for d in range(Dp + 1) for _ in layers[d]])

# 對稱扇區正交基
Q = np.zeros((Ntot, Dp + 1))
for d in range(Dp + 1):
    for i in layers[d]:
        Q[i, d] = 1.0 / np.sqrt(len(layers[d]))
S = Q.T @ Hfull @ Q
Heff = np.diag(V)
Heff[0, 1] = Heff[1, 0] = -np.sqrt(b)
for d in range(1, Dp):
    Heff[d, d + 1] = Heff[d + 1, d] = -np.sqrt(b - 1)
evS, evH = np.sort(np.linalg.eigvalsh(S)), np.sort(np.linalg.eigvalsh(Heff))
dev_ev = float(np.max(np.abs(evS - evH)))
print(f"    樹節點 {Ntot}、層大小 {[len(l) for l in layers]}；矩陣元 ⟨0|H|1⟩ = {S[0,1]:.12f}（−√b = {-np.sqrt(b):.12f}）")
print(f"    ⟨1|H|2⟩ = {S[1,2]:.12f}（−√(b−1) = {-np.sqrt(b-1):.12f}）")
print(f"    本徵值：樹投影 {np.round(evS,10)}")
print(f"            路徑 Heff {np.round(evH,10)}  ⟹ 最大差 {dev_ev:.3e} ⟹ {'✅ PASS' if dev_ev < 1e-10 else '🔴 FAIL'}")
print(f"    ⟹ 歸一化層基之 Hamiltonian **就是** 常數 hopping 路徑 H（生產 OFF 結構）⟹ v_d = Φ_d = √N_d·ψ_d")
print(f"    ⟹ 層機率 = |v_d|² ⟹ w = v²（w = N_d v² 係雙重計數）")

W2 = (abs(-OFF[0] - np.sqrt(B)) < 1e-15) and all(abs(-OFF[d] - np.sqrt(B - 1)) < 1e-15 for d in range(1, L))
print(f"    W2 生產 OFF 逐項 = [−√B, −√(B−1)×{L-1}]：{'✅ PASS' if W2 else '🔴 FAIL'}")

# ---------------- B 重校 ----------------
def amps(a, t, q, mh, ml, conv):
    gcoef = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-gcoef * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag)
    H[np.arange(L), np.arange(1, L + 1)] = OFF
    H[np.arange(1, L + 1), np.arange(L)] = OFF
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = v ** 2 if conv == "v2" else np.array(NODES) * v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    bd = bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX)
    return ALPHA * q * (HBARC / a) * rA * 1000, Cq * (HBARC / a) ** 3 * P0 * 1000, bd, rA, P0, dpk


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


def calib(kappa, conv):
    chd = kappa * G7 * Q_C
    target = A_D - chd
    if target <= 0:
        return None

    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), Q_C, M_C, M_L, conv)[0] - target, np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 300):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, Q_C, M_C, M_L, conv)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))), Q_C, M_C, M_L, conv)[1] - C_D, prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


print("\n[B] 以 w = v² 重校 D 兩數（κ = 0.2150）並與 M7 交接值對照")
sol2 = calib(KAPPA, "v2")
sol1 = calib(KAPPA, "v1")
M7 = {"a": 0.384184, "t": 1.7711, "rB": 0.4812, "P0": 0.3891}
for tag, sol in (("w = v²（新）", sol2), ("w = N_d v²（記錄）", sol1)):
    if sol is None:
        print(f"    {tag}：無解 🔴")
        continue
    a, t = sol
    Ak, ck, bd, rA, P0, dpk = amps(a, t, Q_C, M_C, M_L, "v2" if sol is sol2 else "v1")
    print(f"    {tag}：a = {a:.6f} fm、t = {t*1000:.4f} MeV、⟨B⟩ = {rA:.4f}、P₀ = {P0:.4f}、d_peak = {dpk}、束縛 {bd}")
    if sol is sol2:
        d_a, d_t = abs(a - M7["a"]) / M7["a"], abs(t * 1000 - M7["t"]) / M7["t"]
        print(f"      vs M7 交接 (0.384184, 1.7711)：相對差 a {d_a:.2e}、t {d_t:.2e} ⟹ {'✅ 獨立確認' if max(d_a,d_t) < 1e-3 else '🟡 未吻合'}")
        print(f"      vs M7 ⟨B⟩ 0.4812 / P₀ 0.3891：差 {abs(rA-M7['rB']):.2e} / {abs(P0-M7['P0']):.2e}")

# ---------------- B2 根枚舉（多值性） ----------------
print("\n[B2] 根枚舉（κ = 0.2150，兩種約定各自全部根）—— 校準多值，必須先定位分支")
def all_roots(kappa, conv, ngrid=400):
    chd = kappa * G7 * Q_C
    target = A_D - chd
    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), Q_C, M_C, M_L, conv)[0] - target, np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None
    roots, prev = [], None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, Q_C, M_C, M_L, conv)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))), Q_C, M_C, M_L, conv)[1] - C_D, prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                tt = t_solve(aa)
                if tt and not any(abs(aa - q[0]) < 1e-4 and abs(tt - q[1]) < 1e-6 for q in roots):
                    roots.append((aa, tt))
        prev = (x, k)
    return roots

r1, r2 = all_roots(KAPPA, "v1"), all_roots(KAPPA, "v2")
print(f"    w = N_d v²（記錄）：{len(r1)} 個根 ⟹ " + "；".join(f"({a:.6f} fm, {t*1000:.4f} MeV)" for a, t in r1))
print(f"    w = v²（新）    ：{len(r2)} 個根 ⟹ " + "；".join(f"({a:.6f} fm, {t*1000:.4f} MeV)" for a, t in r2))
hit = [q for q in r2 if abs(q[0] - M7["a"]) / M7["a"] < 1e-3 and abs(q[1] * 1000 - M7["t"]) / M7["t"] < 1e-3]
print(f"    與 M7 交接 (0.384184 fm, 1.7711 MeV) 吻合之 w=v² 根：{len(hit)} 個" +
      ("" if not hit else f" ⟹ ({hit[0][0]:.6f} fm, {hit[0][1]*1000:.4f} MeV) ✅ 獨立確認 M7"))
print(f"    記錄值 (0.377422 fm, 0.7884 MeV) 落在：w=v² 根 {[i for i,q in enumerate(r2) if abs(q[0]-0.377422)<1e-3]}；w=N_d v² 根 {[i for i,q in enumerate(r1) if abs(q[0]-0.377422)<1e-3]}")

# ---------------- C 下游重算 ----------------
print("\n[C] 下游重算（w = v²；同一 κ 網格、同一守門）")
rows = []
for kp in np.arange(0.150, 0.3201, 0.0025):
    s = calib(float(kp), "v2")
    if s is None:
        continue
    a, t = s
    AK, cK, bdK, _, _, dpK = amps(a, t, Q_S, M_S, M_L, "v2")
    AD, cD, bdD, _, _, dpD = amps(a, t, Q_C, M_C, M_L, "v2")
    rK, rD = (AK + kp * G0 * Q_S) / cK, (AD + kp * G7 * Q_C) / cD
    rows.append({"kappa": float(kp), "a0": a, "t": t * 1000, "rK": rK, "rD": rD, "R": rK / rD,
                 "P0K": amps(a, t, Q_S, M_S, M_L, "v2")[4], "P0D": amps(a, t, Q_C, M_C, M_L, "v2")[4],
                 "bdK": bdK, "bdD": bdD, "dpK": dpK, "dpD": dpD})
print(f"    掃描點 {len(rows)} 個（κ ∈ [0.150, 0.320]）")
bothb = [r for r in rows if r["bdK"] and r["bdD"]]
print(f"    K 與 D 同時束縛：{len(bothb)}／{len(rows)}")
if bothb:
    mx = max(bothb, key=lambda r: r["R"])
    ceil = mx["R"] * R_DATA_D
    print(f"    R_max = {mx['R']:.4f}（κ = {mx['kappa']:.4f}、d_peak K/D = {mx['dpK']}/{mx['dpD']}）")
    print(f"    ⟹ r_K 天花板 = {ceil:.4f}；數據 98.7778 ⟹ 短缺 {R_DATA_K/ceil:.2f}×")
    print(f"    P₀(K)/P₀(D) = {mx['P0K']/mx['P0D']:.5f} = 1/{mx['P0D']/mx['P0K']:.3f}")
    dyn = mx["P0K"] / mx["P0D"]
    data_c = 0.0135 / 0.3525
    need = -np.log(data_c / dyn) / np.log(B)
    print(f"    所需 Δd = log₄(數據/動力學) = {need:.4f}；Δd = 1 覆蓋 {dyn*0.25/data_c*100:.2f}%")
    print(f"    A3（Δd = −7）過強 {16384*dyn/data_c:.1f}×")
else:
    mx, ceil, dyn, need = None, None, None, None

out = {"structural": {"tree_vs_path_max_dev": dev_ev, "W2_off_match": bool(W2),
                      "conclusion": "層基 ⟹ w = |v|²；w = N_d v² 雙重計數（M7 正確）"},
       "calib_v2": None if sol2 is None else {"a0_fm": sol2[0], "t_MeV": sol2[1] * 1000},
       "calib_v1": None if sol1 is None else {"a0_fm": sol1[0], "t_MeV": sol1[1] * 1000},
       "M7_cross_check": M7,
       "roots_v1": [[a, t*1000] for a, t in r1], "roots_v2": [[a, t*1000] for a, t in r2], "M7_root_found": bool(hit),
       "scan": rows,
       "ceiling": None if mx is None else {"R_max": mx["R"], "kappa": mx["kappa"], "rK_ceiling": ceil,
                                           "shortfall": R_DATA_K / ceil, "P0_ratio": dyn,
                                           "required_dd": float(need),
                                           "d1_coverage": dyn * 0.25 / (0.0135 / 0.3525),
                                           "a3_over": 16384 * dyn / (0.0135 / 0.3525)}}
with open("data/29_weight_convention.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/29_weight_convention.json")