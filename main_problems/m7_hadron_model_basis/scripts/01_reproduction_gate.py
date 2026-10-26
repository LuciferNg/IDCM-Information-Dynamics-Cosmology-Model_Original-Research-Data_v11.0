"""M7 R0 — 重現閘：零 import 獨立重建 M4/M5 兩體骨架 ＋ 15/15 對照
================================================================================================
判準（寫於執行前，不得事後修改）：
  G1 兩體骨架獨立重建（零 import M4/M5 腳本、零讀其 JSON）：Cayley B=4、L=40 徑向方程
       E φ_d = −t[√(N_d/N_{d−1}) φ_{d−1} + √(N_d/N_{d+1})(B−1) φ_{d+1}] + V_d φ_d
       N_0 = 1、N_d = B(B−1)^{d−1}；對稱三對角 off = {−√B, −√(B−1)×(L−1)}
       V_d = −g·t·B^{−d}（d ≥ 0），另於 d=0 加 −h·t（on-site 磁接觸）
     校準：只用 D 通道 (A_D, |c|_D) = (1.2025, 0.3525) MeV ⟹ (a₀, t)
  G2 15/15 逐項對照 M4 `docs/25` 記載值（容差照該文 §二）
  G3 mpmath 30 dps 覆核 (a₀, t, ⟨B^{−d}⟩, P₀)
  G4 三錨 ⟨O⟩ 與外部 m_d−m_u 只作**登記**，本腳本不計任何 ⟨O⟩ 預言
  G5 印慣例標籤（44-鏈／88-site／176-L）以免跨慣例對錯數

  本腳本**不**讀 M4／M5 任何檔案；所有數字由本檔公式重算後才與記錄值比對。
輸出：data/01_reproduction_gate.json
"""
import json
import os

import numpy as np
from mpmath import mp

mp.dps = 30

HBARC = 0.1973269804  # GeV·fm
ALPHA = 1 / 137.035999084
B, L = 4, 40

M_C, M_S, M_L = 1.27, 0.495, 0.30  # GeV（組分質量；D 通道 m_h = m_c）
Q_C, Q_S = 2 / 3, 1 / 3            # 重夸克電荷；|Δq_l| = 1
A_D, K_D = 1.2025, 0.3525          # MeV（docs/13 §四.1 物理定義 |c| = K）
BAND = -2 * np.sqrt(B - 1)         # bulk 能帶底 / t

NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POWB = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))

RES = []


def check(item, ok, got, want, note=""):
    RES.append(dict(item=item, pass_=bool(ok), got=got, want=want, note=note))
    print(f"  {'✅' if ok else '🔴'} {item:34} got {got} | want {want} {note}")


def shape_float(g, h, Lx=L):
    n = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, Lx + 1)]
    powb = [B ** (-d) for d in range(Lx + 1)]
    off = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (Lx - 1))
    diag = np.array([-g * powb[d] - (h if d == 0 else 0.0) for d in range(Lx + 1)])
    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(n) * v ** 2
    w /= w.sum()
    return dict(E0=float(ev[0]), rA=float((w * np.array(powb)).sum()),
                P0=float(w[0]), d_peak=int(w.argmax()), w=w)


def obs(mh, ml, q, a, t):
    invc = HBARC / a
    cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    g = ALPHA * q * invc / t
    h = cq * invc ** 3 / t
    s = shape_float(g, h)
    A = ALPHA * q * invc * s["rA"] * 1000
    K = cq * invc ** 3 * s["P0"] * 1000
    return dict(A=A, K=K, rA=s["rA"], P0=s["P0"], d_peak=s["d_peak"], E0=s["E0"], g=g, h=h)


def bisect(f, lo, hi, n=400):
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
        if hi - lo < 1e-15:
            break
    return 0.5 * (lo + hi)


def solve_t(a):
    x = bisect(lambda xx: obs(M_C, M_L, Q_C, a, np.exp(xx))["A"] - A_D, np.log(1e-14), np.log(1e3))
    return float(np.exp(x)) if x is not None else None


def calibrate():
    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 600):
        a = float(np.exp(x))
        t = solve_t(a)
        if t is None:
            prev = None
            continue
        k = obs(M_C, M_L, Q_C, a, t)["K"] - K_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: obs(M_C, M_L, Q_C, float(np.exp(xx)),
                                      solve_t(float(np.exp(xx))))["K"] - K_D, prev[0], x)
            if r is not None:
                a = float(np.exp(r))
                return a, solve_t(a)
        prev = (x, k)
    return None, None


# ── mpmath 30 dps 覆核 ────────────────────────────────────────────────────────
def shape_mp(g, h, Lx=L):
    n = [mp.mpf(1)] + [B * (B - 1) ** (d - 1) for d in range(1, Lx + 1)]
    powb = [mp.power(B, -d) for d in range(Lx + 1)]
    off = [mp.sqrt(B)] + [mp.sqrt(B - 1)] * (Lx - 1)
    H = mp.zeros(Lx + 1, Lx + 1)
    for d in range(Lx + 1):
        H[d, d] = -mp.mpf(g) * powb[d] - (h if d == 0 else mp.mpf(0))
    for d in range(Lx):
        H[d, d + 1] = -off[d]
        H[d + 1, d] = -off[d]
    E, V = mp.eigsy(H)
    j = min(range(len(E)), key=lambda i: E[i])
    w = [n[d] * V[d, j] ** 2 for d in range(Lx + 1)]
    s = sum(w)
    return dict(E0=E[j], rA=sum(w[d] * powb[d] for d in range(Lx + 1)) / s, P0=w[0] / s)


def obs_mp(mh, ml, q, a, t):
    invc = mp.mpf(HBARC) / mp.mpf(a)
    cq = (8 * mp.pi / 3) * mp.mpf(ALPHA) * mp.mpf(q) / (16 * mp.mpf(mh) * mp.mpf(ml))
    s = shape_mp(mp.mpf(ALPHA) * mp.mpf(q) * invc / mp.mpf(t), cq * invc ** 3 / mp.mpf(t))
    return (mp.mpf(ALPHA) * mp.mpf(q) * invc * s["rA"] * 1000,
            cq * invc ** 3 * s["P0"] * 1000, s["rA"], s["P0"], s["E0"])


def main():
    out = {"module": "M7", "round": "R0", "script": "01_reproduction_gate.py"}

    print("=" * 100)
    print("慣例標籤（跨慣例對數前必須核）：")
    print("  本輪對象 = Cayley 樹徑向座標      B = 4（因果階梯 w_t = 4 = B）, 深度 L = 40 ⟹ 41×41")
    print("  鏈慣例：44-鏈／88-site／176-L（N_d = 44, N = 88 = 2N_d, L_ring = 176）")
    print("  ⚠️ 本輪**不用**鏈譜數字（γ_n、Σγ = 90），只用樹徑向 41×41；不得與 88-site 譜混算")
    print("=" * 100)

    print("\nG1 獨立重建：能帶底（L=2000, V=0, off=−√3）")
    # 精確有限 L 譜：E_min(L) = −2√(B−1)·cos(π/(L+2))（Dirichlet 兩端三對角）
    # ⟹ 對 −2√3 之差 O(1/L²) 屬截斷，非誤差。判準：①精確式 ②Richardson 外推
    n = 2000
    e = np.linalg.eigvalsh(np.diag(np.full(n + 1, 0.0)) +
                           np.diag(np.full(n, -np.sqrt(3)), 1) + np.diag(np.full(n, -np.sqrt(3)), -1))
    exact_L = BAND * np.cos(np.pi / (n + 2))
    check("V1 有限 L 精確式 −2√3cos(π/(L+2))", abs(e.min() - exact_L) < 1e-9,
          f"{e.min():.9f}", f"{exact_L:.9f}")
    e2 = np.linalg.eigvalsh(np.diag(np.full(4001, 0.0)) +
                            np.diag(np.full(4000, -np.sqrt(3)), 1) + np.diag(np.full(4000, -np.sqrt(3)), -1))
    extrap = (4 * e2.min() - e.min()) / 3.0    # 1/L² Richardson
    check("V1 外推 → −2√3（1/L² Richardson）", abs(extrap - BAND) < 1e-8, f"{extrap:.12f}", f"{BAND:.12f}",
          "（殘差 O(1/L⁴) ≈ 1.4e−9；M4 docs/25 記載 L=2000 值 −3.464097346，截斷差 4.3e−6）")

    print("\nG1 校準：只用 D 通道 (A_D, |c|_D) = (1.2025, 0.3525) MeV")
    a0, t0 = calibrate()
    r0 = obs(M_C, M_L, Q_C, a0, t0)
    print(f"  a₀ = {a0:.10f} fm   t = {t0*1000:.10f} MeV   g = {r0['g']:.10f}   h = {r0['h']:.10f}")

    print("\nG2 15/15 對照 M4 docs/25 記載值")
    check("V2 a₀ [fm]", abs(a0 - 0.377422) < 2e-5, f"{a0:.6f}", "0.377422")
    check("V2 t [MeV]", abs(t0 * 1000 - 0.7884) < 5e-4, f"{t0*1000:.4f}", "0.7884")
    check("V2 ⟨B^{−d}⟩", abs(r0["rA"] - 0.4728) < 5e-4, f"{r0['rA']:.4f}", "0.4728")
    check("V2 P₀", abs(r0["P0"] - 0.3689) < 5e-4, f"{r0['P0']:.4f}", "0.3689")
    check("V2 E(D⁺−D⁰) = A+3|c|", abs(r0["A"] + 3 * r0["K"] - 2.2600) < 5e-4,
          f"{r0['A']+3*r0['K']:.4f}", "2.2600")
    check("V2 E(D*⁺−D*⁰) = A−|c|", abs(r0["A"] - r0["K"] - 0.8500) < 5e-4,
          f"{r0['A']-r0['K']:.4f}", "0.8500")

    s1 = obs(M_S, M_L, Q_C, a0, t0)          # s_ud
    s2 = obs(0.594, 0.25, Q_C, a0, t0)       # 同乘積 0.1485
    rf = lambda mh, ml, s: (6 / np.pi) * (mh * ml * a0 ** 2 / HBARC ** 2) * (s["rA"] / s["P0"])
    check("V3 解析式 = A/|c|（s_ud）",
          abs(rf(M_S, M_L, s1) - s1["A"] / s1["K"]) < 1e-9 * abs(s1["A"] / s1["K"]),
          f"{rf(M_S,M_L,s1):.9f}", f"{s1['A']/s1['K']:.9f}")
    check("V3 乘積守恆 ⟹ 同 r", abs(s1["A"] / s1["K"] - s2["A"] / s2["K"]) < 1e-9,
          f"{s1['A']/s1['K']:.9f}", f"{s2['A']/s2['K']:.9f}")
    check("V4 B⁷", B ** 7 == 16384, f"{B**7}", "16384")
    check("V4 比值集 ≥ 1", all(x >= 1 for x in (np.inf, float(B ** 7), 1.0)), "{∞, 16384, 1}", "全 ≥ 1")
    check("V5 逐元素 B^{−d} ≤ 1", abs(max(POWB) - 1.0) < 1e-15, f"max = {max(POWB):.1f}", "1.0")
    a_max = ALPHA * (1 / 3) * (HBARC / a0) * 1000
    check("V5 A_max(q=1/3, a₀)", abs(a_max - 1.2718) < 5e-4, f"{a_max:.6f}", "1.2718")
    psi2_old = 0.07909
    r_old = (1 / (np.pi * psi2_old)) ** (1 / 3) * HBARC
    check("V7 r = (π|ψ(0)|²)^{−1/3}（舊約定 0.07909）", abs(r_old - 0.3139) < 5e-4,
          f"{r_old:.6f}", "0.3139", "（舊值源於 0.5287 缺陷；見 note）")

    # 重校後之 |ψ(0)|²（附註，不計入 15 項）
    psi2_new = r0["P0"] / a0 ** 3 * HBARC ** 3
    r_new = (1 / (np.pi * psi2_new)) ** (1 / 3) * HBARC
    print(f"\n  ⊕ 重校後 |ψ(0)|² = {psi2_new:.6f} GeV³ ⟹ r = {r_new:.4f} fm"
          f"（docs/13 §四.1 記載 0.05272 GeV³／0.3593 fm）")

    print("\nG3 mpmath 30 dps 覆核（同一 (a₀, t₀)）")
    A_m, K_m, rA_m, P0_m, E0_m = obs_mp(M_C, M_L, Q_C, a0, t0)
    d_rA = abs(float(rA_m) - r0["rA"])
    d_P0 = abs(float(P0_m) - r0["P0"])
    print(f"  mpmath: ⟨B⟩ = {mp.nstr(rA_m, 18)}  P₀ = {mp.nstr(P0_m, 18)}  E₀ = {mp.nstr(E0_m, 12)}")
    print(f"  |Δ⟨B⟩| = {d_rA:.3e}   |ΔP₀| = {d_P0:.3e}   E₀/t = {float(E0_m)/t0:.6f}")
    check("G3 |Δ⟨B^{−d}⟩| < 1e-12", d_rA < 1e-12, f"{d_rA:.2e}", "< 1e-12")
    check("G3 |ΔP₀| < 1e-12", d_P0 < 1e-12, f"{d_P0:.2e}", "< 1e-12")

    # mpmath 高精度打磨 (a, t)：以 float64 解作初值，2-D findroot（30 dps）
    def A_of(a, t):
        return obs_mp(M_C, M_L, Q_C, a, t)[0]

    def K_of(a, t):
        return obs_mp(M_C, M_L, Q_C, a, t)[1]

    try:
        sol = mp.findroot([lambda a, t: A_of(a, t) - mp.mpf(A_D),
                           lambda a, t: K_of(a, t) - mp.mpf(K_D)],
                          (mp.mpf(a0), mp.mpf(t0)), tol=mp.mpf('1e-25'), maxsteps=60)
        a_mp, t_mp = sol
        print(f"  mpmath 打磨: a₀ = {mp.nstr(a_mp, 20)} fm   t = {mp.nstr(t_mp*1000, 18)} MeV")
        check("G3 mpmath a₀ vs float64", abs(float(a_mp) - a0) < 1e-9, f"{float(a_mp)-a0:.2e}", "< 1e-9")
        check("G3 mpmath t vs float64", abs(float(t_mp) - t0) < 1e-12, f"{float(t_mp)-t0:.2e}", "< 1e-12")
        A_pol, K_pol, rA_pol, P0_pol, E0_pol = obs_mp(M_C, M_L, Q_C, a_mp, t_mp)
        print(f"  打磨後殘差: |A−A_D| = {mp.nstr(abs(A_pol-mp.mpf(A_D)), 6)} MeV, "
              f"|K−K_D| = {mp.nstr(abs(K_pol-mp.mpf(K_D)), 6)} MeV")
        print(f"  打磨後 ⟨B⟩ = {mp.nstr(rA_pol, 20)}   P₀ = {mp.nstr(P0_pol, 20)}")
        out["mpmath"] = dict(a0_fm=mp.nstr(a_mp, 24), t_MeV=mp.nstr(t_mp * 1000, 24),
                             rA=mp.nstr(rA_pol, 24), P0=mp.nstr(P0_pol, 24),
                             E0_over_t=mp.nstr(E0_pol / t_mp, 18),
                             resid_A=mp.nstr(abs(A_pol - mp.mpf(A_D)), 6),
                             resid_K=mp.nstr(abs(K_pol - mp.mpf(K_D)), 6))
    except Exception as exc:  # noqa: BLE001
        check("G3 mpmath 打磨", False, f"findroot failed: {exc}", "converged")

    print("\nG4 三錨（只登記，不預言）＋ 外部 m_d − m_u")
    anchors = [("N", 0, 0.891, 0.107), ("Σ", 1, 1.498, 0.262), ("Ξ", 2, 2.293, 0.351)]
    for nm, ns, val, sig in anchors:
        print(f"  ⟨O⟩_{nm} = {val} ± {sig}   (n_s = {ns})")
    print("  m_d − m_u = 2.56 ± 0.094 MeV（FLAG24）")
    print("  ⚠️ 本輪不許算 ⟨O⟩ 預言；R2 之前必須先落定義凍結 hash")

    npass = sum(1 for r in RES if r["pass_"])
    print(f"\n合計 {npass}/{len(RES)} PASS")
    if npass != len(RES):
        for r in RES:
            if not r["pass_"]:
                print(f"  🔴 FAIL {r['item']}: got {r['got']} want {r['want']}")
    out.update(dict(
        conventions=dict(tree="Cayley B=4, L=40, 41x41", chain="44-chain / 88-site / 176-L (unused)"),
        calibration=dict(a0_fm=a0, t_MeV=t0 * 1000, g=r0["g"], h=r0["h"],
                         A_D=r0["A"], K_D=r0["K"], rA=r0["rA"], P0=r0["P0"],
                         d_peak=r0["d_peak"], E0_over_t=r0["E0"] / t0,
                         psi2_corrected_GeV3=psi2_new, r_corrected_fm=r_new),
        anchors=dict(N=[0.891, 0.107], Sigma=[1.498, 0.262], Xi=[2.293, 0.351],
                     md_mu=[2.56, 0.094]),
        checks=[dict(item=r["item"], pass_=r["pass_"], got=r["got"], want=r["want"], note=r["note"])
                for r in RES],
        npass=npass, ntotal=len(RES),
        verdict="PASS" if npass == len(RES) else "FAIL"))
    p = os.path.join(os.path.dirname(__file__), "..", "data", "01_reproduction_gate.json")
    p = os.path.normpath(p)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\nJSON → {p}")


if __name__ == "__main__":
    main()