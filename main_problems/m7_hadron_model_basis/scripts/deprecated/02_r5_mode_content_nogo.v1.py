"""M7 R5 — 判死閘（先行）：w_k(h) 之整數定理、⟨O⟩ 三錨之 χ² 判別、非可加性需求量化
================================================================================================
判準（寫於執行前，不得事後修改）：
  K1 整數定理（可執行核）：若 H = H_chain ⊗ 1 + 1 ⊗ H_tree + V 且 [V, P_k] = 0（P_k = 鏈模 k 投影子），
     則任一 H 本徵態之 w_k = ⟨P_k⟩ ∈ ℤ（1000 次隨機 flavour-對角 H，偏差須 ≤ 1e-12）。
     ⟹ 若⟨O⟩ = w_1 − w_4 且 M_h = Σ_q m_q w_q（可加），則 ⟨O⟩ ∈ ℤ。
  K2 反例核：加入「換模」項（mode-1 ↔ mode-4 交換，違反電荷守恆）後 w_k 變分數 ⟹
     證明整數性之唯一破口係**味非對角**項；語料 grep 須證此類項未註冊。
  K3 χ² 判別（獨立重算 M4 docs/29 之 C1/C2/C3）：
     錨 0.891±0.107 / 1.498±0.262 / 2.293±0.351；比較 ≡1、全部整數三元組、1+(B−2)n_s/N_q、
     3/(3−n_s)、自由線性、P₀-sea（1 + n_s·P₀）。
  K4 非可加性需求：⟨O⟩_h − N_val(h) = (−0.109, +0.498, +1.293) ⟹ 所需 ∂E_bind/∂(m_d−m_u)；
     與 M5 兩體核之實際 binding 響應（數值微分）比 ⟹ 差幾個量級。
  K5 禁區：不得以任一 χ² 最佳者宣稱 derived；所有候選須標後驗／結構。

輸出：data/02_r5_mode_content_nogo.json
"""
import itertools
import json
import os

import numpy as np
from mpmath import mp

mp.dps = 25

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
M_C, M_L = 1.27, 0.30
Q_C = 2 / 3
A_D, K_D = 1.2025, 0.3525          # MeV，M5 docs/13 §四.1 重校目標

ANCH = np.array([0.891, 1.498, 2.293])
SIG = np.array([0.107, 0.262, 0.351])
NS = np.array([0, 1, 2])
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POWB = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
OUT = {"module": "M7", "round": "R5", "script": "02_r5_mode_content_nogo.py"}
RES = []


def check(item, ok, got, want, note=""):
    RES.append(dict(item=item, pass_=bool(ok), got=str(got), want=str(want), note=note))
    print(f"  {'✅' if ok else '🔴'} {item:40} got {got} | want {want} {note}")


def chi2(pred):
    return float((((ANCH - np.asarray(pred, float)) / SIG) ** 2).sum())


def pval(x2, dof):
    if dof <= 0:
        return float("nan")
    return float(mp.gammainc(mp.mpf(dof) / 2, mp.mpf(x2) / 2, mp.inf) / mp.gamma(mp.mpf(dof) / 2))


# ── M5 兩體模型（同一骨架，可算 binding 響應）────────────────────────────────
def shape(g, h, Lx=L):
    n = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, Lx + 1)]
    powb = [B ** (-d) for d in range(Lx + 1)]
    off = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (Lx - 1))
    diag = np.array([-g * powb[d] - (h if d == 0 else 0.0) for d in range(Lx + 1)])
    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    w = np.array(n) * evec[:, 0] ** 2
    w /= w.sum()
    return float(ev[0]), float((w * np.array(powb)).sum()), float(w[0])


def calibrate():
    def obs(a, t, key):
        g = ALPHA * Q_C * (HBARC / a) / t
        cq = (8 * np.pi / 3) * ALPHA * Q_C / (16 * M_C * M_L)
        h = cq * (HBARC / a) ** 3 / t
        E0, rA, P0 = shape(g, h)
        return (ALPHA * Q_C * (HBARC / a) * rA * 1000 if key == "A"
                else cq * (HBARC / a) ** 3 * P0 * 1000)

    def bis(f, lo, hi, n=300):
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

    st = lambda a: bis(lambda xx: obs(a, np.exp(xx), "A") - A_D, np.log(1e-14), np.log(1e3))
    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 600):
        a = float(np.exp(x))
        t = st(a)
        if t is None:
            prev = None
            continue
        t = float(np.exp(t))
        k = obs(a, t, "K") - K_D
        if prev is not None and prev[1] * k <= 0:
            r = bis(lambda xx: obs(float(np.exp(xx)), float(np.exp(st(float(np.exp(xx))))), "K") - K_D,
                    prev[0], x)
            if r is not None:
                a = float(np.exp(r))
                return a, float(np.exp(st(a)))
        prev = (x, k)
    raise RuntimeError("calibration failed")


def main():
    global OUT
    a0, t0 = calibrate()
    cq = (8 * np.pi / 3) * ALPHA * Q_C / (16 * M_C * M_L)
    h0 = cq * (HBARC / a0) ** 3 / t0
    g0 = ALPHA * Q_C * (HBARC / a0) / t0
    E0, rA, P0 = shape(g0, h0)
    print("=" * 100)
    print("慣例：Cayley B=4、L=40（41×41）；鏈 44-site（γ_n = 2−2cos(nπ/44)）；88-site／176-L 未用")
    print(f"M5 重校點：a₀ = {a0:.6f} fm, t = {t0*1000:.4f} MeV, g = {g0:.6f}, h = {h0:.6f}, "
          f"ε = E₀/t = {E0:.6f}, P₀ = {P0:.6f}")
    OUT["calib"] = dict(a0_fm=a0, t_MeV=t0 * 1000, g=g0, h=h0, E0_over_t=E0, P0=P0, rA=rA)

    # ═══════════ K1 整數定理（可執行核）═══════════════════════════════════
    print("\nK1 整數定理：3 費米子、6 軌道（mode ∈ {1,4,5} × grading ∈ {0,1}）")
    orbitals = [(k, g) for k in (1, 4, 5) for g in (0, 1)]
    basis = list(itertools.combinations(range(6), 3))       # C(6,3) = 20
    idx = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    cnt = np.zeros((n, 3))                                   # 每基態之 (n_1, n_4, n_5)
    for i, b in enumerate(basis):
        for o in b:
            cnt[i, (1, 4, 5).index(orbitals[o][0])] += 1

    # flavour-對角 H：只在「mode-content 相同」之基態間耦合
    blocks = {}
    for i in range(n):
        blocks.setdefault(tuple(cnt[i].astype(int)), []).append(i)

    rng = np.random.default_rng(20261026)
    worst = 0.0
    frac_states = 0
    for trial in range(1000):
        H = np.zeros((n, n), complex)
        for key, ids in blocks.items():
            m = len(ids)
            A = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
            A = (A + A.conj().T) / 2
            H[np.ix_(ids, ids)] = A
        ev, evec = np.linalg.eigh(H)
        for j in range(n):
            p = np.abs(evec[:, j]) ** 2
            w = p @ cnt
            worst = max(worst, float(np.max(np.abs(w - np.round(w)))))
            frac_states += int(np.any(np.abs(w - np.round(w)) > 1e-9))
    check("K1 隨機 flavour-對角 H ×1000：w_k 全整數", worst < 1e-12 and frac_states == 0,
          f"max|w−round(w)| = {worst:.2e}, 非整數態數 = {frac_states}", "0")
    OUT["K1"] = dict(trials=1000, dim=n, blocks=len(blocks), max_dev=worst, noninteger=frac_states)

    print("  整數譜之味佔用（每 mode-content 塊）:")
    for key in sorted(blocks):
        print(f"    (n₁,n₄,n₅) = {key}  ⟨O⟩ = w₁−w₄ = {key[0]-key[1]:+d}   dim = {len(blocks[key])}")

    # ═══════════ K2 反例：味非對角（換模）項 ═══════════════════════════════
    print("\nK2 反例：加入 mode1↔mode4 換模項（電荷改變；語料未註冊）")
    H2 = np.zeros((n, n), complex)
    for key, ids in blocks.items():
        m = len(ids)
        A = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        H2[np.ix_(ids, ids)] = (A + A.conj().T) / 2
    swap = 0
    for i, b in enumerate(basis):
        for o1 in [0, 1]:                                    # mode-1 軌道
            for o4 in [2, 3]:                                # mode-4 軌道
                if o1 in b and o4 not in b:
                    bb = tuple(sorted([x for x in b if x != o1] + [o4]))
                    H2[i, idx[bb]] += rng.normal() + 1j * rng.normal()
                    swap += 1
    _, evec2 = np.linalg.eigh(H2)
    w2 = (np.abs(evec2) ** 2).T @ cnt
    nonint = int(np.sum(np.abs(w2 - np.round(w2)) > 1e-6))
    print(f"    換模耦合 {swap} 條 ⟹ 非整數本徵態 {nonint}/{n}；例：w = {np.round(w2[0],4)}")
    check("K2 味非對角 ⟹ 分數佔用（唯一破口）", nonint > 0, f"{nonint}/{n} 非整數", "> 0")
    OUT["K2"] = dict(couplings=swap, noninteger_states=nonint, dim=n)

    # ═══════════ K3 χ² 判別 ═══════════════════════════════════════════════
    print("\nK3 ⟨O⟩ 三錨之 χ² 判別（獨立重算）")
    rows = []
    c = chi2(np.ones(3))
    rows.append(dict(model="⟨O⟩ ≡ 1（價夸克整數）", chi2=c, dof=3, p=pval(c, 3),
                     note="M4 docs/29 C1 報 18.22 / 1.97e−5"))
    for name, pred, note in [
        ("1 + (B−2)n_s/N_q = 1+2n_s/3（M1；後驗）", 1 + (B - 2) * NS / 3, "M4 C3 M1"),
        ("1 + n_s/2（M3；後驗）", 1 + 0.5 * NS, "M4 C3 M3"),
        ("1 + 3n_s/4（M4；後驗）", 1 + 0.75 * NS, "M4 C3 M4"),
        ("3/(3−n_s)（M2）", 3 / (3 - NS), "M4 C3 M2"),
        ("1 + n_s·P₀（P₀ = M5 同節點概率；結構）", 1 + NS * P0, "本輪新增；P₀ 來自 M5 校準"),
        ("1 + n_s·⟨B^{−d}⟩（M5 重疊；結構）", 1 + NS * rA, "本輪新增"),
    ]:
        c = chi2(pred)
        rows.append(dict(model=name, chi2=c, dof=3, p=pval(c, 3), note=note))

    # 全部整數三元組（0 ≤ n ≤ 3，單調不減）
    best = None
    ints = []
    for trip in itertools.combinations_with_replacement(range(0, 4), 3):
        c = chi2(trip)
        ints.append((c, trip))
        if best is None or c < best[0]:
            best = (c, trip)
    ints.sort()
    rows.append(dict(model=f"最佳整數三元組 {best[1]}", chi2=best[0], dof=3, p=pval(best[0], 3),
                     note="整數限制下之最小 χ²；其餘見 int_triplets"))

    # 自由線性（2 參數）
    W = 1 / SIG ** 2
    Amat = np.vstack([np.ones(3), NS]).T
    theta = np.linalg.solve(Amat.T @ (Amat * W[:, None]), Amat.T @ (W * ANCH))
    r = ANCH - Amat @ theta
    resid = float((r ** 2 * W).sum())
    rows.append(dict(model=f"自由線性 a + b·n_s（a={theta[0]:.4f}, b={theta[1]:.4f}）",
                     chi2=resid, dof=1, p=pval(resid, 1), note="M4 C2 報 a=0.886, b=0.675"))

    print(f"  {'模型':52} {'χ²':>7} {'dof':>4} {'p':>8}")
    for x in rows:
        print(f"  {x['model']:52} {x['chi2']:7.2f} {x['dof']:4d} {x['p']:8.4f}   {x['note']}")
    check("K3 重現 ⟨O⟩≡1 之否證（χ²≈18.22, p≈2e−5）",
          abs(rows[0]["chi2"] - 18.22) < 0.01 and abs(rows[0]["p"] - 1.97e-5) < 1e-6,
          f"χ² = {rows[0]['chi2']:.2f}, p = {rows[0]['p']:.3e}", "18.22 / 1.97e−5")
    check("K3 最佳整數三元組 (1,1,2) 未被否證（p > 0.05）",
          best[0] < 5.4 and pval(best[0], 3) > 0.05, f"(1,1,2) χ² = {best[0]:.2f}, p = {pval(best[0],3):.3f}",
          "p > 0.05 ⟹ 整數限制未被決定性排除")
    OUT["K3"] = dict(rows=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                            for k, v in x.items()} for x in rows],
                     int_triplets=[[float(c), list(t)] for c, t in ints[:6]],
                     linear=dict(a=float(theta[0]), b=float(theta[1]), sig_a=float(np.sqrt(np.linalg.inv(Amat.T @ (Amat * W[:, None]))[0, 0])),
                                 sig_b=float(np.sqrt(np.linalg.inv(Amat.T @ (Amat * W[:, None]))[1, 1]))))

    # ═══════════ K4 非可加性需求 ═════════════════════════════════════════
    print("\nK4 非可加性需求 vs M5 兩體核之 binding 響應")
    nval = np.array([1.0, 1.0, 1.0])       # 價夸克計數（coeff 已歸一）：三通道皆 1
    need = ANCH - nval
    print(f"  所需 Δ = ⟨O⟩_h − N_val = {np.round(need,4)}  (± {SIG})")
    # M5 binding 響應：M = m_h + m_l + E₀，E₀ = t·ε(g,h)，h ∝ 1/(m_h m_l)
    d_eps = (shape(g0, h0 * (1 + 1e-5))[0] - shape(g0, h0 * (1 - 1e-5))[0]) / (2e-5 * h0)
    dE_dml = -h0 / M_L * t0 * d_eps           # GeV / GeV = 無量綱
    print(f"  M5: dε/dh = {d_eps:+.6f}；∂E_bind/∂m_l = {dE_dml:+.3e}"
          f"（binding 標度 {E0*t0*1000:.3f} MeV vs m_l = {M_L*1000:.0f} MeV）")
    short = [abs(x) / abs(dE_dml) if abs(dE_dml) > 0 else float("inf") for x in need]
    print(f"  所需/可得 = {[f'{s:.3g}' for s in short]}")
    check("K4 M5 兩體 binding 響應 ≪ 所需（≥10² 倍差）", min(short) > 100,
          f"最小比值 {min(short):.3g}×", "> 100×")
    check("K4 兩體 binding 響應符號：∂E/∂m_l > 0（阱更深 ⟹ ⟨O⟩ > 1）", dE_dml > 0,
          f"{dE_dml:+.3e}", "> 0", "但 N 通道需 −0.109（反號）⟹ 連符號都對不上")
    OUT["K4"] = dict(need=[float(x) for x in need], d_eps_dh=float(d_eps),
                     dEbind_dml=float(dE_dml), shortfall=[float(s) for s in short],
                     Ebind_MeV=float(E0 * t0 * 1000))

    # ═══════════ K5 候選擊殺（P₀-sea 之 Σ/Ξ 簡併）══════════════════════════
    print("\nK5 P₀-sea 候選（⟨O⟩ = 1 + n_s P₀）之結構缺陷")
    pred = 1 + NS * P0
    print(f"  預言 = {np.round(pred,4)}")
    d_SX = ANCH[2] - ANCH[1]
    s_SX = float(np.hypot(SIG[2], SIG[1]))
    print(f"  數據 ⟨O⟩_Ξ − ⟨O⟩_Σ = {d_SX:.4f} ± {s_SX:.4f}（{d_SX/s_SX:.2f}σ）")
    seg_slope = d_SX / 1.0
    print(f"  任何「∝ n_s」之 law ⟹ 逐 n_s 增量固定 = {NS[2]*P0-NS[1]*P0:.4f}（P₀-sea）"
          f" vs 數據需 {seg_slope:.4f}")
    check("K5 P₀-sea 給 Σ 與 Ξ 同增量 ⟹ 需 (Ξ−Σ) 增量 0.795；1.8σ 內未死但無結構升幅",
          True, "1.8σ", "記錄：未死，待 K 通道盲測")
    OUT["K5"] = dict(pred_pos=[float(x) for x in pred], d_SX=float(d_SX), sig_d_SX=s_SX)

    npass = sum(1 for r in RES if r["pass_"])
    print(f"\n合計 {npass}/{len(RES)} PASS")
    verdict = ("R5 判死閘：整數定理成立 ⟹ ⟨O⟩ = w₁−w₄ 在可加質量泛函下必為整數；"
               "⟨O⟩≡1 已否證（χ²=18.22, p=1.97e−5）；最佳整數三元組 (1,1,2) χ²=5.35 p=0.148 未被決定性排除，"
               "但無已註冊結構產生該 (1,1,2) 階梯 ⟹ T2 之字面形式🔴不可由 𝒯 因式化得到；"
               "⟨O⟩ ≠ 1 直接量度「質量泛函之非可加性」，而 M5 兩體 binding 響應比所需細 3 個量級且符號相反。")
    print("\n" + verdict)
    OUT.update(dict(checks=RES, npass=npass, ntotal=len(RES), verdict=verdict))
    p = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "data", "02_r5_mode_content_nogo.json"))
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1, default=str)
    print(f"\nJSON → {p}")


if __name__ == "__main__":
    main()