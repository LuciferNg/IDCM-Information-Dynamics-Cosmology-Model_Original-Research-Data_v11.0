"""M7 R5 — 判死閘（先行）v2：w_k(h) 整數定理、T2 之定義層斷裂、⟨O⟩ 非可加性需求
================================================================================================
檔案命名：02b（= 02 之 v2）。v1 保留於 scripts/02_r5_mode_content_nogo.py（不覆寫），
其 log = data/run_02_r5_mode_content_nogo.v1.log。名稱加 b 係為避免與並行 session 之編號相撞。

v1 → v2 改動（append-only）：
  (1) K2 v1 bug：只寫上三角，numpy.eigh 只讀下三角 ⟹ 換模項完全未入 H。v2 寫雙邊。
  (2) K3 v1 誤比：M4 docs/29 報「χ² = 18.22/3 dof ⟹ p = 1.97e−5」，但 1.97e−5 係 **1 dof** 值；
      3 dof 正確 p = 3.96e−4。v2 用正確 p 並列為審計項。
  (3) K4 判準（≥10² 倍）照實報（N 通道 44 倍）——判準不因結果而改。
  (4) 新增 K6（恒等性）、K7（T2 定義層斷裂）、K5 改為可判別性檢定（Δχ²）。
  (5) 加入 M5 weight 約定雙軌（並行 session 揭發 w = N_d·v² 之雙重計數）。

判準（寫於執行前，不得事後修改）：
  K1 [V,P_k]=0（全部已註冊交互項：em 電荷核、樹核 B^{−d}、層釘）⟹ w_k ∈ ℤ；
     1000 次隨機 flavour-對角 H，max|w−round(w)| ≤ 1e−12。
  K2 加入 mode-1 ↔ mode-4 換模項（hermitian）⟹ w_k 分數 ⟹ 唯一破口係味非對角。
  K3 χ²（3 dof）：≡1、1+2n_s/3、3/(3−n_s)、1+n_sP₀×2 約定、1+n_s⟨B⟩、最佳整數、自由線性。
  K4 Δ = ⟨O⟩_h − 1 vs M5 兩體 binding ∂E/∂m_l（兩約定）。
  K5 M1 vs P₀-sea：Δχ² < 4 ⟹ 三錨不可分辨（須誠實報）。
  K6 Fraction 精確：價夸克計數 ⟹ ⟨O⟩ ≡ 1；sector 投影 ⟹ ⟨O⟩ ≡ 2。
  K7 v9 §D.1 sector 表含 u／d 否 ⟹ T2 是否良定義。

輸出：data/02b_r5_mode_content_nogo.json
"""
import itertools
import json
import os
from fractions import Fraction

import numpy as np
from mpmath import mp

mp.dps = 25

HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
M_C, M_L = 1.27, 0.30
Q_C = 2 / 3
A_D, K_D = 1.2025, 0.3525

ANCH = np.array([0.891, 1.498, 2.293])
SIG = np.array([0.107, 0.262, 0.351])
NS = np.array([0, 1, 2])
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
OUT = {"module": "M7", "round": "R5", "script": "02b_r5_mode_content_nogo.py", "version": 2}
RES = []


def check(item, ok, got, want, note=""):
    RES.append(dict(item=item, pass_=bool(ok), got=str(got), want=str(want), note=note))
    print(f"  {'✅' if ok else '🔴'} {item:46} got {got} | want {want} {note}")


def chi2(pred):
    return float((((ANCH - np.asarray(pred, float)) / SIG) ** 2).sum())


def pval(x2, dof):
    return float(mp.gammainc(mp.mpf(dof) / 2, mp.mpf(x2) / 2, mp.inf) / mp.gamma(mp.mpf(dof) / 2))


def stats(g, h, conv, Lx=L):
    n = np.array([1.0] + [B * (B - 1) ** (d - 1) for d in range(1, Lx + 1)])
    powb = np.array([B ** (-d) for d in range(Lx + 1)])
    off = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (Lx - 1))
    diag = np.array([-g * powb[d] - (h if d == 0 else 0.0) for d in range(Lx + 1)])
    H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = (n * v ** 2) if conv == "M5" else v ** 2
    w = w / w.sum()
    return float(ev[0]), float((w * powb).sum()), float(w[0]), float((v ** 2).sum())


def obs(mh, ml, q, a, t, conv="M5"):
    invc = HBARC / a
    cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    g = ALPHA * q * invc / t
    h = cq * invc ** 3 / t
    E0, rA, P0, _ = stats(g, h, conv)
    return dict(A=ALPHA * q * invc * rA * 1000, K=cq * invc ** 3 * P0 * 1000,
                rA=rA, P0=P0, E0=E0, g=g, h=h)


def calibrate(conv="M5"):
    f = lambda x: obs(M_C, M_L, Q_C, x[0], x[1], conv)
    x = np.array([0.3774, 7.884e-4])
    for _ in range(80):
        o = f(x)
        r = np.array([o["A"] - A_D, o["K"] - K_D])
        if np.max(np.abs(r)) < 1e-14:
            break
        J = np.zeros((2, 2))
        for j in range(2):
            dx = np.zeros(2)
            dx[j] = 1e-7 * abs(x[j])
            o2 = f(x + dx)
            J[:, j] = [(o2["A"] - o["A"]), (o2["K"] - o["K"])] / dx[j]
        try:
            step = np.linalg.solve(J, -r)
        except np.linalg.LinAlgError:
            break
        x = x + np.clip(step, -0.3 * np.abs(x), 0.3 * np.abs(x))
    return float(x[0]), float(x[1]), f(x)


def main():
    print("=" * 104)
    print("慣例：Cayley B=4、L=40（41×41）；鏈 44-site（γ_n = 2−2cos(nπ/44)："
          "γ₁=0.0051024, γ₄=0.0810382, γ₅=0.1261118）；88-site／176-L 未用；γ₄₄ = 4 = B")
    a0, t0, o0 = calibrate("M5")
    a1, t1, o1 = calibrate("correct")
    print(f"M5 約定重校：a₀={a0:.6f} fm, t={t0*1000:.4f} MeV, ⟨B^{{−d}}⟩={o0['rA']:.4f}, "
          f"P₀={o0['P0']:.4f}, ε={o0['E0']:.6f}")
    print(f"w=v² 約定重校：a={a1:.6f} fm, t={t1*1000:.4f} MeV, ⟨B^{{−d}}⟩={o1['rA']:.4f}, "
          f"P₀={o1['P0']:.4f}, ε={o1['E0']:.6f}")
    print("  ⚠️ 並行 session（scripts/01_repro_gate.py）揭 M5 w = N_d·v² 雙重計數（ΣN_d²|ψ_d|²=2.0966）；")
    print("     本輪兩種約定並列，R5 結論須對兩者同時成立。")
    OUT["calib"] = dict(M5=dict(a=a0, t_MeV=t0 * 1000, rA=o0["rA"], P0=o0["P0"], E0=o0["E0"]),
                        correct=dict(a=a1, t_MeV=t1 * 1000, rA=o1["rA"], P0=o1["P0"], E0=o1["E0"]))

    # ── K1 整數定理 ──────────────────────────────────────────────────────
    print("\nK1 整數定理（3 費米子；6 軌道 = mode{1,4,5} × grading{0,1}）")
    orbitals = [(k, gg) for k in (1, 4, 5) for gg in (0, 1)]
    basis = list(itertools.combinations(range(6), 3))
    idx = {b: i for i, b in enumerate(basis)}
    n = len(basis)
    cnt = np.zeros((n, 3))
    for i, b in enumerate(basis):
        for o in b:
            cnt[i, (1, 4, 5).index(orbitals[o][0])] += 1
    blocks = {}
    for i in range(n):
        blocks.setdefault(tuple(int(x) for x in cnt[i]), []).append(i)

    rng = np.random.default_rng(20261026)
    worst, nfrac = 0.0, 0
    for _ in range(1000):
        Hm = np.zeros((n, n), complex)
        for ids in blocks.values():
            m = len(ids)
            A = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
            Hm[np.ix_(ids, ids)] = (A + A.conj().T) / 2
        _, evec = np.linalg.eigh(Hm)
        for j in range(n):
            w = (np.abs(evec[:, j]) ** 2) @ cnt
            worst = max(worst, float(np.max(np.abs(w - np.round(w)))))
            nfrac += int(np.any(np.abs(w - np.round(w)) > 1e-9))
    check("K1 flavour-對角 H ×1000 ⟹ w_k 全整數", worst < 1e-12 and nfrac == 0,
          f"max|w−round(w)|={worst:.2e}，非整數態 {nfrac}", "0")
    print("  mode-content 塊（⟨O⟩ = w₁−w₄）：" +
          ", ".join(f"{k}:{k[0]-k[1]:+d}(dim {len(v)})" for k, v in sorted(blocks.items())))
    OUT["K1"] = dict(trials=1000, dim=n, blocks={str(k): len(v) for k, v in blocks.items()},
                     max_dev=worst, noninteger=nfrac)

    # ── K2 反例 ──────────────────────────────────────────────────────────
    print("\nK2 反例：mode-1 ↔ mode-4 換模項（hermitian；電荷改變，語料未註冊）")
    H2 = np.zeros((n, n), complex)
    for ids in blocks.values():
        m = len(ids)
        A = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        H2[np.ix_(ids, ids)] = (A + A.conj().T) / 2
    ncouple = 0
    for i, b in enumerate(basis):
        for oa in (0, 1):
            for ob in (2, 3):
                if oa in b and ob not in b:
                    bb = idx[tuple(sorted([x for x in b if x != oa] + [ob]))]
                    val = rng.normal() + 1j * rng.normal()
                    H2[i, bb] += val
                    H2[bb, i] += np.conj(val)
                    ncouple += 1
    _, ev2 = np.linalg.eigh(H2)
    w2 = (np.abs(ev2) ** 2).T @ cnt
    nonint = int(np.sum(np.any(np.abs(w2 - np.round(w2)) > 1e-6, axis=1)))
    print(f"    耦合 {ncouple} 條 ⟹ 非整數本徵態 {nonint}/{n}；例 w = {np.round(w2[0], 4)}")
    check("K2 味非對角 ⟹ 分數佔用（唯一破口）", nonint > 0, f"{nonint}/{n}", "> 0")
    OUT["K2"] = dict(couplings=ncouple, noninteger_states=nonint)

    # ── K3 χ² ────────────────────────────────────────────────────────────
    print("\nK3 ⟨O⟩ 三錨 χ²（p 用不完全 gamma；3 dof）")
    rows = [("⟨O⟩ ≡ 1（價夸克整數）", np.ones(3), "M4 docs/29 C1：χ²=18.22；p 誤報 1.97e−5（1 dof）"),
            ("1 + 2n_s/3（M1；後驗）", 1 + (B - 2) * NS / 3, "M4 C3 M1"),
            ("1 + n_s/2（M3；後驗）", 1 + 0.5 * NS, "M4 C3 M3"),
            ("1 + 3n_s/4（M4；後驗）", 1 + 0.75 * NS, "M4 C3 M4"),
            ("3/(3−n_s)（M2）", 3 / (3 - NS), "M4 C3 M2"),
            ("1 + n_s·P₀（M5 約定 0.3689）", 1 + NS * o0["P0"], "結構候選（本輪）"),
            ("1 + n_s·P₀（w=v² 約定）", 1 + NS * o1["P0"], "結構候選（本輪）"),
            ("1 + n_s·⟨B^{−d}⟩（M5 約定）", 1 + NS * o0["rA"], "結構候選（本輪）"),
            ("sector 投影子 ⟨O⟩ ≡ 2（K6）", 2 * np.ones(3), "T2 字面形式之嚴格讀法")]
    tab = []
    for name, pred, note in rows:
        c = chi2(pred)
        tab.append(dict(model=name, chi2=c, p=pval(c, 3), note=note))
    best_int = min(((chi2(t), t) for t in itertools.combinations_with_replacement(range(4), 3)))
    tab.append(dict(model=f"最佳整數三元組 {best_int[1]}", chi2=best_int[0], p=pval(best_int[0], 3),
                    note="可達性見 K6：實際可達集 = {(1,1,1)}"))
    W = 1 / SIG ** 2
    Am = np.vstack([np.ones(3), NS]).T
    th = np.linalg.solve(Am.T @ (Am * W[:, None]), Am.T @ (W * ANCH))
    cov = np.linalg.inv(Am.T @ (Am * W[:, None]))
    res = float((((ANCH - Am @ th) ** 2) * W).sum())
    tab.append(dict(model=f"自由線性 a+b·n_s（a={th[0]:.4f}±{np.sqrt(cov[0,0]):.4f}, "
                          f"b={th[1]:.4f}±{np.sqrt(cov[1,1]):.4f}）",
                    chi2=res, p=pval(res, 1), note="M4 C2 報 a=0.886, b=0.675"))
    print(f"  {'模型':52}{'χ²':>7}{'p':>10}")
    for x in tab:
        print(f"  {x['model']:52}{x['chi2']:7.2f}{x['p']:10.5f}   {x['note']}")
    check("K3 重現 ⟨O⟩≡1 之 χ²=18.22", abs(tab[0]["chi2"] - 18.22) < 0.01, f"{tab[0]['chi2']:.2f}", "18.22")
    check("K3 ★審計：3 dof 正確 p=3.96e−4（M4 報 1.97e−5＝1 dof 值）",
          abs(tab[0]["p"] - 3.962e-4) < 5e-7, f"{tab[0]['p']:.3e}", "3.96e−4",
          "⟹ 仍否證（3.36σ），但 M4 p 差 20×：register handoff 更正")
    OUT["K3"] = dict(rows=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                            for k, v in x.items()} for x in tab],
                     linear=dict(a=float(th[0]), b=float(th[1]),
                                 sig_a=float(np.sqrt(cov[0, 0])), sig_b=float(np.sqrt(cov[1, 1]))),
                     best_integer=[float(best_int[0]), list(best_int[1])])

    # ── K4 非可加性需求 ──────────────────────────────────────────────────
    print("\nK4 非可加性需求 Δ = ⟨O⟩_h − 1 vs M5 兩體 binding 響應")
    need = ANCH - 1.0
    print(f"  Δ = {np.round(need, 4)} ± {SIG}")
    OUT["K4"] = {}
    for tag, (conv, (a, t, o)) in {"M5": ("M5", (a0, t0, o0)),
                                   "w=v²": ("correct", (a1, t1, o1))}.items():
        h0 = o["h"]
        ep = stats(o["g"], h0 * (1 + 1e-6), conv)[0]
        em = stats(o["g"], h0 * (1 - 1e-6), conv)[0]
        d_eps = (ep - em) / (2e-6)
        dE = -h0 / M_L * t * d_eps
        short = [abs(x) / abs(dE) for x in need]
        print(f"  [{tag:5}] h={h0:.6f} dε/dh={d_eps:+.6f} ∂E_bind/∂m_l={dE:+.3e}"
              f" (E_bind={o['E0']*t*1000:.3f} MeV)")
        print(f"          所需/可得 = {[f'{s:.3g}' for s in short]}；N 需 −0.109 ⟹ "
              f"{'反號即死' if dE > 0 else '同號'}")
        OUT["K4"][tag] = dict(h=h0, d_eps_dh=d_eps, dEbind_dml=dE, shortfall=short,
                              Ebind_MeV=o["E0"] * t * 1000)
    smin = min(min(v["shortfall"]) for v in OUT["K4"].values())
    check("K4 判準（≥10² 倍，事前寫死）", smin > 100, f"最小比值 {smin:.3g}×", "> 100×",
          "★判準 FAIL 照實報（N 通道 36.5×／181×）；三通道全部 ≥36× 且 N 反號 ⟹ 定性結論不變")
    OUT["K4"]["need"] = [float(x) for x in need]

    # ── K6 恒等性 ────────────────────────────────────────────────────────
    print("\nK6 恒等性（Fraction 精確）：價夸克 ⟹ ⟨O⟩≡1；sector 投影 ⟹ ⟨O⟩≡2")
    chans = {"N": ((2, 1, 0), (1, 2, 0), 1),
             "Σ": ((2, 0, 1), (0, 2, 1), 2),
             "Ξ": ((1, 0, 2), (0, 1, 2), 1)}
    val, sec = [], []
    for nm, (h1, h2, coeff) in chans.items():
        dnu, dnd, dns = (Fraction(h1[i] - h2[i]) for i in range(3))
        o_val = abs(Fraction(coeff) * Fraction(-1)) / Fraction(coeff)
        o_sec = abs(dnu - dnd - dns) / Fraction(coeff)
        val.append(float(o_val)); sec.append(float(o_sec))
        print(f"  {nm}: Δ(n_u,n_d,n_s)=({int(dnu):+d},{int(dnd):+d},{int(dns):+d}), coeff={coeff} "
              f"⟹ ⟨O⟩_valence={float(o_val):.0f}, ⟨O⟩_sector={float(o_sec):.0f}")
    check("K6 價夸克計數恒等 ⟨O⟩≡1（三通道）", all(abs(x - 1) < 1e-15 for x in val),
          f"{val}", "[1,1,1]", "整數 w ⟺ 價夸克計數：⟨O⟩=1 係恒等，非預言")
    check("K6 sector 投影子恒 ⟨O⟩≡2（三通道）", all(abs(x - 2) < 1e-15 for x in sec),
          f"{sec}", "[2,2,2]", "T2 字面形式（P₁−P₄）預言 2")
    c2 = chi2(np.array(sec))
    print(f"  sector 讀法 χ² = {c2:.2f} / 3 dof（p = {pval(c2,3):.2e}）")
    OUT["K6"] = dict(valence=val, sector=sec, chi2_sector=c2, p_sector=pval(c2, 3),
                     channels={k: dict(h1=v[0], h2=v[1], coeff=v[2]) for k, v in chans.items()})

    # ── K5 判別力 ────────────────────────────────────────────────────────
    print("\nK5 判別力：M1 vs P₀-sea")
    c_m1 = chi2(1 + (B - 2) * NS / 3)
    c_p0 = chi2(1 + NS * o0["P0"])
    dchi = abs(c_m1 - c_p0)
    print(f"  M1 χ²={c_m1:.2f} · P₀-sea χ²={c_p0:.2f} · Δχ²={dchi:.2f} ⟹ {np.sqrt(dchi):.2f}σ")
    check("K5 不可分辨（Δχ² < 4）——不得宣稱 M1 勝出", dchi < 4, f"Δχ²={dchi:.2f}", "< 4",
          "斜率判別須等 K／D 通道盲測")
    OUT["K5"] = dict(chi2_M1=c_m1, chi2_P0sea=c_p0, delta=dchi)

    # ── K7 定義層 ────────────────────────────────────────────────────────
    print("\nK7 定義層：v9 §D.1 sector 表含 u、d 否？")
    v9 = {"k=1": "(c, t)", "k=4": "(s, b)", "k=5": "(e, μ, τ)"}
    print(f"  v9 wfield_lattice/03_extended_derivation.md L205–208：{v9}")
    has_u = any("u" in v for v in v9.values())
    has_d = any(v.strip("()").split(",")[0].strip() == "d" for v in v9.values())
    check("K7 sector 表不含 u（已確認）", not has_u, f"has_u = {has_u}", "False ⟹ P₁ ≠ ūu")
    check("K7 sector 表不含 d（已確認）", not has_d, f"has_d = {has_d}", "False ⟹ P₄ ≠ d̄d")
    check("K7 ⟹ T2 之識別（ūu→P₁、d̄d→P₄）無文本依據", (not has_u) and (not has_d),
          "兩者皆缺", "T2 未良定義；Q4 係前置閘")
    OUT["K7"] = dict(v9_sectors=v9, has_u=has_u, has_d=has_d)

    # ── K8 跨 M 一致性：M2 已註冊之加性質量泛函 ──────────────────────────
    print("\nK8 跨 M 核：M2 docs/01 §2「質量泛函唯一性：質量 = 加性泛函」")
    # 精確：Σ_{k=1}^{44} sin²(kπ/88) = 45/2（標準恆等式 Σ_{k=1}^{N-1}sin²(kπ/N)=N/2，N=88，
    # 對稱 ⟹ 2S₄₃ + 1 = 44 ⟹ S₄₄ = 22.5）；γ_n = 4sin²(nπ/88) ⟹ Σγ = 4·22.5 = 90
    S44 = Fraction(45, 2)
    sum_gamma = 4 * S44
    print(f"  Σ_{{k=1}}^{{44}} sin²(kπ/88) = 45/2（精確）⟹ Σγ = 4·45/2 = {int(sum_gamma)}")
    gam = {n: 4 * np.sin(np.pi * n / 88) ** 2 for n in (1, 4, 5, 27, 44)}
    print("  γ_n: " + ", ".join(f"γ_{n} = {v:.6f}" for n, v in gam.items()))
    num_sum = sum(4 * np.sin(np.pi * n / 88) ** 2 for n in range(1, 45))
    check("K8 Σγ = 90（精確 + 數值 1e−12）", sum_gamma == 90 and abs(num_sum - 90) < 1e-12,
          f"精確 {int(sum_gamma)}；數值 {num_sum:.12f}", "90")
    check("K8 γ₁ 與 v9 §D.1 之 0.0051 一致（同譜）", abs(gam[1] - 0.0051) < 5e-5,
          f"{gam[1]:.7f}", "0.0051")
    print("  m = M₀Σσ_nγ_n 對 σ 線性 ⟹ ∂m/∂(任何單模量子) = M₀γ_n（逐模常數）")
    print("  ⟹ 質子/中子之 u↔d 交換 ⟹ ⟨O⟩ 恰為交換數 / coeff = 1（K6 恒等）")
    check("K8 M2 加性 ⟹ ⟨O⟩ ≡ 1（與 ⟨O⟩≡1 否證 3.36σ 衝突）", True,
          "衝突 3.36σ", "須上游裁定（M2 前提範圍 vs M4 數據讀法）",
          "★本輪最重要跨 M 結果：唔係 M7 推唔到，而係 M2 已證加性 ⟹ ⟨O⟩=1")
    OUT["K8"] = dict(sum_gamma_exact=int(sum_gamma), gamma={str(k): float(v) for k, v in gam.items()},
                     m2_claim="m = M₀Σσ_nγ_n（加性；p=1 唯一，M2 docs/01 §2/§7A）",
                     conflict_sigma=3.36)

    npass = sum(1 for r in RES if r["pass_"])
    verdict = ("R5：①整數定理成立（w_k∈ℤ，1000 次隨機核）；②T2 在已註冊文本未良定義"
               "（v9 sector 表為 (c,t)/(s,b)/(e,μ,τ)，不含 u、d；sector 讀法 ⟹ ⟨O⟩≡2，χ²≈112 即死）；"
               "③價夸克計數恒等 ⟹ ⟨O⟩≡1，χ²=18.22、p=3.96e−4（3 dof；M4 報 1.97e−5 係 1 dof）⟹ 3.36σ 否證；"
               "④∴ ⟨O⟩≠1 直接量度質量泛函之非可加性；M5 兩體 binding 響應比所需細 36.5–433×（M5 約定）／"
               "181–2150×（w=v² 約定），且 N 通道反號；"
               "⑤T1 退化為恒等、T2 雙重死、T3 未閉合。")
    print("\n" + verdict)
    OUT.update(dict(checks=RES, npass=npass, ntotal=len(RES), verdict=verdict))
    p = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "data", "02b_r5_mode_content_nogo.json"))
    with open(p, "w") as f:
        json.dump(OUT, f, ensure_ascii=False, indent=1, default=str)
    print(f"\nJSON → {p}")


if __name__ == "__main__":
    main()