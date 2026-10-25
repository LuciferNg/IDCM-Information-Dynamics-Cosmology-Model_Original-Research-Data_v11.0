"""M4 引用前獨立驗證：M5 可搬運結論逐項重算（2026-10-25 調律者指示「先讀取並驗證過正確才引用」）
=================================================================================================
方法：**不調用 M5 任何腳本**，只按其文檔陳述之模型重建（對稱三對角 41×41、B=4、numpy eigh），
      逐項比對 M5 記載值；凡不能重現者列為 🔴 並記明差異來源。

驗證清單（M4 計劃引用者）：
  V1  能帶底：bulk 譜 [−2√(B−1), 2√(B−1)]·t（B=4 ⟹ ±2√3 t）
  V2  docs/13 重校點：由 (A_D, |c|_D) = (1.2025, 0.3525) 獨立解 (a, t) ⟹ a₀、t、⟨B^{−d}⟩、P₀
      ＋ 回代其自身公式 E(D⁺−D⁰) = A+3|c|、E(D*⁺−D*⁰) = A−|c|
  V3  C27「r 只依 m_h·m_l」之**精確範圍**：解析式 r = (6/π)(m_h m_l a²/(ħc)²)(rA/P₀)
  V4  幾何符號定理算術：B^7 = 16384、比值集 {∞, 16384, 1} 全部 ≥ 1
  V5  封頂不等式：⟨B^{−d}⟩ ≤ 1（逐元素 B^{−d} ≤ 1 ＋ 隨機態抽樣）⟹ A ≤ α q ħc/a
  V6  束縛判別：能帶底判據 ＋ **L 收斂判據**（⟨B⟩、P₀ 隨 L 是否收斂）於兩個校準點
  V7  docs/09 接觸半徑之**約定**：r = (π|ψ(0)|²)^{−1/3}（非 3/(4π) top-hat）

單位紀律（本輪揪出之陷阱）：M5 腳本以 ħc = 0.1973 GeV·fm 計 ⟹ t 以 **GeV** 入式，
  文檔寫 0.7884 MeV ⟺ 7.884e−4 GeV。混用會令 t 差 1000×。
輸出：`data/25_import_verification.json`
"""
import json

import numpy as np

B, L = 4, 40
HBARC, ALPHA = 0.1973269804, 1 / 137.035999084          # GeV·fm
MC, MS, ML, QC, QS = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
EBAND = -2 * np.sqrt(B - 1)                              # 單位 t
RES = []


def check(name, ok, got, want=""):
    RES.append({"item": name, "pass": bool(ok), "got": got, "doc": want})
    print(f"  {'✅' if ok else '🔴'} {name}: {got}" + (f"   （M5 文檔 {want}）" if want else ""))


def build(Lx, a, t, q, mh, ml):
    nodes = np.array([1.0] + [B * (B - 1) ** (d - 1) for d in range(1, Lx + 1)])
    powb = np.array([B ** (-d) for d in range(Lx + 1)])
    off = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (Lx - 1))
    g = ALPHA * q * HBARC / a / t                        # t 以 GeV
    cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = cq * (HBARC / a) ** 3 / t
    dg = -g * powb.copy()
    dg[0] -= h
    return np.diag(dg) + np.diag(off, 1) + np.diag(off, -1), nodes, powb


def shape(a, t, q, mh, ml, Lx=L):
    H, nodes, powb = build(Lx, a, t, q, mh, ml)
    ev, evec = np.linalg.eigh(H)
    w = nodes * evec[:, 0] ** 2
    w /= w.sum()
    rA, p0 = float((w * powb).sum()), float(w[0])
    cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    return dict(E0=float(ev[0]), rA=rA, P0=p0, d_peak=int(w.argmax()),
                A=ALPHA * q * (HBARC / a) * rA * 1000, K=cq * (HBARC / a) ** 3 * p0 * 1000,
                w_profile=[float(x) for x in w[:6]])


print("=" * 96)
print("V1 bulk 能帶底（V = 0，均勻 off-diag = −√3，L = 2000）")
n = 2000
offu = np.full(n - 1, -np.sqrt(B - 1))
evb = np.linalg.eigvalsh(np.diag(offu, 1) + np.diag(offu, -1))
check("bulk 下緣 = −2√3", abs(evb[0] + 2 * np.sqrt(3)) < 1e-3, f"{evb[0]:.9f}", f"{-2*np.sqrt(3):.9f}")

print("\nV2 docs/13 重校點（獨立二分；A_D = 1.2025、|c|_D = 0.3525、q = 2/3、m_c = 1.27、m_l = 0.30）")


def t_for_a(a):
    lo, hi = 1e-6, 0.08                                  # GeV
    for _ in range(300):
        m_ = 0.5 * (lo + hi)
        if shape(a, m_, QC, MC, ML)["A"] > 1.2025:
            lo = m_
        else:
            hi = m_
    return 0.5 * (lo + hi)


lo, hi = 0.02, 3.0
for _ in range(300):
    m_ = 0.5 * (lo + hi)
    if shape(m_, t_for_a(m_), QC, MC, ML)["K"] > 0.3525:
        lo = m_
    else:
        hi = m_
a0 = 0.5 * (lo + hi)
t0 = t_for_a(a0)
r0 = shape(a0, t0, QC, MC, ML)
check("a₀ [fm]", abs(a0 - 0.377422) < 2e-5, f"{a0:.6f}", "0.377422")
check("t [MeV]（GeV 入式）", abs(t0 * 1000 - 0.7884) < 5e-4, f"{t0*1000:.4f} MeV = {t0:.3e} GeV", "0.7884")
check("⟨B^{−d}⟩", abs(r0["rA"] - 0.4728) < 5e-4, f"{r0['rA']:.4f}", "0.4728")
check("P₀", abs(r0["P0"] - 0.3689) < 5e-4, f"{r0['P0']:.4f}", "0.3689")
check("E(D⁺−D⁰) = A + 3|c|", abs(r0["A"] + 3 * r0["K"] - 2.2600) < 5e-4, f"{r0['A']+3*r0['K']:.4f}", "2.2600")
check("E(D*⁺−D*⁰) = A − |c|", abs(r0["A"] - r0["K"] - 0.8500) < 5e-4, f"{r0['A']-r0['K']:.4f}", "0.8500")

print("\nV3 C27 之精確範圍：r = (6/π)(m_h m_l a²/(ħc)²)(rA/P₀)")
s1 = shape(a0, t0, QC, MS, ML)
s2 = shape(a0, t0, QC, 0.594, 0.25)                       # 同乘積 0.1485
rf = lambda mh, ml, s: (6 / np.pi) * (mh * ml * a0 ** 2 / HBARC ** 2) * (s["rA"] / s["P0"])
check("解析式 = A/|c|（s_ud）", abs(rf(MS, ML, s1) - s1["A"] / s1["K"]) < 1e-9 * abs(s1["A"] / s1["K"]),
      f"{rf(MS,ML,s1):.9f} vs {s1['A']/s1['K']:.9f}")
check("乘積守恆 ⟹ r 相同（固定 a,t,q）", abs(s1["A"] / s1["K"] - s2["A"] / s2["K"]) < 1e-9,
      f"{s1['A']/s1['K']:.9f} vs {s2['A']/s2['K']:.9f}", "13.400680（docs/23，於其自身 (a,t)）")
q_sw = shape(a0, t0, QS, MS, ML)
print(f"     ⚠️ 範圍限制：固定 (a,t)、改 q（2/3→1/3）⟹ r = {q_sw['A']/q_sw['K']:.6f}（≠ 原值）")
print(f"        ⟹「與電荷無關」只成立於**乘子層**（α、q 精確相消）；形狀 rA/P₀ 仍依 g = αqħc/(at)")

print("\nV4 幾何符號定理算術")
check("B^7", B ** 7 == 16384, f"{B**7}", "16384")
check("比值集全部 ≥ 1", all(x >= 1 for x in (np.inf, 16384.0, 1.0)), "{∞, 16384, 1}")

print("\nV5 封頂不等式 ⟨B^{−d}⟩ ≤ 1")
powb = np.array([B ** (-d) for d in range(L + 1)])
rng = np.random.default_rng(7)
nodes = np.array([1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)])
rand_max = 0.0
for _ in range(20000):
    v = rng.normal(size=L + 1)
    w = nodes * v ** 2
    w /= w.sum()
    rand_max = max(rand_max, float((w * powb).sum()))
check("逐元素 B^{−d} ≤ 1", bool((powb <= 1).all()), f"max = {powb.max()}")
check("隨機態 max⟨B⟩ ≤ 1", rand_max <= 1.0, f"{rand_max:.6f}")
amax = ALPHA * QS * (HBARC / a0) * 1000
check("A_max(q=1/3, a₀) [MeV]", abs(amax - 1.2718) < 2e-3, f"{amax:.4f}", "1.2718")

print("\nV6 束縛判別（能帶底 ＋ L 收斂）於兩校準點")
pts = [("κ=0 主校準（docs/13）", a0, t0), ("κ=0.215 天花板點（docs/19/21/23）", 0.243274, 2.5743e-3)]
v6 = {}
for tag, a, t in pts:
    print(f"  ── {tag}: a = {a} fm, t = {t*1000:.4f} MeV")
    for nm, (q, mh, ml) in {"s_ud": (QS, MS, ML), "c_ud": (QC, MC, ML), "c_s": (QC, MC, MS)}.items():
        s40 = shape(a, t, q, mh, ml, 40)
        s160 = shape(a, t, q, mh, ml, 160)
        conv = abs(s40["rA"] - s160["rA"]) / max(s160["rA"], 1e-30) < 1e-3
        v6[f"{tag}|{nm}"] = {"E0_over_t": s40["E0"], "below_band": bool(s40["E0"] < EBAND),
                             "rA_L40": s40["rA"], "rA_L160": s160["rA"], "converged": bool(conv),
                             "d_peak": s40["d_peak"], "P0": s40["P0"]}
        print(f"     {nm:5} E₀/t={s40['E0']:+.5f}（能帶底 {EBAND:+.5f}）⟹ {'低於能帶' if s40['E0']<EBAND else '帶內'}"
              f" | ⟨B⟩: L40 {s40['rA']:.5f} → L160 {s160['rA']:.5f} ⟹ {'收斂' if conv else '🔴未收斂'}"
              f" | d_peak {s40['d_peak']:2d} | P₀ {s40['P0']:.5f}")

print("\nV7 docs/09 接觸半徑之約定")
psi2 = 0.07909
r_doc = (1 / (np.pi * psi2)) ** (1 / 3) * HBARC
r_alt = (3 / (4 * np.pi * psi2)) ** (1 / 3) * HBARC
check("r = (π|ψ(0)|²)^{−1/3} [fm]", abs(r_doc - 0.314) < 5e-3, f"{r_doc:.4f}", "0.314")
print(f"     ⚠️ 若改用 top-hat (3/4π)^{{1/3}} ⟹ r = {r_alt:.4f} fm ⟹ 引用時必須連約定一齊寫")

print("\n" + "=" * 96)
npass = sum(1 for r in RES if r["pass"])
print(f"合計 {npass}/{len(RES)} PASS")
for r in RES:
    if not r["pass"]:
        print(f"  🔴 FAIL {r['item']}: got {r['got']} want {r['doc']}")
json.dump({"checks": RES, "V6_binding": v6,
           "constants": {"a0_fm": a0, "t_GeV": t0, "t_MeV": t0 * 1000, "rA": r0["rA"], "P0": r0["P0"],
                         "A_D": r0["A"], "K_D": r0["K"], "E_band_over_t": EBAND}},
          open("data/25_import_verification.json", "w"), ensure_ascii=False, indent=1)
print("→ data/25_import_verification.json")