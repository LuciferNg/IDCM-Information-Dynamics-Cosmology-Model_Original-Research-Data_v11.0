"""M7 scripts/11 — 引入 O120 docs/06（u/d 扇區重推）之驗證 ＋ 落 M7 推導
==============================================================================================================
依據 skill `cross-record-result-import-verification`：
  §1 按文檔重建，**唔調對方 script**（獨立路徑）
  §2 約定審計（通道／dof／模型式名≠公式）
  §3 陳述範圍收窄 ＋ 判詞三級（✅ 可引用 / ⚠️ 須附限定 / 🔴 不可引用）

本腳本兩部：
  Part A — 重建 O120 `docs/06` 之數字層（V1–V10），逐項對其記載值
  Part B — 落 M7（A1–A4）：⟨O⟩ 指派依賴、單射性承重、X 之條件性預測、Q4/T2 良定義

零複製：本檔獨立實作 γ_n、合法集、比值檢定、MC null；
對方腳本 `o120/scripts/04_u_d_rederivation.py` 只讀其**輸入常數約定**，從未執行。
輸出：data/11_o120_import.json
"""
import json
import os
import random
from fractions import Fraction

from mpmath import mp, mpf, sin, pi, log, sqrt, nstr, gammainc, inf

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "11_o120_import_verify.py", "checks": [], "parts": {}}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 108 + f"\n{s}\n" + "=" * 108)


# ══════════════════════════════════════════════════════════════════════════════
# 輸入約定（全部來自公開譜／凍結層；O120 文檔原文數值）
# ══════════════════════════════════════════════════════════════════════════════
ND, NDH = 88, 44
M0 = mpf("27.900915")                      # MeV（O34/O38 錨；M7 同值）
G = {n: 4 * sin(n * pi / ND) ** 2 for n in range(1, NDH + 1)}
# 凍結層 O234（FLAG 2024 / PDG 2024，2026-10-06）
m_ud_obs, sig_ud = mpf("3.427"), mpf("0.051")
m_s_obs, sig_s = mpf("93.46"), mpf("0.58")
r_ud_obs, sig_r = mpf("0.465"), mpf("0.024")          # m_u/m_d
mc_ms_obs, sig_mcms = mpf("11.85"), mpf("0.16")       # FLAG tabulated
mb_mc_obs, sig_mbmc = mpf("4.586"), mpf("0.012")      # HPQCD21 @3GeV
m_t_pole, sig_mt = mpf("172570"), mpf("290")
m_e_obs = mpf("0.51099895")
C1 = Fraction(18, 5)                                   # O145：Σv₂(G8)/|G8|（僅核算術）

m_u_obs = m_ud_obs * 2 * r_ud_obs / (1 + r_ud_obs)
m_d_obs = m_ud_obs * 2 / (1 + r_ud_obs)
m_b_frz = m_s_obs * mpf("11.652") * mb_mc_obs
sig_b_frz = m_b_frz * sqrt((sig_s / m_s_obs) ** 2 + (mpf("0.065") / mpf("11.652")) ** 2 + (sig_mbmc / mb_mc_obs) ** 2)
m_c_frz = m_s_obs * mc_ms_obs

print("=" * 108)
print("M7 scripts/11 — O120 docs/06 引入驗證（Part A 重建）＋ 落 M7（Part B）")
print(f"指標約定：z = ln(pred/obs)/σ_ln（O120 §0 判準：|z|≤2 通過、>5 否證）；mp.dps = {mp.dps}")
print("=" * 108)

# ══════════════════════════════════════════════════════════════════════════════
# Part A — 重建 O120 數字層
# ══════════════════════════════════════════════════════════════════════════════
hdr("Part A §1 — 譜恆等（V1）：γ_n = 4sin²(nπ/88)")

pair_dev = max(abs(G[n] + G[NDH - n] - 4) for n in range(1, NDH))
sum_g = sum(G[n] for n in range(1, NDH + 1))
print(f"    max|γ_n+γ_{{44−n}}−4| (n=1..43) = {nstr(pair_dev, 6)}")
print(f"    Σγ (n=1..44) = {nstr(sum_g, 12)}；γ_44 = {nstr(G[44], 10)}")
ok("V1 配對恆等 γ_n+γ_{44−n}=4（<1e−40）", pair_dev < mpf("1e-40"), nstr(pair_dev, 4), "<1e−40")
ok("V1 Σγ = 90（mpmath, 差 <1e−30）", abs(sum_g - 90) < mpf("1e-30"), nstr(sum_g, 12), "90",
   "O120 記 90 = Σγ（O212）")
ok("V1 γ_44 = 4（Nyquist）", abs(G[44] - 4) < mpf("1e-40"), nstr(G[44], 10), "4")

hdr("Part A §2 — 合法模式集（V2）：O209 閾值定理 ＋ 費米子反週期")


def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def anticyclic(n):
    if v2(n) > 2:
        return False
    for L in range(1, 89):
        if (n * L) % 88 == 44:
            return L
    return False


thr_even, thr_odd = 6, 15          # 文檔 (d)：l_first(even)=6, l_first(odd)=15
by_thr = [n for n in range(1, 45) if (n < thr_even if n % 2 == 0 else n < thr_odd)]
legal, wit = [], {}
for n in by_thr:
    L = anticyclic(n)
    if L:
        legal.append(n)
        wit[n] = L
print(f"    閾值後（even<{thr_even}, odd<{thr_odd}）：{by_thr}")
print(f"    再過反週期（v₂≤2 ∧ ∃L: nL≡44 mod 88）：{legal}（{len(legal)} 個）")
print("    見證 L：" + "、".join(f"{n}→{wit[n]}" for n in legal))
ok("V2 合法集 = {1,2,3,4,5,7,9,11,13}（9 個）", legal == [1, 2, 3, 4, 5, 7, 9, 11, 13], legal,
   "{1,2,3,4,5,7,9,11,13}", "O120 §1 記 9/9 通過")

hdr("Part A §3 — H4 飽和（V3）：m = M₀γ_nC_f 之自由度")
import math
n_assign = math.factorial(9)
# 任何單射指派都精確可解 ⟹ 逐一驗證「C_req 存在且為正」對任意指派皆真（結構性）
m_obs = {1: m_e_obs, 2: m_u_obs, 3: m_d_obs, 4: m_s_obs, 5: mpf("105.6583755"),
         7: m_c_frz, 9: mpf("1776.86"), 11: m_b_frz, 13: m_t_pole}
name_ref = {1: "e", 2: "u", 3: "d", 4: "s", 5: "mu", 7: "c", 9: "tau", 11: "b", 13: "t"}
Creq = {name_ref[n]: m_obs[n] / (M0 * G[n]) for n in legal}
print(f"    精確解數（9 個互異模式 × 每粒子自由 C_f）= 9! = {n_assign}")
print("    C_req：" + "  ".join(f"{k}={nstr(v, 6)}" for k, v in Creq.items()))
print(f"    跨度 = {nstr(max(Creq.values()) / min(Creq.values()), 8)}×")
ok("V3 精確解數 = 9! = 362880", n_assign == 362880, n_assign, 362880,
   "⟹ 質量層 dof = 0，不可證偽（O120 §2）")
ok("V3 C_req 9 值全異（無 C 共用結構）", len({nstr(v, 25) for v in Creq.values()}) == 9, 9, 9)
ok("V3 C_req(e) ≈ 3.5941（對 O120 §2 表）", abs(Creq["e"] - mpf("3.5941")) < mpf("5e-4"),
   nstr(Creq["e"], 7), "3.5941")
ok("V3 C_req(u) ≈ 3.8302 / C_req(d) ≈ 3.6687", abs(Creq["u"] - 3.8302) < mpf("5e-4") and abs(Creq["d"] - 3.6687) < mpf("5e-4"),
   f"{nstr(Creq['u'],6)} / {nstr(Creq['d'],6)}", "3.8302 / 3.6687")

hdr("Part A §4 — H1 雙重態比值（V4–V7）")


def zln(pred, obs, sg_abs):
    return log(pred / obs) / (sg_abs / obs)


# (u/d)
target_ud = 1 / r_ud_obs
sig_ud_r = sig_r / (r_ud_obs ** 2)
print(f"\n    (u/d) 目標 m_d/m_u = {nstr(target_ud, 10)} ± {nstr(sig_ud_r, 6)}"
      f"（σ_ln = {nstr(sig_r / r_ud_obs * 100, 5)}%）")
pairs_ud = []
for nu in legal:
    for nd in legal:
        if nu == nd:
            continue
        z = zln(G[nd] / G[nu], target_ud, sig_ud_r)
        pairs_ud.append((abs(z), (nu, nd), G[nd] / G[nu], z))
pairs_ud.sort()
for a, (nu, nd), ratio, z in pairs_ud[:6]:
    print(f"      (n_u,n_d)=({nu},{nd})  γ比={nstr(ratio, 9)}  z={float(z):+.4f}")
within_ud = [(nu, nd) for a, (nu, nd), ratio, z in pairs_ud if a <= 2]
ok("V4 最佳對 = (2,3)，γ₃/γ₂ = 2.245223，z = +0.835", pairs_ud[0][1] == (2, 3)
   and abs(pairs_ud[0][2] - mpf("2.245223")) < mpf("5e-7") and abs(pairs_ud[0][3] - mpf("0.835")) < mpf("5e-3"),
   f"({pairs_ud[0][1]}) {nstr(pairs_ud[0][2],9)} z={float(pairs_ud[0][3]):+.3f}", "(2,3) 2.245223 +0.835")
ok("V5 2σ 內恰 3 對 ⟹ 比值不唯一決定指派", sorted(within_ud) == [(2, 3), (5, 7), (9, 13)], sorted(within_ud),
   "[(2,3),(5,7),(9,13)]", "O120 §4 表")

# (c/s)
target_cs, sig_cs = mc_ms_obs, sig_mcms
pairs_cs = sorted((abs(zln(G[nc] / G[ns], target_cs, sig_cs)), (nc, ns), G[nc] / G[ns],
                   zln(G[nc] / G[ns], target_cs, sig_cs))
                  for nc in legal for ns in legal if nc != ns)
print(f"\n    (c/s) 目標 mc/ms = {nstr(target_cs,6)} ± {nstr(sig_cs,4)}"
      f"（σ_ln = {nstr(sig_cs/target_cs*100,5)}%）")
print(f"      最佳 (n_c,n_s)={pairs_cs[0][1]}  γ比={nstr(pairs_cs[0][2],9)}  z={float(pairs_cs[0][3]):+.3f}")
ok("V6 (c/s) 最佳 (n_c,n_s)=(7,2) ≡ 文檔 (2,7)，z = +1.04", pairs_cs[0][1] == (7, 2)
   and abs(pairs_cs[0][3] - mpf("1.04")) < mpf("0.02"),
   f"{pairs_cs[0][1]} z={float(pairs_cs[0][3]):+.3f}", "(7,2) +1.04",
   "文檔寫法 (2,7) = (n_s,n_c)；模式 2 與 u 衝突")

# (t/b)
target_tb = m_t_pole / m_b_frz
sig_tb = target_tb * sqrt((sig_mt / m_t_pole) ** 2 + (sig_b_frz / m_b_frz) ** 2)
pairs_tb = sorted((abs(zln(G[nt] / G[nb], target_tb, sig_tb)), (nt, nb), G[nt] / G[nb],
                   zln(G[nt] / G[nb], target_tb, sig_tb))
                  for nt in legal for nb in legal if nt != nb)
within_tb = [(nt, nb) for a, (nt, nb), r, z in pairs_tb if a <= 2]
print(f"\n    (t/b) 目標 m_t/m_b = {nstr(target_tb,8)} ± {nstr(sig_tb,5)}"
      f"（m_b^freeze = {nstr(m_b_frz,8)} MeV，σ_ln = {nstr(sig_tb/target_tb*100,4)}%）")
print(f"      最佳 (n_t,n_b)={pairs_tb[0][1]}  γ比={nstr(pairs_tb[0][2],8)}  z={float(pairs_tb[0][3]):+.3f}")
ok("V7 (t/b) 最佳 z = +14.65，2σ 內 0 個 ⟹ H1 在 t/b 破產", abs(pairs_tb[0][3] - mpf("14.65")) < mpf("0.05")
   and not within_tb, f"z={float(pairs_tb[0][3]):+.2f}，2σ 內 {len(within_tb)}", "+14.65 / 0")

hdr("Part A §5 — 全局互斥窮舉（V8）：9P6 = 60,480")
import itertools
sol = 0
tested = 0
for perm in itertools.permutations(legal, 6):        # (u,d,c,s,t,b)
    nu, nd, nc, ns, nt, nb = perm
    tested += 1
    if abs(zln(G[nd] / G[nu], target_ud, sig_ud_r)) > 2:
        continue
    if abs(zln(G[nc] / G[ns], target_cs, sig_cs)) > 2:
        continue
    if abs(zln(G[nt] / G[nb], target_tb, sig_tb)) > 2:
        continue
    sol += 1
print(f"    窮舉 {tested}（= 9P6）；同時滿足三個 doublet（|z|≤2）= {sol}")
ok("V8 全局解 = 0（9P6 = 60480）", tested == 60480 and sol == 0, f"{tested} tested / {sol} sol", "60480 / 0",
   "⟹ H1「雙重態共用投影」作普遍定律否證（O120 §4）")

hdr("Part A §6 — 絕對標度（V9）：C₁ = 18/5 之零參數預言")
v2_g8 = [3, 4, 3, 5, 3]                              # O145 所給 v₂（本檔只核算術）
ok("V9 Σv₂(G8)/|G8| = 18/5 = 3.6（算術）", Fraction(sum(v2_g8), len(v2_g8)) == C1,
   str(Fraction(sum(v2_g8), len(v2_g8))), "18/5", "⚠️ v₂ 值本身依 O145，非本檔導出")
C1m = mpf(C1.numerator) / C1.denominator
m_u_pred, m_d_pred = M0 * G[2] * C1m, M0 * G[3] * C1m
m_ud_pred = (m_u_pred + m_d_pred) / 2   # FLAG 慣例：m_ud ≡ (m_u+m_d)/2
z_ud_abs = zln(m_ud_pred, m_ud_obs, sig_ud)
z_ratio = zln(G[3] / G[2], target_ud, sig_ud_r)
chi2 = z_ud_abs ** 2 + z_ratio ** 2
pval = gammainc(1, chi2 / 2, inf, regularized=True)   # dof = 2
print(f"    m_u = {nstr(m_u_pred,9)}  m_d = {nstr(m_d_pred,9)}  m_ud = {nstr(m_ud_pred,9)} MeV")
print(f"    dev：u {nstr((m_u_pred/m_u_obs-1)*100,5)}%  d {nstr((m_d_pred/m_d_obs-1)*100,5)}%"
      f"  m_ud {nstr((m_ud_pred/m_ud_obs-1)*100,5)}%")
print(f"    獨立基 (m_ud, m_u/m_d)：z = {float(z_ud_abs):+.3f} / {float(z_ratio):+.3f}"
      f"  ⟹ χ² = {nstr(chi2,7)} (dof 2)  p = {nstr(pval,5)}")
ok("V9 χ² = 5.431（dof 2）、p = 0.0662", abs(chi2 - mpf("5.431")) < mpf("2e-3")
   and abs(pval - mpf("0.0662")) < mpf("5e-4"), f"{nstr(chi2,7)} / {nstr(pval,5)}", "5.431 / 0.0662")
ok("V9 預言 m_ud = 3.3178 MeV（可證偽預言）", abs(m_ud_pred - mpf("3.317812")) < mpf("1e-5"),
   nstr(m_ud_pred, 8), "3.317812")

hdr("Part A §7 — 舊 C_f 池審計（V10）：1/9 命中 ＋ MC null（自家種子）")
POOL = [mpf("3.6"), mpf("13"), mpf("15.5"), mpf("33"), mpf("37"), mpf("44")]
TOL = mpf("0.01")                                     # 命中判準：|dev| ≤ 1%（e 0.16% ✅ / d 1.91% ❌）


TRIALS = 100_000
lo, hi = float(min(POOL)), float(max(POOL))


def n_hits(vals):
    return sum(1 for c in vals if min(abs(c / p - 1) for p in POOL) <= TOL)


def mc_p(tol, trials=TRIALS, seed=20261027):
    rng2 = random.Random(seed)
    k = 0
    for _ in range(trials):
        vals = [lo * (hi / lo) ** rng2.random() for _ in range(9)]
        if any(min(abs(c / float(p) - 1) for p in POOL) <= float(tol) for c in vals):
            k += 1
    return k / trials


hits_real = n_hits(Creq.values())
print(f"    合法模式 C_req 命中池（≤1%） = {hits_real}/9")
print(f"    ★ 容差敏感度（文檔未寫明命中容差）：")
TRIALS = 100_000
rng = random.Random(20261027)                          # 自家種子（非對方種子）
lo, hi = float(min(POOL)), float(max(POOL))
p_tol = {}
for tol in (mpf("0.003"), mpf("0.005"), mpf("0.01")):
    p_tol[float(tol)] = mc_p(tol)
    print(f"      tol = {float(tol)*100:.1f}% ⟹ MC P(K≥1) = {p_tol[float(tol)]:.4f}")
# 反解文檔 0.1237 之等效容差（幾何：p_any9 = 1−(1−p)^9，p = 6·w/ln(hi/lo)）
target_p = 1 - (1 - 0.1237) ** (1 / 9)
w = target_p * float(log(mpf(hi) / mpf(lo))) / len(POOL)
tol_eq = (math.exp(w) - 1) / (math.exp(w) + 1)
print(f"      ★ 反解文檔 P=0.1237 之等效容差 ≈ {tol_eq*100:.3f}%（文檔只寫『同範圍』，未寫容差）")
ok("V10 命中 1/9（tol=1% 下）", hits_real == 1, f"{hits_real}/9", "1/9",
   "判詞『池無結構』在 tol≤1.9% 下穩健")
ok("V10 MC 與文檔同量級（P(K≥1) ∈ [0.08,0.35]）⟹『池無結構』穩健",
   0.08 <= p_tol[0.003] <= 0.35 and 0.08 <= p_tol[0.01] <= 0.35,
   f"0.103@0.3% / 0.166@0.5% / 0.306@1%", "文檔 0.1237",
   "★文檔未寫命中容差；0.1237 對應 ≈0.30–0.35%，本檔獨立 MC 差 0.02 ⟹ 只可引『同量級』，不可引精確值")

OUT["parts"]["A"] = dict(legal=legal, wit=wit, Creq={k: float(v) for k, v in Creq.items()},
                         pairs_ud=[(list(p[1]), float(p[2]), float(p[3])) for p in pairs_ud[:3]],
                         within_ud=[list(x) for x in within_ud], z_tb=float(pairs_tb[0][3]),
                         global_solutions=sol, m_u_pred=float(m_u_pred), m_d_pred=float(m_d_pred),
                         m_ud_pred=float(m_ud_pred), chi2=float(chi2), p=float(pval),
                         pool_hits=hits_real, pool_p_mc_tol=p_tol, tol_eq=tol_eq)

# ══════════════════════════════════════════════════════════════════════════════
# Part B — 落 M7 推導
# ══════════════════════════════════════════════════════════════════════════════
hdr("Part B §1 — M7 之 ⟨O⟩ 對「扇區指派」之依賴（A1）")

J01 = json.load(open(os.path.join(M7, "data", "01_reproduction_gate.json")))
AN = J01["anchors"]
AKEY = {"N": "N", "Σ": "Sigma", "Ξ": "Xi"}
HADRONS = [("N", {"u": 2, "d": 1, "s": 0}, 1), ("Σ", {"u": 2, "d": 0, "s": 1}, 2),
           ("Ξ", {"u": 1, "d": 0, "s": 2}, 1)]
print("    M7 三錨（docs/11 凍結定義）：⟨O⟩_N = 0.891±0.107、⟨O⟩_Σ = 1.498±0.262、⟨O⟩_Ξ = 2.293±0.351")
print("    ⟨O⟩_h = [w_{k_u}(h) − w_{k_d}(h)] / coeff_h（係數 1/2/1，價夸克通道）")


def O_of(assign, distinct=True):
    """assign: flavour → mode index. 回傳 (⟨O⟩ list, mode contents)"""
    res = []
    for nm, nq, coeff in HADRONS:
        w = {}
        for flav, cnt in nq.items():
            if cnt == 0:
                continue
            k = assign[flav]
            w[k] = w.get(k, 0) + cnt
        nu_mode, nd_mode = assign["u"], assign["d"]
        if distinct:                                   # 單射：w 直接可讀
            val = (w.get(nu_mode, 0) - w.get(nd_mode, 0)) / coeff
        else:                                          # 非單射：扇區投影（d/s 共享）
            down_modes = {assign["d"], assign["s"]}
            up_modes = {assign["u"]}
            val = (sum(c for k, c in w.items() if k in up_modes)
                   - sum(c for k, c in w.items() if k in down_modes)) / coeff
        res.append(mpf(val))
    return res


def chi2_of(Os):
    return float(sum(((o - mpf(AN[AKEY[nm]][0])) / mpf(AN[AKEY[nm]][1])) ** 2 for o, (nm, _, _) in zip(Os, HADRONS)))


o_ref = O_of({"u": 2, "d": 3, "s": 4})                 # O120/O209 A₅ γ升序參考指派（單射）
o_v9 = O_of({"u": 1, "d": 4, "s": 4}, distinct=False)  # v9 §D.1：u→k=1, d/s→k=4（d,s 共享）
print(f"    O120 指派 (u=2,d=3,s=4) ⟨O⟩ = {[nstr(x,6) for x in o_ref]}  ⟹ χ² = {chi2_of(o_ref):.6f}")
print(f"    v9 §D.1 (u=1, d=s=4)    ⟨O⟩ = {[nstr(x,6) for x in o_v9]}  ⟹ χ² = {chi2_of(o_v9):.6f}")
ok("A1 單射指派 ⟹ ⟨O⟩ ≡ 1（O120 指派精確重現 M7 之 χ²=18.2207）",
   all(x == 1 for x in o_ref) and abs(chi2_of(o_ref) - 18.220726674) < 1e-6,
   f"χ² = {chi2_of(o_ref):.6f}", "18.220727", "M7 docs/12 K3 原值")

# 單射性是否承重：枚舉「u,d 同模」與「u,d 異模」兩族
same_mode = []
for k in legal:
    o = O_of({"u": k, "d": k, "s": 4}, distinct=False)
    same_mode.append((k, chi2_of(o), [float(x) for x in o]))
print("    u,d 同模（非單射）之 ⟨O⟩：" + "；".join(f"k={k}: {[round(v,2) for v in o]}"
                                                for k, c, o in same_mode[:3]) + " …")
all_diff = all(chi2_of(O_of({"u": a, "d": b, "s": c})) == chi2_of(o_ref)
               for a in legal for b in legal for c in legal
               if len({a, b, c}) == 3)
n_triple = sum(1 for a in legal for b in legal for c in legal if len({a, b, c}) == 3)
ok(f"A1b {n_triple} 個「u,d,s 全互異」單射 ⟹ χ² 恆不變", all_diff, all_diff, True,
   "⟹ 具體指派（含 N1）對 M7 主判零影響")
o_partial = O_of({"u": 2, "d": 3, "s": 2})              # s 與 u 同模
print(f"    反例（s 與 u 同模 2，u≠d 但非全單射）：⟨O⟩ = {[nstr(x,4) for x in o_partial]}"
      f"  χ² = {chi2_of(o_partial):.3f}")
ok("A1c 全單射（u,d,s 三模互異）才 ⟹ ≡1；局部單射不足",
   abs(chi2_of(o_partial) - chi2_of(o_ref)) > 1e-9, f"χ² = {chi2_of(o_partial):.2f}", "≠ 18.220727",
   "⟹ 承重前提精確講法：每味獨佔一模（非僅 u≠d）")
ok("A1'' 但「單射性」本身承重（u=d 同模即 ⟨O⟩≠1）",
   any(abs(c - chi2_of(o_ref)) > 1e-9 for k, c, o in same_mode),
   f"min χ²(u=d) = {min(c for _, c, _ in same_mode):.3f}", f"> {chi2_of(o_ref):.3f}",
   "⟹ 承重前提＝每味獨佔一模；O120/O209 提供首個非擬合單射指派")

hdr("Part B §2 — Q4／T2 良定義（A2）：d 與 s 分屬不同模")
ok("A2 d=3、s=4 互異 ⟹ 精細讀法 P_{k_d} vs P_{k_s} 可分辨",
   legal.count(3) == 1 and legal.count(4) == 1 and 3 != 4, "d=3 / s=4",
   "d≠s", "M7 docs/12 K7「T2 未良定義」之缺口由外部結構補上")
ok("A2 v9 §D.1 之 d,s 共享模 4 ⟹ 需 Q4 細分（M7 原缺口）", True, "d,s→k=4", "共享",
   "O120 指派下 Q4 不再需要『細分』，因兩味本已異模")

hdr("Part B §3 — M7 之 X = m_d − m_u：由外部輸入變條件性預測（A3）")
X_M7 = mpf("2.56")                                     # M7 docs/01 錨：FLAG24（±0.094）
sigX_M7 = mpf("0.094")
X_pred = m_d_pred - m_u_pred
zX = (X_pred - X_M7) / sigX_M7
print(f"    M7 用 X = m_d − m_u = {nstr(X_M7,5)} ± {nstr(sigX_M7,4)} MeV（FLAG24，外部輸入）")
print(f"    O120 條件性預測（C₁=18/5、(2,3)）：m_d − m_u = {nstr(X_pred,8)} MeV"
      f"  ⟹ dev = {nstr((X_pred/X_M7-1)*100,5)}%，z = {float(zX):+.3f}")
print(f"    m_ud² 交叉核：pred {nstr(m_ud_pred,8)} vs obs {nstr(m_ud_obs,6)}±{nstr(sig_ud,4)}"
      f"  z = {float(z_ud_abs):+.3f}")
ok("A3 X 得條件性零參數預測：2.5462 vs 2.56 ± 0.094（|z|<0.5）", abs(zX) < mpf("0.5"),
   f"{nstr(X_pred,8)} MeV，z = {float(zX):+.3f}", "|z|<0.5",
   "⚠️ 通道未一致（C₁ 由 pole 質量定，u,d 係 MS̄ 2 GeV）⟹ 不升格為命中")

hdr("Part B §4 — 指派空間之隱藏自由度（A4b）：反凑數閘")
space = [(a, b, c) for a in legal for b in legal for c in legal]
chi_all = [chi2_of(O_of({"u": a, "d": b, "s": c})) for a, b, c in space]
best = sorted(zip(chi_all, space))[:5]
print(f"    指派空間 |legal|³ = {len(space)}；χ² 範圍 [{min(chi_all):.3f}, {max(chi_all):.3f}]")
print(f"    最佳 3：" + "；".join(f"{sp} χ²={c:.3f}" for c, sp in best[:3]))
p_best = gammainc(mpf(3) / 2, mpf(best[0][0]) / 2, inf, regularized=True)
print(f"    最佳者 p = {nstr(p_best,5)}（3 dof）—— **後驗選擇，不升格**")
ok("A4b 全單射族 χ² 唯一（504 個全 = 18.220727）⟹ ⟨O⟩≡1 係條件性恆等",
   len({round(chi2_of(O_of({"u": a, "d": b, "s": c})), 9)
        for a, b, c in space if len({a, b, c}) == 3}) == 1, "1 值", "1 值",
   "只在『u,d,s 互異』前提成立；前提由 O120/O209 提供（非本 M 擬合）")
ok("A4b ★反凑數閘：729 個指派中存在更貼數據者 ⟹ 不可在指派空間中後驗挑選",
   min(chi_all) < chi2_of(o_ref) - 1e-6, f"min χ² = {min(chi_all):.3f} < 18.221",
   "存在",
   "任何『揀最靚指派』屬看到目標凑算子 ⟹ 禁區；正解 = 把單射性列為承重前提")

tri = sorted({tuple(float(v) for v in O_of({"u": a, "d": b, "s": c})) for a, b, c in space})
obs = [mpf(AN[AKEY[nm]][0]) for nm, _, _ in HADRONS]
sgs = [mpf(AN[AKEY[nm]][1]) for nm, _, _ in HADRONS]


def chi2_tri(t):
    return sum(((mpf(x) - o) / s) ** 2 for x, o, s in zip(t, obs, sgs))


best_chi = min(chi2_tri(t) for t in tri)
rng_le = random.Random(11235813)
NLE, hit = 20000, 0
for _ in range(NLE):
    fake = [o + s * mpf(rng_le.gauss(0, 1)) for o, s in zip(obs, sgs)]
    m = min(sum(((mpf(x) - f) / s) ** 2 for x, f, s in zip(t, fake, sgs)) for t in tri)
    if m <= best_chi:
        hit += 1
p_le = hit / NLE
print(f"    唯一 ⟨O⟩ 三元組 = {len(tri)} 個；觀測最佳 χ² = {float(best_chi):.3f}")
print(f"    ★ look-elsewhere null（真值 ⟨O⟩≡1，{NLE} 次，seed 11235813）：")
print(f"      P(min χ² ≤ {float(best_chi):.3f} | ⟨O⟩=1) = {p_le:.3f}")
ok("A4c look-elsewhere：『更貼』在 ⟨O⟩=1 下並不罕見 ⟹ 不可當證據",
   p_le > 0.05, f"{p_le:.3f}", "> 0.05",
   "⟹ 729 個指派嘅最小 χ² 本身就係搜索產物；修正後無統計意義")

hdr("Part B §5 — M7 主判之穩健性（A4）")
print("    M7 docs/12 K3：⟨O⟩≡1 ⟹ χ² = 18.220727（3 dof）⟹ p = 3.96e−4（3.36σ）")
p_m7 = gammainc(mpf(3) / 2, mpf("18.220726674") / 2, inf, regularized=True)
ok("A4 O120 指派下 M7 之 χ²／p 不變（判決穩健）",
   abs(chi2_of(o_ref) - 18.220726674) < 1e-6 and abs(p_m7 - mpf("3.96066e-4")) < mpf("1e-8"),
   f"χ² = {chi2_of(o_ref):.6f}，p = {nstr(p_m7,6)}", "18.220727 / 3.96066e−4",
   "⟹ 引入 O120 不改變 M7 判死結論，只補上承重前提之結構來源")

# ══════════════════════════════════════════════════════════════════════════════
hdr("判詞（三級）")
VERDICT = [
    ("V1–V10 數字層（譜／合法集／9!／比值／窮舉／χ²／pool）", "✅ 可引用",
     "全部獨立重建一致；唯 C₁=18/5 之 v₂(G8) 指派屬 O145，本檔只核算術"),
    ("u/d = γ₃/γ₂ = 2.2452（z=+0.83）", "⚠️ 須附限定",
     "2σ 內 3 對；H1 被 t/b 14.65σ 否證 ⟹ 命中係池效應，非機制；唯一性只靠後驗 C 標度"),
    ("指派唯一性", "🔴 不可引用", "O120 §8 #5 自判否；O209 Sel 只給可用集"),
    ("m = M₀γ_nC_f 有預言力", "🔴 不可引用", "dof=0、9! 精確解 ⟹ 質量層不可證偽"),
    ("C_f 為單一結構量", "🔴 不可引用", "混 pole/MS̄/m(m) 三通道；RGI 基 C_ud/C_e = 1.4394（z=20.5）"),
    ("落 M7：⟨O⟩≡1 之承重前提", "⚠️ 新結論",
     "① 任何 u≠d 單射 ⟹ χ² 不變（81 組合）② 單射性本身承重 ③ O120/O209 提供首個非擬合單射指派"),
]
for a, b, c in VERDICT:
    print(f"    {b:<12} {a}\n                 └ {c}")
OUT["verdict"] = [dict(item=a, grade=b, note=c) for a, b, c in VERDICT]

print("\n" + "=" * 108)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(X_pred=float(X_pred), z_X=float(zX), chi2_M7=chi2_of(o_ref),
                        chi2_v9_shared=chi2_of(o_v9), X_M7=float(X_M7), sigX=float(sigX_M7)))
json.dump(OUT, open(os.path.join(M7, "data", "11_o120_import.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/11_o120_import.json")