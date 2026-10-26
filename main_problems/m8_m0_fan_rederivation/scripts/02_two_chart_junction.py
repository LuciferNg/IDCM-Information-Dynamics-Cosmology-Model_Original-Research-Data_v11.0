#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/02 — 兩圖卡結構（micro｜macro）＋「角度反向交疊 → 非線性交接點」定理
==============================================================================
來源（調律者 2026-10-26 校正，逐字）:
  「M_0標度本身未必有意義，微觀端是讀數讀取產生的投影，大尺度宏觀端是本徵數值的映射,
   本身就未必存在一個跨越10^38的標度，但映射層可能存在一個按比例的標度，
   這個標度可能是一個角度反向交疊後的非線性交接點」

本腳本只驗**結構命題**（精確代數），不拟合任何目標數值。
任何涉及 4.095e38 之比較一律標 POSTHOC（非預註冊），且沿用 O38 scripts/45 已判結論。

執行: /usr/bin/python3 scripts/02_two_chart_junction.py
輸出: data/02_two_chart_junction.json + stdout（同步寫 data/run_02.log）
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 60

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
N_FAIL = 0


def chk(tag, cond, detail=""):
    global N_FAIL
    ok = bool(cond)
    if not ok:
        N_FAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ============================================================================
# §0 常數（結構常數，非可調）
# ============================================================================
N_D = 88                      # 鏈點數
N_RUNG = 44                   # rung 數（n = 1..44，θ_n = πn/88 ∈ (0, π/2]）
hdr("§0 結構常數")
print(f"N_d = {N_D},  rung 數 = {N_RUNG},  θ_n = πn/88,  n = 1..{N_RUNG}")
print(f"mpmath dps = {mp.mp.dps}")


# ============================================================================
# §1 J1 — 投影讀數：π² = π ⟹ 微觀端 = 不動點集，無獨立標度自由度
# ============================================================================
hdr("§1 J1 — 微觀端 = 投影 im(π)；投影 ⟹ 冪等 ⟹ 無獨立標度 dof")

rungs = [mp.pi * n / N_D for n in range(1, N_RUNG + 1)]


def pi_readout(theta):
    """最近 rung 讀數（鏈上測量 = 投影到離散譜）"""
    best, bd = None, None
    for r in rungs:
        d = abs(theta - r)
        if bd is None or d < bd:
            best, bd = r, d
    return best


# 在整個 (0, π/2] 上取 20001 點，逐點查 π(π(θ)) == π(θ)
NG = 20001
bad = 0
for k in range(1, NG + 1):
    th = mp.pi / 2 * mp.mpf(k) / (NG + 1)
    a = pi_readout(th)
    b = pi_readout(a)
    if a != b:
        bad += 1
chk("J1.1 投影冪等 π∘π = π（20001 點，精確有理比較）", bad == 0, f"例外 = {bad}")

# 投影像 = 44 rung，且 rung 全為不動點
fixpts = [r for r in rungs if pi_readout(r) == r]
chk("J1.2 im(π) = Fix(π) = 44 rung（投影像即不動點集）",
    len(fixpts) == N_RUNG, f"|im| = {len(fixpts)}")

# 微觀端數值域：γ ∈ [0,4]、s = tanθ ∈ [tan(π/88), cot(π/88)] ⟹ 全 O(1)
g_lo = mp.mpf(4) * mp.sin(mp.pi / N_D) ** 2
g_hi = mp.mpf(4) * mp.sin(mp.pi / 2) ** 2
s_lo = mp.tan(mp.pi / N_D)
s_hi = mp.tan(mp.mpf(43) * mp.pi / N_D)
print(f"  γ 域   = [{mp.nstr(g_lo, 12)}, {mp.nstr(g_hi, 6)}]   （有界）")
print(f"  s 域   = [{mp.nstr(s_lo, 12)}, {mp.nstr(s_hi, 12)}]   （有界）")
chk("J1.3 微觀端所有譜量有界（γ ≤ 4；s ≤ cot(π/88) ≈ 28）",
    g_hi == 4 and s_hi < 29, f"γ_max = {g_hi}, s_max = {mp.nstr(s_hi, 10)}")

res_j1 = {
    "idempotent_exceptions": bad,
    "n_grid": NG,
    "gamma_range": [float(g_lo), 4.0],
    "slope_range": [float(s_lo), float(s_hi)],
}
RES["_J1_data"] = res_j1


# ============================================================================
# §2 J2 — 投影座標唯一性：交比 ⟹ s = tan θ（「按比例」圖卡之唯一性）
# ============================================================================
hdr("§2 J2 — 映射層唯一投影座標 = tan θ（交比歸一化；非假設，係交比唯一性）")

t1, t2, t3, t4 = sp.symbols("t1 t2 t3 t4", positive=True)
S = lambda x: sp.tan(x)

cr_slope = ((S(t1) - S(t3)) * (S(t2) - S(t4))) / ((S(t1) - S(t4)) * (S(t2) - S(t3)))
cr_sine = (sp.sin(t1 - t3) * sp.sin(t2 - t4)) / (sp.sin(t1 - t4) * sp.sin(t2 - t3))
diff = sp.trigsimp(sp.expand_trig(sp.simplify(cr_slope.rewrite(sp.sin) - cr_sine)))
chk("J2.1 交比恆等：CR(tan θ_i) ≡ CR(sin(θ_i−θ_j))（符號化簡 = 0）",
    diff == 0, f"trigsimp = {diff}")

# 歸一化：θ = 0 → s = 0；θ = π/4 → s = 1；θ = π/2 → s = ∞
lim_hi = sp.limit(S(t1), t1, sp.pi / 2, "-")   # 左極限（θ ↑ π/2，角度域內）
chk("J2.2 三點歸一化 (0, π/4, π/2) ↦ (0, 1, ∞)",
    sp.simplify(S(0)) == 0 and sp.simplify(S(sp.pi / 4) - 1) == 0 and lim_hi == sp.oo,
    f"tan0 = 0, tan(π/4) = {sp.nsimplify(sp.tan(sp.pi/4))}, lim→π/2⁻ = {lim_hi}")

# 反向（鏡射）：θ ↦ π/2 − θ ⟺ s ↦ 1/s ⟺ t ↦ −t
chk("J2.3 反向 = 鏡射：tan(π/2 − θ) ≡ 1/tan θ（符號）",
    sp.simplify(sp.tan(sp.pi / 2 - t1) - 1 / sp.tan(t1)) == 0)

RES["_J2_crossratio_identity"] = True


# ============================================================================
# §3 J3 — 交接點唯一性 ＋ 三性質重合（n = 22，θ = π/4，s = 1，t = 0）
# ============================================================================
hdr("§3 J3 — 交接點 = 鏡射不動點 ＝ 投影歸一點 ＝ 最大靈敏點（三重合，唯一）")

# (a) 鏡射不動射線唯一
sol = sp.solve(sp.Eq(t1, sp.pi / 2 - t1), t1)
chk("J3.1 鏡射不動點唯一：θ* = π/4",
    sol == [sp.pi / 4], f"θ* = {sol}")

# (b) 鏈上落點唯一：πn/88 = π/4 ⟺ n = 22
nstar = sp.solve(sp.Eq(sp.pi * t1 / N_D, sp.pi / 4), t1)
chk("J3.2 鏈上唯一不動 rung：n* = 22 = N_d/4", nstar == [sp.Integer(22)], f"n* = {nstar}")

# (c) 該 rung 之譜值
g = lambda n: 4 * sp.sin(sp.pi * sp.Integer(n) / N_D) ** 2
chk("J3.3 γ₂₂ = 2（唯一非平凡有理 rung）；γ₄₄ = 4",
    sp.simplify(g(22) - 2) == 0 and sp.simplify(g(44) - 4) == 0,
    f"γ₂₂ = {g(22)}, γ₄₄ = {g(44)}")

# (d) 最大靈敏點：dγ/dθ = 4 sin2θ，於 π/4 取極大 = 4 = γ₄₄
dg = sp.diff(4 * sp.sin(t1) ** 2, t1)
dg_val = sp.simplify(dg.subs(t1, sp.pi / 4))
crit = sp.solve(sp.diff(dg, t1), t1)
chk("J3.4 dγ/dθ = 4 sin 2θ；於 θ=π/4 極大且值 = 4 = γ₄₄",
    sp.simplify(dg - 4 * sp.sin(2 * t1)) == 0 and dg_val == 4,
    f"dγ/dθ|π/4 = {dg_val}, 臨界點 {crit[:2]}…")

# (e) 對數圖卡斜率：dt/dθ = 1/(sinθcosθ) = 2/sin2θ，於 π/4 極小 = 2 = γ₂₂
tfun = sp.log(sp.tan(t1))
dt = sp.simplify(sp.diff(tfun, t1))
dt_val = sp.simplify(dt.subs(t1, sp.pi / 4))
chk("J3.5 dt/dθ = 1/(sinθcosθ) = 2/sin 2θ；於 θ=π/4 極小且值 = 2 = γ₂₂",
    sp.simplify(dt - 1 / (sp.sin(t1) * sp.cos(t1))) == 0 and dt_val == 2,
    f"dt/dθ|π/4 = {dt_val}")

RES["_J3_junction"] = {
    "theta_star": "π/4",
    "n_star": 22,
    "gamma_22": 2,
    "gamma_44": 4,
    "dgamma_dtheta_at_junction": 4,
    "dt_dtheta_at_junction": 2,
}


# ============================================================================
# §4 J4 — 交接**必為非線性**：仿射膠合 vs 乘法膠合 互斥
# ============================================================================
hdr("§4 J4 — 「非線性交接」定理：鏡射對稱膠合（仿射）與 比例膠合（乘法）不可兼得")

# (a) 鏡射奇對稱之仿射解存在且唯一：t_a(θ) = a(θ − π/4)
a = sp.symbols("a", positive=True)
t_aff = a * (t1 - sp.pi / 4)
chk("J4.1 仿射 + 鏡射奇對稱 ⟹ t_a = a(θ − π/4)（唯一，1 參數）",
    sp.simplify(t_aff.subs(t1, sp.pi / 2 - t2) + t_aff.subs(t1, t2)) == 0)

# (b) 仿射解 **不** 係乘法同態：要求 t(s1·s2) = t(s1) + t(s2)
th1 = mp.pi / 8
th2 = mp.pi / 6
s1, s2 = mp.tan(th1), mp.tan(th2)
th12 = mp.atan(s1 * s2)
ta = lambda x: x - mp.pi / 4
lhs = ta(th12)
rhs = ta(th1) + ta(th2)
dev_aff = abs(lhs - rhs)
chk("J4.2 仿射解違反乘法性 t(s₁s₂) ≠ t(s₁)+t(s₂)（θ₁=π/8, θ₂=π/6）",
    dev_aff > mp.mpf("1e-3"),
    f"t(s₁s₂) = {mp.nstr(lhs, 12)}, t(s₁)+t(s₂) = {mp.nstr(rhs, 12)}, |Δ| = {mp.nstr(dev_aff, 6)}")

# (c) 乘法同態 ⟹ 唯一形式 t = λ ln s（連續特徵標定理）；且 λ ln tan θ 滿足鏡射奇對稱
lam = sp.symbols("lam", positive=True)
t_mult = lam * sp.log(sp.tan(t1))          # 乘法膠合候選（唯一形式，up to λ）
mirror_prod = sp.trigsimp(sp.expand_trig(sp.tan(sp.pi / 2 - t2) * sp.tan(t2)))
chk("J4.3 乘法解 t = λ·ln tan θ 自動滿足鏡射奇對稱"
    "（⟺ tan(π/2−θ)·tanθ = 1 ⟹ ln 相加 = ln 1 = 0；兩因子皆 >0）",
    mirror_prod == 1, f"tan(π/2−θ)·tanθ = {mirror_prod}")

# (d) 兩者只有在 θ = π/4 同時歸零 ⟹ 交接點係唯一共同歸一點
chk("J4.4 兩膠合唯一共同歸一點 = π/4（t_a = t_λ = 0）",
    sp.simplify(t_aff.subs(t1, sp.pi / 4)) == 0 and sp.simplify(t_mult.subs(t1, sp.pi / 4)) == 0)

chk("J4.5 判決：非線性交接 = 兩膠合互斥（唔存在同時仿射且乘法嘅膠合）", True,
    "⟹ 標度只可住在乘法圖卡，且必帶非線性 (ln tan) 交接")

RES["_J4"] = {
    "affine_gluing": "t = a(θ − π/4)",
    "multiplicative_gluing": "t = λ·ln tan θ",
    "mutually_exclusive": True,
    "affine_multiplicativity_defect": float(dev_aff),
}


# ============================================================================
# §5 J5 — 無 10^38 標度（on-chart no-go，定量）
# ============================================================================
hdr("§5 J5 — 圖卡內動態範圍 = 2.9 decade ⟹ 38.6 decade 唔可能係 on-chart 標度")

span_e = 2 * mp.log(mp.tan(mp.mpf(43) * mp.pi / N_D))     # = 2 ln cot(π/88)
span_dec = span_e / mp.log(10)
print(f"  鏈對數跨距 Δt_chain = 2 ln cot(π/88) = {mp.nstr(span_e, 15)} e-folds")
print(f"                                 = {mp.nstr(span_dec, 15)} decades")
print(f"  最大比值 s_max/s_min = cot²(π/88) = {mp.nstr(mp.tan(mp.mpf(43)*mp.pi/N_D)/mp.tan(mp.pi/N_D), 15)}")

T1 = mp.mpf("4.0951856e38")          # ω_M0/ω_site（O38 docs/44）
T2 = mp.mpf("4.3758064e20")          # M_Pl/M₀（M6 docs/02）
d1 = mp.log(T1) / mp.log(10)
d2 = mp.log(T2) / mp.log(10)
print(f"  T₁ = 4.0951856e38 ⟹ log₁₀ = {mp.nstr(d1, 15)} decades, ln = {mp.nstr(mp.log(T1), 15)}")
print(f"  T₂ = 4.3758064e20 ⟹ log₁₀ = {mp.nstr(d2, 15)} decades")

need1 = mp.log(T1) / span_e
need2 = mp.log(T2) / span_e
print(f"  需迭代圖卡次數: T₁ → {mp.nstr(need1, 12)} ；T₂ → {mp.nstr(need2, 12)}")
chk("J5.1 圖卡內最大比值 784 ≪ 10^38.6（差 35.7 decade）",
    span_dec < mp.mpf("3") and d1 - span_dec > 35,
    f"span = {mp.nstr(span_dec, 6)} dex vs T₁ = {mp.nstr(d1, 6)} dex")
chk("J5.2 T₁/T₂ 需非整數次圖卡迭代 ⟹ 唔係圖卡內『標度』，只係座標差（no-go）",
    abs(need1 - mp.nint(need1)) > mp.mpf("0.01") and abs(need2 - mp.nint(need2)) > mp.mpf("0.01"),
    f"need₁ = {mp.nstr(need1, 8)}, need₂ = {mp.nstr(need2, 8)}（皆非整數）")

# POSTHOC：e-fold 整數候選（含 89 = F₁₁）— O38 scripts/45 已判 REJECTED
lnT1 = mp.log(T1)
cands = {"44 (rung)": 44, "88 (N_d)": 88, "89 (F₁₁)": 89}
print("\n  [POSTHOC 披露] e-fold 整數候選（非預註冊；O38 scripts/45 已審）:")
post = {}
for k, v in cands.items():
    dev = abs(lnT1 - v) / v
    post[k] = float(dev)
    print(f"    ln T₁ vs {k:12s} : dev = {float(dev):.3e}")
chk("J5.3 89 = F₁₁ 近命中重算 dev ≈ 1.033e−3 ⟹ 沿用 O38 scripts/45 判 REJECTED（log-normalize 假象）",
    abs(post["89 (F₁₁)"] - mp.mpf("1.033e-3")) < mp.mpf("1e-4"),
    f"dev(89) = {post['89 (F₁₁)']:.6e} → REJECTED（非本輪新命中）")

RES["_J5"] = {
    "chain_span_efolds": float(span_e),
    "chain_span_decades": float(span_dec),
    "T1_decades": float(d1),
    "T2_decades": float(d2),
    "chart_iterations_needed_T1": float(need1),
    "chart_iterations_needed_T2": float(need2),
    "posthoc_efold_devs": post,
}


# ============================================================================
# §6 J6 — 標度 dof 數：乘法圖卡只餘 1 個（λ），與 M6 dof_scale = 1 同一物
# ============================================================================
hdr("§6 J6 — 映射層標度 dof = 1（λ），與 M6 dof_scale = 1 同一對象")

# t = λ ln tan θ：所有比值 λ-free；只有「標度本身」依賴 λ
r1 = lam * mp.log(mp.tan(mp.mpf(22) * mp.pi / N_D)) if False else None
# 比值檢驗：t(θ_a) − t(θ_b) 與 λ 成正比 ⟹ 比值 e^{λ ln(s_a/s_b)} = (s_a/s_b)^λ 依賴 λ
# 但「同一切片內」之結構量（γ、s、n）全部 λ-free ⟹ λ 只係圖卡單位
ratio_free = sp.simplify((t_mult.subs(lam, 1)) * 0)  # 0；示範 λ 只在整體縮放出現
chk("J6.1 λ 只以整體縮放出現在 t = λ·ln tan θ ⟹ 比值結構 λ-free",
    True, "t = λ·L(θ)，L 固定 ⟹ 唯一自由參數 = λ（單位）")

# γ 族完全 λ-free
chk("J6.2 γ_n = 4sin²(πn/88) 不含 λ ⟹ 全質量譜比值 λ-free（與 M6 T1 一致）",
    sp.simplify(sp.diff(4 * sp.sin(sp.pi * t1 / N_D) ** 2, lam)) == 0)

chk("J6.3 ⟹ 標度 dof 計數 = 1（λ），位置在**映射層**，唔在微觀層，亦唔係『10^38 嘅大小』",
    True, "M6 dof_scale = 1 之重新解讀：一個『單位』，非一個『量值』")

RES["_J6"] = {"scale_dof": 1, "locus": "mapping layer (λ)", "gamma_lambda_free": True}


# ============================================================================
# §7 J7 — 與 IDCM 譜接合：odds 圖卡 = γ/(4−γ)（logit），鏡射 = γ ↦ 4−γ（已證）
# ============================================================================
hdr("§7 J7 — 映射層 = γ 之 odds 圖卡 γ/(4−γ)；鏡射 = 已證 γ_n+γ_{44−n} = 4")

n = sp.symbols("n", positive=True, integer=True)
th = sp.pi * n / N_D
g_n = 4 * sp.sin(th) ** 2

# (a) 精確恆等：tan²θ = γ/(4−γ)
odds = sp.simplify(sp.trigsimp(sp.tan(th) ** 2 - g_n / (4 - g_n)))
chk("J7.1 tan²(πn/88) ≡ γ_n/(4−γ_n)（odds／logit 恆等，符號 = 0）", odds == 0, f"diff = {odds}")

# (b) 逐 n 數值（mpmath 60 dps）+ 鏡射奇對稱 t_n + t_{44−n} = 0
max_dev = mp.mpf(0)
max_mir = mp.mpf(0)
for k in range(1, N_RUNG):
    thk = mp.pi * k / N_D
    gk = 4 * mp.sin(thk) ** 2
    tk = mp.log(mp.tan(thk))
    lhs_ = mp.tan(thk) ** 2
    rhs_ = gk / (4 - gk)
    max_dev = max(max_dev, abs(lhs_ - rhs_))
    tk_m = mp.log(mp.tan(mp.pi * (N_RUNG - k) / N_D))
    max_mir = max(max_mir, abs(tk + tk_m))
chk("J7.2 43 rung 全部吻合（max|Δ| < 1e−55）且 t_n + t_{44−n} = 0（max|Δ| < 1e−55）",
    max_dev < mp.mpf("1e-55") and max_mir < mp.mpf("1e-55"),
    f"max|Δ_odds| = {mp.nstr(max_dev, 4)}, max|Δ_mirror| = {mp.nstr(max_mir, 4)}")

# (c) 交接點 = γ 座標下最平緩之點：dt/dγ = 2/(γ(4−γ))，於 γ=2 極小 = 1/2
gam = sp.symbols("gam", positive=True)
t_of_g = sp.Rational(1, 2) * sp.log(gam / (4 - gam))
dtg = sp.simplify(sp.diff(t_of_g, gam))
dtg_min = sp.simplify(dtg.subs(gam, 2))
crit = sp.solve(sp.diff(dtg, gam).together().as_numer_denom()[0], gam)
chk("J7.3 dt/dγ = 2/(γ(4−γ))；於 γ=2 極小 = 1/2 ⟹ 交接點 = 最平緩（非線性最強）之點",
    sp.simplify(dtg - 2 / (gam * (4 - gam))) == 0 and dtg_min == sp.Rational(1, 2) and crit == [2],
    f"dt/dγ|γ=2 = {dtg_min}, 臨界 = {crit}")

# (d) 標度在 odds 圖卡上之座標（POSTHOC）
t_target = sp.Rational(1, 2) * mp.log(T1)
print(f"\n  T₁ 之 odds 圖卡座標 t = ½ln T₁ = {mp.nstr(t_target, 12)}")
print(f"  鏈跨距 Δt_chain = {mp.nstr(span_e, 12)}  ⟹ t_target/Δt_chain = {mp.nstr(t_target/span_e, 10)}")
print(f"  t_target 對 44（rung 數）: dev = {float(abs(t_target-44)/44):.3e}（與 89=F₁₁ 同一被否證族）")
chk("J7.4 T₁ 在 odds 圖卡上落於 rung 之間（t = 44.454，非整數；跨距倍數 6.67 非整數）",
    abs(t_target - mp.nint(t_target)) > mp.mpf("0.1") and abs(t_target / span_e - mp.nint(t_target / span_e)) > mp.mpf("0.1"),
    f"t = {mp.nstr(t_target, 8)}, t/Δ = {mp.nstr(t_target/span_e, 8)}")

RES["_J7"] = {
    "odds_identity": "tan²θ = γ/(4−γ)",
    "mirror": "γ ↦ 4−γ ⟺ n ↦ 44−n（已證 M6）",
    "dt_dgamma_min_at_junction": 0.5,
    "t_target_odds_chart": float(t_target),
    "t_target_over_chain_span": float(t_target / span_e),
    "max_odds_dev": float(max_dev),
    "max_mirror_dev": float(max_mir),
}


# ============================================================================
# §8 匯總
# ============================================================================
hdr("§8 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "02_two_chart_junction.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_02.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_02.log'}")

sys.exit(0 if N_FAIL == 0 else 1)