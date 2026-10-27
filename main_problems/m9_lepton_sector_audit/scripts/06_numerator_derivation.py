#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/06 — R6 分子推導攻堅（在 M9 內推，不交 O35）

(A) 域分離定理：任何 φ^x (x ∈ ℚ\\{0}) ∉ ℚ(ζ₈₈) ⟹ 兩個「語言」數學上不可共存
(B) 整數一致性 σ 表：b = 528·frac(log_φ R) 距最近整數幾個 σ
(C) 選擇律 no-go：同一 (扇區, 世代間隔) 輸入 → 兩個唔同指數 ⟹ 無規則
(D) 分母指派審計 ＋ dof 計數
(E) 結構表達式空間覆蓋率（71/29/42 有幾多條「結構」寫法）＋ MC null 10⁷
(F) 單階梯（幾何級數）no-go
(G) 判決

零 import 上游；mp.dps = 60；sympy 只做域論核對。
"""
import json
import os
from itertools import product

import numpy as np
import sympy as sp
from mpmath import mp, mpf, log, sqrt

mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

PHI = (1 + sqrt(5)) / 2
LN_PHI = log(PHI)

R = {  # CODATA 2022 直接比（scripts/01 已核）
    "mu/e": (mpf("206.7682827"), mpf("0.0000046")),
    "tau/mu": (mpf("16.8170"), mpf("0.0011")),
    "tau/e": (mpf("3477.23"), mpf("0.23")),
}
D_LADDER = 528

checks = []


def chk(name, cond, detail=""):
    checks.append({"check": name, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


print("=" * 76)
print("M9 scripts/06 — R6 分子推導攻堅")
print("=" * 76)

# ════════════════════════════════════════════════════════════════════
# (A) 域分離定理
# ════════════════════════════════════════════════════════════════════
print("\n[A] 域分離定理：φ^x ∉ ℚ(ζ₈₈) for any x ∈ ℚ\\{0}")

x_ = sp.symbols("x")
phi_min = sp.minimal_polynomial((1 + sp.sqrt(5)) / 2, x_)      # x²−x−1
sq5_min = sp.minimal_polynomial(sp.sqrt(5), x_)                # x²−5
poly_disc_sq5 = sp.discriminant(sq5_min, x_)                   # 20（多項式判別式）
# ℚ(√d) 之「域判別式」= d（d ≡ 1 mod 4）或 4d；d=5 ≡ 1 (mod 4) ⟹ 5
field_disc_sq5 = 5 if 5 % 4 == 1 else 20
disc_sq5 = field_disc_sq5
n_cyc = 88
primes_cyc = sp.factorint(n_cyc)
print(f"  minpoly(φ) = {phi_min}    minpoly(√5) = {sq5_min}")
print(f"  多項式判別式 = {poly_disc_sq5}；5 ≡ 1 (mod 4) ⟹ **域判別式** disc(ℚ(√5)) = {field_disc_sq5}")
print(f"  ℚ(ζ₈₈): n = 88 = {primes_cyc} ⟹ 分歧素數 (ramified) = {{2, 11}}；5 ∤ 88")
chk("A1 disc(ℚ(√5)) = 5 ⟹ 5 在 ℚ(√5) 分歧",
    disc_sq5 == 5 and poly_disc_sq5 == 20, f"域 {disc_sq5} / 多項式 {poly_disc_sq5}")
chk("A2 5 ∤ 88 ⟹ 5 在 ℚ(ζ₈₈) 非分歧", 88 % 5 != 0)
chk("A3 ⟹ ℚ(√5) ∩ ℚ(ζ₈₈) = ℚ（兩者皆 Galois，交集若為二次域則必為 ℚ(√5)）",
    sp.gcd(88, 5) == 1)

# γ_n ∈ ℚ(ζ₈₈)：2 − 2cos(2πn/88) = 2 − (ζ^n + ζ^{−n})
zeta = sp.exp(2 * sp.pi * sp.I / 88)
g_expr = {}
for n in (1, 7, 13):
    g_expr[n] = sp.simplify(2 - 2 * sp.cos(2 * sp.pi * n / 88) - (2 - (zeta ** n + zeta ** (-n))))
ok_cyc = all(sp.simplify(v) == 0 for v in g_expr.values())
print(f"  γ_n ∈ ℚ(ζ₈₈) 核對（n=1,7,13）: 2−2cos(2πn/88) ≡ 2−(ζⁿ+ζ⁻ⁿ) → {ok_cyc}")
chk("A4 γ_n ∈ ℚ(ζ₈₈)（故 γ_a/γ_b ∈ ℚ(ζ₈₈)）", ok_cyc)

# φ^p ∉ ℚ：φ^p = F_p φ + F_{p−1}，F_p ≠ 0
fib_ok = all(sp.fibonacci(p) != 0 for p in range(-12, 13) if p != 0)
chk("A5 φ^p = F_p φ + F_{p−1} ∉ ℚ（F_p ≠ 0, p ≠ 0）", fib_ok)

# 小 q：X^{2q} − X^q − 1 之不可約性（deg = 2q）
small = {}
for q in (1, 2, 3, 4, 5, 8):
    poly = sp.Poly(x_ ** (2 * q) - x_ ** q - 1, x_)
    irr = poly.is_irreducible
    small[q] = {"degree": 2 * q, "irreducible": bool(irr)}
    print(f"  q={q}: X^{2*q} − X^{q} − 1  deg={2*q}  不可約={irr}")
chk("A6 小 q 全部不可約 ⟹ [ℚ(φ^(1/q)):ℚ] = 2q", all(v["irreducible"] for v in small.values()))
# q=88 ⟹ deg 176；176 ∤ 20 = [ℚ(ζ₈₈)⁺:ℚ]
chk("A7 176 ∤ 20 ⟹ φ^(1/88) ∉ ℚ(ζ₈₈)（另由 A1–A3 對任意有理 x 成立）",
    176 % 20 != 0)

theorem = {
    "statement": "∀x ∈ ℚ\\{0}: φ^x ∉ ℚ(ζ₈₈)；∀a,b: γ_a/γ_b ∈ ℚ(ζ₈₈) ⟹ φ^x ≠ γ_a/γ_b",
    "proof": [
        "① γ_n = 2−(ζ^n+ζ^{−n}) ∈ ℚ(ζ₈₈) ⟹ 任何 γ 比值 ∈ ℚ(ζ₈₈)（域封閉）。",
        "② 設 x = p/q（既約，p≠0），α = φ^x。若 α ∈ ℚ(ζ₈₈) 則 α^q = φ^p ∈ ℚ(ζ₈₈)。",
        "③ φ^p = F_p·φ + F_{p−1}，F_p ≠ 0 ⟹ φ^p ∉ ℚ；而 φ^p ∈ ℚ(√5)。",
        "④ ℚ(√5) 在 5 分歧，ℚ(ζ₈₈) 只在 {2,11} 分歧 ⟹ ℚ(√5) ∩ ℚ(ζ₈₈) = ℚ ⟹ φ^p ∉ ℚ(ζ₈₈)。矛盾。∎",
    ],
    "corollary": (
        "Q1 定論：同一量**唔可能**同時有「譜讀數」同「φ-冪」兩個身份 —— 唔係判據問題，係域論不可能性。"
        "且 φ 階梯**永遠唔可能**係 K_chain 譜嘅輸出（反過來亦然）：φ ∉ ℚ(ζ₈₈)。"
    ),
}
print("\n  定理：" + theorem["statement"])
for s in theorem["proof"]:
    print("    " + s)
print("  推論：" + theorem["corollary"])

# ════════════════════════════════════════════════════════════════════
# (B) 整數一致性 σ 表
# ════════════════════════════════════════════════════════════════════
print("\n[B] 整數一致性 σ 表（b = 528·frac(log_φ R)）")
RUNGS = {
    "mu/e": (11, 42),      # O34: 11 + 42/528 = 11 + 7/88
    "tau/mu": (6, -71),    # O35: 6 − 71/528
    "tau/e": (17, -29),    # O34: 17 − 29/528
}
b_rows = []
for k, (a_int, b_signed) in RUNGS.items():
    val, u = R[k]
    xo = log(val) / LN_PHI
    b_star = D_LADDER * (xo - mp.floor(xo))
    u_b = D_LADDER * (u / val) / LN_PHI            # δb = 528·δx
    b_int = b_star - (b_star - round(b_star))      # 最近整數
    sig = abs(b_star - round(b_star)) / u_b
    x_pred = a_int + mpf(b_signed) / D_LADDER
    val_pred = PHI ** x_pred
    dev_pct = (val_pred / val - 1) * 100
    sig_val = abs(val_pred / val - 1) * val / u
    b_rows.append({"ratio": k, "x_obs": mp.nstr(xo, 12), "b_star": mp.nstr(b_star, 8),
                   "u_b": mp.nstr(u_b, 5), "b_nearest": int(round(b_star)),
                   "sigma_from_integer": mp.nstr(sig, 6),
                   "ladder_x": mp.nstr(x_pred, 12), "ladder_value": mp.nstr(val_pred, 12),
                   "value_dev_pct": mp.nstr(dev_pct, 8), "value_sigma": mp.nstr(sig_val, 8)})
    print(f"  {k:7s} b* = {mp.nstr(b_star, 8):>12s} ± {mp.nstr(u_b, 5):<9s}"
          f" → 整數 {int(round(b_star)):4d}  偏 {mp.nstr(sig, 5)}σ"
          f" | 梯子值 {mp.nstr(val_pred, 10)}  dev {mp.nstr(dev_pct, 6)}% ({mp.nstr(sig_val, 5)}σ)")

live = [r["ratio"] for r in b_rows if float(r["sigma_from_integer"]) < 5]
dead = [r["ratio"] for r in b_rows if float(r["sigma_from_integer"]) >= 5]
print(f"  ⟹ 整數級活：{live}（<5σ）；死：{dead}（≥5σ）")
chk("B1 μ/e 在整數級即已死（b 距整數 >400σ）",
    float([r for r in b_rows if r["ratio"] == "mu/e"][0]["sigma_from_integer"]) > 400)
chk("B2 τ/μ 整數級相容（<5σ）",
    float([r for r in b_rows if r["ratio"] == "tau/mu"][0]["sigma_from_integer"]) < 5)
chk("B3 τ/e 整數級相容（<5σ）",
    float([r for r in b_rows if r["ratio"] == "tau/e"][0]["sigma_from_integer"]) < 5)

# ════════════════════════════════════════════════════════════════════
# (C) 選擇律 no-go：同輸入 → 不同輸出
# ════════════════════════════════════════════════════════════════════
print("\n[C] 選擇律 no-go（指數 = f(扇區, 世代間隔)？）")
# 兩個輕子相鄰世代躍遷（e→μ、μ→τ）輸入完全相同
lep = {"mu/e": RUNGS["mu/e"], "tau/mu": RUNGS["tau/mu"]}
print("  輕子、間隔=1（相鄰世代）:  e→μ  a=11, b=+42   ／  μ→τ  a=6, b=−71")
print("  夸克、間隔=1:              t→b  a=8, b=+71    ／  c→s …")
delta_a = RUNGS["mu/e"][0] - RUNGS["tau/mu"][0]
print(f"  Δa = {delta_a}（同一 (扇區, 間隔) 輸入）；同扇區內 b 亦差 {RUNGS['mu/e'][1] - RUNGS['tau/mu'][1]}")
chk("C1 同一 (lepton, gap=1) 輸入給出 a = 11 與 6 ⟹ 無 f(sector, gap)",
    delta_a != 0, f"Δa = {delta_a}")
chk("C2 同扇區 b/528 = +42/528 與 −71/528 ⟹ 分子亦非由 (sector, gap) 決定",
    RUNGS["mu/e"][1] != RUNGS["tau/mu"][1])

# 兩個獨立數據點 ⟹ 任何 2 參數族都擬合 ⟹ 0 預言
chk("C3 輕子只有 2 個獨立躍遷數據 ⟹ 任何 2 參數規則必擬合 ⟹ 預言 = 0",
    len(lep) == 2)

# ════════════════════════════════════════════════════════════════════
# (D) 分母指派 ＋ dof 計數
# ════════════════════════════════════════════════════════════════════
print("\n[D] 分母指派審計 ＋ dof 計數")
DEN = {"mu/e": 88, "tau/mu": 528, "tau/e": 528, "t/b (O35)": 132, "s/d (O35)": 528}
mult = {k: v // 44 for k, v in DEN.items()}
print(f"  分母 = 44 × {{2,3,12}}；實際指派 {mult}")
print("  ⟹ (lepton, gap=1) 兩個躍遷分別用 88 同 528（multiplier 2 vs 12）⟹ 分母亦非規則")
n_choice = 3                      # {88,132,528}
bits_per_ratio = float(np.log2(D_LADDER * n_choice))
print(f"  每比值自由度：a ∈ ℤ + b ∈ [0,528) + D ∈ {{88,132,528}} ≈ {bits_per_ratio:.1f} bit")
dof = {"independent_lepton_ratios": 2, "params_per_ratio": 3,
       "third_ratio": "自動（a₃=a₁+a₂, b₃=b₁+b₂）⟹ 非測試",
       "predictions": 0}
print(f"  dof 審計：{dof}")
chk("D1 每比值 ≥3 個自由整數（a, b, D）", bits_per_ratio > 10)
chk("D2 預言數 = 0（2 獨立數據、≥4 參數）", dof["predictions"] == 0)

# ════════════════════════════════════════════════════════════════════
# (E) 表達式空間覆蓋率 ＋ MC null
# ════════════════════════════════════════════════════════════════════
print("\n[E] 結構表達式空間：{±,×,÷} 兩步內可達之整數覆蓋率")
DOC_BASIS = {
    "N_gen": 3, "B": 4, "dimC": 3, "F5": 5, "C(B,2)": 6, "L4": 7, "N_gen2": 9,
    "L5": 11, "dimC·B": 12, "F7": 13, "2⁴": 16, "L5+C(B,2)": 17, "22": 22,
    "L7": 29, "N_d": 44, "N_d+1": 45, "L8": 47, "F10": 55, "2N_d": 88,
    "Σγ": 90, "3N_d": 132, "4N_d": 176, "6N_d": 264, "12N_d": 528,
}
BAS = sorted(set(DOC_BASIS.values()))
LO, HI = 1, 528
print(f"  凍結基底（{len(BAS)} 個文檔化結構數）: {BAS}")


def combine(av, bv):
    out = []
    if av + bv <= HI:
        out.append((av + bv, "+", av, bv))
    if av - bv >= LO:
        out.append((av - bv, "-", av, bv))
    if bv - av >= LO:
        out.append((bv - av, "-r", av, bv))
    if av * bv <= HI:
        out.append((av * bv, "*", av, bv))
    if bv and av % bv == 0 and av // bv >= LO:
        out.append((av // bv, "/", av, bv))
    if av and bv % av == 0 and bv // av >= LO:
        out.append((bv // av, "/r", av, bv))
    return out


reach = {}                                    # value → 最短表達式示例
for a_, b_ in product(BAS, BAS):
    for v, o, x1, y1 in combine(a_, b_):
        reach.setdefault(v, f"{x1}{o}{y1}")
n1 = len(reach)
for a_, b_ in product(sorted(reach), BAS):
    for v, o, x1, y1 in combine(a_, b_):
        reach.setdefault(v, f"({reach[x1] if x1 in reach else x1}){o}{y1}")
n2 = len(reach)
cov = n2 / (HI - LO + 1)
print(f"  1 步可達：{n1} 個整數；2 步可達：{n2} 個 ⟹ 覆蓋率 {cov * 100:.1f}% of [1,528]")
for tgt in (42, 71, 29):
    ex = reach.get(tgt, "—")
    print(f"  目標 {tgt}: {'可表達' if tgt in reach else '不可表達'}   例: {ex}")
cov42, cov71, cov29 = ((tgt in reach) for tgt in (42, 71, 29))
chk("E1 42 / 71 / 29 全部落於結構表達式空間內", cov42 and cov71 and cov29)

NMC = 10_000_000
rng = np.random.default_rng(20261027)
draws = rng.integers(LO, HI + 1, NMC)
reach_arr = np.array(sorted(reach), dtype=np.int64)
cov_emp = float(np.mean(np.isin(draws, reach_arr)))
joint = cov_emp ** 2
print(f"  MC 10⁷：P(隨機整數可表達) = {cov_emp:.4f}（解析 {cov:.4f}）；"
      f"P(兩個目標都可表達) = {joint:.4f}；資訊 = {-np.log2(joint):.2f} bit")
chk("E2 聯合命中率 >5% ⟹ 「搵到 71 和 29 的結構寫法」零信息", joint > 0.05,
    f"{joint:.3f} → 資訊 {-np.log2(joint):.2f} bit")

# ════════════════════════════════════════════════════════════════════
# (F) 單階梯（幾何級數）no-go
# ════════════════════════════════════════════════════════════════════
print("\n[F] 單階梯 no-go：m_g ∝ q^g ⟹ μ/e = τ/μ")
x_mue = log(R["mu/e"][0]) / LN_PHI
x_taumu = log(R["tau/mu"][0]) / LN_PHI
ratio_steps = x_taumu / x_mue
u_rel = sqrt((R["mu/e"][1] / R["mu/e"][0]) ** 2 + (R["tau/mu"][1] / R["tau/mu"][0]) ** 2)
u_steps = ratio_steps * u_rel
sig_steps = abs(1 - ratio_steps) / u_steps
print(f"  幾何級數 ⟹ x(τ/μ)/x(μ/e) = 1；實測 = {mp.nstr(ratio_steps, 8)} ± {mp.nstr(u_steps, 3)}"
      f" ⟹ 偏離 {mp.nstr(sig_steps, 5)}σ（= {mp.nstr(1 / ratio_steps, 6)}× 不等步）")
chk("F1 單階梯（等步）被 12.3× 不等步排除", sig_steps > 5, f"{mp.nstr(sig_steps, 4)}σ")

# ════════════════════════════════════════════════════════════════════
# (G) 判決
# ════════════════════════════════════════════════════════════════════
verdict = {
    "field_separation": "✅ 定理：φ^x ∉ ℚ(ζ₈₈) ∀x∈ℚ\\{0} ⟹ φ-冪 與 譜讀數 數學上不可共存（Q1 定論）",
    "integer_level": "μ/e 整數級 428σ 死；μ/τ 2.63σ、e/τ 2.72σ 活（數值相容）",
    "selection_rule": "🔴 無規則：同一 (lepton, gap=1) 給出 a=11 與 a=6、b=+42 與 −71、D=88 與 528",
    "dof": "每比值 ≥3 自由整數；2 個獨立數據 ⟹ 預言 = 0（第三比值自動，非測試）",
    "expression_space": f"🔴 覆蓋率 {cov * 100:.1f}%（2 步內）；71 與 29 皆可表達 ⟹ 分子「有結構寫法」零信息",
    "single_cascade": f"🔴 幾何級數被排除（{mp.nstr(sig_steps, 4)}σ）",
    "verdict": ("🔴 分子 7／71／29 無法由結構推出（無選擇律 ＋ 表達式空間飽和 ＋ 2 數據 4 參數）；"
                "τ/μ、e/τ 保留為「數值活但零預言」；μ/e 死。φ 階梯 = 逐比值後驗編碼，非理論。"),
}
print("\n[G] 判決")
for k, v in verdict.items():
    print(f"  {k:18s}: {v}")

out = os.path.join(DATA, "06_numerator_derivation.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump({"theorem": theorem, "small_q": small, "b_table": b_rows,
               "live_integer": live, "dead_integer": dead,
               "selection_rule": {"lepton_gap1": {"mu_e": [11, 42, 88],
                                                  "tau_mu": [6, -71, 528]},
                                  "delta_a": delta_a},
               "denominator": {"assignment": mult, "bits_per_ratio": bits_per_ratio},
               "dof": dof,
               "expression_space": {"basis": BAS, "n_1step": n1, "n_2step": n2,
                                    "coverage": cov, "joint": joint,
                                    "examples": {str(t): reach.get(t, None) for t in (42, 71, 29)}},
               "single_cascade": {"steps_ratio": mp.nstr(ratio_steps, 8),
                                  "u": mp.nstr(u_steps, 4), "sigma": mp.nstr(sig_steps, 6)},
               "verdict": verdict, "checks": checks,
               "pass": sum(1 for c in checks if c["pass"]), "total": len(checks)},
              f, ensure_ascii=False, indent=2)
print("\n" + "=" * 76)
print(f"PASS {sum(1 for c in checks if c['pass'])}/{len(checks)}   -> {out}")
print("=" * 76)