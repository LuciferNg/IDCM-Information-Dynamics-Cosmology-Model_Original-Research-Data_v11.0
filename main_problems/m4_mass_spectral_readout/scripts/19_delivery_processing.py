"""M4 docs/19 配套 — 處理 O140／P6 交付：核實 ℤ[γ₁] 約束、形式框架可辨識性、A 身份再測（2026-10-24）
====================================================================================================
輸入交付：
  O140 `docs/19_m4_handoff_processing.md`（2026-10-24 13:14 register handed_off O140 → M4）：
    ① R2 核實通過；② 三約束：m_k=M₀γ_k（44 級嚴格遞增）、γ₄₄=4 ⟹ m₄₄=4M₀、
       **單模可達比上限 = 4/γ₁ = 784.9646644423913**、多項式可達集 m = M₀p(γ₁), p ∈ ℤ[x], deg ≤ 43；
    ③ 🔴 flavour→mode 映射**不存在**（三期 null）；禁向 O140 索取該表
  P6 `01_problem_record.md`（13:13 register handed_off P6 → M4）：A 譜身份交回 M4；
    P6 鏈 m_ν = A·γ_n²，m₁/₂/₃ = 0.0209452/0.0226848/0.0536815 eV，Σm = 0.0973116 eV，盲測 0.80%/0.37%

本輪（真計算）：
  A. 獨立核實 O140 五項（γ₁、γ₄₄、1/γ₁、4/γ₁、遞推、γ₂=4γ₁−γ₁²）— sympy 精確
  B. 應用單模上限：已知質量比邊界分類（哪些必須多模/多項式）
  C. **可辨識性審計**（關鍵）：
     C1 deg-0 空洞性：ℤ ⊂ ℤ[γ₁] ⟹ 任何「近整數」比值都被 deg-0 多項式「解釋」
     C2 deg≥1 密度：族 {p(γ₁): deg ≤ 4, |c_i| ≤ 10} 對隨機目標嘅 1% 覆蓋率
  D. A 身份（P6 交回）：量級測試 ＋ {n·γ₁^k} 族掃描 ＋ ℤ[γ₁] 整係數可達性
判準（事前）：核實誤差 ≤ 1e-12 ⟹ ✅；C2 覆蓋率 > 5% ⟹ 「無鑑別力」；D 命中須 <0.1% 且密度 <5% 才可升格
"""
import json
import math
import numpy as np
import sympy as sp

OUT = "data/19_delivery_processing.json"
L = 88
g1 = 4 * math.sin(math.pi / L) ** 2

print("=" * 84)
print("M4 docs/19 — 處理 O140／P6 交付（核實 + 可辨識性 + A 再測）")
print("=" * 84)

# ---------- A. 核實 O140 ----------
print("\n[A] 獨立核實 O140 交付數值／代數")
x = sp.symbols("x")
gam = lambda n: 4 * sp.sin(sp.pi * n / L) ** 2
vals = {"γ₁": float(gam(1)), "γ₄₄": float(gam(44)), "1/γ₁": float(1 / gam(1)),
        "4/γ₁": float(4 / gam(1))}
ref = {"γ₁": 0.005095770779492917, "γ₄₄": 4.0, "1/γ₁": 196.24116611059782,
       "4/γ₁": 784.9646644423913}
for k in vals:
    print(f"  {k:<6} 重算 {vals[k]:.16g}  O140 {ref[k]:.16g}  |Δ| {abs(vals[k]-ref[k]):.2e}")

# 遞推 γ_(n+1) = 2cγ_n − γ_(n−1) + 2γ₁, c = cos(π/44)
c = math.cos(math.pi / 44)
gm = lambda n: 4 * math.sin(n * math.pi / L) ** 2
rec_err = max(abs(gm(n + 1) - (2 * c * gm(n) - gm(n - 1) + 2 * g1)) for n in range(2, 44))
print(f"  遞推 γ_(n+1)=2cγ_n−γ_(n−1)+2γ₁：max|err| = {rec_err:.2e}（n=2..43）")
# γ₂ = 4γ₁ − γ₁² 精確（sympy）
cs = sp.symbols("c", real=True)
g1s, g2s = 2 - 2 * cs, 2 - 2 * (2 * cs ** 2 - 1)      # γ_n = 2 − 2T_n(c)
simp = sp.simplify(g2s - (4 * g1s - g1s ** 2))
print(f"  γ₂ − (4γ₁ − γ₁²)（c 代入，sympy）= {simp} ⟹ {'✅ EXACT' if simp == 0 else '🔴 非零'}")
# c = 1 − γ₁/2
print(f"  c − (1 − γ₁/2) = {c - (1 - float(gam(1))/2):.2e}")
CAP = 4 / g1
print(f"  單模可達比上限 CAP = 4/γ₁ = {CAP:.10f}")

# ---------- B. 單模上限 vs 已知質量比 ----------
print("\n[B] 單模上限 CAP = 784.965 對已知質量比（PDG 2024 引用值，MeV）")
M = {"e": 0.51099895, "μ": 105.6583755, "τ": 1776.86, "u": 2.16, "d": 4.70,
     "s": 93.46, "c": 1273.0, "b": 4180.0, "t": 172760.0, "p": 938.27208816}
pairs = [("τ/e", "τ", "e"), ("t/e", "t", "e"), ("t/c", "t", "c"), ("μ/e", "μ", "e"),
         ("b/s", "b", "s"), ("c/s", "c", "s"), ("s/d", "s", "d"), ("d/u", "d", "u"),
         ("τ/μ", "τ", "μ"), ("p/e", "p", "e")]
rowsB = []
for lab, a, b in pairs:
    R = M[a] / M[b]
    over = R > CAP
    rowsB.append({"pair": lab, "ratio": R, "exceeds_cap": over})
    print(f"  {lab:<5} = {R:>12.4f}  {'🔴 > CAP ⟹ 必須多模/多項式' if over else '🟡 單模可達'}")

# ---------- C1. deg-0 空洞性 ----------
print("\n[C1] deg-0 空洞性：ℤ ⊂ ℤ[γ₁] ⟹ 近整數比值被「解釋」")
rowsC1 = []
for lab, a, b in [("τ/e", "τ", "e"), ("p/e", "p", "e"), ("μ/e", "μ", "e"), ("τ/μ", "τ", "μ")]:
    R = M[a] / M[b]
    n = round(R)
    dev = abs(R / n - 1) * 100
    rowsC1.append({"pair": lab, "ratio": R, "nearest_int": n, "dev_pct": dev})
    print(f"  {lab:<5} R = {R:>10.4f}  nearest int {n:>6}  dev {dev:>7.4f}%  {'✅ <1%' if dev < 1 else '—'}")

# ---------- C2. deg≥1 密度 ----------
print("\n[C2] deg ≥ 1 族密度：{p(γ₁): deg ≤ 4, |c_i| ≤ 10}")
Cs = list(range(-10, 11))
arr = np.array([0.0])
for k in range(5):
    arr = np.add.outer(arr, np.array([cc * g1 ** k for cc in Cs])).ravel()
arr = np.sort(arr[(arr > 0)])
print(f"  族大小（正元素）= {arr.size}")
rng = np.random.default_rng(20261033)
tgt = 10 ** rng.uniform(0, 5, 200000)        # 1 … 1e5
idx = np.clip(np.searchsorted(arr, tgt), 1, arr.size - 1)
dev = np.minimum(np.abs(arr[idx - 1] / tgt - 1), np.abs(arr[idx] / tgt - 1)) * 100
cov = {thr: float((dev < thr).mean()) for thr in (0.1, 1.0)}
for thr, v in cov.items():
    print(f"  隨機目標落 {thr:>4}% 帶嘅機率 = {v*100:6.2f}%")

# ---------- D. A 身份（P6 交回） ----------
print("\n[D] A 譜身份（A = 1.64337 eV，P6 交回）")
M0 = 27.900915e6        # eV（O38 錨，MeV→eV）
A = 1.64337
ratio = A / M0
print(f"  A/M₀ = {ratio:.6e}   （γ₁ = {g1:.6e}，γ₁⁴ = {g1**4:.6e}）")
best = None
for k in range(1, 13):
    n = ratio / g1 ** k
    ni = round(n)
    if ni == 0:
        continue
    dev = abs(ni * g1 ** k / ratio - 1) * 100
    if best is None or dev < best[3]:
        best = (ni, k, ni * g1 ** k * M0, dev)
    print(f"  k={k:>2}: n = ratio/γ₁^k = {n:>12.4f} ⟹ 最近整數 {ni:>4} ⟹ A_pred = {ni*g1**k*M0:>10.4f} eV  dev {dev:>8.3f}%")
print(f"  最佳（無界限）：n={best[0]}, k={best[1]} ⟹ A_pred = {best[2]:.5f} eV, dev {best[3]:.3f}%（n 巨大 ⟹ 無結構意義）")
best_small = None
for k in range(1, 9):
    for n in range(1, 201):
        d_ = abs(n * g1 ** k * M0 / A - 1) * 100
        if best_small is None or d_ < best_small[3]:
            best_small = (n, k, n * g1 ** k * M0, d_)
print(f"  最佳（係數界限 n ≤ 200）：n={best_small[0]}, k={best_small[1]} ⟹ A_pred = {best_small[2]:.5f} eV, dev {best_small[3]:.3f}%")
print(f"    結構備註：n=87 = L₈−1（O85 亦有 Λ = M_Pl/87）⟹ 只可當線索（docs/15 身份冗餘教訓）")
# 密度：同 (n ≤ 200, k ≤ 12) 族中落 1% 嘅組合數
nhit = 0
for k in range(1, 13):
    for n in range(1, 201):
        if abs(n * g1 ** k / ratio - 1) < 0.01:
            nhit += 1
print(f"  族 {{n·γ₁^k: n≤200, k≤12}} = 2400 個值；落 A/M₀ 1% 帶嘅組合數 = {nhit} ⟹ 密度 {nhit/2400*100:.2f}%")
# ℤ[γ₁] 整係數可達性：需要 p(γ₁) = ratio，而最小正元素 ≈ γ₁
print(f"  量級障：p(γ₁) 最小正值 = γ₁ = {g1:.3e} ≫ A/M₀ = {ratio:.3e}（差 {g1/ratio:.0f}×）")
print(f"  ⟹ A 唔可能係 ℤ[γ₁] 低次整係數多項式值（除非高次近抵消；見上 k 掃描最佳 dev {best[3]:.2f}%）")

res = {"meta": {"script": "19_delivery_processing", "date": "2026-10-24", "seed": 20261033},
       "A_verify": {"gamma1": g1, "gamma44": float(gm(44)), "inv_g1": 1 / g1, "cap_4_over_g1": CAP,
                    "recurrence_max_err": rec_err, "gamma2_identity_sympy": str(simp),
                    "c_identity_err": c - (1 - g1 / 2)},
       "B": rowsB, "C1": rowsC1,
       "C2": {"family_size": int(arr.size), "p_within_1pct": cov[1.0], "p_within_0p1pct": cov[0.1]},
       "D": {"A_over_M0": ratio, "best_unbounded": {"n": best[0], "k": best[1], "dev_pct": best[3]},
             "best_n_le_200": {"n": best_small[0], "k": best_small[1], "A_pred_eV": best_small[2], "dev_pct": best_small[3]},
             "family_hits": nhit, "family_size": 2400}}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")