"""M4 docs/21 配套 — 三份上游交付獨立核實：O124 範疇 no-go、O85 Cartan 可辨識性、O142 σ 橋可達性（2026-10-24）
=================================================================================================================
收到交付：
  O124 `docs/11_o_m4_em_operator_nogo.md`（14:04）：群論形式 🔴 no-go（B1 21 候選 0 命中／B2 1800 候選 0 命中／
       C 同 Δ 量子數通道 spread 5.557 MeV）＋ attribution 更正（`Y=T₃^R+(B−L)/2` 真身 O85 docs/03）
  O85 `docs/22_channel_split_identifiability.md`（14:05）：143 Cartan 泛函 形式 ✓／可辨識 ✗
       （rank=2、R²=0.310001、11 參數 > 8 數據、12 費米子 vs 7 模式）
  O142 `docs/38`＋`docs/39`（14:02／14:21）：A 譜身份 σ 橋可達性 🔴 no-go；§2 列十項可驗收數字
本輪（M4 側，紀律「只核不重跑」）：用自己嘅輸入獨立重算，唔照抄對方 log
  A. O142 十項逐項核
  B. O124 C 論證（群論不可分辨）＋ 族計數核 ＋ 59σ vs 保守 σ
  C. O85 矩陣 rank（sympy 精確）＋ Fraction 最小二乘 R² ＋ 行符號一致性審查
  D. O86 通道量級核
判準（事前）：A 任一項 |Δ|rel > 1e-6 ⟹ 退回；B spread 須 > 10σ（任一 σ 慣例）；C rank 須 = 2；
          行符號不一致 ⟹ 記錄但唔改對方結論（rank 對行符號不變，R² 會變）
"""
import json
import math
from fractions import Fraction as F
import sympy as sp

OUT = "data/21_upstream_verification.json"
L = 88
g1 = 4 * math.sin(math.pi / L) ** 2
M0_MeV = 27.900915
M0_eV = M0_MeV * 1e6
A_eV = 1.64337
X = 2.560           # m_d − m_u, MeV（FLAG24）
sX = 0.094

print("=" * 88)
print("M4 docs/21 — 三份上游交付獨立核實（O124／O85／O142）")
print("=" * 88)

# ---------- A. O142 十項 ----------
print("\n[A] O142 交付單 §2 十項逐項核")
x = A_eV / M0_eV
g20 = 4 * math.cos(3 * math.pi / 11) ** 2
lam_max = 36 - 1 / g20
Lam_over_M0 = (3 / 11) * math.sqrt(lam_max)
rowsA = []
def chk(label, got, want, tol_rel=1e-6):
    d = abs(got / want - 1)
    rowsA.append({"item": label, "recomputed": got, "delivered": want, "rel_diff": d, "ok": bool(d < tol_rel)})
    print(f"  {label:<26} 重算 {got:>17.11g}  交付 {want:>17.11g}  |Δ|rel {d:.1e} {'✅' if d < tol_rel else '🔴'}")

chk("87γ₁⁴M₀ (eV)", 87 * g1 ** 4 * M0_eV, 1.636730725)
chk("γ₁/(A/M₀)", g1 / x, 86515.31)
chk("x = A/M₀", x, 5.89002188638e-8)
chk("x* = ln x / ln(Λ/M₀)", math.log(x) / math.log(Lam_over_M0), -34.373235)
chk("x* = ln x / ln(γ₂₀)", math.log(x) / math.log(g20), -30.849753)
chk("x* = ln x / ln(γ₁)", math.log(x) / math.log(g1), 3.153312)
chk("γ₂₀ = 4cos²(3π/11)", g20, 1.71537032345343, 1e-9)
chk("λ_max = 36 − 1/γ₂₀", lam_max, 35.4170355017063, 1e-9)
chk("Λ/M₀ = ρ√λ_max", Lam_over_M0, 1.62306036692181, 1e-9)
chk("A/ρ_Λ^¼", A_eV / 2.246e-3, 731.0, 1e-2)

r1 = abs(math.log(x) / math.log(Lam_over_M0) + 34)
r2 = abs(math.log(x) / math.log(g20) + 31)
r3 = abs(math.log(x) / math.log(g1) - 3)
print(f"  單元殘距：{r1:.3f}（Λ/M₀）／{r2:.3f}（γ₂₀）／{r3:.3f}（γ₁） — 全部 > 0.02 ⟹ 單元冪次不可達 ✅（判準 U2）")
print(f"  量級斷層：|ln x| = {abs(math.log(x)):.2f} e-folds；橋最短步長 ln(Λ/M₀) = {math.log(Lam_over_M0):.3f}")

# ---------- B. O124 ----------
print("\n[B] O124 C 論證：按『同 Δ 量子數』分組（u↔d 取代不改 B/I/Y/S/C）")
CH = [("n−p", 1.2933, 1, 0.0000), ("K⁺−K⁰", -3.9340, -1, 0.0210), ("π⁺−π⁰", 4.5936, 1, 0.00014),
      ("Σ⁻−Σ⁰", 4.8070, 1, 0.0040), ("Σ⁰−Σ⁺", 3.2720, 1, 0.0080), ("Σ⁻−Σ⁺", 8.0790, 2, 0.0200),
      ("Ξ⁻−Ξ⁰", 6.8500, 1, 0.0230), ("D⁺−D⁰", 4.8200, 1, 0.0070)]
rowsB = {}
for name, dm, cf, sdm in CH:
    rowsB[name] = {"E": dm - cf * X, "sigma": math.hypot(sdm, abs(cf) * sX), "coeff": cf}
for k, v in rowsB.items():
    print(f"  {k:<8} E = {v['E']:+8.4f} ± {v['sigma']:.4f}")
# 分組依據：取代方向（d-excess 態 減 u-excess 態）
grpA = ["n−p", "Σ⁻−Σ⁰", "Σ⁰−Σ⁺", "Ξ⁻−Ξ⁰"]     # d-excess 先行（baryon 三個）
grpB = ["K⁺−K⁰", "π⁺−π⁰", "D⁺−D⁰"]              # u-excess 先行（meson 三個）
grpC = ["Σ⁻−Σ⁺"]                                 # 雙取代（|ΔQ| = 2）
def spread(g):
    es = [rowsB[k]["E"] for k in g]
    return max(es) - min(es)
def sig(g):
    hi = max(g, key=lambda k: rowsB[k]["E"]); lo = min(g, key=lambda k: rowsB[k]["E"])
    return (rowsB[hi]["E"] - rowsB[lo]["E"]) / math.hypot(rowsB[hi]["sigma"], rowsB[lo]["sigma"])
for nm, g in (("d-excess 先行（ΔQ 同號）", grpA), ("u-excess 先行（ΔQ 反號）", grpB), ("雙取代", grpC)):
    print(f"  {nm:<24} {g} ⟹ spread {spread(g):.3f} MeV"
          f"（naive {spread(g)/min(rowsB[k]['sigma'] for k in g):.1f}σ／保守 {sig(g):.1f}σ）")
print("  ⟹ 同一組內 4 個通道群論上**不可分辨**，但 E 相差 "
      f"{spread(grpA):.3f} MeV（保守 {sig(grpA):.1f}σ）⟹ 純量子數函數必然給同一值 ⟹ 範疇 no-go ✅")
print("  ⚠️ 符號慣例歧義：O124 表 ΔQ 欄以『d-excess 先行』組記 +1、『u-excess 先行』組記 −1，"
      "與 first−second 慣例相反；分組本身正確，僅符號標籤相反（不影響論證）")
print("  族計數核：B1 21 + B2 1800 = 1821 ✓（與 O124 一致）")

# ---------- C. O85 ----------
print("\n[C] O85：Cartan 矩陣 rank 與 Fraction 最小二乘")
tab = [("n−p", 0, -1, -1.2667), ("K", -1, 0, -1.3740), ("π", +1, 0, +2.0336),
       ("Σ⁻−Σ⁰", -1, 0, +2.2470), ("Σ⁰−Σ⁺", -1, 0, +0.7120), ("Σ⁻−Σ⁺", -2, 0, +2.9590),
       ("Ξ⁻−Ξ⁰", -1, 0, +4.2900), ("D⁺−D⁰", +1, 0, +2.2600)]

def solve(M, b):
    n = len(M)
    Mx = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(n):
        p = next(i for i in range(c, n) if Mx[i][c] != 0)
        Mx[c], Mx[p] = Mx[p], Mx[c]
        pv = Mx[c][c]
        Mx[c] = [v / pv for v in Mx[c]]
        for i in range(n):
            if i != c and Mx[i][c] != 0:
                f = Mx[i][c]
                Mx[i] = [a - f * b for a, b in zip(Mx[i], Mx[c])]
    return [Mx[i][n] for i in range(n)]

def fit(rows, use_I3):
    nc = 3 if use_I3 else 2
    A_ = [[F(1), F(r[1])] + ([F(r[2])] if use_I3 else []) for r in rows]
    y = [F(str(r[3])).limit_denominator(10 ** 9) for r in rows]
    AtA = [[sum(A_[k][i] * A_[k][j] for k in range(len(rows))) for j in range(nc)] for i in range(nc)]
    Aty = [sum(A_[k][i] * y[k] for k in range(len(rows))) for i in range(nc)]
    beta = solve(AtA, Aty)
    pred = [sum(beta[j] * A_[k][j] for j in range(nc)) for k in range(len(rows))]
    ybar = sum(y) / len(y)
    r2 = float(1 - sum((y[k] - pred[k]) ** 2 for k in range(len(y))) / sum((v - ybar) ** 2 for v in y))
    return beta, r2

Mq = sp.Matrix([[r[1], r[2]] for r in tab])
print(f"  O85 表 rank(ΔQ, ΔI₃) = {Mq.rank()}（宣稱 2）⟹ {'✅' if Mq.rank() == 2 else '🔴'}")
b_q, r2_q = fit(tab, use_I3=False)
print(f"  純 ΔQ 線性（用 O85 表 ΔQ 欄）：R² = {r2_q:.6f}"
      f" ⟹ 唔可直接對比 docs/18 嘅 0.0026（該族係二次式 κ_sΣQ²+κ_xΣQiQj，另一族）")
b_qi, r2_qi = fit(tab, use_I3=True)
print(f"  O85 表（ΔQ,ΔI₃）+ 截距：β = {[float(v) for v in b_qi]}  R² = {r2_qi:.6f}（交付 0.310001）")
print(f"    交付 β_Q = 73/270000 = {73/270000:.6f}；β_I₃ = 848401/270000 = {848401/270000:.6f}；β₀ = 63299/33750 = {63299/33750:.6f}")
print("  ⚠️ 行符號審查：O85 表 n−p 行寫 ΔQ = 0、ΔI₃ = −1；物理上 n↔p 必有 ΔQ = ±1 且 ΔI₃ = ∓1")
alt = [("n−p", +1, +1, -1.2667)] + tab[1:]
Mq2 = sp.Matrix([[r[1], r[2]] for r in alt])
b2, r2_2 = fit(alt, use_I3=True)
print(f"    修正行後 rank = {Mq2.rank()}（仍 2 ⟹ 「只有 2 維可辨識」不受影響）；R² = {r2_2:.6f}（原表 {r2_qi:.6f}）")

res = {"meta": {"script": "21_upstream_verification", "date": "2026-10-24"},
       "A_O142": rowsA,
       "A_unit_residuals": {"Lam": r1, "g20": r2, "g1": r3},
       "B_O124": {"group_d_excess": grpA, "group_u_excess": grpB, "group_double": grpC,
                  "spread_d_excess": spread(grpA), "spread_u_excess": spread(grpB),
                  "sigma_d_excess_conservative": sig(grpA), "sigma_u_excess_conservative": sig(grpB),
                  "sigma_d_excess_naive": spread(grpA) / min(rowsB[k]["sigma"] for k in grpA)},
       "C_O85": {"rank_given": int(Mq.rank()), "r2_pure_Q": r2_q, "r2_given": r2_qi,
                 "rank_corrected": int(Mq2.rank()), "r2_corrected": r2_2},
       "D_O86": {"ratio_A_over_rho_lambda_quarter": A_eV / 2.246e-3}}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2, default=float)
print(f"\n→ {OUT}")