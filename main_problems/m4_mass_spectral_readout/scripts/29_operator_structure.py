"""M4 docs/29 配套 — ⟨O⟩ 之結構審計：⟨O⟩≡1 之 χ² 否證、n_s 線性檢定、零參數候選與盲測設計
=========================================================================================
輸入（`docs/28`，由 BMW PRL 111, 252001 TABLE I 反推）：
    ⟨O⟩_N = 0.891 ± 0.107（n_s = 0）
    ⟨O⟩_Σ = 1.498 ± 0.262（n_s = 1）
    ⟨O⟩_Ξ = 2.293 ± 0.351（n_s = 2）
判準（事前寫死）：
    C1 `⟨O⟩ ≡ 1` 之 χ²（3 dof）→ 若 p < 1% 則否證。
    C2 線性 `a + b·n_s` 之加權擬合：報 (a, b) 及協方差；檢定 a = 1。
    C3 候選模型**必須**分兩類標記：**零參數（事前）** vs **後驗選取（睇完 3 錨才選）**。
        本輪全部候選皆屬**後驗**（因 3 錨已在手），故一律不得升格；只可作盲測設計之輸入。
    C4 盲測：候選模型對 (K, D) 之預言寫死並儲存，待獨立 lattice QED 值落刀。
    C5 唔准用 `docs/18` 之 `E` 值反推 ⟨O⟩ 再「驗證」候選（循環）。
輸出：`data/29_operator_structure.json`
"""
import json

import numpy as np

OBS = [(0, 0.891, 0.107, "N"), (1, 1.498, 0.262, "Σ"), (2, 2.293, 0.351, "Ξ")]
X, SX = 2.56, 0.094


def chi2(pred, obs=OBS):
    return sum(((o - pred(n)) / s) ** 2 for n, o, s, _ in obs)


print("=" * 100)
print("⟨O⟩ 結構審計：⟨O⟩ ≡ 1 否證 ＋ n_s 線性檢定 ＋ 零參數候選（全部後驗，付盲測）")
print("=" * 100)
print("輸入：", "；".join(f"{lbl}(n_s={n}) = {o:.3f} ± {s:.3f}" for n, o, s, lbl in OBS))

from math import exp, erfc, gamma as _g, sqrt  # noqa: E402


def chi2_sf(x, k):
    """上尾機率（regularized upper incomplete gamma），k 為自由度"""
    if k % 2 == 0:
        s, t, i = 0.0, 1.0, 0
        while i < k // 2:
            s += t
            i += 1
            t *= x / 2 / i
        return exp(-x / 2) * s
    return erfc(sqrt(x / 2)) + exp(-x / 2) * sum(
        (x / 2) ** (j + 0.5) / _g(j + 1.5) for j in range(0))


# C1
c1 = chi2(lambda n: 1.0)
p1 = chi2_sf(c1, 3)
print(f"\n[C1] ⟨O⟩ ≡ 1：χ² = {c1:.2f} / 3 dof", end="")
print(f" ⟹ p = {p1:.2e} ⟹ {'🔴 否證' if p1 < 0.01 else '✅ 未否證'}")

# C2 加權線性擬合
A = np.array([[1.0, n] for n, _, _, _ in OBS])
y = np.array([o for _, o, _, _ in OBS])
w = np.array([1 / s ** 2 for _, _, s, _ in OBS])
W = np.diag(w)
cov = np.linalg.inv(A.T @ W @ A)
beta = cov @ (A.T @ W @ y)
a, b = beta
sa, sb = np.sqrt(np.diag(cov))
res = y - A @ beta
c2 = float(res @ W @ res)
print(f"\n[C2] 加權線性 a + b·n_s：a = {a:.3f} ± {sa:.3f}（距 1 有 {abs(a-1)/sa:.2f}σ）；b = {b:.3f} ± {sb:.3f}")
print(f"     χ² = {c2:.2f} / 1 dof ⟹ p = {chi2_sf(c2,1):.3f} ⟹ "
      f"{'✅ 線性相容' if chi2_sf(c2,1) > 0.05 else '🟡 線性有張力'}")

# C3 候選（全部後驗）
MODELS = {
    "M1: 1 + (2/3)·n_s = 1 + (B−2)n_s/N_q": lambda n: 1 + (2 / 3) * n,
    "M2: 3/(3 − n_s)": lambda n: 3 / (3 - n),
    "M3: 1 + 0.5·n_s": lambda n: 1 + 0.5 * n,
    "M4: 1 + 0.75·n_s": lambda n: 1 + 0.75 * n,
    "M5: 自由線性（2 參數；參考）": lambda n: a + b * n,
}
print("\n[C3] 候選模型 χ²（**全部屬後驗**；M5 有 2 參數，其餘近似 0 參數）")
best = None
for name, f in MODELS.items():
    c = chi2(f)
    dof = 1 if name.startswith("M5") else 3
    p = chi2_sf(c, dof)
    print(f"     {name:<44} χ² = {c:>6.2f} / {dof} dof   p = {p:.3f}")
    if best is None or c < best[1]:
        best = (name, c)
print(f"     ⟹ 最佳：{best[0]}（χ² = {best[1]:.2f}）")

# C4 盲測預言（候選 M1；取向：em_M4 = Δm − coeff·X·⟨O⟩）
print("\n[C4] 盲測預言（用 M1：⟨O⟩ = 1 + (2/3)n_s；寫死待外部核）")
# coeff 取向鐵律（docs/18）：u→d 方向（末態較重）⟹ coeff = +1；d→u（末態較輕）⟹ coeff = −1。
#   K⁺(u s̄) → K⁰(d s̄)：u→d ⟹ coeff = −1（故 docs/18 之 E = −3.934 + 2.56 = −1.374 ✓）
#   D⁺(c d̄) → D⁰(c ū)：夸克質量 d̄ > ū ⟹ coeff = +1（docs/18 之 E = +4.82 − 2.56 = +2.26 ✓）
#   其餘通道皆 d→u 或含 d→u 兩次 ⟹ coeff = +1 / +2
CAND = [
    ("K⁺−K⁰", -1, -3.9340, 1, 0.0963, "n_s = 1（s̄ 旁觀）；coeff = −1"),
    ("D⁺−D⁰", +1, +4.8200, 0, 0.1172, "n_s = 0（c 旁觀）；coeff = +1"),
    ("Σ⁰−Σ⁺", +1, +3.2720, 1, 0.1196, "n_s = 1；coeff = +1"),
    ("Σ⁻−Σ⁰", +1, +4.8070, 1, 0.1010, "n_s = 1；coeff = +1"),
    ("Σ⁻−Σ⁺", +2, +8.0790, 1, 0.2028, "n_s = 1；coeff = +2"),
]
preds = []
for name, nf, dm, ns, sE, note in CAND:
    o = 1 + (2 / 3) * ns
    em = dm - nf * X * o
    preds.append({"channel": name, "O_pred": o, "em_pred": em, "n_s": ns, "sigma_E": sE, "note": note})
    print(f"     {name:<9} ⟨O⟩_pred = {o:.3f} ⟹ em_pred = {em:+.4f} MeV   （{note}）")

print("\n[C4b] 判別力：要在 n_s = 1 之新錨上分辨 b = 2/3（0.667）vs 0.5 / 0.75")
for btest in (0.5, 0.75):
    need = abs(2 / 3 - btest) * 1  # Δ⟨O⟩ at n_s = 1
    print(f"     b = {btest} ⟹ ⟨O⟩(n_s=1) 差 {need:.3f} ⟹ 要 3σ 分辨需 σ_⟨O⟩ ≤ {need/3:.3f}"
          f"（現 Σ 之 σ = 0.262 ⟹ {'不足' if 0.262 > need/3 else '足夠'}）")

json.dump({"chi2_O_eq_1": c1, "p_O_eq_1": p1,
           "linear_fit": {"a": a, "sa": sa, "b": b, "sb": sb, "chi2": c2},
           "models": {k: chi2(v) for k, v in MODELS.items()},
           "blind_predictions_M1": preds},
          open("data/29_operator_structure.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/29_operator_structure.json")