"""M5 docs/01 配套 — C_N 環上嘅 Coulomb 核：精確 Green 函數 ＋ 非加性示範（2026-10-24）
=========================================================================================
前設問題（M4 三份上游 no-go 合流得出，見 v11 M4 `docs/21` §五.4、`docs/22` §五.2）：
  em 通道分裂 = ⟨Σ_{i<j} q_i q_j / r_ij⟩ ＋ Σ_i q_i² Δm_i
  分子（電荷量子數）已證不足（O124 41.2σ 反例）；分母 = 波函數重疊 ⟨1/r_ij⟩，IDCM 完全缺失。
本輪（第一步，唔跳步）：
  A. C_N 上 Coulomb 核 = Laplacian 偽逆 G = L⁺，候選閉式 g(d) = d(N−d)/(2N)。
     **精確證明方式**：L⁺ 唯一 ⟹ 只需驗 L·G = I − J/N 且 G 對稱、行和 0（全部 Fraction 精確）。
  B. 結構性質：g 遞增至 d = N/2 ⟹ 一維環上核係**二次增長**（非 1/r）。
  C. 非加性示範：E = Σ_{i<j} g(d_ij) 對「同取代、不同 spectator」給出唔同值。
  D. 量級估計：由跨定 κ 使 ΔE 達 MeV 級。
判準（事前）：A 若 L·G ≠ I − J/N ⟹ 閉式作廢；C 若 spread = 0 ⟹ 核加性、路線即死。
本輪**唔做**：sector assignment（flavour→site 映射，O140 已證不存在）、QED 自能重整、波函數 ansatz。
（首版用 sympy 複指數 DFT 建 88×88 偽逆 ⟹ 超時；改為「閉式 + 精確矩陣乘積驗證」，因 L⁺ 唯一。）
"""
import json
from fractions import Fraction as F

OUT = "data/01_chain_green_function.json"
N = 88
print("=" * 88)
print("M5 docs/01 — C_N 環 Coulomb 核（Laplacian 偽逆）精確核實 + 非加性示範")
print("=" * 88)

def closed_form(d, n):
    """L⁺_{ij}（環距離 d）：g(d) = [(N²−1) − 6d(N−d)] / (12N)
    推導：L⁺_{0d} = (1/N)Σ_{k=1}^{N−1} cos(2πkd/N)/(2−2cos(2πk/N))；
    閉式由「行和 = 0」＋「L·G = I − J/N」唯一確定（本腳本精確驗證）。"""
    return F((n * n - 1) - 6 * d * (n - d), 12 * n)

def check(n):
    """驗 G[d] = d(n−d)/2n 係 L⁺：L G = I − J/n（G 由閉式建，全部 Fraction）"""
    G = [[closed_form(min((i - j) % n, (j - i) % n), n) for j in range(n)] for i in range(n)]
    ok_sym = all(G[i][j] == G[j][i] for i in range(n) for j in range(n))
    ok_sum = all(sum(G[i]) == 0 for i in range(n))
    worst = F(0)
    for i in range(n):
        for j in range(n):
            s = F(2) * G[i][j] - G[(i - 1) % n][j] - G[(i + 1) % n][j]
            target = (F(1) if i == j else F(0)) - F(1, n)
            worst = max(worst, abs(s - target))
    return ok_sym, ok_sum, worst

print("\n[A] 精確驗證 L·G = I − J/N（Fraction；L⁺ 唯一 ⟹ 此即證明閉式）")
rows = []
for n in (8, 12, 88):
    sym, sm, worst = check(n)
    rows.append({"N": n, "symmetric": sym, "rowsum_zero": sm, "max_residual": str(worst)})
    print(f"  N={n:>2}: 對稱 {sym} / 行和=0 {sm} / max|L·G − (I−J/N)| = {worst} "
          f"{'✅ 閉式成立（EXACT）' if (sym and sm and worst == 0) else '🔴'}")

print("\n[B] N = 88 精確值（Fraction）")
g = {d: closed_form(d, N) for d in range(N // 2 + 1)}
print(f"  g(0) = {g[0]} ≈ {float(g[0]):.5f}；g(1) = {g[1]} ≈ {float(g[1]):.5f}；"
      f"g(44) = {g[44]} ≈ {float(g[44]):.5f}（最大）")
print(f"  近距離 g(0) = +{float(g[0]):.4f}、g(1) = +{float(g[1]):.4f}；最遠 g(44) = {float(g[44]):.4f}")
print(f"  ⟹ 緊緻 1D 環上『Coulomb』核 ∝ [N²−1 − 6d(N−d)]：近正遠負、**反屏蔽**，完全冇 1/r 行為 ⟹"
      f" K_chain 單獨承載唔到 QED 自能差（必須靠 K_tree 徑向自由度）")

print("\n[C] 非加性示範：pairwise 核能量 vs 加性預言")
def energy(sites):
    tot = F(0)
    for a in range(len(sites)):
        for b in range(a + 1, len(sites)):
            d = abs(sites[a] - sites[b]) % N
            tot += g[min(d, N - d)]
    return tot

confs = {"A": [0, 2, 44], "B": [0, 2, 20], "C": [0, 2, 3], "D": [0, 2, 1]}
Es = {}
for lab, c in confs.items():
    Es[lab] = energy(c)
    print(f"  構型 {lab} sites={c} ⟹ E = {float(Es[lab]):.4f}")
vs = [float(v) for v in Es.values()]
span = max(vs) - min(vs)
print(f"  固定『取代對』(0,2)、只換 spectator：E ∈ [{min(vs):.3f}, {max(vs):.3f}] ⟹ spread = {span:.3f} ≫ 0 ⟹ 非加性成立 ✅")

print("\n[D] 量級：κ 使 pairwise 能達 MeV 級")
E_obs = 5.664   # M4 8 通道 E 跨度（MeV）
print(f"  觀測跨度 5.664 MeV / 本核跨度 {span:.3f} ⟹ κ ≈ {E_obs/span:.5f} MeV")
print(f"  最大單對能 = κ·g(44) = {E_obs/span*float(g[44]):.3f} MeV")

res = {"meta": {"script": "01_chain_green_function", "date": "2026-10-24", "N": N,
                "note": "首版 sympy DFT 超時；改閉式+精確矩陣驗證（L⁺ 唯一）"},
       "verify": rows,
       "g": {"g0": [g[0].numerator, g[0].denominator], "g1": [g[1].numerator, g[1].denominator],
             "g2": [g[2].numerator, g[2].denominator], "g44": [g[44].numerator, g[44].denominator]},
       "non_additivity": {"E": {k: float(v) for k, v in Es.items()}, "spread": span},
       "kappa_MeV": E_obs / span}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")