"""M4 docs/31 配套 — 味插入算符：列舉式首實現 ＋ MC null（防掃描假陽性）
=========================================================================================
`docs/30` §五 之規格要求零參數計 `⟨O⟩ = N_q·P_flip`。首試（rank = depth、per-site 歸一、κ = 1）
得 ⟨O⟩ ~ 0.005–0.06（細兩個數量級）且**方向相反**（Ξ < N）。本輪按紀律改為**列舉 + null**：

列舉空間（全部事前可列，非事後挑）：
  E1 味→深度賦值：`{u, d, s} → {0, 1, 2}` 之 **6 個置換**（heavy-deeper 與其反向皆含在內）。
  E2 位點密度歸一：`per-site |v_d|²` 抑或 `shell-summed N_d|v_d|²`（2 種）。
  E3 耦合 `κ ∈ {0.5, 1, 2, 3, 5, 10, 30, 100, 300}`（9 種）。
  ⟹ 總 6×2×9 = **108 個模型**（全部零新常數；`B = 4`、`N_q = 3`、核 = `B^{−|Δd|}`）。
判準：
  J1 逐模型 χ²（3 錨）。報最佳者，但**不得**只報最佳而不報掃描次數。
  J2 **MC null**：抽 3 個 iid log-uniform `⟨O⟩ ∈ [0.2, 4]`（配同一組 σ），在 108 個模型上取 min χ²；
      重複 2000 次 ⟹ `p = P(min χ²_null ≤ min χ²_obs)`。**若 p ≳ 0.05 則命中不可採信**。
  J3 方向性：命中之模型須給 `N < Σ < Ξ` 單調（否則與資料結構相反，即使 χ² 低亦不採）。
輸出：`data/31_flavour_insertion_enum.json`
"""
import itertools
import json

import numpy as np

B, NQ = 4, 3
ANCH = {"N": (0.891, 0.107), "Σ": (1.498, 0.262), "Ξ": (2.293, 0.351)}
FLAV = ("u", "d", "s")
# 通道之環境（探針 = 翻轉夸克 d）
ENV = {"N": ("u", "d"), "Σ": ("u", "s"), "Ξ": ("s", "s")}
KAPPA = (0.5, 1, 2, 3, 5, 10, 30, 100, 300)
L = 12


def density(dep, kappa, env, res="site"):
    d = np.arange(L + 1)
    V = np.zeros(L + 1)
    for p in env:
        V -= kappa * float(B) ** (-np.abs(d - dep[p]))
    off = np.full(L, -np.sqrt(B - 1))
    H = np.diag(V) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    nodes = np.array([1.0] + [B * (B - 1) ** (k - 1) for k in range(1, L + 1)])
    v = v / np.sqrt((nodes * v ** 2).sum())
    rho = v ** 2
    return rho if res == "site" else nodes * rho


def chi2_of(triple):
    return sum(((triple[i] - ANCH[k][0]) / ANCH[k][1]) ** 2 for i, k in enumerate(("N", "Σ", "Ξ")))


print("=" * 100)
print(f"味插入算符列舉：6 置換 × 2 歸一 × {len(KAPPA)} κ = {6*2*len(KAPPA)} 模型（L = {L}）")
print("=" * 100)
models = []
for perm in itertools.permutations(range(3)):
    dep = {FLAV[i]: perm[i] for i in range(3)}
    for res in ("site", "shell"):
        for k in KAPPA:
            trip = tuple(float(NQ * density(dep, k, ENV[s], res)[dep["d"]]) for s in ("N", "Σ", "Ξ"))
            models.append({"dep": {f: int(v) for f, v in dep.items()}, "norm": res, "kappa": k,
                           "triple": trip, "chi2": chi2_of(trip)})

models.sort(key=lambda m: m["chi2"])
print(f"{'排名':>4} {'u,d,s 深度':>12} {'歸一':>7}{'κ':>7}{'⟨O⟩(N,Σ,Ξ)':>28}{'χ²':>9}{'單調':>6}")
for i, m in enumerate(models[:8], 1):
    t = m["triple"]
    mono = t[0] < t[1] < t[2]
    print(f"{i:>4} {str(tuple(m['dep'][f] for f in FLAV)):>12} {m['norm']:>7}{m['kappa']:>7}"
          f"{t[0]:>10.3f}{t[1]:>8.3f}{t[2]:>9.3f}{m['chi2']:>9.2f}{('✅' if mono else '—'):>6}")

best = models[0]
print(f"\n[J1] 最佳：{best['dep']} × {best['norm']} × κ = {best['kappa']} ⟹ χ² = {best['chi2']:.2f}，"
      f"單調 = {best['triple'][0] < best['triple'][1] < best['triple'][2]}")

# J3 方向性
mono_models = [m for m in models if m["triple"][0] < m["triple"][1] < m["triple"][2]]
print(f"[J3] 單調升（N<Σ<Ξ）之模型：{len(mono_models)}/{len(models)}")
if mono_models:
    bm = min(mono_models, key=lambda m: m["chi2"])
    print(f"     其中最佳：{bm['dep']} × {bm['norm']} × κ = {bm['kappa']} ⟹ χ² = {bm['chi2']:.2f}，"
          f"⟨O⟩ = ({bm['triple'][0]:.3f}, {bm['triple'][1]:.3f}, {bm['triple'][2]:.3f})")

# J2 MC null
rng = np.random.default_rng(20261025)
sig = np.array([ANCH[k][1] for k in ("N", "Σ", "Ξ")])
grid = np.array([m["triple"] for m in models])            # (M, 3)
NMC = 2000
mins = np.empty(NMC)
for i in range(NMC):
    fake = np.exp(rng.uniform(np.log(0.2), np.log(4.0), 3))
    mins[i] = (((grid - fake[None, :]) / sig[None, :]) ** 2).sum(axis=1).min()
p_null = float((mins <= best["chi2"]).mean())
print(f"\n[J2] MC null（{NMC} 次，{len(models)} 模型）：min χ² 之中位 = {np.median(mins):.2f}，"
      f"5% 分位 = {np.percentile(mins,5):.2f}")
print(f"     p = P(min χ²_null ≤ {best['chi2']:.2f}) = {p_null:.3f} ⟹ "
      f"{'🔴 掃描即可得，不可採信' if p_null > 0.05 else '✅ 通過 null'}")

json.dump({"L": L, "n_models": len(models), "best": best,
           "best_monotone": (min(mono_models, key=lambda m: m["chi2"]) if mono_models else None),
           "n_monotone": len(mono_models), "p_null": p_null,
           "null_median": float(np.median(mins))},
          open("data/31_flavour_insertion_enum.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/31_flavour_insertion_enum.json")