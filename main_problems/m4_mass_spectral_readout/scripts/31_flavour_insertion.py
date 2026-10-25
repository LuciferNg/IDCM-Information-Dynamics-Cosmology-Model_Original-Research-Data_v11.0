"""M4 docs/31 配套 — 味插入算符之首次實現：樹上局部插入 ⟹ ⟨O⟩ = N_q·|ψ(d_flip)|²
=========================================================================================
目的（`docs/29` §五 V1–V5、`docs/30` §五 之規格）：在 `K_chain ⊗ K_tree` 骨架上**零參數**計
    ⟨O⟩_ch = N_q · P_flip ,   P_flip = 探針夸克在翻轉位點之歸一密度
並與三錨（`⟨O⟩_N = 0.891 ± 0.107`、`_Σ = 1.498 ± 0.262`、`_Ξ = 2.293 ± 0.351`）比對。

模型（全部參數皆為已註冊整數或 ∈ {0,1}，**無新常數**）：
  T1 樹：分支 `B = 4`，深度 `d = 0..L`，殼層簡併 `N_d = B(B−1)^{d−1}`（`N_0 = 1`）。
  T2 徑向哈密頓（沿用 M5 對稱三對角骨架）：`H = diag(V_d) + offdiag(−√(B−1))`，
     歸一化 `Σ_d N_d |v_d|² = 1` ⟹ 單位點密度 `|v_d|²`。
  T3 味位點：味序 `u < d < s` ⟹ 深度 `d_u = 0, d_d = 1, d_s = 2`（**rank = depth**，事前寫死）。
  T4 環境位能：探針（翻轉夸克）感受**其餘兩粒夸克**之核位能
     `V_d = −κ · Σ_{p ∈ passengers} B^{−|d − d_p|}`（核 = 樹核，`κ` = 耦合；零參數取 `κ = 1`）。
  T5 翻轉通道：
       n−p  : 探針 = d（深度 1），環境 = {u(0), d(1)}
       Σ    : 探針 = d（深度 1），環境 = {u(0), s(2)}
       Ξ    : 探針 = d（深度 1），環境 = {s(2), s(2)}
  T6 `⟨O⟩ = N_q · |v_{d_probe}|²`，`N_q = 3`。

紀律：
  Q1 `κ = 1` 為**事前**選擇（核單位 = 跳躍單位）；`κ` 掃描只作**診斷**，若只有特定 `κ` 命中則記為「需調參」= negative。
  Q2 不得改味位點賦值（`rank = depth`）以求命中。
  Q3 不得改 `N_q = 3` 或 `⟨O⟩ = N_q·|v|²` 之形式以求命中。
  Q4 深度 `L` 須夠大且結果對 `L` 收斂（報 L = 6, 8, 10）。
輸出：`data/31_flavour_insertion.json`
"""
import json

import numpy as np

B, NQ = 4, 3
ANCH = {"N": (0.891, 0.107), "Σ": (1.498, 0.262), "Ξ": (2.293, 0.351)}
DEPTH = {"u": 0, "d": 1, "s": 2}
ENV = {"N": ["u", "d"], "Σ": ["u", "s"], "Ξ": ["s", "s"]}


def densities(L, kappa, env):
    """返回 (v_d 向量, 各深度單點密度)"""
    d = np.arange(L + 1)
    V = np.zeros(L + 1)
    for p in env:
        V -= kappa * float(B) ** (-np.abs(d - DEPTH[p]))
    off = np.full(L, -np.sqrt(B - 1))
    H = np.diag(V) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    nodes = np.array([1.0] + [B * (B - 1) ** (k - 1) for k in range(1, L + 1)])
    v = v / np.sqrt((nodes * v ** 2).sum())          # Σ N_d v_d² = 1
    return v, v ** 2


print("=" * 100)
print("味插入算符（零參數首實現）：⟨O⟩ = N_q·|v_{d_flip}|²；κ = 1（核單位 = 跳躍單位）")
print("=" * 100)

print("\n[L 收斂檢查] κ = 1")
out = {}
for L in (6, 8, 10, 12):
    line = f"  L = {L:2d}："
    for sec in ("N", "Σ", "Ξ"):
        _, rho = densities(L, 1.0, ENV[sec])
        O = NQ * rho[DEPTH["d"]]
        out.setdefault(sec, []).append((L, O))
        line += f"  {sec}: ⟨O⟩ = {O:7.4f}"
    print(line)

print("\n[主結果 κ = 1, L = 10]")
res = {}
for sec in ("N", "Σ", "Ξ"):
    _, rho = densities(10, 1.0, ENV[sec])
    O = float(NQ * rho[DEPTH["d"]])
    o, s = ANCH[sec]
    z = (O - o) / s
    res[sec] = {"O": O, "anchor": o, "sigma": s, "z": z}
    print(f"  {sec:>2}：⟨O⟩ = {O:7.4f}   錨 = {o:.3f} ± {s:.3f}   ⟹ {z:+6.2f}σ")

chi2 = sum((res[k]["z"]) ** 2 for k in res)
print(f"\n  χ² = {chi2:.2f} / 3 dof ⟹ {'✅ 命中（<3）' if chi2 < 3 else '🔴 未命中'}")

print("\n[診斷：κ 掃描（只作診斷；若只有特定 κ 命中 ⟹ 需調參 = negative）]")
print(f"{'κ':>8}" + "".join(f"{k:>12}" for k in ("N", "Σ", "Ξ")) + f"{'χ²':>9}")
best = None
for k in (0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0):
    vals, c = [], 0.0
    for sec in ("N", "Σ", "Ξ"):
        _, rho = densities(10, k, ENV[sec])
        O = float(NQ * rho[DEPTH["d"]])
        o, s = ANCH[sec]
        vals.append(O)
        c += ((O - o) / s) ** 2
    print(f"{k:>8.1f}" + "".join(f"{v:>12.4f}" for v in vals) + f"{c:>9.2f}")
    if best is None or c < best[1]:
        best = (k, c)
print(f"  ⟹ 最佳 κ = {best[0]}（χ² = {best[1]:.2f}）")

print("\n[方向性檢查] 若模型正確，⟨O⟩ 應隨 n_s 單調升（N < Σ < Ξ）")
ds = [res[k]["O"] for k in ("N", "Σ", "Ξ")]
print(f"  ⟨O⟩(N, Σ, Ξ) = {ds[0]:.4f}, {ds[1]:.4f}, {ds[2]:.4f} ⟹ "
      f"{'✅ 單調升' if ds[0] < ds[1] < ds[2] else '🔴 非單調'}")

json.dump({"kappa": 1.0, "L": 10, "result": res, "chi2": chi2,
           "kappa_scan_best": {"kappa": best[0], "chi2": best[1]},
           "L_scan": {k: v for k, v in out.items()},
           "model": "H = diag(-κ Σ B^{-|d-d_p|}) + offdiag(-√3); <O> = N_q |v_1|^2"},
          open("data/31_flavour_insertion.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/31_flavour_insertion.json")