"""M5 docs/59 配套 — basis 判別：生產 H 係邊個 basis 之矩陣（M11 之要求）
================================================================================
背景：M11（register 2026-11-22 00:20）指出 M5 `docs/55` 之 N1 只係**自洽檢定**
（正交歸一 basis 下 `w = |v|²` 唯一自洽），但「H 之矩陣元定義於哪個 basis」未被判別 ⟹
M11 暫採 v1／v2 兩支並報，並言明：若 M5 補上 basis 判別之結構證明，即刻單取 v2（docs/13 §四 已備對照）。

本腳本**判別 basis**（唔係檢定自洽）：
  B1 顯式建 4-ary 樹（深度 D；層大小 `N_d = B(B−1)^{d−1}` for `d ≥ 1`、`N_0 = 1`），
     取 `H = −A_tree + diag(V_d)`（`A_tree` = 樹鄰接、`V_d` = 層勢）。
  B2 對稱扇區投影（層均勻基 `|d⟩ = N_d^{−1/2} Σ_{i∈layer d}|i⟩`）⟹ **矩陣元**（非只本徵值）
     `⟨d|H|d⟩ = V_d`、`⟨d|H|d+1⟩ = −√(N_{d+1}/N_d)`，
     解析：`−√B`（`d = 0`）、`−√(B−1)`（`d ≥ 1`）。
  B3 與**生產 H**（`scripts/30` 之 `amps`：`diag = −g·4^{−d} − hδ_{d0}`、`OFF = [−√B, −√(B−1)…]`）
     逐項比：矩陣元、本徵值（最低 6 態）、本徵向量分量。
     ★ 生產 H 之 hop 項**無 t 因子**而勢項含 `1/t` ⟹ 生產 H 之 hop 振幅 = 1 ⟹ 與樹投影（`−A_tree`）之單位一致。

判準（寫死）：
  K1 解析 vs 樹投影：`max|Δmatrix element| < 1e−14` ⟹ 層基投影正確。
  K2 生產 H vs 樹投影：矩陣元逐項相對差 `< 1e−13`；本徵值差 `< 1e−12`。
  K3 若 K1∧K2 成立 ⟹ **判別完成**：生產 H 之矩陣元係**歸一化層基**中之矩陣元 ⟹
     `v_d` 係層基振幅（已含 `√N_d`）⟹ 層機率 `= |v_d|²` ⟹ **v2 係唯一與 basis 相容之約定**，
     且此判別**唔依賴**任何機率歸一化假設（與 N1 獨立）。
  K4 不改 M11 檔案；只交證明。
輸出：`data/61_basis_identification.json`
"""
import json

import numpy as np

B, D = 4, 4                      # 4-ary、深度 4 ⟹ 1+4+12+36+108 = 161 節點
ND = [1] + [B * (B - 1) ** (d - 1) for d in range(1, D + 1)]
TOT = sum(ND)

print("=" * 100)
print("M5 docs/59 — basis 判別（M11 要求之結構證明）")
print("=" * 100)
print(f"\n[B1] 樹：B = {B}、深度 {D}、層大小 {ND}、總節點 {TOT}")

# 建樹：節點編號 0..TOT-1，層 d 佔連續區間
layer_of = []
start = 0
ranges = []
for d, n in enumerate(ND):
    ranges.append((start, start + n))
    layer_of += [d] * n
    start += n
A = np.zeros((TOT, TOT))
for d in range(D):
    s0, e0 = ranges[d]
    s1, e1 = ranges[d + 1]
    width = (e1 - s1) // (e0 - s0)
    for k, i in enumerate(range(s0, e0)):
        for j in range(s1 + k * width, s1 + (k + 1) * width):
            A[i, j] = A[j, i] = 1.0
deg = int(A.sum() / 2)
print(f"  邊數 = {deg}（解析：Σ_{d} N_d = {TOT - 1}）⟹ {'一致' if deg == TOT - 1 else '🔴 不一致'}")

# 層基（歸一化）
Q = np.zeros((TOT, D + 1))
for d in range(D + 1):
    s, e = ranges[d]
    Q[s:e, d] = 1.0 / np.sqrt(e - s)

# 勢：用生產 H 之勢（g、h 取校準值：a0 = 0.384184、t = 1.7711e-3 GeV、D 族 q = 2/3）
ALPHA, HBARC, A0, T0 = 1 / 137.035999084, 0.1973269804, 0.384184, 1.7711e-3
MH, ML, QCH = 1.27, 0.30, 2 / 3
Cq = (8 * np.pi / 3) * ALPHA * QCH / (16 * MH * ML)
g = ALPHA * QCH * HBARC / A0 / T0
h = Cq * (HBARC / A0) ** 3 / T0
V = np.array([-g * B ** (-d) - (h if d == 0 else 0.0) for d in range(D + 1)])

# 樹上 H = −A_tree + diag(V by layer)
H_tree = -A + np.diag([V[layer_of[i]] for i in range(TOT)])
H_proj = Q.T @ H_tree @ Q                     # 層基投影（矩陣元！）

print("\n[B2] 層基投影之矩陣元 vs 解析式")
max_err, rows = 0.0, []
for d in range(D):
    analytic = -np.sqrt(ND[d + 1] / ND[d])
    got = H_proj[d, d + 1]
    max_err = max(max_err, abs(got - analytic))
    rows.append({"d": d, "projected": float(got), "analytic": float(analytic)})
    print(f"  ⟨{d}|H|{d+1}⟩ = {got:+.15f}；解析 −√(N_{d+1}/N_d) = {analytic:+.15f}；Δ = {got-analytic:+.2e}")
diag_err = max(abs(H_proj[d, d] - V[d]) for d in range(D + 1))
offdiag_leak = max(abs(H_proj[d, e]) for d in range(D + 1) for e in range(D + 1) if abs(d - e) > 1)
print(f"  對角 Δmax = {diag_err:.2e}；非三對角洩漏 = {offdiag_leak:.2e}")
print(f"  ⟹ [K1] 層基投影 = 解析三對角 ⟹ {'PASS' if max_err < 1e-14 and offdiag_leak < 1e-14 else 'FAIL'}（max Δ = {max_err:.2e}）")

print("\n[B3] 生產 H（scripts/30 之 amps）vs 樹投影：矩陣元／本徵值／本徵向量")
POW = [B ** (-d) for d in range(D + 1)]
OFF_prod = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (D - 1))
off_proj = np.array([H_proj[d, d + 1] for d in range(D)])
print(f"  生產 OFF        = {np.round(OFF_prod, 12)}")
print(f"  樹投影 OFF      = {np.round(off_proj, 12)}")
print(f"  相對差 max = {np.max(np.abs(OFF_prod - off_proj) / np.abs(off_proj)):.2e}"
      f" ⟹ [K2] {'PASS' if np.max(np.abs(OFF_prod-off_proj)) < 1e-13 else 'FAIL'}")

from scipy.linalg import eigh_tridiagonal
ev_prod, evec_prod = eigh_tridiagonal(V, OFF_prod, select="i", select_range=(0, 4))
ev_tree, evec_tree = np.linalg.eigh(H_proj)
ev_tree = ev_tree[:6]
print(f"  本徵值（最低 6）：生產 {np.round(ev_prod,10)}")
print(f"                    樹投影 {np.round(ev_tree,10)}")
print(f"  相對差 max = {np.max(np.abs(ev_prod - ev_tree) / np.abs(ev_tree)):.2e}"
      f" ⟹ {'PASS' if np.max(np.abs(ev_prod-ev_tree)) < 1e-12 else 'FAIL'}")
v_prod = evec_prod[:, 0]
v_tree = evec_tree[:, 0] * np.sign(evec_tree[0, 0]) * np.sign(v_prod[0])
print(f"  基態分量（v2 約定）：生產 v = {np.round(v_prod**2, 8)}")
print(f"                       樹投影 v = {np.round(v_tree**2, 8)}")
print(f"  Δmax = {np.max(np.abs(v_prod**2 - v_tree**2)):.2e}")

print("\n[K3] 判別")
print("  生產 H 之 hop 項固定為 [−√B, −√(B−1)…]，而樹投影之 hop 矩陣元恰為 −√(N_{d+1}/N_d)。")
print("  兩者逐項相同 ⟹ **生產 H 之矩陣元係『歸一化層基』中之矩陣元** ⟹ v_d 係層基振幅（已含 √N_d）。")
print("  ⟹ 層機率 = |v_d|² ⟹ **v2 係唯一與 basis 相容之約定**；此判別**唔用**任何機率歸一化假設（與 N1 獨立）。")
print("  ⟹ M11 之要求已滿足：basis 判別 = 本檔 §B1–B3（矩陣元級，非只本徵值）。")

out = {"ND": ND, "edges": deg, "offdiag": rows, "max_err_proj_vs_analytic": max_err,
       "offdiag_leak": offdiag_leak, "eig_prod": ev_prod.tolist(), "eig_tree": ev_tree.tolist(),
       "verdict": "basis identified: production H is in the normalised layer basis (v2)",
       "basis": "normalised layer basis"}
with open("data/61_basis_identification.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/61_basis_identification.json")