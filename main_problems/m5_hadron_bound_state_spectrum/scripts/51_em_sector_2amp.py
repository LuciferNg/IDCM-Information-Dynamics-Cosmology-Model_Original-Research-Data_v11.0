"""M5 docs/49 配套 — 格點 EM 扇區之完整 2 振幅空間：可算、且已算死
================================================================================
問題：`Δ0` 能否由 M5 **現有**格點核之 EM 扇區（2 個振幅）重現？
即：`Δ0(f) = κ_s·g(0)·Δ(Σq²)(f) + κ_c·g(d_pair(f))·Δ(q_h q_l)(f)`

  · 自能通道（一體）：權重 `g(0)`（**格點上有限**，連續極限才發散 ⟹ 「硬」項之結構來源候選）
  · 對偶通道（兩體）：權重 `g(d_pair)`（精確 Poisson 核，`L·G = I − J/N` 已於 docs/01 精確證）
  · 層指派（已註冊，P10 `docs/00` §6）：`ℓ(u) = ℓ(c) = 0`、`ℓ(d) = ℓ(s) = 7`

判準（事前寫死）：
  E1 電荷不變量於 `e²/9` 單位：`Σq²`、`q_h q_l`；`π⁰` 取理想混合 `(uū − dd̄)/√2`。
  E2 核：`g(d) = [(N²−1) − 6d(N−d)]/(12N)`，`N = 88`（`Fraction` 精確）。
  E3 2 參數最小平方擬合 π/K/D；報殘差與 pull。若 max|pull| > 3 ⟹ EM 扇區**計算上封閉、無解**。
  E4 B：層釘係 **sector-based**（`docs/12` P1：Up(u,c) → 0；Down(d,s) → 7）⟹ `b` 屬 down-type
     ⟹ **`ℓ(b) = 7`（按規則延伸，非 15 成員表所列值）** ⟹ `d_pair(B) = |7−7| = 0` ⟹ B **可表示**。
     故 B 之抑制**不可能係幾何**（`d_pair` 被釘死為 0）⟹ 幾何路線由結構死亡，非由擬合否定。
  E5 若 `g(d)` 可為負（`d ≥ 21`），則符號可翻 ⟹ 紀錄此自由度（但已由 `docs/04` 全枚舉否證）。
  E6 全程 Fraction 核；擬合用 float 但報精確核值。
輸出：`data/51_em_sector_2amp.json`
"""
import json
import math
from fractions import Fraction as F

N = 88


def gx(d):
    return F((N * N - 1) - 6 * d * (N - d), 12 * N)


LAY = {"u": 0, "c": 0, "d": 7, "s": 7}
print("=" * 100)
print("M5 docs/49 — 格點 EM 扇區（2 振幅）可算性")
print("=" * 100)
print(f"\n[E2] g(0) = {gx(0)} = {float(gx(0)):.6f}；g(7) = {gx(7)} = {float(gx(7)):.6f}；"
      f"g(20) = {gx(20)} = {float(gx(20)):.6f}；g(44) = {gx(44)} = {float(gx(44)):.6f}")

# 電荷（e/3 單位）：u=+2, d=-1, s=-1, c=+2, b=-1；反粒子取負
# (族, 中性態 [q_a, q_b], 電荷態 [q_a, q_b], q_h, d_pair 用之兩味)
CASES = [
    ("π", None, None, "u", "d"),      # 中性取理想混合，另算
    ("K", (-1, +1), (+2, +1), "s", "d"),   # K⁰ = d s̄ ; K⁺ = u s̄
    ("D", (+2, -2), (+2, +1), "c", "d"),   # D⁰ = c ū ; D⁺ = c d̄
]
SQ2 = {"π": 0.0}   # π 之 Δ(Σq²) = 5 − 5 = 0（理想混合）
CROSS = {"π": -4.5}  # π⁰ ⟨q_a q_b⟩ = (−4 + −1)/2 = −2.5；π⁺ = +2 ⟹ Δ = −4.5
rows = [{"fam": "π", "dsq": 0.0, "dcross": -4.5, "d_pair": abs(LAY["u"] - LAY["d"]),
         "g_pair": float(gx(abs(LAY["u"] - LAY["d"]))), "d0": -4.5936, "sig": 0.0005}]
for nm, qN, qC, h, l in CASES[1:]:
    dsq = (qN[0] ** 2 + qN[1] ** 2) - (qC[0] ** 2 + qC[1] ** 2)
    dcr = (qN[0] * qN[1]) - (qC[0] * qC[1])
    dpr = abs(LAY[h] - LAY[l])
    rows.append({"fam": nm, "dsq": dsq / 9.0, "dcross": dcr / 9.0, "d_pair": dpr,
                 "g_pair": float(gx(dpr)), "d0": 0.0, "sig": 0.02})

rows[1]["d0"], rows[1]["sig"] = 3.9340, 0.0198
rows[2]["d0"], rows[2]["sig"] = -4.8200, 0.0707

print("\n[E1] 不變量表（e²/9 單位）")
for r in rows:
    print(f"  {r['fam']:2s} ΔΣq² = {r['dsq']:+7.4f}；Δq_hq_l = {r['dcross']:+7.4f}；"
          f"d_pair = {r['d_pair']:2d} ⟹ g = {r['g_pair']:+8.4f}；Δ0 = {r['d0']:+8.4f} ± {r['sig']:.4f}")

# 設計矩陣：Δ0 = κ_s·g(0)·ΔΣq² + κ_c·g_pair·Δq_hq_l
g0 = float(gx(0))
A = [[g0 * r["dsq"], r["g_pair"] * r["dcross"]] for r in rows]
y = [r["d0"] for r in rows]
# 正規方程
S11 = sum(a[0] * a[0] for a in A); S12 = sum(a[0] * a[1] for a in A); S22 = sum(a[1] * a[1] for a in A)
b1 = sum(a[0] * yy for a, yy in zip(A, y)); b2 = sum(a[1] * yy for a, yy in zip(A, y))
det = S11 * S22 - S12 * S12
ks = (b1 * S22 - b2 * S12) / det
kc = (b2 * S11 - b1 * S12) / det
print(f"\n[E3] 2 參數最小平方（π/K/D）：κ_s = {ks:.6f}、κ_c = {kc:.6f}")
print(f"  {'族':2s} {'預言':>10s} {'觀測':>10s} {'殘差':>10s} {'pull':>7s}")
worst = 0.0
for r, a in zip(rows, A):
    pv = ks * a[0] + kc * a[1]
    res = pv - r["d0"]
    pull = abs(res) / r["sig"]
    worst = max(worst, pull)
    print(f"  {r['fam']:2s} {pv:+10.4f} {r['d0']:+10.4f} {res:+10.4f} {pull:7.1f}")
print(f"  ⟹ max|pull| = {worst:.1f} ⟹ "
      f"{'🔴 EM 扇區計算上封閉、無解' if worst > 3 else '相容'}")

print("\n[E4] B：sector pin ⟹ ℓ(b) = 7 ⟹ d_pair(B) = 0（與 K 同幾何）")
LAY_B = {"b": 7}
d_pair_B = abs(LAY_B["b"] - LAY["d"])
print(f"  層釘（docs/12 P1）：Up(u,c) → 0；Down(d,s) → 7 ⟹ b 屬 down-type ⟹ ℓ(b) = 7")
print(f"  ⟹ d_pair(B) = |7 − 7| = {d_pair_B} ⟹ B **可表示**，且鏈幾何與 K 完全相同（q_h = −1/3、|Δq_l| = 1、g(0)）")
print(f"  ⟹ B 之抑制**唔可能係幾何**（d_pair 被釘死 0）⟹ 幾何路線由結構死亡")
print(f"  ⟹ 模型對 B 之鏈項 = 對 K 之鏈項 ⟹ 若 K 需 4.4492 MeV 之一體項，B 何以缺之？此為真正問題。")
target = 0.3100 / 3.9340
print(f"  （附：若硬要幾何解釋，需 g(d)/g(0) = {target:.5f} ⟹ 連續解 d ≈ 16.67，與 pin 之 d = 0 相斥 ⟹ 不可能）")

print("\n[E5] 核符號自由度")
zs = [d for d in range(45) if float(gx(d)) < 0]
print(f"  g(d) < 0 之格距：d ∈ [{min(zs)}, {max(zs)}] ⟹ 符號可翻（但 `docs/04` 已全枚舉否證含號核類）")

out = {"g0": g0, "g7": float(gx(7)), "kappa_s": ks, "kappa_c": kc,
       "max_pull": worst, "b_representable": True, "b_d_pair": d_pair_B,
       "b_need_ratio": target, "rows": rows}
with open("data/51_em_sector_2amp.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/51_em_sector_2amp.json")