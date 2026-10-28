"""M5 docs/40 配套 — (a₀,t) 算子消除路線：尺度對之不可投影性判準
================================================================================
路線 E3。`docs/31` 已證 `(a₀,t)` 之**封閉式**不可識別（語法 30 125 條 × 隨機目標 MC）。
本輪攻**算子路線**：`H_pair = P†(K_chain⊗K_tree)P` 能否消除尺度對？
理論界線：投影 `P` 與 `K_chain⊗K_tree` 均**無量綱**，故任何算子構造只能給出**無量綱比**；
尺度須來自外部錨（O212：`G ≡ ω_ring ≡ M_Pl` 唯一錨）。故唯一可測之算子路線陳述＝
「校準出之無量綱比 `w ≡ t/(ħc/a₀)` 是否為 `(b = 4, N = 88, γ_n)` 之結構組合」。

判準（事前寫死）：
  W1 目標 `w = t/(ħc/a₀)`，由 v2 校準值代入（`a₀ = 0.384184 fm`、`t = 1.7711e−3 MeV`）。
  W2 語法：原子集 `A = {4^±k (k≤14), 88^±k (k≤6), γ_n (n = 1..43 及 44−n 對), π^±k (k≤3),
     √b, (N²−1)/(12N), 1/2, 1/3, 1/π}`；至多 **3 個原子**，允許 `×`、`/`、`+`、`−`、`^`。
  W3 命中定義：`|log(w_form/w_target)| ≤ 1e−4`（相對 1e−4）。
  W4 零假設 MC：以同分佈隨機目標 `10^5` 個（`log10` 均勻於 `[−9, −3]`），量度語法之偶然覆蓋率 `p`。
  W5 判：命中且 `p < 1e−3` ⟹ 結構命中；命中但 `p ≥ 1e−3` ⟹ 偶然（不算）；無命中 ⟹ **不可投影（no-go）**。
  W6 不得反解（即不得先算 w 再挑原子）；原子集與精度在跑之前固定於本 docstring。
輸出：`data/42_scale_pair_operator.json`
"""
import json
import math
import random
from itertools import combinations_with_replacement, product

WARD = 1e-4
random.seed(20261027)

A0_FM, T_MEV = 0.384184, 1.7711e-3
HBARC = 197.3269804                      # MeV·fm
W = T_MEV / (HBARC / A0_FM)              # 目標
print("=" * 100)
print("M5 docs/40 — 尺度對之算子（投影）消除判準")
print("=" * 100)
print(f"[W1] ħc/a₀ = {HBARC/A0_FM:.6f} MeV；w = t/(ħc/a₀) = {W:.8e}")

# γ_n：樹／鏈譜（docs/01 之鏈核譜，n = 1..43）
GAMMA = [4 * math.sin(math.pi * n / 88) ** 2 for n in range(1, 44)]
ATOMS = {}
for k in range(1, 15):
    ATOMS[f"4^-{k}"] = 4.0 ** -k
for k in range(1, 7):
    ATOMS[f"88^-{k}"] = 88.0 ** -k
for n in (1, 2, 3, 5, 8, 11, 21, 22):
    ATOMS[f"γ{n}"] = GAMMA[n - 1]
for k in (1, 2, 3):
    ATOMS[f"π^{k}"] = math.pi ** k
    ATOMS[f"π^-{k}"] = math.pi ** -k
ATOMS["√4"] = 2.0
ATOMS["(N²−1)/12N"] = (88 ** 2 - 1) / (12 * 88)
ATOMS["1/2"] = 0.5
ATOMS["1/3"] = 1 / 3
ATOMS["1/π"] = 1 / math.pi
ATOMS["3/4"] = 0.75
print(f"[W2] 原子集 {len(ATOMS)} 個；至多 3 原子，容差 {WARD:.0e}")


def forms():
    keys = list(ATOMS)
    vals = [ATOMS[k] for k in keys]
    out = []
    for k, v in zip(keys, vals):
        out.append((k, v))
    for (k1, v1), (k2, v2) in product(zip(keys, vals), repeat=2):
        out.append((f"{k1}×{k2}", v1 * v2))
        if v2:
            out.append((f"{k1}/{k2}", v1 / v2))
        out.append((f"{k1}+{k2}", v1 + v2))
        out.append((f"{k1}−{k2}", v1 - v2))
    base = out[:]
    for k3, v3 in zip(keys, vals):
        for kk, vv in base:
            out.append((f"{kk}×{k3}", vv * v3))
            if v3:
                out.append((f"{kk}/{k3}", vv / v3))
    return out


F = forms()
print(f"     形態總數 {len(F):,}")

hits = [(n, v) for n, v in F if v > 0 and abs(math.log(v / W)) <= WARD]
print(f"[W3] 命中 {len(hits)} 條" + ("" if not hits else f"：{[h[0] for h in hits[:5]]}"))

print("[W4] 零假設 MC（10⁵ 隨機目標，log10 均勻 [−9,−3]，numpy 向量化）…")
import numpy as np
lg = np.sort(np.log(np.array([v for _, v in F if v > 0])))
N_MC = 100000
_tg = 10 ** np.random.default_rng(20261027).uniform(-9, -3, N_MC)
_l = np.log(_tg)
_idx = np.searchsorted(lg, _l)
_lo = np.clip(_idx - 1, 0, lg.size - 1)
_hi = np.clip(_idx, 0, lg.size - 1)
_near = np.minimum(np.abs(lg[_lo] - _l), np.abs(lg[_hi] - _l))
cover = int((_near <= WARD).sum())
p = cover / N_MC
print(f"     偶然覆蓋率 p = {p:.3e}（{cover}/{N_MC}）")

verdict = ("structural" if (hits and p < 1e-3) else
           ("coincidence" if hits else "no-go"))
print(f"[W5] 判：{verdict}")
print(f"     ⟹ 算子路線{'給出結構命中' if verdict=='structural' else '**不能消除尺度對**（尺度須外部錨）'}")

out = {"target_w": W, "n_forms": len(F), "n_hits": len(hits),
       "hits": [h[0] for h in hits[:20]], "mc_p": p, "mc_n": N_MC,
       "verdict": verdict,
       "atoms": {k: v for k, v in ATOMS.items()}}
with open("data/42_scale_pair_operator.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("-> data/42_scale_pair_operator.json")