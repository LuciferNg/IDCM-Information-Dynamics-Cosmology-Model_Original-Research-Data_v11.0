"""M5 docs/47 配套 — B 異常：正確歸一化（sign 而非 |q_h|）＋ 抑制機制候選檢
================================================================================
自糾：`docs/46` 用 `r_f ≡ Δ0/(−q_h)`，把 `|q_h| = 2/3`（π、D）與 `1/3`（K、B）兩族混歸一化，
令 K 之 `r_f` 被人為放大 2×。正確歸一化係按**符號**：`A_f ≡ Δ0/(−sign q_h)`。

判準（事前寫死）：
  B1 `A_f = Δ0/(−sign q_h)`；`q_h`：`c = +2/3`、`s = b = −1/3`、π 取 `+2/3`。
  B2 若 π/K/D 之 `A_f` 相對離散 ≤ 15% ⟹ 通用標度成立（`A*`）；B 另判。
  B3 `B` 之偏離以 `z = |A_B − A*|/σ(A_B)` 報（`σ(A_B) = σ(Δ0_B)`，PDG）。
  B4 抑制因子 `s ≡ A_B/A*`；候選機制逐一比對，命中需 `|ln(候選/s)| ≤ ln(1.10)`（10% 內）。
  B5 候選集**先固定**（不得反解後再加）：`m_c/m_b`、`(m_c/m_b)²`、`m_s/m_b`、`m_s/m_c`、`(m_c/m_b)³`、
     `ħc/(m_b a₀)`、`[ħc/(m_b a₀)]²`、`(m_s/m_b)²`、`1/√N_d`、`1/N_d`。
  B6 不得以擬合係數收工；本輪只判「有無候選」。
輸出：`data/49_B_anomaly.json`
"""
import json
import math
from pathlib import Path

HBARC = 0.1973269804
A0 = 0.384184
M_ = {}
ERR = {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    try:
        ev = abs(float(line[52:60]))
    except ValueError:
        ev = 0.0
    for k in ids:
        if k:
            M_[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0

FAM = [("π", 111, 211, +1), ("K", 311, 321, -1), ("D", 421, 411, +1), ("B", 511, 521, -1)]
print("=" * 100)
print("M5 docs/47 — B 異常：sign 歸一化 ＋ 抑制候選")
print("=" * 100)

print("\n[B1] A_f = Δ0/(−sign q_h)（sign 歸一化）")
rows = []
for nm, i0, i1, sq in FAM:
    d0 = M_[i0] - M_[i1]
    sd = math.hypot(ERR[i0], ERR[i1])
    A = d0 / (-sq)
    rows.append({"fam": nm, "d0": d0, "sigma_d0": sd, "sign_qh": sq, "A": A})
    print(f"  {nm:2s} Δ0 = {d0:+8.4f} ± {sd:.4f}；sign(q_h) = {sq:+d} ⟹ A_f = {A:+8.4f} MeV")

pi_k_d = [r["A"] for r in rows if r["fam"] in ("π", "K", "D")]
mean = sum(pi_k_d) / 3
sd3 = math.sqrt(sum((v - mean) ** 2 for v in pi_k_d) / 2)
print(f"\n[B2] π/K/D：A* = {mean:.4f} ± {sd3:.4f} MeV（相對離散 {100*sd3/mean:.2f}%）"
      f" ⟹ {'✅ 通用標度成立' if sd3 / mean <= 0.15 else '🔴 不成立'}")
for r in rows:
    if r["fam"] in ("π", "K", "D"):
        print(f"     {r['fam']:2s} 偏離 A*：{100*(r['A']/mean-1):+6.2f}%")

B = [r for r in rows if r["fam"] == "B"][0]
z = abs(B["A"] - mean) / B["sigma_d0"]
print(f"\n[B3] B：A_B = {B['A']:.4f} ± {B['sigma_d0']:.4f} vs A* = {mean:.4f} ± {sd3:.4f}"
      f" ⟹ 以 A* 預言 {mean:.4f}，偏離 {z:.1f}σ")

s_obs = B["A"] / mean
print(f"\n[B4] 抑制因子 s ≡ A_B/A* = {s_obs:.5f}（1/{1/s_obs:.2f}）")

print("\n[B5] 候選機制（先固定集，10% 容差）")
CAND = {
    "m_c/m_b": 1.27 / 4.18,
    "(m_c/m_b)²": (1.27 / 4.18) ** 2,
    "m_s/m_b": 0.495 / 4.18,
    "m_s/m_c": 0.495 / 1.27,
    "(m_c/m_b)³": (1.27 / 4.18) ** 3,
    "ħc/(m_b·a₀)": HBARC / (4.18 * A0),
    "[ħc/(m_b·a₀)]²": (HBARC / (4.18 * A0)) ** 2,
    "(m_s/m_b)²": (0.495 / 4.18) ** 2,
    "1/√N_d": 1 / math.sqrt(88),
    "1/N_d": 1 / 88,
}
hits = []
for k, v in CAND.items():
    dev = abs(math.log(v / s_obs))
    ok = dev <= math.log(1.10)
    if ok:
        hits.append(k)
    print(f"  {k:18s} = {v:0.5f}（vs s = {s_obs:.5f}）⟹ 偏離 {100*(v/s_obs-1):+8.1f}%"
          f" {'✅ 命中' if ok else ''}")
print(f"\n[B6] 命中候選：{hits if hits else '無'}")

print("\n[B7] 逆推律：由 D 預言 B／K（`ΔM ∝ q_h/m_h^b`，b = 1（HQET）與 b = 2）")
mh = {"K": 0.495, "D": 1.27, "B": 4.18}
qhd = {"K": -1 / 3, "D": +2 / 3, "B": -1 / 3}
dD = [r for r in rows if r["fam"] == "D"][0]["d0"]
pred = {}
for b in (1, 2):
    print(f"  b = {b}：")
    for f in ("K", "B"):
        pv = dD * (qhd[f] / qhd["D"]) * (mh["D"] / mh[f]) ** b
        obs = [r for r in rows if r["fam"] == f][0]["d0"]
        sg = [r for r in rows if r["fam"] == f][0]["sigma_d0"]
        zz = abs(pv - obs) / sg
        pred[f"b{b}_{f}"] = {"pred": pv, "obs": obs, "z": zz}
        print(f"     {f}: 預言 {pv:+7.4f}；觀測 {obs:+7.4f} ± {sg:.4f} ⟹ {zz:5.1f}σ")

out = {"rows": rows, "hqet_pred": pred, "A_star": mean, "sigma_A_star": sd3, "B_A": B["A"], "B_z": z,
       "suppression": s_obs, "candidates": CAND, "hits": hits, "hqet_pred": pred}
with open("data/49_B_anomaly.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/49_B_anomaly.json")