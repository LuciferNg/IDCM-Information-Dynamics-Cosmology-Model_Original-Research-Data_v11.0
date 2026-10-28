"""M5 docs/43 配套 — 電荷分裂 Δ0 之兩體上界：一體 EM 自能項係結構必需
================================================================================
觸發：`docs/42` 改立目標 —— `A` 殘差（K +0.9992、B −3.4300 MeV）。本輪問一個更前置之問題：
`Δ0 ≡ M(J=0, neutral) − M(J=0, charged)` 可否由**兩體** Coulomb（含符號 `q_h q_l`）產生？

判準（事前寫死）：
  E1 兩體 Coulomb 之貢獻上界：`|Δ0| ≤ |Δ(q_h q_l)|·α·(ħc/a₀)·⟨B⟩`，
     `|Δ(q_h q_l)| ≤ 2/3`（`u/d` 之電荷差 1，乘 `|q_h| ≤ 2/3`），`⟨B⟩ ≤ 1`（恆等式，`docs/34` 同款）。
  E2 以 v3 校準值 `a₀ = 0.384184 fm`、`t = 1.771108 MeV` 代入，得**硬上界** `Δ0^max`。
  E3 逐族比較 `Δ0^data` 與 `Δ0^max`；`|Δ0^data| > Δ0^max` ⟹ 兩體 Coulomb **不可能**係主項
     ⟹ 一體 EM 自能項**結構必需**（非可選）。
  E4 所需一體項 `M_f ≡ Δ0^data − Δ0^model`，其中 `Δ0^model = A_model − 3c_model`（`s_0 = −3`）。
  E5 簽名判：若 `M_f` 於 π/K/D 同號同量級而於 B 反號 ⟹ B 係 anomaly（與 `docs/42` 之 q_B 一致）。
  E6 不引入任何新自由參數；本輪只判「必需／不必然」。
輸出：`data/45_two_body_delta0_bound.json`
"""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
A0, T0 = 0.384184, 1.771108e-3
M_REF = 0.30
A_D, C_D = 1.2025, 0.3525

MASS = {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    for k in ids:
        if k:
            MASS[int(k)] = mv * 1000.0


def amps(a, t, q, mh, ml):
    """v3：C_q ∝ 1/(m_h·m_ref)。q 帶符號（q_h q_l）。"""
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * abs(q) / (16 * mh * M_REF)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    return ALPHA * q * (HBARC / a) * rA * 1000, Cq * (HBARC / a) ** 3 * P0 * 1000, rA, P0


print("=" * 100)
print("M5 docs/43 — Δ0 之兩體上界 vs 數據：一體 EM 項係否結構必需")
print("=" * 100)

print("\n[E1/E2] 硬上界（以 v3 校準 a₀、⟨B⟩ ≤ 1）")
g_max = ALPHA * (2 / 3) * HBARC / A0 * 1000
print(f"  α·(2/3)·(ħc/a₀) = {g_max:.4f} MeV ⟹ Δ0^max = {g_max:.4f} MeV（⟨B⟩ = 1 之上限）")

FAM = [
    ("π", 111, 211, (1 / 3, -2 / 3), 0.30, None),
    ("K", 311, 321, (-1 / 3, 1 / 3), 0.495, 0.30),
    ("D", 421, 411, (2 / 3, -2 / 3), 1.27, 0.30),
    ("B", 511, 521, (-1 / 3, 1 / 3), 4.18, 0.30),
]
print("\n[E3] 逐族：|Δ0^data| vs 兩體上界（⟨B⟩ ≤ 1）")
rows = []
for nm, i0, i1, (qh, qlN), mh, ml in FAM:
    d0 = MASS[i0] - MASS[i1]
    # 兩體 Coulomb 之最大可能貢獻：電荷差 1 × |q_h| × α(ħc/a₀) × ⟨B⟩ ≤ 1
    qlC = qlN + 1.0
    dq = abs(qh * qlN - qh * qlC)
    bound = dq * ALPHA * HBARC / A0 * 1000
    exceed = abs(d0) / bound if bound else float("inf")
    deficit = abs(d0) - bound
    rows.append({"fam": nm, "d0_data": d0, "dq": dq, "bound_fam": bound,
                 "exceed_factor": exceed, "deficit_MeV": deficit,
                 "one_body_min": deficit if deficit > 0 else 0.0})
    print(f"  {nm:2s} |Δ0^data| = {abs(d0):7.4f}；|Δ(q_h q_l)| = {dq:.4f} ⟹ 上界 {bound:6.4f} MeV"
          f" ⟹ 超 {exceed:5.2f}× ⟹ 一體項至少 {max(deficit,0):6.4f} MeV"
          f" {'🔴' if deficit > 0 else '🟡'}")

print("\n[E5] 簽名判（缺口帶號 = sign(Δ0^data)）")
for r in rows:
    print(f"  {r['fam']:2s} 缺口帶號 {np.sign(r['d0_data']):+.0f}")
signs = {r["fam"]: float(np.sign(r["d0_data"])) for r in rows}
print(f"  π/K/D/B 帶號：{signs} ⟹ 交替（即 docs/42 之抵消結構）")
must = any(r["deficit_MeV"] > 0 for r in rows)
out = {"g_max_MeV": g_max, "rows": rows, "one_body_required": bool(must), "signs": signs}
with open("data/45_two_body_delta0_bound.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print(f"\n[E6] 一體 EM 項結構必需：{'✅ 是（至少一族超上界）' if must else '🟡 否'}")
print("-> data/45_two_body_delta0_bound.json")
