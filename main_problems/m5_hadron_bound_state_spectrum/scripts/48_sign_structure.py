"""M5 docs/46 配套 — Δ0 之符號結構：兩體 Coulomb 通道 vs 數據（4 族）
================================================================================
`docs/43` 已證量級不足（π/K/D 超上界 1.9–3.7×）。本輪問更硬嘅問題：**符號**係否一致？
`V ∝ α q_h q_l / r` ⟹ 較吸引者較輕 ⟹ `Δ0 = M(neutral) − M(charged)` 之符號可由電荷乘積決定。

判準（事前寫死）：
  S1 電荷指派：`q_h, q_l` 以 `e/3` 為單位（`u = +2, d = −1, s = −1, c = +2, b = −1`；反粒子取負）。
     `π⁰` 取 `(uū − dd̄)/√2` ⟹ 有效乘積 `⟨q_h q_l⟩ = (−4/9 − 1/9)/2 = −5/18`。
  S2 模型符號：`sign(Δ0^model) = sign(⟨q_hq_l⟩_neutral − ⟨q_hq_l⟩_charged)`（越吸引 ⟹ 越低 ⟹ 越輕）。
  S3 數據符號：`sign(Δ0^data) = sign(M(neutral) − M(charged))`，PDG 2024 直讀。
  S4 判：4 族中任一符號相反 ⟹ 兩體 Coulomb **連符號都錯** ⟹ 該通道徹底出局（不只量級不足）。
  S5 另報「歸一化分裂」`r_f ≡ Δ0/(−q_h)`，檢係否存在 `−q_h` 之通用標度（正則 `q_h` 取重夸克電荷）。
  S6 不引入擬合；本輪只作符號與標度判定。
輸出：`data/48_sign_structure.json`
"""
import json
from pathlib import Path

M = {}
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
            M[int(k)] = mv * 1000.0

print("=" * 100)
print("M5 docs/46 — Δ0 符號結構：兩體 Coulomb vs 數據")
print("=" * 100)

# (族, 中性 id, 電荷 id, 中性之 ⟨q_h q_l⟩, 電荷之 ⟨q_h q_l⟩, q_h)
FAM = [
    ("π", 111, 211, (-4 / 9 - 1 / 9) / 2, (+2 / 3) * (+1 / 3), 2 / 3),
    ("K", 311, 321, (-1 / 3) * (+1 / 3), (-1 / 3) * (-2 / 3), -1 / 3),
    ("D", 421, 411, (+2 / 3) * (-2 / 3), (+2 / 3) * (+1 / 3), 2 / 3),
    ("B", 511, 521, (-1 / 3) * (+1 / 3), (-1 / 3) * (-2 / 3), -1 / 3),
]
print("\n[S1/S2/S3] 電荷乘積 → 模型符號；數據符號")
rows = []
for nm, i0, i1, pN, pC, qh in FAM:
    d0 = M[i0] - M[i1]
    s_model = -1 if pN < pC else (+1 if pN > pC else 0)   # 越吸引越低 ⟹ Δ0 = E_N − E_C 之符號
    s_data = 1 if d0 > 0 else -1
    agree = (s_model == s_data)
    r = d0 / (-qh)
    rows.append({"fam": nm, "d0": d0, "pN": pN, "pC": pC,
                 "s_model": s_model, "s_data": s_data, "agree": bool(agree), "r_reduced": r})
    print(f"  {nm:2s} ⟨q_hq_l⟩: N {pN:+.4f} vs C {pC:+.4f} ⟹ 模型 sign {s_model:+d}；"
          f"數據 Δ0 = {d0:+8.4f} ⟹ sign {s_data:+d} ⟹ {'✅ 同' if agree else '🔴 反'}"
          f"；r ≡ Δ0/(−q_h) = {r:+8.4f}")

print("\n[S4] 判")
bad = [r["fam"] for r in rows if not r["agree"]]
print(f"  符號相反之族：{bad if bad else '無'}")
print(f"  ⟹ 兩體 Coulomb {'連符號都錯（徹底出局）' if bad else '符號上可用'}")

print("\n[S5] 歸一化分裂 r_f ≡ Δ0/(−q_h)（檢 −q_h 通用標度）")
import statistics
for r in rows:
    print(f"  {r['fam']:2s} r = {r['r_reduced']:+8.4f} MeV")
pi_D = [r["r_reduced"] for r in rows if r["fam"] in ("π", "D")]
print(f"  π 與 D 之 r：{pi_D[0]:.4f} vs {pi_D[1]:.4f} ⟹ 差 "
      f"{100*abs(pi_D[0]/pi_D[1]-1):.2f}% ⟹ {'一致' if abs(pi_D[0]/pi_D[1]-1) < 0.1 else '不一致'}")
K = [r["r_reduced"] for r in rows if r["fam"] == "K"][0]
Bv = [r["r_reduced"] for r in rows if r["fam"] == "B"][0]
print(f"  K 偏離 π/D 之均值 {100*(K/statistics.mean(pi_D)-1):+.1f}%；B 偏離 {100*(Bv/statistics.mean(pi_D)-1):+.1f}%")

out = {"rows": rows, "sign_fail": bad,
       "r_pi_D_mean": statistics.mean(pi_D), "r_K": K, "r_B": Bv}
with open("data/48_sign_structure.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/48_sign_structure.json")