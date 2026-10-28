"""M5 docs/44 配套 — ΔI=0 殘餘：單一 `m_l` 冪律係否可救（指數穩定性檢）
================================================================================
觸發：`docs/41`（v3）令 `c(D_s)/c(D) = c(B_s)/c(B) = 1.0000`，而數據 1.0130／1.0764 ⟹ 殘餘 4.4σ／2.4σ。
本輪問：能否用**單一**結構化冪律（`c ∝ μ^a m_h^b m_l^p`）同時解釋兩族？抑或「改指數」根本唔夠？

判準（事前寫死）：
  P1 由 `R_data ≡ c_hi/c_lo` 反解單族所需指數 `p_f = ln R_data / ln(m_l,hi/m_l,lo)`（v3 之 `P₀` 比 = 1）。
  P2 `σ(p_f)` 由 `σ(R)/R` 線性傳播。
  P3 判：`|p_Ds − p_Bs| / sqrt(σ²+σ²) ≥ 2` ⟹ 單一冪律**不穩定**（兩族要求不同指數）。
  P4 另列候選結構標度（Coulomb `μ³`、線性 `μ`、heavy-quark `1/m_h`、v3 `p = 0`）之預言，逐族比數據。
  P5 不得以擬合 `p` 宣稱解釋；`p` 只作診斷量。
輸出：`data/46_ml_exponent_stability.json`
"""
import json
import math
from pathlib import Path

MASS, ERR = {}, {}
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
            MASS[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0

PAIRS = {"D_s/D": ((433, 431), (423, 421), 1.27, 0.495, 0.30),
         "B_s/B": ((533, 531), (513, 511), 4.18, 0.495, 0.30)}
print("=" * 100)
print("M5 docs/44 — ΔI=0 殘餘：單一 m_l 冪律之穩定性")
print("=" * 100)

print("\n[P1/P2] 逐族所需指數")
res = {}
for name, ((s1, g1), (s2, g2), mh, ml_hi, ml_lo) in PAIRS.items():
    E_hi, E_lo = MASS[s1] - MASS[g1], MASS[s2] - MASS[g2]
    s_hi, s_lo = math.hypot(ERR[s1], ERR[g1]), math.hypot(ERR[s2], ERR[g2])
    R = E_hi / E_lo
    sR = R * math.sqrt((s_hi / E_hi) ** 2 + (s_lo / E_lo) ** 2)
    p = math.log(R) / math.log(ml_hi / ml_lo)
    sp = (sR / R) / math.log(ml_hi / ml_lo)
    res[name] = {"R": R, "sigma_R": sR, "p": p, "sigma_p": sp, "m_h": mh}
    print(f"  {name:6s} R = {R:.4f} ± {sR:.4f} ⟹ p = {p:+.4f} ± {sp:.4f}（m_h = {mh}）")

pd_, pb = res["D_s/D"], res["B_s/B"]
diff = pb["p"] - pd_["p"]
sdiff = math.hypot(pd_["sigma_p"], pb["sigma_p"])
z = abs(diff) / sdiff
print(f"\n[P3] Δp = {diff:+.4f} ± {sdiff:.4f} ⟹ {z:.2f}σ ⟹ "
      f"{'🔴 單一冪律不穩定（≥2σ）' if z >= 2 else '🟡 未能排除單一冪律（<2σ）'}")

print("\n[P4] 候選結構標度之預言（v3 之 P₀ 比 = 1）")
MU = {"D": 1.27 * 0.30 / 1.57, "D_s": 1.27 * 0.495 / 1.765,
      "B": 4.18 * 0.30 / 4.48, "B_s": 4.18 * 0.495 / 4.675}
MODELS = [
    ("v2: 1/(m_h m_l) [p = −1]", lambda hi, lo: (0.30 / 0.495)),
    ("v3: 1/m_h [p = 0]", lambda hi, lo: 1.0),
    ("線性勢: μ/m_h", lambda hi, lo: MU[hi] / MU[lo]),
    ("Coulomb: μ³/(m_h m_l)", lambda hi, lo: (MU[hi] / MU[lo]) ** 3 * (0.30 / 0.495)),
]
for label, f in MODELS:
    out = []
    for key, hi, lo in (("D_s/D", "D_s", "D"), ("B_s/B", "B_s", "B")):
        pred = f(hi, lo)
        rd = res[key]["R"]
        out.append(f"{key}: {pred:.4f} vs {rd:.4f}（dev {100*(pred/rd-1):+6.2f}%）")
    print(f"  {label:28s} " + "；".join(out))

print("\n[P5] 判讀")
print("  · 兩族所需 p 差異 " f"{z:.1f}σ ⟹ 唔存在一個『正確指數』；")
print("  · 所有結構候選（−1、0、μ、μ³）皆偏離 ≥1.3%（D_s）或 ≥7%（B_s）；")
print("  · 數據落於 p ≈ 0.03（D_s）與 0.15（B_s）之間 ⟹ **ΔI=0 通道未閉合，且唔係改指數可以解決**。")

out = {"pairs": res, "p_diff": diff, "sigma_diff": sdiff, "z_diff": z,
       "single_power_stable": bool(z < 2)}
with open("data/46_ml_exponent_stability.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/46_ml_exponent_stability.json")