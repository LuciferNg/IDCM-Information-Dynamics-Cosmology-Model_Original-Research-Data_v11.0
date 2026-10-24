"""M4 docs/20 — 可辨識性曲線：IDCM 質量／em 預言在現時精度下「不可檢定」嘅定量界（2026-10-24）
==============================================================================================
背景（本 O 十輪經驗）：docs/09–18 全部死於「密度殺」（34.73%／27.90%／5.59%／4.39%／1.97%…）。
本輪把呢個現象升級為**可計算嘅界**：對每個 candidate 族，計 coverage(δ) = 隨機目標落 δ 帶嘅機率，
反解「要使單一命中具鑑別力（coverage ≤ 1e-3）需要幾精確」，再對比現時數據精度。

族（全部事前定義，唔准見結果再改）：
  F1 = {γ_a/γ_b : 1 ≤ a,b ≤ 44}                       （單模比，1936 值）
  F2 = {p(γ₁) : deg ≤ 2, |c_i| ≤ 10}                   （O140 形式框架，低次）
  F3 = {p(γ₁) : deg ≤ 4, |c_i| ≤ 10}                   （O140 形式框架，docs/19 用過）
  F4 = {n·γ₁^k : n ≤ 100, 1 ≤ k ≤ 12}                  （指數階梯族，docs/15 用過）

判準（事前）：
  coverage(δ) 用 2×10⁵ 個 log-uniform 目標（域 1e-2 … 1e4）估計；
  要求 coverage(δ*) ≤ 1e-3 ⟹ δ* = 「單一命中可辨識」所需相對精度；
  對比：quark ratio 精度（FLAG ~0.1–1%）、X = m_d−m_u 精度 3.7%、E_ch 精度 4.6–7.4%
  若 δ* ≪ 現時精度 ⟹ 「現階段不可檢定」係結構事實，唔係努力不足
"""
import json
import math
import numpy as np

OUT = "data/20_identifiability.json"
L = 88
g1 = 4 * math.sin(math.pi / L) ** 2
NMAX = 44
gam = np.array([4 * math.sin(n * math.pi / L) ** 2 for n in range(1, NMAX + 1)])

print("=" * 84)
print("M4 docs/20 — 可辨識性曲線（族密度 → 所需精度）")
print("=" * 84)
print(f"\nγ₁ = {g1:.12e}")

# ---------- 族 ----------
fam = {}
fam["F1 γ_a/γ_b"] = np.array([gam[i] / gam[j] for i in range(NMAX) for j in range(NMAX) if i != j])
def polys(deg, C):
    arr = np.array([0.0])
    for k in range(deg + 1):
        arr = np.add.outer(arr, np.array([c * g1 ** k for c in range(-C, C + 1)])).ravel()
    return np.unique(arr[arr > 0])
fam["F2 p(γ₁) deg≤2,|c|≤10"] = polys(2, 10)
fam["F3 p(γ₁) deg≤4,|c|≤10"] = polys(4, 10)
fam["F4 n·γ₁^k (n≤100,k≤12)"] = np.array([n * g1 ** k for n in range(1, 101) for k in range(1, 13)])

rng = np.random.default_rng(20261034)
NT = 200_000
tgt = 10 ** rng.uniform(-2, 4, NT)

def coverage(values, delta):
    v = np.sort(values)
    idx = np.clip(np.searchsorted(v, tgt), 1, v.size - 1)
    dev = np.abs(v[idx] / tgt - 1)
    dev = np.minimum(dev, np.abs(v[idx - 1] / tgt - 1))
    return float((dev < delta).mean())

def required_delta(values, target=1e-3):
    """coverage 隨 δ 單調遞增 ⟹ 二分：cov(mid) > target ⟹ hi = mid（收窄 δ）"""
    lo, hi = 1e-12, 1e-1
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if coverage(values, mid) > target:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)

print(f"\n{'族':<24}{'大小':>10}{'cov(0.1%)':>12}{'cov(1%)':>10}{'δ*(cov≤1e-3)':>16}{'δ*(cov≤1e-4)':>16}")
rows = []
for name, vals in fam.items():
    c01 = coverage(vals, 1e-3)
    c1 = coverage(vals, 1e-2)
    d3 = required_delta(vals, 1e-3)
    d4 = required_delta(vals, 1e-4)
    rows.append({"family": name, "size": int(len(vals)), "cov_0p1pct": c01, "cov_1pct": c1,
                 "delta_for_1e-3": d3, "delta_for_1e-4": d4})
    print(f"{name:<24}{len(vals):>10}{c01*100:>11.2f}%{c1*100:>9.2f}%{d3*100:>15.4f}%{d4*100:>15.5f}%")

# ---------- 對比現時精度 ----------
print("\n現時數據精度（相對）：")
prec = {"quark mass ratio (FLAG)": 0.001, "X = m_d−m_u (FLAG)": 0.037,
        "E_ch 最差 (docs/18)": 0.074, "em ratio Δπ (PDG)": 3e-5}
for k, v in prec.items():
    print(f"  {k:<26} {v*100:>7.3f}%")

loose = max(rows, key=lambda r: r["delta_for_1e-3"])   # δ* 最大 = 最寬鬆
strict = min(rows, key=lambda r: r["delta_for_1e-3"])
print(f"\n最寬鬆族：{loose['family']} ⟹ δ* = {loose['delta_for_1e-3']*100:.4f}%"
      f"（最嚴：{strict['family']} ⟹ {strict['delta_for_1e-3']*100:.5f}%）")
print("\n所需精度 vs 現有精度（短欠倍數 = 現有精度 / δ*）：")
shortfall = {}
for r in rows:
    line = f"  {r['family']:<24}"
    for k, v in prec.items():
        f_ = v / r["delta_for_1e-3"]
        line += f"  {k.split(' (')[0][:14]}: {f_:>9.0f}×"
    shortfall[r["family"]] = {k: prec[k] / r["delta_for_1e-3"] for k in prec}
    print(line)
print("\n判讀：短欠倍數 < 1 ⟹ 現有精度**足夠**；> 1 ⟹ 精度不足（需要更精確量度）")
for r in rows:
    ok = [k for k, v in shortfall[r["family"]].items() if v < 1]
    print(f"  {r['family']:<24} 精度足夠嘅數據源：{ok if ok else '無'}")

res = {"meta": {"script": "20_identifiability", "date": "2026-10-24", "seed": 20261034,
                "targets": NT, "domain_log10": [-2, 4]},
       "families": rows, "data_precision": prec, "loosest_family": loose, "strictest_family": strict,
       "shortfall": shortfall}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")