"""M5 docs/52 配套 — A_eff 與 B 抑制因子之閉式掃描（預註冊語法 ＋ MC 零假設）
================================================================================
兩項未導出量：
  ② `A_eff = 4.4492 ± 0.4603 MeV` ⟹ 無量綱目標 `T₂ ≡ A_eff/(α ħc/a₀) = 1.18707`
  ③ B 之抑制 `s = 1/14.35` ⟹ 目標 `T₃ ≡ 1/s = 14.352`

判準（**先預註冊，見下 `ATOMS`／語法；不得事後增刪**）：
  C1 原子集只取**已登記**結構常數（IDCM 核、群論階、簡單整數、φ、√n、π、e、ln2）。
  C2 語法：`v = a^i · b^j`，`a,b ∈ ATOMS`、`a ≠ b`、`i,j ∈ {−2,−1,1,2}`（`i = j = 0` 排除）。
  C3 命中 = `|v/T − 1| ≤ 1e−4`（相對）。
  C4 MC 零假設：`T` 在語法值域之對數區間內均勻抽 10⁷ 個 ⟹ 得命中數分佈；
     `p = P(N_hit ≥ N_obs)`。若 `p > 0.01` ⟹ **命中數與隨機無異，不得宣稱結構式**。
  C5 只報「命中清單 ＋ p 值」，不作任何「最接近」之挑選（禁事後凑）。
輸出：`data/54_closed_form_scan.json`
"""
import itertools
import json
import math
import random

ALPHA = 1 / 137.035999084
HBARC = 197.3269804
A0 = 0.384184

ATOMS = {
    "α^-1": 1 / ALPHA,
    "g(0)": 2581 / 352,
    "g(7)": 1447 / 352,
    "N_d": 88.0,
    "b": 4.0,
    "N_gen": 3.0,
    "Σγ": 90.0,
    "N_d+1": 45.0,
    "φ": (1 + 5 ** 0.5) / 2,
    "√5": 5 ** 0.5,
    "√2": 2 ** 0.5,
    "√3": 3 ** 0.5,
    "π": math.pi,
    "e": math.e,
    "ln2": math.log(2),
    "2": 2.0, "3": 3.0, "4": 4.0, "5": 5.0, "6": 6.0, "7": 7.0, "8": 8.0,
    "9": 9.0, "12": 12.0, "19": 19.0, "41": 41.0, "60": 60.0,
}
EXPS = (-2, -1, 1, 2)

print("=" * 100)
print("M5 docs/52 — 預註冊閉式掃描（A_eff 與 B 抑制）")
print("=" * 100)

# 建語法值表
vals = []
names = list(ATOMS)
for a, b in itertools.permutations(names, 2):
    for i in EXPS:
        for j in EXPS:
            v = ATOMS[a] ** i * ATOMS[b] ** j
            if v > 0 and math.isfinite(v):
                vals.append(v)
vals = sorted(vals)
vmin, vmax = vals[0], vals[-1]
loglo, loghi = math.log(vmin), math.log(vmax)
print(f"\n[C1/C2] 原子 {len(names)} 個；語法表達式 {len(vals)} 個；"
      f"值域 [{vmin:.4g}, {vmax:.4g}]（對數寬 {loghi-loglo:.3f}）")

T2 = 4.4492 / (ALPHA * HBARC / A0)
T3 = 1 / (0.3100 / 4.4492)
print(f"\n目標：T₂ = A_eff/(α ħc/a₀) = {T2:.6f}；T₃ = 1/s = {T3:.6f}")

TOL = 1e-4


def scan(T):
    hits = []
    lo, hi = T * (1 - TOL), T * (1 + TOL)
    import bisect
    i = bisect.bisect_left(vals, lo)
    while i < len(vals) and vals[i] <= hi:
        hits.append(vals[i])
        i += 1
    return hits


print("\n[C3] 掃描")
h2, h3 = scan(T2), scan(T3)
print(f"  T₂ 命中 {len(h2)} 個：{h2 if h2 else '無'}")
print(f"  T₃ 命中 {len(h3)} 個：{h3 if h3 else '無'}")

print("\n[C4] MC 零假設（10⁷ 個均勻對數目標）")
random.seed(20261028)
NS = 10_000_000
counts = {0: 0, 1: 0, 2: 0}
more = 0
import bisect
for _ in range(NS):
    T = math.exp(random.uniform(loglo, loghi))
    lo, hi = T * (1 - TOL), T * (1 + TOL)
    i = bisect.bisect_left(vals, lo)
    n = 0
    while i < len(vals) and vals[i] <= hi:
        n += 1
        i += 1
    if n >= 2:
        more += 1
    else:
        counts[n] += 1
p0 = counts[0] / NS
p1 = counts[1] / NS
p2 = more / NS
print(f"  P(N=0) = {p0:.6f}；P(N=1) = {p1:.6f}；P(N≥2) = {p2:.6f}")
print(f"  合理檢查：P(N=0) 應 ≈ 1 − {2*TOL:.2e}·{len(vals)}/({loghi-loglo:.2f})")


def pval(n_obs):
    if n_obs == 0:
        return p0
    if n_obs == 1:
        return p1 + p2
    return p2


for nm, hits, T in (("T₂", h2, T2), ("T₃", h3, T3)):
    n = len(hits)
    p = pval(n)
    print(f"\n  {nm}：N_obs = {n} ⟹ p = {p:.4f} ⟹ "
          f"{'🟡 命中數與隨機無異（不得宣稱結構式）' if p > 0.01 else '★ 超出隨機'}")

print("\n[C5] 紀律：即使命中亦只列值，不挑「最接近」；且 T₂ 之 1.187 與 T₃ 之 14.35 皆由數據定出 ⟹ 屬單點反解。")

print("\n[附] B 之抑制所需指數（三種歸一化）")
mC, mB = 1.27, 4.18
for lab, r in (("m_b/m_c", mB / mC), ("(r_h/a₀)_D/(r_h/a₀)_B", (HBARC / (mC * A0)) / (HBARC / (mB * A0)))):
    k = math.log(T3) / math.log(r)
    print(f"  ∝ ({lab})^−k ⟹ k = {k:.4f}（√5 = {5**0.5:.4f}，ln10 = {math.log(10):.4f}；"
          f"偏離 {100*(k/5**0.5-1):+.2f}% / {100*(k/math.log(10)-1):+.2f}%）")

out = {"n_expr": len(vals), "log_range": [loglo, loghi], "T2": T2, "T3": T3,
       "hits_T2": h2, "hits_T3": h3, "p0": p0, "p1": p1, "p2": p2,
       "p_T2": pval(len(h2)), "p_T3": pval(len(h3)), "N_samples": NS}
with open("data/54_closed_form_scan.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/54_closed_form_scan.json")