#!/usr/bin/python3
"""
M20 — scripts/02_single_chord_exhaustive.py
===========================================
R2／R3：**單鍵 ±1 缺陷類之窮盡掃描**（全部位置 × 兩種耦合約定 × ±1）。

目的（三問）：
  A. 符號定律：AF 鍵 vs FM 鍵對 ρ_d 之方向效應（結構性，非統計）
  B. 目標是否可由「單位強度單鍵」到達（若否則為**類內 no-go**）
  C. look-elsewhere：N_trials = 全部成員數，MC null 定 p

精確性：
  · 全部 G、δ、ρ 以 Fraction 精確（Woodbury 秩 1／秩 2，G₀ 閉式）
  · 只計 4×4 相關塊（k,m ∈ {2,3,4}）—— 因 δ 只涉此三站
  · 浮點只用於 MC null（明示 Monte-Carlo）與顯示
"""

import sys
import numpy as np
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 30
N, H = 44, 45
KU, KD, KS = 2, 3, 4

T_D, S_D = F(40350365, 10**8), F(351, 10000)
T_U, S_U = F(-2588837, 10**8), F(3292, 100000)

def g0(k, m):
    return F(min(k, m) * (H - max(k, m)), H)

SITES = (2, 3, 4)                                   # 讀出式所涉站

def block_of(Gf):
    return {(k, m): Gf(k, m) for k in SITES for m in SITES}

def ratios_from_block(B):
    d_d = B[(3, 3)] - B[(2, 3)]
    d_s = B[(3, 4)] - B[(2, 4)]
    d_u = B[(3, 2)] - B[(2, 2)]
    if d_s == 0:
        return None, None, (d_u, d_d, d_s)
    return d_d / d_s, d_u / d_s, (d_u, d_d, d_s)

# ---------------------------------------------------------------- 兩個約定之候選
def eval_L(i, j, lam):
    """Φ = λ·L_(i,j) = λ·u uᵀ，u = e_i − e_j，λ=±1。
    G = G₀ − (G₀u)(G₀u)ᵀ/(λ + uᵀG₀u)。"""
    u = {k: g0(k, i) - g0(k, j) for k in range(1, N + 1)}
    uTu = g0(i, i) - 2 * g0(i, j) + g0(j, j)
    den = lam + uTu                      # C⁻¹ = 1/λ = λ（λ = ±1）
    if den == 0:
        return None
    def Gf(k, m):
        return g0(k, m) - u[k] * u[m] / den
    return Gf, ("L", lam, i, j, den)

def eval_M(i, j, lam):
    """Φ = λ·(e_i e_jᵀ + e_j e_iᵀ)。U=[e_i,e_j]，C=[[0,λ],[λ,0]]，C⁻¹=C。
    M2 = C⁻¹ + UᵀG₀U = [[G₀(i,i), λ+G₀(i,j)],[λ+G₀(i,j), G₀(j,j)]]；G = G₀ − (G₀U) M2⁻¹ (G₀U)ᵀ。"""
    a, b, c = g0(i, i), g0(i, j), g0(j, j)
    M2 = [[a, lam + b], [lam + b, c]]
    det = M2[0][0] * M2[1][1] - M2[0][1] * M2[1][0]
    if det == 0:
        return None
    Minv = [[M2[1][1] / det, -M2[0][1] / det], [-M2[1][0] / det, M2[0][0] / det]]
    ci = [g0(k, i) for k in range(1, N + 1)]
    cj = [g0(k, j) for k in range(1, N + 1)]
    def Gf(k, m):
        p = [ci[k - 1], cj[k - 1]]
        q = [ci[m - 1], cj[m - 1]]
        corr = sum(p[r] * Minv[r][s] * q[s] for r in range(2) for s in range(2))
        return g0(k, m) - corr
    return Gf, ("M", lam, i, j, det)

print("=" * 92)
print("M20 / scripts/02 — 單鍵 ±1 缺陷類：窮盡掃描（全部位置 × {L,M} × {±1}）")
print("=" * 92)

base_d, base_u, _ = ratios_from_block(block_of(g0))
print(f"\n基線 ρ_d = {base_d} = {float(base_d):.8f} ｜ ρ_u = {base_u} = {float(base_u):.8f}")

rows = []
skipped = 0
for i in range(1, N + 1):
    for j in range(i + 2, N + 1):
        for lam in (+1, -1):
            for conv in ("L", "M"):
                r = eval_L(i, j, lam) if conv == "L" else eval_M(i, j, lam)
                if r is None:
                    skipped += 1
                    continue
                Gf, meta = r
                rd, ru, dl = ratios_from_block(block_of(Gf))
                if rd is None:
                    skipped += 1
                    continue
                rows.append((conv, lam, i, j, rd, ru, dl, meta[4]))

print(f"成員數 = {len(rows)}（奇異/退化略去 {skipped}；全對 (i,j), |i−j|≥2 ）")

rd_f = np.array([float(r[4]) for r in rows])
ru_f = np.array([float(r[5]) for r in rows])
td, tu = float(T_D), float(T_U)
sd, su = float(S_D), float(S_U)
zd = (rd_f - td) / sd
zu = (ru_f - tu) / su
chi = np.sqrt(zd ** 2 + zu ** 2)

# ---------------------------------------------------------------- A. 符號定律
print("\n" + "-" * 92)
print("A. 符號定律（結構性）：耦合符號 ⟹ ρ_d 之方向")
print("-" * 92)
for conv in ("L", "M"):
    for lam in (+1, -1):
        m = np.array([r[0] == conv and r[1] == lam for r in rows])
        if m.sum() == 0:
            continue
        tag = ("AF" if lam > 0 else "FM") if conv == "L" else ("FM" if lam > 0 else "AF")
        print(f"  {conv}-form λ={lam:+d} （{tag}）: n={m.sum():4d}  "
              f"ρ_d ∈ [{rd_f[m].min():+.6f}, {rd_f[m].max():+.6f}]  "
              f">42/41 之比例 = {np.mean(rd_f[m] > float(base_d))*100:5.1f}%")
print("  （L-form: 非對角 = −λ ⟹ λ=+1 = AF、λ=−1 = FM；M-form: 非對角 = +λ ⟹ 相反）")

# ---------------------------------------------------------------- B. 目標可達性
print("\n" + "-" * 92)
print("B. 目標可達性（㈡：1/ρ_d ≈ 2.4783 ⟺ ρ_d ≈ 0.40350）")
print("-" * 92)
best_d = int(np.argmin(np.abs(zd)))
print(f"  min |z_d| = {abs(zd[best_d]):.4f} 於 {rows[best_d][:4]}  ρ_d = {float(rows[best_d][4]):.8f}")
for thr in (1.0, 2.0, 3.0):
    print(f"  |z_d| < {thr}: {int((np.abs(zd) < thr).sum())} / {len(rows)}"
          f"   （同 |z_u| < 2 且 |z_d| < {thr}: "
          f"{int(((np.abs(zd) < thr) & (np.abs(zu) < 2)).sum())}）")
best_c = int(np.argmin(chi))
print(f"  min χ = {chi[best_c]:.4f} 於 {rows[best_c][:4]}  "
      f"ρ_d={float(rows[best_c][4]):.8f} (z={zd[best_c]:+.3f})  "
      f"ρ_u={float(rows[best_c][5]):.8f} (z={zu[best_c]:+.3f})")

# ---------------------------------------------------------------- C. MC null
print("\n" + "-" * 92)
print("C. look-elsewhere MC null（10^4 次；零假設 = 目標真值落於其誤差雲內）")
print("-" * 92)
rng = np.random.default_rng(20261119)
R = 10000
mins = np.empty(R)
for r in range(R):
    t1 = rng.normal(td, sd)
    t2 = rng.normal(tu, su)
    mins[r] = np.min(np.sqrt(((rd_f - t1) / sd) ** 2 + ((ru_f - t2) / su) ** 2))
p = float(np.mean(mins <= chi[best_c]))
print(f"  N_trials = {len(rows)}；min χ_obs = {chi[best_c]:.4f}")
print(f"  null min χ：中位 {np.median(mins):.4f}｜5% 分位 {np.percentile(mins,5):.4f}"
      f"｜min {mins.min():.4f}")
print(f"  p = P(min χ_null ≤ min χ_obs) = {p:.4f}   ⟹ {'信號' if p < 0.05 else '零命中（look-elsewhere 吸收）'}")

# ---------------------------------------------------------------- 精確值（前 3 名 + 錨點）
print("\n" + "-" * 92)
print("D. 精確值（Fraction）：min χ 前 3 名 ＋ 味錨定候選 (2,4)")
print("-" * 92)
order = np.argsort(chi)[:3]
for idx in order:
    conv, lam, i, j, rd, ru, dl, det = rows[idx]
    print(f"  {conv} λ={lam:+d} ({i},{j})  δ_u={dl[0]} δ_d={dl[1]} δ_s={dl[2]}")
    print(f"      ρ_d = {rd} = {float(rd):.8f}  z={zd[idx]:+.3f}   1/ρ_d = {1/rd}")
    print(f"      ρ_u = {ru} = {float(ru):.8f}  z={zu[idx]:+.3f}")
for r in rows:
    if r[0] == "L" and r[2] == 2 and r[3] == 4:
        print(f"  [錨定] L λ={r[1]:+d} (2,4)  ρ_d = {r[4]}  ρ_u = {r[5]}")
    if r[0] == "M" and r[2] == 2 and r[3] == 4:
        print(f"  [錨定] M λ={r[1]:+d} (2,4)  ρ_d = {r[4]}  ρ_u = {r[5]}")

print("\n" + "=" * 92)
print("結論欄（自動生成）")
print("=" * 92)
n_hit = int(((np.abs(zd) < 1) & (np.abs(zu) < 1)).sum())
print(f"  · 同時 |z_d|<1 且 |z_u|<1 之成員數 = {n_hit} / {len(rows)}")
verdict = "有信號" if p < 0.05 else "無信號（N_trials 密度吸收）"
print(f"  · min χ 之 p = {p:.4f} ⟹ {verdict}")

# ---------------------------------------------------------------- E. 分母結構
print("\n" + "-" * 92)
print("E. 分母結構：ρ ∈ ℚ/? （類級精度界限）")
print("-" * 92)
dens_d = {}
dens_u = {}
for r in rows:
    dens_d[r[4].denominator] = dens_d.get(r[4].denominator, 0) + 1
    dens_u[r[5].denominator] = dens_u.get(r[5].denominator, 0) + 1
print(f"  ρ_d 之分母分佈 = {dict(sorted(dens_d.items()))}")
print(f"  ρ_u 之分母分佈 = {dict(sorted(dens_u.items()))}")
all41 = all(d in (1, 41) for d in dens_d) and all(d in (1, 41) for d in dens_u)
print(f"  ⟹ 全類 ρ ∈ ℚ/41（分母 ∈ {{1,41}}）：{'✅ 是' if all41 else '🔴 否'}")
if all41:
    print("     結構原因：δ_s = G(3,4) − G(2,4)，站 4 之 Dirichlet 距離 45−4 = 41")
    print("     ⟹ M12 R25 之「格點吻合機率 100%」由單點升級為**類級**：本類任何成員皆不可判別")

# ---------------------------------------------------------------- F. 條件可達集
print("\n" + "-" * 92)
print("F. 條件可達集：ρ_u ≈ 目標（|z_u|<1）之成員，其 ρ_d 為何？")
print("-" * 92)
sel = np.abs(zu) < 1
print(f"  |z_u|<1 之成員數 = {int(sel.sum())} / {len(rows)}")
vals = sorted({rows[k][4] for k in range(len(rows)) if sel[k]})
print(f"  其 ρ_d 取值集（精確，共 {len(vals)} 個）= {vals}")
print(f"  其中 |z_d|<3 者 = {[str(v) for v in vals if abs(float(v)-td)/sd < 3]}")
print("=" * 92)
sys.exit(0)
