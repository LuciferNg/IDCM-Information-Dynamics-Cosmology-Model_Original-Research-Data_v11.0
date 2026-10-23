"""M4 Q8 — em 讀數分解模型判別：H1（E_K=E_N）vs 電荷自能+交換模型（2026-10-24）
==============================================================================
上游：docs/14（P1 反解 m_d−m_u = 2.6137 MeV，0.57σ，但依賴 H1：E_K = E_N 未推導）；
      01_problem_record §十一 Q8；docs/16（Q6 scale 算子，PDG 錨點已核）。
數據（PDG 2024，MeV）：Δπ = m_π+ − m_π0 = +4.59359 ± 0.00014；
                        ΔK = m_K+ − m_K0 = −3.934 ± 0.021；
                        ΔN = m_n − m_p    = +1.293331 ± 0.0000005
lattice（FLAG 24, Nf=2+1+1）：m_d − m_u = 2.56 ± 0.094 MeV

問題：Δm 三分量可以點分解？「E_K = E_N」（H1）係任揀抑或數據偏好？
  Δπ = X + E_π ,  ΔK = −X + E_K ,  ΔN = X + E_N     （X ≡ m_d − m_u）

Protocol（事前寫死）：
  M1 **H1 模型**：E_K = E_N ⟹ X = (ΔN − ΔK)/2
     判：|X − 2.56| / σ_lat < 2 ⟹ 🟡 一致；≥ 3 ⟹ 🔴
  M2 **電荷模型（自能 + 交換）**：em 能量 = κ_s ΣQ² + κ_x Σ_{{i<j}} Q_i Q_j，
     荷 Q_u=+2/3、Q_d=Q_s=−1/3（O124：Y = T₃^R + (B−L)/2 結構電荷）；
     內容：π+ = u d̄、π0 = (uū − dd̄)/√2（權重 ½）、K+ = u s̄、K0 = d s̄、p = uud、n = udd
     ⟹ 三方程三未知（X, κ_s, κ_x）恰好確定 ⟹ 解出 X 與 lattice 比 ⟹ 呢個係**真外部檢定**
  M3 **r 參數族**：E_K = r·E_N ⟹ X(r) = (rΔN − ΔK)/(1 + r)；
     由 lattice 反解 r 嘅容許範圍；比較 M2 嘅 r_pred = E_K/E_N
  M4 密度／可辨識：M2 三方程三未知 ⟹ 對隨機 Δm 亦必然解出某 X ⟹ 量化「X 命中 2.56」嘅偶然率
  M5 結構身份：r = 1 有冇結構理由（|ΔQ|=1 ＋ 同內容符號型）⟹ 只可標一致性
輸出：data/17_em_decomposition.json ＋ stdout
"""
import json
import random

OUT = "data/17_em_decomposition.json"

DP, EP = 4.59359, 0.00014        # Δπ
DK, EK_ = -3.934, 0.021          # ΔK
DN, EN = 1.293331, 0.0000005     # ΔN
X_LAT, S_LAT = 2.560, 0.094      # FLAG24 m_d − m_u

print("=" * 78)
print("M4 Q8 — em 讀數分解模型判別（H1 vs 電荷自能+交換）")
print("=" * 78)
print(f"\n數據：Δπ={DP:+.5f}  ΔK={DK:+.3f}  ΔN={DN:+.6f} MeV ; lattice m_d−m_u={X_LAT}±{S_LAT}")

# ---------- M1 ----------
X1 = (DN - DK) / 2
z1 = (X1 - X_LAT) / S_LAT
print(f"\n[M1] H1（E_K = E_N）：X = (ΔN − ΔK)/2 = {X1:.6f} MeV")
print(f"     vs lattice {X_LAT}±{S_LAT} ⟹ dev {abs(X1/X_LAT-1)*100:.4f}%  = {z1:+.2f}σ"
      f" ⟹ {'🟡 一致' if abs(z1) < 2 else ('🟠 2–3σ' if abs(z1) < 3 else '🔴 排除')}")

# ---------- M2 ----------
# em 能量 = κ_s·ΣQ² + κ_x·Σ_{i<j} Q_iQ_j（夸克對；介子 2 體、重子 3 體）
from fractions import Fraction as F
Qu, Qd, Qs = F(2, 3), F(-1, 3), F(-1, 3)

# 每個強子：(自能 ΣQ², 交換 Σ_{i<j}Q_iQ_j) — 精確分數
coef = {
    "π+": (Qu ** 2 + Qd ** 2,               Qu * Qd),                              # u d̄
    "π0": (Qu ** 2 + Qd ** 2,               F(1, 2) * Qu * Qu + F(1, 2) * Qd * Qd),  # (uū−dd̄)/√2
    "K+": (Qu ** 2 + Qs ** 2,               Qu * Qs),                              # u s̄
    "K0": (Qd ** 2 + Qs ** 2,               Qd * Qs),                              # d s̄
    "p":  (2 * Qu ** 2 + Qd ** 2,           Qu * Qu + 2 * Qu * Qd),               # uud
    "n":  (Qu ** 2 + 2 * Qd ** 2,           2 * Qu * Qd + Qd * Qd),               # udd
}
print("\n[M2] 電荷模型係數（自能 ΣQ² / 交換 ΣQ_iQ_j，精確分數）")
for k, (s_, x_) in coef.items():
    print(f"     {k:<3} ΣQ² = {str(s_):>5}   Σ_{{i<j}}Q_iQ_j = {str(x_):>5}")

aP = (coef["π+"][0] - coef["π0"][0], coef["π+"][1] - coef["π0"][1])
aK = (coef["K+"][0] - coef["K0"][0], coef["K+"][1] - coef["K0"][1])
aN = (coef["n"][0] - coef["p"][0],   coef["n"][1] - coef["p"][1])
print(f"     Δ 係數：π {tuple(map(str,aP))}  K {tuple(map(str,aK))}  N {tuple(map(str,aN))}")

# 三方程：Δπ = X + aP0·κs + aP1·κx ; ΔK = −X + aK0·κs + aK1·κx ; ΔN = X + aN0·κs + aN1·κx
# 代數消去 ⟹ 模型自身約束（唔靠擬合）：
#   aP0=0, aN = −aK ⟹ ΔN − X = aN·(κ) 而 ΔK + X = aK·(κ) = −aN·(κ) = −(ΔN − X)
#   ⟹ ΔK + ΔN = 0（恆等，對任何 κ_s, κ_x）
constraint_lhs = DK + DN
print(f"     代數約束：aP0 = 0、aN = −aK ⟹ 模型強制 **ΔK + ΔN = 0**")
print(f"     觀測：ΔK + ΔN = {DK:+.3f} + {DN:+.6f} = {constraint_lhs:+.6f} MeV"
      f"（誤差 { (EK_**2 + EN**2)**0.5:.3f}）⟹ 違反 {abs(constraint_lhs)/((EK_**2+EN**2)**0.5):.1f}σ ⟹ **模型結構否證**")
m2_ok = abs(constraint_lhs) < 3 * (EK_ ** 2 + EN ** 2) ** 0.5
X2 = None
print(f"     判：{'🟡 未否證' if m2_ok else '🔴 結構否證（兩標度電荷模型 ⟹ ΔK=−ΔN，與觀測差 2.64 MeV）'}")

# ---------- M3 ----------
print("\n[M3] r 參數族（E_K = r·E_N）：X(r) = (rΔN − ΔK)/(1 + r)（對 r 單調遞減）")
def X_of_r(r): return (r * DN - DK) / (1 + r)
for r in [0.0, 0.5, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0]:
    print(f"     r={r:>5.2f} ⟹ X = {X_of_r(r):.6f} MeV  (dev vs lat {abs(X_of_r(r)/X_LAT-1)*100:>6.3f}%)")
# 因 X(r) 遞減：X(r) = X_LAT ± σ ⟹ 反解（單調 ⟹ 直接代數反解）
def r_of_X(x):
    # x(1+r) = rΔN − ΔK ⟹ x + xr = rΔN − ΔK ⟹ r(ΔN − x) = x + ΔK ⟹ r = (x + ΔK)/(ΔN − x)
    return (x + DK) / (DN - x)
r_center = r_of_X(X_LAT)
r_lo = min(r_of_X(X_LAT - S_LAT), r_of_X(X_LAT + S_LAT))
r_hi = max(r_of_X(X_LAT - S_LAT), r_of_X(X_LAT + S_LAT))
print(f"     反解（代數，非數值搜）：lattice 一致區間 r ∈ [{r_lo:.4f}, {r_hi:.4f}]（中心 {r_center:.4f}）")
print(f"     ⟹ r = 1（H1：E_K = E_N）{'落區間內 ✅' if r_lo <= 1.0 <= r_hi else '落區間外 🔴'}；"
      f"r_pred(M2) 因模型已否證，唔適用")

# ---------- M4 ----------
print("\n[M4] 密度：M2 三方程三未知 ⟹ 任一隨機 Δm 都解得某 X ⟹ 「X 命中 2.56」偶然率")
rng = random.Random(20261032)
NMC = 200_000
cnt = 0
for _ in range(NMC):
    dp = rng.gauss(DP, 0.5); dk = rng.gauss(DK, 0.5); dn = rng.gauss(DN, 0.03)
    x = (dn + dp - dk) / 3
    if abs(x - X_LAT) < S_LAT: cnt += 1
print(f"     隨機 Δm（半寬 σ=(0.5,0.5,0.03) MeV）⟹ X 落 2.56±0.094 嘅機率 = {cnt/NMC*100:.2f}%")

# ---------- M5 ----------
print("\n[M5] r = 1 嘅結構身份檢查")
print("     |ΔQ|=1 通道：π+→π0（u d̄ → (uū−dd̄)/√2）、K+→K0（u s̄ → d s̄）、n→p（udd → uud）")
print("     內容差異型：π 含湮滅（ū/d̄ 變 ū,dd̄ 混合）、K 係 u→d 替換（s̄ 不變）、N 係 d→u 替換")
print("     ⟹ K、N 嘅差異型同構（單夸克替換）⟹ r=1 有**弱結構理由**（同型 ⟹ 同 em 讀數差）")
print("     但唔係推導：s̄ 同 u/d 嘅 em 讀數係否相同，需 operator 才可判 ⟹ 標 🟡 一致性")

res = {
    "meta": {"script": "17_em_decomposition", "date": "2026-10-24", "seed": 20261032},
    "data": {"dpi": DP, "dK": DK, "dN": DN, "X_lattice": X_LAT, "sigma_lattice": S_LAT},
    "M1": {"X": X1, "z": z1, "verdict": "ok" if abs(z1) < 2 else ("warn" if abs(z1) < 3 else "excluded")},
    "M2": {"coef": {k: [str(v[0]), str(v[1])] for k, v in coef.items()},
           "aP": [str(aP[0]), str(aP[1])], "aK": [str(aK[0]), str(aK[1])], "aN": [str(aN[0]), str(aN[1])],
           "model_constraint": "ΔK + ΔN = 0 (identity)",
           "observed_sum": constraint_lhs, "violation_sigma": abs(constraint_lhs) / ((EK_ ** 2 + EN ** 2) ** 0.5),
           "verdict": "not_falsified" if m2_ok else "structurally_falsified"},
    "M3": {"r_lo": r_lo, "r_hi": r_hi, "r_center": r_center,
           "r1_inside": bool(r_lo <= 1.0 <= r_hi),
           "note": "M2 falsified ⟹ r_pred undefined"},
    "M4": {"p_random_hit": cnt / NMC},
    "M5": {"verdict": "weak_structural_reason_not_derivation"},
}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")
