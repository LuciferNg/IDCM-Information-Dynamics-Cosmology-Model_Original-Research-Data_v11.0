"""M4 Q7 — `M_EW = M_Pl·γ₁^k·N_d` 嘅指數 k = 8 可否事前推導（唔准反選）（2026-10-24）
==============================================================================
上游：docs/10 §四（M_EW = M_Pl γ₁⁸ N_d = 244.2372173 GeV vs v_SM 246.2196 ⟹ dev 0.8051%；
          族 multiplicity thr 0.5%/0.8% → 0、1.0% → 1；指數 8 由 Λ²(12)=8 **反選** ⟹ 維持 🟡）；
      docs/13 §六 Q7；01_problem_record §十一 Q7。
問題：k = 8 有冇**事前**結構推導？抑或只係「見到 v 之後揀指數」？（反選 = 禁區）

Protocol（事前寫死）：
  G-family: F = { M_Pl · γ₁^k · S }，k ∈ [−45,45]（91），S ∈ 18 個結構常數（共 13104，同 docs/12 族）
  anchors (8)：1 eV, m_e, Λ_IDCM, M₀, v_SM, √(M_Pl·Λ), M_Pl, Λ_IDCM²/M₀（單位全 MeV）

  T1 跨標度指數一致性：對 8 個 anchor 各自求 best (k,S,dev)；問「best k 係咪全部 = 8？」
     判：若 ≥3 個 anchor 嘅 best k = 8 ⟹ 指數有跨實例支持；若只有 v ⟹ v-specific（🔴）
  T2 k=8 分支 multiplicity：固定 k=8，列 S 使 dev<1%；再算 density（隨機目標域 = [10,1000] GeV）
  T3 機制 M1「γ₁ = 每跨一層 K_chain 同步邊界」：k_mech = log(v/M_Pl)/log(γ₁)
     判：|k_mech − 8| < 0.05 ⟹ M1 支持；否則 M1 否證（γ₁ 非「每層」因子）
  T4 機制 M2「指數 = dim su(N_gen) = N_gen²−1 = 8」：
     判：若 k=8 只係 v 嘅解、其他 colour-dressed 標度（Λ, M₀）嘅 best k ≠ 8 ⟹ M2 唔被支持
  T5 identity grammar multiplicity：「8 有身份」係咪罕有？
     grammar（事前固定）：atoms A = {1,2,3,4,5,7,8,11,12,13,18,21,29,44,76,88}；
       ops: a+b, a−b(>0), a×b(≤200), a/b(整數), a², a(a−1)/2, 2a, a/2, a−1, a+1
     計每個 n ∈ [1,30] 嘅相異表達式數 → 8 嘅排名／百分位
     判：若 8 唔喺 top-3（即身份唔罕有）⟹ 「8 = Λ²(12)／N_gen²−1／L/L₅」屬**多樣性冗餘**，無資訊
  T6 leave-one-out：用其他 anchor 嘅 best-k 規律預測 v 嘅 k（線性／最近鄰），問係咪 8
輸出：data/15_exponent8.json ＋ stdout
"""
import json
import numpy as np
import random

OUT = "data/15_exponent8.json"

# ---------- 譜 ----------
L = 88
N = np.arange(1, L // 2 + 1)
GAM = 4 * np.sin(N * np.pi / L) ** 2
G1 = GAM[0]
LOGG1 = np.log(G1)

M_PL = 1.220910e22       # MeV
ANCH = {"1 eV": 1e-6, "m_e": 0.51099895, "Λ_IDCM": 45.28487, "M0": 27.900915,
        "v_SM": 246219.6, "sqrt(MPl·Λ)": np.sqrt(M_PL * 45.28487),
        "M_Pl": M_PL, "Λ²/M0": 45.28487 ** 2 / 27.900915}
S18 = {"1": 1, "B": 4, "N_d": 44, "L": 88, "N_gen": 3, "F5": 5, "F7": 13, "F8": 21,
       "L5": 11, "L7": 29, "L9": 76, "N_gen²": 9, "L6": 18, "B·N_gen": 12,
       "L+N_d": 132, "N_d/2": 22, "L/2": 44, "√5": np.sqrt(5)}
KS = np.arange(-45, 46)

print("=" * 78)
print("M4 Q7 — M_EW 指數 k=8 可否事前推導")
print("=" * 78)
print(f"\nγ₁ = {G1:.12f}   log(γ₁) = {LOGG1:.9f}   M_Pl = {M_PL:.6g} MeV")
print("族 F = { M_Pl·γ₁^k·S }，每 anchor |F| = %d（8 anchor 共 %d）" % (len(KS) * len(S18), len(KS) * len(S18) * len(ANCH)))

# ---------- T1 ----------
print("\n[T1] 跨標度指數一致性（8 anchor × 91 k × 18 S）")
print(f"  {'anchor':<14}{'best k':>7}{'S':>8}{'value':>14}{'dev%':>10}")
T1 = {}
for an, av in ANCH.items():
    best = None
    for k in KS:
        base = M_PL * G1 ** k
        for sn, sv in S18.items():
            v = base * sv
            dev = abs(v / av - 1) * 100
            if best is None or dev < best[4]:
                best = (an, int(k), sn, v, dev)
    a, k, sn, v, dev = best
    T1[an] = {"k": k, "S": sn, "value": v, "dev_pct": dev}
    print(f"  {an:<14}{k:>7}{sn:>8}{v:>14.6g}{dev:>10.4f}")
ks = [T1[a]["k"] for a in ANCH]
print(f"  best-k 分佈: {ks}")
n8 = sum(1 for k in ks if k == 8)
print(f"  best k = 8 嘅 anchor 數 = {n8}/8 ⟹ {'🟢 跨實例支持' if n8 >= 3 else '🔴 指數 8 係 v-specific'}")

# ---------- T2 ----------
print("\n[T2] k=8 分支 multiplicity（目標 v_SM，固定 k=8）")
V_T = ANCH["v_SM"]
rows = []
base8 = M_PL * G1 ** 8
for sn, sv in S18.items():
    v = base8 * sv
    rows.append((sn, v, abs(v / V_T - 1) * 100))
rows.sort(key=lambda r: r[2])
for sn, v, dev in rows[:5]:
    print(f"  M_Pl γ₁⁸ × {sn:<8} = {v:>14.6g} MeV  dev {dev:.4f}%")
n1 = sum(1 for _, _, d in rows if d < 1.0)
print(f"  dev<1% 嘅 S 數 = {n1}/18（注意 S18 內 N_d = L/2 = 44 係同一數值兩個標籤）")

# density：隨機目標域 [10, 1000] GeV（事前固定），thr = N_d 命中自身 dev
thr = rows[0][2] if rows[0][0] == "N_d" else dict((s, d) for s, _, d in rows)["N_d"]
rng = random.Random(20261031)
NMC = 500_000
samp = 10 ** np.array([rng.uniform(np.log10(10e3), np.log10(1000e3)) for _ in range(NMC)])  # MeV
vals8 = np.sort(np.array([base8 * sv for sv in S18.values()]))
idx = np.clip(np.searchsorted(vals8, samp), 1, len(vals8) - 1)
d = np.minimum(np.abs(vals8[idx - 1] / samp - 1), np.abs(vals8[idx] / samp - 1)) * 100
p1 = float((d < thr).sum()) / NMC
print(f"  單點 density（50 萬隨機目標，域 10–1000 GeV，thr {thr:.4f}%）：P1 = {p1*100:.2f}%")

# ---------- T3 ----------
print("\n[T3] 機制 M1「γ₁ = 每跨一層同步邊界」")
k_mech = np.log(V_T / M_PL) / LOGG1
print(f"  k_mech = log(v/M_Pl)/log(γ₁) = {k_mech:.6f}")
print(f"  |k_mech − 8| = {abs(k_mech-8):.4f} ⟹ {'🟢 M1 支持' if abs(k_mech-8) < 0.05 else '🔴 M1 否證（γ₁ 唔係每層因子）'}")
k_mech_l = np.log(V_T / M_PL) / np.log(1 / G1)   # 若用倒數
print(f"  參考：log(v/M_Pl)/log(1/γ₁) = {k_mech_l:.6f}（負值方向）")

# ---------- T4 ----------
print("\n[T4] 機制 M2「指數 = dim su(N_gen) = N_gen²−1 = 8」")
t4 = {a: T1[a]["k"] for a in ["Λ_IDCM", "M0", "m_e", "v_SM"]}
print(f"  best k：Λ_IDCM {t4['Λ_IDCM']}, M0 {t4['M0']}, m_e {t4['m_e']}, v_SM {t4['v_SM']}")
support = sum(1 for a in ["Λ_IDCM", "M0"] if t4[a] == 8)
print(f"  其他 colour-dressed 標度（Λ, M0）best k = 8 嘅數 = {support}/2 ⟹ "
      f"{'🟡 部分支持' if support else '🔴 M2 唔被支持（k=8 只屬 v）'}")

# ---------- T5 ----------
print("\n[T5] identity grammar multiplicity（8 嘅身份罕有嗎）")
ATOMS = [1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 18, 21, 29, 44, 76, 88]
exprs = {n: set() for n in range(1, 31)}
def add(n, s):
    if 1 <= n <= 30:
        exprs[int(n)].add(s)
for a in ATOMS:
    add(a, str(a))
    add(a + 1, f"{a}+1"); add(a - 1, f"{a}-1"); add(2 * a, f"2·{a}")
    add(a ** 2, f"{a}²"); add(a * (a - 1) // 2, f"C({a},2)")
    if a % 2 == 0:
        add(a // 2, f"{a}/2")
for a in ATOMS:
    for b in ATOMS:
        add(a + b, f"{a}+{b}")
        if a > b:
            add(a - b, f"{a}−{b}")
        if a * b <= 200:
            add(a * b, f"{a}·{b}")
        if b and a % b == 0 and a != b:
            add(a // b, f"{a}/{b}")
counts = {n: len(exprs[n]) for n in range(1, 31)}
rank = sorted(counts, key=lambda n: -counts[n])
print(f"  每個整數嘅相異 identity 數（1–30）：")
for n in range(1, 31):
    print(f"    {n:>2}: {counts[n]:>4}" + ("   ← 8" if n == 8 else ""))
pos8 = rank.index(8) + 1
med = np.median(list(counts.values()))
print(f"  8 嘅排名 = {pos8}/30（count {counts[8]}）；中位數 {med:.0f}")
print(f"  判：{'🔴 身份多樣性冗餘（8 嘅 identity 數 ≥ 中位數 ⟹ 「8 = X」冇資訊）' if counts[8] >= med else '🟡 8 嘅身份相對罕有'}")
print(f"  8 嘅 identity（前 12 個）：{sorted(exprs[8])[:12]}")

# ---------- T6 ----------
print("\n[T6] leave-one-out：用其他 anchor 預測 v 嘅 k")
others = [(a, T1[a]["k"]) for a in ANCH if a != "v_SM"]
# 最近鄰：以 log10(anchor/M_Pl) 為 x，best-k 為 y
xs = np.array([np.log10(ANCH[a] / M_PL) for a, _ in others])
ys = np.array([k for _, k in others], float)
xv = np.log10(ANCH["v_SM"] / M_PL)
# 線性 fit
A_ = np.vstack([xs, np.ones_like(xs)]).T
coef, *_ = np.linalg.lstsq(A_, ys, rcond=None)
k_pred_lin = coef[0] * xv + coef[1]
# 最近兩點平均
order = np.argsort(np.abs(xs - xv))[:2]
k_pred_nn = float(np.mean(ys[order]))
print(f"  線性 fit：k_pred = {k_pred_lin:.3f}（observed best k = {T1['v_SM']['k']}）")
print(f"  最近兩 anchor 平均：k_pred = {k_pred_nn:.3f}")
scat = float(np.std(ys))
print(f"  best-k 離散度 std = {scat:.3f}（8 anchor）")
print(f"  判：{'🔴 LOO 無判別力（離散度過大，k 分佈 0–13）' if scat > 1.0 else ('🟡 LOO 預測接近 8' if abs(k_pred_lin-8) < 1.0 else '🔴 LOO 預測唔支持 8')}")

res = {"meta": {"script": "15_exponent8", "date": "2026-10-24", "seed": 20261031,
                "family_size": int(len(KS) * len(S18)), "gamma1": float(G1)},
       "T1_cross_scale": {a: {k: (int(v) if k == "k" else (float(v) if isinstance(v, (int, float, np.floating)) else v))
                              for k, v in d.items()} for a, d in T1.items()},
       "T1_n_bestk8": n8,
       "T2": {"n_S_under_1pct": n1, "thr_pct": float(thr), "density_P1": p1,
              "top5": [{"S": s, "value": float(v), "dev_pct": float(dd)} for s, v, dd in rows[:5]]},
       "T3": {"k_mech": float(k_mech), "abs_diff": float(abs(k_mech - 8)),
              "verdict": "green" if abs(k_mech - 8) < 0.05 else "refuted"},
       "T4": {"best_k": t4, "support": support, "verdict": "partial" if support else "unsupported"},
       "T5": {"counts": counts, "rank8": pos8, "median": float(med),
              "identities8": sorted(exprs[8])[:20],
              "verdict": "vacuous" if counts[8] >= med else "rare"},
       "T6": {"k_pred_linear": float(k_pred_lin), "k_pred_nn": k_pred_nn,
              "scatter_std": float(np.std(ys)),
              "verdict": "no_discriminating_power" if float(np.std(ys)) > 1.0 else ("partial" if abs(k_pred_lin - 8) < 1.0 else "unsupported")}}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")
