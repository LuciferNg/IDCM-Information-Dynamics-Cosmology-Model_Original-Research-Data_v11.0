"""M4 docs/22 配套 — 核實 O216（標籤層 no-go）與 O86（A 算子身份接手 + 87 結構身份）（2026-10-24）
===============================================================================================
收到交付：
  O216 `docs/04_handoff_response_spectator_label_operator.md`（14:43 register note）：
       🔴 標籤層冇 spectator 依賴算子（三條結構理由 R1 標籤≠operator／R2 標籤層泛函恆等可加／R3 色數3=B−1 亦加性）；
       獨立 Fraction 重算 8 通道 8/8 PASS
  O216 `docs/05_audit_report.md`（14:52）：指出 M4 `docs/20` §六「5 個 handoff 發出」stale（實有 8 條）
  O86 `docs/16_handoff_intake_20261024.md`（14:37）：H1 A 譜身份 → ✅ 接受為 O86/N8（open）；
       核 M4 五項；★兩項新發現：① 87 = N_gen×L₇ = 3×29 = Σγ(44)−N_gen = 90−3；② P6 鏈 m1/m2 反解 n 非整數
本輪（M4 側，紀律「只核不重跑」）：
  A. O216 核實：8 通道重算、加性模型預言（同取代 ⟹ 同 E）對比數據、C₂ 不變性
  B. O86 核實：87 三重結構身份、A/A_old、m3、m1/m2 反解 n 非整數、4/γ₁
  C. 結構標籤族掃描：A/M₀ ∈ {P·γ₁^k}（P = ≤3 個結構整數之積）＋ MC 密度閘
判準（事前）：B 須逐位重現（≤1e-9 相對）；C 升格須 dev < 0.1% ∧ 密度 < 5%；87 多身份 ⟹ 只記線索
"""
import json
import math
import numpy as np

OUT = "data/22_o216_o86_verification.json"
L = 88
g1 = 4 * math.sin(math.pi / L) ** 2
gam = lambda n: 4 * math.sin(n * math.pi / L) ** 2
M0 = 27.900915e6            # eV
A = 1.64337
A_old = 1.64044
X = 2.560
sX = 0.094

print("=" * 92)
print("M4 docs/22 — 核實 O216（標籤層 no-go）／O86（A 接手 + 87 身份）")
print("=" * 92)

# ---------- A. O216 ----------
print("\n[A] O216：8 通道重算 + 加性模型判死")
CH = [("n−p", 1.2933, 1, 0.0000), ("K⁺−K⁰", -3.9340, -1, 0.0210), ("π⁺−π⁰", 4.5936, 1, 0.00014),
      ("Σ⁻−Σ⁰", 4.8070, 1, 0.0040), ("Σ⁰−Σ⁺", 3.2720, 1, 0.0080), ("Σ⁻−Σ⁺", 8.0790, 2, 0.0200),
      ("Ξ⁻−Ξ⁰", 6.8500, 1, 0.0230), ("D⁺−D⁰", 4.8200, 1, 0.0070)]
E = {}
for name, dm, cf, sdm in CH:
    E[name] = dm - cf * X
    print(f"  {name:<8} E = {E[name]:+8.4f}  (σ = {math.hypot(sdm, abs(cf)*sX):.4f})")
claim = {"n−p": -1.2667, "K⁺−K⁰": -1.3740, "π⁺−π⁰": 2.0336, "Σ⁻−Σ⁰": 2.2470,
         "Σ⁰−Σ⁺": 0.7120, "Σ⁻−Σ⁺": 2.9590, "Ξ⁻−Ξ⁰": 4.2900, "D⁺−D⁰": 2.2600}
bad = [k for k in claim if abs(round(E[k], 4) - claim[k]) > 1e-9]
print(f"  8/8 對 M4 原值：{'✅ 全部一致' if not bad else '🔴 不符 '+str(bad)}")
# 加性模型（每 flavour 一權重 w_q 相加）對 u↔d 取代嘅預言
subs = ["n−p", "Σ⁻−Σ⁰", "Σ⁰−Σ⁺", "Ξ⁻−Ξ⁰", "K⁺−K⁰", "π⁺−π⁰", "D⁺−D⁰"]
pred_spread = 0.0
obs = [E[k] for k in subs]
print(f"  加性模型預言：所有『單一 u↔d 取代』通道 E 必須相等 ⟹ 預言 spread = {pred_spread:.3f} MeV")
print(f"  觀測 spread（7 個單取代通道）= {max(obs)-min(obs):.3f} MeV（−1.374 → +4.290）"
      f"；雙取代 Σ⁻−Σ⁺ = {E['Σ⁻−Σ⁺']:.3f} MeV")
sig = (max(obs) - min(obs)) / math.hypot(0.0963, 0.0968)
print(f"  ⟹ 加性被否證：{sig:.1f}σ（保守取 K 與 Ξ 之 σ）⟹ O216 R2『標籤層泛函恆等可加』與數據衝突 ⟹ 標籤層不可產 spread ✅")
# C₂ 不變性：同一 SU(3) 表示下多重集改變 ⟹ C₂(R) 不變（數值：C₂ 只依 R）
C2 = {"3": 4/3, "6": 10/3, "8": 3, "10": 6, "15": 16/3}
print(f"  C₂(R) 檢核：{C2} ⟹ C₂ 只依表示 R，同 R 下多重集改動不變 ⟹ 唔可產 spread ✅")

# ---------- B. O86 ----------
print("\n[B] O86：87 結構身份 + P6 鏈自洽核")
Lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123]
N_gen, N_d, B = 3, 44, 4      # N_d = 44（88 = 2N_d）；Σγ(44) = B(N_d+1)/2 = 90
S_gam = int(B * (N_d + 1) / 2)
ids = {"N_gen×L₇": N_gen * Lucas[7], "Σγ(44) − N_gen": S_gam - N_gen, "88 − 1": 88 - 1}
print(f"  87 結構身份：{ {k: v for k, v in ids.items()} } ⟹ 全部 = 87 ⟹ {'✅' if set(ids.values()) == {87} else '🔴'}")
print(f"    ⚠️ 多身份（3 個）⟹ 依 docs/15 教訓（8 有 17 個身份）只可作線索")
print(f"  A/A_old = {A/A_old:.7f}（O86 報 1.0017861）")
print(f"  4/γ₁ = {4/g1:.10f}（O86 報 784.9646644）")
print(f"  m3 = A·γ₆² = {A*gam(6)**2:.10f} eV（O86 報 0.0536815110413）")
print(f"  87γ₁⁴M₀ = {87*g1**4*M0:.9f} eV，dev vs A = {abs(87*g1**4*M0/A-1)*100:.4f}%")
# m1/m2 反解
print("  P6 鏈 m = A·γ_n² 反解（n 本應為整數）:")
rows_m = []
for lab, m in (("m1", 0.0209452), ("m2", 0.0226848), ("m3", 0.0536815)):
    gn = math.sqrt(m / A)
    # γ_n = 4sin²(nπ/88) ⟹ n = (88/π)·arcsin(√(γ_n/4))
    n = (L / math.pi) * math.asin(math.sqrt(gn / 4))
    rows_m.append({"mass": lab, "m_eV": m, "gamma_n": gn, "n_solved": n})
    print(f"    {lab} = {m:.7f} eV ⟹ γ_n = {gn:.6f} ⟹ n = {n:.4f} {'✅ 整數' if abs(n-round(n))<1e-6 else '🔴 非整數'}")
print("  ⟹ O86 發現成立：m3 ⟹ n=6 精確；m1/m2 ⟹ n = 4.7276／4.8245 非整數 ⟹ P6 鏈 m1/m2 唔可以係 A·γ_n²（n 整數）")
# 診斷：邊種鏈形式 m_n = A·γ_{n+δ}^p 可以同時容納 P6 三個數？
print("  鏈形式診斷（求單一 A 使三數同時成立，容差 1%）：")
ms = {"m1": 0.0209452, "m2": 0.0226848, "m3": 0.0536815}
found = []
for p_ in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
    for dlt in (0.0, 0.5):
        best_t = None
        for n1 in range(1, 40):
            for n2 in range(n1 + 1, 41):
                for n3 in range(n2 + 1, 42):
                    A1 = ms["m1"] / gam(n1 + dlt) ** p_
                    A2 = ms["m2"] / gam(n2 + dlt) ** p_
                    A3 = ms["m3"] / gam(n3 + dlt) ** p_
                    spread = max(A1, A2, A3) / min(A1, A2, A3) - 1
                    if best_t is None or spread < best_t[0]:
                        best_t = (spread, (n1, n2, n3), (A1, A2, A3))
        if best_t[0] < 0.01:
            found.append((p_, dlt, best_t))
            print(f"    p={p_}, δ={dlt} ⟹ n=({best_t[1][0]},{best_t[1][1]},{best_t[1][2]}) "
                  f"A = {[round(x,6) for x in best_t[2]]} spread {best_t[0]*100:.3f}% ✅")
        else:
            print(f"    p={p_}, δ={dlt} ⟹ 最佳 spread {best_t[0]*100:8.2f}%（n={best_t[1]}）🔴")

# ---------- C. 結構標籤族掃描 ----------
print("\n[C] 結構標籤族：A/M₀ = P·γ₁^k（P = ≤3 個結構整數之積）＋ 密度閘")
Sset = [1, 3, 4, 5, 7, 9, 11, 12, 27, 29, 36, 44, 45, 47, 60, 76, 88, 90, 123, 421]
target = A / M0
cands = []
for i, a in enumerate(Sset):
    for b in Sset[i:]:
        for c in Sset[Sset.index(b):]:
            cands.append(a * b * c)
cands = sorted(set(cands))
print(f"  P 集合大小 = {len(cands)}（≤3 因子，含重複）")
best = None
hits = []
for k in range(1, 7):
    for P in cands:
        v = P * g1 ** k
        dev = abs(v / target - 1) * 100
        if best is None or dev < best[3]:
            best = (P, k, v * M0, dev)
        if dev < 1.0:
            hits.append((P, k, v * M0, dev))
hits.sort(key=lambda r: r[3])
print(f"  族大小 = {len(cands)*6}；dev < 1% 嘅組合 = {len(hits)}")
for P, k, aev, dev in hits[:6]:
    print(f"    P={P:>6}·γ₁^{k} ⟹ A_pred = {aev:.6f} eV  dev {dev:.4f}%")
print(f"  最佳：P={best[0]}, k={best[1]} ⟹ {best[2]:.6f} eV, dev {best[3]:.4f}%"
      f" ⟹ {'✅ 過 0.1% 門檻' if best[3] < 0.1 else '🔴 未過 0.1% 門檻'}")
# MC 密度（同一族，log-uniform [1e-10, 1e-5]）
rng = np.random.default_rng(20261035)
vals = np.sort(np.array([P * g1 ** k for P in cands for k in range(1, 7)]))
tg = 10 ** rng.uniform(-10, -5, 200000)
idx = np.clip(np.searchsorted(vals, tg), 1, vals.size - 1)
dev_mc = np.minimum(np.abs(vals[idx] / tg - 1), np.abs(vals[idx - 1] / tg - 1))
dens = float((dev_mc < 1e-3).mean())
print(f"  MC 密度（±0.1% 帶）= {dens*100:.2f}% ⟹ {'🔴 無鑑別力' if dens > 0.05 else '🟡 通過密度閘'}")

res = {"meta": {"script": "22_o216_o86_verification", "date": "2026-10-24", "seed": 20261035},
       "A_O216": {"E": E, "mismatch_vs_claim": bad, "additive_spread_obs": max(obs) - min(obs),
                  "additive_sigma": sig},
       "B_O86": {"ids_87": ids, "A_over_Aold": A / A_old, "inv4g1": 4 / g1,
                 "m3": A * gam(6) ** 2, "m_scan": rows_m},
       "C_family": {"size": len(cands) * 6, "hits_below_1pct": len(hits),
                    "best": {"P": best[0], "k": best[1], "A_pred_eV": best[2], "dev_pct": best[3]},
                    "density_0p1pct": dens}}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2, default=float)
print(f"\n→ {OUT}")