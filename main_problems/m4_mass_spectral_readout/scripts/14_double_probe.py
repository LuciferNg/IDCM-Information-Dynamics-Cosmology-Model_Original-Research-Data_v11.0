"""M4 G5 — 雙 probe 讀數差（em/strong 反相）算子化第一步（2026-10-24）
=====================================================================
上游：docs/05 §四 [D]（π⁺−π⁰ > 0 但 K⁺−K⁰ < 0，em 與 strong-isospin 反相）、
      docs/06 §四 [C]（|K/(n−p)| = 3.04415，hop 比 3.0 dev 1.472%）、docs/13 §六 gap G5。

問題：反相係要一個「雙通道混合算子」，抑或只係「單一 em 讀數 + 夸克內容符號」？

Protocol（事前寫死，唔准見結果改）：
  P0 數據（PDG 2024，pole/charged-neutral 質量差；一律引用值 + 誤差）：
       Δπ ≡ m(π⁺) − m(π⁰) = +4.59359 ± 0.00014 MeV   （純 em 通道）
       ΔK ≡ m(K⁺) − m(K⁰) = −3.934 ± 0.021 MeV
       ΔN ≡ m(n)  − m(p)  = +1.293331 ± 0.0000005 MeV
       lattice 對照：m_d − m_u = 2.56 ± 0.094 MeV（FLAG24 Nf=2+1+1：4.70(5) − 2.14(8)）
  P1 假設 H1：K 與 N 屬同一 |ΔQ| = 1 通道 ⟹ em 讀數相同 E_K = E_N ≡ E
       ⟹ E = (ΔN + ΔK)/2 ；m_d−m_u = ΔN − E = (ΔN − ΔK)/2
       ⟹ 判準：dev vs lattice < 1% 且 < 2σ ⟹ 一致（🟡，後驗 premise）
  P2 T2：同一 ΔQ² = 1 但 π 通道 em 讀數 E_π = Δπ ⟹ 比值 R_em = E_π/|E|
       ⟹ 判準：R_em 落 pre-registered 族 <0.1% 且 density PN < 5% 才可升格
  P3 T3：反相結構條件 —— m_d−m_u > |E| ⟹ ΔN > 0 且 ΔK < 0（自動反相）；
       判：若 |E| < m_d−m_u 且 E < 0 則「單一 em 讀數 + 夸克內容符號」足夠（唔需要雙通道算子）
  P4 T4：E 嘅 IDCM 譜形式 —— 族 F = { anchor × γ₁^k × S }（8 anchors × 91 k × 18 S = 13104）
       ＋ density（同族對隨機目標，thr = 命中自身 dev）
輸出：data/14_double_probe.json ＋ stdout
"""
import json
import numpy as np
import random

OUT = "data/14_double_probe.json"

# ---------- 譜 ----------
L = 88
N = np.arange(1, L // 2 + 1)                      # 1..44
GAM = 4 * np.sin(N * np.pi / L) ** 2
GAM2 = GAM ** 2
G1 = GAM[0]

# ---------- P0 數據 ----------
D_PI = 4.59359;   S_PI = 0.00014
D_K  = -3.934;    S_K  = 0.021
D_N  = 1.293331;  S_N  = 0.0000005
MD_MU_LAT = 2.56; S_MD_MU_LAT = 0.094

print("=" * 78)
print("M4 G5 — 雙 probe 讀數差（em / strong-isospin 反相）")
print("=" * 78)
print("\n[P0] 數據（PDG 2024 + FLAG24）")
print(f"  Δπ = m(π⁺)−m(π⁰) = {D_PI:+.5f} ± {S_PI:.5f} MeV   (純 em 通道)")
print(f"  ΔK = m(K⁺)−m(K⁰) = {D_K:+.5f} ± {S_K:.5f} MeV")
print(f"  ΔN = m(n)−m(p)   = {D_N:+.6f} ± {S_N:.7f} MeV")
print(f"  lattice  m_d−m_u = {MD_MU_LAT:.4f} ± {S_MD_MU_LAT:.4f} MeV (FLAG24 2+1+1)")

# ---------- P1 ----------
E = (D_N + D_K) / 2
MDMU_DER = (D_N - D_K) / 2
dev1 = abs(MDMU_DER / MD_MU_LAT - 1) * 100
sig1 = np.hypot(np.hypot(S_N, S_K) / 2, S_MD_MU_LAT)
z1 = (MDMU_DER - MD_MU_LAT) / sig1
print("\n[P1] H1（K/N 同 |ΔQ|=1 通道 ⟹ E_K = E_N = E）")
print(f"  E = (ΔN + ΔK)/2 = {E:+.6f} MeV")
print(f"  反解 m_d−m_u = (ΔN − ΔK)/2 = {MDMU_DER:+.6f} MeV")
print(f"  vs lattice {MD_MU_LAT:.4f} ± {S_MD_MU_LAT:.4f}: dev {dev1:.4f}%  = {z1:+.2f}σ")
print(f"  判：{'🟡 一致（1% 內、2σ 內）' if (dev1 < 1 and abs(z1) < 2) else '🔴 不一致'}")

# ---------- P2 ----------
R_EM = D_PI / abs(E)
print("\n[P2] T2：同 ΔQ²=1 通道嘅 em 讀數比 R_em = Δπ/|E|")
print(f"  Δπ = {D_PI:+.5f} MeV ; |E| = {abs(E):.6f} MeV ⟹ R_em = {R_EM:.6f}")
print(f"  ⟹ π 通道 em 讀數係 K/N 通道嘅 {R_EM:.4f}×（正負亦相反）")

# 族 1：γ 比池（docs/09 protocol，1936 值）
POOL_G = np.sort(np.array([a / b for a in GAM for b in GAM]))
# 族 2：結構整數比池（docs/11 T4 同族）
SI = [1, 2, 3, 4, 5, 7, 8, 11, 13, 21, 29, 44, 76, 88]
POOL_S = np.array(sorted({p / q for p in SI for q in SI if p != q}))
best_g = POOL_G[np.argmin(np.abs(POOL_G / R_EM - 1))]
best_s = POOL_S[np.argmin(np.abs(POOL_S / R_EM - 1))]
dev_g = abs(best_g / R_EM - 1) * 100
dev_s = abs(best_s / R_EM - 1) * 100
print(f"  γ 比池最佳 {best_g:.6f} dev {dev_g:.4f}% | 結構整數比池最佳 {best_s:.6f} dev {dev_s:.4f}%")

rng = random.Random(20261030)
NMC = 1_000_000
dom = np.log10([1.0, 10.0])
samp = 10 ** np.array([rng.uniform(*dom) for _ in range(NMC)])

def density(pool, thr, targets, n_ratios=1):
    idx = np.searchsorted(pool, targets)
    idx = np.clip(idx, 1, len(pool) - 1)
    lo, hi = pool[idx - 1], pool[idx]
    d = np.minimum(np.abs(lo / targets - 1), np.abs(hi / targets - 1)) * 100
    p1 = float((d < thr).sum()) / len(targets)
    pn = 1 - (1 - p1) ** n_ratios
    return p1, pn

p1g, png = density(POOL_G, dev_g, samp)
p1s, pns = density(POOL_S, dev_s, samp)
print(f"  density（10⁶，域 1–10，thr = 命中自身 dev）：")
print(f"    γ 比池 1936 值  : P1 {p1g*100:.2f}%  PN {png*100:.2f}%")
print(f"    結構整數比 158 值: P1 {p1s*100:.2f}%  PN {pns*100:.2f}%")
p2_verdict = "🟢 可辨識" if (min(dev_g, dev_s) < 0.1 and max(p1g, p1s) < 0.05) else "🔴 唔可辨識（密度殺）"
print(f"  判：{p2_verdict}")

# ---------- P3 ----------
print("\n[P3] T3：反相結構條件")
cond = (E < 0) and (MDMU_DER > abs(E))
print(f"  E = {E:+.6f} MeV (<0) ; m_d−m_u = {MDMU_DER:.6f} > |E| = {abs(E):.6f} ⟹ {cond}")
print(f"  ΔN = (m_d−m_u) + E = {MDMU_DER + E:+.6f} MeV (觀察 {D_N:+.6f})")
print(f"  ΔK = −(m_d−m_u) + E = {-(MDMU_DER) + E:+.6f} MeV (觀察 {D_K:+.5f})")
print(f"  判：{'✅ 反相 = 單一 em 讀數 (E<0) + 夸克內容符號翻轉；唔需要雙通道算子' if cond else '🔴 反相需要額外通道'}")

# ---------- P4 ----------
print("\n[P4] T4：E 嘅 IDCM 譜形式（族 {anchor × γ₁^k × S}，全部單位 MeV）")
# 全部統一 MeV（1 eV = 1e-6 MeV；m_e = 0.51099895 MeV；Λ/M0/v/M_Pl 皆 MeV）
ANCH = {"1 eV": 1e-6, "m_e": 0.51099895, "Λ_IDCM": 45.28487, "M0": 27.900915,
        "v_SM": 246219.6, "sqrt(MPl·Λ)": np.sqrt(1.220910e22 * 45.28487),
        "M_Pl": 1.220910e22, "Λ²/M0": 45.28487 ** 2 / 27.900915}
S18 = {"1": 1, "B": 4, "N_d": 44, "L": 88, "N_gen": 3, "F5": 5, "F7": 13, "F8": 21,
       "L5": 11, "L7": 29, "L9": 76, "N_gen²": 9, "L6": 18, "B·N_gen": 12,
       "L+N_d": 132, "N_d/2": 22, "L/2": 44, "√5": np.sqrt(5)}
KS = range(-45, 46)
TARGET = abs(E)
entries = []
for an, av in ANCH.items():
    for k in KS:
        base = av * G1 ** k
        for sn, sv in S18.items():
            entries.append((base * sv, an, k, sn))
vals = np.array([e[0] for e in entries])
devs = np.abs(vals / TARGET - 1) * 100
order = np.argsort(devs)
print(f"  族大小 {len(entries)} ; 目標 |E| = {TARGET:.6f} MeV")
for i in order[:5]:
    v, an, k, sn = entries[i]
    print(f"    {an} × γ₁^{k} × {sn} = {v:.6g} MeV  dev {devs[i]:.4f}%")
hits = int((devs < 0.1).sum())
print(f"  <0.1% 命中數 = {hits}")

# density：同族對隨機目標（log-均勻 [min,max] of entries）
lm = np.log10([vals.min(), vals.max()])
samp2 = 10 ** np.array([rng.uniform(*lm) for _ in range(200_000)])
logv = np.sort(np.log10(vals))
ls = np.log10(samp2)
idx = np.searchsorted(logv, ls)
idx = np.clip(idx, 1, len(logv) - 1)
d1 = np.abs(10 ** logv[idx - 1] / samp2 - 1) * 100
d2 = np.abs(10 ** logv[idx] / samp2 - 1) * 100
dd = np.minimum(d1, d2)
p4_p1 = float((dd < 0.1).sum()) / len(samp2)
print(f"  density（10⁵，域 10^[{lm[0]:.2f},{lm[1]:.2f}]，thr 0.1%）：P1 = {p4_p1*100:.2f}%")
p4_verdict = "🟢 可辨識" if (hits and p4_p1 < 0.05) else "🔴 唔可辨識"
print(f"  判：{p4_verdict}")

res = {
    "meta": {"script": "14_double_probe", "date": "2026-10-24", "seed": 20261030,
             "NMC_density": NMC, "NMC_spectral": 200000},
    "data": {"d_pi": D_PI, "s_pi": S_PI, "d_k": D_K, "s_k": S_K, "d_n": D_N, "s_n": S_N,
             "md_mu_lattice": MD_MU_LAT, "s_md_mu_lattice": S_MD_MU_LAT},
    "P1": {"E_MeV": E, "md_mu_derived": MDMU_DER, "dev_pct": dev1, "sigma": sig1, "z": z1,
           "verdict": "partial_yellow" if (dev1 < 1 and abs(z1) < 2) else "red"},
    "P2": {"R_em": R_EM, "best_gamma_pool": best_g, "dev_gamma_pool_pct": dev_g,
           "best_struct_pool": best_s, "dev_struct_pool_pct": dev_s,
           "P1_gamma": p1g, "PN_gamma": png, "P1_struct": p1s, "PN_struct": pns,
           "verdict": p2_verdict},
    "P3": {"E_negative": bool(E < 0), "m_dm_u_gt_absE": bool(MDMU_DER > abs(E)),
           "pred_dN": MDMU_DER + E, "pred_dK": -MDMU_DER + E, "verdict": "structural_green" if cond else "red"},
    "P4": {"family_size": len(entries), "hits_0p1pct": hits, "P1_density": p4_p1,
           "top5": [{"value": float(entries[i][0]), "anchor": entries[i][1], "k": int(entries[i][2]),
                     "S": entries[i][3], "dev_pct": float(devs[i])} for i in order[:5]],
           "verdict": p4_verdict},
}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")
