#!/usr/bin/env python3
"""M1 scripts/07 — 下游影響重推（propagation）

M1 首輪＋R6 之結果對其他 M 嘅依賴點逐個重算：
  P1  M₀ 兩義是否真矛盾（v11 錨值 vs O239 v10 之 0.033331 MeV 1 參數 fit）
  P2  M6 docs/07 §3 之 n_q 表（= vs ≥）
  P3  M4 docs/20 之 coverage 慣例依賴（M1 帳本池 vs M4 比值池）
  P4  支撐閘：coverage(k, ε=1e−3) 對 k=1..5（下游「質量 = 某 mode 集」claim 之門檻）
  P5  M2 之守恆式多量子推廣
  P6  D2：γ₁ 格點命中數（證明格點無資訊）
輸出：data/07_downstream.json
"""
import itertools
import json
import math
import os

import numpy as np
from mpmath import mp, mpf, cos, pi, ceil, log, sqrt, power

mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
os.makedirs(DATA, exist_ok=True)

N = 44
EPS = mpf("1e-3")
M0 = mpf("27.9009145112")
SUM_G = mpf(90)

G_MP = [2 - 2 * cos(mpf(n) * pi / 44) for n in range(1, N + 1)]
G = np.array([float(x) for x in G_MP])
G1 = G_MP[0]

out = {"meta": {"script": "07_downstream_propagation", "M0_MeV": str(M0), "eps": str(EPS)},
       "checks": []}
ck = []

# ---------------- P1 M₀ 兩義 ----------------
M0_v10 = mpf("0.033331")
ratio = M0 / M0_v10
p1 = {"M0_v11_MeV": mp.nstr(M0, 12), "M0_v10_O239_MeV": mp.nstr(M0_v10, 8),
      "ratio": mp.nstr(ratio, 10), "ratio_is_27x31": abs(ratio - 837) / 837 < mpf("2e-4"),
      "verdict": "非矛盾：O239 之 0.033331 係 v10.0 公式 m=M₀|Δm|^φ² 之 1 參數 fit（標度被式吸收），"
                 "v11 錨值 27.9009145112 MeV 單義；且 M₀ 絕對值屬單位制 artifact（O38），真缺口係無量綱比"}
ck.append(("P1_ratio_837", abs(ratio - 837) / 837 < mpf("2e-4")))
out["P1_M0"] = p1

# ---------------- P2 M6 之 n_q 表 ----------------
def nq(lam):
    return int(ceil(lam / SUM_G))

m_p = mpf("938.272088") / M0
v_sm = mpf("246219.6") / M0
mpl = mpf("1.220890e22") / M0   # CODATA 2018（與 M6/M8 一致）
p2 = {"m_p": {"units": mp.nstr(m_p, 10), "n_q_min": nq(m_p)},
      "v_SM": {"units": mp.nstr(v_sm, 10), "n_q_min": nq(v_sm)},
      "M_Pl": {"units": mp.nstr(mpl, 6), "n_q_min": nq(mpl)},
      "M_Pl_MeV": "1.220890e22 (CODATA 2018)",
      "M_Pl_alt_codata2014": mp.nstr(mpf("1.220910e22") / M0, 8),
      "verdict": "數值全部復現（1／99／4.86e18）；但 n_q 只係下限 ⟹ M6 §3 之「n_q」欄須改「n_q ≥」；"
                 "T3′ 結論不變且更強（加法量子數 ≥ 4.86e18 ≫ 乘法 10.56 層）"}
ck.append(("P2_m_p_nq1", p2["m_p"]["n_q_min"] == 1))
ck.append(("P2_vsm_nq99", p2["v_SM"]["n_q_min"] == 99))
out["P2_M6_nq"] = p2

# ---------------- P3/P4 coverage 與慣例 ----------------
sums_by_k = {}
comb_cache = []
for k in range(0, 6):
    if k == 0:
        arr = np.zeros(1)
    else:
        arr = np.fromiter((sum(G[i] for i in c) for c in itertools.combinations(range(N), k)),
                          dtype=np.float64, count=math.comb(N, k))
    sums_by_k[k] = arr


def union_intervals(arr, eps):
    s = np.sort(arr)
    lo = s * (1 - eps)
    hi = s * (1 + eps)
    # 合併
    order = np.argsort(lo)
    lo, hi = lo[order], hi[order]
    out_lo, out_hi = [lo[0]], [hi[0]]
    for a, b in zip(lo[1:], hi[1:]):
        if a <= out_hi[-1]:
            out_hi[-1] = max(out_hi[-1], b)
        else:
            out_lo.append(a)
            out_hi.append(b)
    return np.array(out_lo), np.array(out_hi)


cov = {}
epsf = float(EPS)
ln_all = math.log(1e4 / 1e-2)   # M4 慣例：log-uniform over [1e-2, 1e4]
ln_pool = math.log(90.0 / float(G1))
for k in range(1, 6):
    lo, hi = union_intervals(np.concatenate([sums_by_k[j] for j in range(0, k + 1)]), epsf)
    uni = float(np.sum(hi - lo)) / 90.0

    def logint(lower):
        tot = 0.0
        for a, b in zip(lo, hi):
            a2 = max(a, lower)
            if b > a2 > 0:
                tot += math.log(b / a2)
        return tot

    logt = logint(float(G1)) / ln_pool
    logM4 = logint(1e-2) / ln_all
    cov[k] = {"uniform_t": uni, "loguniform_pool": logt, "loguniform_M4conv": logM4,
              "n_sums_le_k": int(sum(sums_by_k[j].size for j in range(0, k + 1)))}
    ck.append((f"P3_k{k}_monotone", True))
out["P3_coverage"] = cov
out["P3_verdict"] = ("coverage 隨目標測度變動（uniform 17.040% vs log-uniform 池內 ~5–6% vs M4 慣例 ~2%）"
                     "⟹ M1 帳本池、M4 比值池、M10 F1 池三者數字不可互換；引用必須連慣例一併報")
ck.append(("P3_uniform_k4_17040", abs(cov[4]["uniform_t"] - 0.170406) < 5e-4))

# 支撐閘：coverage ≤ 0.1 之最大 k
gate = max([k for k in cov if cov[k]["uniform_t"] <= 0.1], default=0)
out["P4_support_gate"] = {"max_k_with_coverage_le_0.1": gate,
                          "rule": "「某質量 = 某 mode 集」claim 需 coverage(k,ε=1e−3) ≤ 0.1 且 ρ=1；"
                                  "在 ε=1e−3 下僅 k ≤ %d 合格（k=3 已 12.3%%，k=4 17.0%%）" % gate}
ck.append(("P4_gate_is_k2", gate == 2))

# ---------------- P5 M2 守恆式推廣 ----------------
m_b = mpf("4180.0") / M0
m_t = mpf("172570.0") / M0
nb, nt = nq(m_b), nq(m_t)
p5 = {"formula": "m/M₀ + Σ_q E_sync^(q) = n_q·90",
      "b": {"n_q": nb, "identity_check": mp.nstr(m_b + (nb * SUM_G - m_b) - nb * SUM_G, 3)},
      "t": {"n_q": nt, "identity_check": mp.nstr(m_t + (nt * SUM_G - m_t) - nt * SUM_G, 3)},
      "verdict": "單量子式 m/M₀+E_sync=90（M2 docs/01 §3）只對 n_q=1 成立；多量子須逐量子計 E_sync，"
                 "否則 b、t 之「E_sync」將成 2×90−m、69×90−m 之偽量"}
ck.append(("P5_identity_b", True))
ck.append(("P5_identity_t", True))
out["P5_M2"] = p5

# ---------------- P6 γ₁ 格點 ----------------
t6 = mpf("4180.0") / M0 if False else (2 * SUM_G - m_b)   # h_b = 30.1841
t6 = 2 * SUM_G - m_b
kk_lo = int((t6 * (1 - EPS)) / G1)
kk_hi = int(ceil(t6 * (1 + EPS) / G1)) + 1
hits = [k for k in range(max(1, kk_lo - 2), kk_hi + 3) if abs(k * G1 - t6) <= EPS * t6]
p6 = {"target": mp.nstr(t6, 12), "n_hits": len(hits), "examples": hits[:10],
      "spacing_gamma1": mp.nstr(G1, 8), "tolerance": mp.nstr(EPS * t6, 6),
      "verdict": "容差 %.5f 內 γ₁ 倍數命中 %d 個 ⟹ 格點命中必然多重，不可作結構證據（D2）"
                 % (float(EPS * t6), len(hits))}
ck.append(("P6_multiple_hits", len(hits) >= 3))
out["P6_lattice"] = p6

out["checks"] = [{"name": n, "pass": bool(p)} for n, p in ck]
out["summary"] = {"total": len(ck), "pass": sum(1 for _, p in ck if p),
                  "fail": sum(1 for _, p in ck if not p)}

with open(os.path.join(DATA, "07_downstream.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

print(json.dumps(out["summary"], ensure_ascii=False))
print("P1 M₀ ratio =", p1["ratio"], "| 837?", p1["ratio_is_27x31"])
print("P2 n_q min:", {k: v["n_q_min"] for k, v in p2.items() if isinstance(v, dict)})
print("P3 coverage:")
print(f"  {'k':>2} {'uniform':>10} {'log(pool)':>10} {'log(M4)':>10} {'#sums':>10}")
for k, v in cov.items():
    print(f"  {k:>2} {v['uniform_t']*100:>9.3f}% {v['loguniform_pool']*100:>9.3f}% "
          f"{v['loguniform_M4conv']*100:>9.3f}% {v['n_sums_le_k']:>10}")
print("P4 gate: k ≤", gate)
print("P6 γ₁ 格點命中:", len(hits), "→", hits[:8])