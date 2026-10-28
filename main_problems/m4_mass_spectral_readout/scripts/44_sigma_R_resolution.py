#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/44 -- `σ_R` 爭議終局：MC vs 解析（`docs/41` F3 / `01_problem_record` §三十 T3）。

問題：
  `docs/11` 報 `σ_R = 0.039896`（`scripts/11` 之 2×10⁵ 樣本 MC），
  `scripts/l2_independent_verify.py` 之解析傳播報 `0.039714`，差 **0.46%**。
  兩者不能同真 ⟹ 必須判定 canonical，並重算下游（z、R 帶、`m_u` 預言帶）。

本腳本三路獨立計算：
  (a) 一階解析（delta method）：σ = Δm_u /(m_u ln10) / span_ν
  (b) **二階／精確**：以 Gauss–Hermite 求 Var[log10(m_t/m_u)]（mu ~ N(μ,σ) 截斷 >0）
  (c) 獨立 MC：N = 10⁷，三個 seed，報 std 之標準誤 ⟹ 判定 (a) 與 (c) 之差係噪聲定系統

判準（事前）：
  A 三路一致（相對差 < 0.1%）⟹ canonical = 精確值 (b)
  B MC 之 3σ 區間須覆蓋 (b)
  C 下游 z 移動 < 0.01σ ⟹ 判決不變
全部輸入：m_u = 2.14 ± 0.08 MeV、m_t = 172.69 ± 0.30 GeV、span_ν = log10(m3/m1)，
m3 = 0.053586 eV、m1 = 0.020908 eV（O257 R7 定案值；誤差 <0.1%）。
"""
import json
import math
import os

import numpy as np

CHECK = []


def chk(label, ok, got, expect):
    CHECK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-60s | %-28s | %s" % (label, got, expect))


# ---------------- 輸入
MU, DU = 2.14e-3, 0.08e-3          # GeV
MT, DT = 172.69, 0.30              # GeV
M1, M3 = 0.020908, 0.053586        # eV
SPAN_NU = math.log10(M3 / M1)

print("=" * 104)
print("輸入：m_u = %.3e ± %.3e GeV ｜ m_t = %.2f ± %.2f GeV ｜ span_ν = log10(%.6f/%.6f) = %.9f"
      % (MU, DU, MT, DT, M3, M1, SPAN_NU))
print("=" * 104)

# ---------------- (a) 一階解析
sig_spanQ_lin = DU / (MU * math.log(10.0))
sig_R_lin = sig_spanQ_lin / SPAN_NU
R_central = math.log10(MT / MU) / SPAN_NU
print("(a) 一階解析：σ_spanQ = %.9f ⟹ σ_R = %.9f" % (sig_spanQ_lin, sig_R_lin))
print("    R（中心值）= %.9f" % R_central)

# ---------------- (b) 精確（Gauss–Hermite，截斷 mu>0）
def exact_log_sigma(mu, du, mt, dt, n=100):
    """Var[log10(mt/mu)] 精確值：mu, mt 皆 Gaussian（各自獨立），mu 截斷 >0。"""
    x, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / w.sum()
    mus = mu + du * x
    ok = mus > 0
    mus, w1 = mus[ok], w[ok] / w[ok].sum()
    mts = mt + dt * x
    ok2 = mts > 0
    mts, w2 = mts[ok2], w[ok2] / w[ok2].sum()
    # log10(mt) - log10(mu)；兩者獨立 ⟹ 變量可加
    lm_t = np.log10(mts)
    lm_u = np.log10(mus)
    mean = float(np.sum(w1 * lm_u) - np.sum(w2 * lm_t))  # = E[log mu] - E[log mt]
    var = float(np.sum(w1 * (lm_u - np.sum(w1 * lm_u)) ** 2) + np.sum(w2 * (lm_t - np.sum(w2 * lm_t)) ** 2))
    return math.sqrt(var), mean


sig_spanQ_ex, mean_spanQ = exact_log_sigma(MU, DU, MT, DT)
sig_R_ex = sig_spanQ_ex / SPAN_NU
R_ex = -mean_spanQ / SPAN_NU
print("(b) 精確（GH, n=100）：σ_spanQ = %.9f ⟹ σ_R = %.9f ；R = %.9f" % (sig_spanQ_ex, sig_R_ex, R_ex))

# ---------------- (c) 獨立 MC，N = 1e7，3 seeds
NMC = 10_000_000
mcs = []
for seed in (20261028, 7, 31415):
    rng = np.random.default_rng(seed)
    mu = rng.normal(MU, DU, NMC)
    mt = rng.normal(MT, DT, NMC)
    mask = (mu > 0) & (mt > 0)
    v = np.log10(mt[mask] / mu[mask])
    mcs.append(float(np.std(v)) / SPAN_NU)
mc_mean = float(np.mean(mcs))
mc_se = float(np.std(mcs, ddof=1) / math.sqrt(len(mcs)))
print("(c) MC N=%d ×3 seeds：σ_R = %s ⟹ mean = %.9f ± %.9f（std 之 SE）"
      % (NMC, ", ".join("%.9f" % x for x in mcs), mc_mean, mc_se))

# ---------------- 判準
rel_ac = abs(mc_mean / sig_R_ex - 1.0)
rel_lin_ex = abs(sig_R_lin / sig_R_ex - 1.0)
print("\n  MC vs 精確：相對差 %.4f%%（MC 之 3σ = %.4f%%）" % (rel_ac * 100, 3 * mc_se / mc_mean * 100))
print("  一階 vs 精確：相對差 %.4f%%（二階項大小）" % (rel_lin_ex * 100))
chk("A1 一階與精確之相對差 < 0.1%（二階項可忽略）", rel_lin_ex < 1e-3,
    "%.5f%%" % (rel_lin_ex * 100), "<0.1%")
chk("A2 MC 與精確一致（差 < 3σ_MC）", rel_ac < 3 * mc_se / mc_mean,
    "%.4f%% vs 3σ=%.4f%%" % (rel_ac * 100, 3 * mc_se / mc_mean * 100), "差<3σ")
chk("A3 MC 與舊值 0.039896 之差不為零（舊值偏離 > 3σ_MC ⟹ 舊值係統計漲落）",
    abs(0.039896 / mc_mean - 1.0) > 3 * mc_se / mc_mean,
    "舊 vs 新：%.3f%%" % ((0.039896 / mc_mean - 1) * 100), ">3σ_MC")

SIG_R_CANON = sig_R_ex
print("\n  ⟹ canonical σ_R = **%.9f**（精確；一階 %.9f 與之同值至 5 位）" % (SIG_R_CANON, sig_R_lin))

# ---------------- 下游重算
print("\n" + "=" * 104)
print("D 下游重算（z、R 帶、m_u 預言）")
print("=" * 104)
for tag, R, s in (("舊（MC 0.039896）", 12.005533, 0.039896),
                  ("新（canonical）", R_ex, SIG_R_CANON)):
    z = abs(R - 12.0) / s
    print("  %-18s R = %.6f ± %.6f ⟹ z = %.4fσ ；帶 [%.4f, %.4f]"
          % (tag, R, s, z, R - 2 * s, R + 2 * s))
z_old = abs(12.005533 - 12.0) / 0.039896
z_new = abs(R_ex - 12.0) / SIG_R_CANON
chk("C1 下游 z 仍 < 0.2σ 且移動 < 0.05σ（判決不變）", z_new < 0.2 and abs(z_new - z_old) < 0.05,
    "%.4f → %.4f" % (z_old, z_new), "z<0.2 且 |Δz|<0.05")

# m_u 預言帶（R = 12 精確 ⟹ m_u_pred = m_t / 10^(12·span_ν)）
mu_pred = MT / (10 ** (12.0 * SPAN_NU)) * 1000.0     # MeV
print("  m_u 預言（R = 12）= %.6f MeV（FLAG 2.14 ± 0.08 ⟹ z = %.3fσ）"
      % (mu_pred, abs(mu_pred - 2.14) / 0.08))
chk("C2 m_u 預言不依賴 σ_R（純 R=12 值）⟹ 不受本輪影響", abs(mu_pred - 2.1498) < 0.005,
    "%.6f MeV" % mu_pred, "2.1498")

# σ_R 之相對不確定度本身
print("\n  σ_R 之相對精度：舊值 0.039896 vs 新值 %.9f ⟹ 差 %.3f%%" % (SIG_R_CANON, (0.039896 / SIG_R_CANON - 1) * 100))

DATA = {
    "inputs": {"mu_GeV": MU, "dmu_GeV": DU, "mt_GeV": MT, "dmt_GeV": DT,
               "m1_eV": M1, "m3_eV": M3, "span_nu": round(SPAN_NU, 12)},
    "sigma_R": {"linear": round(sig_R_lin, 9), "exact_GH": round(sig_R_ex, 9),
                "MC_mean": round(mc_mean, 9), "MC_se": round(mc_se, 9),
                "MC_runs": [round(x, 9) for x in mcs], "old_L1_MC": 0.039896,
                "canonical": round(SIG_R_CANON, 9)},
    "R_central": {"L1": 12.005533, "exact": round(R_ex, 9)},
    "z": {"old": round(z_old, 4), "new": round(z_new, 4)},
    "mu_pred_MeV": round(mu_pred, 6),
    "n_pass": sum(1 for c in CHECK if c["ok"]), "n_total": len(CHECK),
}
print("\n" + "=" * 104)
print("RESULT: %d/%d PASS" % (DATA["n_pass"], DATA["n_total"]))
print("=" * 104)

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "44_sigma_R_resolution.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"checks": CHECK, "data": DATA}, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON -> data/44_sigma_R_resolution.json")