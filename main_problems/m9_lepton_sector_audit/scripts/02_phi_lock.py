#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/02 — R2 φ-lock 閘（P3）

(A) K_φ 凍結集 → 逐指數判成員資格（事前指定 vs 反解）
(B) 格點殘差指紋：必需分子 b* = D·frac(x*)；觀測殘差 / 格距
(C) 族密度 E：{a + b/D}（D ∈ 88/132/528）→ 半格值容差、隨機命中率
(D) MC null 10^7：log-uniform 目標 → P(|rel dev| ≤ τ) 實測
(E) 上游 primary-source 審計（唯讀）：o34 scripts/21 之構造次序
(F) 判決

零 import 上游；零寫入他 M/O 目錄。mp.dps = 50；MC 用 numpy。
"""
import hashlib
import json
import os
import re

import numpy as np
from mpmath import mp, mpf, log, sqrt

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

PHI = (1 + sqrt(5)) / 2
LN_PHI = log(PHI)
R_MUE = mpf("206.7682827")      # CODATA 2022 direct entry（scripts/01 已核）
U_MUE = mpf("0.0000046")
R_TAUMU = mpf("16.8170")
R_TAUE = mpf("3477.23")

checks = []


def chk(name, cond, detail=""):
    checks.append({"check": name, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


# ════════════════════════════════════════════════════════════════════
# (A) K_φ 凍結集（O38 scripts/33 慣例；事前指定）
#     K_φ = {−2,1,2,3,4,5,7,11,18,29,47,76,123,199,322,521}
#           ∪ {291.5 − ℓ | ℓ = 0..33}
# ════════════════════════════════════════════════════════════════════
LUCAS_CLOSURE = [-2, 1, 2, 3, 4, 5, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521]
K_TREE = [mpf(291) + mpf(1) / 2 - ell for ell in range(34)]
K_PHI = [mpf(v) for v in LUCAS_CLOSURE] + K_TREE

print("=" * 74)
print("M9 scripts/02 — R2 φ-lock 閘")
print("=" * 74)
print(f"|K_φ| = {len(K_PHI)}  (Lucas 閉合 16 ＋ K_tree μ_ℓ 層 34)")

X_CAND = {
    "C1  μ/e   φ^(11+7/88)": 11 + mpf(7) / 88,
    "L34 μ/τ   φ^(6−71/528)": 6 - mpf(71) / 528,
    "L34 e/τ   φ^(17−29/528)": 17 - mpf(29) / 528,
}
X_OBS = {
    "C1  μ/e   φ^(11+7/88)": log(R_MUE) / LN_PHI,
    "L34 μ/τ   φ^(6−71/528)": log(R_TAUMU) / LN_PHI,
    "L34 e/τ   φ^(17−29/528)": log(R_TAUE) / LN_PHI,
}

print("\n[A] K_φ 成員資格（指數必須事前落在凍結集）")
kphi_rows = []
for k, x in X_CAND.items():
    inbox = any(abs(x - kp) < mpf("1e-30") for kp in K_PHI)
    nearest = min(K_PHI, key=lambda kp: abs(x - kp))
    print(f"  {k:26s} x = {mp.nstr(x, 12)}   ∈K_φ = {inbox}   nearest K_φ = "
          f"{mp.nstr(nearest, 8)} (Δ = {mp.nstr(abs(x - nearest), 4)})")
    kphi_rows.append({"candidate": k, "x": mp.nstr(x, 12), "in_K_phi": bool(inbox),
                      "nearest_K_phi": mp.nstr(nearest, 8),
                      "delta_to_nearest": mp.nstr(abs(x - nearest), 4)})
chk("A1 11+7/88 ∉ K_φ", not kphi_rows[0]["in_K_phi"])
chk("A2 6−71/528 ∉ K_φ", not kphi_rows[1]["in_K_phi"])
chk("A3 17−29/528 ∉ K_φ", not kphi_rows[2]["in_K_phi"])

# ════════════════════════════════════════════════════════════════════
# (B) 格點殘差指紋
#     若分子係事前結構整數 → 殘差應 ~U[0, 1/2] 格；殘差 ≪ 格 = 先取 log 後凑
# ════════════════════════════════════════════════════════════════════
print("\n[B] 格點殘差指紋（b* = D·frac(x_obs)，D ∈ {88,132,528}）")
grid_rows = []
for d in (88, 132, 528):
    print(f"  ── D = {d}  (格距 Δx = 1/{d} = {mp.nstr(mpf(1) / d, 6)}；"
          f"值容差半格 = {mp.nstr((PHI ** (mpf(1) / (2 * d)) - 1) * 100, 5)} %)")
    for k in X_CAND:
        xo = X_OBS[k]
        b_star = d * (xo - mp.floor(xo))
        resid = abs(b_star - round(b_star))
        cell_frac = resid / mpf(1)          # 以格為單位（b 單位 = 1 格）
        p_rand = min(mpf(1), 2 * cell_frac)
        print(f"     {k:26s} b* = {mp.nstr(b_star, 8):>12s}  → 取整 {round(b_star):3d}"
              f"  殘差 = {mp.nstr(resid, 4)} 格 ({mp.nstr(cell_frac * 100, 3)}% 格)"
              f"  P_rand = {mp.nstr(p_rand, 4)}")
        grid_rows.append({"D": d, "candidate": k, "b_star": mp.nstr(b_star, 8),
                          "b_int": int(round(b_star)), "resid_cell": mp.nstr(resid, 6),
                          "p_rand": mp.nstr(p_rand, 6)})

# ════════════════════════════════════════════════════════════════════
# (C) 族密度 E：{a + b/D}，a ∈ 1..30，b ∈ 0..D−1
#     三分母之並集 = 全 1/528 格（因 1/88 = 6/528、1/132 = 4/528）
# ════════════════════════════════════════════════════════════════════
print("\n[C] 族密度 E（隨機目標命中率）")
dens_rows = []
for d in (88, 132, 528):
    cell_x = 1.0 / d
    cell_rel = float(PHI ** (mpf(1) / d) - 1)      # 相鄰格之相對值差
    half_rel = float(PHI ** (mpf(1) / (2 * d)) - 1)
    row = {"D": d, "n_points": 30 * d, "cell_rel": cell_rel, "half_cell_rel": half_rel}
    for tau in (2.2e-8, 1e-5, 1e-4, 1e-3):
        row[f"p_(|dev|<={tau:g})"] = min(1.0, 2 * tau / cell_rel)
    dens_rows.append(row)
    print(f"  D={d:4d}  點數 {30 * d:6d}  格距(值) {cell_rel * 100:8.4f}%  "
          f"半格 {half_rel * 100:8.4f}%   P(1e-5)={row['p_(|dev|<=1e-05)']:.4f}  "
          f"P(1e-4)={row['p_(|dev|<=0.0001)']:.4f}  P(σ=2.2e-8)={row['p_(|dev|<=2.2e-08)']:.6f}")

# ════════════════════════════════════════════════════════════════════
# (D) MC null 10^7（log-uniform 目標；seed 固定）
# ════════════════════════════════════════════════════════════════════
print("\n[D] MC null 10^7（x = log_φ(target) ~ U[1,30]）")
NMC = 10_000_000
rng = np.random.default_rng(20261027)
x_mc = rng.uniform(1.0, 30.0, NMC)
mc_rows = []
for d in (88, 132, 528):
    for tau in (2.2e-8, 1e-5, 1e-4):
        delta_x = float(np.log1p(tau) / float(LN_PHI))
        res = np.abs(x_mc * d - np.round(x_mc * d)) / d     # 距最近格（x 單位）
        hit = float(np.mean(res <= delta_x))
        analytic = min(1.0, 2 * tau / float(PHI ** (mpf(1) / d) - 1))
        mc_rows.append({"D": d, "tau": tau, "p_mc": hit, "p_analytic": analytic})
        print(f"  D={d:4d} τ={tau:<9g}  P_MC = {hit:.6f}   P_analytic = {analytic:.6f}")
chk("D1 528 格於 τ=1e-5 之命中率 >1% (無鑑別力)",
    [r for r in mc_rows if r["D"] == 528 and r["tau"] == 1e-5][0]["p_mc"] > 0.01,
    f"{[r for r in mc_rows if r['D'] == 528 and r['tau'] == 1e-5][0]['p_mc']:.4f}")
chk("D2 528 格於 σ=2.2e-8 之命中率 <1% (σ 級全滅)",
    [r for r in mc_rows if r["D"] == 528 and r["tau"] == 2.2e-8][0]["p_mc"] < 0.01,
    f"{[r for r in mc_rows if r['D'] == 528 and r['tau'] == 2.2e-8][0]['p_mc']:.6f}")

# ════════════════════════════════════════════════════════════════════
# (E) 上游 primary-source 審計（唯讀）
# ════════════════════════════════════════════════════════════════════
print("\n[E] 上游 primary-source 審計：O34 scripts/21_generation_lucas.py（唯讀）")
UP = os.path.expanduser(
    "~/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/open_problems_fresh/"
    "o34_observer_introspection/scripts/21_generation_lucas.py")
UP_DOC = os.path.expanduser(
    "~/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/open_problems_fresh/"
    "o35_beats_mass_linkup/docs/02_transfer_matrix_closure.md")
audit = {"path": UP}
with open(UP, "rb") as f:
    raw = f.read()
audit["sha256"] = hashlib.sha256(raw).hexdigest()
src = raw.decode("utf-8")

# E1：先算目標 k12，再硬編 (11, 7)
e1 = re.search(r'k12 = log_phi\(m_mu / m_e\)', src) is not None
e2 = re.search(r'"k\(e→μ\)": \(k12, 11, 7, "L₅ \+ L₄/88"\)', src) is not None
# E2：顯式 (a,b) 網格搜尋 1800 點
e3 = ("for a in range(1, 31):" in src) and ("for b in range(1, 61):" in src)
# E3：反解打印
e4 = "71 = 528×(6-k23)" in src
# E4：損壞目標 m_e = 0.511
e5 = re.search(r"^m_e\s*=\s*0\.511\s*$", src, re.M) is not None
print(f"  sha256 = {audit['sha256']}")
print(f"  E1 先算目標 k12 = log_φ(m_μ/m_e)                    : {e1}")
print(f"  E2 之後硬編 (a=11, b=7) 標籤 'L₅ + L₄/88'           : {e2}")
print(f"  E3 顯式 (a∈1..30, b∈1..60) 網格搜尋（1800 點反搜）   : {e3}")
print(f"  E4 反解打印 '71 = 528×(6-k23)'                      : {e4}")
print(f"  E5 目標用 m_e = 0.511（非 0.51099895）              : {e5}")

# E5 後果：用 m_e = 0.511 重現文件記載之 0.0003% / 0.0015%
m_e_bad = mpf("0.511")
m_mu_doc = mpf("105.658")
k12_bad = log(m_mu_doc / m_e_bad) / LN_PHI
cand = 11 + mpf(7) / 88
err_k_pct = abs(cand - k12_bad) / k12_bad * 100
val_bad = m_mu_doc / m_e_bad
err_val_pct = abs(PHI ** cand / val_bad - 1) * 100
print(f"  → 用 m_e=0.511 重現：k 空間偏差 = {mp.nstr(err_k_pct, 6)} %（文件記 0.0003%）")
print(f"  → 用 m_e=0.511 重現：值偏差     = {mp.nstr(err_val_pct, 6)} %（文件記 0.0015%）")
k12_true = log(R_MUE) / LN_PHI
err_k_true = abs(cand - k12_true) / k12_true * 100
err_val_true = abs(PHI ** cand / R_MUE - 1) * 100
print(f"  → 權威值：k 空間 {mp.nstr(err_k_true, 6)} %，值 {mp.nstr(err_val_true, 6)} %")

audit.update({"E1_target_first": bool(e1), "E2_hardcode_after": bool(e2),
              "E3_grid_search_1800": bool(e3), "E4_reverse_print": bool(e4),
              "E5_corrupt_m_e": bool(e5),
              "repro_k_pct_with_0511": mp.nstr(err_k_pct, 6),
              "repro_val_pct_with_0511": mp.nstr(err_val_pct, 6),
              "true_k_pct": mp.nstr(err_k_true, 6),
              "true_val_pct": mp.nstr(err_val_true, 6)})

chk("E1 上游先算目標後硬編指數（後驗指紋）", e1 and e2)
chk("E2 上游含 1800 點 (a,b) 網格反搜", e3)
chk("E3 上游含反解打印 71 = 528×(6−k23)", e4)
chk("E4 上游目標用 m_e = 0.511（損壞 ~2e−6）", e5)
chk("E5 記載 0.0003%（k 空間）可由 m_e=0.511 重現",
    abs(err_k_pct - mpf("0.0003")) < mpf("0.00008"), mp.nstr(err_k_pct, 6))
chk("E6 記載 0.0015%（值空間）可由 m_e=0.511 重現",
    abs(err_val_pct - mpf("0.0015")) < mpf("0.0004"), mp.nstr(err_val_pct, 6))

# E7：O35 「分子代數」表（Step 2）缺 7
with open(UP_DOC, "r", encoding="utf-8") as f:
    doc = f.read()
num_table = doc[doc.find("### Step 2"):doc.find("### Step 3")]
has7 = bool(re.search(r"\|\s*\*?\*?7\*?\*?\s*\|", num_table))
print(f"  E7 O35 docs/02 Step-2「分子代數」表含 7 作為 μ/e 分子 : {has7}")
audit["E7_numerator_7_absent_in_O35_table"] = not has7
chk("E7 7 未出現於 O35 分子代數表（7 無上游推導）", not has7)

# ════════════════════════════════════════════════════════════════════
# (F) 判決
# ════════════════════════════════════════════════════════════════════
sigma_C1 = abs(PHI ** (11 + mpf(7) / 88) - R_MUE) / U_MUE
verdict = {
    "K_phi_membership": "FAIL（三個候選指數全部 ∉ K_φ）",
    "pre_specification": "FAIL（E1/E2/E3/E4：先算目標、後硬編、1800 點網格反搜、反解打印）",
    "numeric": f"C1 值偏差 {mp.nstr(sigma_C1, 6)}σ → 精確身份已否證",
    "resolution_mismatch": (
        f"ladder 最細格（528）半格 = {mp.nstr((PHI ** (mpf(1) / (2 * 528)) - 1) * 100, 5)}%"
        f" vs 觀測相對精度 {mp.nstr(U_MUE / R_MUE * 100, 6)}%"
        f" → 粗 {mp.nstr((PHI ** (mpf(1) / (2 * 528)) - 1) / (U_MUE / R_MUE), 4)}×"),
    "verdict": "🔴 REFUTED（後驗 fit ＋ 數值否證；K_φ FAIL）",
}
print("\n[F] 判決")
for k, v in verdict.items():
    print(f"  {k:22s}: {v}")

out = os.path.join(DATA, "02_phi_lock.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump({"k_phi": {"lucas": LUCAS_CLOSURE,
                         "k_tree": [str(x) for x in K_TREE],
                         "size": len(K_PHI)},
               "membership": kphi_rows, "grid_residual": grid_rows,
               "density": dens_rows, "mc": mc_rows, "upstream_audit": audit,
               "verdict": verdict, "checks": checks,
               "pass": sum(1 for c in checks if c["pass"]), "total": len(checks)},
              f, ensure_ascii=False, indent=2)
print("\n" + "=" * 74)
print(f"PASS {sum(1 for c in checks if c['pass'])}/{len(checks)}   -> {out}")
print("=" * 74)