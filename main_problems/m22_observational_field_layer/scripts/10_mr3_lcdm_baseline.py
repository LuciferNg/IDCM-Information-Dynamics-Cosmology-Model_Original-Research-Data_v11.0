#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/10 — MR3 攻堅：ΛCDM 實際譜作基準（CAMB，Planck 2018 best-fit）
=========================================================================
MR3 原文：`ratio16` 之對照基準依賴「填回 ℓ≤5 之尺度不變延拓」。
本檔改用**真 ΛCDM 譜**（CAMB 2.0.1，Planck 2018 TT,TE,EE+lowE+lensing best-fit）：

  (a) ΛCDM 之 σ²_16 缺損 D = 1 − σ²_16(ℓ=6..47)/σ²_16(ℓ=2..47)（healpy pixwin(16)）
      = 模型「ℓ≤5 恰零」在 Nside=16 解析度上嘅等價缺損（ν 一個純 ΛCDM 數，M17 可直接用）
  (b) 用**真 ΛCDM 權重**重算 f 之 GRF null（10⁷）⟹ 修正 docs/03 之 SW／平 C_ℓ 代理門檻
  (c) ΛCDM 之 ratio16 = σ²_16(ℓ≤47)/σ²_total(ℓ≤2500) 作對照尺度
"""
import json
import os
import sys

import camb
import healpy as hp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ── (0) CAMB：Planck 2018 best-fit ─────────────────────────────
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, omk=0, tau=0.0544)
pars.InitPower.set_params(As=np.exp(3.044e0) / 1e10, ns=0.9649)
pars.set_for_lmax(2500, lens_potential_accuracy=2)
pars.Want_CMB_lensing = True
res = camb.get_results(pars)
tot = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)
ClTT = tot["total"][:, 0]
ell = np.arange(ClTT.size)
print(f"   CAMB ok：ℓ_max = {ClTT.size-1}；C_2 = {ClTT[2]:.1f} μK²；C_10 = {ClTT[10]:.1f} μK²", flush=True)
ck("MR3-CAMB", ClTT.size > 100 and ClTT[2] > 100,
   f"CAMB ΛCDM 譜（Planck 2018 best-fit）：C₂ = {ClTT[2]:.1f}、C₁₀ = {ClTT[10]:.1f} μK²")

# ── (a) σ²_16 缺損 ────────────────────────────────────────────
W16 = hp.pixwin(16)
LM = 47


def s16(lo, hi):
    return float(sum((2 * l + 1) / (4 * np.pi) * ClTT[l] * W16[l] ** 2 for l in range(lo, hi + 1)))


full = s16(2, LM)
drop = 1.0 - s16(6, LM) / full
band = s16(2, 5) / full
print(f"   ΛCDM：σ²_16(ℓ2..47) = {full:.3f} μK² ｜ ℓ≤5 佔比 = {band*100:.2f}% ｜ 缺損 D = {drop*100:.2f}%",
      flush=True)
ck("MR3-缺損", 0.05 < drop < 0.6,
   f"ΛCDM 之 Nside-16 variance 缺損（ℓ≤5 歸零）= {drop*100:.2f}% ⟹ 模型預言之等價幅度")

# 不同 Nside 之穩健性
rows = {}
for ns in (8, 16, 32):
    W = hp.pixwin(ns)
    lm = min(3 * ns - 1, ClTT.size - 1)
    f_ = float(sum((2 * l + 1) / (4 * np.pi) * ClTT[l] * W[l] ** 2 for l in range(2, lm + 1)))
    d_ = 1.0 - float(sum((2 * l + 1) / (4 * np.pi) * ClTT[l] * W[l] ** 2 for l in range(6, lm + 1))) / f_
    rows[ns] = (f_, d_)
    print(f"   Nside={ns:2d}（ℓ≤{lm}）：σ² = {f_:.3f} ｜ 缺損 = {d_*100:.2f}%", flush=True)
ck("MR3-Nside穩健", rows[16][1] > rows[32][1] > 0,
   "缺損隨解析度單調減（低解析度更受 ℓ≤5 支配）")

# ── (b) 真 ΛCDM 權重之 f null（10⁷）──────────────────────────
print("\n(b) 真 ΛCDM 權重之 f null（N = 10⁷）", flush=True)
TOT = list(range(2, 28))
LOW = list(range(2, 6))
w = {l: ClTT[l] for l in TOT}
NMC = 10_000_000
rng = np.random.default_rng(20261201)
out = np.empty(NMC)
done, chunk = 0, 2_000_000
while done < NMC:
    n = min(chunk, NMC - done)
    X = np.zeros(n)
    Y = np.zeros(n)
    for l in TOT:
        chi = rng.chisquare(2 * l + 1, size=n)
        if l in LOW:
            X += w[l] * chi
        Y += w[l] * chi
    out[done:done + n] = X / Y
    done += n
q = lambda p: float(np.quantile(out, p))
print(f"   f(LCDM權重): mean = {out.mean():.4f} ｜ q0.1% = {q(0.001):.5f} ｜ q1% = {q(0.01):.5f} "
      f"｜ q5% = {q(0.05):.5f}", flush=True)
a_bonf = 0.05 / 22
fc_b, fc_5 = q(a_bonf), q(5.7e-7)
print(f"   f_crit(Bonf α={a_bonf:.2e}) = {fc_b:.5f} ｜ f_crit(5σ) = {fc_5:.5f}", flush=True)
ck("MR3-f門檻", fc_5 < fc_b < 1.0, f"真 ΛCDM 權重門檻：{fc_b:.5f}（Bonf）／{fc_5:.5f}（5σ）")

# ── (c) 幅度之 convention 依賴：comb 譜 vs ΛCDM 譜 ─────────────
print("\n(c) 「ℓ≤5 歸零」缺損之幅度：comb 譜 vs ΛCDM 譜", flush=True)
J03 = json.load(open(os.path.join(DATA, "03_comb_field_stats.json"), encoding="utf-8"))
supp = J03["bands"]["suppression_sigma16"]
print(f"   模型 comb 譜（Reading F，24 種子）：缺損 = {supp[0]*100:.2f}% … {supp[1]*100:.2f}%", flush=True)
print(f"   ΛCDM 譜（本檔）：缺損 = {drop*100:.2f}%", flush=True)
spread = (drop*100) / (supp[1]*100)
ck("MR3-幅度依賴", spread > 3.0,
   f"缺損幅度跨 {supp[0]*100:.2f}%–{drop*100:.2f}%（ΛCDM/comb 比 = {spread:.1f}×）"
   " ⟹ 幅度**不可預言**（只零位結構可）；docs/03 之「帶窄」係 comb 藍譜之 artifact")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/10：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
outj = {"C2": float(ClTT[2]), "C10": float(ClTT[10]),
        "sigma16_full": full, "drop_low5": drop, "band_low5": band,
        "per_nside": {str(k): [v[0], v[1]] for k, v in rows.items()},
        "f_null_lcdm": {"mean": float(out.mean()), "q001": q(0.001), "q01": q(0.01),
                        "q05": q(0.05), "f_crit_bonf": fc_b, "f_crit_5sigma": fc_5},
        "suppression_model_comb": supp, "deficit_spread_ratio": spread, "N_mc": NMC,
        "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
        "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "10_mr3_lcdm_baseline.json"), "w", encoding="utf-8") as f:
    json.dump(outj, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
