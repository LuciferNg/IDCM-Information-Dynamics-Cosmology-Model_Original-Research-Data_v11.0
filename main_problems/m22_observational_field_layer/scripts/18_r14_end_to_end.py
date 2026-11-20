#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/18 — R14：通道 2 之端到端演練（HEALPix 離散化系統誤差 ＋ null 校準 ＋ 提取步驟）
=============================================================================================
R13 證 κ₃ 可用「單 ℓ 地圖立方積分」實作（免 3j）。本檔補最後一格（R13 §四 之 🟡）：
  (1) HEALPix 像素立方積分之離散化誤差 vs Nside（ℓ = 6, 12, 26）
  (2) null 校準：像素估計器之 κ₃ 分佈 vs R9 之解析 sd（200 實現，ℓ = 6）
  (3) 提取步驟：多 ℓ 地圖 → 投影出單 ℓ → κ₃ 復原（M17 真實管線必做之一步）
全程**自建 Y 基底於像素中心**（不經 healpy 之 alm 公約）⟹ 只測離散化，不混公約。
"""
import json
import os
import sys

import numpy as np
from scipy.special import gammaln, lpmv

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


try:
    import healpy as hp
except Exception as e:                                     # pragma: no cover
    print("healpy 不可用：" + str(e))
    sys.exit(1)

J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
J13 = json.load(open(os.path.join(DATA, "13_shape_invariants.json"), encoding="utf-8"))
J17 = json.load(open(os.path.join(DATA, "17_r13_estimator_robustness.json"), encoding="utf-8"))
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items()}


def norm_lm(l, m):
    am = abs(m)
    return np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))


def Yrow(l, m, ct, ph):
    am = abs(m)
    n = norm_lm(l, m) * (np.sqrt(2) if m != 0 else 1.0)
    P = lpmv(am, l, ct)
    return n * P * (1.0 if m == 0 else (np.cos(am * ph) if m > 0 else np.sin(am * ph)))


def pixels(nside):
    npix = hp.nside2npix(nside)
    th, ph = hp.pix2ang(nside, np.arange(npix))   # healpy 回 θ（弧度），非 cos θ
    ct = np.clip(np.cos(th), -1 + 1e-12, 1 - 1e-12)   # lpmv 於 |x| = 1 會出 NaN（極點）
    return ct, ph, npix


def fmap(coef, l, ct, ph):
    out = np.zeros_like(ct)
    for m in range(-l, l + 1):
        c = coef[l + m]
        if abs(c) > 1e-16:
            out += c * Yrow(l, m, ct, ph)
    return out


print("(1) 像素立方積分之離散化誤差（vs NG=200 求積之收斂值）", flush=True)
exact = {}
for l in J17["ratios"]:
    v = J17["ratios"][l]
    if v["quad"] is not None:
        exact[int(l)] = v["quad"]
errs = {}
for l in [6, 12, 26]:
    print(f"   ℓ = {l}：精確 ∫f³ = {exact[l]:+.9f}", flush=True)
    for ns in (64, 128, 256, 512):
        ct, ph, npix = pixels(ns)
        f = fmap(INV[l], l, ct, ph)
        q = float(np.sum(f ** 3) * (4 * np.pi / npix))
        e = abs(q - exact[l]) / abs(exact[l])
        errs[f"{l},{ns}"] = e
        print(f"      Nside = {ns:4d}（npix = {npix:8d}）：∫f³ = {q:+.9f}｜相對誤差 = {e:.2e}", flush=True)
best = {l: min(errs[f"{l},{ns}"] for ns in (64, 128, 256, 512)) for l in (6, 12, 26)}
ck("R14-收斂", all(best[l] < 1e-2 for l in best),
   "最細誤差（Nside=512）：" + "｜".join(f"ℓ={l} → {best[l]:.2e}" for l in best))
mono = all(errs[f"{l},64"] > errs[f"{l},512"] for l in (6, 12, 26))
ck("R14-誤差隨 Nside 下降", mono, "Nside 64 → 512 之誤差單調下降（3 個 ℓ）")

print("\n(2) null 校準：像素估計器 vs R9 解析 sd（ℓ = 6，200 實現，Nside = 256）", flush=True)
NS_NULL = 256
ct, ph, npix = pixels(NS_NULL)
Omega = 4 * np.pi / npix
l = 6
# 由求積反推 κ₃ 之常數 C_ℓ（R13 恆等式：∫f³ = C_ℓ·I₃）
C6 = J17["ratios"]["6"]["pred"]
rng = np.random.default_rng(20270111)
ks = []
Nsamp = 200
bases = np.array([Yrow(l, m, ct, ph) for m in range(-l, l + 1)])   # (13, npix)
for i in range(Nsamp):
    re = rng.normal(size=2 * l + 1)
    im = rng.normal(size=2 * l + 1)
    a = (re + 1j * im) / np.sqrt(2)
    for m in range(1, l + 1):
        a[l - m] = (-1) ** m * np.conj(a[l + m])
    # 實場：f = Σ_{m>0} √2 Re(a_m) cos + ... —— 直接用 a 之實/虛部構成實係數
    coef = np.zeros(2 * l + 1)
    coef[l] = np.real(a[l])
    for m in range(1, l + 1):
        coef[l + m] = np.sqrt(2) * np.real(a[l + m])
        coef[l - m] = -np.sqrt(2) * np.imag(a[l + m])
    f = fmap(coef, l, ct, ph)
    I3_pix = float(np.sum(f ** 3) * Omega) / C6
    I1 = float(np.sum(coef ** 2))
    ks.append(I3_pix / I1 ** 1.5)
ks = np.array(ks)
sd_pix = float(ks.std())
sd_ana = float(J13["null"]["6"]["sd"])
print(f"   像素估計器：sd = {sd_pix:.5f}（均值 {ks.mean():+.5f}）", flush=True)
print(f"   R9 解析：   sd = {sd_ana:.5f}", flush=True)
print(f"   比 = {sd_pix/sd_ana:.4f}（MC 誤差 ±{1/np.sqrt(2*Nsamp-2):.4f}）", flush=True)
ck("R14-null校準", abs(sd_pix / sd_ana - 1) < 0.15,
   f"sd 比 = {sd_pix/sd_ana:.4f}（200 實現之抽樣誤差 ~7%）⟹ 像素估計器之 null 校準合格")

print("\n(3) 提取步驟：多 ℓ 地圖 → 投影出單 ℓ → κ₃ 復原", flush=True)
ct, ph, npix = pixels(128)
Omega = 4 * np.pi / npix
ls = [6, 12]
cs = {6: 1.0, 12: 0.37}
f_multi = np.zeros(npix)
for l in ls:
    f_multi += cs[l] * fmap(INV[l], l, ct, ph)
recovered = {}
for l in ls:
    B = np.array([Yrow(l, m, ct, ph) for m in range(-l, l + 1)])   # (2ℓ+1, npix)
    G = B @ B.T * Omega
    a = np.linalg.solve(G, B @ (f_multi * Omega))
    C = J17["ratios"][str(l)]["pred"]
    f_est = B.T @ a
    I3 = float(np.sum(f_est ** 3) * Omega) / C
    I1 = float(np.sum(a ** 2))
    k3 = I3 / I1 ** 1.5
    k3_exact = J13["model"][str(l)]["kappa3_real"]
    recovered[l] = {"k3": k3, "k3_exact": k3_exact, "dev": abs(k3 - k3_exact)}
    print(f"   ℓ = {l:2d}（多 ℓ 混合，c = {cs[l]}）：κ₃ 復原 = {k3:+.6f}｜"
          f"精確 = {k3_exact:+.6f}｜差 = {abs(k3-k3_exact):.2e}", flush=True)
ck("R14-提取復原", all(v["dev"] < 5e-3 for v in recovered.values()),
   "多 ℓ 混合下提取單 ℓ 後 κ₃ 復原：" +
   "｜".join(f"ℓ={l} 差 {v['dev']:.1e}" for l, v in recovered.items()))

print("\n(4) 交付：M17 建議參數", flush=True)
ck("R14-建議", best[6] < 1e-3 and best[12] < 1e-3,
   f"Nside = 512 時 ℓ=6/12 之離散化誤差 < 1e−3 ⟹ 建議 Nside ≥ 512（ℓ ≤ 26）；"
   f"ℓ=26 之誤差 {best[26]:.2e}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/18：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"discretization": errs, "best": {str(k): v for k, v in best.items()},
       "null_pix_sd": sd_pix, "null_ana_sd": sd_ana, "null_ratio": sd_pix / sd_ana,
       "extraction": {str(k): v for k, v in recovered.items()},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "18_r14_end_to_end.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
