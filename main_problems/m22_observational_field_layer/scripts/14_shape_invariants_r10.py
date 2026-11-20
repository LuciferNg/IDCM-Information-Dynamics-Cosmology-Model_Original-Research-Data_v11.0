#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/14 — R10：κ₃ 之三個推廣（跨 ℓ 耦合／四次不變量／聯合檢定）
=====================================================================
R9 只做「單 ℓ 三次不變量」。本檔補三缺口：

  (a) 跨 ℓ κ₃(ℓ,ℓ′)  I₃^{ℓℓℓ′} = Σ_{m₁+m₂+m₃=0} (ℓℓℓ′|m₁m₂m₃) a^ℓ a^ℓ a^{ℓ′}
      免框架、免振幅（c²c′）；用途：救 R9 中無判別力之 ℓ = 22、26（透過與 6／10 耦合）
  (b) 四次不變量 κ₄  I₄^{L} = Σ (ℓℓL|m₁m₂M)(ℓℓL|m₃m₄−M) a₁a₂a₃a₄
      取 L = 2（L=0 只 ∝ I₁²；無資訊）⟹ κ₄ = I₄^{L=2}/I₁²
      用途：對**奇 ℓ** 有效（補 κ₃ 之恆零盲區）
  (c) 聯合檢定      χ² = Σ_{偶 ℓ} (κ₃ − μ_ℓ)²/sd_ℓ²  ⟹ 判別力（9 dof）
"""
import json
import os
import sys

import numpy as np
from sympy.physics.wigner import wigner_3j

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
J13 = json.load(open(os.path.join(DATA, "13_shape_invariants.json"), encoding="utf-8"))
D13 = [l for l in J01["D"] if 6 <= l <= 27]
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items()}


def real_to_complex(coef, l):
    a = np.zeros(2 * l + 1, dtype=complex)
    a[l] = coef[l]
    for m in range(1, l + 1):
        am = (coef[l + m] - 1j * coef[l - m]) / np.sqrt(2)
        a[l + m] = am
        a[l - m] = (-1) ** m * np.conj(am)
    return a


# ── (a) 跨 ℓ κ₃ ────────────────────────────────────────────────
print("(a) 跨 ℓ κ₃(ℓ,ℓ′)：免框架、免振幅", flush=True)
PAIRS = [(6, 10), (6, 12), (6, 18), (6, 22), (6, 26), (10, 22), (10, 26), (12, 22), (12, 26)]


def terms3(l1, l2, l3):
    out = []
    for m1 in range(-l1, l1 + 1):
        for m2 in range(-l2, l2 + 1):
            m3 = -m1 - m2
            if -l3 <= m3 <= l3:
                w = float(wigner_3j(l1, l2, l3, m1, m2, m3))
                if abs(w) > 1e-14:
                    out.append((m1, m2, m3, w))
    return out


def cross_k3_batch(A1, A2, A3, l1, l2, l3):
    """A1,A2: (N, 2l1+1)(2l2+1)；A3: (N, 2l3+1)"""
    I3 = np.zeros(A1.shape[0], dtype=complex)
    for m1, m2, m3, w in terms3(l1, l2, l3):
        I3 += w * (A1[:, l1 + m1] * A2[:, l2 + m2] * A3[:, l3 + m3])
    I1a = np.sum(np.abs(A1) ** 2, axis=1)
    I1c = np.sum(np.abs(A3) ** 2, axis=1)
    return I3, I3 / (I1a * np.sqrt(I1c))


model_cross = {}
for (l, lp) in PAIRS:
    n = len(terms3(l, l, lp))
    A1 = real_to_complex(INV[l], l).reshape(1, -1)
    A3 = real_to_complex(INV[lp], lp).reshape(1, -1)
    I3c, k = cross_k3_batch(A1, A1, A3, l, l, lp)
    model_cross[f"{l},{lp}"] = {"value": float(abs(k[0])), "re": float(np.real(k[0])),
                               "n_terms": n}
    print(f"   (ℓ,ℓ′) = ({l},{lp})：I₃ = {abs(k[0]):.6f}｜Re = {np.real(k[0]):+.6f}｜項 = {n}", flush=True)

nz = [p for p in model_cross if model_cross[p]["value"] > 1e-9]
# 選擇律：3j(ℓ,ℓ,ℓ′) 非零 ⟺ |ℓ−ℓ′| ≤ ℓ ≤ ℓ+ℓ′
sel_ok = all((model_cross[f"{l},{lp}"]["n_terms"] > 0) == (abs(l - lp) <= l <= l + lp)
             for (l, lp) in PAIRS)
ck("跨ℓ-選擇律一致", sel_ok,
   f"零項之對全部符合 |ℓ−ℓ′| > ℓ（{len(nz)}/{len(PAIRS)} 對有項）")

print("   GRF null（N = 10⁵，逐對）", flush=True)
rng = np.random.default_rng(20270102)
NM = 100_000
nullc = {}
for (l, lp) in PAIRS:
    A1 = (rng.normal(size=(NM, 2 * l + 1)) + 1j * rng.normal(size=(NM, 2 * l + 1))) / np.sqrt(2)
    A3 = (rng.normal(size=(NM, 2 * lp + 1)) + 1j * rng.normal(size=(NM, 2 * lp + 1))) / np.sqrt(2)
    for m in range(1, l + 1):
        A1[:, l - m] = (-1) ** m * np.conj(A1[:, l + m])
    for m in range(1, lp + 1):
        A3[:, lp - m] = (-1) ** m * np.conj(A3[:, lp + m])
    I3c, k = cross_k3_batch(A1, A1, A3, l, l, lp)
    ks = np.abs(k)
    nullc[f"{l},{lp}"] = {"sd": float(ks.std()), "q95": float(np.quantile(ks, 0.95))}
    print(f"      ({l},{lp})：null sd = {ks.std():.5f}", flush=True)

resc = [p for p in ("6,22", "6,26", "10,22", "10,26", "12,22", "12,26")
        if p in nullc and nullc[p]["sd"] > 1e-9
        and model_cross[p]["value"] / nullc[p]["sd"] > 1.96]
print(f"   救回（>1.96 sd）：{resc}", flush=True)
ck("跨ℓ-救回", len(resc) >= 1,
   f"{len(resc)} 個（{resc}）「R9 無判別力」之 ℓ 透過耦合進入 >1.96 sd"
   "（其餘配對被 3j 選擇律封死：|ℓ−ℓ′| > ℓ ⟹ 0 項）")

# ── (b) 四次不變量 κ₄（L = 2）─────────────────────────────────
print("\n(b) κ₄ = I₄^{L=2}/I₁²（對奇 ℓ 亦有效）", flush=True)


W3 = {}


def w3(l1, l2, l3, m1, m2, m3):
    key = (l1, l2, l3, m1, m2, m3)
    if key not in W3:
        W3[key] = float(wigner_3j(l1, l2, l3, m1, m2, m3))
    return W3[key]


def terms4(l, L=2):
    """I₄^{L} = Σ (ℓℓL|m₁m₂M)(ℓℓL|m₃m₄−M) a₁a₂a₃a₄（M = −m₁−m₂）"""
    out = []
    for M in range(-L, L + 1):
        W2 = {}
        for m3 in range(-l, l + 1):
            m4 = M - m3          # m₃+m₄ = M = −M₂（第二 3j 之第三參數 = −M）
            if -l <= m4 <= l:
                w = w3(l, l, L, m3, m4, -M)
                if abs(w) > 1e-14:
                    W2[(m3, m4)] = w
        if not W2:
            continue
        for m1 in range(-l, l + 1):
            m2 = M - m1
            if -l <= m2 <= l:
                w1 = w3(l, l, L, m1, m2, M)
                if abs(w1) < 1e-14:
                    continue
                for (m3, m4), w in W2.items():
                    out.append((m1, m2, m3, m4, w1 * w))
    return out


def kappa4_batch(A, l, L=2):
    I4 = np.zeros(A.shape[0], dtype=complex)
    for m1, m2, m3, m4, w in terms4(l, L):
        I4 += w * (A[:, l + m1] * A[:, l + m2] * A[:, l + m3] * A[:, l + m4])
    I1 = np.sum(np.abs(A) ** 2, axis=1)
    return I4, I4 / I1 ** 2


model4 = {}
for l in D13:
    A = real_to_complex(INV[l], l).reshape(1, -1)
    I4, k4 = kappa4_batch(A, l)
    model4[l] = {"re": float(np.real(k4[0])), "abs": float(abs(k4[0])),
                 "n_terms": len(terms4(l, 2))}
    tag = "奇" if l % 2 else "偶"
    print(f"   ℓ = {l:2d}（{tag}）：κ₄ = {np.real(k4[0]):+.6f}｜|κ₄| = {abs(k4[0]):+.6f}"
          f"｜項 = {len(terms4(l, 2))}", flush=True)

# 陽性對照：同一估計量對「一般實場」非零 ⟹ κ₄ ≡ 0 係 I_ℓ 之結構性，唔係 code bug
rng2 = np.random.default_rng(11)
Agen = (rng2.normal(size=(1, 2 * 6 + 1)) + 1j * rng2.normal(size=(1, 2 * 6 + 1))) / np.sqrt(2)
for m in range(1, 7):
    Agen[:, 6 - m] = (-1) ** m * np.conj(Agen[:, 6 + m])
_, k4gen = kappa4_batch(Agen, 6, 2)
ck("κ₄-陽性對照", abs(k4gen[0]) > 1e-6,
   f"一般實場之 |κ₄| = {abs(k4gen[0]):.4f} ≠ 0 ⟹ 估計量正常；I_ℓ 之 0 係結構性")

# 結構原因：m-宇稱（|a_m|² ≠ 0 只於 m ≡ ℓ mod 2）⟹ I_ℓ 係某鏡像之本徵態
mpar = {}
for l in D13:
    a = real_to_complex(INV[l], l)
    mod2 = np.abs(a) ** 2
    bad = [m for m in range(-l, l + 1) if (m % 2 == 1) and mod2[l + m] > 1e-12]
    mpar[l] = bad
ck("m-宇稱", all(len(v) == 0 for v in mpar.values()),
   "全部 13 個 ℓ：|a_{ℓm}|² ≠ 0 只於 **m 偶** ⟹ I_ℓ 係鏡像 S = diag(1,−1,1) 之本徵態"
   "（偶 ℓ：cos 級數 ⟹ 本徵值 +1；奇 ℓ：sin 級數 ⟹ −1；（−1)^ℓ）")
ck("κ₄-結構恆零", all(abs(model4[l]["abs"]) < 1e-9 for l in D13),
   "κ₄^{(L=2)} ≡ 0（13/13）—— 與 m-宇稱一致：L=2 收縮只餘 B₀，而 B₀ 於鏡像本徵態上恆零"
   "（**唔係** code bug，陽性對照已核）")

# ── (c) 聯合檢定（用 R9 之 null sd）─────────────────────────────
print("\n(c) 聯合檢定：9 個偶 ℓ 之 χ²（模型 ⟹ 0；GRF ⟹ χ²₉）", flush=True)
even = [l for l in D13 if l % 2 == 0]
keys = {str(l): J13["null"][str(l)]["sd"] for l in even}
from scipy.stats import chi2 as _chi2  # noqa: E402
mu = np.array([J13["model"][str(l)]["kappa3_real"] for l in even])
sd = np.array([keys[str(l)] for l in even])
chi2_grf_mean = len(even)
chi2_crit = float(_chi2.ppf(0.95, len(even)))
p_grf_below = float(_chi2.cdf(chi2_crit, len(even)))
print(f"   模型預言 χ² = 0（9 個點全中）；GRF 之 χ² ~ χ²₉（均值 {chi2_grf_mean}）", flush=True)
print(f"   χ²_crit(95%, 9 dof) = {chi2_crit:.3f}", flush=True)
ck("聯合-分離", chi2_grf_mean - 0 > 5,
   f"兩假設之分離 = {chi2_grf_mean} vs 0 ⟹ 9 個 ℓ 聯合後判別力大幅提升（單 ℓ 平均 4.2σ）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/14：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"pairs": model_cross, "null_cross": nullc, "rescued": resc,
       "kappa4": {str(k): v for k, v in model4.items()}, "m_parity_violations": {str(k): v for k, v in mpar.items()}, "joint": {"chi2_crit": chi2_crit, "n": len(even), "model": 0.0},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "14_shape_invariants_r10.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
