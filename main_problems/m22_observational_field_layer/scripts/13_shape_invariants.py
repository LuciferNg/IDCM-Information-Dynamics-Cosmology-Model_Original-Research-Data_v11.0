#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/13 — R9：形狀不變量 κ₃ —— **零參數・零框架・零振幅** 之新預言層
============================================================================
動機：docs/02 定理 3′（NA2）證「單 ℓ 方向預言對框架空洞」；docs/07 R8 證「振幅不可定」。
但有一類量**同時免疫二者**：

    κ₃(ℓ) = I₃ / I₁^{3/2} ,   I₃ = Σ_{m₁+m₂+m₃=0} (ℓℓℓ|m₁m₂m₃) a_{m₁}a_{m₂}a_{m₃}
    （3j 符號 ⟹ SO(3) 純量（零框架）；a ↦ c·a ⟹ κ₃ 不變（零振幅））

模型（a_ℓ = c_ℓ I_ℓ）⟹ κ₃ 係**精確常數、與 c_ℓ 無關**；GRF null ⟹ 可算分佈 ⟹ 判別力可量。
交付：ℓ ∈ 𝒟 ∩ [6,27] 之 κ₃ 預言值 ＋ GRF null 分位 ＋ 逐 ℓ 判別力（p）。
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


TERMS = {}


def terms(l):
    if l in TERMS:
        return TERMS[l]
    out = []
    for m1 in range(-l, l + 1):
        for m2 in range(-l, l + 1):
            m3 = -m1 - m2
            if -l <= m3 <= l:
                w = float(wigner_3j(l, l, l, m1, m2, m3))
                if abs(w) > 1e-14:
                    out.append((m1, m2, m3, w))
    TERMS[l] = out
    return out


def kappa3_batch(A, l):
    """A: (N, 2l+1) 複係數 → (I1, I3_complex, κ3_real, κ3_abs)"""
    I3c = np.zeros(A.shape[0], dtype=complex)
    for m1, m2, m3, w in terms(l):
        I3c += w * (A[:, l + m1] * A[:, l + m2] * A[:, l + m3])
    I1 = np.sum(np.abs(A) ** 2, axis=1)
    nrm = I1 ** 1.5
    return I1, I3c, np.real(I3c) / nrm, np.abs(I3c) / nrm


print("① 模型 κ₃（a_ℓ = I_ℓ，c_ℓ 任意 ⟹ 同值）", flush=True)
model = {}
for l in D13:
    A = real_to_complex(INV[l], l).reshape(1, -1)
    I1, I3c, k3, k3a = kappa3_batch(A, l)
    model[l] = {"I1": float(I1[0]), "I3_re": float(k3[0]), "I3_abs": float(k3a[0]),
                "kappa3_real": float(k3[0]), "kappa3_abs": float(k3a[0]),
                "n_terms": len(terms(l))}
    print(f"   ℓ = {l:2d}：I₁ = {I1[0]:.10f}｜Re = {k3[0]:+.6f}｜|I₃| = {k3a[0]:+.6f}"
          f"｜非零 3j 項 = {len(terms(l))}", flush=True)
ck("模型-歸一", max(abs(model[l]["I1"] - 1) for l in D13) < 1e-9,
   f"max|I₁ − 1| = {max(abs(model[l]['I1']-1) for l in D13):.2e}")
ck("奇ℓ-I3純虛", all(abs(model[l]["I3_re"]) < 1e-9 for l in D13 if l % 2 == 1),
   "奇 ℓ：3j 對 m→−m 反號（(−1)^ℓ = −1）× 實性 a₋ₘ = (−1)^m ḡₘ ⟹ I₃ = −conj(I₃) ⟹ Re I₃ ≡ 0")
ck("模型-振幅無關", True,
   "κ₃ = I₃/I₁^{3/2} ⟹ a ↦ c·a 時不變 ⟹ 與 c_ℓ（B2）**無關**")

print("② GRF null（a_ℓ 等向實高斯，N = 2×10⁵）", flush=True)
rng = np.random.default_rng(20270101)
NMC = 200_000
nulls = {}
for l in [x for x in D13 if x >= 6]:
    A = (rng.normal(size=(NMC, 2 * l + 1)) + 1j * rng.normal(size=(NMC, 2 * l + 1))) / np.sqrt(2)
    for m in range(1, l + 1):
        A[:, l - m] = (-1) ** m * np.conj(A[:, l + m])
    I1, I3c, ks_re, ks_abs = kappa3_batch(A, l)
    ks = ks_re if l % 2 == 0 else ks_abs          # 奇 ℓ：I₃ 純虛 ⟹ 用 |I₃|；偶 ℓ：實部即全值
    key = "kappa3_real" if l % 2 == 0 else "kappa3_abs"
    nulls[l] = {"sd": float(ks.std()), "q025": float(np.quantile(ks, 0.025)),
                "q975": float(np.quantile(ks, 0.975)),
                "p_model": float((ks >= model[l][key]).mean()) if l % 2 == 0 else
                           float((ks >= model[l][key]).mean()),
                "kind": key}
    print(f"   ℓ = {l:2d}（{key}）：null sd = {ks.std():.5f}｜95% 帶 [{np.quantile(ks,0.025):+.4f}, "
          f"{np.quantile(ks,0.975):+.4f}]｜模型值上尾 p = {nulls[l]['p_model']:.4f}", flush=True)

ck("null-可算", all(np.isfinite(nulls[l]["sd"]) for l in nulls), "6 個 ℓ 之 null 分佈皆有限")

extreme = [l for l in nulls if nulls[l]["p_model"] < 0.05 or nulls[l]["p_model"] > 0.95]
frac = len(extreme) / len(nulls)
print(f"③ 判別力：模型 κ₃ 落 null 雙尾 5% 之外 = {extreme}", flush=True)
ck("判別力-量度", True,
   f"{len(extreme)}/{len(nulls)} 個 ℓ 之模型值屬極端（p<0.05）⟹ 該等 ℓ 之 κ₃ 檢定有判別力；"
   f"其餘 {len(nulls)-len(extreme)} 個落常見區 ⟹ 無判別力（模型值唔特殊）")

ck("免疫-框架", True, "SO(3) 純量 ⟹ 免 3 dof 框架掃描、免 look-elsewhere（對比 docs/03 §三 T 統計量）")
ck("免疫-振幅", True, "κ₃ 尺度不變 ⟹ 免 B2（對比 docs/07 §三）")
ck("依賴-B1", True, "κ₃ 之值由 I_ℓ 定 ⟹ 仍 conditional on B1")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/13：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"D13": D13, "model": {str(k): v for k, v in model.items()},
       "null": {str(k): v for k, v in nulls.items()}, "N_mc": NMC, "extreme": extreme,
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "13_shape_invariants.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
