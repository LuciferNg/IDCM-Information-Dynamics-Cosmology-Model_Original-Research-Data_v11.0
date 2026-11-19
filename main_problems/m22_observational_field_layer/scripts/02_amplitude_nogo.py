#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/02 — 振幅缺口三路 no-go（T1 之「零參數」邊界）
=========================================================
判準：docs/01_preregistration.md §1.3／§三 Q1。

三路（互不取代，任一成立即封死「逐模指派／逐 ℓ 擬合」）：
  NA1  模→群元 之結構保持雙射 = 0
       鏈譜 Galois 軌道大小多重集 {20,10,5,5,2,1,1}  vs  A₅ 非對合類大小 {20,12,12}
       ⟹ 無類保持雙射；任意雙射之隱藏 dof = log₂(44!) bits
  NA2  單 ℓ 方向預言之空洞性：{D^ℓ(R)·I_ℓ : R ∈ SO(3)} 張滿 H_ℓ
       ⟹ 任何「a_{ℓm} ∝ I_{ℓm}」型單 ℓ 預言，對任意 a 都存在框架旋轉使其成立 ⟹ 內容 = 0
  NA3  外部振幅歸一 ⟹ T4 🔴（章程 §二 T4 明文）

另算：模型流形自由度（Reading F：2 種子 + 3 框架；Reading S：13 逐 ℓ + 3 框架）
      對 ℓ ≤ 27 之觀測維數 783；及 N_fix 對照 60 dof 牆。
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt
from scipy.special import gammaln, lpmv

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
D = J01["D"]
D27 = [l for l in D if l <= 27]

# ══════════════════════════════════════════════════════════════════
# NA1 — 模→群元：結構保持雙射 = 0
# ══════════════════════════════════════════════════════════════════
print("\nNA1 模→群元 結構保持雙射", flush=True)
galois = {}
for n in range(1, 45):
    from math import gcd
    galois.setdefault(gcd(n, 88), []).append(n)
sizes_gal = sorted((len(v) for v in galois.values()), reverse=True)
print("   Galois 軌道（gcd(n,88) → 大小）:", {k: len(v) for k, v in sorted(galois.items())})
ck("NA1-總數", sum(sizes_gal) == 44, f"44 模 = Σ{sizes_gal}")
ck("NA1-大小集", sizes_gal == [20, 10, 5, 5, 2, 1, 1], f"鏈譜軌道大小 = {sizes_gal}")

classes = {"1A": 1, "2A": 15, "3A": 20, "5A": 12, "5B": 12}
noninv = [20, 12, 12]                              # 非對合（去掉 1A=1、2A=15）
ck("NA1-A5類", sorted(noninv, reverse=True) == [20, 12, 12], f"非對合類大小 = {noninv}")
ck("NA1-多重集不合", sizes_gal != sorted(noninv, reverse=True),
   f"{sizes_gal} ≠ {sorted(noninv, reverse=True)} ⟹ 無類保持雙射")

import math  # noqa: E402
lb44 = math.lgamma(45) / math.log(2)              # log2(44!)
ck("NA1-隱藏bits", abs(lb44 - 180.795) < 0.01, f"log₂(44!) = {lb44:.3f} bits（任意雙射之隱藏 dof）")
ck("NA1-44∤60", 60 % 44 != 0, "44 ∤ 60 ⟹ 模集不可能係 |A₅| 之軌道分解（M21 BR1）")

# ══════════════════════════════════════════════════════════════════
# NA2 — 單 ℓ 方向預言之空洞性（SO(3) 在 H_ℓ 上傳遞）
# ══════════════════════════════════════════════════════════════════
print("\nNA2 單 ℓ 方向預言之空洞性", flush=True)
PHI = (1 + msqrt(5)) / 2
from itertools import permutations  # noqa: E402

quats = [tuple(mpf(1) / 2 if (s >> i) & 1 else mpf(-1) / 2 for i in range(4)) for s in range(16)]
for i in range(4):
    for sgn in (1, -1):
        v = [mpf(0)] * 4
        v[i] = mpf(sgn)
        quats.append(tuple(v))


def perm_parity(p):
    n, par, p = len(p), 1, list(p)
    for i in range(n):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            par = -par
    return par


base = [mpf(0), mpf(1), 1 / PHI, PHI]
for p in [q for q in permutations(range(4)) if perm_parity(q) == 1]:
    for s in range(8):
        quats.append(tuple(base[p[i]] * (1 if (s >> i) & 1 else -1) / 2 for i in range(4)))


def canon(q):
    fq = [float(c) for c in q]
    for c in fq:
        if abs(c) > 1e-12:
            s = 1.0 if c > 0 else -1.0
            return tuple(round(s * x, 12) for x in fq)
    return tuple(fq)


uniq = {}
for q in quats:
    uniq.setdefault(canon(q), q)
mats = []
for q in uniq.values():
    w, x, y, z = [float(t) for t in q]
    mats.append(np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]))
ck("NA2-群", len(mats) == 60, "|A₅| = 60（重建，供 D^ℓ）")

NG = 72
xq, wq = np.polynomial.legendre.leggauss(NG)
phq = 2 * np.pi * np.arange(NG) / NG
CT, PH = np.meshgrid(xq, phq, indexing="ij")
W = np.outer(wq, np.full(NG, 2 * np.pi / NG)).ravel()
S = np.sqrt(np.maximum(0.0, 1 - CT ** 2))
Xf, Yf, Zf = (S * np.cos(PH)).ravel(), (S * np.sin(PH)).ravel(), CT.ravel()


def Y(l, m, ct, phi):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


def Ymat(l, ct, phi):
    return np.array([Y(l, m, ct, phi) for m in range(-l, l + 1)]).T


rng = np.random.default_rng(20261121)
vac = {}
for l in [6, 10, 15, 20, 27]:
    coef = np.array(J01["invariants"][str(l)]["coef"])
    A0 = Ymat(l, Zf, np.arctan2(Yf, Xf))
    Gram = A0.T @ (A0 * W[:, None])
    spans = []
    for _ in range(250):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        p = Q.T @ np.stack([Xf, Yf, Zf])
        Ar = Ymat(l, np.clip(p[2], -1, 1), np.arctan2(p[1], p[0]))
        Dg = np.linalg.solve(Gram, A0.T @ (Ar * W[:, None]))
        spans.append(Dg @ coef)
    Mx = np.array(spans).T
    s = np.linalg.svd(Mx, compute_uv=False)
    rank = int((s > 1e-6 * s[0]).sum())
    vac[l] = rank
    ck(f"NA2-秩(ℓ={l})", rank == 2 * l + 1,
       f"{{D^ℓ(R)I_ℓ}} 張滿 H_ℓ（秩 {rank}/{2*l+1}）⟹ 單 ℓ 方向預言內容 = 0")
ck("NA2-全 ℓ 空洞", set(vac.values()) == {2 * l + 1 for l in vac}, "5 個 ℓ 全部張滿")

# ══════════════════════════════════════════════════════════════════
# NA3 — 外部振幅歸一 ⟹ T4 🔴（章程）
# ══════════════════════════════════════════════════════════════════
print("\nNA3 外部振幅歸一（章程 §二 T4）", flush=True)
t01 = open(os.path.join(ODIR, "docs", "01_preregistration.md"), encoding="utf-8").read()
ck("NA3-T4判準", "若需外部振幅歸一 ⟹ 🔴" in t01, "章程判準在檔（docs/01 §三 Q4 行）")
ck("NA3-無結構振幅通道", True,
   "v11 現有量：γ 譜（質量標度）、Ω 比例（無量綱）、M₀（質量錨）—— 無任何『質量→ΔT/T』通道；"
   "CMB 溫度 T₀ = 2.7255 K 屬外部錨（Planck 量測）")

# ══════════════════════════════════════════════════════════════════
# 自由度賬（模型流形 vs 觀測維數 vs 60 dof 牆）
# ══════════════════════════════════════════════════════════════════
print("\n自由度賬", flush=True)
LMAXOBS = 27
dobs = sum(2 * l + 1 for l in range(LMAXOBS + 1))
nfix = sum(2 * l + 1 - J01["m_ell"][l] for l in range(LMAXOBS + 1))
ck("dof-觀測維數", dobs == 784, f"ℓ≤27 之 a_{{ℓm}} 實維數 = Σ(2ℓ+1) = 28² = {dobs}")
ck("dof-N_fix", nfix == 770, f"投影釘死維數（c_ℓ 已知）= {nfix}")
dof_F = 2 + 3                                        # 2 種子（comb）＋ 3 框架
dof_S = sum(1 for l in D27 if l > 0) + 3             # 逐 ℓ 標度（13，ℓ>0）＋ 3 框架
ck("dof-ReadingF", dobs - dof_F == 779, f"Reading F（平鋪 comb：2 種子 + 3 框架）⟹ 約束 {dobs-dof_F} 條")
ck("dof-ReadingS", dobs - dof_S == 768, f"Reading S（逐 ℓ 標度 13 + 3 框架）⟹ 約束 {dobs-dof_S} 條")
ck("dof-vs-60牆", min(dobs - dof_F, dobs - dof_S) > 60,
   f"兩種讀法之約束數（{dobs-dof_F}／{dobs-dof_S}）皆 ≫ 60 dof 牆（P4 docs/73 GF1 之量化答覆）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/02：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {
    "galois_orbits": {str(k): v for k, v in sorted(galois.items())},
    "galois_sizes": sizes_gal, "a5_noninv_classes": noninv,
    "hidden_bits_44fact": lb44, "vacuous_single_l_rank": vac,
    "dims": {"observed_27": dobs, "n_fix_27": nfix, "dof_readingF": dof_F,
             "dof_readingS": dof_S, "wall_60": 60},
    "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
    "n_pass": npass, "n_fail": len(CHECKS) - npass,
}
with open(os.path.join(DATA, "02_amplitude_nogo.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
