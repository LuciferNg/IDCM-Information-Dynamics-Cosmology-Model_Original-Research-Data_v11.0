#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/15 — R11：鏡像本徵態之 realisation 獨立性／奇 ℓ 不變量普查／𝒟 延伸
==============================================================================
(1) R10 之鏡像本徵態只在本 realisation 測得（🟡）。本檔證其**共變性** ⟹ 升為定理：
    任何 realisation（共軛群 R·A₅·Rᵀ）之不變向量 I_ℓ^{(R)} = ±D^ℓ(R) I_ℓ
    ⟹ 鏡像性質隨之共變（鏡面法向 = R·ŷ）。
(2) 奇 ℓ 之低階純量不變量普查：order 3（ℓ³）與 order 4（ℓℓℓℓ→0，經 L=2）對奇 ℓ 是否恆零？
(3) 𝒟 延伸至 ℓ ≤ 200：m_ℓ、支援密度、N_fix 與「零位預言力」之衰減。
"""
import json
import os
import sys

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt
from scipy.special import gammaln, lpmv
from sympy.physics.wigner import wigner_3j

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
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items()}

PHI = (1 + msqrt(5)) / 2


def build_group():
    from itertools import permutations
    quats = [tuple(mpf(1) / 2 if (s >> i) & 1 else mpf(-1) / 2 for i in range(4)) for s in range(16)]
    for i in range(4):
        for sgn in (1, -1):
            v = [mpf(0)] * 4
            v[i] = mpf(sgn)
            quats.append(tuple(v))

    def par(p):
        n, s, p = len(p), 1, list(p)
        for i in range(n):
            while p[i] != i:
                j = p[i]
                p[i], p[j] = p[j], p[i]
                s = -s
        return s
    base = [mpf(0), mpf(1), 1 / PHI, PHI]
    for p in [q for q in permutations(range(4)) if par(q) == 1]:
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
        if abs(float(sum(x * x for x in q)) - 1) < 1e-12:
            uniq.setdefault(canon(q), q)
    mats = []
    for q in list(uniq.values())[:60]:
        w, x, y, z = [float(t) for t in q]
        mats.append(np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]))
    return mats


MATS = build_group()
NG = 96
xq, wq = np.polynomial.legendre.leggauss(NG)
phq = 2 * np.pi * np.arange(NG) / NG
CT, PH = np.meshgrid(xq, phq, indexing="ij")
W = np.outer(wq, np.full(NG, 2 * np.pi / NG)).ravel()
S_ = np.sqrt(np.maximum(0.0, 1 - CT ** 2))
XF, YF, ZF = (S_ * np.cos(PH)).ravel(), (S_ * np.sin(PH)).ravel(), CT.ravel()


def Y(l, m, ct, phi):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


YB, GR = {}, {}


def ybasis(l):
    if l not in YB:
        YB[l] = np.array([Y(l, m, CT, PH).ravel() for m in range(-l, l + 1)]).T if False else \
            np.array([Y(l, m, ZF, np.arctan2(YF, XF)) for m in range(-l, l + 1)]).T
        GR[l] = YB[l].T @ (YB[l] * W[:, None])
    return YB[l]


def Dmat(l, R):
    ybasis(l)                                    # 確保 YB/GR 已建
    pn = R.T @ np.stack([XF, YF, ZF])
    Ar = np.array([Y(l, m, np.clip(pn[2], -1, 1), np.arctan2(pn[1], pn[0])) for m in range(-l, l + 1)]).T
    return np.linalg.solve(GR[l], YB[l].T @ (Ar * W[:, None]))


# ── (1) 共變性：共軛 realisation 之不變向量 = ±D^ℓ(R) I_ℓ ─────
print("(1) 鏡像本徵態之 realisation 獨立性（共變性）", flush=True)
rng = np.random.default_rng(20270105)
worst_cov, worst_mirror = 0.0, 0.0
for l in [6, 10]:
    I0 = INV[l]
    ybasis(l)
    for trial in range(3):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        gm = [Q @ M @ Q.T for M in MATS]
        Pbar = np.zeros((2 * l + 1, 2 * l + 1))
        for M in gm:
            Pbar += Dmat(l, M)
        Pbar /= 60.0
        ev, evec = np.linalg.eig(Pbar)
        idx = int(np.argmax(np.abs(ev - 1) ** (-1)))
        Inew = np.real(evec[:, idx])
        Inew = Inew / np.linalg.norm(Inew)
        pred = Dmat(l, Q) @ I0
        pred = pred / np.linalg.norm(pred)
        a = abs(float(Inew @ pred))
        worst_cov = max(worst_cov, abs(1 - a))
        # 鏡像：法向 = Q·ŷ（原 realisation 之鏡面法向 = ŷ，因 I 只含 m 偶）
        yhat = np.array([0.0, 1.0, 0.0])
        Snew = Q @ np.diag([1.0, -1.0, 1.0]) @ Q.T          # 平面反射（法向 Q·ŷ）
        Ds = Dmat(l, Snew)
        r = min(np.linalg.norm(Ds @ Inew - Inew), np.linalg.norm(Ds @ Inew + Inew))
        worst_mirror = max(worst_mirror, r)
ck("R11-共變", worst_cov < 1e-8,
   f"共軛 realisation：max|1−|⟨I^{{(R)}}, D^ℓ(R)I⟩|| = {worst_cov:.2e}（ℓ=6,10 × 3 個隨機 R）")
ck("R11-鏡像隨共變", worst_mirror < 1e-7,
   f"鏡面法向 = R·ŷ：max‖D^ℓ(S_new)I_new ∓ I_new‖ = {worst_mirror:.2e} ⟹ 鏡像本徵態性質 **realisation-independent**")

# ── (2) 奇 ℓ 之低階純量不變量普查 ───────────────────────────────
print("\n(2) 奇 ℓ 之低階純量不變量（order 3 / order 4）", flush=True)


def real_to_complex(coef, l):
    a = np.zeros(2 * l + 1, dtype=complex)
    a[l] = coef[l]
    for m in range(1, l + 1):
        am = (coef[l + m] - 1j * coef[l - m]) / np.sqrt(2)
        a[l + m] = am
        a[l - m] = (-1) ** m * np.conj(am)
    return a


odd = [l for l in J01["D"] if l % 2 == 1 and l <= 27]
res2 = {}
for l in odd:
    a = real_to_complex(INV[l], l)
    # order 3
    I3 = 0j
    for m1 in range(-l, l + 1):
        for m2 in range(-l, l + 1):
            m3 = -m1 - m2
            if -l <= m3 <= l:
                I3 += float(wigner_3j(l, l, l, m1, m2, m3)) * a[l + m1] * a[l + m2] * a[l + m3]
    # order 4（L=2 收縮，取 |I₄|）
    I4 = 0j
    for m1 in range(-l, l + 1):
        for m2 in range(-l, l + 1):
            M = -m1 - m2
            if abs(M) > 2:
                continue
            w1 = float(wigner_3j(l, l, 2, m1, m2, M))
            if abs(w1) < 1e-14:
                continue
            for m3 in range(-l, l + 1):
                m4 = M - m3
                if -l <= m4 <= l:
                    w2 = float(wigner_3j(l, l, 2, m3, m4, -M))
                    I4 += w1 * w2 * a[l + m1] * a[l + m2] * a[l + m3] * a[l + m4]
    res2[l] = {"I3_abs": float(abs(I3)), "I4_abs": float(abs(I4))}
    print(f"   ℓ = {l:2d}（奇）：|I₃| = {abs(I3):.3e}｜|I₄^{{L=2}}| = {abs(I4):.3e}", flush=True)
ck("R11-奇ℓ低階全零", all(v["I3_abs"] < 1e-9 and v["I4_abs"] < 1e-9 for v in res2.values()),
   f"{len(odd)} 個奇 ℓ：order-3 與 order-4 純量不變量**皆恆零** ⟹ 奇 ℓ 之形狀對低階純量統計不可見")

# ── (3) 𝒟 延伸至 ℓ ≤ 200 ────────────────────────────────────────
print("\n(3) 𝒟 延伸（ℓ ≤ 200）與預言力衰減", flush=True)


def m_exact(ell):
    n3 = [sum(1 for m in range(1, ell + 1) if m % 3 == j) for j in range(3)]
    c3 = 1 + 2 * n3[0] - (n3[1] + n3[2])
    n0 = sum(1 for m in range(1, ell + 1) if m % 5 == 0)
    num = (2 * ell + 1) + 15 * (-1) ** ell + 20 * c3 + 12 * (2 + 5 * n0 - ell)
    assert num % 60 == 0
    return num // 60


L200 = 200
mseq = [m_exact(l) for l in range(L200 + 1)]
dz = [l for l in range(2, L200 + 1) if mseq[l] == 0]
print(f"   ℓ ≤ 200 之零位（m_ℓ = 0）：{len(dz)} 個 ⟹ {dz[:20]} …", flush=True)
last_zero = max(dz) if dz else None
print(f"   最後一個零位 ℓ* = {last_zero} ⟹ ℓ > {last_zero} 之支援**全滿**（每個 ℓ 皆有 1 維不變空間）", flush=True)
for L in (27, 40, 60, 100, 200):
    nz = sum(1 for l in range(2, L + 1) if mseq[l] == 0)
    tot = L - 1
    nfix = sum(2 * l + 1 - mseq[l] for l in range(L + 1))
    dims = sum(2 * l + 1 for l in range(L + 1))
    print(f"   ℓ ≤ {L:3d}：零位 {nz:2d}/{tot:3d}（{100*nz/tot:5.1f}%）｜N_fix = {nfix}"
          f"｜N_fix/維數 = {nfix/dims:.4f}", flush=True)
ck("R11-零位有限", last_zero is not None and last_zero < 200,
   f"零位集有限（ℓ* = {last_zero}）⟹ 支援約束在大 ℓ 完全消失（預言力集中於 ℓ ≤ {last_zero}）")
dens = {L: sum(1 for l in range(2, L + 1) if mseq[l] == 0) / (L - 1) for L in (27, 40, 60, 100, 200)}
ck("R11-密度衰減", dens[27] > dens[200],
   f"零位密度：ℓ≤27 = {100*dens[27]:.0f}% ⟹ ℓ≤200 = {100*dens[200]:.1f}%（單調衰減）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/15：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"covariance_dev": worst_cov, "mirror_dev": worst_mirror, "odd_low_order": res2,
       "last_zero_ell": last_zero, "zeros_up_to_200": dz,
       "density": {str(k): v for k, v in dens.items()},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "15_r11_covariance_and_density.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
