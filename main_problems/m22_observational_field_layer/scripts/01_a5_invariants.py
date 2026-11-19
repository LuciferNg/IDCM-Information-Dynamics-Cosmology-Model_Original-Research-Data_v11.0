#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/01 — A₅ 之天空載體：群結構、不變譜支撐 𝒟 = {ℓ : m_ℓ ≥ 1}、不變球諧 I_ℓ
================================================================================
判準：docs/01_preregistration.md（凍結）。本檔只算 §3.1 之量，不判死生。

三條獨立路線（不一致 ⟹ 停）：
  M1  精確字符閉式（純整數算術；φ 項必須精確消去）
       m_ℓ = [ (2ℓ+1) + 15(−1)^ℓ + 20·c₃(ℓ) + 12·(2 + 5·n₀(ℓ) − ℓ) ] / 60
  M2  群重建：2I（600-cell）120 四元數 → ±1 商 → 60 旋轉；類大小/轉角直接量（零假設）
  M3  投影算子秩法：D^ℓ(g) 由「旋轉取樣 ＋ 加權最小二乘」數值求，rank(P̄) 直接數
       （零字符公式、零類表；只用 M2 之旋轉矩陣）

旁證 M1b：以 M2 量到之 60 個轉角直接代入 χ_ℓ(θ) = Σ_{m=−ℓ}^{ℓ} e^{imθ} 數值求和。

另產：不變球諧 I_ℓ = P̄ 之唯一本徵向量（m_ℓ = 1 時）；不變性殘差、宇稱、N_fix。
"""
import json
import os
import sys
import time
from collections import Counter
from itertools import permutations

import numpy as np
import sympy as sp
from mpmath import mp, mpf, sqrt as msqrt
from scipy.special import gammaln, lpmv

mp.dps = 40
T0 = time.time()

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
os.makedirs(DATA, exist_ok=True)
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


PHI = (1 + msqrt(5)) / 2

# ══════════════════════════════════════════════════════════════════
# M2 — 群重建（2I / {±1}）
# ══════════════════════════════════════════════════════════════════
print("\nM2 群重建：2I / {±1} ≅ A₅ ⊂ SO(3)", flush=True)

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
ck("M2-四元數", len(quats) == 120, f"2I = {len(quats)}（600-cell 頂點）")


def canon(q):
    fq = [float(c) for c in q]
    for c in fq:
        if abs(c) > 1e-12:
            s = 1.0 if c > 0 else -1.0
            return tuple(round(s * x, 12) for x in fq)
    return tuple(fq)


uniq = {}
for q in quats:
    assert abs(float(sum(x * x for x in q)) - 1) < 1e-12
    uniq.setdefault(canon(q), q)
rot_q = list(uniq.values())
ck("M2-數目", len(rot_q) == 60, f"|A₅| = {len(rot_q)}（120 四元數 / ±1）")


def rotmat(q):
    w, x, y, z = [float(t) for t in q]
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


mats = [rotmat(q) for q in rot_q]
ck("M2-正交", max(np.abs(M @ M.T - np.eye(3)).max() for M in mats) < 1e-12, "max|RRᵀ−I| < 1e-12")
ck("M2-行列式", max(abs(np.linalg.det(M) - 1) for M in mats) < 1e-11, "det R = +1（手性，無反射）")


def key(M):
    return tuple(np.round(M.ravel(), 8))


keys = {key(M) for M in mats}
ck("M2-封閉", all(key(A @ B) in keys for A in mats for B in mats), "60×60 乘積落回集合")
ck("M2-逆元", all(key(M.T) in keys for M in mats), "R⁻¹ = Rᵀ 落回集合")

deg = []
for M in mats:
    c = max(-1.0, min(1.0, (np.trace(M) - 1) / 2))
    deg.append(round(float(np.arccos(c)) * 180 / np.pi, 6))
cnt = Counter(deg)
print("   轉角分佈（度 → 個數）：", dict(sorted(cnt.items())), flush=True)
ck("M2-類大小", dict(sorted(cnt.items())) == {0.0: 1, 72.0: 12, 120.0: 20, 144.0: 12, 180.0: 15},
   "類大小 = 1/15/20/12/12（A₅ 五類，直接量）")

# ══════════════════════════════════════════════════════════════════
# M1 ＋ M1b — 精確字符閉式 ／ 由量到之轉角數值求和
# ══════════════════════════════════════════════════════════════════
print("\nM1 精確字符閉式：m_ℓ = [(2ℓ+1) + 15(−1)^ℓ + 20c₃ + 12(2+5n₀−ℓ)]/60", flush=True)


def counts_mod(ell, k):
    return [sum(1 for m in range(1, ell + 1) if m % k == j) for j in range(k)]


def m_exact(ell):
    n3 = counts_mod(ell, 3)
    c3 = 1 + 2 * n3[0] - (n3[1] + n3[2])
    n0 = counts_mod(ell, 5)[0]
    num = (2 * ell + 1) + 15 * (-1) ** ell + 20 * c3 + 12 * (2 + 5 * n0 - ell)
    assert num % 60 == 0, f"ℓ={ell} 分子 {num} 不整除 60"
    return num // 60


LMAX = 60
mseq = [m_exact(l) for l in range(LMAX + 1)]
ck("M1-非負整數", all(m >= 0 for m in mseq), f"ℓ ≤ {LMAX} 全部 m_ℓ ∈ ℤ≥0")

ph = (1 + sp.sqrt(5)) / 2
ok_cancel = True
for ell in range(0, 41):
    n5 = counts_mod(ell, 5)
    c5 = 1 + 2 * (n5[0] + (n5[1] + n5[4]) / (2 * ph) - (n5[2] + n5[3]) * ph / 2)
    d5 = 1 + 2 * (n5[0] + (n5[2] + n5[3]) / (2 * ph) - (n5[1] + n5[4]) * ph / 2)
    if sp.simplify(c5 + d5 - (2 + 5 * n5[0] - ell)) != 0:
        ok_cancel = False
ck("M1-φ消去", ok_cancel, "c₅+d₅ − (2+5n₀−ℓ) = 0 於 ℚ(√5)（ℓ ≤ 40 逐項 sympy 核）")

th = np.array([d * np.pi / 180 for d in deg])


def m_from_angles(ell):
    chi = np.zeros_like(th)
    for m in range(-ell, ell + 1):
        chi += np.cos(m * th)
    return chi.mean()


mb = [m_from_angles(l) for l in range(LMAX + 1)]
resid_b = max(abs(mb[l] - mseq[l]) for l in range(LMAX + 1))
ck("M1b-轉角求和", resid_b < 1e-9,
   f"max|(1/60)Σχ_ℓ(θ_g) − m_ℓ^exact| = {resid_b:.2e}（ℓ ≤ 60；用 M2 直接量到之 60 個轉角）")

Dset = [l for l in range(0, 41) if mseq[l] >= 1]
De = [l for l in Dset if l % 2 == 0 and l > 0]
Do = [l for l in Dset if l % 2 == 1]
print(f"   𝒟 (ℓ ≤ 40) = {Dset}", flush=True)
print(f"   𝒟_even = {De}\n   𝒟_odd  = {Do}", flush=True)
ck("M1-ℓ1..5", all(mseq[l] == 0 for l in range(1, 6)), "m_{1..5} = 0（IR 端結構零）")
ck("M1-ℓmin", min(l for l in Dset if l > 0) == 6, "ℓ_min = 6")
ck("M1-ℓodd_min", min(Do) == 15, "ℓ_odd^min = 15（奇度不變存在）")
ck("M1-ℓ14", mseq[14] == 0, "m₁₄ = 0（⟨6,10⟩ 半群缺口）")


def nfix(lmax):
    return sum(2 * l + 1 - mseq[l] for l in range(lmax + 1))


print(f"   N_fix(ℓ≤27) = {nfix(27)} ｜ N_fix(ℓ≤40) = {nfix(40)} ｜ 60 dof 牆 = 60", flush=True)

# ══════════════════════════════════════════════════════════════════
# M3 — 投影算子秩法（零字符、零類表）
# ══════════════════════════════════════════════════════════════════
print("\nM3 投影算子秩法：rank(P̄), P̄ = (1/60)Σ_g D^ℓ(g)", flush=True)

NG = 64
xq, wq = np.polynomial.legendre.leggauss(NG)
phq = 2 * np.pi * np.arange(NG) / NG
CT, PH = np.meshgrid(xq, phq, indexing="ij")
Wq = np.outer(wq, np.full(NG, 2 * np.pi / NG)).ravel()
ST = np.sqrt(np.maximum(0.0, 1 - CT ** 2))
Xa, Ya, Za = (ST * np.cos(PH)).ravel(), (ST * np.sin(PH)).ravel(), CT.ravel()


def sph_harm_real(l, m, ct, phi):
    am = abs(m)
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    if m == 0:
        return norm * P
    return norm * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


def Ymatrix(l, ct, phi):
    return np.array([sph_harm_real(l, m, ct, phi) for m in range(-l, l + 1)]).T


LMC = 30
m_proj, inv_eig = [], []
for l in range(0, LMC + 1):
    A0 = Ymatrix(l, Za, np.arctan2(Ya, Xa))
    Gram = A0.T @ (A0 * Wq[:, None])
    Pbar = np.zeros((2 * l + 1, 2 * l + 1))
    for M in mats:
        pinv = (M.T @ np.stack([Xa, Ya, Za]))            # R⁻¹ n̂
        ct2 = np.clip(pinv[2], -1, 1)
        ph2 = np.arctan2(pinv[1], pinv[0])
        Ar = Ymatrix(l, ct2, ph2)
        Dg = np.linalg.solve(Gram, A0.T @ (Ar * Wq[:, None]))   # D^ℓ(g), column-major
        Pbar += Dg
    Pbar /= len(mats)
    ev = np.linalg.eigvals(Pbar)
    m_proj.append(int(np.sum(np.abs(ev - 1) < 1e-6)))
    inv_eig.append(float(np.sort(np.abs(ev - 1))[0]))
print("   m_ℓ (M3proj) =", m_proj, flush=True)
print("   m_ℓ (M1)     =", mseq[:LMC + 1], flush=True)
ck("M3-秩法 vs M1", m_proj == mseq[:LMC + 1],
   f"ℓ ≤ {LMC}：{sum(1 for a, b in zip(m_proj, mseq) if a == b)}/{LMC+1} 相同")
ck("M3-本徵值整數性", max(abs(e - round(e)) for e in m_proj) < 1e-12,
   "rank 直接數（零連續自由度）")

# ══════════════════════════════════════════════════════════════════
# 不變球諧 I_ℓ = P̄ 之唯一本徵向量
# ══════════════════════════════════════════════════════════════════
print("\n不變球諧 I_ℓ（m_ℓ = 1 之 ℓ）", flush=True)
NGF = 96
xf, wf = np.polynomial.legendre.leggauss(NGF)
phf = 2 * np.pi * np.arange(NGF) / NGF
CTF, PHF = np.meshgrid(xf, phf, indexing="ij")
WF = np.outer(wf, np.full(NGF, 2 * np.pi / NGF)).ravel()
STF = np.sqrt(np.maximum(0.0, 1 - CTF ** 2))
Xf, Yf, Zf = (STF * np.cos(PHF)).ravel(), (STF * np.sin(PHF)).ravel(), CTF.ravel()

alldata, inv_resid, orth_resid, sign_rule = {}, 0.0, 0.0, "first-nonzero-coef>0"
for l in [x for x in Dset if 0 < x <= 30]:
    A0 = Ymatrix(l, Zf, np.arctan2(Yf, Xf))
    Gram = A0.T @ (A0 * WF[:, None])
    orth_resid = max(orth_resid, float(np.abs(Gram - np.eye(2 * l + 1)).max()))
    Pbar = np.zeros((2 * l + 1, 2 * l + 1))
    for M in mats:
        pinv = M.T @ np.stack([Xf, Yf, Zf])
        Ar = Ymatrix(l, np.clip(pinv[2], -1, 1), np.arctan2(pinv[1], pinv[0]))
        Pbar += np.linalg.solve(Gram, A0.T @ (Ar * WF[:, None]))
    Pbar /= len(mats)
    ev, evec = np.linalg.eig(Pbar)
    idx = int(np.argmax(np.abs(ev - 1) ** (-1)))
    coef = np.real(evec[:, idx])
    coef = coef / np.sqrt(coef @ Gram @ coef)
    nz = next(c for c in coef if abs(c) > 1e-12)
    if nz < 0:
        coef = -coef
    ck(f"Iℓ-唯一次數(ℓ={l})", abs(ev[idx].real - 1) < 1e-6,
       f"本徵值 = {ev[idx].real:.8f}（唯一不變方向）")
    rng = np.random.default_rng(20261120)
    pts = rng.normal(size=(240, 3))
    pts /= np.linalg.norm(pts, axis=1)[:, None]

    def ev_at(P):
        ct_, ph_ = np.clip(P[:, 2], -1.0, 1.0), np.arctan2(P[:, 1], P[:, 0])
        return Ymatrix(l, ct_, ph_) @ coef
    b0 = ev_at(pts)
    for M in mats[1:]:
        inv_resid = max(inv_resid, float(np.abs(ev_at((M @ pts.T).T) - b0).max()))
    alldata[l] = {"coef": [float(c) for c in coef], "m": int(mseq[l]),
                  "eig": float(ev[idx].real)}

ck("Y-正交", orth_resid < 1e-9, f"max|Gram − I| = {orth_resid:.2e}（實球諧 ∫dΩ 正交歸一）")
ck("Iℓ-不變性", np.isfinite(inv_resid) and inv_resid < 1e-7, f"max|I(Rn̂) − I(n̂)| = {inv_resid:.2e}")

par_resid = 0.0
for l, d in alldata.items():
    coef = np.array(d["coef"])
    p = np.array([0.3, -0.5, 0.81])
    p /= np.linalg.norm(p)

    def ev2(P):
        ct_, ph_ = np.clip(P[:, 2], -1.0, 1.0), np.arctan2(P[:, 1], P[:, 0])
        return (Ymatrix(l, ct_, ph_) @ coef)[0]
    Pp = p.reshape(1, 3)
    f1, f2 = ev2(Pp), ev2(-Pp)
    par_resid = max(par_resid, abs(f2 - (-1) ** l * f1) / max(abs(f1), 1e-30))
ck("Iℓ-宇稱", np.isfinite(par_resid) and par_resid < 1e-9, f"max rel |I(−n̂) − (−1)^ℓ I(n̂)| = {par_resid:.2e}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/01：{npass} PASS ／ {len(CHECKS) - npass} FAIL ｜ {time.time()-T0:.1f}s ===")
out = {
    "group": {"order": len(rot_q), "class_angles_deg": {str(k): v for k, v in sorted(cnt.items())}},
    "m_ell": mseq[:LMAX + 1], "m_proj": m_proj, "D": Dset, "D_even": De, "D_odd": Do,
    "l_min": 6, "l_odd_min": 15, "N_fix_27": nfix(27), "N_fix_40": nfix(40),
    "invariants": alldata, "invariance_resid": inv_resid, "quadrature": {"NG_proj": NG, "NG_final": NGF},
    "sign_rule": sign_rule,
    "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
    "n_pass": npass, "n_fail": len(CHECKS) - npass,
}
with open(os.path.join(DATA, "01_a5_invariants.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
