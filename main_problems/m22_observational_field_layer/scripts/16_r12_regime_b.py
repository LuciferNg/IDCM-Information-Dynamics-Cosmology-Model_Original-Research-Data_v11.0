#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/16 — R12：Regime B（ℓ ≥ 30）—— 不變子空間秩之獨立核 ＋ κ₃ 由「點」退化為「區間」
============================================================================================
R11 §3.2 只定性講「ℓ ≥ 30 形狀唔唯一」。本檔量化：

  (1) 獨立核 m_ℓ 公式於 ℓ = 30..40（用投影算子秩，不依賴 docs/02 之字符閉式）
  (2) 對 m_ℓ ≥ 2 之不變子空間，κ₃ 之可達集 = 區間 [min, max]
      （κ₃ 係 SO(3) 純量 ⟹ 對任意 realisation 同值 ⟹ **仍免框架**）
  (3) 區間寬度 vs GRF null sd ⟹ 剩餘判別力（預言力由「點」退化為「帶」之定量）
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
PHI_F = np.arctan2(YF, XF)
YB, GR = {}, {}


def Y(l, m):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ZF)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * PHI_F) if m > 0 else np.sin(am * PHI_F))


def ybasis(l):
    if l not in YB:
        YB[l] = np.array([Y(l, m) for m in range(-l, l + 1)]).T
        GR[l] = YB[l].T @ (YB[l] * W[:, None])
    return YB[l]


def Dmat(l, R):
    ybasis(l)
    pn = R.T @ np.stack([XF, YF, ZF])
    cols = []
    for m in range(-l, l + 1):
        am = abs(m)
        nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
        P = lpmv(am, l, np.clip(pn[2], -1, 1))
        ph = np.arctan2(pn[1], pn[0])
        cols.append(nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * ph) if m > 0 else np.sin(am * ph)))
    Ar = np.array(cols).T
    return np.linalg.solve(GR[l], YB[l].T @ (Ar * W[:, None]))


def m_formula(ell):
    n3 = [sum(1 for m in range(1, ell + 1) if m % 3 == j) for j in range(3)]
    c3 = 1 + 2 * n3[0] - (n3[1] + n3[2])
    n0 = sum(1 for m in range(1, ell + 1) if m % 5 == 0)
    return ((2 * ell + 1) + 15 * (-1) ** ell + 20 * c3 + 12 * (2 + 5 * n0 - ell)) // 60


print("(1) 不變子空間之秩：投影算子法 vs m_ℓ 公式（ℓ = 30..40）", flush=True)
ranks, subspaces = {}, {}
for l in range(30, 41):
    Pbar = np.zeros((2 * l + 1, 2 * l + 1))
    for M in MATS:
        Pbar += Dmat(l, M)
    Pbar /= 60.0
    ev, evec = np.linalg.eigh((Pbar + Pbar.T) / 2)
    idx = [i for i in range(len(ev)) if abs(ev[i] - 1) < 1e-8]
    ranks[l] = len(idx)
    subspaces[l] = np.real(evec[:, idx])
    print(f"   ℓ = {l:2d}：秩 = {len(idx)}｜m_ℓ 公式 = {m_formula(l)}"
          f"｜{'✅' if len(idx) == m_formula(l) else '❌'}", flush=True)
ck("R12-秩=公式", all(ranks[l] == m_formula(l) for l in range(30, 41)),
   f"ℓ = 30..40：投影算子秩 (= {sorted(set(ranks.values()))}) 全等於 m_ℓ 公式 ⟹ **獨立核**")


def kappa3_vec(a, l, terms):
    """偶 ℓ：Re(I₃)/I₁^{3/2}（I₃ 實）；奇 ℓ：|I₃|/I₁^{3/2}（I₃ 純虛 ⟹ 實部恆零）"""
    I3 = 0j
    for m1, m2, m3, w in terms:
        I3 += w * a[l + m1] * a[l + m2] * a[l + m3]
    I1 = float(np.sum(np.abs(a) ** 2))
    num = float(np.real(I3)) if l % 2 == 0 else float(abs(I3))
    return num / I1 ** 1.5


TERMS = {}


def terms(l):
    if l not in TERMS:
        out = []
        for m1 in range(-l, l + 1):
            for m2 in range(-l, l + 1):
                m3 = -m1 - m2
                if -l <= m3 <= l:
                    w = float(wigner_3j(l, l, l, m1, m2, m3))
                    if abs(w) > 1e-14:
                        out.append((m1, m2, m3, w))
        TERMS[l] = out
    return TERMS[l]


print("\n(2) κ₃ 之可達集（不變子空間內，m_ℓ ≥ 2 ⟹ 由點退化為區間）", flush=True)

J13 = json.load(open(os.path.join(DATA, "13_shape_invariants.json"), encoding="utf-8"))
rng = np.random.default_rng(20270107)
intervals = {}
for l in [l for l in range(30, 41) if m_formula(l) >= 2]:
    Q = subspaces[l]
    mm = Q.shape[1]
    Ns = 4000 if mm <= 3 else 2000
    coef = rng.normal(size=(Ns, mm))
    coef /= np.linalg.norm(coef, axis=1, keepdims=True)
    vals = np.array([kappa3_vec(Q @ coef[i], l, terms(l)) for i in range(Ns)])
    intervals[l] = {"m": int(mm), "min": float(vals.min()), "max": float(vals.max()),
                    "width": float(vals.max() - vals.min())}
    print(f"   ℓ = {l:2d}（m = {mm}）：κ₃ ∈ [{vals.min():+.5f}, {vals.max():+.5f}]"
          f"｜寬 = {vals.max()-vals.min():.5f}", flush=True)

deg = [l for l in range(30, 41) if m_formula(l) >= 2]
pt = [l for l in range(30, 41) if m_formula(l) == 1]
ck("R12-區間-退化判準",
   all(intervals[l]["width"] > 0 for l in intervals) and len(intervals) == len(deg),
   f"m_ℓ ≥ 2 之 ℓ = {deg}（{len(deg)} 個）⟹ κ₃ 區間非退化（寬 "
   + ", ".join(f"{intervals[l]['width']:.5f}" for l in intervals) + "）"
   f"；m_ℓ = 1 之 ℓ（{len(pt)} 個）⟹ 點預言（寬 = 0，如 ℓ=31/33）")

print("\n(3) 退化後之剩餘判別力（對比 R9 之單點 null sd）", flush=True)
for l in [30, 31, 32, 33]:
    sd = J13["null"].get(str(l), {}).get("sd")
    if sd is None:
        print(f"   ℓ = {l:2d}：R9 未算 null（ℓ ∉ R9 樣本）⟹ 剩餘判別力待補", flush=True)
    else:
        print(f"   ℓ = {l:2d}：區間寬 / null sd = {intervals[l]['width']/sd:.2f}", flush=True)

# 由 R9 之 null 尺度外推：sd ∝ 1/√(mode 數)
sd_anchor = J13["null"]["26"]["sd"]          # 偶 ℓ（奇 ℓ 之 sd ≡ 0，不可作錨）
for l in intervals:
    sd_ext = sd_anchor * np.sqrt((2 * 26 + 1) / (2 * l + 1))
    intervals[l]["sd_extrap"] = float(sd_ext)
    intervals[l]["width_over_sd"] = float(intervals[l]["width"] / sd_ext)
    print(f"   ℓ = {l:2d}：外推 sd = {sd_ext:.5f}（∝1/√(2ℓ+1)）⟹ 寬/sd = "
          f"{intervals[l]['width_over_sd']:.2f} ⟹ 剩餘判別力 ≈ "
          f"{'有' if intervals[l]['width_over_sd'] < 3 else '弱'}", flush=True)

ck("R12-剩餘判別力", all(intervals[l]["width_over_sd"] < 3 for l in intervals),
   "m_ℓ ≥ 2 之 ℓ：區間寬 < 3× 外推 sd（錨 ℓ=26）⟹ κ₃ 仍保留判別力（唔係全失）")

# (4) m_ℓ 普查（ℓ ≤ 200）：修正「ℓ ≥ 30 ⟹ m ≥ 2」之誤述
allm = {l: m_formula(l) for l in range(2, 201)}
zero = [l for l, v in allm.items() if v == 0]
one = [l for l, v in allm.items() if v == 1]
multi = [l for l, v in allm.items() if v >= 2]
print(f"\n(4) m_ℓ 普查（ℓ ≤ 200）：m=0 → {len(zero)} 個｜m=1 → {len(one)} 個｜"
      f"m≥2 → {len(multi)} 個 {multi[:12]}", flush=True)
ck("R12-普查", len(one) >= 20 and max(one) > 29 and len(multi) < 199 - 14,
   f"**更正**：m_ℓ ≥ 2 只出現於 {len(multi)} 個 ℓ（{multi[:6]} …），非「ℓ ≥ 30 全部」；"
   f"m_ℓ = 1 有 {len(one)} 個 ℓ ⟹ 形狀定理（a_ℓ ∝ I_ℓ）之有效範圍**遠超 ℓ ≤ 29**")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/16：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"ranks": {str(k): v for k, v in ranks.items()},
       "m_formula": {str(l): m_formula(l) for l in range(30, 41)},
       "intervals": {str(k): v for k, v in intervals.items()},
       "m_census": {"zero": zero, "one": len(one), "multi": multi},
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "16_r12_regime_b.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
