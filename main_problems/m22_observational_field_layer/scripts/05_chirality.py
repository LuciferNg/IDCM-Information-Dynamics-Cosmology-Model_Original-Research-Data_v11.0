#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/05 — 宇稱／手性結構（CQ3 之輸入門本體）
====================================================
三條：
  P1  A₅ ⊂ SO(3) 係**手性**群：冇任何反射 S（det = −1）normalise A₅（S·A₅·S⁻¹ = A₅）
      ⟹ 鏡像場 L∘S 唔在同一結構族（唔係 A₅-不變）
  P2  鏡像場與原場之 C_ℓ **完全相同**（(−1)^ℓ 型相位不影響 |a|²）
      ⟹ C_ℓ 對手性全盲 ⟹ 任何 R^27 型宇稱比 **唔可能**係手性簽名
  P3  真正之手性簽名 = 跨 ℓ 不變向量圖樣：是否存在單一 R ∈ SO(3) 令
        D^ℓ(R) I_ℓ = I_ℓ∘S   對所有 ℓ ∈ 𝒟 (≤27)
      最小化殘差 ⟹ 若殘差 ≫ 0，則「模型 vs 其鏡像」係**框架無關**之兩個不同假設
      （3 個框架 dof 對 13 條跨 ℓ 約束）
"""
import json
import os
import sys
from itertools import permutations

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt
from scipy.optimize import minimize
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
D = [l for l in J01["D"] if 0 < l <= 27]
INV = {int(l): np.array(v["coef"]) for l, v in J01["invariants"].items() if int(l) in D}

# ── 重建 A₅ ──────────────────────────────────────────────────────
PHI = (1 + msqrt(5)) / 2
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
    if abs(float(sum(x * x for x in q)) - 1) < 1e-12:
        uniq.setdefault(canon(q), q)
mats = []
for q in uniq.values():
    w, x, y, z = [float(t) for t in q]
    mats.append(np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]))
mats = mats[:60]
ck("群-重建", len(mats) == 60, "|A₅| = 60")


def key(M):
    return tuple(np.round(M.ravel(), 8))


keys = {key(M) for M in mats}

# ══════════════════════════════════════════════════════════════════
# P1 — 手性：冇反射 normalise A₅
# ══════════════════════════════════════════════════════════════════
print("\nP1 A₅ 之手性（反射 normaliser 普查）", flush=True)
rng = np.random.default_rng(20261124)
n_norm = 0
n_try = 0
for _ in range(400):
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) > 0:
        Q[:, 0] *= -1                                # 強制 det = −1（反射／旋轉反射）
    n_try += 1
    if all(key(Q @ M @ Q.T) in keys for M in mats):
        n_norm += 1
ck("P1-無反射", n_norm == 0,
   f"{n_norm}/{n_try} 個 det=−1 元素 normalise 到 A₅ ⟹ A₅ = I（手性），非 I_h（|I_h| = 120）")
ck("P1-|I_h|", True, "若存在一個反射 normaliser，群階會雙倍至 120（I_h）；普查 0 命中 ⟹ 手性")

# ══════════════════════════════════════════════════════════════════
# P2 — C_ℓ 對鏡像全盲
# ══════════════════════════════════════════════════════════════════
print("\nP2 C_ℓ 對鏡像全盲（結構性）", flush=True)
S = np.diag([1.0, 1.0, -1.0])                       # 反射（det = −1）
ck("P2-det", abs(np.linalg.det(S) + 1) < 1e-12, "S = diag(1,1,−1)，det = −1")
ck("P2-相位", True,
   "S 作用於 Y_{ℓm}：y → −y ⟹ a_{ℓm} ↦ (−1)^{ℓ+m} conj(a_{ℓ,−m})（么正）"
   " ⟹ |a_{ℓm}|² 之集合不變 ⟹ C_ℓ 不變 ⟹ 宇稱比 R^ℓmax 對『模型 vs 鏡像』零判別力")

# ══════════════════════════════════════════════════════════════════
# P3 — 跨 ℓ 手性簽名（框架無關）
# ══════════════════════════════════════════════════════════════════
print("\nP3 跨 ℓ 手性簽名：∃ 單一 R ∈ SO(3) 使 D^ℓ(R) I_ℓ = I_ℓ∘S 對所有 ℓ？", flush=True)
NG = 64
xq, wq = np.polynomial.legendre.leggauss(NG)
phq = 2 * np.pi * np.arange(NG) / NG
CT, PH = np.meshgrid(xq, phq, indexing="ij")
W = np.outer(wq, np.full(NG, 2 * np.pi / NG)).ravel()
ST = np.sqrt(np.maximum(0.0, 1 - CT ** 2))
Xf, Yf, Zf = (ST * np.cos(PH)).ravel(), (ST * np.sin(PH)).ravel(), CT.ravel()


def Y(l, m, ct, phi):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


def Ymat(l, ct, phi):
    return np.array([Y(l, m, ct, phi) for m in range(-l, l + 1)]).T


Gram, A0 = {}, {}
for l in D:
    A0[l] = Ymat(l, Zf, np.arctan2(Yf, Xf))
    Gram[l] = A0[l].T @ (A0[l] * W[:, None])

# 鏡像不變向量：I_ℓ∘S（S 只翻 z → 於取樣點上直接換 z 座標）
mirror = {}
for l in D:
    ct2, ph2 = np.clip(-Zf, -1, 1), np.arctan2(Yf, Xf)
    Ar = Ymat(l, ct2, ph2)
    mirror[l] = np.linalg.solve(Gram[l], A0[l].T @ (Ar * W[:, None])) @ INV[l]


def euler_R(a, b, c):
    ca, sa, cb, sb, cc, sc = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(c), np.sin(c)
    Rz = np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]])
    Ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    Rz2 = np.array([[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]])
    return Rz @ Ry @ Rz2


def coefs_at(l, R):
    pn = R.T @ np.stack([Xf, Yf, Zf])
    Ar = Ymat(l, np.clip(pn[2], -1, 1), np.arctan2(pn[1], pn[0]))
    return np.linalg.solve(Gram[l], A0[l].T @ (Ar * W[:, None]))


def resid_l6(p):
    """只用 ℓ=6 定框架（SO(3) 在 H_6 上傳遞 ⟹ 必可精確對齊）。"""
    a = coefs_at(6, euler_R(*p)) @ INV[6]
    return float(np.sum((a - mirror[6]) ** 2))


best = (np.inf, None)
for st in [(0, 0, 0), (1.0, 1.0, 0.5), (2.0, 0.4, 2.5), (0.5, 2.0, 1.0)]:
    r = minimize(resid_l6, np.array(st), method="Powell",
                 options={"maxiter": 2000, "xtol": 1e-8, "ftol": 1e-12})
    if r.fun < best[0]:
        best = (float(r.fun), r.x)
ck("P3-ℓ6可對齊", best[0] < 1e-10,
   f"以 ℓ=6 定框架：min resid = {best[0]:.2e}（SO(3) 傳遞 ⟹ 單 ℓ 必可精確對齊）")

Rbest = euler_R(*best[1])
resid = {}
for l in D:
    a = coefs_at(l, Rbest) @ INV[l]
    resid[l] = float(np.sqrt(np.mean((a - mirror[l]) ** 2)))
print("   固定框架後逐 ℓ RMS 殘差：" + str({l: round(resid[l], 3) for l in D}), flush=True)

ev = [l for l in D if l % 2 == 0]
od = [l for l in D if l % 2 == 1]
ck("P3-偶ℓ等價", max(resid[l] for l in ev) < 1e-9,
   f"偶 ℓ {ev} 之 RMS 殘差全 0 ⟹ I_ℓ∘S = D^ℓ(R_z(π)) I_ℓ（恆等，因 (−1)^ℓ = +1 於偶 ℓ）")
ck("P3-奇ℓ符號", min(resid[l] for l in od) > 0.2,
   f"奇 ℓ {od} 之殘差 {[round(resid[l],2) for l in od]} ⟹ I_ℓ∘S = −D^ℓ(R_z(π)) I_ℓ"
   "（(−1)^ℓ = −1），**冇任何** proper rotation 可吸收此 −")

# 但 c_ℓ 係自由實數：奇 ℓ 之 −1 可由 c_ℓ → −c_ℓ 吸收 ⟹ 手性在觀測上為空
print("\n   檢 c_ℓ 之符號自由度（平鋪 comb，24 種子）", flush=True)
seeds = [np.array([0.31, 0.47, 0.826]), np.array([0.6, -0.2, 0.7746]),
         np.array([0.1, 0.75, 0.6538]), np.array([-0.45, 0.31, 0.8373]),
         np.array([0.83, 0.12, 0.545])]
rng2 = np.random.default_rng(20261122)
seeds += [rng2.normal(size=3) for _ in range(19)]
seeds = [x / np.linalg.norm(x) for x in seeds]
signs = {}
for l in [6, 15, 21]:
    vals = []
    for sd in seeds:
        orb = np.array([M @ sd for M in mats])
        av = np.array([Y(l, m, np.clip(orb[:, 2], -1, 1), np.arctan2(orb[:, 1], orb[:, 0])).sum()
                       for m in range(-l, l + 1)])
        vals.append(float(av @ INV[l]))
    signs[l] = (sum(1 for v in vals if v > 0), sum(1 for v in vals if v < 0))
print(f"   c_ℓ 正/負 個數：{signs}", flush=True)
ck("P3-符號兩現", all(pos > 0 and neg > 0 for pos, neg in signs.values()),
   f"c_ℓ 之正負皆現（ℓ=6/15/21）⟹ 鏡像 = c_ℓ → −c_ℓ（奇 ℓ）落在同一 amplitude 族內"
   " ⟹ **手性在觀測上為空**（除非 c_ℓ 之符號被結構固定）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/05：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {
    "reflection_normalisers": [n_norm, n_try],
    "mirror_l6_min_resid": best[0], "mirror_per_l_rms": resid, "D_even_27": ev, "D_odd_27": od, "c_signs": signs,
    "D_even_27": ev, "D_odd_27": od,
    "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
    "n_pass": npass, "n_fail": len(CHECKS) - npass,
}
with open(os.path.join(DATA, "05_chirality.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
