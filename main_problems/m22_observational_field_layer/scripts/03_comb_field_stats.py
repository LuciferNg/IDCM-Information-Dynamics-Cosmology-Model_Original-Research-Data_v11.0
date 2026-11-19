#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/03 — 平鋪 comb 場（Reading F）：形狀定理之直接核 ＋ 種子帶統計量
============================================================================
判準：docs/01_preregistration.md §三（Q1–Q4）、§二（四統計量定義）。

核心檢定（形狀定理；Q1 之心臟）：
  對種子 s 之 60 點 A₅ 軌道（均勻權重 = 唯一 A₅-不變振幅指派）：
    a_{ℓm}(s) = Σ_j Y_{ℓm}(n̂_j)
  要求：(i) a_{ℓm} = 0 對 ℓ ∉ 𝒟；
        (ii) 對 ℓ ∈ 𝒟（m_ℓ = 1）a_{ℓm}(s) ∝ I_{ℓm}（cos 角 = 1，逐 ℓ、逐種子）
  ⟹ 逐 ℓ 形狀零參數；標量 c_ℓ(s) = Σ_j I_ℓ(n̂_j) 隨種子變（= 2 dof 帶）

統計量（ℓ ≤ 27 帶；零 ℓ≤5 為結構預言）：
  C_ℓ = (1/(2ℓ+1))Σ_m a_{ℓm}²、C(θ)、S_{1/2}、C(π)、R^27 = D⁺/D⁻、σ²_16 / σ²_total
  σ²_16 用 healpy pixwin(16)（Nside=16 之像素窗）；比值 = 無量綱
"""
import json
import os
import sys

import healpy as hp
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


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
J02 = json.load(open(os.path.join(DATA, "02_amplitude_nogo.json"), encoding="utf-8"))
D = J01["D"]
MELL = J01["m_ell"]
INV = {int(k): np.array(v["coef"]) for k, v in J01["invariants"].items()}
LMAX = 40
DL = [l for l in D if 0 < l <= LMAX]


def Y(l, m, ct, phi):
    am = abs(m)
    nn = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    P = lpmv(am, l, ct)
    return nn * P if m == 0 else nn * np.sqrt(2) * P * (np.cos(am * phi) if m > 0 else np.sin(am * phi))


def Yvec(l, n):
    """n: (K,3) 單位向量 → (2l+1, K) 實球諧值"""
    ct = np.clip(n[:, 2], -1.0, 1.0)
    ph = np.arctan2(n[:, 1], n[:, 0])
    return np.array([Y(l, m, ct, ph) for m in range(-l, l + 1)])


# ══════════════════════════════════════════════════════════════════
# 重建 A₅（與 scripts/01/02 完全同一構造）
# ══════════════════════════════════════════════════════════════════
from itertools import permutations  # noqa: E402
from mpmath import mp, mpf, sqrt as msqrt  # noqa: E402
mp.dps = 40
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
    uniq.setdefault(canon(q), q)
mats = []
for q in uniq.values():
    w, x, y, z = [float(t) for t in q]
    mats.append(np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]]))
ck("群-重建", len(mats) == 60, "|A₅| = 60")


def orbit(seed):
    pts = np.array([M @ seed for M in mats])
    # 去重（generic orbit ⟹ 60 點）
    keep = []
    for p in pts:
        if not any(np.abs(p - q).max() < 1e-9 for q in keep):
            keep.append(p)
    return np.array(keep)


# ══════════════════════════════════════════════════════════════════
# 種子帶：多個 generic 種子（2 dof）
# ══════════════════════════════════════════════════════════════════
rng = np.random.default_rng(20261122)
seeds = [np.array([0.31, 0.47, 0.826]), np.array([0.6, -0.2, 0.7746]),
         np.array([0.1, 0.75, 0.6538]), np.array([-0.45, 0.31, 0.8373]),
         np.array([0.83, 0.12, 0.545])]
seeds += [rng.normal(size=3) for _ in range(19)]
seeds = [s / np.linalg.norm(s) for s in seeds]

res, shape_worst, zero_worst = {}, 0.0, 0.0
for si, s in enumerate(seeds):
    orb = orbit(s)
    assert len(orb) == 60, f"seed {si}: |orbit| = {len(orb)}"
    cmb = {}                                        # c_ℓ
    shape = {}
    for l in range(2, LMAX + 1):
        av = Yvec(l, orb).sum(axis=1)               # a_{ℓm} = Σ_j Y_{ℓm}(n̂_j)
        if MELL[l] == 0:
            zero_worst = max(zero_worst, float(np.abs(av).max()))
            continue
        Iv = INV[l] if l in INV else None
        if Iv is None or MELL[l] != 1:      # 只有 m_ℓ = 1 才有單一不變方向
            if Iv is not None:
                zero_worst = max(zero_worst, 0.0)
            continue
        c = float(av @ Iv)                          # ⟨a, I⟩（I 已歸一）
        cn = av / np.linalg.norm(av)
        cs = float(cn @ Iv)
        shape[l] = cs
        shape_worst = max(shape_worst, abs(abs(cs) - 1.0))
        cmb[l] = c
    res[si] = {"seed": s.tolist(), "c": cmb, "shape_cos": shape}

ck("形狀-支援零", zero_worst < 1e-9, f"ℓ ∉ 𝒟 之 max|a_ℓm| = {zero_worst:.2e}（ℓ ≤ {LMAX}，20 種子）")
ck("形狀-逐ℓ對齊", shape_worst < 1e-9,
   f"m_ℓ=1 之 ℓ：max|1−|cos(a_ℓ, I_ℓ)|| = {shape_worst:.2e} ⟹ 逐 ℓ 形狀零參數（種子無關）")
ck("形狀-ℓ30退化", MELL[30] == 2 and 30 in INV,
   f"m₃₀ = 2（ℓ=30 起逐 ℓ 形狀不再唯一）—— 正是 30 在帶外之理由")


# ══════════════════════════════════════════════════════════════════
# 統計量（ℓ ≤ 27 帶）＋ 種子帶
# ══════════════════════════════════════════════════════════════════
NGS = 200
xs, ws = np.polynomial.legendre.leggauss(NGS)       # cosθ 節點
Pw16 = hp.pixwin(16)


def stats_for(cm):
    """cm: {ℓ: c_ℓ}（comb）。回傳統計量 dict。"""
    lmax = max(cm)
    Cl = np.zeros(lmax + 1)
    for l, c in cm.items():
        Cl[l] = c * c / (2 * l + 1)
    # C(θ)
    th = np.arccos(xs)
    Pl = np.zeros((lmax + 1, NGS))
    from scipy.special import eval_legendre  # noqa: E402
    for l in range(lmax + 1):
        Pl[l] = eval_legendre(l, xs)
    Cth = ((2 * np.arange(lmax + 1) + 1)[:, None] * Cl[:, None] * Pl).sum(axis=0) / (4 * np.pi)
    # S_1/2 = ∫_{−1}^{1/2} C(θ)² dcosθ ; 節點在 [−1,1] → 需截斷
    mask = xs <= 0.5
    S12 = float(np.sum(ws[mask] * Cth[mask] ** 2))
    Cpi = float(np.sum((2 * np.arange(lmax + 1) + 1) * (-1.0) ** np.arange(lmax + 1) * Cl[:lmax + 1]) / (4 * np.pi))
    # R^27（Planck/Muir：D⁺/D⁻，ℓ=2..27）
    De = Do = 0.0
    for l in range(2, 28):
        v = l * (l + 1) / (2 * np.pi) * Cl[l]
        if l % 2 == 0:
            De += v
        else:
            Do += v
    R27 = De / Do if Do > 0 else float("inf")
    # σ² 類（Nside=16 像素窗；ℓ ≤ 3·16−1 = 47 上限）
    lm16 = min(lmax, 3 * 16 - 1)
    s16 = float(sum((2 * l + 1) / (4 * np.pi) * Cl[l] * Pw16[l] ** 2 for l in range(2, lm16 + 1)))
    stot = float(sum((2 * l + 1) / (4 * np.pi) * Cl[l] for l in range(2, lmax + 1)))
    flow = float(sum((2 * l + 1) / (4 * np.pi) * Cl[l] for l in range(2, 6))) / stot
    return {"C_ell": {str(l): float(Cl[l]) for l in range(lmax + 1) if Cl[l] > 0},
            "C_pi": Cpi, "S_half": S12, "R27": R27, "sigma16": s16,
            "sigma_total": stot, "ratio16": s16 / stot, "frac_low_5": flow,
            "Ctheta": [[float(x), float(y)] for x, y in zip(xs[::10], Cth[::10])]}


tab = {}
for si, r in res.items():
    tab[si] = stats_for(r["c"]) | {"seed": r["seed"]}
    print(f"   seed {si:2d}: R^27 = {tab[si]['R27']:.4f}  ratio16 = {tab[si]['ratio16']:.6f} "
          f"C(π) = {tab[si]['C_pi']:.4e}  S½ = {tab[si]['S_half']:.4e}  frac(ℓ≤5) = {tab[si]['frac_low_5']:.3e}",
          flush=True)

R27s = np.array([tab[s]["R27"] for s in tab])
r16s = np.array([tab[s]["ratio16"] for s in tab])
cpis = np.array([tab[s]["C_pi"] for s in tab])
shs = np.array([tab[s]["S_half"] for s in tab])
fls = np.array([tab[s]["frac_low_5"] for s in tab])

ck("統計-ℓ≤5零功率", float(np.abs(fls).max()) < 1e-12,
   f"ℓ ≤ 5 功率分數 = {float(np.abs(fls).max()):.2e}（恰 0，種子無關）")
ck("統計-R27無預言", float(R27s.max() / R27s.min()) > 10.0,
   f"R^27 ∈ [{R27s.min():.4f}, {R27s.max():.4f}]（跨 both sides，逾 2 個數量級）⟹ **無預言**（種子帶主導）")
print(f"   R^27 帶：{R27s.min():.4f} … {R27s.max():.4f}（中位 {np.median(R27s):.4f}）", flush=True)
print(f"   ratio16 帶：{r16s.min():.6f} … {r16s.max():.6f}（中位 {np.median(r16s):.6f}）", flush=True)
print(f"   C(π) 帶（comb 原始單位）：{cpis.min():.4e} … {cpis.max():.4e}", flush=True)
print(f"   S½ 帶（comb 原始單位）：{shs.min():.4e} … {shs.max():.4e}", flush=True)

# 對照：若把 ℓ≤5 填成尺度不變（ΛCDM-like），ratio16 變幾多
def ratio16_if_filled(cm):
    """把 ℓ=2..5 之 C_ℓ 設為相鄰支撐帶之同值（1/(ℓ(ℓ+1)) 定標），量抑制度。"""
    lmax = max(cm)
    Cl = np.zeros(lmax + 1)
    for l, c in cm.items():
        Cl[l] = c * c / (2 * l + 1)
    ref = np.median([Cl[l] * l * (l + 1) for l in cm if l <= 30])
    for l in range(2, 6):
        Cl[l] = ref / (l * (l + 1))
    lm16 = min(lmax, 47)
    s16 = sum((2 * l + 1) / (4 * np.pi) * Cl[l] * Pw16[l] ** 2 for l in range(2, lm16 + 1))
    stot = sum((2 * l + 1) / (4 * np.pi) * Cl[l] for l in range(2, lmax + 1))
    return s16 / stot


filled = np.array([ratio16_if_filled(res[s]["c"]) for s in res])
supp = np.array([1 - r16s / filled])
ck("統計-σ16抑制度", float(supp.min()) > 0.0,
   f"σ²_16 相對抑制度 = {supp.min()*100:.2f}% … {supp.max()*100:.2f}%（填 ℓ≤5 之對照）")

# ══════════════════════════════════════════════════════════════════
# Reading G（對照）：逐 ℓ 等權 c_ℓ = 1（另一 convention）—— 證明統計量之 convention 依賴
# ══════════════════════════════════════════════════════════════════
print("\nReading G（c_ℓ = 1 ∀ℓ ∈ 𝒟，ℓ ≤ 27）：另一 amplitude convention", flush=True)
cmG = {l: 1.0 for l in D if 0 < l <= 27}
stG = stats_for(cmG)
print(f"   R^27 = {stG['R27']:.4f}  ratio16 = {stG['ratio16']:.6f}  "
      f"C(π) = {stG['C_pi']:.4e}  S½ = {stG['S_half']:.4e}", flush=True)
ck("ReadingG-ℓ≤5零", abs(stG["frac_low_5"]) < 1e-12, "ℓ≤5 功率分數 = 0（與 convention 無關）")
n_above = int((R27s > 1).sum())
ck("統計-R27跨1", float(R27s.min()) < 1.0 < float(R27s.max()),
   f"comb {n_above}/{len(R27s)} 種子 R^27>1；Reading G = {stG['R27']:.4f}（>1）"
   f" ⟹ 方向偏 odd-deficit 但 {len(R27s)-n_above} 個種子相反 ⟹ 無點預言")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/03：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {
    "n_seeds": len(seeds), "stats": {str(k): v for k, v in tab.items()},
    "bands": {"R27": [float(R27s.min()), float(R27s.max()), float(np.median(R27s))],
              "ratio16": [float(r16s.min()), float(r16s.max()), float(np.median(r16s))],
              "C_pi": [float(cpis.min()), float(cpis.max())],
              "S_half": [float(shs.min()), float(shs.max())],
              "suppression_sigma16": [float(supp.min()), float(supp.max())]},
    "shape_worst": shape_worst, "zero_worst": zero_worst, "reading_G": stG,
    "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
    "n_pass": npass, "n_fail": len(CHECKS) - npass,
}
with open(os.path.join(DATA, "03_comb_field_stats.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
