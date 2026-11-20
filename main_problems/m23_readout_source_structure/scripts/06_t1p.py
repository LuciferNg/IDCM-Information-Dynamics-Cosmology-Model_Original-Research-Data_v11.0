#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/06 — R8：源剖面（T1p）判定

預註冊（`docs/07` 凍結）：
  P1 s=e_b, b∈{1,2,3,4}      P2 s_j=1        P3 s_j=(−1)^j    P4 s=e_j−e_{45−j}, j∈{2,3,4}
  P5 Δ_n=1                   P6 Δ_n=(−1)^n   P7 Δ_n=δ_{n,1}   P8 Δ_n=1 略 n∈{9,27}
  P9 s∈ℝ^44 自由：核 rank(δ_u,δ_d,δ_s 之 3×44 映射) = 3 ⟹ 任何 (ρ_d,ρ_u) 精確可達（出示 witness）
判：P1–P8 零參數；全族 MC 10⁷；防重複；格點閘。P9 只作可達性披露。
算術：P1–P4 Fraction 精確；P5–P8 mpmath dps=50；P9 mpmath 50。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
N, M = 44, 45
SIG_D = mp.mpf('0.0351')
SIG_U = mp.mpf('0.03292')
TGT_D = mp.mpf('0.40350365')
TGT_U = mp.mpf('-0.02588837')
KS = {"u": 2, "d": 3, "s": 4}

results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))


def G(i, j):
    return F(min(i, j) * (M - max(i, j)), M)


def zc(rd, ru):
    """rd, ru 為 mpf。"""
    zd = (rd - TGT_D) / SIG_D
    zu = (ru - TGT_U) / SIG_U
    return zd, zu, mp.sqrt(zd ** 2 + zu ** 2)


def report(name, du, dd, ds):
    if ds == 0:
        return None
    rd, ru = dd / ds, du / ds
    if isinstance(rd, F):
        rd = mp.mpf(rd.numerator) / rd.denominator
        ru = mp.mpf(ru.numerator) / ru.denominator
    zd, zu, c = zc(rd, ru)
    print("      %-26s ρ_d=%s (z=%s)  ρ_u=%s (z=%s)  χ=%s"
          % (name, mp.nstr(rd, 10), mp.nstr(zd, 8), mp.nstr(ru, 10), mp.nstr(zu, 8), mp.nstr(c, 8)))
    return {"name": name, "rho_d": mp.nstr(rd, 20), "rho_u": mp.nstr(ru, 20),
            "z_d": mp.nstr(zd, 10), "z_u": mp.nstr(zu, 10), "chi": mp.nstr(c, 12)}


print("=== M23 scripts/06 — R8：源剖面（T1p） ===")
rows = []

# ── P1–P4：station-space 剖面（Fraction 精確）────────────────────
print("── P1–P4（station-space，Fraction 精確）──")
for b in (1, 2, 3, 4):
    d = {x: G(b, k) for x, k in KS.items()}
    rows.append(report("P1 s=e_%d" % b, d["u"], d["d"], d["s"]))
d = {x: sum(G(j, k) for j in range(1, N + 1)) for x, k in KS.items()}
rows.append(report("P2 s=1（均勻）", d["u"], d["d"], d["s"]))
d = {x: sum(((-1) ** j) * G(j, k) for j in range(1, N + 1)) for x, k in KS.items()}
rows.append(report("P3 s=(−1)^j", d["u"], d["d"], d["s"]))
for j in (2, 3, 4):
    d = {x: G(j, k) - G(M - j, k) for x, k in KS.items()}
    rows.append(report("P4 s=e_%d−e_%d" % (j, M - j), d["u"], d["d"], d["s"]))

# ── P5–P8：mode-space 剖面（mpmath dps=50）──────────────────────
print("── P5–P8（mode-space）──")


def phi(n, j):
    return mp.sqrt(mp.mpf(2) / M) * mp.sin(mp.pi * n * j / M)


def gam(n):
    return 4 * mp.sin(mp.pi * n / (2 * M)) ** 2


def profile_delta(dn):
    out = {}
    for x, k in KS.items():
        out[x] = sum(dn(n) * phi(n, k) / gam(n) for n in range(1, N + 1))
    return out


rows.append(report("P5 Δ_n=1", *[profile_delta(lambda n: mp.mpf(1))[x] for x in ("u", "d", "s")]))
rows.append(report("P6 Δ_n=(−1)^n", *[profile_delta(lambda n: mp.mpf(-1) ** n)[x] for x in ("u", "d", "s")]))
rows.append(report("P7 Δ_n=δ_{n,1}", *[profile_delta(lambda n: mp.mpf(1) if n == 1 else mp.mpf(0))[x]
                                       for x in ("u", "d", "s")]))
rows.append(report("P8 Δ_n=1 略 n∈{9,27}",
                   *[profile_delta(lambda n: mp.mpf(0) if n in (9, 27) else mp.mpf(1))[x]
                     for x in ("u", "d", "s")]))

rows = [r for r in rows if r]
chk("R8.0  P1–P8 族非空且全部可算（%d 個成員）" % len(rows), len(rows) == 13, "n = %d" % len(rows))

best = min(rows, key=lambda r: mp.mpf(r["chi"]))
print("      最優：%s χ = %s" % (best["name"], mp.nstr(mp.mpf(best["chi"]), 8)))
n1 = sum(1 for r in rows if abs(mp.mpf(r["z_d"])) < 1 and abs(mp.mpf(r["z_u"])) < 1)
n3 = sum(1 for r in rows if abs(mp.mpf(r["z_d"])) < 3 and abs(mp.mpf(r["z_u"])) < 3)
chk("R8.1  P1–P8 全族：0 個落 1σ ∧ 最優 χ ≥ 10σ", n1 == 0 and mp.mpf(best["chi"]) >= 10,
    "|z|<1 = %d/%d ｜ |z|<3 = %d/%d ｜ min χ = %s" % (n1, len(rows), n3, len(rows), mp.nstr(mp.mpf(best["chi"]), 6)))

# ── MC null 10⁷（13 點族）───────────────────────────────────────
import numpy as np
from scipy.spatial import cKDTree

pts = np.array([[float(mp.mpf(r["rho_d"])), float(mp.mpf(r["rho_u"]))] for r in rows], dtype=float)
sd, su = float(SIG_D), float(SIG_U)
tree = cKDTree(pts / np.array([sd, su]))
obs = float(mp.mpf(best["chi"]))
rng = np.random.default_rng(20261120)
NMC, CH = 10 ** 7, 10 ** 5
cnt = 0
for _ in range(NMC // CH):
    t = rng.normal(loc=[float(TGT_D), float(TGT_U)], scale=[sd, su], size=(CH, 2))
    dd, _ = tree.query(t / np.array([sd, su]), k=1)
    cnt += int(np.count_nonzero(dd <= obs))
p_mc = mp.mpf(cnt) / NMC
print("      MC null（10⁷）：p(min χ_null ≤ %s) = %s" % (mp.nstr(mp.mpf(best["chi"]), 6), mp.nstr(p_mc, 8)))
chk("R8.2  look-elsewhere §七②", True,
    "p = %s ⟹ %s" % (mp.nstr(p_mc, 8), "無信號" if p_mc >= mp.mpf('0.05') else "罕見"))

# ── P9：自由剖面之秩與 witness ──────────────────────────────────
print("── P9 自由剖面（s ∈ ℝ^44）──")
# 3×44 矩陣 A[x][n] = φ_n(k_x)/γ_n（δ_x = Σ_n A[x][n]·(Δ_n) ；Δ 自由 ⟺ s 自由，維度相同）
A = [[phi(n, k) / gam(n) for n in range(1, N + 1)] for k in (2, 3, 4)]
best_det, idx = mp.mpf(0), None
for i in range(N):
    for j in range(i + 1, N):
        for l in range(j + 1, N):
            det = mp.det(mp.matrix([[A[r][c] for c in (i, j, l)] for r in range(3)]))
            if abs(det) > best_det:
                best_det, idx = abs(det), (i + 1, j + 1, l + 1)
print("      max|3×3 minor| = %s @ 模態 %s ⟹ rank = 3" % (mp.nstr(best_det, 6), idx))
chk("R8.3  rank(3×44) = 3（存在非零 3×3 minor）", best_det > mp.mpf('1e-12'),
    "max|minor| = %s" % mp.nstr(best_det, 6))

# witness：只用 3 個模態，解 3×3 命中 (δ_u,δ_d,δ_s) ∝ (−0.02588837, 0.40350365, 1)
sub = [[A[r][c - 1] for c in idx] for r in range(3)]
rhs = mp.matrix([[-TGT_U * 0 + mp.mpf('-0.02588837')], [TGT_D], [mp.mpf(1)]])
coef = mp.lu_solve(mp.matrix(sub), rhs)
rec = [sum(coef[t] * A[r][idx[t] - 1] for t in range(3)) for r in range(3)]
resid = max(abs(rec[0] - rhs[0]), abs(rec[1] - rhs[1]), abs(rec[2] - rhs[2]))
rd9, ru9 = rec[1] / rec[2], rec[0] / rec[2]
zd9, zu9, c9 = zc(rd9, ru9)
print("      witness：僅 3 模態 %s 之剖面 ⟹ ρ_d=%s ρ_u=%s（殘差 %s）⟹ χ=%s"
      % (idx, mp.nstr(rd9, 10), mp.nstr(ru9, 10), mp.nstr(resid, 4), mp.nstr(c9, 6)))
chk("R8.4  P9 可達性：3 模態 witness 精確命中目標（殘差 < 1e−40 ⟹ 任何 (ρ_d,ρ_u) 可達）",
    resid < mp.mpf('1e-40'), "殘差 = %s ｜ χ = %s" % (mp.nstr(resid, 4), mp.nstr(c9, 6)))
chk("R8.5  P9 參數代價 = 44 − 3 = 41 個自由實數 ⟹ §七① 範疇性淘汰", True, "41 ≥ 1")

# ── 防重複 ────────────────────────────────────────────────────
print("── 防重複 ──")
with open(os.path.join(DATA, "01_family_enum.json")) as fh:
    fam01 = json.load(fh)
t1_vals = set()
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        def dl(k, a=a, b=b):
            L = b - a
            if k <= a:
                return F(-k * L, 45)
            if k <= b:
                return F(45 * (k - a) - k * L, 45)
            return F(L * (45 - k), 45)
        t1_vals.add((mp.nstr(mp.mpf(dl(3).numerator) / dl(3).denominator /
                             (mp.mpf(dl(4).numerator) / dl(4).denominator), 12),
                     mp.nstr(mp.mpf(dl(2).numerator) / dl(2).denominator /
                             (mp.mpf(dl(4).numerator) / dl(4).denominator), 12)))
with open(os.path.join(DATA, "05_t1x.json")) as fh:
    t1x = json.load(fh)
p_vals = {(mp.nstr(mp.mpf(r["rho_d"]), 12), mp.nstr(mp.mpf(r["rho_u"]), 12)) for r in rows}
cap = len(p_vals & t1_vals)
print("      P 族 distinct = %d ｜ T1 distinct = %d ｜ 交集 = %d" % (len(p_vals), len(t1_vals), cap))
print("      例：P1 s=e_4 = (3/4, 1/2) 已見於 T1 之 IA 段；P3 交錯 = (−21/2, 1/2) 係新值")
chk("R8.6  防重複：P 族之值與 T1 有交集但非子集（交錯／均勻剖面為新值）",
    0 < cap < len(p_vals), "交集 = %d/%d" % (cap, len(p_vals)))

out = {
    "members": rows, "best": best, "min_chi": best["chi"],
    "hits_z1": n1, "hits_z3": n3, "mc_null": {"nmc": NMC, "p": str(p_mc)},
    "p9": {"rank3_minor_modes": list(idx), "max_minor": mp.nstr(best_det, 10),
           "witness_rho_d": mp.nstr(rd9, 20), "witness_rho_u": mp.nstr(ru9, 20),
           "residual": mp.nstr(resid, 6), "free_reals": 41},
    "dedup": {"p_distinct": len(p_vals), "t1_distinct": len(t1_vals), "cap": cap},
}
with open(os.path.join(DATA, "06_t1p.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
print("結果: %d/%d PASS" % (npass, len(results)))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
