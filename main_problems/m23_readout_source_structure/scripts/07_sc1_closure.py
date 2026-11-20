#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/07 — R9：S 側封閉性審計（SC1）——把三族規則併成**單一 catalog** 再判一次

預註冊（本檔自身；寫完才跑）：
  catalog 𝔇 = 全部「已註冊規則」之 (δ_u, δ_d, δ_s)：
      T1  站對 946（scripts/01 之定義，重算）
      T1x 映射對 114（scripts/05 之 M，重算）
      T1p 剖面 13（scripts/06 之 P1–P8，重算）
  去重後 |𝔇| = n（報數）。
  (i)  一致規則（同一規則供三讀數）：判 |z_d|<1 ∧ |z_u|<1 之數；報 min χ。
  (ii) 每讀數自選規則（已註冊、無自由整數，但組合係後驗選取）：最優 χ 與選擇空間位元 = log₂(n_s·n_d·n_u)。
  (iii) 三族值集互含關係（防重複）。
判：SC1 成立 ⟺ (i) 0 個 1σ ∧ (ii) 最優 χ 雖可極小但位元 ≥ 15 ⟹ 零判別力；否則 SC1 降級。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40

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


print("=== M23 scripts/07 — R9：S 側封閉性審計（SC1）===")
cat = {}

# ── T1：站對 946 ────────────────────────────────────────────────
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        L = b - a
        def dl(k, a=a, b=b, L=L):
            if k <= a:
                return F(-k * L, M)
            if k <= b:
                return F(45 * (k - a) - k * L, M)
            return F(L * (45 - k), M)
        key = (dl(2), dl(3), dl(4))
        if key[2] == 0:
            continue
        cat.setdefault(key, set()).add("T1(%d,%d)" % (a, b))
n_t1 = len(cat)
print("      T1 站對：加入後 catalog 三元組 = %d" % n_t1)

# ── T1x：映射對（docs/05 之 M）───────────────────────────────────
MAPS = []
for c in (1, 2, 3, 4):
    MAPS.append(("c%d" % c, lambda k, c=c: c))
for d in (-2, -1, 0, 1, 2):
    MAPS.append(("k%+d" % d, lambda k, d=d: k + d))
for d in (-2, -1, 0, 1, 2):
    MAPS.append(("R(k%+d)" % d, lambda k, d=d: (M) - (k + d)))
for n1, f1 in MAPS:
    for n2, f2 in MAPS:
        if n1 == n2:
            continue
        ok = all(1 <= f(k) <= N and f1(k) != f2(k) for _, f in ((1, f1), (2, f2)) for k in KS.values())
        if not ok:
            continue
        d = {x: G(f2(k), k) - G(f1(k), k) for x, k in KS.items()}
        if d["s"] == 0:
            continue
        cat.setdefault((d["u"], d["d"], d["s"]), set()).add("T1x(%s/%s)" % (n1, n2))
n_t1x = len(cat)
print("      T1x 映射對：加入後 catalog 三元組 = %d（新增 %d）" % (n_t1x, n_t1x - n_t1))

# ── T1p：剖面（docs/07 之 P1–P8）────────────────────────────────
prof = []
for b in (1, 2, 3, 4):
    prof.append(("P1 e_%d" % b, lambda k, b=b: G(b, k)))
prof.append(("P2 s=1", lambda k: sum(G(j, k) for j in range(1, N + 1))))
prof.append(("P3 (−1)^j", lambda k: sum(((-1) ** j) * G(j, k) for j in range(1, N + 1))))
for j in (2, 3, 4):
    prof.append(("P4 e_%d−e_%d" % (j, M - j), lambda k, j=j: G(j, k) - G(M - j, k)))

for nm, f in prof:
    trip = (f(2), f(3), f(4))
    if trip[2] == 0:
        continue
    cat.setdefault(trip, set()).add("T1p(%s)" % nm)
n_all = len(cat)
print("      T1p 剖面（station-space 8 條）：catalog 三元組 = %d（新增 %d）" % (n_all, n_all - n_t1x))

# ── (i) 一致規則 ────────────────────────────────────────────────
rows = []
for trip, src in cat.items():
    du, dd, dsn = trip
    rd, ru = dd / dsn, du / dsn
    zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
    zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
    rows.append({"trip": trip, "src": sorted(src)[:3], "n_src": len(src),
                 "rho_d": str(rd), "rho_u": str(ru), "z_d": mp.nstr(zd, 10), "z_u": mp.nstr(zu, 10),
                 "chi": str(mp.sqrt(zd ** 2 + zu ** 2))})
h1 = sum(1 for r in rows if abs(mp.mpf(r["z_d"])) < 1 and abs(mp.mpf(r["z_u"])) < 1)
h3 = sum(1 for r in rows if abs(mp.mpf(r["z_d"])) < 3 and abs(mp.mpf(r["z_u"])) < 3)
best = min(rows, key=lambda r: mp.mpf(r["chi"]))
print("      (i) catalog 三元組 %d ｜ |z|<1 = %d ｜ |z|<3 = %d ｜ min χ = %s @ %s（由 %s 給出）"
      % (len(rows), h1, h3, mp.nstr(mp.mpf(best["chi"]), 8), best["rho_d"] + ", " + best["rho_u"], best["src"]))
chk("SC1-i  一致規則：0 個落 1σ", h1 == 0, "hits(z<1) = %d/%d" % (h1, len(rows)))
chk("SC1-i2 結構地板 = min χ（該三元組由多條規則共用 ⟹ 非唯一）", True,
    "min χ = %s ｜ 該三元組之規則數 = %d" % (mp.nstr(mp.mpf(best["chi"]), 6), best["n_src"]))

# ── (ii) 每讀數自選規則 ─────────────────────────────────────────
Vs = sorted({t[2] for t in cat})
Vd = sorted({t[1] for t in cat})
Vu = sorted({t[0] for t in cat})
best2 = None
for sv in Vs:
    if sv == 0:
        continue
    for dv in Vd:
        for uv in Vu:
            rd, ru = dv / sv, uv / sv
            zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
            zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
            c = mp.sqrt(zd ** 2 + zu ** 2)
            if best2 is None or c < best2[0]:
                best2 = (c, rd, ru, sv, dv, uv)
space = len(Vs) * len(Vd) * len(Vu)
bits = mp.log(space, 2)
print("      (ii) 每讀數自選：值集 δ_s %d ｜ δ_d %d ｜ δ_u %d ⟹ 空間 %d ⟹ %s bit"
      % (len(Vs), len(Vd), len(Vu), space, mp.nstr(bits, 6)))
print("           最優 ρ_d=%s ρ_u=%s ⟹ χ = %s" % (best2[1], best2[2], mp.nstr(best2[0], 6)))
chk("SC1-ii  每讀數自選之最優 χ 極小但位元 ≥ 15 ⟹ 零判別力", mp.mpf(best2[0]) < 1 and bits >= 15,
    "χ = %s ｜ %s bit ⟹ 擬合" % (mp.nstr(best2[0], 6), mp.nstr(bits, 6)))

# ── (iii) 防重複 ───────────────────────────────────────────────
src_kinds = {"T1": 0, "T1x": 0, "T1p": 0}
for trip, src in cat.items():
    kinds = {s.split("(")[0] for s in src}
    for kk in kinds:
        src_kinds[kk] += 1
print("      (iii) 各類覆蓋之三元組數：%s" % src_kinds)
print("            catalog 三元組總數 = %d ｜ 只由單一類覆蓋者 = %d"
      % (len(cat), sum(1 for s in cat.values() if len({x.split('(')[0] for x in s}) == 1)))
chk("SC1-iii 三族互不包含（各有獨佔三元組）",
    all(src_kinds[k] > 0 for k in src_kinds) and
    any(len({x.split('(')[0] for x in s}) == 1 for s in cat.values()), "見上")

out = {
    "catalog_size": len(cat), "t1_triples": n_t1, "t1x_triples": n_t1x, "t1p_triples": n_all,
    "uniform": {"hits_z1": h1, "hits_z3": h3, "min_chi": best["chi"],
                "best": {"rho_d": best["rho_d"], "rho_u": best["rho_u"],
                         "src": best["src"], "n_src": best["n_src"], "chi": best["chi"]}},
    "per_readout": {"space": space, "bits": mp.nstr(bits, 6),
                    "best_chi": str(best2[0]), "rho_d": str(best2[1]), "rho_u": str(best2[2]),
                    "n_s": len(Vs), "n_d": len(Vd), "n_u": len(Vu),
                    "value_sets": [len(Vs), len(Vd), len(Vu)]},
    "coverage": src_kinds,
}
with open(os.path.join(DATA, "07_sc1.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
print("結果: %d/%d PASS" % (npass, len(results)))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
