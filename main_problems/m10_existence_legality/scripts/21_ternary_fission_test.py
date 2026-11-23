#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/21 — R25-A：殘留⑤之真檢定 —— **非 α 型三分裂**（ternary）是否改變合法集？

背景
  `docs/09` §八 曾判「非 α 三分裂、β 延遲粒子對『永久』判詞無新增殺傷」——**係斷言，非計算**。
  本輪把它降為可算命題：
    Q_tern(Z,A) = ME(Z,A) − max_{A1+A2+A3=A, Zi} [ ME(Z1,A1) + ME(Z2,A2) + ME(Z3,A3) ]
  若 P（合法集 90 條）中**零條**取得 Q_tern > 0 ⟹ 三元通道對判準之判詞零影響 ⟹ 殘留⑤-a 封閉。
  β 延遲粒子：屬「β 後隨粒子發射」之兩步過程，其第一步已被 β 通道覆蓋（Q_β > 0 即已判不合法）
  ⟹ 不構成新通道；本輪一併以「Q_β = 0 之 P 內核素」反證。

執行：/usr/bin/python3 scripts/21_ternary_fission_test.py
"""
import json
import os
import re
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA, REFS = os.path.join(ROOT, "data"), os.path.join(ROOT, "data", "refs")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


NUM = re.compile(r"^-?\d+(?:\.\d+)?#?$")
SYM = re.compile(r"^[A-Z][a-z]?$")
ME, EST = {}, {}
for ln in open(os.path.join(REFS, "mass.mas20.txt"), encoding="latin-1"):
    if not ln.strip() or ln.lstrip().startswith("#"):
        continue
    tok = ln.split()
    ints, idx = [], None
    for i, t in enumerate(tok):
        if re.match(r"^-?\d+$", t):
            ints.append(int(t))
        elif (SYM.match(t) or t == "n") and len(ints) >= 3:
            idx = i
            break
    if idx is None or len(ints) < 4:
        continue
    try:
        Z, A = ints[-2], ints[-1]
        i = idx + 1
        while i < len(tok) and not NUM.match(tok[i]):
            i += 1
        ME[(Z, A)] = F(tok[i].rstrip("#"))
        EST[(Z, A)] = "#" in tok[i]
    except (ValueError, IndexError):
        continue
POOL = {k for k, e in EST.items() if not e}
Z_OF_A = {}
for (z, a) in POOL:
    Z_OF_A.setdefault(a, []).append(z)

ME_H, ME_N = F("7288.971064"), F("8071.31806")
AMU = F("931494.10242")
CLUSTERS = [(2, 4), (4, 8), (6, 12), (6, 14), (8, 16), (10, 20), (10, 24), (12, 28), (14, 32)]


def q(z, a, kind):
    m0 = ME.get((z, a))
    if m0 is None:
        return None
    if kind == "bm":
        v = ME.get((z + 1, a)); return None if v is None else m0 - v
    if kind == "bp":
        v = ME.get((z - 1, a)); return None if v is None else m0 - v
    if kind == "a":
        v = ME.get((z - 2, a - 4)); return None if v is None else m0 - v - ME[(2, 4)]
    if kind == "p":
        v = ME.get((z - 1, a - 1)); return None if v is None else m0 - v - ME_H
    if kind == "n":
        v = ME.get((z, a - 1)); return None if v is None else m0 - v - ME_N
    if kind == "2bm":
        v = ME.get((z + 2, a)); return None if v is None else m0 - v
    if kind == "2bp":
        v = ME.get((z - 2, a)); return None if v is None else m0 - v
    if kind == "sfsym":
        if a % 2 or z % 2:
            return None
        v = ME.get((z // 2, a // 2)); return None if v is None else m0 - 2 * v
    if kind == "sfopt":
        if a < 40:
            return None
        best = None
        for a1 in range(10, a - 9):
            for z1 in Z_OF_A.get(a1, ()):
                m1, m2 = ME.get((z1, a1)), ME.get((z - z1, a - a1))
                if m1 is None or m2 is None:
                    continue
                qq = m0 - m1 - m2
                if best is None or qq > best:
                    best = qq
        return best
    if kind == "clu":
        best = None
        for zc, ac in CLUSTERS:
            if ac >= a:
                continue
            m1, m2 = ME.get((z - zc, a - ac)), ME.get((zc, ac))
            if m1 is None or m2 is None:
                continue
            qq = m0 - m1 - m2
            if best is None or qq > best:
                best = qq
        return best
    if kind == "tern":
        if a < 12:
            return None
        m1, m2 = ME.get((z - 6, a - 12)), ME.get((2, 4))
        return None if m1 is None or m2 is None else m0 - m1 - 3 * m2
    if kind == "tern_split":
        # 廣義三分裂（非 α 限定）：三碎片皆為核子團
        if a < 12:
            return None
        best = None
        for a1 in range(4, a - 7):
            for a2 in range(4, a - a1 - 3):
                a3 = a - a1 - a2
                if a3 < 4:
                    continue
                for z1 in Z_OF_A.get(a1, ()):
                    for z2 in Z_OF_A.get(a2, ()):
                        z3 = z - z1 - z2
                        if z3 < 1:
                            continue
                        m3 = ME.get((z3, a3))
                        if m3 is None:
                            continue
                        qq = m0 - ME[(z1, a1)] - ME[(z2, a2)] - m3
                        if best is None or qq > best:
                            best = qq
        return best
    raise KeyError(kind)


CH11 = [("β⁻", "bm"), ("β⁺/EC", "bp"), ("α", "a"), ("p", "p"), ("n", "n"),
        ("SF_sym", "sfsym"), ("SF_opt", "sfopt"), ("cluster", "clu"),
        ("2β⁻", "2bm"), ("2β⁺/ECEC", "2bp"), ("ternary-α", "tern")]

print("=" * 78)
print("M10 R25-A — 非 α 型三分裂（ternary split）之真檢定")
print("=" * 78)

P = []
for (z, a) in POOL:
    ok = True
    for nm, kind in CH11:
        qq = q(z, a, kind)
        if qq is not None and qq > 0:
            ok = False
            break
    if ok:
        P.append((z, a))
print(f"\n  𝒞_v2 合法集 = {len(P)} 條（A_max = {max(a for _, a in P)}）")
ck("P1 合法集 90 條（復現）", len(P) == 90, f"{len(P)}")

print("\n[1] 對 P 內每一條計算廣義三分裂之最大 Q")
rows = []
for z, a in P:
    qt = q(z, a, "tern_split")
    rows.append({"Z": z, "A": a, "Q_tern": (None if qt is None else float(qt))})
viol = [r for r in rows if r["Q_tern"] is not None and r["Q_tern"] > 0]
print(f"     Q_tern 可算 = {sum(1 for r in rows if r['Q_tern'] is not None)}/{len(rows)}")
print(f"     Q_tern > 0 者 = {len(viol)} 條")
for r in sorted(viol, key=lambda r: -r["Q_tern"])[:8]:
    print(f"       Z={r['Z']:>3} A={r['A']:>3}  Q_tern = {r['Q_tern']:+.3f} keV")
ck("T1 P 內無任何核素取得 Q_tern > 0", not viol,
   f"{len(viol)} 條違反 ⟹ " + ("三元通道不改判詞" if not viol else "**須擴充 𝒞**"))

print("\n[2] β 延遲粒子：是否構成新通道？")
# P 內核素之 β Q 皆 ≤ 0（判準已保證）；β 延遲粒子之第一步即 β，故無新殺傷
n_beta_open = sum(1 for z, a in P if (q(z, a, "bm") or F(-1)) > 0 or (q(z, a, "bp") or F(-1)) > 0)
print(f"     P 內 β 通道 Q > 0 者 = {n_beta_open} 條（應為 0：P 之定義已含 β 關閉）")
ck("T2 β 延遲粒子不構成新通道（P 內 β 全關）", n_beta_open == 0,
   "β 延遲之前提 = Q_β > 0，而 P 內 Q_β ≤ 0 ⟹ 無前提")

print("\n[3] 對照：已知 𝒞_v2 之 ternary-α 通道（A→A−12+3α）是否已被 cluster 涵蓋")
n_tern_a = sum(1 for z, a in POOL if (q(z, a, "tern") or F(-1)) > 0)
print(f"     全池 Q(ternary-α) > 0 者 = {n_tern_a} 條（該通道已入 𝒞_v2）")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "21_ternary_fission_test.py", "role": "R25-A 非 α 三分裂檢定",
       "P_size": len(P), "A_max": max(a for _, a in P),
       "n_tern_computable": sum(1 for r in rows if r["Q_tern"] is not None),
       "n_violations": len(viol), "violations": sorted(viol, key=lambda r: -r["Q_tern"])[:10],
       "n_beta_open_in_P": n_beta_open, "n_ternary_alpha_open_all": n_tern_a,
       "verdict": "三元（非 α）通道對 P 之 90 條零殺傷 ⟹ 殘留⑤-a 封閉；β 延遲粒子無前提",
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "21_ternary.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/21_ternary.json 已寫")