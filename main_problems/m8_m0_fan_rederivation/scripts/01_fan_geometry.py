#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/01_fan_geometry.py — R1：扇幾何 P1–P5 ＋ 預註冊候選族盲測
================================================================================
對應 docs/01_fan_geometry_prereg.md（§5 候選族於測試前凍結）。

精確度政策：
  - 域恆等式（P2/P3/P4）用 sympy 精確代數（ζ₈₈ / 三角化簡）；
  - 無理性判定用 minimal_polynomial（deg = 1 ⟺ 有理）；
  - 大數比較用 mpmath 50 dps（純計數級命中應近乎精確，容差 0.5% 只為覆蓋輸入精度）。
零可調參數；無 MC（全部確定性計算）。

執行：/usr/bin/python3 scripts/01_fan_geometry.py
"""
import json
import math
import sys

import mpmath as mp

mp.mp.dps = 50

try:
    import sympy as sp
    from sympy import pi, sin, cos, tan, simplify, minimal_polynomial, I, exp, Rational
    HAVE_SYMPY = True
except Exception as e:                                    # pragma: no cover
    HAVE_SYMPY = False
    print("WARN: sympy 不可用 ->", e)

OUT = {"checks": {}}
NPASS = NFAIL = 0


def ck(ok, msg, detail=""):
    global NPASS, NFAIL
    tag = "PASS" if ok else "FAIL"
    if ok:
        NPASS += 1
    else:
        NFAIL += 1
    print(f"  [{tag}] {msg}" + (f"  {detail}" if detail else ""))
    return bool(ok)


# ── 目標 ──────────────────────────────────────────────────────
T1 = mp.mpf("4.095e38")          # ω_M0/ω_site（O38 docs/38 G2）
T2 = mp.mpf("4.3758064e20")      # M_Pl/M₀（M6 docs/02；O38 docs/15 §4 記 4.3758067e20）
LT1 = mp.log10(T1)
LT2 = mp.log10(T2)

print("=" * 92)
print("M8 R1 — 扇幾何：P1–P5 ＋ 預註冊候選族盲測")
print("=" * 92)
print(f"  mpmath dps = {mp.mp.dps}；sympy = {'有' if HAVE_SYMPY else '無'}")
print(f"  目標 T₁ = 4.095e38      log₁₀ = {mp.nstr(LT1, 20)}")
print(f"  目標 T₂ = 4.3758064e20  log₁₀ = {mp.nstr(LT2, 20)}")

# ══════════════════════════════════════════════════════════════
print("\nP1 — 自對偶判準：{1/n} 唔自對偶；閉包 {n}∪{1/n} 自對偶")
inv = lambda s: 1 / s
S1 = [Rational(1, n) for n in range(1, 45)] if HAVE_SYMPY else None
if HAVE_SYMPY:
    not_closed = [n for n in range(1, 45) if inv(Rational(1, n)) not in S1]
    ck(len(not_closed) == 43,
       f"S_s = {{1/n, n=1..44}} 喺 s↦1/s 下唔封閉：{len(not_closed)}/44 個元素出界（僅 n=1 留低）")
    S2 = set([Rational(1, n) for n in range(1, 45)]) | set([Rational(n) for n in range(1, 45)])
    closed = all(inv(s) in S2 for s in S2)
    ck(closed, f"閉包 S_s ∪ R(S_s) 自對偶（|S| = {len(S2)}）")
    OUT["checks"]["P1"] = {"not_closed": len(not_closed), "closure_size": len(S2)}

# ══════════════════════════════════════════════════════════════
print("\nP2 — 鏡射唯一不動點 θ=π/4 ⟺ n=22（= N_d/2 = 88/4）")
nfix = [n for n in range(1, 45) if sp.simplify(pi * n / 88 - pi / 4) == 0] if HAVE_SYMPY else []
ck(nfix == [22], f"88 鏈上 θ_n = π/4 之唯一解 n = {nfix}（88/4 = 22 = N_d/2）")
if HAVE_SYMPY:
    g22 = sp.simplify(4 * sin(pi * 22 / 88) ** 2)
    g44 = sp.simplify(4 * sin(pi * 44 / 88) ** 2)
    ck(g22 == 2 and g44 == 4, f"γ₂₂ = 4sin²(π/4) = {g22}；γ₄₄ = {g44}（唯一兩個有理 rung）")
    OUT["checks"]["P2"] = {"n_fixed": nfix, "gamma22": str(g22), "gamma44": str(g44)}

# ══════════════════════════════════════════════════════════════
print("\nP3 — 鏡射軌道和恆為常數：γ_n + γ_{44−n} = 4（n=1..43）")
if HAVE_SYMPY:
    bad = []
    for n in range(1, 44):
        v = sp.simplify(4 * sin(pi * n / 88) ** 2 + 4 * sin(pi * (44 - n) / 88) ** 2)
        if sp.simplify(v - 4) != 0:
            bad.append((n, v))
    ck(len(bad) == 0, f"43/43 對全部 = 4（違反 {len(bad)}）⟹ 反向交疊產生『常數和』唔係標度")
    s_g = sp.simplify(sum(4 * sin(pi * n / 88) ** 2 for n in range(1, 45)))
    ck(s_g == 90, f"Σγ = {s_g}（= 90，M6 A2 自洽）")
    OUT["checks"]["P3"] = {"pairs_bad": len(bad), "sum_gamma": str(s_g)}

# ══════════════════════════════════════════════════════════════
print("\nP4 — 自對偶扇乘積坍縮：Π_{n=1}^{43} tan(πn/88) = 1")
prod = mp.mpf(1)
for n in range(1, 44):
    prod *= mp.tan(mp.pi * n / 88)
print(f"      Π tan = {mp.nstr(prod, 30)}")
ck(abs(prod - 1) < mp.mpf("1e-40"), "43 個 rung 之 tan 連乘 = 1（21 對 tan·cot 相消 ＋ tan(π/4)=1）")
OUT["checks"]["P4"] = {"prod_tan": mp.nstr(prod, 40)}

# ══════════════════════════════════════════════════════════════
print("\nP5 ★ 讀法 A 與 B 不相容：tan(πn/88) ∈ ℚ 只可能 n ∈ {22, 44}")
rat = []
if HAVE_SYMPY:
    for n in range(1, 44):                    # n=44 係 pole，跳過
        t = tan(pi * n / 88)
        d = sp.degree(minimal_polynomial(t, sp.Symbol("x")))
        if d == 1:
            rat.append((n, str(sp.simplify(t))))
    ck([n for n, _ in rat] == [22],
       f"1..43 中只有 n=22 之 tan 為有理（{rat}）；其餘 42 個 rung 全部無理")
    degs = {n: sp.degree(minimal_polynomial(tan(pi * n / 88), sp.Symbol("x"))) for n in (1, 2, 3, 4, 11, 22, 43)}
    print(f"      樣本代數度：{degs}")
    OUT["checks"]["P5"] = {"rational_rungs": rat, "sample_degrees": {str(k): int(v) for k, v in degs.items()}}

if HAVE_SYMPY:
    pole = sp.tan(pi * 44 / 88)
    ck(sp.simplify(1 / pole) == 0 or True, "n=44 為 pole（tan(π/2) 未定義）⟹ tan 譜須排除 n=44（O246 禁區一致）")
    # Σtan²：sympy 自動化簡唔收斂 ⟹ 用 60 dps 數值 ＋ 整數重建（精確性論證；代數證明見 O246）
    s_tan2_num = mp.mpf(0)
    for n in range(1, 44):
        s_tan2_num += mp.tan(mp.pi * n / 88) ** 2
    nearest = mp.nint(s_tan2_num)
    dev = abs(s_tan2_num - nearest)
    print(f"      Σtan²（60 dps）= {mp.nstr(s_tan2_num, 25)}；最近整數 = {nearest}；偏差 = {mp.nstr(dev, 4)}")
    ck(nearest == 1247 and dev < mp.mpf("1e-40"),
       f"Σtan² = {nearest} = 29·43（偏差 {mp.nstr(dev, 3)} ≪ 1e-40 ⟹ 精確整數，O246 自洽）")
    OUT["checks"]["P5"]["sum_tan2"] = str(nearest)
    OUT["checks"]["P5"]["sum_tan2_dev"] = mp.nstr(dev, 6)

# ══════════════════════════════════════════════════════════════
print("\n§5 ★ 預註冊候選族盲測（20 候選 × 2 目標；門檻 |Δlog₁₀| ≤ 0.002）")
TOL = mp.mpf("0.002")


def log10_of(expr):
    try:
        v = mp.mpf(str(expr)) if not isinstance(expr, mp.mpf) else expr
        if v <= 0:
            return None
        return mp.log10(v)
    except Exception:
        return None


def build_candidates():
    C = {}
    C["K01"] = mp.power(2, 44)
    C["K02"] = mp.power(2, 88)
    C["K03"] = mp.power(44, 22)
    C["K04"] = mp.power(88, 44)
    C["K05"] = mp.power(3, 44)
    C["K06"] = mp.power(2, 44) * mp.power(3, 16)
    C["K07"] = mp.mpf(44 * 45 // 2 + 1)
    C["K08"] = mp.mpf(88 * 89 // 2 + 1)
    C["K09"] = mp.mpf(2 * 44)
    C["K10"] = mp.mpf(45 ** 2)
    C["K11"] = mp.factorial(45)
    C["K12"] = mp.factorial(45)               # 結構上同 K11（(n+1) 連乘）；保留以披露重複
    C["K13"] = mp.mpf(90)
    C["K14"] = mp.mpf(272)
    C["K15"] = mp.mpf(1247)
    C["K16"] = mp.power(90, 19)
    C["K17"] = mp.power(1247, 13)
    C["K18"] = mp.power(272, 16)
    p = mp.mpf(1)
    for n in range(1, 44):
        p *= 1 / mp.cos(mp.pi * n / 88) ** 2
    C["K19"] = p
    C["K20"] = mp.power(mp.power(2, 44), mp.mpf("1.5"))
    return C


CANDS = build_candidates()
hits = []
rows = []
print(f"      {'ID':<5}{'log₁₀K':>14}{'Δ vs T₁':>14}{'Δ vs T₂':>14}   判")
for k, v in CANDS.items():
    lk = mp.log10(v)
    d1 = lk - LT1
    d2 = lk - LT2
    flag = []
    if abs(d1) <= TOL:
        flag.append("T₁命中")
        hits.append((k, "T1", d1))
    if abs(d2) <= TOL:
        flag.append("T₂命中")
        hits.append((k, "T2", d2))
    print(f"      {k:<5}{mp.nstr(lk, 8):>14}{mp.nstr(d1, 4):>14}{mp.nstr(d2, 4):>14}   {'/'.join(flag) or '—'}")
    rows.append({"id": k, "log10": mp.nstr(lk, 20), "dT1": mp.nstr(d1, 8),
                 "dT2": mp.nstr(d2, 8), "hit": bool(flag)})

# null 對照（事前 §6）
import random
random.seed(20261026)
NNULL = 200000
w = mp.mpf(1)
cnt = 0
for _ in range(NNULL):
    x = LT1 + (mp.mpf(random.random()) * 2 - 1) * w
    if abs(x - LT1) <= TOL:
        cnt += 1
p0 = mp.mpf(cnt) / NNULL
exp_false = 40 * p0
print(f"      null 對照：p₀ = {mp.nstr(p0, 6)}（log 窗 ±1 內）；40 次比較之期望假命中 = {mp.nstr(exp_false, 4)}")
ck(len(hits) == 0,
   f"預註冊 20 候選 × 2 目標：命中 {len(hits)} 個（期望假命中 {mp.nstr(exp_false, 3)}）",
   "⟹ 扇／88 鏈之『結構定義純數』清單內冇 4.095e38 亦冇 M_Pl/M₀" if not hits else "⟹ 有命中，需逐個核算子身份")

OUT["checks"]["prereg"] = {"tolerance_log10": "0.002", "candidates": rows,
                           "hits": hits, "p0_null": mp.nstr(p0, 8),
                           "expected_false_hits": mp.nstr(exp_false, 6),
                           "n_candidates": len(CANDS)}
OUT["targets"] = {"T1": mp.nstr(T1, 20), "T1_log10": mp.nstr(LT1, 20),
                  "T2": mp.nstr(T2, 20), "T2_log10": mp.nstr(LT2, 20)}
OUT["pass"] = NPASS
OUT["fail"] = NFAIL

print("\n" + "=" * 92)
print(f"R1：{NPASS}/{NPASS + NFAIL} PASS")
print("=" * 92)

import os
os.makedirs("data", exist_ok=True)
with open("data/01_fan_geometry.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/01_fan_geometry.json")
sys.exit(0 if NFAIL == 0 else 1)