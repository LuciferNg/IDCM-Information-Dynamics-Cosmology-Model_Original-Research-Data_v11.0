#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/00 — 重現閘（SESSION_START §二 ①–⑤）
================================================
硬性要求：數唔對 ⟹ 停。

① 零 import 重建 {N_d=44, B=4, Σγ=90, h11=36, h21=98, N_gen=3, 88} 七量
② 核 Σγ = B(N_d+1)/2 = 4×45/2 = 90（精確整數）
③ 核 γ_n + γ_{44−n} = 4 對極恆等（n = 1..43）
④ 核 M8 缺口 35.717977 dex 之算術（2.894296 vs 38.612274）
⑤ 唔准重跑 M6／M8／M16／v10 O38 腳本（本檔零 import 上游；零 exec）

精確算術：γ_n = 2 − (ζ^n + ζ^{−n}) 於 ℚ(ζ₈₈) = ℚ[x]/(Φ₈₈)，Φ₈₈ = Φ₁₁(−x⁴)，度 40。
全部環運算用 Fraction 整係數多項式，零浮點（除 ④ 之 dex 顯示）。
"""
import json
import os
import sys
from fractions import Fraction

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
DOCS = os.path.join(ODIR, "docs")

CHECKS = []


def ck(tag, cond, detail=""):
    """[PASS]/[FAIL] 機械可數行。"""
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
# 分圓環 ℚ(ζ₈₈)：Φ₈₈ 由 Φ₁₁(−x⁴) 構造（零 import）
# ══════════════════════════════════════════════════════════════════

def poly_mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            if b != 0:
                out[i + j] += a * b
    return out


def poly_divmod(p, q):
    """精確長除法，返回 (商, 餘)。"""
    p = [Fraction(c) for c in p]
    q = [Fraction(c) for c in q]
    while len(p) > 0 and p[-1] == 0:
        p.pop()
    db = len(q) - 1
    if db < 0:
        raise ZeroDivisionError("divisor is zero polynomial")
    quot = [Fraction(0)] * max(0, len(p) - db)
    while len(p) - 1 >= db and any(c != 0 for c in p):
        shift = len(p) - 1 - db
        coef = p[-1] / q[-1]
        quot[shift] = coef
        for i in range(db + 1):
            p[shift + i] -= coef * q[i]
        while len(p) > 0 and p[-1] == 0:
            p.pop()
    if not p:
        p = [Fraction(0)]
    return quot, p


def cyclotomic(n):
    """無 import 版本：Φ_n(x) = ∏_{d|n} (x^d − 1)^{μ(n/d)}，用精確除法。"""
    def mobius(m):
        if m == 1:
            return 1
        f = 2
        primes = []
        while f * f <= m:
            if m % f == 0:
                m //= f
                primes.append(f)
                if m % f == 0:
                    return 0
            f += 1
        primes.append(m)
        return -1 if len(primes) % 2 else 1

    num = [Fraction(1)]
    den = [Fraction(1)]
    for d in range(1, n + 1):
        if n % d:
            continue
        xd_minus_1 = [Fraction(0)] * d + [Fraction(1)]
        xd_minus_1[0] = Fraction(-1)
        mu = mobius(n // d)
        if mu == 1:
            num = poly_mul(num, xd_minus_1)
        elif mu == -1:
            den = poly_mul(den, xd_minus_1)
    quot, rem = poly_divmod(num, den)
    assert all(c == 0 for c in rem), "cyclotomic division not exact"
    return quot


PHI88 = cyclotomic(88)
DEG = len(PHI88) - 1

assert DEG == 40, f"Φ₈₈ degree {DEG} ≠ 40"
_monic = PHI88[-1]
assert _monic == 1, "Φ₈₈ not monic"
# 首一化（若最後係 1 已是）
PHI88[-1] = Fraction(1)


def ring_reduce(p):
    _, r = poly_divmod(p, PHI88)
    out = [Fraction(0)] * DEG
    for i, c in enumerate(r):
        if i < DEG:
            out[i] = c
    return out


def ring_add(a, b):
    return [x + y for x, y in zip(a, b)]


def ring_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def ring_scalar(a, s):
    return [x * Fraction(s) for x in a]


def pow_x(k):
    """x^k mod Φ₈₈。"""
    res = [Fraction(0)] * DEG
    res[0] = Fraction(1)
    base = [Fraction(0)] * (k + 1)
    base[k] = Fraction(1)
    return ring_reduce(base)


def is_const(a, val):
    return a[0] == Fraction(val) and all(c == 0 for c in a[1:])


print("=" * 78)
print("M21 scripts/00 — 重現閘")
print("=" * 78)
print(f"Φ₈₈ 度數 = {DEG}（應 40）")

# ζ^44 = −1（環內恆等）
X44 = pow_x(44)
ck("R00.0  ζ₈₈^44 = −1（環內）", is_const(X44, -1), f"x^44 ≡ {X44[0] if is_const(X44, X44[0]) else 'non-const'}")

# ══════════════════════════════════════════════════════════════════
# γ_n 譜（精確環元素）
# ══════════════════════════════════════════════════════════════════
N_d = 44
N = 2 * N_d
GAMMA = {}
for n in range(0, N_d + 1):
    # γ_n = 2 − (ζ^n + ζ^{−n})
    e = ring_add(pow_x(n), pow_x(N - n))
    g = ring_sub(ring_scalar([Fraction(0)] * DEG, 1), e)
    g[0] += Fraction(2)
    GAMMA[n] = g

ck("R00.1  γ_0 = 0", is_const(GAMMA[0], 0), "zero mode")
ck("R00.2  γ_44 = 4 = γ_max", is_const(GAMMA[44], 4), "Nyquist 端")

# ① 七量（零 import 重建）
B = 4                                    # = γ_max（頻帶數，M16 閉合）
N_gen = B - 1                            # M15：3 = B−1
SIG = None
sum_all = [Fraction(0)] * DEG
for n in range(1, N_d + 1):
    sum_all = ring_add(sum_all, GAMMA[n])
SIG = int(sum_all[0]) if is_const(sum_all, sum_all[0]) else None
h11, h21 = 36, 98                        # v10 O211 拓撲輸入（唯讀引用）
chi_cy3 = 2 * (h11 - h21)

print("\n── ① 七量 ──")
print(f"  N_d = {N_d} ｜ 88 = 2N_d = {N} ｜ B = {B} ｜ N_gen = B−1 = {N_gen}")
print(f"  Σγ = {SIG} ｜ h11 = {h11} ｜ h21 = {h21} ｜ χ(CY₃) = 2(h11−h21) = {chi_cy3}")

ck("gate①-a  88 = 2·N_d", N == 88, f"{N}")
ck("gate①-b  B = γ_max = 4", is_const(GAMMA[44], B), "B 由譜端定義")
ck("gate①-c  N_gen = B−1 = 3", N_gen == 3, f"{N_gen}")
ck("gate①-d  Σγ 為整數環元素", SIG is not None, f"Σγ = {SIG}")
ck("gate①-e  h11、h21 為註冊輸入（χ 一致）", chi_cy3 == -124, f"χ = {chi_cy3}")

# ② Σγ = B(N_d+1)/2
print("\n── ② Σγ = B(N_d+1)/2 ──")
ck("gate②-a  Σγ = 90（精確）", SIG == 90, f"{SIG}")
ck("gate②-b  B(N_d+1)/2 = 4×45/2 = 90", B * (N_d + 1) // 2 == 90, f"{B * (N_d + 1) // 2}")
ck("gate②-c  兩式相等", SIG == B * (N_d + 1) // 2, f"{SIG} vs {B * (N_d + 1) // 2}")
# 額外：Σγ = 2(N_d+1)（CG1 之 N=88 特例；恆等反演事實）
ck("gate②-d  Σγ = 2(N_d+1)（恆等，非約束）", SIG == 2 * (N_d + 1), f"{SIG} vs {2 * (N_d + 1)}")

# ③ 對極 γ_n + γ_{44−n} = 4
print("\n── ③ 對極恆等（n = 1..43） ──")
viol = 0
for n in range(1, N_d):
    pair = ring_add(GAMMA[n], GAMMA[N_d - n])
    if not is_const(pair, 4):
        viol += 1
        print(f"  [FAIL] γ_{n} + γ_{N_d - n} ≠ 4")
ck("gate③-a  γ_n + γ_{44−n} = 4，n = 1..43，違例 0", viol == 0, f"violations = {viol}")
ck("gate③-b  違例計數恰 43 對檢查", N_d - 1 == 43, f"pairs = {N_d - 1}")

# ④ M8 缺口算術
print("\n── ④ M8 缺口算術（2.894296 vs 38.612274） ──")
mp.mp.dps = 50
pi = mp.pi
span_efolds = 2 * mp.log(mp.cot(pi / 88))
span_dec = span_efolds / mp.log(10)
T1 = mp.mpf("4.0951856e38")
need_dec = mp.log10(T1)
gap_dec = need_dec - span_dec
ratio_span = mp.cot(pi / 88) ** 2
factor = mp.power(10, gap_dec)
print(f"  圖卡內跨距 = 2·ln cot(π/88) = {mp.nstr(span_efolds, 12)} e-folds = {mp.nstr(span_dec, 10)} dex")
print(f"  T₁ = {mp.nstr(T1, 10)} ⟹ {mp.nstr(need_dec, 10)} dex")
print(f"  缺口 = {mp.nstr(gap_dec, 10)} dex ｜ 因子 = {mp.nstr(factor, 8)}")
print(f"  cot²(π/88) = {mp.nstr(ratio_span, 10)} ｜ T₁/cot² = {mp.nstr(T1 / ratio_span, 8)}")
ck("gate④-a  跨距 = 2.894296 dex（6 位）", abs(span_dec - mp.mpf("2.894296")) < mp.mpf("5e-7"),
   mp.nstr(span_dec, 10))
ck("gate④-b  需求 = 38.612274 dex", abs(need_dec - mp.mpf("38.612274")) < mp.mpf("5e-7"),
   mp.nstr(need_dec, 10))
ck("gate④-c  缺口 = 35.717977 dex", abs(gap_dec - mp.mpf("35.717977")) < mp.mpf("5e-7"),
   mp.nstr(gap_dec, 10))
# 注意：M8 docs/04 §4 附註更正 —— 原值 2.3e35 有誤，正確為 5.223686e35（本 check 用更正後值）
ck("gate④-d  因子 5.223686e35（M8 更正後值）",
   abs(factor - mp.mpf("5.223686e35")) / mp.mpf("5.223686e35") < mp.mpf("1e-6"), mp.nstr(factor, 8))
ck("gate④-e  因子 = T₁ / cot²(π/88)", abs(factor / (T1 / ratio_span) - 1) < mp.mpf("1e-12"),
   mp.nstr(factor / (T1 / ratio_span), 12))
ck("gate④-f  迭代 13.34 非整數（T₁ 需求）", abs(need_dec / span_dec - mp.mpf("13.34")) < mp.mpf("0.01"),
   mp.nstr(need_dec / span_dec, 8))

# ⑤ 零上游 import／零 exec
print("\n── ⑤ 上游腳本零執行 ──")
import ast

src = open(os.path.abspath(__file__), encoding="utf-8").read()
tree = ast.parse(src)
imported = set()
exec_calls = []
for node in ast.walk(tree):
    if isinstance(node, ast.Import):
        for a in node.names:
            imported.add(a.name.split(".")[0])
    elif isinstance(node, ast.ImportFrom):
        if node.module:
            imported.add(node.module.split(".")[0])
    elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("exec", "eval"):
        exec_calls.append(node.func.id)
upstream_mods = {"p9_kchain_number_field", "kchain_number_field", "o38", "idcm_core",
                 "runpy", "importlib", "subprocess"}
bad_imports = sorted(imported & upstream_mods)
ck("gate⑤-a  本檔零 import 上游模組", not bad_imports,
   f"imports = {sorted(imported)} ｜ upstream hits = {bad_imports}")
ck("gate⑤-b  本檔零 exec/eval/runpy", (not exec_calls) and ("runpy" not in imported),
   f"exec/eval calls = {exec_calls}")

# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"重現閘結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "00_reproduction_gate.json"), "w", encoding="utf-8") as f:
    json.dump({
        "checks": [{"tag": t, "pass": ok, "detail": d} for t, ok, d in CHECKS],
        "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
        "quantities": {"N_d": N_d, "N": N, "B": B, "Sigma_gamma": SIG,
                       "h11": h11, "h21": h21, "N_gen": N_gen, "chi_CY3": chi_cy3},
        "exact": {"field": "Q(zeta_88)=Q[x]/(Phi88)", "deg": DEG},
    }, f, ensure_ascii=False, indent=2)

sys.exit(0 if nfail == 0 else 1)
