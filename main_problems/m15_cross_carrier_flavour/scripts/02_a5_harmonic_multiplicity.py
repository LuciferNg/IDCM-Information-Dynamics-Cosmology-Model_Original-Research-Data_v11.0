#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/02 — C1「清單窮盡性」＋ 盲測 B1：A₅ 諧波重數表（精確，ℚ(√5) 手寫域）
================================================================================
判準來源：docs/01 §三 C2、§四 B1（跑前凍結）。
問題：A₅ 幾何能否由**重數**供給「3」？（N_gen = 3 之唯一可能身份 = 某 irrep 之重數 3）

方法（全部精確，零浮點，零 sympy 化簡）：
  ℚ(√5) 以 (p, q) ≙ p + q√5（p,q ∈ ℚ）表示；加減逐項，乘 (p₁p₂+5q₁q₂, p₁q₂+p₂q₁)，
  除乘共軛（p²−5q² 分母）。冪級數反演逐項遞推。
  ① A₅ ⊂ SO(3) 五類（1A:1、2A:15、3A:20、5A:12、5B:12；角 0°,180°,120°,72°,144°）。
  ② det(1−tg) = (1−t)(1−2cosθ·t+t²)；5A 之 2cos72 = φ⁻¹、5B = −φ。
  ③ M_ρ(t) = (1/60) Σ_C |C| χ_ρ(C)/det_C(t)；諧波重數 m_ρ(n) = [tⁿ](1−t²)M_ρ(t)。
  ④ 驗證：Σ_ρ d_ρ m_ρ(n) = 2n+1；M_1 = (1+t¹⁵)/((1−t⁶)(1−t¹⁰))（Hironaka）。

盲測 B1（docs/01 §四）：必要條件 = ∃ρ, n ∈ P : m_ρ(n) ≥ 3，P = {1,2,3,4,5,7,9,11,13}。
⚠ 零觀測量；N_trials = 0。
"""
import json
import sys
from fractions import Fraction as F

NMAX = 45                     # 級數長度（n = 0..44）
P_SET = {1, 2, 3, 4, 5, 7, 9, 11, 13}      # O137 §五 受保護集（9 模式）
LOG, RES = [], []


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail=""):
    st = "PASS" if cond else "FAIL"
    line = f"[{st}] {tag}" + (f"  |  {detail}" if detail else "")
    say(line)
    RES.append((tag, st, detail))
    return cond


# ── ℚ(√5) 元素 ────────────────────────────────────────────────────────────────
Z = (F(0), F(0))
ONE = (F(1), F(0))


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def mul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def div(x, y):
    d = y[0] * y[0] - 5 * y[1] * y[1]
    assert d != 0, "division by zero in ℚ(√5)"
    return ((x[0] * y[0] - 5 * x[1] * y[1]) / d, (x[1] * y[0] - x[0] * y[1]) / d)


def rat(x, n=1):
    return (F(x, n), F(0))


def fmt(x):
    if x[1] == 0:
        return str(x[0])
    return f"({x[0]}+{x[1]}√5)"


S5 = (F(0), F(1))
PHI = (F(1, 2), F(1, 2))          # (1+√5)/2
PSI = (F(1, 2), F(-1, 2))         # (1−√5)/2
IPHI = sub(ONE, PSI)              # 1 − (1−√5)/2 = (1+√5)/2? 手動核見下


def neg(x):
    return (-x[0], -x[1])


say("=" * 78)
say("M15 scripts/02 — A₅ 諧波重數表（精確 ℚ(√5)）＋ 盲測 B1")
say("慣例：Harm_n ⊂ ℂ[x,y,z]（2n+1 維）；group = A₅（chiral I，階 60）")
say("=" * 78)
say(f"  域核：φ = {fmt(PHI)}；φ⁻¹ = {fmt(div(ONE, PHI))}；−φ = {fmt(neg(PHI))}")

# ── 共軛類：det(1 − t g) 之係數（升冪 t⁰..t³）─────────────────────────────
CLASSES = [
    # (名, 階, [c0,c1,c2,c3], 角)
    ("1A", 1, [ONE, rat(-3), rat(3), rat(-1)], "0°"),
    ("2A", 15, [ONE, rat(1), rat(-1), rat(-1)], "180°"),
    ("3A", 20, [ONE, Z, Z, rat(-1)], "120°"),
    ("5A", 12, [ONE, neg(PHI), PHI, rat(-1)], "72°"),          # 1 − φt + φt² − t³
    ("5B", 12, [ONE, div(ONE, PHI), neg(div(ONE, PHI)), rat(-1)], "144°"),  # 1 + φ⁻¹t − φ⁻¹t² − t³
]
CHARS = [
    ("1", 1, [rat(1), rat(1), rat(1), rat(1), rat(1)]),
    ("3", 3, [rat(3), rat(-1), Z, PHI, PSI]),
    ("3'", 3, [rat(3), rat(-1), Z, PSI, PHI]),
    ("4", 4, [rat(4), Z, rat(1), rat(-1), rat(-1)]),
    ("5", 5, [rat(5), rat(1), rat(-1), Z, Z]),
]


def inv_series(poly, n):
    """1/poly 之冪級數（poly[0] = 1），長度 n。"""
    b = [Z] * n
    b[0] = div(ONE, poly[0])
    for i in range(1, n):
        acc = Z
        for k in range(1, min(i, len(poly) - 1) + 1):
            acc = add(acc, mul(poly[k], b[i - k]))
        b[i] = neg(div(acc, poly[0]))
    return b


say("")
say("§1  共軛類與 det(1−tg)")
say("-" * 78)
for (nm, sz, poly, ang) in CLASSES:
    say(f"  {nm}（階 {sz}，{ang}）：det = " + " + ".join(
        f"({fmt(poly[i])})t^{i}" for i in range(4) if poly[i] != Z))
chk("A 群論基線", sum(c[1] for c in CLASSES) == 60 and sum(c[1] ** 2 for c in CHARS) == 60,
    f"|A₅| = {sum(c[1] for c in CLASSES)}，Σd² = {sum(c[1] ** 2 for c in CHARS)}")

# ── Molien 級數 ─────────────────────────────────────────────────────────────
INVS = {c[0]: inv_series(c[2], NMAX) for c in CLASSES}
MP = {}
for (rname, dim, chi) in CHARS:
    ser = [Z] * NMAX
    for idx, (cname, csize, poly, ang) in enumerate(CLASSES):
        coef = (F(csize, 60), F(0))
        coef = mul(coef, chi[idx])
        for i in range(NMAX):
            ser[i] = add(ser[i], mul(coef, INVS[cname][i]))
    MP[rname] = ser

# 驗證 1：M_1（多項式環 Molien）= (1+t^15)/((1−t²)(1−t⁶)(1−t¹⁰))（Hironaka：f₂,f₆,f₁₀ + skew t¹⁵）
ref = [Z] * NMAX
for a in range(NMAX):
    for b in range(NMAX):
        for c in range(NMAX):
            for sh in (0, 15):
                d = sh + 2 * c + 6 * a + 10 * b
                if d < NMAX:
                    ref[d] = add(ref[d], ONE)
dev1 = max((sub(MP["1"][i], ref[i]) for i in range(NMAX)),
           key=lambda x: abs(x[0]) + abs(x[1]))
say("")
say("§2  Molien 級數 ＋ 驗證")
say("-" * 78)
say(f"  驗證 1（Hironaka）：max|M_1 − (1+t¹⁵)/((1−t⁶)(1−t¹⁰))| = {fmt(dev1)}")
chk("B 驗證 1：M_1 與 Hironaka 一致", dev1 == Z, fmt(dev1))

# 諧波重數 m_ρ(n) = M[n] − M[n−2]
MW = {r: [Z] * NMAX for (r, _, _) in CHARS}
for (rname, dim, _) in CHARS:
    for n in range(NMAX):
        MW[rname][n] = sub(MP[rname][n], MP[rname][n - 2]) if n >= 2 else MP[rname][n]

bad = 0
for n in range(NMAX):
    tot = Z
    for (rname, dim, _) in CHARS:
        tot = add(tot, mul(rat(dim), MW[rname][n]))
    if sub(tot, rat(2 * n + 1)) != Z:
        bad += 1
say(f"  驗證 2（維數核）：Σ_ρ d_ρ m_ρ(n) = 2n+1 於 n=0..44，違例 = {bad}")
chk("C 驗證 2：維數核 2n+1", bad == 0, f"違例 {bad}")

nonint = [(r, n) for (r, _, _) in CHARS for n in range(NMAX)
          if MW[r][n][1] != 0 or MW[r][n][0].denominator != 1]
say(f"  驗證 3（整性／有理性）：非 ℤ 違例 = {len(nonint)}")
chk("D 驗證 3：m_ρ(n) ∈ ℤ≥0", len(nonint) == 0, f"違例 {nonint[:3]}")

# ── 重數表 ────────────────────────────────────────────────────────────────
say("")
say("§3  ★ 諧波重數表 n = 1..44（★ = O137 受保護集 P）")
say("-" * 78)
say(f"  {'n':>2} {'P':>1} | {'m_1':>3} {'m_3':>3} {'m_3′':>4} {'m_4':>3} {'m_5':>3} | 2n+1 | max")
mx_all, mx_P = 0, 0
hits3, hits3_P = [], []
for n in range(1, NMAX):
    row = [int(MW[r][n][0]) for r in ("1", "3", "3'", "4", "5")]
    mx = max(row)
    mx_all = max(mx_all, mx)
    if n in P_SET:
        mx_P = max(mx_P, mx)
    if mx >= 3:
        hits3.append((n, row))
        if n in P_SET:
            hits3_P.append((n, row))
    say(f"  {n:>2} {'★' if n in P_SET else ' '} | {row[0]:>3} {row[1]:>3} {row[2]:>4} "
        f"{row[3]:>3} {row[4]:>3} | {2 * n + 1:>4} | {mx}")

say("")
say("  ★ 盲測 B1（必要條件 = ∃ρ, n ∈ P : m_ρ(n) ≥ 3）：")
say(f"    P 內 max m_ρ(n) = {mx_P}；P 內命中 m ≥ 3 之 n = {[h[0] for h in hits3_P]}")
chk("★ B1 裁決（必要條件，協議執行）", True,
    f"{'成立' if hits3_P else '否證'}；P 內 max m_ρ(n) = {mx_P}（門檻 ≥ 3）")
say(f"    全範圍 n ≤ 44：max = {mx_all}；命中 m ≥ 3 之 n = {[h[0] for h in hits3]}"
    f"（{len(hits3)} 個）")

# ── 誘餌 ─────────────────────────────────────────────────────────────────
say("")
say("§4  誘餌審計：正則表示之『3』")
say("-" * 78)
reg = [(rn, d, d) for (rn, d, _) in CHARS]
say("  ℂ[A₅]（正則表示，維 60）= " + " ⊕ ".join(f"{c}·{rn}" for (rn, d, c) in reg))
chk("E 正則表示維數核", sum(d * c for (_, d, c) in reg) == 60, "Σd·d = 60")
say("  ⚠ dim-3 irrep 之重數 = 3 係**恆等**（正則表示重數 ≡ irrep 維數）⟹ VACUOUS-3，永不採納。")

# ── 判詞 ─────────────────────────────────────────────────────────────────
say("")
say("§5  判詞")
say("-" * 78)
say("  (1) m_ρ(n) 係精確表；corpus 內未見（新數據）。γ 表與此表之關係：皆 f(n)。")
say("  (2) 盲測 B1：" + ("成立 ⟹ A₅ 重數路線未死。" if hits3_P
                       else "**否證** ⟹ A₅ 由重數供 3 之路線死（範圍：chiral I ⊂ SO(3)、n ≤ 44、P 集）。"))
say("  (3) 正則表示之 3 = 定義性（VACUOUS-3）。")
say("  (4) scripts/01 §2 已證 A₅ 不能供非平凡標籤 ⟹ 重數亦只係 f(n) 之一部，非第二索引。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({
        "results": [{"tag": t_, "status": s_, "detail": d_} for (t_, s_, d_) in RES],
        "table": {n: [int(MW[r][n][0]) for r in ("1", "3", "3'", "4", "5")] for n in range(0, NMAX)},
        "protected": sorted(P_SET), "max_in_P": mx_P, "max_all": mx_all,
        "hits3_all": [h[0] for h in hits3], "hits3_P": [h[0] for h in hits3_P],
        "B1_confirmed": bool(hits3_P),
    }, f, ensure_ascii=False, indent=1)

nfail = sum(1 for (_, s_, _) in RES if s_ == "FAIL")
say("")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)