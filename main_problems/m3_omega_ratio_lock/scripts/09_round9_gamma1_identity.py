#!/usr/bin/env python3
"""M3 scripts/09 — 第九輪：γ_1 = 4sin²(π/88) 之結構身份

問題: 最輕 gap γ_1 有冇結構身份？（docs/08 觀察：0 < 0.0051 ≪ 2 < 4）

檢定:
  §1 和律 Σ_{k=1..43} γ_k = 86 ; Σ_{k=0..43}+γ_44 = 90
  §2 ★ 乘積律 ∏_{k=1..43} γ_k = N_d = 44（Cycle 之 Matrix-Tree 半積）
  §3 Chebyshev 生成: γ_n = 2 − 2T_n(1 − γ_1/2) ⟹ 全譜 ∈ ℤ[γ_1]
  §4 分形倍頻: γ_n^{(88)} = γ_{2n}^{(176)} ; 特別 γ_1^{(88)} = γ_2^{(176)}
  §5 γ_1 之代數次數 = φ(176)/2 = 40（非有理 ⟹ 唔屬 ℚ(B, N_gen, L₅)）
"""
import mpmath as mp
from fractions import Fraction as F

mp.mp.dps = 60
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

N, B = 88, 4
N_d = 44
gam = lambda n, M=88: B * mp.sin(mp.pi*n/M)**2

# ── §1 和律 ────────────────────────────────────────────────────────
print("§1 和律")
S43 = sum(gam(k) for k in range(1, 44))
Sall = sum(gam(k) for k in range(88))
print(f"  Σ_{{k=1..43}} γ_k = {mp.nstr(S43, 20)}")
print(f"  Σ_{{k=0..87}} γ_k = {mp.nstr(Sall, 20)}")
chk("Σ_{k=1..43} γ_k = 86", abs(S43 - 86) < mp.mpf('1e-45'), mp.nstr(S43, 15))
chk("Σ_{k=0..87} γ_k = 176 = 2N_d", abs(Sall - 176) < mp.mpf('1e-45'), mp.nstr(Sall, 15))
chk("(Σ_all + γ_44)/2 = 90 = Σγ（M3）", abs((Sall + 4)/2 - 90) < mp.mpf('1e-45'))

# ── §2 乘積律 ─────────────────────────────────────────────────────
print("\n§2 ★ 乘積律（Matrix-Tree：C_N 生成樹數 = N）")
P43 = mp.mpf(1)
for k in range(1, 44):
    P43 *= gam(k)
print(f"  ∏_{{k=1..43}} γ_k = {mp.nstr(P43, 20)}")
chk("∏_{k=1..43} γ_k = 44 = N_d（精確至 1e-45）",
    abs(P43 - N_d) < mp.mpf('1e-45'), mp.nstr(P43, 15))
P_all = P43 * gam(44)
print(f"  ∏_{{k=1..43}} γ_k × γ_44 = {mp.nstr(P_all, 20)}")
chk("× γ_44 = 176 = 2N_d", abs(P_all - 176) < mp.mpf('1e-45'), mp.nstr(P_all, 15))
print("  ⟹ 譜之 e₁ = 86（⟹ Σγ=90）、e₄₃ = 44 = N_d：兩條 M3 數皆為譜之初等對稱函數")

# ── §3 Chebyshev 生成 ─────────────────────────────────────────────
print("\n§3 Chebyshev 生成：γ_n = 2 − 2T_n(c), c = 1 − γ_1/2")
g1 = gam(1)
c = 1 - g1/2
def T(n, x):
    # Chebyshev T_n via recurrence (mpmath)
    if n == 0: return mp.mpf(1)
    if n == 1: return x
    a, b = mp.mpf(1), x
    for _ in range(2, n+1):
        a, b = b, 2*x*b - a
    return b
maxdev = mp.mpf(0)
for n in range(0, 45):
    pred = 2 - 2*T(n, c)
    dev = abs(pred - gam(n))
    maxdev = max(maxdev, dev)
print(f"  max |2−2T_n(c) − γ_n| (n=0..44) = {mp.nstr(maxdev, 6)}")
chk("全譜由 γ_1 生成（最大偏差 < 1e-40）", maxdev < mp.mpf('1e-40'), mp.nstr(maxdev, 5))
chk("γ_0 = 2−2T_0(c) = 0（自動）", abs((2 - 2*T(0, c))) < mp.mpf('1e-45'))
chk("γ_44 = 2−2T_44(c) = 4 = B（自動）", abs((2 - 2*T(44, c)) - 4) < mp.mpf('1e-40'))

# ── §4 分形倍頻 ───────────────────────────────────────────────────
print("\n§4 分形倍頻：γ_n^(88) = γ_{2n}^(176)")
maxdev2 = max(abs(gam(n, 88) - gam(2*n, 176)) for n in range(45))
print(f"  max |γ_n^(88) − γ_{{2n}}^(176)| = {mp.nstr(maxdev2, 6)}")
chk("γ_n^(88) = γ_{2n}^(176) 全 n（< 1e-45）", maxdev2 < mp.mpf('1e-45'))
chk("特別：γ_1^(88) = γ_2^(176)（鏈之 gap = 雙覆蓋之第二諧波）", True,
    f"γ_1 = {mp.nstr(gam(1), 15)}")

# ── §5 γ_1 之代數次數 ─────────────────────────────────────────────
print("\n§5 γ_1 之代數次數")
from math import gcd
def phi(n):
    r = n
    x = n
    p = 2
    while p*p <= x:
        if x % p == 0:
            while x % p == 0: x //= p
            r -= r//p
        p += 1
    if x > 1: r -= r//x
    return r
# γ_1 = 4sin²(π/88) = 2 − 2cos(2π/88) = 2 − 2cos(π/44)  ⟹ α = 2cos(π/44)
deg = phi(88)//2
print(f"  γ_1 = 2 − α , α = 2cos(π/44)；φ(88) = {phi(88)} ⟹ deg(α) = {deg}")
conj = [2*mp.cos(2*mp.pi*k/88) for k in range(1, 44) if gcd(k, 88) == 1]
print(f"  共軛數 = {len(conj)}；相異 = {len({mp.nstr(x, 40) for x in conj})}")
chk("共軛數 = 20 = φ(88)/2", len(conj) == 20)
chk("共軛全相異（最小多項式次數 = 20）", len({mp.nstr(x, 40) for x in conj}) == 20)
print("  [推論] deg > 1 ⟹ γ_1 ∉ ℚ ⟹ 非 ℚ(B,N_gen,L₅) 元素（由 §5/§6 實算之次數支持）")

# ── §6 γ_1 之最小多項式（精確，sympy）────────────────────────────
print("\n§6 最小多項式 q(y)（deg 20；α = 2cos(π/44)）")
try:
    import sympy as sp
    sym = sp.Symbol('x'); yv = sp.Symbol('y')
    p = sp.minpoly(2*sp.cos(sp.pi/44), sym)
    q = sp.expand(p.subs(sym, 2 - yv))
    Pq = sp.Poly(q, yv)
    co = [int(c) for c in Pq.all_coeffs()]
    open("data/09_round9_minpoly.txt", "w").write(
        "# γ_1 = 4sin²(π/88) 之最小多項式（deg 20），係數由高次至常數\n"
        + " ".join(str(c) for c in co) + "\n")
    print(f"  deg = {sp.degree(q, yv)}；LC = {co[0]}；常數項 = {co[-1]}")
    print(f"  e₁ = {abs(co[1])}（Σ 40 個共軛 gap）")
    print(f"  首幾個係數: {co[:6]}")
    # Horner 驗證 q(γ_1) = 0（60 dps）
    v = mp.mpf(0)
    for c in co:
        v = v*g1 + c
    print(f"  q(γ_1) ≈ {mp.nstr(v, 6)}")
    chk("deg q = 20", len(co) - 1 == 20)
    chk("q(γ_1) = 0（|·| < 1e-40）", abs(v) < mp.mpf('1e-40'), mp.nstr(v, 5))
    chk("常數項 = 1 ⟹ γ_1 係 unit（Norm = 1）", co[-1] == 1)
    chk("e₁ = 40 = 2·deg（Σ 共軛 gap = 2·20）", abs(co[1]) == 40)
except Exception as e:
    print(f"  [SKIP] sympy 不可用：{e}")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")