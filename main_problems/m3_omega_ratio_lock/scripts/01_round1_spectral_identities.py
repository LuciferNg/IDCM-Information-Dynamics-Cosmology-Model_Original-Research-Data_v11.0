#!/usr/bin/env python3
"""M3 scripts/01 — 第一輪：譜恆等式 + 分割枚舉 null (T1) + 結構剛性 (T2)

譜定義（權威）: 88-site dual chain，44-mode 半譜
    γ_n = 4 sin²(π n / 88),  n = 1..44          (等價 2−2cos(πn/44))
    88-mode 全譜: γ_k = 4 sin²(πk/176), k=1..88 (Σ=178, odd=88, even=90)
本輪全部用 44-mode（O215/O163 的預算約定）；兩約定之換算見 docs/01 §0。

跑法: /usr/bin/python3 scripts/01_round1_spectral_identities.py
輸出: stdout（正式記錄見 data/01_round1.out）
"""
from mpmath import mp, mpf, sin, pi
mp.dps = 50

PASS = FAIL = 0

def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:
        FAIL += 1; print(f"  [FAIL] {name}  {detail}")

# ── §0 譜量（44-mode） ────────────────────────────────────────────────
g = [4*sin(mpf(n)*pi/88)**2 for n in range(1, 45)]
S_all  = sum(g)
S_odd  = sum(g[0::2])      # n = 1,3,...,43  (22 項)
S_even = sum(g[1::2])      # n = 2,4,...,44  (22 項)
gmax   = g[-1]
B, N_d, N_gen, K_tree = 4, 44, 3, 33

print("§0 44-mode 譜量 (mpmath dps=50)")
chk("Σγ = 90",            abs(S_all - 90)  < mpf('1e-40'), f"|Δ|={mp.nstr(abs(S_all-90),4)}")
chk("Σγ_odd = 44 = N_d",  abs(S_odd - 44)  < mpf('1e-40'), f"|Δ|={mp.nstr(abs(S_odd-44),4)}")
chk("Σγ_even = 46",       abs(S_even - 46) < mpf('1e-40'), f"|Δ|={mp.nstr(abs(S_even-46),4)}")
chk("γ₄₄ = 4 = B (Nyquist)", abs(gmax - 4) < mpf('1e-40'), f"|Δ|={mp.nstr(abs(gmax-4),4)}")

# ── §0b 核心身份（EXACT） ────────────────────────────────────────────
print("\n§0b 60 份預算的身份鏈 (EXACT, |Δ|<1e-40)")
ident = [
    ("3  = γ_max − 1",              gmax - 1,                        N_gen),
    ("16 = γ_max²",                 gmax**2,                         B**2),
    ("41 = Σγ_odd − γ_max + 1",     S_odd - gmax + 1,                N_d - N_gen),
    ("60 = Σγ_odd + γ_max²",        S_odd + gmax**2,                 60),
    ("19 = 3 + 16",                 N_gen + B**2,                    19),
    ("123 = K_tree + Σγ",           K_tree + S_all,                  123),
]
for name, lhs, rhs in ident:
    chk(name, abs(lhs - rhs) < mpf('1e-40'), f"({mp.nstr(lhs,12)} vs {rhs})")

# ── §1 一般 M 定理（半譜，γ_n = 4sin²(πn/(2M))） ─────────────────────
print("\n§1 一般 M 定理：半譜求和（有理數精確，sympy）")
import sympy as sp
for m in (4, 6, 8, 10, 12, 20, 44):
    tot  = sp.simplify(sum(sp.sin(sp.pi*k/(2*m))**2 for k in range(1, m+1)))
    odd  = sp.simplify(sum(sp.sin(sp.pi*k/(2*m))**2 for k in range(1, m+1) if k % 2))
    even = sp.simplify(sum(sp.sin(sp.pi*k/(2*m))**2 for k in range(1, m+1) if not k % 2))
    chk(f"M={m:2d}: Σsin²=(M+1)/2, odd=M/4, even=M/4+½",
        tot == sp.Rational(m+1, 2) and odd == sp.Rational(m, 4)
        and even == sp.Rational(m, 4) + sp.Rational(1, 2),
        f"Σγ={4*tot}, Σγ_odd={4*odd}, Σγ_even={4*even}")
print("  ⟹ 定理: Σγ = 2M+2 = 2N_d+2 ; Σγ_odd = M = N_d ; Σγ_even = M+2 ; γ_M = 4")
print("  ⟹ 44-mode 之 Σγ_odd = N_d 係半譜定理（非 88 特例）")

# ── §2 T1: 分割枚舉 + χ²（O215 觀測協議） ────────────────────────────
print("\n§2 T1 分割枚舉 (a<b<c, a+b+c=60; O215 觀測協議)")
ops = sorted({1,2,3,4,5,7,8,9,11,13,16,17,18,21,22,25,29,33,41,44,88,90,123})
obs = (0.0493, 0.2645, 0.6847)      # Ω_b, Ω_dm, Ω_Λ (Planck 2018, O215 docs/04)
sig = (0.0003, 0.0050, 0.0074)

def chi2(t):
    return sum(((v/60 - o)/s)**2 for v, o, s in zip(t, obs, sig))

scored = sorted((chi2(t), t) for t in
                ((a, b, c) for a in ops for b in ops for c in ops
                 if a < b < c and a + b + c == 60))
ch2_id = chi2((3, 16, 41))
chk("3:16:41 喺 operator 字典內", (3, 16, 41) in [t for _, t in scored])
print(f"  總分割數 {len(scored)}；3:16:41 χ² = {ch2_id:.3f}")
print(f"  χ² ≤ χ²(3:16:41) 之分割: {sum(1 for c,_ in scored if c <= ch2_id)} / {len(scored)}")
print("  top-5:")
for chi, t in scored[:5]:
    z = [(v/60 - o)/s for v, o, s in zip(t, obs, sig)]
    print(f"    χ²={chi:7.2f}  {t[0]:>3}:{t[1]:>3}:{t[2]:>3}  z=({z[0]:+.2f},{z[1]:+.2f},{z[2]:+.2f})")
print(f"  rank(3:16:41) = {[t for _,t in scored].index((3,16,41))+1} / {len(scored)}")

# ── §2b T1b: 全整數字典對照（測試 χ² 本身是否可辨識） ────────────────
print("\n§2b T1b 全整數字典 1..59（χ² 可辨識性對照）")
allint = sorted((chi2(t), t) for t in
                ((a, b, c) for a in range(1,60) for b in range(1,60) for c in range(1,60)
                 if a < b < c and a + b + c == 60))
print(f"  總分割數 {len(allint)}；3:16:41 χ² = {ch2_id:.3f}")
print(f"  χ² ≤ χ²(3:16:41): {sum(1 for c,_ in allint if c <= ch2_id)} / {len(allint)}")
print(f"  rank = {[t for _,t in allint].index((3,16,41))+1} / {len(allint)}")
print("  ⟹ 全整數下仍唯一，但成因 = Ω_b σ 極細把 a 釘死 3；χ² 反映觀測，唔係結構證據")

# ── §3 T2: 結構剛性 ──────────────────────────────────────────────────
print("\n§3 T2 結構剛性")
print("  (a) (B−1)+B²+(N_d−B+1) ≡ N_d+B² 對所有 B 成立（代數恆等 → 無剛性）：")
for b in range(2, 9):
    print(f"      B={b}: 總和 {b-1+b*b+44-b+1} = N_d+B² {44+b*b}  ✓")
from math import ceil, log
print("  (b) κ=B⁻² 可見模式數 ⌈ln B²⌉ vs B−1：")
for b in range(2, 9):
    nv = ceil(log(b*b))
    print(f"      B={b}: ⌈ln B²⌉={nv}, B−1={b-1}  {'MATCH' if nv==b-1 else 'MISMATCH'}")
print("      ⟹ 只在 B=4,5 命中；IDCM 之 B 由 Nyquist γ_max=4 釘死 ⟹ 屬「條件恆等」，非普適定理")
print("  (c) Σγ_odd = N_d 為半譜定理（與 B 無關）⟹ 分母 60 的鏈條特異性全部來自 M=44")

# ── §4 60 身份普查（防雙重計數 |A₅|） ────────────────────────────────
print("\n§4 60 身份普查")
cands = {
    "N_d + B² = 44+16 (sector 份額之和)": 44 + 16,
    "|A₅| = 60 (icosahedral)":            60,
    "5!/2":                               120 // 2,
    "B·N_gen·F₅ = 4·3·5":                 4*3*5,
    "3×16+12":                            3*16+12,
    "2N_d − 28":                          2*44-28,
}
comp = "N_d + B² = 44+16 (sector 份額之和)"
for k, v in cands.items():
    tag = "★ 組合恆等（= 三份額之和）" if k == comp else "數值巧合（非組合來源）"
    print(f"  {k:38} = {v:3d}   {tag}")
print("  註: 組合恆等 (B−1)+B²+(N_d−B+1) ≡ N_d+B² 對任意 B 成立（§3a）")
print("      ⟹ 分母並非獨立假設，而係三份額之和；|A₅|=60 屬另一幾何來源，")
print("         若同時用作預算分母 = 雙重計數陷阱（O223 之 60 與此 60 必須分開引用）")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")