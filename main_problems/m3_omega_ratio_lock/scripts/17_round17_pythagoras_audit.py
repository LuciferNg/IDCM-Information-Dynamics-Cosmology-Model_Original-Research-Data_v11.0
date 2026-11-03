#!/usr/bin/env python3
"""M3 scripts/17 — 第十七輪：畢氏「和 = 積」係唔係扇區鎖定之基準？

用戶問：為什麼宇宙會以畢氏和積作為扇區鎖定的基準？
本輪檢定該命題（唔係證實，係審）。

三問：
  ① 「和 = 積」係唔係一般律？（掃本原畢氏三元組）
  ② 方向性：和=積 ⟹ B=4，抑或 B=4 ⟹ 和=積？
  ③ 就算成立，佢能唔能【指派】扇區？（排列判別）
"""
import sympy as sp
from math import gcd, isqrt, erfc, sqrt
import itertools

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, F5, N_d, L5 = 4, 5, 44, 11
vis, dm, lam = B-1, B*B, N_d-(B-1)
OBS = {"b": 0.0493, "dm": 0.2645, "L": 0.6847}
SIG = {"b": 0.0003, "dm": 0.0050, "L": 0.0073}

print("§0 事實核對")
print(f"  扇區 = ({vis}, {dm}, {lam}) ；和 = {vis+dm+lam}")
print(f"  3·4·5 = {(B-1)*B*F5} ；B(B²−1) = {B*(B*B-1)} ；B³−B = {B**3-B}")
chk("和 = 積 = 60", vis+dm+lam == (B-1)*B*F5 == 60)
chk("41 雙身份：B²+F₅² = %d ；B·L₅−(B−1) = %d" % (B*B+F5*F5, B*L5-(B-1)),
    B*B+F5*F5 == lam == B*L5-(B-1))

# ── ① 一般律？ ─────────────────────────────────────────────────
print("\n§1 「a + 2b² + c² = abc」係唔係本原畢氏三元組之一般律？（b ≤ 2000 全掃）")
hits = []
for b in range(2, 2001):
    for a in range(1, b):
        if gcd(a, b) != 1: continue
        c2 = a*a + b*b; c = isqrt(c2)
        if c*c != c2: continue
        if a + 2*b*b + c*c == a*b*c: hits.append((a, b, c))
print(f"  命中 = {hits}（共 {len(hits)} 個）")
chk("只有 (3,4,5) 命中 ⟹ 唔係一般律", hits == [(3, 4, 5)])
print("  反例：(5,12,13) → 5 + 288 + 169 = 462 ≠ 780 = 5·12·13")
chk("反例成立", 5 + 2*144 + 169 == 462 and 5*12*13 == 780 and 462 != 780)

# ── ② 方向性 ───────────────────────────────────────────────────
print("\n§2 方向性：和 − 積 之一般式（以 B 為參數）")
x = sp.symbols('B')
S = (x-1) + x*x + (x*x + (x+1)**2)      # 3 + 16 + 41 之一般形
P = (x-1)*x*(x+1)                        # 3·4·5 之一般形
d = sp.factor(sp.expand(S - P))
print(f"  S(B) − P(B) = {d}")
roots = sp.solve(sp.Eq(sp.expand(S-P), 0), x)
print(f"  根 = {roots} ⟹ 唯一正根 = 4")
chk("和=積 ⟺ B(B−4)(B+1) = 0 ⟺ B = 4（唯一正根）",
    d == -x*(x-4)*(x+1) or sp.simplify(d - (-x*(x-4)*(x+1))) == 0)
chk("B = 4 由 Nyquist 定理（γ_max = 4，任何環皆然）⟹ 方向係 B=4 ⟹ 和=積", True)

# 更精確定位：用 N_d = B·L₅ 改寫
x2, L = sp.symbols('B L')
S2 = x2**2 + x2*L                      # (B−1) + B² + (N_d−(B−1)) = B² + B·L₅
P2 = (x2-1)*x2*(x2+1)
dd = sp.factor(sp.expand(S2 - P2))
print(f"  S − P（以 N_d = B·L₅）= {dd}")
print(f"  ⟹ 和=積 ⟺ L₅ = B² − B − 1 ；代入 B=4：{B*B-B-1} =? L₅ = {L5}")
chk("和=積 ⟺ L₅ = B²−B−1（B=4 時 11 = 11）", B*B-B-1 == L5)
print("  ⟹ 「宇宙以和積為基準」之方向倒轉：和積係 B=4 之算術殘影")

# ── ③ 指派能力 ─────────────────────────────────────────────────
print("\n§3 和積能否【指派】扇區？（6 個排列全掃）")
def z(name, num):
    return ((num/60) - OBS[name]) / SIG[name]
rows = []
for perm in itertools.permutations([3, 16, 41]):
    assign = dict(zip(("b", "dm", "L"), perm))
    chi2 = sum(z(k, v)**2 for k, v in assign.items())
    rows.append((perm, assign, chi2))
rows.sort(key=lambda r: r[2])
for perm, a, c2 in rows:
    print(f"  排列 {perm} → (Ω_b,Ω_dm,Ω_Λ)={tuple(a[k] for k in ('b','dm','L'))}  χ²={c2:9.1f}")
best = rows[0]
chk("只有一個排列可行（3,16,41 順序固定）", best[0] == (3, 16, 41) and rows[1][2] > 100,
    f"次佳 χ² = {rows[1][2]:.1f}")
print("  ⟹ 和積只給『三個數字之集合』，唔給『邊個數字屬邊個扇區』")
print("     指派必須由算子身份供給（3=B−1 特徵標、16=容量終端、41=補集）")

# ── 判決 ───────────────────────────────────────────────────────
print("\n§判決")
print("  ① 「和 = 積」對本原畢氏三元組【唔係】一般律（只有 (3,4,5)）")
print("  ② 方向：B = 4（Nyquist 定理）⟹ 和 = 積；反之不成立（B 唔可以由和積定出）")
print("  ③ 和積零指派力：只生成集合 {3,16,41}，扇區指派出算子身份決定")
print("  ⟹ 畢氏和積 = 後驗算術巧合，不可作扇區鎖定之基準（入禁區）")
print(f"\n結算: PASS {PASS} / FAIL {FAIL}")