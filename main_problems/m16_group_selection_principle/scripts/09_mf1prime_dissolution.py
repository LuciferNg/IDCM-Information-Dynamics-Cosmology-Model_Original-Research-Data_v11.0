#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/09_mf1prime_dissolution.py — M16 追加輪 4：攻 MF1′（「三帶連續 / 等 γ 間距」是否為必要前提？）

結論預告：**唔必要**。{3,4,5} 已由「Σ = 12」＋ Burnside（L2）唯一逼出，
「等距（AP）」與「三帶連續」皆係事後自動成立之性質 ⟹ MF1′ 撤回。

§A 窮盡 Σ = 12 之相異整數集（任何大小）＋ 非平凡 ⟹ 8 集
§B 套 L2（∏ > Σd²） ⟹ 唯一 {3,4,5}
§C 帶號連續性與 γ-AP 係「結果」而非「前提」：對 {3,4,5} 自動成立；對被殺者無意義
§D 群層重核（92,570 群）：唯一 A₅
§E look-elsewhere 記賬：56 → 8 → 1
"""
import os, sys, re, itertools
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCAN = os.path.join(ROOT, "data", "run_02_scan_all_500.txt")

NPASS = NFAIL = 0
OUT = []


def chk(n, c, d=""):
    global NPASS, NFAIL
    if c:
        NPASS += 1
        OUT.append(f"PASS {n}" + (f"  [{d}]" if d else ""))
    else:
        NFAIL += 1
        OUT.append(f"FAIL {n}" + (f"  [{d}]" if d else ""))


# ── §A Σ=12 之全部相異整數集（元素 ≥ 2，任何大小）
OUT.append("§A 窮盡 Σ = 12 之相異整數集（元素 ≥ 2，任何大小；「非平凡」已內建）")
CAND = []
for r in range(1, 13):
    for c in itertools.combinations(range(2, 13), r):
        if sum(c) == 12:
            CAND.append(c)
OUT.append(f"  共 {len(CAND)} 集：{CAND}")
chk("A:count-8", len(CAND) == 8)

# ── §B L2：∏ > Σd²
OUT.append("§B L2（Burnside 必要條件：∏_{d} d > Σ_{d} d²，且 |G| = ∏）：")
SURV = []
for c in CAND:
    prod = 1
    for d in c:
        prod *= d
    sq = sum(d * d for d in c)
    ok = prod > sq
    OUT.append(f"  {str(list(c)):<12} ∏={prod:3d}  Σd²={sq:3d}  {'✓ 生還' if ok else '✗'}")
    if ok:
        SURV.append(c)
OUT.append(f"  ⟹ 生還 = {SURV}")
chk("B:unique-survivor", SURV == [(3, 4, 5)])

# ── §C 連續性 / γ-AP 係結果
OUT.append("§C 「三帶連續」與「γ-AP」係結果，唔係前提")
gam = lambda m: sp.simplify(2 * (1 - sp.cos(m * sp.pi / 8)))
g = [3, 4, 5]
cons = (g[1] - g[0] == 1) and (g[2] - g[1] == 1)
ap = sp.simplify(gam(g[0]) + gam(g[2]) - 2 * gam(g[1])) == 0
OUT.append(f"  {g}: 連續={cons}；γ_3+γ_5−2γ_4 = {sp.simplify(gam(3)+gam(5)-2*gam(4))} ⟹ AP={ap}")
chk("C:cons-and-ap-hold", cons and ap)
OUT.append("  註：另 2 組 AP 三元組之死亡路徑各異 ——")
OUT.append("      (1,4,7)：含帶 1 ⟹ 度數 = 1 = trivial ⟹ 由「非平凡」剔除（唔入 CAND）")
OUT.append("      (2,4,6)：入 CAND（Σ=12）但 ∏=48 < Σd²=56 ⟹ 由 L2 殺")
chk("C:ap-others-killed", (2, 4, 6) in CAND and (1, 4, 7) not in CAND,
    "(2,4,6) 由 L2 殺；(1,4,7) 由非平凡殺")

# ── §D 群層重核
OUT.append("§D 群層重核（階 ≤ 500 之 92,570 群）")
LINE = re.compile(r"^(\d+)\s+(\d+)\s+\[(.*)\]\s*$")
hits = []
N = 0
for ln in open(SCAN):
    m = LINE.match(ln.strip())
    if not m:
        continue
    N += 1
    degs = []
    for d, k in re.findall(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", m.group(3)):
        degs.extend([int(d)] * int(k))
    S = sorted({d for d in degs if d > 1})
    if not S:
        continue
    prod = 1
    for d in S:
        prod *= d
    if sum(S) == 12 and prod == int(m.group(1)):
        hits.append((int(m.group(1)), int(m.group(2)), S))
OUT.append(f"  掃描群數 = {N}；命中（Σ=12 ∧ ∏=|G|）= {hits}")
chk("D:unique-A5", hits == [(60, 5, [3, 4, 5])])

# ── §E look-elsewhere 記賬
OUT.append("§E look-elsewhere 記賬（表面扇區選擇之自由度）")
n_all = len(list(itertools.combinations(range(1, 9), 3)))
OUT.append(f"  ① 無任何約束：8 帶取 3 = C(8,3) = {n_all} 組 ⟹ 1/{n_all} = {1/n_all:.4f}（O171 docs/09 R1b 之 1.79%）")
OUT.append(f"  ② 加 Σ = 12（12 = B(B−1) 已導出，M15 docs/17）：{len(CAND)} 集 ⟹ 1/{len(CAND)} = {1/len(CAND):.4f}")
OUT.append(f"  ③ 再加 L2（Burnside 必要條件）：1 集 ⟹ 自由度 0")
OUT.append("  ⟹ 理論側剩餘代價全部落在「Σ = 12」呢一個數字；而 c 掃描（docs/02 §4.1）顯示")
OUT.append("     c = 1..40 中只有 c = 12 會唯一命中 A₅ ⟹ 該數字已被測過，唔係新自由度。")
chk("E:cascade", 1 / n_all < 1 / len(CAND))

print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)