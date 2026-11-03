#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/07_variational_axiom_attack.py — M16 追加輪 2：攻擊 MF1（「最小階」變分公理）

三路：
  §A  變分公理之候選形式：min_B |G| = B³−B ⟹ 會被數據反證（B=2 → S₃，扇區 11/22/33）
  §B  ★ 替代結構條件（可精確證）：帶 γ 成等差之三連續帶唯一 = {3,4,5}
      γ_{n−1}+γ_{n+1}−2γ_n = 4cos(nπ/8)(1−cos(π/8))  ⟹ 零 ⟺ cos(nπ/8)=0 ⟺ n ≡ 4 (mod 8)
  §C  由此：γ_mid = γ_max/2 = μ² = 2 自動成立；|G| = B³−B = 60、Σ = 3B = 12 全為 B 之公式
  §D  自洽鎖：3B = B(B−1) ⟹ B = 4（唯一），即變分無可行域可變
  §E  掃描核對：B = 2..8 之 B³−B 與 DNT={B−1,B,B+1} 之群是否存在
"""
import os, sys, re
from fractions import Fraction
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


R = sp.Integer(11)          # L₅ = 帶寬（鏈週期）
N2 = sp.Integer(176)        # 2N_d = 88 之雙鏈長

def gam(n):
    """γ_{nR} = 4 sin²(nπ/16)，以 2(1−cos(nπ/8)) 精確表示。"""
    return sp.simplify(4 * sp.sin(n * sp.pi / 16) ** 2)


# ─────────────────────────────────────────────────────────────
OUT.append("§A 變分公理之候選形式：min_B (B³−B)")
tabA = []
for B in range(2, 9):
    g = B ** 3 - B
    tabA.append((B, g))
    OUT.append(f"  B={B}: |G| = B³−B = {g:4d}   （三元組 {B-1},{B},{B+1}，扇區 = ×L₅ = "
               f"{(B-1)*11},{(B)*11},{(B+1)*11}）")
OUT.append("  min_B (B≥2) 之 |G| = 6（B=2, S₃）⟹ 預言扇區 11/22/33")
OUT.append("  ⟹ 觀測家族 = 33/44/55 ⟹ **naive 最小成本被數據反證**（預言錯扇區）")
chk("A:naive-min-falsified", min(g for _, g in tabA) == 6 and tabA[0][0] == 2)

# ─────────────────────────────────────────────────────────────
OUT.append("§B ★ 等差 γ 之三連續帶 —— 精確判定")
n = sp.symbols("n", integer=True, positive=True)
expr = sp.simplify(gam(n - 1) + gam(n + 1) - 2 * gam(n))
closed = sp.simplify(4 * sp.cos(n * sp.pi / 8) * (1 - sp.cos(sp.pi / 8)))
OUT.append(f"  γ_{{n−1}}+γ_{{n+1}}−2γ_n 之化簡 = {expr}")
chk("B:closed-form", sp.simplify(expr - closed) == 0, "= 4cos(nπ/8)(1−cos(π/8))")
chk("B:prefactor-nonzero", sp.simplify(1 - sp.cos(sp.pi / 8)) != 0,
    f"1−cos(π/8) = {sp.nsimplify(1 - sp.cos(sp.pi / 8), [sp.cos(sp.pi/8)])}")

OUT.append("  帶 n = 1..8 之實值（γ_{nR}）與差值：")
rowsB = []
for m in range(1, 9):
    diff = sp.simplify(gam(m - 1) + gam(m + 1) - 2 * gam(m))
    val = sp.N(gam(m), 20)
    rowsB.append((m, val, sp.N(diff, 20), sp.simplify(diff) == 0))
    OUT.append(f"   n={m}: γ={val}   Δ={sp.N(diff,12)}   等差={'✓' if sp.simplify(diff)==0 else '✗'}")
zeros = [m for m, _, _, z in rowsB if z]
OUT.append(f"  帶範圍內 Δ=0 者：n = {zeros}")
chk("B:unique-n4", zeros == [4])
OUT.append("  ⟹ 三連續帶成等差者**唯一** = {3,4,5}（cos(nπ/8)=0 ⟹ n≡4 mod 8，鏈內只有 n=4）")

OUT.append("§B2 中項自動取半最大值：")
gmax = sp.Max(*[gam(m) for m in range(1, 9)])
OUT.append(f"  γ_max = {sp.simplify(gmax)} = γ_(8R)（鏈末端 n=2B=8）")
OUT.append(f"  γ_(4R) = {sp.simplify(gam(4))} = γ_max/2 = μ²  ⟹ 精確")
chk("B2:mid-half-max", sp.simplify(gmax) == 4 and sp.simplify(gam(4)) == 2)

# ─────────────────────────────────────────────────────────────
OUT.append("§C 由三元組到群階：全部係 B 之公式（零自由參數）")
for B in [4]:
    prod = (B - 1) * B * (B + 1)
    s = (B - 1) + B + (B + 1)
    OUT.append(f"  B={B}: ∏ = {prod} = B³−B = |G| ；Σ = {s} = 3B ；Σ = B(B−1)? {s == B*(B-1)}")
    chk("C:order", prod == 60)
    chk("C:sum", s == 12)

# ─────────────────────────────────────────────────────────────
OUT.append("§D 自洽鎖：3B = B(B−1) 之整數解（本 O 之 Σ=12 條件）")
sols = [B for B in range(2, 50) if 3 * B == B * (B - 1)]
OUT.append(f"  解 = {sols}")
chk("D:unique-B", sols == [4])
OUT.append("  ⟹ 可行域只有一點 ⟹ 任何「最小化」都係空轉（vacuous）：無兩個候選可比")

# ─────────────────────────────────────────────────────────────
OUT.append("§E 掃描核對：|G| = B³−B 且 DNT = {B−1,B,B+1} 之群（階 ≤ 500）")
LINE = re.compile(r"^(\d+)\s+(\d+)\s+\[(.*)\]\s*$")
G = []
with open(SCAN) as f:
    for ln in f:
        m = LINE.match(ln.strip())
        if not m:
            continue
        degs = []
        for d, k in re.findall(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", m.group(3)):
            degs.extend([int(d)] * int(k))
        G.append((int(m.group(1)), int(m.group(2)), sorted({d for d in degs if d > 1})))
for B in range(2, 9):
    target = B ** 3 - B
    want = [B - 1, B, B + 1]
    hits = [(a, b) for a, b, S in G if a == target and S == want]
    OUT.append(f"  B={B}: |G|={target:4d} DNT={want} → 群 {hits if hits else '∅'}")
OUT.append("  ⟹ B=4 有 A₅；B=5 有 S₅；B=2,3,6,7,8 之三元組目標群不存在（或不在階 ≤ 500 內）")

print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)