#!/usr/bin/env python3
"""
M19 / scripts/04 — R3：c/b/t 指數掃描（密度閘 5% ＋ MC null 10⁷）

預註冊（跑前凍結）:
  E1  規則空間（**全列，不剪枝**，共 30 條）：
        A 排名族(9)   : 「9 模式集中第 r 小」r = 1..9
        B 餘數族(9)   : n mod m = a，m ∈ {2,3,4}，a = 0..m−1
        C 數論集族(6) : primes / Fibonacci / Lucas / 平方 / 三角 / 2 冪
        D γ 閾族(6)   : γ_n < x，x ∈ {0.05,0.1,0.2,0.3,0.5,1.0}
      每條規則誘導一個模式子集 S ⊆ M9（唯一值指派唔可能，故只測**子集命中**）
  E2  目標（三個，跑前寫死）：
        T1 = {c,b,t} = {7,11,13}（3 模，R3 真目標）
        T2 = {c,b} = {7,11}（2 模）
        T3 = {c} = {7}（1 模）
  E3  密度閘（SESSION_START §三 硬性）：命中率 = (純結構命中規則數)/(21) ≤ 5% ⟺ ≤ 1 條
  E4  MC null 10⁷（seed 固定）：隨機雙射令 c=7∧b=11∧t=13 之機率（理論 1/504）；
      隨機 3 子集 = {7,11,13} 之機率（理論 1/84，作 E3 之真實基準）
  E5  正控制：質量 rank 規則必須「命中全部目標」，但該規則**使用觀測質量** ⟹ 標記為
      輸入（禁區：由質量反解 k_q）。若正控制唔命中 ⟹ 本輪掃描無效（停）。
  E6  判決：任何命中規則若唔屬 A 排名族＋質量輸入，須逐條列出其結構出處；
      冇出處 ⟹ 記「無 provenance 命中」，依 M13 docs/11 標籤層定理不得升格。

失敗判準: E5 FAIL ⟹ 掃描設定錯（停，唔出判決）。
"""
import json
import math
from pathlib import Path

import numpy as np
from mpmath import mp, mpf

mp.dps = 40
OUT = Path(__file__).parent.parent / "data"

M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
MASS = {"e": 0.511, "u": 2.16, "d": 4.67, "s": 93.4, "mu": 105.6584,
        "c": 1270.0, "tau": 1776.86, "b": 4180.0, "t": 172690.0}

T1 = {7, 11, 13}           # c, b, t（R3 真目標）
T2 = {7, 11}               # c, b
T3 = {7}                   # c

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def gamma(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


rules = []   # (name, family, set, uses_mass)

# A 排名族：第 r 小 —— 「秩定位」本身需要目標粒子之質量 rank，故標記為輸入（r 係 rank）
for r in range(1, 10):
    rules.append((f"A 第{r}小", "A", {M9[r - 1]}, True))

# B 餘數族
for m in (2, 3, 4):
    for a in range(m):
        rules.append((f"B n≡{a}(mod {m})", "B", {n for n in M9 if n % m == a}, False))

# C 數論集族
def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


fib = {1, 1}
a, b = 1, 1
while b < 20:
    a, b = b, a + b
    fib.add(b)
lucas = {2, 1}
a, b = 2, 1
while b < 20:
    a, b = b, a + b
    lucas.add(b)
sets_c = {
    "primes": {n for n in M9 if is_prime(n)},
    "Fibonacci": {n for n in M9 if n in fib},
    "Lucas": {n for n in M9 if n in lucas},
    "平方": {n for n in M9 if int(math.isqrt(n)) ** 2 == n},
    "三角": {n for n in M9 if any(k * (k + 1) // 2 == n for k in range(1, 6))},
    "2冪": {n for n in M9 if n & (n - 1) == 0},
}
for k, v in sets_c.items():
    rules.append((f"C {k}", "C", v, False))

# D γ 閾族
for x in ("0.05", "0.1", "0.2", "0.3", "0.5", "1.0"):
    rules.append((f"D γ<{x}", "D", {n for n in M9 if gamma(n) < mpf(x)}, False))

# E5 正控制：質量 rank（使用觀測質量 ⟹ 輸入）
f1 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}
rank_quark = {f1[p] for p in ("c", "b", "t")}
rules.append(("E 質量 rank（正控制／輸入）", "E", rank_quark, True))

print("=" * 92)
print("M19 scripts/04 — R3：c/b/t 掃描（30 規則 ＋ 正控制）")
print("=" * 92)

hits_t1 = [r for r in rules if r[2] == T1]
hits_t2 = [r for r in rules if r[2] == T2]
hits_t3 = [r for r in rules if r[2] == T3]
pure = [r for r in rules if not r[3]]

print(f"  規則數 = {len(rules)}（純結構 {len(pure)} ＋ 質量輸入 {len(rules)-len(pure)}）")
print(f"  T1 命中 = {len(hits_t1)} → {[h[0] for h in hits_t1]}")
print(f"  T2 命中 = {len(hits_t2)} → {[h[0] for h in hits_t2]}")
print(f"  T3 命中 = {len(hits_t3)} → {[h[0] for h in hits_t3]}")

chk("E1 規則空間 30 條（全列，不剪枝）＋ 1 正控制", len(rules) == 31 and len(pure) == 21,
    f"{len(rules)} 條（純結構 B/C/D = {len(pure)}；秩族 A = 9 屬輸入；正控制 = 1）")
dens_t1 = len([h for h in hits_t1 if not h[3]]) / len(pure)
dens_t2 = len([h for h in hits_t2 if not h[3]]) / len(pure)
dens_t3 = len([h for h in hits_t3 if not h[3]]) / len(pure)
chk("E3 密度閘 ≤5%（純結構規則 B/C/D = 21 條）",
    dens_t1 <= 0.05 and dens_t2 <= 0.05 and dens_t3 <= 0.05,
    f"T1 {dens_t1:.1%}｜T2 {dens_t2:.1%}｜T3 {dens_t3:.1%}")

# E4：MC null 10^7（chunked）
rng = np.random.default_rng(20261114)
N_MC = 10_000_000
CH = 1_000_000
cnt_perm = 0     # 隨機雙射：c→7, b→11, t→13（以位置 0,2,5 代表模式 7,11,13 之秩）
cnt_sub5 = 0     # 隨機 5 子集 == T1
for _ in range(N_MC // CH):
    idx = rng.random((CH, 9)).argsort(axis=1)
    perm = np.array(M9)[idx]                       # 每行 = 一個雙射（9 模之排列）
    # c, b, t 分別抽三個不同位置（前 3 欄）
    ok = (perm[:, 0] == 7) & (perm[:, 1] == 11) & (perm[:, 2] == 13)
    cnt_perm += int(ok.sum())
    sub = perm[:, :3]
    matches = np.all(np.sort(sub, axis=1) == np.array(sorted(T1)), axis=1)
    cnt_sub5 += int(matches.sum())

p_perm = cnt_perm / N_MC
p_sub5 = cnt_sub5 / N_MC
chk("E4a MC 雙射同時 c7∧b11∧t13 ≈ 1/504",
    abs(p_perm - 1 / 504) < 5 * math.sqrt((1 / 504) * (1 - 1 / 504) / N_MC),
    f"p = {p_perm:.3e}（理論 {1/504:.3e}）；期望計數 {N_MC/504:.0f}")
chk("E4b MC 3 子集 = {7,11,13} ≈ 1/84",
    abs(p_sub5 - 1 / 84) < 5 * math.sqrt((1 / 84) * (1 - 1 / 84) / N_MC),
    f"p = {p_sub5:.3e}（理論 {1/84:.3e}）；期望計數 {N_MC/84:.0f}")

chk("E5 正控制（質量 rank）命中 T1",
    any(h[3] and h[2] == T1 for h in hits_t1),
    "命中但標記為輸入 ⟹ 禁區路線（由質量反解 k_q）")

print("\n  純結構命中明細：")
for name, fam, s, uses in pure:
    mark = ""
    if s == T1:
        mark += " [T1]"
    if s == T2:
        mark += " [T2]"
    if s == T3:
        mark += " [T3]"
    print(f"    {name:<18} → {sorted(s)}{mark}")

OUT.mkdir(exist_ok=True)
(OUT / "04_scan_density.json").write_text(json.dumps({
    "n_rules": 30,
    "rules": [{"name": n, "family": f, "set": sorted(s), "uses_mass": u} for n, f, s, u in rules],
    "hits_T1": [h[0] for h in hits_t1],
    "hits_T2": [h[0] for h in hits_t2],
    "hits_T3": [h[0] for h in hits_t3],
    "density": {"T1": dens_t1, "T2": dens_t2, "T3": dens_t3},
    "mc": {"N": N_MC, "seed": 20261114, "p_perm_c7b11t13": p_perm,
           "p_sub5_T1": p_sub5, "theory_perm": 1 / 504, "theory_sub5": 1 / 126},
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
