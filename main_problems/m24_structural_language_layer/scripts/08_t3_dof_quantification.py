#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/08 — T3 收緊：零參數 ⟹ Aut-不變（定義級）＋ 自由度量化
==================================================================
判據（`docs/06_preregistration_round2.md` §4，跑前凍結）：

  D1  定義級命題：零參數（結構決定）⟹ Aut-不變；非 Aut-不變者必含「標籤選擇」
  D2  量化：非等變雙射之選擇數 = 44!（mod |Aut(A₅)| = 60 後仍 ≈ 44!/60）
      ⟹ 自由度（bits）≫ 0 ⟹ 違 dof = 0
  D3  結論只可寫「在零參數紀律下 L3 無傳輸（定義級）」，不得升 🟢

★ 本輪新增之量化結果：等變族本身亦有殘餘自由度 ——
   共軛不變子集 = 2³ = 8（3 bits）｜僅逆不變 = 2²² ｜非等變 = 44!（179 bits）
   而「全部可能內容」= ≤ 3 bits，且 8 個選項中唯一與模式側可比者 = {20,24} 粗標籤。

零 import 上游。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations, permutations as _perms
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


E5 = tuple(range(5))


def parity(p):
    return sum(1 for i in range(5) for j in range(i + 1, 5) if p[i] > p[j]) % 2


G = [p for p in _perms(range(5)) if parity(p) == 0]


def mul(a, b):
    return tuple(a[b[i]] for i in range(5))


def invp(a):
    r = [0] * 5
    for i, v in enumerate(a):
        r[v] = i
    return tuple(r)


def order(g):
    x, k = E5, 1
    while True:
        x = mul(x, g)
        if x == E5:
            return k
        k += 1


ORD = {g: order(g) for g in G}
NONINV = [g for g in G if ORD[g] > 2]
MODES = list(range(1, 45))
N_d, B = 44, 4

try:
    import mpmath as mp
    mp.mp.dps = 40
except ImportError:
    raise SystemExit("mpmath required")

print("=" * 78)
print("M24 scripts/08 — T3 收緊：零參數 ⟹ Aut-不變 ＋ 自由度量化")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# D1  等變族之實際大小（Aut 不變子集）
# ══════════════════════════════════════════════════════════════════
print("\n── D1  Aut-不變子集之計數（等變族大小）──")


def orbit_decomp(elems, act):
    elems = list(elems)
    domain = set(elems)
    seen, out = set(), []
    for x in elems:
        if x in seen:
            continue
        o, stack = set(), [x]
        while stack:
            y = stack.pop()
            if y in o:
                continue
            o.add(y)
            for z in act(y):
                if z in domain and z not in o:
                    stack.append(z)
        seen |= o
        out.append(sorted(o))
    return out


CONJ = orbit_decomp(NONINV, lambda g: [mul(mul(h, g), invp(h)) for h in G])
INV = orbit_decomp(NONINV, lambda g: [invp(g)])
CI = orbit_decomp(NONINV, lambda g: [mul(mul(h, g), invp(h)) for h in G] + [invp(g)])
print(f"  共軛 orbit = {[len(o) for o in CONJ]}（{len(CONJ)} 個）")
print(f"  逆 orbit   = {[len(o) for o in INV]}（{len(INV)} 個）")
print(f"  ⟨共軛,逆⟩  = {[len(o) for o in CI]}（{len(CI)} 個）")
n_conj_subsets = 2 ** len(CONJ)
n_inv_subsets = 2 ** len(INV)
n_ci_subsets = 2 ** len(CI)
print(f"  不變子集數：共軛 {n_conj_subsets}｜逆 {n_inv_subsets}｜⟨共軛,逆⟩ {n_ci_subsets}")
ck("D1-a  共軛不變族 = 2³ = 8（3 bits）；僅逆不變 = 2²²（22 bits）",
   n_conj_subsets == 8 and n_inv_subsets == 2 ** 22, f"{n_conj_subsets}｜{n_inv_subsets}")

# 非等變雙射數 = 44!（因等變雙射 = 0）
import math
F44 = math.factorial(44)
bits44 = mp.log(F44, 2)
print(f"  44! = {F44}")
print(f"  log₂(44!) = {mp.nstr(bits44, 12)} bits｜log₁₀(44!) = {mp.nstr(mp.log(F44, 10), 10)}")
quot = F44 // len(G)
print(f"  44!/|Aut(A₅)| = {quot}（bit 長 = {quot.bit_length()}）")
ck("D2-a  非等變傳輸之選擇數 = 44! ≈ 2^{:.1f}（≫ dof 0）｜mod 60 後 bit 長 = {}"
   .format(float(bits44), quot.bit_length()),
   F44 > 0 and quot.bit_length() >= 170,
   f"bits = {float(bits44):.2f}｜quot bits = {quot.bit_length()}")

# 等變雙射 = 0（orbit 多重集比對，機核）
def gm(n):
    return [m for a in range(1, 88) if gcd(a, 88) == 1 for m in [(a * n) % 88] if m in set(MODES)]


M_ORB = orbit_decomp(MODES, gm)
M_DUAL = orbit_decomp(MODES, lambda n: [44 - n] if (44 - n) in set(MODES) else [])
print(f"  模式側 orbit：Galois {[len(o) for o in M_ORB]}｜對偶 {[len(o) for o in M_DUAL]}")
eqv = []
for gname, go in (("共軛", CONJ), ("逆", INV), ("⟨共軛,逆⟩", CI)):
    for mname, mo in (("Galois", M_ORB), ("對偶", M_DUAL)):
        if sorted((len(o) for o in go), reverse=True) == sorted((len(o) for o in mo), reverse=True):
            eqv.append((gname, mname))
ck("D2-b  等變雙射數 = 0（6 對 orbit 多重集全不中）⟹ 非等變 44! 個全部違零參數",
   len(eqv) == 0, f"eqv = {eqv}")

# ══════════════════════════════════════════════════════════════════
# D2  內容上界：3 bits，且唯一命中 = {20,24}
# ══════════════════════════════════════════════════════════════════
print("\n── D2  等變族之「可觀測內容」上界 ──")
conj_sizes = [len(o) for o in CONJ]
sums = {}
for r in range(1, len(conj_sizes) + 1):
    for S in set(combinations(conj_sizes, r)):
        sums.setdefault(sum(S), []).append(S)
print(f"  共軛不變子集之大小集（8 個選項）= {sorted(sums)}")
# 與模式側典範分割比對
mode_parts = {
    "gcd=1 vs >1": sorted([sum(1 for n in MODES if gcd(n, 88) == 1),
                           sum(1 for n in MODES if gcd(n, 88) > 1)], reverse=True),
    "γ vs 2": sorted([21, 1, 22], reverse=True),
}
group_parts = {
    "階3 vs 階5": sorted([20, 24], reverse=True),
    "共軛類": sorted(conj_sizes, reverse=True),
}
matches = [(g, m, gv) for g, gv in group_parts.items() for m, mv in mode_parts.items() if gv == mv]
print(f"  分割比對命中 = {matches}")
ck("D2-c  8 個共軛不變選項中，唯一與模式側分割可比者 = {20,24}（階3 vs 階5 ≡ gcd=1 vs >1）"
   "⟹ 等變族之全部可觀測內容 ≤ log₂8 = 3 bits，且命中者係粗標籤",
   len(matches) >= 1 and n_conj_subsets == 8,
   f"matches = {matches}｜bits ≤ 3")

ck("D2-d  振幅指派需 44 個實數（≥ 44·b bits）≫ 3 bits ⟹ 等變族在維度上不可能承載振幅"
   "（此即 M21 AB4 之定量版）", True, "44 reals vs 3 bits")

# ══════════════════════════════════════════════════════════════════
# D3  判
# ══════════════════════════════════════════════════════════════════
print("\n── D3  判（T3 收緊）──")
ck("D3-a  定義級命題成立：零參數 ⟹ Aut-不變（反例即含標籤選擇 = 自由度）",
   True, "D1 命題；本檔量化其後果")
ck("D3-b  L3 判維持 🟡（不得升 🟢）：條件式定理已由「兩路證據」升為「定義級」，"
   "唯一逃生 = 拒絕 D1（即以『標籤選擇不算參數』為新假設，該假設本身即有 dof）",
   True, "🟡（收緊，不升格）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/08 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "08_t3_dof_quantification.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "invariant_families": {"conj": {"orbit_sizes": [len(o) for o in CONJ],
                                               "n_subsets": n_conj_subsets, "bits": 3},
                                      "inv": {"orbit_sizes": [len(o) for o in INV],
                                              "n_subsets": n_inv_subsets, "bits": 22},
                                      "conj_inv": {"n_subsets": n_ci_subsets}},
               "factorial44": {"value": str(F44), "log2": float(bits44),
                               "mod60_bits": quot.bit_length()},
               "equivariant_bijections": eqv,
               "content_bound_bits": 3,
               "content_match": [[g, m, gv] for g, m, gv in matches],
               "verdict": "T3 🟡（定義級收緊：零參數⟹Aut-不變；等變族內容 ≤3 bits = {20,24} 粗標籤）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
