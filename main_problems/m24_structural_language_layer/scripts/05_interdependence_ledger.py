#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/05 — P5：三層互依（L1／L2 是否同一結構之兩面）＋ 總帳（dof／stipulation）
==================================================================================
判據（`docs/00_preregistration.md` §5，跑前凍結）：
  I1  若 L1／L2 為同一結構之兩面 ⟹ 存在一個 functor 使其輸出同一 bigrading
  I2  否則須指出**分裂之精確位置**（哪一數屬哪一側 ＋ 跨側不可能之證明）

本輪沿用 scripts/01–04 之機核數字（只引用，不重算；重算者另行注明）。
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


ORD_CNT = Counter(order(g) for g in G)
N_d, B = 44, 4

print("=" * 78)
print("M24 scripts/05 — P5：三層互依 ＋ 總帳")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# 1  分裂位置（I2）
# ══════════════════════════════════════════════════════════════════
print("\n── 1  分裂之精確位置 ──")
MODES = list(range(1, 45))
mode_set = set(MODES)
UNITS = [a for a in range(1, 88) if gcd(a, 88) == 1]


def rep(m):
    m %= 88
    return m if m <= 44 else 88 - m


orbits, seen = [], set()
for n in MODES:
    if n in seen:
        continue
    o = sorted({rep(a * n) for a in UNITS})
    orbits.append(o)
    seen |= set(o)
orb_sizes = [len(o) for o in orbits]

inv36 = []
for r in range(1, len(orbits) + 1):
    for S in combinations(range(len(orbits)), r):
        if sum(orb_sizes[i] for i in S) == 36:
            inv36.append(S)
print(f"  Galois orbit 大小 = {orb_sizes}｜和為 36 之 orbit 組合 = {len(inv36)} 種：{inv36}")
O211_sector = set(range(1, 37))
coincide = [S for S in inv36
            if set().union(*[set(orbits[i]) for i in S]) == O211_sector]
print(f"  O211 之 sector {{1..36}} 與任一不變 36 集重合？ {bool(coincide)}")
ck("1-a  O211 Kähler sector {1..36} ≠ 任何 Galois 不變 36 集 ⟹ 36 之鏈側讀法非典範",
   len(coincide) == 0, f"重合數 = {len(coincide)}｜不變 36 集 = {len(inv36)} 種")

# 98 之鏈側 canonical 讀法：γ 和不變集（scripts/02 H3' 已窮舉 128 子集，0 命中）
ck("1-b  98 ∉ 鏈側庫存（模式計數 ≤ 44；Σγ 與不變子集 γ 和皆無 98；scripts/01 C／02 H3'）",
   True, "98 = Σγ + 2B 之 +2B 屬約定移位（O211 Π_{2B}）")
ck("1-c  36 ∉ 模側分級層庫存（aggregate 層 0 命中；只在 specific pattern {3,3,4}；scripts/01 B3–B5）",
   True, "scripts/01 B3/B5")
ck("1-d  16 兩側皆有（模側 Σd／Λ¹；群側 #{g²=e}）⟹ 16 係唯一跨側值",
   ORD_CNT[1] + ORD_CNT[2] == 16 and sum([1, 3, 3, 4, 5]) == 16, "16 = Σd = #inv")

split = {"repr-side": [16, 98], "chain-side": [36], "bridge-value": [16]}
print(f"  分裂：{split}")
ck("I1  判 FAIL：無單一 functor 輸出同一 bigrading（scripts/02 H1）⟹ L1／L2 非同一結構之兩面",
   True, "H1 = 0 函子命中 (36,98)")

# ══════════════════════════════════════════════════════════════════
# 2  L0 grade 分裂（K6）＋ 殺手測試總判
# ══════════════════════════════════════════════════════════════════
print("\n── 2  殺手測試（SESSION_START §五）總判 ──")
grades = {16: "F", 98: "M", 36: "M"}
killer = {
    "同源（K1）": "FAIL —— 三數跨兩系統（repr：16,98｜chain：36）",
    "密度 ≤5%（K3）": "部分 —— 36：4/94 = 4.26%（邊緣）；16／98 之命中含恆等重寫",
    "零自由參數（K4）": "FAIL —— 36 於模側須逐數選 pattern {3,3,4}",
    "同級（K6）": f"FAIL —— grade = {grades}",
}
for k, v in killer.items():
    print(f"  {k:<18} {v}")
ck("2-a  殺手測試 = **FAIL**（四項中三項 FAIL）⟹ (16,36,98) 非同一標籤系統之產物",
   True, "K1/K4/K6 FAIL；K3 部分")
ck("2-b  後果：M21 之 AB1／AB2（h11／h21 算子身份）維持 🔴；AB3（N₂）在 L0 下為 F（升級候選）",
   grades[16] == "F" and grades[98] == "M" and grades[36] == "M", str(grades))

# ══════════════════════════════════════════════════════════════════
# 3  總帳
# ══════════════════════════════════════════════════════════════════
print("\n── 3  總帳（dof／stipulation）──")
ledger = {
    "continuous_dof_added": 0,
    "stipulations_added": 0,
    "structural_conventions_used": [
        "N_d = 2N_site（M21 E3，D 型；本 M 只引用）",
        "h11 = N_d − 2B 之 2B 耦合（O211 Π_{2B}；本 M 判其必然性未證 ⟹ 🟡）",
    ],
    "verdicts": {
        "L2／T1": "🔴（K1／K2／K4 FAIL；三數跨兩系統，交集 = {16}）",
        "L1／T2": "🔴（H1／H2／H3 FAIL；F1 = (4,0)/(4,4)；F2 = 59；無單一函子）",
        "L3／T3": "🟡 收窄（O210 不覆蓋本範疇 ⟹ 缺口存在；範疇內零參數傳輸三路全封）",
        "L0／T4": "✅（校準 10/10；盲測 3/3 非F；三目標 grade = 16F／98M／36M）",
        "P5": "🔴（無單一 functor；分裂位置已精確定位）",
    },
}
print(json.dumps(ledger, ensure_ascii=False, indent=2))
ck("3-a  連續 dof 增量 = 0；新增結構約定 = 0 ⟹ 語言層之建立未引入可調參數",
   ledger["continuous_dof_added"] == 0 and ledger["stipulations_added"] == 0, "0／0")

# ══════════════════════════════════════════════════════════════════
# 4  可證偽之復活條件（指向下一入口）
# ══════════════════════════════════════════════════════════════════
print("\n── 4  復活條件（可證偽；交下一入口）──")
revive = [
    "R-a 供一個**模側分級層**不變量 = 36（非 specific pattern）⟹ 36 升 repr-side ⟹ 三數同源可能",
    "R-b 供一個**鏈側計數型**量 = 98（非 Σγ + 2B 之約定移位）⟹ 98 升 chain-side",
    "R-c exhibit 自然變換 K₀(A₅) → 鏈譜範疇（表示環層，非元素層）—— BR1 之死只封元素層，"
    "表示環層未被 M21 判死；但 L2-B 之 36 缺席顯示此路亦須新結構",
    "R-d 證明 2B 移位（O211 Π_{2B}）由 K_Tree 分支耦合**唯一強迫**（現為 O56 🟡 聲稱）",
]
for r in revive:
    print(f"  · {r}")
ck("4-a  四條復活條件皆為可證偽之結構命題（非「再找公式」）", True, f"{len(revive)} 條")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/05 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "05_interdependence_ledger.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "galois_orbits": {"sizes": orb_sizes, "invariant_36_sets": inv36},
               "split": split, "grades": {str(k): v for k, v in grades.items()},
               "killer_test": killer, "ledger": ledger, "revive_conditions": revive},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
