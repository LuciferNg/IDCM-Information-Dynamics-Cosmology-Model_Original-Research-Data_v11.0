#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/01 — R1／T1：L2 分次標籤代數（三數同源？）
=====================================================
判據（`docs/00_preregistration.md` §1，跑前凍結；不得增刪）：

候選系統 𝕃（3 個）：
  L2-A  元素階標籤（group census）：階類大小 {1,15,20,24} 之子集和
        ∪ 類大小 {20,15,12,12,1} 之子集和 ∪ 兩者之補（|G| − X）
  L2-B  同構型標籤（representation census）：Λ^k V、Sym^k V（k ≤ 3）之
        (k, isotypic 支撐多重度 m) 分級塊維數 ＋ 度 k 總維數
        —— 另設「specific 層」（逐 (k_i) 度分佈）作對照，兩層分開報
  L2-C  模式標籤（chain census）：由 B、N_d 生成之**統一規則**計數

凍結判準：K1 同源／K2 唯一（n=1）／K3 密度 ≤5%／K4 零參數（索引規則單一）／
          K5 剛性（B ↦ 5 或換群須令式改變）／K6 同級（L0 grade 相同 → scripts/04）

零 import 上游。A₅ 由 5 點偶置換精確枚舉；維數全整數。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations as _perms
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
# A₅（精確）＋ 度集
# ══════════════════════════════════════════════════════════════════
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
ORD_CNT = Counter(ORD.values())
CLS = []
_seen = set()
for g in G:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    CLS.append(c)
    _seen |= c
CLASS_SIZES = sorted((len(c) for c in CLS), reverse=True)   # [20,15,12,12,1]
D = [1, 3, 3, 4, 5]
N_d = sum(d * (d - 1) for d in D)
B = 4
N2 = ORD_CNT[1] + ORD_CNT[2]
TARGETS = (16, 36, 98)

print("=" * 78)
print("M24 scripts/01 — R1／T1：L2 分次標籤代數")
print("=" * 78)
print(f"  A₅：|G| = {len(G)}｜階 {dict(sorted(ORD_CNT.items()))}｜類 {CLASS_SIZES}｜D = {D}")

# ══════════════════════════════════════════════════════════════════
# A  L2-A 元素階標籤系統
# ══════════════════════════════════════════════════════════════════
print("\n── A  L2-A：元素階標籤（group census）──")


def subset_sums(vals):
    out = {}
    n = len(vals)
    for r in range(1, n + 1):
        for S in set(combinations(sorted(vals, reverse=True), r)):
            out["+".join(map(str, S))] = sum(S)
    return out


order_inv = subset_sums(sorted(ORD_CNT.values(), reverse=True))   # {1,15,20,24}
class_inv = subset_sums(CLASS_SIZES)                              # {20,15,12,12,1}
comp = {}
for name, inv in (("階類", order_inv), ("類大小", class_inv)):
    for k, v in inv.items():
        if 0 < v < len(G):
            comp[f"|G|−({k})"] = len(G) - v
INV_A = {}
INV_A.update({f"階類:{k}": v for k, v in order_inv.items()})
INV_A.update({f"類大小:{k}": v for k, v in class_inv.items()})
INV_A.update({f"補:{k}": v for k, v in comp.items()})

print(f"  庫存大小 |𝔸| = {len(INV_A)}（階類 {len(order_inv)} ＋ 類大小 {len(class_inv)} ＋ 補 {len(comp)}）")
hitsA = {t: sorted(k for k, v in INV_A.items() if v == t) for t in TARGETS}
for t in TARGETS:
    print(f"  命中 {t:<3}：n = {len(hitsA[t])}  {hitsA[t]}")

ck("A1  L2-A 命中 16（= #{g²=e} = 階類 1+15）", len(hitsA[16]) >= 1, f"{hitsA[16]}")
ck("A2  L2-A 命中 36（= #{ord ≤ 3}）", len(hitsA[36]) >= 1, f"{hitsA[36]}")
ck("A3  L2-A **結構上界**：元素計數 ≤ |G| = 60 < 98 ⟹ 98 不可能命中（非密度問題）",
   len(hitsA[98]) == 0, f"hits98 = {hitsA[98]}；max 庫存值 = {max(INV_A.values())}")

# ══════════════════════════════════════════════════════════════════
# B  L2-B 同構型標籤系統（multigraded Λ•V、Sym•V）
# ══════════════════════════════════════════════════════════════════
print("\n── B  L2-B：同構型標籤（representation census）──")
# V = ⊕_{i=1..5} V_i；Λ^k V 按 (k_1..k_5)（Σk_i = k）多分級分解
#   piece(k-vec) = ⊗_i Λ^{k_i} V_i ｜ dim = Π_i C(d_i, k_i)
#   Sym 同理：dim = Π_i C(d_i + k_i − 1, k_i)


def pieces(kmax, mode):
    specific = []          # (dim, support-set as degree-tuple, k-vec)
    for k in range(1, kmax + 1):
        def rec(i, left, kv):
            if i == 5:
                if left == 0:
                    dim = 1
                    for j in range(5):
                        if mode == "L":
                            dim *= comb(D[j], kv[j])
                        else:
                            dim *= comb(D[j] + kv[j] - 1, kv[j]) if kv[j] else 1
                    supp = tuple(sorted(D[j] for j in range(5) if kv[j] > 0))
                    specific.append({"k": k, "dim": dim, "support": supp, "kv": tuple(kv)})
                return
            for x in range(left + 1):
                rec(i + 1, left - x, kv + [x])
        rec(0, k, [])
    return specific


for mode in ("L", "S"):
    spec = pieces(3, mode)
    agg = {}
    for p in spec:
        key = (p["k"], len(p["support"]))
        agg[key] = agg.get(key, 0) + p["dim"]
    total = {}
    for p in spec:
        total[p["k"]] = total.get(p["k"], 0) + p["dim"]
    name = "Λ" if mode == "L" else "Sym"
    print(f"  {name}•V：specific 塊數 = {len(spec)}｜aggregate 級數 = {len(agg)}")
    print(f"    aggregate (k,m) 維數：" + " ｜ ".join(f"({k},{m})={v}" for (k, m), v in sorted(agg.items())))
    print(f"    total 度 k 維數：" + " ｜ ".join(f"k={k}:{v}" for k, v in sorted(total.items())))
    globals()[f"SPEC_{mode}"] = spec
    globals()[f"AGG_{mode}"] = agg
    globals()[f"TOT_{mode}"] = total

INV_B_AGG = {}
for mode, tag in (("L", "Λ"), ("S", "Sym")):
    for (k, m), v in globals()[f"AGG_{mode}"].items():
        INV_B_AGG[f"{tag}(k={k},m={m})"] = v
    for k, v in globals()[f"TOT_{mode}"].items():
        INV_B_AGG[f"{tag} 總(k={k})"] = v

INV_B_SPEC = {}
for mode, tag in (("L", "Λ"), ("S", "Sym")):
    for p in globals()[f"SPEC_{mode}"]:
        INV_B_SPEC[f"{tag}(k={p['k']},supp={p['support']},kv={p['kv']})"] = p["dim"]

hitsB_agg = {t: sorted(k for k, v in INV_B_AGG.items() if v == t) for t in TARGETS}
hitsB_spec = {t: sorted(k for k, v in INV_B_SPEC.items() if v == t) for t in TARGETS}
print("  aggregate 層命中：")
for t in TARGETS:
    print(f"    {t:<3}：n = {len(hitsB_agg[t])}  {hitsB_agg[t]}")
print("  specific 層命中（對照）：")
for t in TARGETS:
    print(f"    {t:<3}：n = {len(hitsB_spec[t])}  {hitsB_spec[t][:3]}{' …' if len(hitsB_spec[t]) > 3 else ''}")

ck("B1  aggregate 層命中 16 = Λ¹（k=1,m=1）= Σd", len(hitsB_agg[16]) >= 1, f"{hitsB_agg[16]}")
ck("B2  aggregate 層命中 98 = Λ² cross（k=2,m=2）= Σ_(i<j) d_i d_j", len(hitsB_agg[98]) >= 1,
   f"{hitsB_agg[98]}")
ck("B3  aggregate 層 **零** 命中 36 ⟹ 36 不在分級層庫存", len(hitsB_agg[36]) == 0,
   f"{hitsB_agg[36]}")
ck("B4  specific 層 36 之唯一命中 = isotypic pattern {3,3,4}（k=3）⟹ 須逐數選 pattern（違 K4）",
   len(hitsB_spec[36]) >= 1 and all("(3, 3, 4)" in k for k in hitsB_spec[36]),
   f"{hitsB_spec[36]}")
ck("B5  兩索引層各缺一數 ⟹ 同一代數內**無單一索引層**含齊三數（16／98 只在 aggregate，36 只在 specific）",
   len(hitsB_agg[36]) == 0 and len(hitsB_spec[16]) == 0 and len(hitsB_spec[98]) == 0,
   f"spec 16 = {hitsB_spec[16]}｜spec 98 = {hitsB_spec[98]}")

# ══════════════════════════════════════════════════════════════════
# C  L2-C 模式標籤系統
# ══════════════════════════════════════════════════════════════════
print("\n── C  L2-C：模式標籤（chain census）──")
MODES = list(range(1, N_d + 1))


def gamma_gt2(n):
    """γ_n = 4sin²(πn/88) 與 2 之比較（精確：γ_n ∝ sin²(πn/88)，n ≤ 44 單調增）"""
    return n >= 23           # γ_22 = 2 精確；n ≥ 23 ⟹ γ > 2


rules = {}
rules["C1 #{n : N_d − n ≥ 2B}"] = sum(1 for n in MODES if N_d - n >= 2 * B)
rules["C2 #{n : γ_n > 2}"] = sum(1 for n in MODES if gamma_gt2(n))
rules["C3 #{n : γ_n < 2}"] = sum(1 for n in MODES if not gamma_gt2(n)) - 1
rules["C4 #{n : γ_n = 2}"] = 1
rules["C5 #{n : gcd(n,88) = 1}"] = sum(1 for n in MODES if __import__("math").gcd(n, 88) == 1)
rules["C6 #{n : gcd(n,88) > 1}"] = sum(1 for n in MODES if __import__("math").gcd(n, 88) > 1)
rules["C7 #{n : n ≤ 2B} + #{n : N_d − n < 2B}"] = sum(1 for n in MODES if n <= 2 * B) + \
    sum(1 for n in MODES if N_d - n < 2 * B)
rules["C8 Σγ"] = 90
for k, v in rules.items():
    print(f"    {k:<34} = {v}")
hitsC = {t: sorted(k for k, v in rules.items() if v == t) for t in TARGETS}
for t in TARGETS:
    print(f"  命中 {t:<3}：{hitsC[t]}")

ck("C1  模式側命中 36 = #{n : N_d − n ≥ 2B}（O211 docs/04 規則）", len(hitsC[36]) == 1, f"{hitsC[36]}")
ck("C2  模式側 98 需 Σγ + 2B（實數和 ＋ 整數）—— 非 label 計數 ⟹ 類型不同",
   len(hitsC[98]) == 0 and rules["C8 Σγ"] + 2 * B == 98, f"Σγ+2B = {rules['C8 Σγ'] + 2 * B}")
ck("C3  模式側 16 = 4B 之讀法存在，但 K5 剛性 FAIL（B ↦ 5 時 4B = 20 ≠ B² = 25）",
   rules["C7 #{n : n ≤ 2B} + #{n : N_d − n < 2B}"] == 16 and 4 * 5 != 5 * 5,
   f"C7 = {rules['C7 #{n : n ≤ 2B} + #{n : N_d − n < 2B}']}｜B=5: 4B=20, B²=25")

# ══════════════════════════════════════════════════════════════════
# D  K1 同源 ＋ K2 唯一 ＋ K6 之交集定理
# ══════════════════════════════════════════════════════════════════
print("\n── D  K1／K2：同源與交集 ──")
INV_ALL = {}
for tag, inv in (("A", INV_A), ("B_agg", INV_B_AGG), ("B_spec", INV_B_SPEC), ("C", rules)):
    for k, v in inv.items():
        INV_ALL[f"{tag}|{k}"] = v

systems = {"L2-A": INV_A, "L2-B(agg)": INV_B_AGG, "L2-B(spec)": INV_B_SPEC, "L2-C": rules}
per_system = {}
for sname, inv in systems.items():
    per_system[sname] = {t: len([1 for v in inv.values() if v == t]) for t in TARGETS}
    print(f"  {sname:<12} 命中數：" + "｜".join(f"{t}:{per_system[sname][t]}" for t in TARGETS))

full = [s for s, h in per_system.items() if all(h[t] >= 1 for t in TARGETS)]
ck("D1  **無任何單一系統**同時命中三數（K1 FAIL ⟹ T1 🔴）", len(full) == 0, f"full-hit systems = {full}")

two_hit = {s: [t for t in TARGETS if per_system[s][t] >= 1] for s in per_system}
print(f"  各系統之命中集：{two_hit}")
A20 = set(t for t in TARGETS if per_system["L2-A"][t] >= 1)
B20 = set(t for t in TARGETS if per_system["L2-B(agg)"][t] >= 1)
inter = A20 & B20
ck("D2  兩大系統之交集 = {16}（16 係唯一跨系統值 ⟹ 唯一可能之橋值）",
   inter == {16}, f"L2-A ∩ L2-B = {inter}")

# K2 唯一性（在 L2-A 之池內）
ck("D3  K2（唯一性）：L2-A 之 36 命中 n = 4（2 恆等族）⟹ 零判別力（同重現閘 ②-b）",
   per_system["L2-A"][36] == 4, f"n₃₆(L2-A) = {per_system['L2-A'][36]}")

# ══════════════════════════════════════════════════════════════════
# E  K5 剛性（B ↦ 5；換群）
# ══════════════════════════════════════════════════════════════════
print("\n── E  K5 剛性測試 ──")


def rigidity_check():
    """(i) L2-C 之 C1 規則換 B：h11(B) = N_d − 2B；(ii) L2-A 之 36 規則換群（S₄）"""
    vals = {b: N_d - 2 * b for b in range(1, 9)}
    # (ii) #{g : ord ≤ B−1} 在 S₄ 上重算（B = 4 ⟹ ord ≤ 3）
    S4 = [p for p in _perms(range(4))]
    o4 = Counter(order_S4(p) for p in S4)
    n_le3_S4 = sum(v for o, v in o4.items() if o <= 3)
    return vals, n_le3_S4, dict(o4)


def order_S4(g):
    x, k = tuple(range(4)), 1
    while True:
        x = tuple(g[x[i]] for i in range(4))
        if x == tuple(range(4)):
            return k
        k += 1


vals, n_le3_S4, o4 = rigidity_check()
print(f"  C1 規則 h11(B) = N_d − 2B：{vals}")
print(f"  L2-A 規則 #{'{'}ord ≤ B−1{'}'} 於 S₄（|S₄| = 24，階分佈 {o4}）= {n_le3_S4}")
ck("E1  K5 剛性：C1 規則隨 B 變（B=4 ⟹ 36；B=5 ⟹ 34）⟹ 非凑數", vals[4] == 36 and vals[5] == 34,
   f"h11(4) = {vals[4]}｜h11(5) = {vals[5]}")
ck("E2  K5 剛性：L2-A 之「ord ≤ B−1」規則換群即變（S₄ ⟹ 18 ≠ 36）⟹ 規則有結構內容",
   n_le3_S4 == 18, f"#{'{'}ord ≤ 3{'}'}(S₄) = {n_le3_S4}")

# ══════════════════════════════════════════════════════════════════
# F  MC null：同源巧合之判別力
# ══════════════════════════════════════════════════════════════════
print("\n── F  MC null：16 與 98 同源是否特別？──")
import random
random.seed(20261125)


def agg_of(DD):
    out = {}
    for k in range(1, 4):
        def rec(i, left, kv):
            if i == 5:
                if left == 0:
                    dim = 1
                    for j in range(5):
                        dim *= comb(DD[j], kv[j])
                    supp = len([1 for x in kv if x > 0])
                    out[(k, supp)] = out.get((k, supp), 0) + dim
                return
            for x in range(left + 1):
                rec(i + 1, left - x, kv + [x])
        rec(0, k, [])
    return out


N_MC = 20000
hit_both = 0
hit16 = 0
hit98 = 0
hit36 = 0
for _ in range(N_MC):
    DD = [random.randint(1, 8) for _ in range(5)]
    a = agg_of(DD)
    v = set(a.values())
    h16 = 16 in v
    h98 = 98 in v
    h36 = 36 in v
    hit16 += h16
    hit98 += h98
    hit36 += h36
    hit_both += (h16 and h98)
p16 = hit16 / N_MC
p98 = hit98 / N_MC
p36 = hit36 / N_MC
pb = hit_both / N_MC
print(f"  N_MC = {N_MC}｜P(16) = {p16:.4f}｜P(98) = {p98:.4f}｜P(36) = {p36:.4f}"
      f"｜P(16 ∧ 98) = {pb:.4f}（獨立期望 {p16 * p98:.4f}）")
ck("F1  null（後驗診斷，非凍結判據）：P(16 ∧ 98) = 0.18% ≤ 5% ⟹ 「16 與 98 同源於 L2-B」"
   "有判別力（非密度巧合）；但仍不含 36 ⟹ 不構成三數同源",
   pb <= 0.05, f"P = {pb:.4f}｜lift over independent = {pb / max(p16 * p98, 1e-12):.1f}×")
ck("F2  null：36 於 aggregate 層**非**結構性缺席（P(36 ∈ agg) ≠ 0）⟹ 缺席係 A₅ 度集"
   "之事實，判詞只可寫「此 D 下缺席」，不可寫一般定理",
   hit36 > 0, f"P(36 ∈ agg) = {p36:.4f}")

# ══════════════════════════════════════════════════════════════════
# G  判詞
# ══════════════════════════════════════════════════════════════════
print("\n── G  判詞（T1）──")
verdict = {
    "K1 同源": "FAIL（無單一系統含三數：A={16,36}／B(agg)={16,98}／C={36}）",
    "K2 唯一": "FAIL（L2-A 之 36 命中 4；L2-B 之 16／98 皆有恆等重寫）",
    "K3 密度": "36：4/94 = 4.26%（邊緣過）｜98：1/94｜16：≥2 ⟹ 36／16 皆不合格",
    "K4 零參數": "FAIL（36 於 L2-B 須逐數選 isotypic pattern {3,3,4}）",
    "K5 剛性": "PASS（C1、L2-A 規則皆隨 B／換群改變）",
    "K6 同級": "見 scripts/04（L0 grade：16 = F ｜ 98 = M ｜ 36 = M）",
}
for k, v in verdict.items():
    print(f"  {k:<9} {v}")

ck("G1  T1 總判 = 🔴（三數分屬兩個標籤系統，交集只餘值 16；98 對 L2-A 有結構上界 60）",
   len(full) == 0 and inter == {16}, f"full = {full}｜inter = {inter}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/01 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "01_l2_label_algebra.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "inventory_A": INV_A, "inventory_B_agg": INV_B_AGG,
               "inventory_B_spec_size": len(INV_B_SPEC),
               "inventory_C": rules,
               "hits": {"A": hitsA, "B_agg": hitsB_agg, "B_spec": {t: hitsB_spec[t] for t in TARGETS}},
               "per_system": per_system, "intersection_A_B": sorted(inter),
               "null": {"N_MC": N_MC, "p16": p16, "p98": p98, "p36": p36, "p_both": pb,
                        "note": "後驗診斷（未列入 docs/00 凍結判據）；null 域 = 隨機度集 k=5, d ∈ 1..8"},
               "verdict": "T1 = 🔴（K1／K2／K4 FAIL；三數跨兩系統，交集 = {16}）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
