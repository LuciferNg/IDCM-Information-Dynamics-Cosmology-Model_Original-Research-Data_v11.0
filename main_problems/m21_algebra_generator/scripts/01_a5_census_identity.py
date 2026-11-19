#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/01 — A₅ 普查恆等式：N_d = B·Σ_χ(χ(1)−1)（R1／T1／P1）
=====================================================================
判準全部照 docs/00_preregistration.md §3（跑前寫死）：

  CA1  度集唯一性：5 件、d≥1、Σd²=60、Σd=16 之解集 |S|=1 且 = {1,3,3,4,5}
  CA2  普查恆等式：Σ_χ (d−1)(d−B) = 0 於 B = 4（精確整數）
  CA3  推論：Σd(d−1) = B·Σ(d−1) = 4·11 = 44
  CA4  對合／FS：|{g: g²=e}| = 16；每元素與其逆共軛（FS 指標全 +1 之先決）
  CA5  函數族密度：q_obs(44) ≤ 5%
  CA6  三元唯一性：(33,44,55) 由 S_k = A + k·B 之生成元對 {A,B} 唯一
  BR1  橋：A₅ → 44 點作用之軌道分解存在性（存在 ≠ 唯一）

零 import 上游；A₅ 用精確置換群枚舉（5 點全 60 偶置換）。
POSTHOC 披露：三元對應之觀察見 docs/00 §2 —— CA5/CA6 屬後驗檢定。
"""
import json
import os
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════
# A₅ 精確構造：5 點偶置換（60 個）
# ══════════════════════════════════════════════════════════════════
def parity(p):
    inv = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv += 1
    return inv % 2


from itertools import permutations as _perms

E = tuple(range(5))
G = [p for p in _perms(range(5)) if parity(p) == 0]
assert len(G) == 60
Gset = set(G)


def compose(a, b):
    """(a∘b)(i) = a[b[i]]"""
    return tuple(a[b[i]] for i in range(5))


def inv(a):
    r = [0] * 5
    for i, v in enumerate(a):
        r[v] = i
    return tuple(r)


ORDER = {}
for g in G:
    x = E
    k = 1
    while True:
        x = compose(x, g)
        if x == E:
            ORDER[g] = k
            break
        k += 1

# 共軛類
classes = []
seen = set()
for g in G:
    if g in seen:
        continue
    cl = set()
    for h in G:
        cl.add(compose(compose(h, g), inv(h)))
    classes.append(sorted(cl))
    seen |= cl
classes.sort(key=lambda c: (len(c), c[0]))

n_invol = sum(1 for g in G if g != E and ORDER[g] == 2)
n_sqrt_e = sum(1 for g in G if compose(g, g) == E)

print("=" * 78)
print("M21 scripts/01 — A₅ 普查恆等式")
print("=" * 78)
print(f"  |A₅| = {len(G)}")
print(f"  #共軛類 k = {len(classes)}   類大小 = {[len(c) for c in classes]}")
print(f"  #對合（order 2）= {n_invol}")
print(f"  #{{g: g² = e}} = {n_sqrt_e}   （= 1 + #對合 = {1 + n_invol}）")
print(f"  元素階分佈 = " + str({k: sum(1 for g in G if ORDER[g] == k) for k in sorted(set(ORDER.values()))}))

ck("A5.0  |G| = 60", len(G) == 60, str(len(G)))
ck("A5.1  k = 5 共軛類", len(classes) == 5, str([len(c) for c in classes]))
ck("A5.2  #對合 = 15", n_invol == 15, str(n_invol))
ck("A5.3  #{g² = e} = 16", n_sqrt_e == 16, str(n_sqrt_e))

# 每元素與其逆共軛 ⟹ 所有特徵標實值 ⟹ FS 指標全 +1
non_selfinv = 0
for g in G:
    gi = inv(g)
    if not any(compose(compose(h, g), inv(h)) == gi for h in G):
        non_selfinv += 1
ck("CA4-a  每元素與其逆共軛（違例 0）⟹ 特徵標全實 ⟹ FS = +1", non_selfinv == 0, f"violations = {non_selfinv}")
ck("CA4-b  Σχ(1) = #{g²=e} = 16（FS 定理，先決已證）", n_sqrt_e == 16, "16")

# 註冊度集（M15 §三 / M16）
D = [1, 3, 3, 4, 5]
k_classes = len(classes)
Sd2 = sum(d * d for d in D)
Sd = sum(D)
print(f"\n  註冊度集 d = {D} ｜ Σd² = {Sd2} ｜ Σd = {Sd} ｜ k = {len(D)}")
ck("CA4-c  Burnside 自檢：Σd² = |G| = 60", Sd2 == len(G), f"{Sd2} vs {len(G)}")
ck("CA4-d  #度 = #共軛類 = 5", len(D) == k_classes, f"{len(D)} vs {k_classes}")
ck("CA4-e  Σd = #{g²=e} = 16", Sd == n_sqrt_e, f"{Sd} vs {n_sqrt_e}")

# ══════════════════════════════════════════════════════════════════
# CA1 度集唯一性（5 件、d ≥ 1、Σd² = 60、Σd = 16）
# ══════════════════════════════════════════════════════════════════
print("\n── CA1 度集唯一性 ──")
sols = []
lo, hi = 1, 16


def rec(prefix, remaining, target_sum, target_sq):
    if remaining == 0:
        if target_sum == 0 and target_sq == 0:
            sols.append(tuple(prefix))
        return
    start = prefix[-1] if prefix else lo
    for d in range(start, hi + 1):
        if d * d > target_sq:
            break
        if d > target_sum:
            break
        rec(prefix + [d], remaining - 1, target_sum - d, target_sq - d * d)


rec([], 5, 16, 60)
sols = sorted(set(sols))
print(f"  解集 S = {sols}  ｜ |S| = {len(sols)}")
ck("CA1  度集唯一 ⟹ {1,3,3,4,5}", len(sols) == 1 and sols[0] == (1, 3, 3, 4, 5),
   f"|S| = {len(sols)}, S = {sols}")

# ══════════════════════════════════════════════════════════════════
# CA2 / CA3：普查恆等式與推論
# ══════════════════════════════════════════════════════════════════
print("\n── CA2 / CA3 普查恆等式 ──")
B = 4
sum_dm1 = sum(d - 1 for d in D)                     # Σ(d−1)
sum_dm1_sq = sum((d - 1) ** 2 for d in D)           # Σ(d−1)²
sum_dd_1 = sum(d * (d - 1) for d in D)              # Σd(d−1)
sum_dsq_1 = sum(d * d - 1 for d in D)               # Σ(d²−1)
centred = sum((d - 1) * (d - B) for d in D)         # Σ(d−1)(d−B)
print(f"  Σ(d−1)     = {sum_dm1}")
print(f"  Σ(d−1)²    = {sum_dm1_sq}")
print(f"  Σd(d−1)    = {sum_dd_1}")
print(f"  Σ(d²−1)    = {sum_dsq_1}")
print(f"  Σ(d−1)(d−B) = {centred}   （B = {B}）")
ck("CA2  加權均值恆等 Σ(d−1)(d−B) = 0（B = 4）", centred == 0, str(centred))
ck("CA3-a  Σd(d−1) = 44", sum_dd_1 == 44, str(sum_dd_1))
ck("CA3-b  B·Σ(d−1) = 4×11 = 44", B * sum_dm1 == 44, f"{B}×{sum_dm1} = {B * sum_dm1}")
ck("CA3-c  兩式精確相等 ⟹ N_d = B·Σ(d−1)", sum_dd_1 == B * sum_dm1, f"{sum_dd_1} = {B * sum_dm1}")
ck("CA3-d  Σ(d−1) = 11 = L₅", sum_dm1 == 11, "11")
ck("CA3-e  Σ(d−1)² = 33 = K_tree", sum_dm1_sq == 33, "33")
ck("CA3-f  Σ(d²−1) = 55 = F₁₀", sum_dsq_1 == 55, "55")
ck("CA3-g  三元 AP：(33, 44, 55) 公差 11", (sum_dm1_sq, sum_dd_1, sum_dsq_1) == (33, 44, 55)
   and sum_dd_1 - sum_dm1_sq == sum_dm1 and sum_dsq_1 - sum_dd_1 == sum_dm1,
   f"{sum_dm1_sq}, {sum_dd_1}, {sum_dsq_1}")

# 等價形式（供文件引用）
ck("CA3-h  Σd(d−1) = |G| − Σd = 60 − 16 = 44", sum_dd_1 == len(G) - Sd, f"{len(G) - Sd}")
ck("CA3-i  33 + 44 + 55 = 132 = 3·44（M16 記賬自檢）", 33 + 44 + 55 == 132 == 3 * 44, "132")

# ══════════════════════════════════════════════════════════════════
# CA5 密度閘：𝔉 = {Σ(d−1)^a d^b (d+1)^c} ∪ {k, |G|, Σd², Σd}
# ══════════════════════════════════════════════════════════════════
print("\n── CA5 函數族密度閘 ──")
fam = {}
for a in range(3):
    for b in range(3):
        for c in range(3):
            if a + b + c == 0:
                continue
            name = f"Σ(d-1)^{a} d^{b} (d+1)^{c}"
            fam[name] = sum(((d - 1) ** a) * (d ** b) * ((d + 1) ** c) for d in D)
fam["k"] = k_classes
fam["|G|"] = len(G)
fam["Σd²"] = Sd2
fam["Σd"] = Sd

vals = {}
for name, v in fam.items():
    vals.setdefault(v, []).append(name)

hits44 = vals.get(44, [])
hits33 = vals.get(33, [])
hits55 = vals.get(55, [])
q44 = len(hits44) / len(fam)
print(f"  族大小 |𝔉| = {len(fam)}")
print(f"  命中 44 者（{len(hits44)} 個）：")
for h in hits44:
    print(f"    · {h}")
print(f"  命中 33 者（{len(hits33)}）：{hits33}")
print(f"  命中 55 者（{len(hits55)}）：{hits55}")
print(f"  q_obs(44) = {len(hits44)}/{len(fam)} = {q44:.4f}")
ck("CA5  密度閘 q_obs(44) ≤ 5%", q44 <= 0.05, f"q_obs = {q44:.4f} ({len(hits44)}/{len(fam)})")

# ══════════════════════════════════════════════════════════════════
# CA6 三元唯一性：S_k = A + k·B, k = 0,1,2 → (33, 44, 55)
# ══════════════════════════════════════════════════════════════════
print("\n── CA6 三元 AP 生成元唯一性 ──")
realized = sorted(set(fam.values()))
pairs = []
for A in realized:
    for Bv in realized:
        if Bv == 0:
            continue
        if (A, A + Bv, A + 2 * Bv) == (33, 44, 55):
            pairs.append((A, Bv))
print(f"  已實現值域（{len(realized)} 個）：{realized}")
print(f"  生成元對 {A,B} 解 = {pairs}")
ck("CA6  (33,44,55) 之生成元對在已實現值域內唯一", len(pairs) == 1, str(pairs))
ck("CA6-b  生成元 = {Σ(d−1)², Σ(d−1)} = {33, 11}", pairs == [(33, 11)], str(pairs))

# ══════════════════════════════════════════════════════════════════
# BR1 橋：A₅ 作用於 44 點之軌道分解存在性
# ══════════════════════════════════════════════════════════════════
print("\n── BR1 群↔晶格橋：44 點 A₅-集 ──")


def closure(gens):
    H = {E}
    frontier = [E]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = compose(x, g)
            if y not in H:
                H.add(y)
                frontier.append(y)
            y2 = compose(g, x)
            if y2 not in H:
                H.add(y2)
                frontier.append(y2)
    return frozenset(H)


subgroups = set()
subgroups.add(frozenset({E}))
for g in G:
    subgroups.add(closure([g]))
for g, h in combinations(G, 2):
    subgroups.add(closure([g, h]))
sub_orders = sorted({len(H) for H in subgroups})
orbit_sizes = sorted({60 // o for o in sub_orders})
print(f"  A₅ 子群階（完整枚舉）= {sub_orders}")
print(f"  ⟹ 可能軌道大小（= 指數）= {orbit_sizes}")

# 44 點集之軌道分解：44 = Σ（可用軌道大小）
decomps = []


def decompose(rem, idx, path):
    if rem == 0:
        decomps.append(tuple(path))
        return
    for j in range(idx, len(orbit_sizes)):
        s = orbit_sizes[j]
        if s <= rem:
            decompose(rem - s, j, path + [s])


decompose(44, 0, [])
decomps = sorted(set(decomps), key=lambda t: (len(t), t))
print(f"  44 之軌道分解數 = {len(decomps)}")
print(f"  例：{decomps[:6]} …")
ck("BR1-a  44 點 A₅-集存在（軌道分解非空）", len(decomps) > 0, f"{len(decomps)} 種分解")
ck("BR1-b  存在 ≠ 唯一（分解數 = 1 才算唯一）", len(decomps) > 1, f"分解數 = {len(decomps)}")
# 注意：44 不在 A₅ 之「傳遞」作用大小內（44 ∤ 60）
ck("BR1-c  44 ∤ 60 ⟹ 無傳遞 A₅-作用於 44 點", (60 % 44) != 0, f"60 mod 44 = {60 % 44}")

# ══════════════════════════════════════════════════════════════════
# 綜合
# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/01 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "01_a5_census.json"), "w", encoding="utf-8") as f:
    json.dump({
        "checks": [{"tag": t, "pass": ok, "detail": d} for t, ok, d in CHECKS],
        "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
        "a5": {"order": len(G), "n_classes": len(classes), "class_sizes": [len(c) for c in classes],
               "n_involutions": n_invol, "n_sqrt_identity": n_sqrt_e,
               "degrees": D, "sum_d2": Sd2, "sum_d": Sd},
        "census": {"sum_dm1": sum_dm1, "sum_dm1_sq": sum_dm1_sq, "sum_dd_1": sum_dd_1,
                   "sum_dsq_1": sum_dsq_1, "centred_at_B": centred, "B": B,
                   "N_d_pred": B * sum_dm1},
        "density": {"family_size": len(fam), "hits_44": hits44, "q_obs_44": q44,
                    "hits_33": hits33, "hits_55": hits55},
        "generator_pair": pairs,
        "bridge": {"subgroup_orders": sub_orders, "orbit_sizes": orbit_sizes,
                   "n_decompositions_44": len(decomps), "examples": [list(x) for x in decomps[:8]]},
    }, f, ensure_ascii=False, indent=2)

sys.exit(0 if nfail == 0 else 1)
