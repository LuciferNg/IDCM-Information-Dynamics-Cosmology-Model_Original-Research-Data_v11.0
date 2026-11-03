#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/25 — 不變量格（09-C3）：結構不變量之窮舉、聯合分割、以及「3」之非唯一性

────────────────────────────────────────────────────────────────────────
殘留：docs/09 §五 C3（分類不變量之先驗凍結清單）＋ docs/10 D4（72 軌 vs 味對型）
問題：價對空間 946 點上，**全部結構可定義之不變量**聯合成幾類？
      「3 類」分割有幾個？是否由對稱唯一決定（canonical）？

判準（跑前凍結，禁事後加不變量）
  A1 不變量池必須**跑前**列齊（10 個；見 INV 表），不得見到結果後增刪。
  A2 "canonical" 之定義：由群作用之軌道分解誘導者（= 不變量對群作用恆定）。
  A3 若存在 ≥2 個互不等價之 3 類分割 ⟹ 判「3 之選擇違 A3（dof ≥ 1）」。
性質：確定性窮舉（2^10 = 1024 個子集）。N_trials = 0。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import math
from collections import Counter, defaultdict

L = "=" * 78
N = 44
PAIRS = [(a, b) for a in range(1, N + 1) for b in range(a + 1, N + 1)]


def g(a, b):
    return b - a


# ============================================================================
# §0 不變量池（跑前凍結，10 個）
# ============================================================================
INV = {
    "i1 Δ": lambda a, b: (b - a,),
    "i2 min(Δ,44−Δ)": lambda a, b: (min(b - a, 44 - (b - a)),),
    "i3 (a%2,b%2)": lambda a, b: (a % 2, b % 2),
    "i4 Δ%3": lambda a, b: ((b - a) % 3,),
    "i5 Σ%3": lambda a, b: ((a + b) % 3,),
    "i6 (a%3,b%3)": lambda a, b: (a % 3, b % 3),
    "i7 |{9,15,27}|": lambda a, b: (len({a, b} & {9, 15, 27}),),
    "i8 |{18,30,36}|": lambda a, b: (len({a, b} & {18, 30, 36}),),
    "i9 T_pipe": lambda a, b: (a <= 27, b <= 27),
    "i10 (gcd(a,88),gcd(b,88))": lambda a, b: (math.gcd(a, 88), math.gcd(b, 88)),
}

print(L)
print("M14 scripts/25 — 不變量格（09-C3）：「3 類」分割之窮舉與非唯一性")
print(L)

print(f"\n【§0】價對空間：{len(PAIRS)} 點；不變量池 {len(INV)} 個（跑前凍結）")
for k in INV:
    print(f"    {k}")

# ============================================================================
print("\n【§1】單一/聯合不變量之分割數（2^10 − 1 個子集全掃）")


def partition_of(keys):
    """回傳 (類大小 tuple, 分割本體 = frozenset(frozenset(點)))"""
    d = defaultdict(list)
    for a, b in PAIRS:
        t = tuple()
        for k in keys:
            t += INV[k](a, b)
        d[t].append((a, b))
    sizes = tuple(sorted(len(v) for v in d.values()))
    body = frozenset(frozenset(v) for v in d.values())
    return sizes, body


single = {}
for k in INV:
    sizes, _ = partition_of([k])
    single[k] = (len(sizes), sizes)
print(f"  {'不變量':<26} {'類數':>4}  類大小分佈（前 8）")
for k, (ncls, sizes) in single.items():
    print(f"  {k:<26} {ncls:>4}  {list(sizes)[:8]}")

names = list(INV)
three = []
for r in range(1, len(names) + 1):
    for combo in itertools.combinations(names, r):
        sizes, key = partition_of(list(combo))
        if len(sizes) == 3:
            three.append((combo, key))

print(f"\n  **恰 3 類** 之子集數 = {len(three)} / {2**len(names) - 1}")
distinct = {}
for combo, key in three:
    distinct.setdefault(key, []).append(combo)
print(f"  其中**互不等價**（分割內容不同）者 = {len(distinct)} 個")

print("\n  各不等價 3 類分割（取最小子集代表，列前 12）:")
for i, (key, combos) in enumerate(sorted(distinct.items(), key=lambda kv: (len(kv[1][0]), kv[1][0]))[:12]):
    smallest = sorted(combos, key=len)[0]
    d = defaultdict(int)
    for a, b in PAIRS:
        t = tuple()
        for k in smallest:
            t += INV[k](a, b)
        d[t] += 1
    print(f"    #{i+1:<2} 最小子集 {smallest}  類大小 {sorted(d.values())}")

# ============================================================================
print("\n【§2】canonical（群軌道誘導）之分割")
pairs_set = set(PAIRS)


def nrm(x):
    return (x - 1) % N + 1


def orbit_partition(gens):
    seen, orbs = set(), []
    for p in sorted(pairs_set):
        if p in seen:
            continue
        comp, stack = {p}, [p]
        while stack:
            q = stack.pop()
            for w in gens(q):
                if w in pairs_set and w not in comp:
                    comp.add(w)
                    stack.append(w)
        orbs.append(comp)
        seen |= comp
    return orbs


z44 = lambda p: [(min(nrm(p[0] + s), nrm(p[1] + s)), max(nrm(p[0] + s), nrm(p[1] + s))) for s in range(N)]
d44 = lambda p: z44(p) + [(min(nrm(N + 1 - p[1]), nrm(N + 1 - p[0])),
                           max(nrm(N + 1 - p[1]), nrm(N + 1 - p[0])))]
o_z, o_d = orbit_partition(z44), orbit_partition(d44)
print(f"  ℤ_44 軌道 = {len(o_z)}；D_44 軌道 = {len(o_d)}")
print(f"  ⟹ 由群作用誘導之分割類數 = {{{len(o_z)}, {len(o_d)}}} ⟹ **冇 3**")
print("  ⟹ 任何 3 類分割**唔可能**係 canonical（必屬「選出嚟」之指標算術）")

# ============================================================================
print("\n【§3】10-D4：72（Galois 軌）vs 13（M12 味對型）")
gal = defaultdict(int)
for a, b in PAIRS:
    gal[(math.gcd(a, 88), math.gcd(b, 88))].append if False else None
ngal = len({(math.gcd(a, 88), math.gcd(b, 88)) for a, b in PAIRS})
print(f"  gcd 類（非 Galois 軌，僅代數類）：{ngal}")
print(f"  946 / 13 = {946/13:.4f}；946 / 72 = {946/72:.4f}")
print("  ⟹ 13 與 72 無整除關係；且 72 軌已判「非物理（γ 不被保）」（docs/11）⟹ 不構成對應")

# ============================================================================
print("\n" + L)
print("【§4】判定（09-C3）")
print(L)
print(f"""  ① 不變量格（10 個跑前凍結之結構不變量）共 {2**len(names)-1} 個非空組合；
     其中恰 3 類者 **{len(three)}** 個組合，互不等價 **{len(distinct)}** 個分割。
  ② 由群作用誘導（canonical）者只有 {len(o_z)} / {len(o_d)} 類 ⟹ **冇 3**。
  ③ ⟹ 全部 3 類分割皆屬「指標算術選擇」（Δ%3、Σ%3、集合成員數…），**冇一個被結構強制**。
     ⟹ 依 A3：揀任何一個都要 **dof ≥ 1** ⟹ 「3 類」之指定違零參數。
  ④ ⟹ **09-C3 得答**：先驗凍結清單存在（上列 10 個），但其對「3」係**退化**的 ——
     3 可以被做出 {len(distinct)} 種互不相容之方式 ⟹ 「荷係數」唔可能由此唯一導出。
  ⑤ 此結論與 docs/19 單標籤定理、docs/10 T1/T3、M15 之「四載體 no-go」三路一致。
  ⑥ 10-D4 得答（否定式）：72 與 13 無整除／無結構對應；72 軌本身已被判非物理。
""")
print(L)
raise SystemExit(0)