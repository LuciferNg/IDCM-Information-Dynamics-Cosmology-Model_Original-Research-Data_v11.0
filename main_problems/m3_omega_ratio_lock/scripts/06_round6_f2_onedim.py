#!/usr/bin/env python3
"""M3 scripts/06 — 第六輪：F2 單鏈 1D 性

命題: 鏈 = 相位圓，唔係空間晶格；鏈維數 = 序參量相位環面之秩。
  (i)  O204: W = W₀e^{iθ}（單一複純量，V=−μ²|W|²+λ|W|⁴, U(1)_W）⟹ 相位空間 S¹ ⟹ rank 1
  (ii) 譜指紋: Coxeter 配對 γ_n + γ_{44−n} = B 對【全部】44 模式成立 ⟺ 1D 環；
       2D 環面 C_m×C_n（88 = 2×44, 4×22, 8×11）只對少數模式成立 ⟹ 可直接判別。
  (iii) 群約束: A₅ 之 abelian 子群全部循環（秩 1）⟹ 二十面體唔提供第二獨立週期。
  (iv) 平移群: Z₈₈ 循環（秩 1）⟹ 單一軌道 ⟹ 一條鏈。
"""
import itertools
import numpy as np
from collections import Counter

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N = 4, 88

# ── §1 譜指紋：1D 環 vs 2D 環面 ──────────────────────────────────────
print("§1 譜指紋（Coxeter 配對 γ_n + γ_{44−n} = B 之覆蓋率）")
def lap_cycle(N):
    A = np.zeros((N, N))
    for i in range(N):
        A[i, (i+1) % N] = 1; A[(i+1) % N, i] = 1
    return np.diag(A.sum(1)) - A

def lap_torus(m, n):
    A = np.zeros((m*n, m*n))
    idx = lambda a, b: (a % m)*n + (b % n)
    for a in range(m):
        for b in range(n):
            for (da, db) in ((1,0), (0,1)):
                A[idx(a,b), idx(a+da, b+db)] = 1
                A[idx(a+da,b+db), idx(a,b)] = 1
    return np.diag(A.sum(1)) - A

def pairing_cover(L, target=B, tol=1e-8):
    ev = np.linalg.eigvalsh(L)
    hit = 0
    for x in ev:
        if np.min(np.abs(ev - (target - x))) < tol:
            hit += 1
    return hit, len(ev), len(set(np.round(ev, 8)))

hit1, tot1, dist1 = pairing_cover(lap_cycle(N))
print(f"  1D C₈₈ : 配對命中 {hit1}/{tot1}（{100*hit1/tot1:.0f}%）；相異特徵值 {dist1}")
chk("1D 環：全部 88 模式皆有 4−λ 配對", hit1 == tot1, f"{hit1}/{tot1}")

for (m, n) in ((2, 44), (4, 22), (8, 11)):
    hit, tot, dist = pairing_cover(lap_torus(m, n))
    print(f"  2D C_{m}×C_{n}: 配對命中 {hit}/{tot}（{100*hit/tot:.0f}%）；相異特徵值 {dist}")
    chk(f"2D C_{m}×C_{n} 之配對覆蓋 < 100%", hit < tot, f"{hit}/{tot}")

# 1D：每個相異特徵值之重數 ≤ 2（單一分支，無偏振）
ev1 = np.round(np.linalg.eigvalsh(lap_cycle(N)), 8)
mult = Counter(ev1)
chk("1D：相異特徵值重數 ≤ 2（單分支，無偏振指標）",
    max(mult.values()) <= 2, f"max mult = {max(mult.values())}")
chk("1D：45 個相異特徵值（k=0..44 各異）", len(mult) == 45, f"{len(mult)}")

# ── §2 A₅ 階結構 ────────────────────────────────────────────────────
print("\n§2 A₅（60 偶置換）階結構")
def comp(a, b): return tuple(a[b[i]] for i in range(5))
def order(a):
    e = tuple(range(5)); x = a; k = 1
    while x != e: x = comp(x, a); k += 1
    return k

from itertools import permutations
A5 = [p for p in permutations(range(5))
      if sum(1 for i in range(5) for j in range(i+1, 5) if p[i] > p[j]) % 2 == 0]
chk("|A₅| = 60", len(A5) == 60, f"{len(A5)}")
ords = Counter(order(g) for g in A5)
print(f"  階分佈: {dict(sorted(ords.items()))}")
chk("階分佈 = {1:1, 2:15, 3:20, 5:24}", dict(ords) == {1:1, 2:15, 3:20, 5:24})
chk("無 6 階（⟹ 無 Z₂×Z₃ ≅ Z₆）", ords.get(6, 0) == 0)
chk("無 4 階", ords.get(4, 0) == 0)

# ── §3 abelian 子群：真結構（含否證）────────────────────────────────
print("\n§3 abelian 子群實際結構（含對「無秩-2」假設之否證）")
def closure_pair(g, h):
    S = {tuple(range(5))}
    frontier = [g, h]
    while frontier:
        x = frontier.pop()
        if x in S: continue
        S.add(x)
        for y in (g, h) + tuple(S):
            z = comp(x, y)
            if z not in S: frontier.append(z)
    return frozenset(S)

pairs = 0; abelians = set(); max_ab = 1
for g, h in itertools.combinations(A5, 2):
    if comp(g, h) != comp(h, g): continue
    pairs += 1
    S = closure_pair(g, h)
    abelians.add(S); max_ab = max(max_ab, len(S))
types = sorted({len(S) for S in abelians})
print(f"  交換對 = {pairs}；abelian 子群 {len(abelians)} 個；階集合 {types}；最大階 {max_ab}")
chk("最大 abelian 子群階 = 5（素數）", max_ab == 5, f"{max_ab}")
chk("★ 存在 V₄（Z₂×Z₂）⟹ 「A₅ 無秩-2」假設【否證】", 4 in types,
    f"階 4 子群 {sum(1 for S in abelians if len(S)==4)} 個")
chk("無 6 階元素（⟹ 無 Z₂×Z₃ ≅ Z₆）", ords.get(6, 0) == 0)
chk("無 4 階元素（⟹ 無 Z₄）", ords.get(4, 0) == 0)
chk("A₅ 冇 8 階 / 11 階 abelian 子群（⟹ 冇 Z₈、Z₁₁）",
    8 not in types and 11 not in types, f"{types}")
chk("88 ∤ 60 ⟹ Z₈₈ ⊄ A₅（鏈唔係群軌道）", 88 % 60 != 0)
print("  ⟹ 群路線【否證】：A₅ 唔可能提供 88 位環面；1D 性唔可以靠群子群論證")

# ── §4 平移群 Z₈₈：秩 1（唯一 2 階元素） ───────────────────────────
print("\n§4 平移群 Z₈₈")
orbit = {0}; x = 0
while True:
    x = (x+1) % N; orbit.add(x)
    if x == 0: break
chk("單一生成元傳遞作用（88 位一條軌道）", len(orbit) == N, f"{len(orbit)}")
divs = [d for d in range(1, N+1) if N % d == 0]
chk("子群數 = τ(88) = 8（全循環）", len(divs) == 8, f"{divs}")
inv2 = [a for a in range(N) if (2*a) % N == 0 and a != 0]
chk("唯一 2 階元素 ⟹ 無 V₄ ⟹ Z₈₈ abelian 秩 = 1", len(inv2) == 1, f"{inv2}")
print("  ⚠️ 註：Z₈₈ 自身秩 1；但此為『鏈之平移群』，唔係由 A₅ 導出（§3 已否證群路線）")

# ── §5 與 3D 空間並存 ───────────────────────────────────────────────
print("\n§5 一致性（空間 3D vs 鏈 1D）")
print("  [CITE] 空間維 = N_gen = 3（O197 docs/06 #11）—— 語料既有, 本輪不重算")
print("  [CITE] K_chain = 形狀/相位扇區（O204 final_result 結論）—— 語料既有")
chk("1D 鏈每模式 1 分量 ⟹ 供唔出偏振 2 態（O221 docs/13）",
    max(mult.values()) <= 2, "重數 ≤ 2 來自 ±n 宇稱，非偏振")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")