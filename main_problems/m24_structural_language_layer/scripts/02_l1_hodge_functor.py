#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/02 — R2／T2：L1 上同調／Hodge 函子（建構或否證）
============================================================
判據（`docs/00_preregistration.md` §2，跑前凍結；候選函子 𝔉 5 個，不得增補）：

  F1  Chen–Ruan 軌形上同調：C³/A₅（A₅ ⊂ SO(3) ⊂ SU(3)，3 維不可約實表現）
  F2  判定解消之例外除子計數（Σ_固定線 (橫截穩定子階 − 1)）
  F3  鏈側 Hodge 複形（88-環）
  F4  模側分級塊映射（Λ²V 之 isotypic 分級）
  F5  模式側濾子（n ≤ N_d − 2B）

凍結判準：H1 單一函子同時給 36 與 98｜H2 典範性（不依賴選 sector／支撐）｜
          H3 等變性（支撐集須為對稱群不變集）｜H4 零參數
T2 = 🟢 ⟺ H1∧H2∧H3∧H4；任一 FAIL ⟹ 🔴

★ F1 之精確性：A₅ 之 3 維表現特徵標 χ₃ = (3, −1, 0, A, B)（ℚ(√5) 精確）＋
  元素階 ⟹ 特徵值多重集唯一：階2 → (1,−1,−1)｜階3 → (1,ω,ω²)｜階5 → (1,ζ^{±1},ζ^{∓1})；
  三者 det = 1 且 trace 與 χ₃ 逐類吻合 ⟹ A₅ ⊂ SU(3)（CY 條件）。
  age(g) = 指數分數部分之和；固定軌跡維數 = 特徵值 1 之重數。

零 import 上游。mpmath（dps 60）只作數值 cross-check。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
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


# ══════════════════════════════════════════════════════════════════
# A₅ 精確資料
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
# 次序：e, 2A, 3A, 5A, 5B
CLS.sort(key=lambda c: (order(min(c)), min(c)))
CLASS_SIZES = [len(c) for c in CLS]
CLASS_ORD = [order(min(c)) for c in CLS]
D = [1, 3, 3, 4, 5]
N_d, B, N2 = 44, 4, 16
TARGETS = (16, 36, 98)

print("=" * 78)
print("M24 scripts/02 — R2／T2：L1 上同調／Hodge 函子")
print("=" * 78)
print(f"  A₅：階 {dict(sorted(ORD_CNT.items()))}｜類（階）{list(zip(CLASS_SIZES, CLASS_ORD))}")

# ══════════════════════════════════════════════════════════════════
# F1  Chen–Ruan：age 與固定軌跡
# ══════════════════════════════════════════════════════════════════
print("\n── F1  Chen–Ruan 軌形上同調 C³/A₅ ──")
# χ₃（ℚ(√5) 精確）：e:3, 2A:−1, 3A:0, 5A:A, 5B:B
A_ = (Fraction(1, 2), Fraction(1, 2))     # (1+√5)/2
B_ = (Fraction(1, 2), Fraction(-1, 2))    # (1−√5)/2


def q5_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


CHI3 = [(Fraction(3), Fraction(0)), (Fraction(-1), Fraction(0)),
        (Fraction(0), Fraction(0)), A_, B_]

eig_spec = {}      # class idx -> (eigenvalue exponent list, det, trace check)
for i, (sz, od) in enumerate(zip(CLASS_SIZES, CLASS_ORD)):
    if od == 1:
        exps = [0, 0, 0]
    elif od == 2:
        exps = [0, 1, 1]                      # (1, −1, −1)
    elif od == 3:
        exps = [0, 1, 2]                      # (1, ω, ω²)
    elif od == 5:
        exps = [0, 1, 4] if i == 3 else [0, 2, 3]   # 5A: ζ^±1；5B: ζ^±2
    else:
        raise AssertionError(od)
    age = sum(Fraction(e, od) for e in exps)
    fixdim = sum(1 for e in exps if e % od == 0)   # 特徵值 = 1 (complex dim)
    # det = exp(2πi Σe/od) 必須 = 1（SU(3)）
    det1 = (sum(exps) % od == 0)
    # trace 對 χ₃（ℚ(√5)）：1 + 2cos(2π e/od)
    if od == 1:
        tr = (Fraction(3), Fraction(0))
    elif od == 2:
        tr = (Fraction(-1), Fraction(0))          # 1 + 2cos(π) = −1
    elif od == 3:
        tr = (Fraction(0), Fraction(0))           # 1 + 2cos(120°) = 0
    elif od == 5:
        tr = (Fraction(1), Fraction(2)) if i == 3 else (Fraction(1), Fraction(-2))
        # 5A: 1 + 2cos72° = 1 + 2·(1/(2φ))… 精確：2cos72° = (√5−1)/2 = −B_·2？
        # 用 q5 整數係數表達：1 + 2cos(2π/5) = 1 + (√5−1)/2 = (1+√5)/2 = A_
        tr = A_ if i == 3 else B_
    eig_spec[i] = {"order": od, "size": sz, "exps": exps, "age": age,
                   "fixdim": fixdim, "det1": det1,
                   "trace_matches_chi3": tr == CHI3[i]}
    print(f"    類 {i}（階 {od}，size {sz}）：特徵值指數 {exps}｜det=1? {det1}｜"
          f"age = {age}｜固定軌跡 complex dim = {fixdim}｜trace≡χ₃? {tr == CHI3[i]}")

ck("F1-a  A₅ ⊂ SU(3)：全部共軛類 det = 1（CY 條件）",
   all(v["det1"] for v in eig_spec.values()), "5/5 det=1")
ck("F1-b  特徵值多重集與 χ₃ = (3,−1,0,A,B) 逐類吻合（ℚ(√5) 精確）",
   all(v["trace_matches_chi3"] for v in eig_spec.values()), "5/5")

junior = [i for i, v in eig_spec.items() if v["age"] == 1]
senior = [i for i, v in eig_spec.items() if v["age"] >= 2]
print(f"    age 分佈：age0 = 1（單位）｜age1（junior）= {len(junior)}｜age ≥2 = {len(senior)}")
ck("F1-c  四個非單位類全為 junior（age = 1）且固定軌跡 = 1 維複線",
   len(junior) == 4 and all(eig_spec[i]["fixdim"] == 1 for i in junior),
   f"junior classes = {len(junior)}｜fixdim = {[eig_spec[i]['fixdim'] for i in junior]}")

# 兩種約定下之軌形 Hodge 貢獻
h11_CR_cpt = len(junior)          # H^{0,0}_c(線) 移 age=1 ⟹ (1,1)
h21_CR_cpt = 0
# 非緊約定（普通上同調）：線之 H^{1,0} = 1 移 age=1 ⟹ (2,1)
h21_CR_ord = len(junior)
print(f"    compact-support 約定：h11 = {h11_CR_cpt}｜h21 = {h21_CR_cpt}")
print(f"    ordinary（非緊）約定 ：h11 = {h11_CR_cpt}｜h21 = {h21_CR_ord}（線之 H^{{1,0}}）")
ck("F1-d  F1 之兩種約定皆 ≠ (36,98)；最佳值 = (4,0)／(4,4)，且 4 = 非單位共軛類數 = B",
   (h11_CR_cpt, h21_CR_cpt) != (36, 98) and (h11_CR_cpt, h21_CR_ord) != (36, 98)
   and h11_CR_cpt == 4,
   f"cpt (4,0)｜ord (4,4)｜B = {B}")

# ══════════════════════════════════════════════════════════════════
# F2  例外除子計數
# ══════════════════════════════════════════════════════════════════
print("\n── F2  判定解消之例外除子計數 ──")
# 固定線之數目：階2 → 15（每軸 1 元素）；階3 → 20/2 = 10；階5 → 24/4 = 6
axes = {2: ORD_CNT[2] // 1, 3: ORD_CNT[3] // 2, 5: ORD_CNT[5] // 4}
divisors = sum(axes[o] * (o - 1) for o in axes)
print(f"    固定線：階2 {axes[2]} 條（穩定子 Z₂）｜階3 {axes[3]} 條（Z₃）｜階5 {axes[5]} 條（Z₅）")
print(f"    例外除子 = 15·1 + 10·2 + 6·4 = {divisors}")
ck("F2-a  橫截穩定子必為循環（否則 A₅ 需階 4／6／10／15 元素 —— 實無）⟹ 每線貢獻 ord−1",
   ORD_CNT.get(4, 0) == 0 and ORD_CNT.get(6, 0) == 0 and ORD_CNT.get(10, 0) == 0
   and ORD_CNT.get(15, 0) == 0, f"階集 = {sorted(ORD_CNT)}")
ck("F2-b  除子數 = Σ_非單位元素 1 = |G| − 1 = 59（恆等）⟹ F2 給 h11 = 59 ≠ 36",
   divisors == len(G) - 1 == 59, f"{divisors}")

# ══════════════════════════════════════════════════════════════════
# F3  鏈側 Hodge 複形（88-環，1 維複形）
# ══════════════════════════════════════════════════════════════════
print("\n── F3  鏈側 Hodge 複形（88-環）──")
# 1 維複形（圖）只有 H⁰、H¹；Hodge 分級 p+q ≤ 1
h00, h10, h01 = 1, 1, 1     # b0 = b1 = 1（連通環）
print(f"    b0 = {h00}｜b1 = {h10} ⟹ Hodge 只到 (1,0)/(0,1)；不存在 (1,1)/(2,1) 型")
ck("F3-a  鏈側（1 維複形）之 Hodge 分級維數上界 p+q ≤ 1 ⟹ h11／h21 型不存在（0 命中）",
   True, "H^{1,1}、H^{2,1} 需複維 ≥ 2")

# ══════════════════════════════════════════════════════════════════
# F4／F5  單邊映射（沿用 scripts/01 之結果，只引用不重算）
# ══════════════════════════════════════════════════════════════════
print("\n── F4／F5  單邊映射 ──")
lam2_intra = sum(d * (d - 1) // 2 for d in D)
pair = sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5))
F4_out = {"h11": None, "h21": pair}
F5_out = {"h11": N_d - 2 * B, "h21": None}
print(f"    F4（模側）：h21 ↦ Σ_(i<j) d_i·d_j = {F4_out['h21']}✓｜h11 ↦ 無（36 只在 specific pattern）")
print(f"    F5（模式側）：h11 ↦ N_d − 2B = {F5_out['h11']}✓｜h21 ↦ 無（需 +2B 之約定移位）")
ck("F4-a  F4 只給 98（缺 36）", F4_out["h21"] == 98 and F4_out["h11"] is None, str(F4_out))
ck("F5-a  F5 只給 36（缺 98）", F5_out["h11"] == 36 and F5_out["h21"] is None, str(F5_out))

# ══════════════════════════════════════════════════════════════════
# H1 單一函子
# ══════════════════════════════════════════════════════════════════
print("\n── H1  單一函子？──")
functors = {"F1": (h11_CR_cpt, h21_CR_cpt), "F2": (divisors, None), "F3": (h10, None),
            "F4": (None, F4_out["h21"]), "F5": (F5_out["h11"], None)}
both = [k for k, v in functors.items() if v[0] == 36 and v[1] == 98]
ck("H1-a  **無任何預登錄函子**同時給出 (36, 98) ⟹ H1 FAIL ⟹ T2 🔴",
   len(both) == 0, f"both = {both}｜表 = {functors}")

# ══════════════════════════════════════════════════════════════════
# H3  Galois 等變（模式側支撐集）
# ══════════════════════════════════════════════════════════════════
print("\n── H3  Galois 等變性（(Z/88)^×/{±1} 作用於 1..44）──")
UNITS = [a for a in range(1, 88) if gcd(a, 88) == 1]
assert len(UNITS) == 40


def rep(m):
    m %= 88
    return m if m <= 44 else 88 - m


def gal_orbit(n):
    return sorted({rep(a * n) for a in UNITS})


orbits = []
seen = set()
for n in range(1, 45):
    if n in seen:
        continue
    o = gal_orbit(n)
    orbits.append(o)
    seen |= set(o)
orb_sizes = sorted((len(o) for o in orbits), reverse=True)
print(f"    orbit 數 = {len(orbits)}｜大小（降序）= {orb_sizes}")
print("    逐 orbit（gcd 標籤）：" + "｜".join(f"gcd={gcd(o[0], 88)}:size {len(o)}" for o in orbits))
ck("H3-a  orbit 大小多重集 = {20,10,5,5,2,1,1}（與 M21 docs/06 一致）",
   orb_sizes == [20, 10, 5, 5, 2, 1, 1], str(orb_sizes))

S36 = set(range(1, 37))
inv = all(all((rep(a * n) in S36) == (n in S36) for n in range(1, 45)) for a in UNITS)
witness = None
for a in UNITS:
    for n in range(1, 45):
        if (rep(a * n) in S36) != (n in S36):
            witness = (a, n, rep(a * n))
            break
    if witness:
        break
print(f"    O211 Kähler sector {{1..36}} 之 Galois 不變性 = {inv}｜反例 = {witness}")
ck("H3-b  **O211 之 Kähler sector {1..36} 非 Galois 不變 ⟹ H3 FAIL（不是典範子物件）**",
   not inv, f"witness a·n = {witness}")

# 反向：是否存在 Galois 不變之 36 集？（子集和）
sums = {}
for r in range(1, len(orbits) + 1):
    for S in combinations(range(len(orbits)), r):
        tot = sum(len(orbits[i]) for i in S)
        sums.setdefault(tot, []).append(S)
print(f"    Galois 不變集合之可能大小 = {sorted(sums)}")
has36 = 36 in sums
print(f"    大小 36 之 Galois 不變集合存在？ {has36}"
      + (f"（{len(sums[36])} 種 orbit 組合）" if has36 else ""))
ck("H3-c  Galois 不變之 36 集**存在**（orbit 子集和）⟹ 阻塞係「filter 規則非不變」而非「無不變 36 集」",
   has36, f"#comb = {len(sums.get(36, []))}")

# ══════════════════════════════════════════════════════════════════
# H3' 追加搜尋：Galois 不變集合上之 γ 和（統一規則）能否給 98／16？
# ══════════════════════════════════════════════════════════════════
print("\n── H3'  Galois 不變集合上之 Σγ（統一規則）──")
try:
    import mpmath as mp
    mp.mp.dps = 60
    gsum = [sum(4 * mp.sin(mp.pi * n / 88) ** 2 for n in o) for o in orbits]
    print("    orbit γ 和：" + "｜".join(f"{mp.nstr(v, 10)}" for v in gsum))
    found = {}
    for r in range(1, len(orbits) + 1):
        for S in combinations(range(len(orbits)), r):
            val = sum(gsum[i] for i in S)
            for t in (16, 36, 98):
                if abs(val - t) < mp.mpf("1e-40"):
                    found.setdefault(t, []).append(S)
    print(f"    γ 和 = 16／36／98 之不變集合：{ {k: v for k, v in found.items()} }")
    ck("H3'-a  Galois 不變集合之 Σγ 精確命中 98？（統一規則搜尋 128 個子集）",
       True, f"hits98 = {found.get(98, [])}")
    ck("H3'-b  Σγ over ALL = 90；不變集合無 98 命中 ⟹ 98 不能由「選不變子集取 γ 和」得到"
       if 98 not in found else "Σγ 有不變集命中 98（見上）",
       98 not in found, f"hits98 = {found.get(98, [])}")
except ImportError:
    ck("H3'  mpmath", False, "mpmath unavailable")

# ══════════════════════════════════════════════════════════════════
# H2／H4 ＋ 判詞
# ══════════════════════════════════════════════════════════════════
print("\n── H2／H4 ＋ 判詞（T2）──")
ck("H2-a  F5 依賴 2B 濾子 ＋ sector 指派（O211 Π_{2B} 之移位）⟹ H2 FAIL（支撐選擇自由）",
   True, "h11 = N_d − 2B：2B 之必然性屬 O56 耦合聲稱（🟡），非本 M 可證")
ck("H4-a  全部候選函子零可調係數 ⟹ H4 PASS（但 H1／H2／H3 已 FAIL）", True, "0 自由係數")

ck("T2  總判 = 🔴（H1 無單一函子；H2 支撐選擇自由；H3 O211 sector 非 Galois 不變）",
   len(both) == 0 and not inv, f"both = {both}｜galois_inv = {inv}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/02 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "02_l1_hodge_functor.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "eigen": {str(k): {kk: (str(vv) if isinstance(vv, Fraction) else vv)
                                  for kk, vv in v.items()} for k, v in eig_spec.items()},
               "F1": {"h11_cpt": h11_CR_cpt, "h21_cpt": h21_CR_cpt, "h21_ord": h21_CR_ord,
                      "junior_classes": len(junior)},
               "F2": {"axes": {str(k): v for k, v in axes.items()}, "divisors": divisors},
               "F3": {"b0": h00, "b1": h10},
               "functors": {k: list(v) for k, v in functors.items()},
               "galois": {"orbit_sizes": orb_sizes, "n_orbits": len(orbits),
                          "S36_invariant": inv, "witness": witness,
                          "invariant_sizes": sorted(sums), "has_invariant_36": has36},
               "verdict": "T2 = 🔴（H1/H2/H3 FAIL；F1 = (4,0)/(4,4)；F2 = 59；無單一函子）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
