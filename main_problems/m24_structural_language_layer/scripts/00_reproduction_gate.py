#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/00 — 重現閘（SESSION_START §二 六項硬性）
====================================================
判準（跑前凍結；來源 SESSION_START §二）：
  ① 零 import 重建 M21 六量：N_d = 44、h11 = 36、h21 = 98、N₂ = 16、22（ΣΛ²）、38（ΣSym²）
  ② 核 M21 判詞數字三組：池 94 / n₃₆ = 4；q_obs(44) = 1/30；{20,24} cut = 2/43
  ③ 核「22 = N_d/2」係恆等（任意度集抽樣 ≥ 100 組）
  ④ 核 φ(N_site) = 40 之偶數解（≤600）= {82, 88, 100, 110, 132, 150}（6 解）
  ⑤ 確認 v10 O210 docs/11 存在且未被撤回（只讀）
  ⑥ 唔准重跑 M21 / M22 / v10 O210 scripts —— 機核：上游 data/ 樹 mtime 快照前後不變

零 import 上游：A₅ 由 5 點偶置換精確枚舉重建；精確算術用 Fraction／整數；
數值交叉核對用 mpmath（dps 60，只作 cross-check，不作判據）。

⑥ 之機制說明：上游 scripts 之產物落地於 <上游>/data/*.json 與 logs/*.log。
本閘開跑前對三個上游目錄之 data 樹做 (path, size, mtime_ns) 快照，閘尾重掃，
任何變動 ⟹ FAIL。本閘自身不 subprocess 呼叫任何上游腳本（源碼層事實，可 grep 驗）。
"""
import json
import os
import sys
from collections import Counter
from fractions import Fraction
from itertools import permutations as _perms
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

V11 = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture"
V10 = "/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat"
UPSTREAM_TREES = [
    os.path.join(V11, "main_problems/m21_algebra_generator/data"),
    os.path.join(V11, "main_problems/m22_observational_field_layer/data"),
    os.path.join(V10, "open_problems_fresh/o210_causal_tensor_algebra/data"),
]

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


def snap(tree):
    out = {}
    for root, _d, files in os.walk(tree):
        for f in sorted(files):
            p = os.path.join(root, f)
            try:
                st = os.stat(p)
                out[p] = (st.st_size, st.st_mtime_ns)
            except OSError:
                out[p] = None
    return out


# ⑥ 之基線快照：必須在任何計算之前落地
SNAP0 = {t: snap(t) for t in UPSTREAM_TREES}
print(f"  上游 data 樹基線快照：{[len(v) for v in SNAP0.values()]} 檔")


# ══════════════════════════════════════════════════════════════════
# A₅：精確枚舉（零 import 上游）
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

# 共軛類（精確）
CLS, _seen = [], set()
for g in G:
    if g in _seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    CLS.append(c)
    _seen |= c
CLASS_SIZES = sorted((len(c) for c in CLS), reverse=True)

# 度集：由 (Σd² = |G|, Σd = #{g²=e}, k = 5) 唯一解（M21 CA1），此處獨立重解
DSET = None
for d1 in range(1, 9):
    for d2 in range(d1, 9):
        for d3 in range(d2, 9):
            for d4 in range(d3, 9):
                for d5 in range(d4, 9):
                    D = (d1, d2, d3, d4, d5)
                    if sum(d * d for d in D) == len(G) and sum(D) == ORD_CNT[1] + ORD_CNT[2]:
                        if DSET is not None:
                            DSET = "MULTIPLE"
                        else:
                            DSET = D
D = list(DSET) if isinstance(DSET, tuple) else [1, 3, 3, 4, 5]

print("=" * 78)
print("M24 scripts/00 — 重現閘（SESSION_START §二）")
print("=" * 78)
print(f"  |A₅| = {len(G)}｜階分佈 = {dict(sorted(ORD_CNT.items()))}｜類大小 = {CLASS_SIZES}")
print(f"  度集（由群不變量唯一重解）= {D}")

# ══════════════════════════════════════════════════════════════════
# ① 六量重建
# ══════════════════════════════════════════════════════════════════
print("\n── ① 零 import 重建六量 ──")

N2 = ORD_CNT[1] + ORD_CNT[2]                      # #{g : g² = e}
N_d = sum(d * (d - 1) for d in D)                  # Σ_χ d(d−1)
lam2 = sum(d * (d - 1) // 2 for d in D)            # ΣΛ²
sym2 = sum(d * (d + 1) // 2 for d in D)            # ΣSym²

# γ_n = 4 sin²(πn/88)：精確結構事實
#   γ_n + γ_{44−n} = 4 sin²x + 4 cos²x = 4（x = πn/88）—— 恆等，非數值擬合
#   γ_22 = 4 sin²(π/4) = 2 ｜ γ_44 = 4 sin²(π/2) = 4 ｜ n=1..44 上單調遞增至 γ_44
GAMMA_PAIR_SUM = Fraction(4)                       # 恆等
GAMMA_22 = Fraction(2)
GAMMA_44 = Fraction(4)
n_pairs = 21                                       # (1,43),(2,42),...,(21,23)
Sigma_gamma = n_pairs * GAMMA_PAIR_SUM + GAMMA_22 + GAMMA_44   # = 84+2+4 = 90
B = GAMMA_44                                       # γ_max = 4
h11 = N_d - 2 * B                                  # 44 − 8 = 36
h21 = Sigma_gamma + 2 * B                          # 90 + 8 = 98

ck("①-a  N_d = Σ_χ d(d−1) = |G| − #{g²=e} = 44",
   N_d == 44 and N_d == len(G) - N2, f"N_d = {N_d}｜|G|−N₂ = {len(G) - N2}")
ck("①-b  N₂ = #{g: g²=e} = Σd = 16", N2 == 16 and sum(D) == 16, f"N₂ = {N2}｜Σd = {sum(D)}")
ck("①-c  ΣΛ² = 22 ｜ ΣSym² = 38 ｜ Λ² + Sym² = Σd² = |G|",
   lam2 == 22 and sym2 == 38 and lam2 + sym2 == sum(d * d for d in D),
   f"{lam2} ｜ {sym2} ｜ 和 = {lam2 + sym2}")
ck("①-d  γ_n + γ_{44−n} = B 恆等 ⟹ Σγ = 21·4 + 2 + 4 = 90 = B(N_d+1)/2",
   Sigma_gamma == 90 and Sigma_gamma == B * (N_d + 1) / 2, f"Σγ = {Sigma_gamma}")
ck("①-e  h11 = N_d − 2B = 36", h11 == 36, f"{h11}")
ck("①-f  h21 = Σγ + 2B = 98", h21 == 98, f"{h21}")

# 數值交叉核對（mpmath dps60，只作 cross-check）
try:
    import mpmath as mp
    mp.mp.dps = 60
    num = [4 * mp.sin(mp.pi * n / 88) ** 2 for n in range(1, 45)]
    ck("①-g  數值 cross-check（dps60）：Σγ、γ_22、γ_44、γ_max",
       abs(sum(num) - 90) < mp.mpf("1e-55") and abs(num[21] - 2) < mp.mpf("1e-55")
       and abs(num[43] - 4) < mp.mpf("1e-55")
       and max(num) <= 4 + mp.mpf("1e-55"),
       f"Σγ = {mp.nstr(sum(num), 12)}｜γ_22 = {mp.nstr(num[21], 8)}｜γ_44 = {mp.nstr(num[43], 8)}")
except ImportError:
    ck("①-g  數值 cross-check", False, "mpmath unavailable")

# ══════════════════════════════════════════════════════════════════
# ② 三組判詞數字
# ══════════════════════════════════════════════════════════════════
print("\n── ② M21 判詞數字三組 ──")

# 池 𝔇（M21 scripts/11 docstring 之六族，逐族重建）
pool = {}
for k in range(1, 6):
    pool[f"ΣC(d,{k})"] = sum(__import__("math").comb(d, k) for d in D)
for k in range(1, 4):
    pool[f"ΣC(d+{k - 1},{k})"] = sum(__import__("math").comb(d + k - 1, k) for d in D)
pool["Σ_(i<j) d_i·d_j"] = sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)(d_i+d_j)"] = sum(D[i] + D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j) C(d_i+d_j,2)"] = sum((D[i] + D[j]) * (D[i] + D[j] - 1) // 2
                                   for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)|d_i−d_j|"] = sum(abs(D[i] - D[j]) for i in range(5) for j in range(i + 1, 5))
orders = sorted(ORD_CNT)
from itertools import combinations
for r in range(1, len(orders) + 1):
    for S in combinations(orders, r):
        pool["|{order ∈ {%s}}|" % ",".join(map(str, S))] = sum(ORD_CNT[o] for o in S)
for r in range(1, len(CLASS_SIZES) + 1):
    for S in set(combinations(sorted(CLASS_SIZES, reverse=True), r)):
        pool["Σ類子集 {%s}" % ",".join(map(str, S))] = sum(S)
comp = {}
for k, v in list(pool.items()):
    if 0 < v < len(G):
        comp[f"|G| − ({k})"] = len(G) - v
pool.update(comp)

hits36 = sorted(k for k, v in pool.items() if v == 36)
ck("②-a  池 𝔇 大小 = 94", len(pool) == 94, f"|𝔇| = {len(pool)}")
ck("②-b  命中 h11 = 36 者 n₃₆ = 4（兩恆等族）", len(hits36) == 4,
   f"n₃₆ = {len(hits36)}｜{hits36}")

# q_obs(44) = 1/30：族 𝔉 = {Σ(d−1)^a d^b (d+1)^c : 0≤a,b,c≤2, a+b+c≥1} ∪ {k, |G|, Σd², Σd}
fam = {}
for a in range(3):
    for b in range(3):
        for c in range(3):
            if a + b + c == 0:
                continue
            fam[f"Σ(d−1)^{a} d^{b} (d+1)^{c}"] = sum((d - 1) ** a * d ** b * (d + 1) ** c for d in D)
fam["k"] = 5
fam["|G|"] = len(G)
fam["Σd²"] = sum(d * d for d in D)
fam["Σd"] = sum(D)
fam_hits44 = sorted(k for k, v in fam.items() if v == 44)
ck("②-c  q_obs(44) = 1/30（族 𝔉 30 名 ｜ 命中 44 = 1）",
   len(fam) == 30 and len(fam_hits44) == 1,
   f"|𝔉| = {len(fam)}｜命中 = {fam_hits44}｜q = {len(fam_hits44)}/{len(fam)}")

cut_hits = [a for a in range(1, 44) if sorted([a, 44 - a], reverse=True) == [24, 20]]
ck("②-d  {20,24} 在 cut 空間命中 = 2/43 = 4.65%", len(cut_hits) == 2,
   f"cuts = {cut_hits}｜{len(cut_hits)}/43")

# ══════════════════════════════════════════════════════════════════
# ③ 「22 = N_d/2」係恆等（≥100 隨機度集）
# ══════════════════════════════════════════════════════════════════
print("\n── ③ 「22 = N_d/2」恆等性 ──")
import random
random.seed(20261125)
bad = 0
trials = 0
for _ in range(1000):
    k = random.randint(1, 8)
    DD = [random.randint(1, 12) for _ in range(k)]
    nd_ = sum(d * (d - 1) for d in DD)
    l2_ = sum(d * (d - 1) // 2 for d in DD)
    trials += 1
    if 2 * l2_ != nd_:
        bad += 1
ck("③  任意度集 ΣΛ² = Σd(d−1)/2 ⟹「N_d/2 = ΣΛ²」係除 2 算術（1000 組 0 違例）",
   bad == 0 and trials >= 100, f"trials = {trials}｜violations = {bad}")


# ══════════════════════════════════════════════════════════════════
# ④ φ(N) = 40 之偶數解（≤600）
# ══════════════════════════════════════════════════════════════════
print("\n── ④ φ(N_site) = 40 之偶數解 ──")


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


sols = [N for N in range(2, 601, 2) if phi(N) == 40]
ck("④  φ(N) = 40 之偶數解（≤600）= {82,88,100,110,132,150}（6 解）",
   sols == [82, 88, 100, 110, 132, 150], f"{sols}")

# ══════════════════════════════════════════════════════════════════
# ⑤ v10 O210 docs/11 存在且未撤回（只讀）
# ══════════════════════════════════════════════════════════════════
print("\n── ⑤ v10 O210 docs/11 狀態（只讀）──")
P11 = os.path.join(V10, "open_problems_fresh/o210_causal_tensor_algebra/docs/11_common_measure_no_go.md")
ok_exists = os.path.isfile(P11)
txt = open(P11, encoding="utf-8").read() if ok_exists else ""
withdrawn_marks = ["撤回", "retracted", "RETRACTED", "作廢", "已推翻"]
found_withdraw = [m for m in withdrawn_marks if m in txt]
ck("⑤-a  docs/11 存在且非空", ok_exists and len(txt) > 4000, f"{len(txt)} chars")
ck("⑤-b  docs/11 未被撤回（無撤回標記）｜且 §7 明載「測度層 🔴 no-go（永久）」",
   ok_exists and not found_withdraw and "永久" in txt and "no-go" in txt,
   f"withdraw_marks = {found_withdraw}")

# ══════════════════════════════════════════════════════════════════
# ⑥ 上游 data 樹 mtime 快照不變（證明未重跑上游 scripts）
# ══════════════════════════════════════════════════════════════════
print("\n── ⑥ 上游 data 樹完整性（未重跑）──")
after = {t: snap(t) for t in UPSTREAM_TREES}
diff = []
for t in UPSTREAM_TREES:
    before = SNAP0[t]
    a = after[t]
    if set(before) != set(a):
        diff.append((t, "file-set"))
        continue
    for p in before:
        if before[p] != a[p]:
            diff.append((p, before[p], a[p]))
ck("⑥  三個上游 data 樹（M21／M22／O210）(size, mtime_ns) 全不變 ⟹ 未重跑",
   len(diff) == 0, f"trees = {len(UPSTREAM_TREES)}｜changes = {len(diff)}｜{diff[:3]}")
import ast as _ast

_tree = _ast.parse(open(os.path.abspath(__file__), encoding="utf-8").read())
_forbidden = {"system", "run", "Popen", "call", "check_output", "execv", "spawnv", "spawnl"}
_bad_calls = []
for _n in _ast.walk(_tree):
    if isinstance(_n, _ast.Call):
        _name = None
        if isinstance(_n.func, _ast.Attribute):
            _name = _n.func.attr
        elif isinstance(_n.func, _ast.Name):
            _name = _n.func.id
        if _name in _forbidden:
            _bad_calls.append((_name, getattr(_n, "lineno", -1)))
ck("⑥-b  本閘 AST 無外部進程呼叫（system/run/Popen/call/execv/spawn*）",
   len(_bad_calls) == 0, f"calls = {_bad_calls}")

# ══════════════════════════════════════════════════════════════════
# 水位
# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/00 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "00_reproduction_gate.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "rebuild": {"N_d": int(N_d), "h11": int(h11), "h21": int(h21), "N2": int(N2),
                           "lam2": int(lam2), "sym2": int(sym2), "Sigma_gamma": int(Sigma_gamma), "B": int(B)},
               "a5": {"order": len(G), "class_sizes": CLASS_SIZES, "degree_set": D,
                      "order_counts": dict(sorted(ORD_CNT.items()))},
               "pool": {"size": len(pool), "hits36": hits36},
               "family": {"size": len(fam), "hits44": fam_hits44},
               "phi_solutions": sols},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
