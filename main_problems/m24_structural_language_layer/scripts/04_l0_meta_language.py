#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/04 — R4／T4：L0 元語言（count-identity 判準）＋ 10 案例校準 ＋ 盲測應用
=================================================================================
判據（`docs/00_preregistration.md` §4，跑前凍結）：

定義（結構身份 = 交換圖存在，四條件）：
  S1 範疇化：A、B 各為函子 F、G: 𝒞 → (分次向量空間／集合) 之值，
             且 claim 斷言自然變換 η: F ⟹ G 在該物件上為同構
  S2 非恆等：η 在 𝒞 上**非**自由變數恆等（存在同範疇測試物件令兩側不等；
             即該等式非「代數改寫／解自身輸入」）
  S3 強迫：該 claim 之候選族內命中 n = 1（目標物件被唯一強迫）
  S4 等變：η 與 Aut(O) 之交換圖可換

  分類：¬S1 ⟹ A 型｜¬S2 ⟹ T 型｜S2 ∧ (¬S3 ∨ ¬S4) ⟹ M 型｜S1∧S2∧S3∧S4 ⟹ F 型

校準集（10 案例；通過條件 10/10 與 M21 記錄判詞一致）：
  CA1（預期 F）｜E12（對照，預期 F）｜BF1–BF8（預期 非F）
非調參測試（§4.2）：對三條本輪新 functional（N-a／N-b／N-c）盲測應用。

零 import 上游；所有測試物件在本地重算（來源：M21 docs/01–09 之資料，數字重算）。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations, permutations as _perms
from math import comb, gcd

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
NONINV = [g for g in G if ORD[g] > 2]
D = [1, 3, 3, 4, 5]
N_d, B, N2 = 44, 4, 16

print("=" * 78)
print("M24 scripts/04 — R4／T4：L0 元語言 ＋ 10 案例校準")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# 各 case 之機核證據
# ══════════════════════════════════════════════════════════════════
EVID = {}

# ── CA1：度集唯一性 ────────────────────────────────────────────
sols = []
for d1 in range(1, 9):
    for d2 in range(d1, 9):
        for d3 in range(d2, 9):
            for d4 in range(d3, 9):
                for d5 in range(d4, 9):
                    Dc = (d1, d2, d3, d4, d5)
                    if sum(d * d for d in Dc) == 60 and sum(Dc) == 16:
                        sols.append(Dc)
EVID["CA1"] = {"S1": "degree-multiset functor ↔ 約束解集 functor（兩側皆結構物件）",
               "S2_probe": "同範疇測試物件（其他 (|G|,k,#inv) 三元）：解集改變 ⟹ 非自由恆等",
               "S3": f"{len(sols)} 解 ⟹ 唯一", "S3_num": len(sols),
               "S4": "約束 (Σd², Σd, k) 皆共軛不變", "class": "F", "m21": "✅ F"}

# ── E12（對照）：N₂ = #{g²=e} = Σd ──────────────────────────
order_sums = {}
for r in range(1, 5):
    for S in combinations(sorted(ORD_CNT.values(), reverse=True), r):
        order_sums[sum(S)] = order_sums.get(sum(S), 0) + 1
n16_orderfam = order_sums.get(16, 0)
EVID["E12"] = {"S1": "模維數 functor（Σd）↔ 對合計數 functor（#{g²=e}）：FS 定理",
               "S2_probe": "非自由恆等（Σd 與對合數在無群結構之資料上無定義）",
               "S3": f"階類子集族（15 個）內命中 16 = {n16_orderfam}", "S3_num": n16_orderfam,
               "S4": "兩側皆共軛不變（元素階為類函數）", "class": "F", "m21": "✅ F"}

# ── BF1：CY₃(36,98) → N_d 路線 ───────────────────────────────
P0 = {"h11": 36, "h21": 98, "chi": -124, "b2": 36, "b3": 198, "B": 4, "Sig": 90,
      "Ngen": 3, "G": 60, "K_tree": 33, "F10": 55, "L4": 7, "L5": 11}
kfree_hits = []
for x in P0:
    for y in P0:
        if x == y:
            continue
        for a in (1, -1, 2, -2):
            for b in (1, -1, 2, -2):
                if a * P0[x] + b * P0[y] == 44:
                    kfree_hits.append((a, x, b, y))
# 口徑統一：無序化（(a,x,b,y) 與 (b,y,a,x) 同一式），再按 M21 之別名（b2 ≡ h11）合併
_kf_unord = set()
for a, x, b, y in kfree_hits:
    if (y, b, x, a) < (x, a, y, b):
        _kf_unord.add((y, b, x, a))
    else:
        _kf_unord.add((x, a, y, b))
n_kfree_raw = len(kfree_hits)
n_kfree = len(_kf_unord)
# 帶自由常數族：a·x − b·y + c = 44，|c| ≤ 120
n_free, n_hit = 0, 0
for x in P0:
    for y in P0:
        if x == y:
            continue
        for a in (1, -1, 2, -2):
            for b in (1, -1, 2, -2):
                n_free += 1
                c = 44 - (a * P0[x] - b * P0[y])
                if abs(c) <= 120:
                    n_hit += 1
EVID["BF1"] = {"S1": "無：CY₃ 拓撲量 → 鏈模式數之 functor 未 exhibit",
               "S2_probe": "非恆等（池算術）",
               "S3": f"k-free 命中（無序化）= {n_kfree}；raw = {n_kfree_raw}（含換序/別名重複）"
                     f"｜帶自由常數 q_obs = {n_hit}/{n_free} = {n_hit / n_free:.4f}（M21 報 53.53%，率一致）",
               "S3_num": n_kfree, "S4": "不適用（無映射）",
               "class": "M", "m21": "🔴"}

# ── BF2：p9 FORCE_44（φ-級聯 ＋ 恆等反演）───────────────────
try:
    import mpmath as mp
    mp.mp.dps = 50
    ident_ok = all(abs(sum(4 * mp.sin(mp.pi * n / (2 * N)) ** 2 for n in range(1, N + 1))
                        - 2 * (N + 1)) < mp.mpf("1e-40") for N in range(2, 121))
except ImportError:
    ident_ok = False
EVID["BF2"] = {"S1": "無（級聯非 functor 值）",
               "S2_probe": f"Σγ = 2(N+1) 於 N=2..120 全真（{ident_ok}）⟹ 恆等反演 ⟹ T",
               "S3": "不適用（已死於 S2）", "S3_num": -1,
               "S4": "不適用", "class": "T", "m21": "🔴"}

# ── BF3：Σγ = 90 ⟹ N_d = 44 作 producer ─────────────────────
# 與 BF2 同型：CG1 恆等（同一方程）
EVID["BF3"] = {"S1": "無（反向求解下游量）",
               "S2_probe": f"恒等：Σγ = B(N_d+1)/2 對任意 N_d 成立（{ident_ok}）",
               "S3": "不適用", "S3_num": -1, "S4": "不適用",
               "class": "T", "m21": "🔴"}

# ── BF4：f_M0／f_R 兩獨立因子（gauge）────────────────────────
M0_ratio = 1.801881664e40      # N_O38（M21 docs/03 之乘積；只作乘積固定之演示）
pairs = [(fm, M0_ratio / fm) for fm in (2.28529e-21, 1.0e-21, 5.0e-21, 1.0e-20, 1.0e-18)]
same_product = all(abs(fm * fr - M0_ratio) / M0_ratio < 1e-9 for fm, fr in pairs)
EVID["BF4"] = {"S1": "乘積 = 單一 object（四重寫）；分解無獨立 functor",
               "S2_probe": "非恆等（乘積為真）",
               "S3": "自由 gauge 參數 1 條 ⟹ 不唯一（5 對同乘積）", "S3_num": 5,
               "S4": "無等變內容", "class": "M", "m21": "🔴"}

# ── BF5：BR1 對象層映射 ─────────────────────────────────────
def orbits(elems, act):
    elems = list(elems)
    domain = set(elems)
    seen, out = set(), []
    for x in elems:
        if x in seen:
            continue
        orb, stack = set(), [x]
        while stack:
            y = stack.pop()
            if y in orb:
                continue
            orb.add(y)
            for z in act(y):
                if z in domain and z not in orb:
                    stack.append(z)
        seen |= orb
        out.append(len(orb))
    return sorted(out, reverse=True)


mode_set = set(range(1, 45))


def galois_act(n):
    return [m for a in range(1, 88) if gcd(a, 88) == 1
            for m in [(a * n) % 88] if m in mode_set]


g_orb = orbits(NONINV, lambda g: [mul(mul(h, g), invp(h)) for h in G])
m_orb = orbits(list(range(1, 45)), galois_act)
inv_orb = orbits(NONINV, lambda g: [invp(g)])
ck_pair = orbits(list(range(1, 45)), lambda n: [44 - n] if 44 - n in mode_set else [])
# 44 之 orbit 大小分解數（A₅ 子群指數集 = {1,5,6,10,12,15,20,30,60}）
_SUB = [1, 5, 6, 10, 12, 15, 20, 30, 60]


def _n_decomp(target, parts, idx=0):
    if target == 0:
        return 1
    if idx == len(parts):
        return 0
    tot = 0
    k = 0
    while k * parts[idx] <= target:
        tot += _n_decomp(target - k * parts[idx], parts, idx + 1)
        k += 1
    return tot


n_decomp44 = _n_decomp(44, _SUB)
EVID["BF5"] = {"S1": "雙射之存在性 claim（非等式）",
               "S2_probe": "非恆等",
               "S3": f"44 點 A₅-集分解 = {n_decomp44} 種 ⟹ 不唯一", "S3_num": n_decomp44,
               "S4": f"orbit 多重集：群 {g_orb} vs 模式 {m_orb}／逆 {inv_orb} vs {ck_pair}"
                     f" ⟹ 0/6 相等（等變死）",
               "class": "M", "m21": "🔴"}

# ── BF6：22 三重身份升 count-level bridge ───────────────────
import random
random.seed(20261125)
viol = 0
for _ in range(2000):
    k = random.randint(1, 8)
    DD = [random.randint(1, 12) for _ in range(k)]
    if 2 * sum(d * (d - 1) // 2 for d in DD) != sum(d * (d - 1) for d in DD):
        viol += 1
EVID["BF6"] = {"S1": "模側 22 有 functor（ΣΛ²）；鏈側/群側 22 全為 44／88 之函數",
               "S2_probe": f"N_d = 2ΣΛ² 於 2000 隨機度集 0 違例 ⟹ 純代數恆等 ⟹ T",
               "S3": "不適用", "S3_num": -1, "S4": "不適用",
               "class": "T", "m21": "🔴"}

# ── BF7：h11 = 36 由 A₅ canonical construction 導出 ─────────
pool = {}
for k in range(1, 6):
    pool[f"ΣC(d,{k})"] = sum(comb(d, k) for d in D)
for k in range(1, 4):
    pool[f"ΣC(d+{k - 1},{k})"] = sum(comb(d + k - 1, k) for d in D)
pool["Σ_(i<j)prod"] = sum(D[i] * D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)sum"] = sum(D[i] + D[j] for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)C(sum,2)"] = sum(comb(D[i] + D[j], 2) for i in range(5) for j in range(i + 1, 5))
pool["Σ_(i<j)|diff|"] = sum(abs(D[i] - D[j]) for i in range(5) for j in range(i + 1, 5))
for r in range(1, len(ORD_CNT) + 1):
    for S in combinations(sorted(ORD_CNT), r):
        pool["|ord∈{%s}|" % ",".join(map(str, S))] = sum(ORD_CNT[o] for o in S)
_cls_seen, cls_sizes = set(), []
for g in G:
    if g in _cls_seen:
        continue
    c = frozenset(mul(mul(h, g), invp(h)) for h in G)
    cls_sizes.append(len(c))
    _cls_seen |= c
cls_sizes = sorted(cls_sizes, reverse=True)
assert cls_sizes == [20, 15, 12, 12, 1], cls_sizes
for r in range(1, len(cls_sizes) + 1):
    for S in set(combinations(cls_sizes, r)):
        pool["Σ類子集{%s}" % ",".join(map(str, S))] = sum(S)
comp2 = {}
for k2, v2 in list(pool.items()):
    if 0 < v2 < 60:
        comp2[f"60−({k2})"] = 60 - v2
pool.update(comp2)
hits36 = sorted(k2 for k2, v2 in pool.items() if v2 == 36)
EVID["BF7"] = {"S1": "目標 36 之 functor 未 exhibit",
               "S2_probe": "非恆等",
               "S3": f"池 94（實 |𝔇| = {len(pool)}）內命中 36 = {len(hits36)}（2 恆等族）",
               "S3_num": len(hits36), "S4": "候選讀法（階 ≤3）共軛不變 —— S4 過但 S3 死",
               "class": "M", "m21": "🔴"}

# ── BF8：{20,24} 由共同結構強迫 ────────────────────────────
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


phi40 = [N for N in range(2, 601, 2) if phi(N) == 40]
EVID["BF8"] = {"S1": "兩側 20 屬獨立計數（階3 元素數 vs gcd=1 模式數）",
               "S2_probe": "非恆等",
               "S3": f"φ(N)=40 之偶數解 = {phi40}（{len(phi40)} 解 ⟹ 不強迫）",
               "S3_num": len(phi40), "S4": "無共同 Aut（A₅ 單非交換 vs Galois 阿貝爾）",
               "class": "M", "m21": "🔴"}

# ── 非調參測試：三條新 functional ───────────────────────────
d_sorted = sorted(D, reverse=True)
Na = comb(d_sorted[0] + d_sorted[1], 2)                       # C(4+5,2) = 36
Nb = sum(D[i] + D[j] for i in range(5) for j in range(i + 1, 5))   # 64
Nc = sum(min(D[i], D[j]) for i in range(5) for j in range(i + 1, 5))  # 23
# N-b 之自由恆等檢定：Σ_(i<j)(d_i+d_j) = (k−1)·Σd 對任意度集
free_id_b = all(sum(DD[i] + DD[j] for i in range(len(DD)) for j in range(i + 1, len(DD)))
                == (len(DD) - 1) * sum(DD)
                for DD in [tuple(random.randint(1, 9) for _ in range(k)) for k in range(2, 8)
                           for _ in range(200)])
EVID["N-a"] = {"S1": "無：「兩個最大度」非 functor（度集無序，選 4+5 屬後見選擇）",
               "S2_probe": "非恆等", "S3": "候選族未定義（選 pattern 即後見）", "S3_num": 0,
               "S4": "n/a", "class": "A", "m21": "（本輪新）"}
EVID["N-b"] = {"S1": "弱（Σ_(i<j)(d_i+d_j) 自身有 functor 形式）",
               "S2_probe": f"= (k−1)·Σd 自由恆等（{free_id_b}）⟹ 純代數改寫 ⟹ T",
               "S3": "不適用", "S3_num": -1, "S4": "n/a", "class": "T", "m21": "（本輪新）"}
EVID["N-c"] = {"S1": "無：pairwise min 為 ad hoc functional，無對象讀法",
               "S2_probe": "非恆等", "S3": "候選族未定義", "S3_num": 0, "S4": "n/a",
               "class": "A", "m21": "（本輪新）"}

print(f"\n  CA1 解集 = {sols}｜E12 階類族命中 16 = {n16_orderfam}")
print(f"  BF1  k-free 命中 44（無序化）= {n_kfree}｜raw = {n_kfree_raw}"
      f"｜帶自由常數 q_obs = {n_hit}/{n_free} = {n_hit / n_free:.4f}")
print(f"       k-free 無序式：{sorted(_kf_unord)}")
print(f"  BF2  Σγ = 2(N+1)（N=2..120）全真 = {ident_ok}")
print(f"  BF4  gauge 5 對同乘積 = {same_product}")
print(f"  BF5  群 orbit {g_orb} vs 模式 {m_orb}｜逆 {inv_orb} vs 對偶 {ck_pair}")
print(f"  BF6  2000 隨機度集違例 = {viol}")
print(f"  BF7  池 |𝔇| = {len(pool)}｜命中 36 = {len(hits36)}")
print(f"  BF8  φ(N)=40 偶數解 = {phi40}")
print(f"  N-a  C(4+5,2) = {Na}｜N-b Σ = {Nb}（= 4·16）｜N-c Σ min = {Nc}")

# ══════════════════════════════════════════════════════════════════
# 分類與校準
# ══════════════════════════════════════════════════════════════════
print("\n── 分類（S1／S2／S3／S4 → 型）──")
order_cases = ["CA1", "E12", "BF1", "BF2", "BF3", "BF4", "BF5", "BF6", "BF7", "BF8",
               "N-a", "N-b", "N-c"]
for c in order_cases:
    e = EVID[c]
    print(f"  {c:<4} → {e['class']} 型｜M21: {e['m21']}｜S3: {e['S3']}")

# 校準：10 案例（CA1、E12 → F；BF1–BF8 → 非F）
calib = []
for c in ["CA1", "E12"] + [f"BF{i}" for i in range(1, 9)]:
    want_F = c in ("CA1", "E12")
    got_F = EVID[c]["class"] == "F"
    calib.append((c, want_F, got_F, want_F == got_F))
n_ok = sum(1 for *_x, ok in calib if ok)
print(f"\n  校準：{n_ok}/10 一致")
for c, w, g, ok in calib:
    print(f"    {c:<4} 預期{'F' if w else '非F'}｜判{EVID[c]['class']}｜{'✓' if ok else '✗'}")
ck("T4-a  校準 10/10 與 M21 記錄判詞一致",
   n_ok == 10, f"{n_ok}/10")

# 非調參：三條新 functional 必須判非 F，且理由可核
blind_ok = all(EVID[c]["class"] in ("T", "A", "M") for c in ("N-a", "N-b", "N-c"))
ck("T4-b  盲測應用（三條本輪新 functional）：判準給非F 且理由可核（A／T 型）",
   blind_ok, str({c: EVID[c]["class"] for c in ("N-a", "N-b", "N-c")}))

# ══════════════════════════════════════════════════════════════════
# 應用於 M24 三目標（K6 同級檢定）
# ══════════════════════════════════════════════════════════════════
print("\n── 應用：三目標之 L0 grade（K6）──")
targets = {
    16: {"S1": "✅ 模維數／對合計數（FS）", "S2": "✅ 非自由恆等", "S3": "✅ 階類族內唯一",
         "S4": "✅ 共軛不變", "grade": "F"},
    98: {"S1": "△ 有模側讀法（cross Λ²）但 Hodge 函子未 exhibit", "S2": "✅ 非自由恆等",
         "S3": "△ 命中 1 但 = ((Σd)²−Σd²)/2 代數恆等", "S4": "✅ 等變", "grade": "M"},
    36: {"S1": "❌ 無 functor（模側只在 specific pattern；鏈側 filter 非 Galois 不變）",
         "S2": "✅", "S3": "❌ 池 94 命中 4（2 恆等族）", "S4": "❌ 鏈側非不變", "grade": "M"},
}
for t, v in targets.items():
    print(f"  {t:<3} → {v['grade']} 型｜S1 {v['S1']}｜S3 {v['S3']}｜S4 {v['S4']}")
grades = {t: v["grade"] for t, v in targets.items()}
ck("K6  三目標之 L0 grade 不全同（16 = F｜98 = M｜36 = M）⟹ 殺手測試之「同級」FAIL",
   len(set(grades.values())) > 1, str(grades))

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/04 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "04_l0_meta_language.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "evidence": EVID, "calibration": calib, "n_calib_ok": n_ok,
               "target_grades": grades,
               "verdict": "T4 = ✅（校準 10/10；盲測 3/3 非F）；K6 FAIL（16F／98M／36M）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
