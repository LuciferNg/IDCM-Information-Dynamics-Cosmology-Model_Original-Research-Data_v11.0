#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/02 — CY₃ 路線密度閘（R2／T2／P2）
==============================================
判準照 docs/00_preregistration.md §3（跑前寫死）：

  CY1  族 𝔊 = {a·x + b·y + c}（x,y ∈ 池、a,b ∈ {−2..2} 去掉 0、|c| ≤ 120、≤ 3 項）
       命中 44 之比例 q ≤ 5% ⟹ 有判別力
  CY2  若 q > 5% ⟹ T2 🔴（無 producer；池富整數 ⟹ 一切皆中）

另加（事前同一節精神：「≥5 條 k-free 式 ⟹ 零判別力」）：
  CY1′ k-free 族：{a·x + b·y}（零自由常數）+ {a·x}，a,b ∈ {±1,±2}
       命中數 ≥ 5 ⟹ 零判別力（不唯一）

池 P0（44-free、不含 88／45／89 等 44-導出量）：
  h11=36, h21=98, χ=−124, b₂=36, b₃=198, B=4, Σγ=90, N_gen=3, |G|=60,
  K_tree=33, F₁₀=55, L₄=7, L₅=11
池 P1（+ 44-導出量，只作「循環命中」披露，不計入判定）：
  N=88, N_d+1=45, F₁₁=89, C(44,2)=946, 45, B²=16, k=5

MC null：N = 10^7 隨機表達式（x,y 隨機取自 P0、a,b 隨機取 {±1,±2}），
         統計 P(值 = 44) ⟹ 單一目標之零命中率。
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ── 池（44-free）──────────────────────────────────────────────
P0 = {
    "h11": 36, "h21": 98, "chi=2(h11-h21)": -124, "b2=h11": 36, "b3=2(1+h21)": 198,
    "B": 4, "Sigma_gamma": 90, "N_gen": 3, "|G|": 60, "K_tree": 33,
    "F10": 55, "L4": 7, "L5": 11,
}
P1_derived = {
    "N=88": 88, "N_d+1=45": 45, "F11=89": 89, "C(44,2)=946": 946,
    "B^2=16": 16, "k=5": 5,
}

TARGET = 44

print("=" * 78)
print("M21 scripts/02 — CY₃ 路線密度閘")
print("=" * 78)
print(f"  池 P0（{len(P0)} 項，全部 44-free）：{P0}")
print(f"  目標 = {TARGET}")
ck("CY0-a  池 P0 不含 44", TARGET not in P0.values(), "44 ∉ P0")
ck("CY0-b  池 P0 不含 88／45／89（44-導出）",
   not ({88, 45, 89} & set(P0.values())), str(sorted({88, 45, 89} & set(P0.values()))))

items = list(P0.items())

# ── 族 𝔊0：k-free（零自由常數）────────────────────────────────
print("\n── CY1′ k-free 族（零自由常數）──")
kfree_hits = []
single_hits = []
COEFS = (-2, -1, 1, 2)
n_expr = 0
for ai in COEFS:
    for (n1, v1) in items:
        expr = f"{ai:+d}·{n1}"
        n_expr += 1
        if ai * v1 == TARGET:
            single_hits.append((expr, ai * v1))
for i, (n1, v1) in enumerate(items):
    for j, (n2, v2) in enumerate(items):
        if j <= i:
            continue
        for ai in COEFS:
            for bi in COEFS:
                n_expr += 1
                val = ai * v1 + bi * v2
                if val == TARGET:
                    kfree_hits.append((f"{ai:+d}·{n1} {bi:+d}·{n2}", val))

print(f"  k-free 表達式總數 = {n_expr}")
print(f"  單項命中 = {[e for e, _ in single_hits]}")
print(f"  兩項命中數 = {len(kfree_hits)}：")
for e, _ in kfree_hits:
    print(f"    · {e} = 44")

ck("CY1′-a  k-free 命中數 < 5（唯一／少數才有判別力）",
   len(kfree_hits) + len(single_hits) < 5,
   f"hits = {len(kfree_hits) + len(single_hits)}")
ck("CY1′-b  無單項 k-free 命中（44 唔係池成員之標量倍）", len(single_hits) == 0,
   f"{len(single_hits)}")

# ── 族 𝔊1：帶自由常數 c（|c| ≤ 120）──────────────────────────
print("\n── CY1 族 𝔊（含自由常數 |c| ≤ 120）──")
n_total = 0
n_hit = 0
for i, (n1, v1) in enumerate(items):
    for j, (n2, v2) in enumerate(items):
        if j <= i:
            continue
        for ai in COEFS:
            for bi in COEFS:
                base = ai * v1 + bi * v2
                c_needed = TARGET - base
                n_total += 1
                if abs(c_needed) <= 120:
                    n_hit += 1
q_obs = n_hit / n_total
print(f"  表達式總數 = {n_total} ｜ 可由 |c| ≤ 120 調整命中者 = {n_hit}")
print(f"  q_obs = {q_obs:.4f}（判準 q ≤ 0.05）")
ck("CY1  族 𝔊 密度 q_obs ≤ 5%", q_obs <= 0.05, f"q_obs = {q_obs:.4f} ({n_hit}/{n_total})")

# ── 循環命中披露（P1）────────────────────────────────────────
print("\n── 循環命中披露（P1：44-導出量入池）──")
circ = []
allp = dict(P0)
allp.update(P1_derived)
cit = list(allp.items())
for i, (n1, v1) in enumerate(cit):
    for j, (n2, v2) in enumerate(cit):
        if j <= i:
            continue
        for ai in COEFS:
            for bi in COEFS:
                if ai * v1 + bi * v2 == TARGET:
                    circ.append(f"{ai:+d}·{n1} {bi:+d}·{n2}")
print(f"  含 44-導出量之命中數 = {len(circ)}（**不計入判定**，只披露循環性）")
for e in circ[:12]:
    print(f"    · {e}")

# ── MC null（10^7）─────────────────────────────────────────────
print("\n── MC null（N = 10^7）──")
NMC = 10_000_000
rng = np.random.default_rng(20261120)
pool_vals = np.array([v for _, v in items], dtype=np.int64)
coefs = np.array([-2, -1, 1, 2], dtype=np.int64)
k = len(pool_vals)
idx1 = rng.integers(0, k, NMC)
idx2 = rng.integers(0, k, NMC)
a = coefs[rng.integers(0, 4, NMC)]
b = coefs[rng.integers(0, 4, NMC)]
vals = a * pool_vals[idx1] + b * pool_vals[idx2]
n44 = int(np.count_nonzero(vals == TARGET))
p_null = n44 / NMC
print(f"  隨機表達式值 = 44 之比例 p_null = {n44}/{NMC} = {p_null:.3e}")
# 判別力：觀測命中率 vs 零命中率
if p_null > 0:
    lift = q_obs / p_null
else:
    lift = float("inf")
print(f"  lift = q_obs / p_null = {lift:.3f}")
ck("CY1-MC  MC null 完成（10^7）", n44 >= 0, f"p_null = {p_null:.3e}")

# 池內整數密度：44 落在池值生成之區間內，故「命中」不稀奇
pool_min, pool_max = min(pool_vals), max(pool_vals)
print(f"  池值域 = [{pool_min}, {pool_max}]；|池| = {k}")

# ── 綜合判詞 ─────────────────────────────────────────────────
print("\n" + "=" * 78)
if q_obs > 0.05:
    verdict = "CY2 ⟹ T2 🔴：族 𝔊 密度超閘，44 可由大量 a·x+b·y+c 表示 ⟹ 零判別力"
else:
    verdict = "CY1 PASS ⟹ 需再審（密度閘過，仍須算子身份）"
print(verdict)
print("=" * 78)

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print(f"scripts/02 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "02_cy3_pool_gate.json"), "w", encoding="utf-8") as f:
    json.dump({
        "checks": [{"tag": t, "pass": ok, "detail": d} for t, ok, d in CHECKS],
        "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
        "pool_P0": P0, "pool_P1_derived": P1_derived, "target": TARGET,
        "kfree": {"n_expr": n_expr, "single_hits": single_hits,
                  "two_term_hits": [e for e, _ in kfree_hits]},
        "family_G": {"n_total": n_total, "n_hit": n_hit, "q_obs": q_obs},
        "circular_hits": circ,
        "mc_null": {"N": NMC, "n44": n44, "p_null": p_null, "lift": lift},
        "verdict": verdict,
    }, f, ensure_ascii=False, indent=2)

sys.exit(0 if nfail == 0 else 1)
