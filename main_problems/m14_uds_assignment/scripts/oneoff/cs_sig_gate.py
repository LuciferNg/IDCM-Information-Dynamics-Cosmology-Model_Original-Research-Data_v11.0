#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""CS／SIG 之結構推導閘
判準（跑前凍結）
 ① c == ∏_j sg_j（8 行，逐行比）
 ② 恆等式 W_table(v) == ½ Σ_{ε∈{±1}⁴} (∏ε)·[ε·v ≡ 0 (mod 90)]
    於 500,000 個隨機 4-向量 v ∈ {1..44}⁴ 逐點比（容差：整數零差）
 ③ 全域翻轉不變性：ε ↦ −ε 使 ∏ε 不變（4 腿）⟹ 16 元素和 = 2 × 8 元素和
 ④ 若 ①②③ 全過 ⟹ CS 非輸入，係「四腿符號乘積 × 環動量守恆」之強制結果
"""
import itertools
import numpy as np

N = 44
MOD = 90
SIG = np.array([(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
                (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)], dtype=np.int64)
CS = np.array([1, 1, -1, -1, -1, -1, 1, 1], dtype=np.int64)

print("【①】逐行核 c == ∏ sg")
ok = 0
for i, (c, sg) in enumerate(zip(CS, SIG)):
    p = int(np.prod(sg))
    ok += (p == c)
    print("  row {}：sg = {}｜∏ = {:+d}｜c = {:+d}｜{}".format(
        i + 1, tuple(int(x) for x in sg), p, int(c), "✓" if p == c else "✗"))
print("  ⟹ 命中 {}/8".format(ok))

E4 = np.array(list(itertools.product([1, -1], repeat=4)), dtype=np.int64)
PE = np.prod(E4, axis=1)

def W_table(T):
    s = np.zeros(T.shape[0], dtype=np.int64)
    for c, sg in zip(CS, SIG):
        s += np.where(np.abs(T @ sg) % MOD == 0, c, 0)
    return s

def W_ident(T):
    s = np.zeros(T.shape[0], dtype=np.int64)
    for pe, e in zip(PE, E4):
        s += np.where(np.abs(T @ e) % MOD == 0, pe, 0)
    return s // 2

rng = np.random.default_rng(20261124)
NS = 500_000
T = rng.integers(1, N + 1, size=(NS, 4), dtype=np.int64)
wt = W_table(T)
wi = W_ident(T)
bad = int((wt != wi).sum())
print("\n【②】恆等式 W_table == ½Σ(∏ε)[ε·v≡0]")
print("  樣本 = {:,}｜不符 = {:,}｜max|差| = {}".format(NS, bad, int(np.abs(wt - wi).max())))
u, c = np.unique(wt, return_counts=True)
print("  W 值分佈（樣本）：", {int(a): int(b) for a, b in zip(u, c)})

print("\n【③】全域翻轉不變性（∏ε 於 4 腿下不變）")
flip = (PE == np.prod(-E4, axis=1)).all()
print("  ∏(−ε) == ∏ε 逐元素：{}".format(bool(flip)))
same = []
for pe, e in zip(PE, E4):
    m = np.abs(T @ e) % MOD == 0
    m2 = np.abs(T @ (-e)) % MOD == 0
    same.append(bool((m == m2).all()))
print("  條件 [|ε·v|≡0] 對 ε 與 −ε 同值：{}（8/8）".format(all(same)))
print("\n【④】判詞")
print("  ① {}/8｜② 不符 {:,}（0 則恆等成立）｜③ {} ⟹ {}".format(
    ok, bad, all(same),
    "CS 非輸入，係四腿符號乘積 × 環動量守恆之強制結果" if (ok == 8 and bad == 0 and all(same)) else "未閉合"))

import json
json.dump({
 "rows_c_equals_prod": int(ok),
 "identity_sample": NS, "identity_bad": bad,
 "W_distribution_sample": {str(int(a)): int(b) for a, b in zip(u, c)},
 "flip_invariance": bool(all(same)),
 "E4_16_terms": [list(map(int, e)) for e in E4],
 "PE": [int(x) for x in PE],
}, open("/home/wsl/idcm_ckpt/m4_cs_sig.json", "w"), indent=1, ensure_ascii=False)
print("\n已寫 /home/wsl/idcm_ckpt/m4_cs_sig.json")
