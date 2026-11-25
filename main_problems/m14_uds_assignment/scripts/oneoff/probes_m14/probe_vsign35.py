#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scripts/48_fifth_order.py — MR1 機器（七）：五階微擾 E_5（全扇區，原有模式）

性質
    零參數、零觀測、無擬合。純代數／組合計算。
    沿用 docs/41–46 之既有機器（MR1），**不引入任何新假設**。

公約（與 docs/41–47 完全一致，勿混用）
    /90 站空間：φ_n(j) = √(2/45) sin(πnj/45)，j = 1..44
    γ_n = 4 sin²(πn/90)，n = 1..44，Σγ = 88，γ_max = γ_44 = 3.99512810052
    屏蔽模態：1↔3（p=1, q=3, 權 6）、3↔3（p=2, q=2, 權 4）、
              3↔5（p=1, q=3, 權 6）、5↔5（p=2, q=2, 權 4）
    模式數變化：Δ|T| ∈ {0, ±2}（原有模式；Δ=±4 之 p=0/4 通道**不在本腳本範圍**，
    見 docs/47 §七 之範圍宣告與 register note）

N_trials
    本腳本為**完整枚舉**（deterministic），無 MC、無隨機抽樣、無誤差棒。
    §二 之對稱核抽樣為**固定 40 個樣本**（不隨機），僅作算子閘。

範圍宣告
    計算 E_2 … E_5 於 n ∈ {1, 17, 43, 44}（代表：最輕／E_4 最大／兩個最重）。
    全部 44 個 n 未算（成本 5→5 ≈ 250–320 s/n）。
    本腳本**不**宣稱五階完整：p=0（四創建）通道未含（凍結，見 register）。

輸出紀律
    stdout 不含任何計時、時間、pid、集合迭代序；耗時一律 stderr。
    兩跑必須 byte-identical（審計 ⑦）。
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import itertools
import math
import sys
import time
from collections import Counter
from math import comb, factorial

import numpy as np
import scipy.sparse as sp

T0 = time.time()
L = "=" * 78
N = 44
ALL = "--all" in sys.argv          # 全 44 n 掃描（預設關；預設路徑輸出不得改變）

SIG = np.array([(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
                (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)], dtype=np.int32)
CS = np.array([1, 1, -1, -1, -1, -1, 1, 1], dtype=np.int32)


def W4v(T):
    """χ⁴ = W/90 之分子 W，對一批 4-index 組合向量化計算"""
    s = np.zeros(T.shape[0], dtype=np.int32)
    for c, sg in zip(CS, SIG):
        s += np.where(np.abs(T @ sg) % 90 == 0, c, 0)
    return s


def nrm(t):
    r = 1
    for v in Counter(t).values():
        r *= factorial(v)
    return float(r)


G = np.array([0.0] + [4 * np.sin(np.pi * n / 90) ** 2 for n in range(1, N + 1)])


def make_rank(k):
    """多重集 (a1≤…≤ak) 於 itertools.combinations_with_replacement 之 lex rank"""
    M = N + k - 1
    CT = np.zeros((M + 2, k + 2), dtype=np.int64)
    for nn in range(M + 2):
        for rr in range(k + 2):
            CT[nn, rr] = comb(nn, rr)
    ARR = np.arange(k)

    def f(a):
        a = np.asarray(a, dtype=np.int64)
        b = a + ARR
        prev = np.zeros(len(a), dtype=np.int64)
        out = np.zeros(len(a), dtype=np.int64)
        for i in range(k):
            r = k - i - 1
            out += CT[M - prev, r + 1] - CT[M - b[:, i] + 1, r + 1]
            prev = b[:, i]
        return out

    return f


rank2 = make_rank(2)
rank3 = make_rank(3)
rank5 = make_rank(5)

# ------------------------------------------------------------------ §一
print(L)
print("【§一】五階 RS 之結構：模式數可達性（原有模式 Δ|T| ∈ {0,±2}）")
print(L)


def reach(k, deltas):
    best = 0
    paths = []

    def dfs(m, step, trail):
        nonlocal best
        if step == k:
            if m == 1:
                paths.append(tuple(trail))
            return
        for d in deltas:
            if m + d >= 1:
                dfs(m + d, step + 1, trail + [m + d])

    dfs(1, 0, [1])
    return (max(max(p) for p in paths) if paths else 0), paths


for k in (2, 3, 4, 5, 6):
    mx, ps = reach(k, [0, 2, -2])
    ex = sorted(set(p for p in ps if max(p) == mx))[:1]
    print(f"  E_{k}（{k} 頂點）：最大可達模數 = {mx}｜鏈例 {ex}")
print("  ⟹ 四階最大 5 模（1→3→5→3→1）⟹ docs/45+46 已完整；五階亦最大 5 模；7 模屬 E_6")

# ------------------------------------------------------------------ §二
print()
print(L)
print("【§二】建 3↔5 算子（p=1, q=3）＋ 對稱核閘")
print(L)
TRI = [tuple(t) for t in itertools.combinations_with_replacement(range(1, N + 1), 3)]
IDX3 = {t: i for i, t in enumerate(TRI)}
N3 = len(TRI)
NT3 = np.array([nrm(t) for t in TRI])
ET3 = np.array([sum(G[a] for a in t) for t in TRI])
QP = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 5)), dtype=np.int32)
NQ = len(QP)
NT5 = np.array([nrm(map(int, t)) for t in QP])
EQ = np.array([float(np.sum(G[t])) for t in QP])
ADD2 = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 2)), dtype=np.int32)
ADD3 = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 3)), dtype=np.int32)
print(f"  三模態 {N3:,}｜五模態 {NQ:,}｜加對 {len(ADD2)}｜加三 {len(ADD3)}")

qidx_l, vals35_l, tri_l = [], [], []
for i, T in enumerate(TRI):
    for pos in range(3):
        surv = [T[q] for q in range(3) if q != pos]
        rem = int(T[pos])
        A = np.column_stack([np.full(len(ADD3), rem, dtype=np.int32), ADD3])
        w = W4v(A)
        nz = np.nonzero(w)[0]
        if len(nz) == 0:
            continue
        Qs = np.sort(np.column_stack([np.repeat(np.array(surv, dtype=np.int32)[None, :], len(nz), axis=0),
                                      ADD3[nz]]), axis=1)
        r = rank5(Qs)
        cs = Counter(surv)
        emb = np.ones(len(nz))
        for m, cnt in cs.items():
            matches = (Qs == m).sum(1)
            emb *= np.array([comb(int(x), cnt) for x in matches], dtype=float)
        amp = 6.0 * (w[nz] / 90.0) * emb / (NT3[i] * NT5[r]) ** 0.5
        qidx_l.append(r.astype(np.int32))
        vals35_l.append(amp)
        tri_l.append(np.full(len(r), i, dtype=np.int32))
qidx = np.concatenate(qidx_l)
vals35 = np.concatenate(vals35_l)
tri = np.concatenate(tri_l)
print(f"  3↔5 非零元 = {len(qidx):,}｜每 trio 均值 {len(qidx) / N3:.1f}")


# ★ 凍結 3↔5 算子之陣列（防止後段 5↔5 LUT 重用同名 tri/qidx 造成污染）
TRI35 = tri.copy(); QIDX35 = qidx.copy(); VALS35F = vals35.copy()
print("  [凍結] 3↔5：tri35 {}｜qidx35 {}｜vals35 {}".format(TRI35.shape, QIDX35.shape, VALS35F.shape))


def apply_35(p3):
    return np.bincount(QIDX35, weights=VALS35F * p3[TRI35], minlength=NQ)


def apply_53(p5):
    return np.bincount(TRI35, weights=VALS35F * p5[QIDX35], minlength=N3)


def mel35_tside(T, Q):
    """由 T 側（p=1: 去 1 添 3）獨立直算 ⟨Q|H/g₄|T⟩"""
    cQ = Counter(Q)
    tot = 0.0
    for pos in range(3):
        rem = [T[pos]]                                  # 被移除之 1 模（p=1）
        sv = [T[q] for q in range(3) if q != pos]       # 倖存 2 模
        cS = Counter(sv)
        if any(cQ[m] < c for m, c in cS.items()):
            continue
        add = list((cQ - cS).elements())
        if len(add) != 3:
            continue
        w = sum(c for c, sg in zip(CS, SIG)
                if sum(x * y for x, y in zip(sg, tuple(rem + add))) % 90 == 0)
        if not w:
            continue
        emb = 1
        for m, c in cS.items():
            emb *= comb(cQ[m], c)
        tot += 6 * (w / 90.0) * emb / (nrm(T) * nrm(Q)) ** 0.5
    return tot


def col35(T):
    """由 3↔5 表取出 ⟨Q|H/g₄|T⟩ 對全部 Q（稀疏）"""
    m = tri == IDX3[T]
    return qidx[m], vals35[m]


SAMP = [((2, 9, 33), (2, 5, 14, 27, 41)), ((7, 7, 12), (7, 12, 20, 20, 33)),
        ((1, 1, 1), (1, 1, 2, 3, 5)), ((5, 17, 44), (5, 17, 18, 29, 40)),
        ((11, 22, 33), (11, 22, 33, 40, 41)), ((3, 15, 27), (3, 15, 27, 27, 27)),
        ((4, 8, 16), (4, 8, 16, 24, 32)), ((6, 6, 18), (6, 6, 18, 18, 30)),
        ((9, 18, 36), (9, 18, 27, 36, 44)), ((13, 26, 39), (13, 26, 39, 39, 44))]
mx35 = 0.0
nz35 = 0
for T, Q in SAMP:
    qi_, vv_ = col35(T)
    got = 0.0
    for j, v in zip(qi_, vv_):
        if tuple(int(x) for x in QP[j]) == Q:
            got += float(v)
    ref = mel35_tside(T, Q)
    mx35 = max(mx35, abs(got - ref))
    nz35 += int(ref != 0)
print(f"  對稱核（10 個固定樣本、{nz35} 個非零）：max|表 − 直算| = {mx35:.3e}")


# ================================================== 3↔5 頂點符號律（vals35 之 sign 結構）
# 目的：判定 sign(V^{3↔5}) 是否只依賴「改變之四模向量」v = rem ⊕ add（mod 90 結構）
#   amp = 6·(w/90)·emb/√(NT3·NT5)，emb = 1 ⟹ w = amp·90·√(NT3 NT5)/6
#   w = Σ_c c·[⟨sg, v⟩ ≡ 0 (mod 90)]
import json as _json
print("\n【§V】3↔5 頂點符號律")
amp = VALS35F
den = np.sqrt(NT3[TRI35] * NT5[QIDX35])
w = amp * 90.0 * den / 6.0
print("  非零元 = {:,}".format(len(amp)))
print("  w 之 min / max = {:+.6e} / {:+.6e}".format(float(w.min()), float(w.max())))
rw = np.rint(w)
print("  max|w − round(w)| = {:.3e} ⟹ w 整數性 {}".format(float(np.abs(w - rw).max()),
      "✓" if np.abs(w - rw).max() < 1e-9 else "✗"))
u, c = np.unique(rw.astype(np.int64), return_counts=True)
print("  w 之相異值 {} 個：".format(len(u)))
for uu, cc in zip(u, c):
    print("     w = {:+3d} ⟹ {:>12,} 元（{:6.3f}%）".format(int(uu), int(cc), 100.0 * cc / len(amp)))
neg = rw < 0
print("\n  負 w 之元數 = {:,}（{:.4f}%）".format(int(neg.sum()), 100.0 * neg.sum() / len(amp)))
print("  sign(amp) == sign(w) 一致率 = {:.6f}".format(float((np.sign(amp) == np.sign(rw)).mean())))
# ★ 是否只依賴 v：抽樣反算 v 並重算 w
rng = np.random.default_rng(20261124)
idx = rng.choice(len(amp), size=200000, replace=False)
import collections
def mult_diff(T, Q):
    cT = collections.Counter(T); cQ = collections.Counter(Q)
    surv = cT & cQ
    rem = list((cT - surv).elements()); add = list((cQ - surv).elements())
    return tuple(sorted(rem)), tuple(sorted(add))
def w_of(v):
    tot = 0
    for cc, sg in zip(CS, SIG):
        if sum(int(x) * int(y) for x, y in zip(sg, v)) % 90 == 0:
            tot += int(cc)
    return tot
bad = 0; nmatch = 0; negcls = collections.Counter(); poscls = collections.Counter()
for k in idx:
    T = TRI[TRI35[k]]; Q = QP[QIDX35[k]]
    rem, add = mult_diff(T, Q)
    v = tuple(rem) + tuple(add)          # ★ 保持 W4v 之引數次序（rem 在前）
    ww = w_of(v)
    # emb = ∏_m C(#Q 中 m 之數, #surv 中 m 之數)（正值）
    cQ = collections.Counter(Q); cS = collections.Counter([T[q] for q in range(3) if T[q] != rem] if False else [])
    survm = collections.Counter(T); 
    for _m in list(T):
        survm[_m] -= 1
    survm = +survm
    embv = 1
    for _m, _cnt in survm.items():
        embv *= int(math.comb(cQ[_m], _cnt))
    if ww * embv != int(rw[k]):
        bad += 1
    else:
        nmatch += 1
    if len(rem) == 1 and len(add) == 3:
        (negcls if ww < 0 else poscls)[(rem, add)] += 1
print("\n  抽樣 {:} 元：w 由 v=rem⊕add 重算命中 = {:,}｜不符 = {:,}".format(len(idx), nmatch, bad))
print("  （不符 = 0 ⟹ 振幅之符號只由改變之四模向量 v 決定，與 (T,Q) 之其餘內容無關）")
print("\n  唯一類別數：負 w 類 = {}｜正 w 類 = {}".format(len(negcls), len(poscls)))
print("  負 w 類之前 12 個（rem, add）：")
for (rem, add), cc in negcls.most_common(12):
    print("     rem={} add={} ⟹ {} 次（w={}）".format(rem, add, cc, w_of(tuple(rem)+tuple(add))))
_json.dump({
 "n_nonzero": int(len(amp)), "w_min": float(w.min()), "w_max": float(w.max()),
 "w_values": {str(int(uu)): int(cc) for uu, cc in zip(u, c)},
 "negative_count": int(neg.sum()), "negative_fraction": float(neg.sum()/len(amp)),
 "sign_consistent": float((np.sign(amp) == np.sign(rw)).mean()),
 "sample": int(len(idx)), "sample_match": int(nmatch), "sample_bad": int(bad),
 "neg_classes": int(len(negcls)), "pos_classes": int(len(poscls)),
 "neg_top": [[list(map(int, rem)), list(map(int, add)), int(cc)] for (rem, add), cc in negcls.most_common(20)],
}, open("/home/wsl/idcm_ckpt/e8/mrB/vsign35.json", "w"), indent=1, ensure_ascii=False)
print("\n已寫 /home/wsl/idcm_ckpt/e8/mrB/vsign35.json")
