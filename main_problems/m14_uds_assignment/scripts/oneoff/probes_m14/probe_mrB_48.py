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


def apply_35(p3):
    return np.bincount(qidx, weights=vals35 * p3[tri], minlength=NQ)


def apply_53(p5):
    return np.bincount(tri, weights=vals35 * p5[qidx], minlength=N3)


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

# ------------------------------------------------------------------ §三
print()
print(L)
print("【§三】查表：χ⁴(990×990)、三元組 (rank, emb) LUT、V33")
print(L)
t0 = time.time()
W = np.zeros((990, 990), dtype=np.int16)
for i in range(990):
    A = np.column_stack([np.full(990, ADD2[i, 0], dtype=np.int32),
                         np.full(990, ADD2[i, 1], dtype=np.int32), ADD2])
    W[i] = W4v(A)
LUTR = np.zeros((N3, 990), dtype=np.int32)
LUTE = np.zeros((N3, 990), dtype=np.uint8)
for i, tr in enumerate(TRI):
    Tp = np.sort(np.concatenate([np.repeat(np.array(tr, dtype=np.int32)[None, :], 990, axis=0),
                                 ADD2], axis=1), axis=1)
    LUTR[i] = rank5(Tp).astype(np.int32)
    emb = np.ones(990, dtype=np.int64)
    for kk in range(3):
        v = tr[kk]
        cnt = (Tp == v).sum(1)
        prior = sum(1 for k2 in range(kk) if tr[k2] == v)
        emb *= (cnt - prior)
    LUTE[i] = emb.astype(np.uint8)
print(f"  χ⁴ 表非零率 {(W != 0).mean():.4f}｜LUT {LUTR.shape}｜耗時 {time.time() - t0:.0f} s", file=sys.stderr)

rows, cols, vals = [], [], []
for i, T in enumerate(TRI):
    for pos in range(3):
        surv = T[pos]
        rem = [T[q] for q in range(3) if q != pos]
        A = np.column_stack([np.full(len(ADD2), rem[0], dtype=np.int32),
                             np.full(len(ADD2), rem[1], dtype=np.int32), ADD2])
        w = W4v(A)
        for kk in np.nonzero(w)[0]:
            Tp = tuple(sorted([surv, int(ADD2[kk, 0]), int(ADD2[kk, 1])]))
            j = IDX3[Tp]
            rows.append(i)
            cols.append(j)
            vals.append(4 * (float(w[kk]) / 90.0) * Counter(Tp)[surv] / (NT3[i] * NT3[j]) ** 0.5)
V33 = sp.csr_matrix((vals, (rows, cols)), shape=(N3, N3))
asym = abs(V33 - V33.T).max()
print(f"  V33：非零 {V33.nnz:,}｜max|V−Vᵀ| = {asym:.3e}")
PAIRS5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
SURV5 = {p: [q for q in range(5) if q not in p] for p in PAIRS5}


def apply_55(p5, chunk=400, only=None):
    """5↔5（p=2, q=2）：查表法；only = (start, stop) 可限段（供 δ 測試）"""
    out = np.zeros(NQ)
    lo, hi = (0, NQ) if only is None else only
    for st in range(lo, hi, chunk):
        Qc = QP[st:st + chunk]
        m = len(Qc)
        pv = p5[st:st + chunk]
        nQin = NT5[st:st + chunk][:, None]
        for pi in PAIRS5:
            rem = np.sort(Qc[:, list(pi)], axis=1)
            pidx = rank2(rem)
            sv = Qc[:, SURV5[pi]]
            tidx = rank3(sv)
            ranks = LUTR[tidx]
            emb = LUTE[tidx].astype(np.float64)
            amp = 4.0 * (W[pidx] / 90.0) * emb / np.sqrt(nQin * NT5[ranks])
            out += np.bincount(ranks.ravel(), weights=(amp * pv[:, None]).ravel(), minlength=NQ)
    return out


def mel55_slow(Qp, Q):
    """5↔5 之獨立直算（不作任何查表）"""
    cQ = Counter(Qp)
    tot = 0.0
    for pi in PAIRS5:
        rem = [Q[i] for i in pi]
        sv = [Q[q] for q in range(5) if q not in pi]
        cS = Counter(sv)
        if any(cQ[m] < c for m, c in cS.items()):
            continue
        add = list((cQ - cS).elements())
        if len(add) != 2:
            continue
        w = sum(c for c, sg in zip(CS, SIG)
                if sum(x * y for x, y in zip(sg, tuple(rem + add))) % 90 == 0)
        if not w:
            continue
        emb = 1
        for m, c in cS.items():
            emb *= comb(cQ[m], c)
        tot += 4 * (w / 90.0) * emb / (nrm(Q) * nrm(Qp)) ** 0.5
    return tot


Q0 = (3, 11, 19, 27, 40)
j0 = int(rank5(np.array([Q0]))[0])
d = np.zeros(NQ)
d[j0] = 1.0
chunk0 = (j0 // 400) * 400
col = apply_55(d, only=(chunk0, chunk0 + 400))
cand = np.argsort(-np.abs(col))[:8]
mx55 = 0.0
for jj in cand:
    Qp = tuple(int(x) for x in QP[jj])
    ref = mel55_slow(Qp, Q0)
    mx55 = max(mx55, abs(ref - float(col[jj])))
print(f"  δ 測試（Q₀ = {Q0}）：LUT vs 直算 max|差| = {mx55:.3e}（前 8 大元）")

# ------------------------------------------------------------------ §四
print()
print(L)
print("【§四】RS 遞歸至五階（n = 1, 17, 43, 44）")
print(L)
print("  φ_1 = R₃b｜E_2 = ⟨b|φ_1⟩｜φ_2 = R V φ_1｜E_3 = ⟨b|φ_2⟩")
print("  φ_3 = R(Vφ_2 − E_2φ_1)｜E_4 = ⟨b|φ_3⟩｜φ_4 = R(Vφ_3 − E_2φ_2 − E_3φ_1)｜E_5 = ⟨b|φ_4⟩")
print()
print("   n ｜ E_2             ｜ E_3              ｜ E_4                ｜ E_5（含 5↔5）        ｜ 5↔5 修正")
RES = {}
MRB = []
for n in [int(x) for x in os.environ.get('MR_B_NS', '30,35,39,41').split(',')]:
    A = np.column_stack([np.full(N3, n, dtype=np.int32), np.array(TRI, dtype=np.int32)])
    b = 6 * (W4v(A) / 90.0) / np.sqrt(NT3)
    inv3 = 1.0 / (G[n] - ET3)
    inv5 = 1.0 / (G[n] - EQ)
    p1 = inv3 * b
    E2 = float(b @ p1)
    p2_3 = inv3 * (V33 @ p1)
    E3 = float(b @ p2_3)
    X = apply_35(p1)
    p2_5 = inv5 * X
    p3_3 = inv3 * (V33 @ p2_3 + apply_53(p2_5) - E2 * p1)
    E4 = float(b @ p3_3)
    v2_5 = apply_35(p2_3)
    p3_5 = inv5 * (v2_5 - E2 * p2_5)
    E5_no = float(b @ (inv3 * (V33 @ p3_3 + apply_53(p3_5) - E2 * p2_3 - E3 * p1)))
    t55 = time.time()
    Bp = apply_55(p2_5)
    dt55 = time.time() - t55
    # ---- MR-B §二 診斷（docs/49 §二 定義）----
    nrm25 = float(np.linalg.norm(p2_5))
    ip_v = float(p2_5 @ v2_5)
    tE2 = E2 * nrm25 * nrm25
    ipB = float(p2_5 @ Bp)
    net = ipB - tE2
    print(f"    [MR-B] ‖p₂⁵‖ = {nrm25:.6e}｜⟨p₂⁵|v₂⁵⟩ = {ip_v:+.6e}｜E_2‖p₂⁵‖² = {tE2:+.6e}"
          f"｜⟨p₂⁵|B p₂⁵⟩ = {ipB:+.6e}｜淨 Δ = {net:+.6e}｜sign·E_2 = {np.sign(net)*np.sign(E2):+.0f}")
    print(f"    [MR-B] min|DEN5| = {float(np.min(np.abs(inv5**-1))):.6e}")
    p3_5b = inv5 * (v2_5 + Bp - E2 * p2_5)
    E5 = float(b @ (inv3 * (V33 @ p3_3 + apply_53(p3_5b) - E2 * p2_3 - E3 * p1)))
    rel = (E5 - E5_no) / E5_no if E5_no else float("nan")
    MRB.append(dict(n=n, E2=E2, E3=E3, E4=E4, p25_norm=nrm25, ip_v2p25=ip_v, E2_norm2=tE2,
                    ip_B=ipB, net=net, rel=rel, sign_prod=float(np.sign(net) * np.sign(E2)),
                    den5_min=float(np.min(np.abs(inv5 ** -1)))))
    RES[n] = (E2, E3, E4, E5_no, E5, rel)
    if ALL:
        print(f"  {n:2d} ｜{E2:+.12e}｜{E3:+.12e}｜{E4:+.12e}｜{E5:+.12e}｜{rel:+.6%}")
    else:
        print(f"  {n:2d} ｜{E2:+16.9f}｜{E3:+17.9f}｜{E4:+19.6f}｜{E5:+20.6f}｜{rel:+.4%}")
    print(f"      （5↔5 未含之 E_5 = {E5_no:+.6f}；5↔5 耗時 {dt55:.0f} s）", file=sys.stderr)


import json
print()
print("=" * 78)
print("MR-B 匯總（docs/49 §二 定義；normalization = scripts/48 之 b = 6·W4v/90/√NT3）")
print("=" * 78)
print(f"  {'n':>3} ｜ {'E_2':>14} ｜ {'‖p₂⁵‖':>14} ｜ {'E_2‖p‖²':>14} ｜ {'⟨p|Bp⟩':>14} ｜ {'淨 Δ':>14} ｜ {'rel':>9}")
for r in MRB:
    print(f"  {r['n']:3d} ｜{r['E2']:+14.6f}｜{r['p25_norm']:14.4e}｜{r['E2_norm2']:+14.4e}"
          f"｜{r['ip_B']:+14.4e}｜{r['net']:+14.4e}｜{r['rel']:+9.3%}")
print()
for r in MRB:
    verdict = "符號例外（Δ 與 −E_2 同號）" if r['sign_prod'] > 0 else "符合律 sign(Δ) = −sign(E_2)"
    dom = "⟨p|Bp⟩ 主導" if abs(r['ip_B']) > abs(r['E2_norm2']) else "E_2‖p‖² 主導"
    print(f"  n={r['n']:2d}：E_2 = {r['E2']:+.6f}（{'細小' if abs(r['E2']) < 7 else '不小'}）⟹ {dom}；{verdict}")
json.dump(MRB, open(os.environ.get("MR_B_OUT", "/home/wsl/idcm_ckpt/e8/mrB/decomp48.json"), "w"),
          indent=1, ensure_ascii=False)
print("\n已寫", os.environ.get("MR_B_OUT", "/home/wsl/idcm_ckpt/e8/mrB/decomp48.json"))
