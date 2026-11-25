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
MRB_TERMS = []

# ============================================================ L3 定律之代數檢定
# 目標：sign(t2) = −sign(E_2) 44/44 之【充分條件檢定】
#   t2 = inv3·inv5·⟨b|apply_53(v2_5 − E2 p2_5)⟩
#      = inv3³·inv5·[ P1 − inv5·R·P2 ]
#   P1 = ⟨b|apply_53(apply_35(V33 b))⟩ > 0
#   P2 = ⟨b|apply_53(apply_35 b)⟩       > 0
#   R  = ‖b‖² > 0
#   ⟹ sign(t2) = sign(inv3)·sign(inv5)·sign(P1 − inv5·R·P2)
# 判準（跑前凍結）
#   ① vals35 ≥ 0（算子全非負）⟹ 逐元素 min ≥ 0
#   ② 代數分解對 t2 之相對差 ≤ 1e−10
#   ③ pred ≡ sign(t2) 之命中數（44 n）
#   ④ pred ≡ −sign(E2) 之命中數（44 n）
import json as _json
print("\n【§L3】非負性與代數分解檢定")
print("  vals35 之 min = {:.6e}｜qidx/tri 長度 = {:,}".format(float(vals35.min()), len(vals35)))
T44 = {r["n"]: r for r in _json.load(open("/home/wsl/idcm_ckpt/e8/mrB/terms44.json"))}
D44 = {r["n"]: r for r in _json.load(open("/home/wsl/idcm_ckpt/e8/mrB/decomp48_all_with_E5no.json"))}
print("\n   n ｜ sign(inv3) ｜ sign(inv5) ｜      R      ｜      P1      ｜      P2      ｜ P1/(R P2) ｜ pred ｜ t2 實 ｜ E2 實 ｜ ① ｜ ②")
rows = []
okA = okB = 0; mxrel = 0.0
for n in range(1, 45):
    A = np.column_stack([np.full(N3, n, dtype=np.int32), np.array(TRI, dtype=np.int32)])
    b = 6 * (W4v(A) / 90.0) / np.sqrt(NT3)
    R = float(b @ b)
    inv3 = 1.0 / (G[n] - ET3)
    inv5 = 1.0 / (G[n] - EQ)
    E2 = float(b @ (inv3 * b))
    P1 = float(b @ apply_53(apply_35(V33 @ b)))
    P2 = float(b @ apply_53(apply_35(b)))
    t2 = inv3**3 * inv5 * (P1 - inv5 * R * P2)
    t2a = T44[n]["t2"]
    rel = abs(t2 - t2a) / max(abs(t2a), 1e-300)
    mxrel = max(mxrel, rel)
    pred = np.sign(inv3) * np.sign(inv5) * np.sign(P1 - inv5 * R * P2)
    c1 = int(pred == np.sign(t2a))
    c2 = int(pred == -np.sign(D44[n]["E2"]))
    okA += c1; okB += c2
    rows.append(dict(n=n, inv3=inv3, inv5=inv5, R=R, P1=P1, P2=P2, kap=P1 / (R * P2),
                     t2_pred=t2, t2_act=t2a, rel=rel, E2=D44[n]["E2"], pred=int(pred)))
    print("  {:2d} ｜ {:+.0f} ｜ {:+.0f} ｜ {:9.3e} ｜ {:9.3e} ｜ {:9.3e} ｜ {:8.4f} ｜ {:+.0f} ｜ {:+.0f} ｜ {:+.4f} ｜ {} ｜ {:.1e}".format(
        n, np.sign(inv3), np.sign(inv5), R, P1, P2, P1 / (R * P2), pred, np.sign(t2a), D44[n]["E2"], c1, rel))
print("\n  ★ 代數分解對 t2 之 max 相對差 = {:.3e}（判準 ≤1e−10）".format(mxrel))
print("  ★ pred == sign(t2) 命中 = {}/44".format(okA))
print("  ★ pred == −sign(E2) 命中 = {}/44".format(okB))
# ⑥ 充分條件：inv5<0 ⟹ 自動；inv5>0 ⟹ 需 inv5·R·P2 > P1 ⟺ γ_n − EQ < P1/(R P2)
auto = [r["n"] for r in rows if r["inv5"] < 0]
need = [r["n"] for r in rows if r["inv5"] > 0]
fail = [r["n"] for r in rows if r["inv5"] > 0 and not (r["inv5"] * r["R"] * r["P2"] > r["P1"])]
print("\n  inv5 < 0（自動成立）之 n = {} 個｜inv5 > 0（需不等式）之 n = {} 個".format(len(auto), len(need)))
print("  inv5 > 0 而不等式失敗之 n = {}".format(fail or "無"))
print("  ⟹ 判詞：{}".format("充分條件 44/44 成立（定律可證）" if not fail else "存在反例，定律不可由此路證"))
_out = {
 "note": "docs/49 §十：sign(t2) = −sign(E_2) 之代數分解與充分條件檢定",
 "vals35_min": float(vals35.min()),
 "max_rel_vs_actual_t2": mxrel,
 "hits_pred_vs_t2": okA, "hits_pred_vs_negE2": okB,
 "inv5_negative_n": auto, "inv5_positive_n": need, "inequality_failures": fail,
 "rows": rows,
}
_json.dump(_out, open("/home/wsl/idcm_ckpt/e8/mrB/law_L3.json", "w"), indent=1, ensure_ascii=False)
print("\n已寫 /home/wsl/idcm_ckpt/e8/mrB/law_L3.json")
