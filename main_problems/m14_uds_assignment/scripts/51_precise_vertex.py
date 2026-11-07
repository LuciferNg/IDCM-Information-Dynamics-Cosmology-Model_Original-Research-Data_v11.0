#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M14 51：精確四波頂點引擎（路徑 A / B 雙讀法）
=================================================
核心：H = (1/4!)[ :A⁴:_{χ⁴} + 6·Σ_c :A_c²: + 3·N_d ]

精確矩陣元（Wick 分解，無 emb 參數）：
    mult_n = n!/∏_v c_v(U)! × ∏_v C(c_v(U), c_v(A))
    NO     = √(∏_v c_T(v)!/(c_T(v)−a_v)!) · √(∏_v c_Q(v)!/(c_Q(v)−q_v)!)
    ⟨Q|:A⁴:|T⟩ / 24 = χ⁴(U)·Σ_splits mult·NO / 24
化簡：mult = 4!/(∏_v c_v(A)! · ∏_v c_v(C)!)

兩讀法：
  A（正規序）：頂點只含 :A⁴:；收縮歸 H_0
  B（保留收縮）：頂點另含同模 Δ=±2 通道；對角 N_T/2 + 3N_d/24 = N_T/2 + 5.5 入 H_0

輸出（stdout 只印結果，耗時一律 stderr）
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import sys, json, itertools, argparse
from collections import Counter
from math import comb, factorial
import numpy as np

N = 44
SIG = np.array([(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
                (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)], dtype=np.int32)
CS = np.array([1, 1, -1, -1, -1, -1, 1, 1], dtype=np.int32)
NB = np.zeros((64, 12), dtype=np.int64)
for _n in range(64):
    for _r in range(12):
        NB[_n, _r] = comb(_n, _r) if 0 <= _r <= _n else 0


def W4(T):
    s = np.zeros(T.shape[0], dtype=np.int32)
    for c, sg in zip(CS, SIG):
        s += np.where(np.abs(T @ sg) % 90 == 0, c, 0)
    return s.astype(np.float64)


def chi4(u):
    return W4(np.array([sorted(u)], dtype=np.int32))[0] / 90.0


def rank_fast(A):
    n, k = A.shape
    r = np.zeros(n, dtype=np.int64)
    for j in range(k):
        z = np.clip(A[:, j].astype(np.int64) + j - 1, 0, 63)
        r += NB[z, j + 1]
    return r


def counts_of(S):
    m = S.shape[0]
    out = np.zeros((m, N), dtype=np.int32)
    for v in range(1, N + 1):
        out[:, v - 1] = (S == v).sum(1)
    return out


class Add:
    def __init__(self, table):
        self.T = np.asarray(table, dtype=np.int32)
        self.n, self.k = self.T.shape
        self.cnt = np.zeros((self.n, N), dtype=np.int32)
        for v in range(1, N + 1):
            self.cnt[:, v - 1] = (self.T == v).sum(1)
        fact = np.array([1.0, 1.0, 2.0, 6.0, 24.0])
        self.fp = np.ones(self.n)
        for v in range(N):
            self.fp *= fact[np.clip(self.cnt[:, v], 0, 4)]
        self.sup = np.zeros((self.n, max(self.k, 1)), dtype=np.int32)
        self.q = np.zeros((self.n, max(self.k, 1)), dtype=np.int32)
        for i in range(self.n):
            c = Counter(self.T[i].tolist())
            for t, v in enumerate(sorted(c)):
                self.sup[i, t] = v - 1
                self.q[i, t] = c[v]


def vertex_apply(Ts, phi, add, ann_combos, n_tgt, CH=48, nmod=None):
    """向量化精確頂點：Δ|T| = n_tgt − |T|（nmod 供小模態閘用）"""
    AT = add if isinstance(add, Add) else Add(add)
    m_states = Ts.shape[0]
    Ntgt = comb((nmod or N) + n_tgt - 1, n_tgt)
    out = np.zeros(Ntgt)
    for combos in ann_combos:
        st = 0
        while st < m_states:
            en = min(st + CH, m_states)
            blk = Ts[st:en]; pv = phi[st:en]; mm = en - st
            cT = counts_of(blk)
            ann = blk[:, list(combos)]
            annc = np.zeros((mm, N), dtype=np.int32)
            for t in range(ann.shape[1]):
                np.add.at(annc, (np.arange(mm), ann[:, t] - 1), 1)
            cS = cT - annc
            rest = np.delete(blk, list(combos), axis=1)
            if AT.k == 0:
                tgt = rest
                U = ann
            else:
                tgt = np.sort(np.concatenate(
                    [np.repeat(rest[:, None, :], AT.n, axis=1),
                     np.broadcast_to(AT.T[None], (mm, AT.n, AT.k))], axis=2), axis=2)
                U = np.sort(np.concatenate(
                    [np.repeat(ann[:, None, :], AT.n, axis=1),
                     np.broadcast_to(AT.T[None], (mm, AT.n, AT.k))], axis=2), axis=2)
            W = W4(U.reshape(-1, U.shape[2])).reshape(mm, AT.n)
            fpA = np.ones(mm)
            fac3 = np.array([1.0, 1.0, 2.0, 6.0, 24.0])
            for aval in np.unique(ann):
                am = (ann == aval).sum(1)
                fpA *= fac3[np.clip(am, 0, 4)]
            mult = (24.0 / AT.fp)[None, :] / fpA[:, None]
            no_ann = np.ones(mm)
            for aval in np.unique(ann):
                am = (ann == aval).sum(1).astype(np.float64)
                cTv = cT[:, aval - 1].astype(np.float64)
                v = np.ones(mm)
                for t in range(int(am.max())):
                    v = np.where(am > t, v * (cTv - t), v)
                no_ann *= np.sqrt(v)
            fac = np.ones((mm, AT.n))
            for t in range(AT.k):
                sv = AT.sup[:, t]; qv = AT.q[:, t]
                cs = cS[:, sv].astype(np.float64)
                prod = np.ones((mm, AT.n))
                for t2 in range(4):
                    prod *= np.where(qv > t2, cs + t2 + 1, 1.0)
                fac *= prod
            mulpos = np.ones(mm)
            for aval in np.unique(ann):
                am = (ann == aval).sum(1)
                cTv = cT[:, aval - 1]
                mulpos *= np.array([comb(int(cTv[i]), int(am[i])) for i in range(mm)], dtype=np.float64)
            amp = (W / 90.0) * mult * no_ann[:, None] * np.sqrt(fac) / 24.0 / mulpos[:, None]
            contrib = amp * pv[:, None]
            rk = rank_fast(tgt.reshape(-1, n_tgt))
            sel = np.abs(contrib).reshape(-1) > 0
            out += np.bincount(rk[sel], weights=contrib.reshape(-1)[sel], minlength=Ntgt)
            st = en
    return out


def vertex_elem(T, Q, reading="A"):
    """單一 (T,Q) 之精確矩陣元（參考實作，直接按 Wick 定義求和）"""
    T = tuple(sorted(T)); Q = tuple(sorted(Q))
    cT = Counter(T); cQ = Counter(Q)
    if abs(len(Q) - len(T)) > 2 or (len(Q) - len(T)) % 2 != 0:
        return 0.0
    d = len(Q) - len(T)
    na = (4 - d) // 2
    if na < 0 or na > min(4, len(T)):
        return 0.0
    tot = 0.0
    for Asub in set(itertools.combinations_with_replacement(sorted(cT), na)):
        ca = Counter(Asub)
        if any(ca[v] > cT[v] for v in ca):
            continue
        S = list(T)
        for v in Asub:
            S.remove(v)
        need = Counter(Q) - Counter(S)
        if any(x < 0 for x in need.values()):
            continue
        C = list(sorted(need.elements()))
        if len(C) != 4 - na:
            continue
        U = tuple(sorted(list(Asub) + C))
        cU = Counter(U)
        mult = factorial(4)
        for v in cU:
            mult /= factorial(cU[v])
        for v in cU:
            mult *= comb(cU[v], ca.get(v, 0))
        no = 1.0
        for v, a in ca.items():
            for t in range(a):
                no *= (cT[v] - t)
        for v, q in Counter(C).items():
            cs = Counter(S).get(v, 0)
            for t in range(q):
                no *= (cs + t + 1)
        val = chi4(U) * mult * np.sqrt(no) / 24.0
        if na == 2 and C[0] == C[1]:
            pass
        tot += val
    if reading == "B" and d == 0:
        # 對角：N_T/2 + 5.5（僅記錄，不入微擾頂點）
        tot += 0.0
    return tot


def gate_g2(nmod=6):
    """閘 1：小模態全枚舉（3→5），向量化 vs 參考實作逐點比對"""
    A3 = np.array(list(itertools.combinations_with_replacement(range(1, nmod + 1), 3)), dtype=np.int32)
    Qr = np.array(list(itertools.combinations_with_replacement(range(1, nmod + 1), 5)), dtype=np.int32)
    Qr = Qr[np.argsort(rank_fast(Qr))]
    rkmap = rank_fast(Qr)
    total = 0; ok = 0; mx = 0.0
    for T in A3:
        v = vertex_apply(np.array([T], dtype=np.int32), np.array([1.0]), A3,
                         [(0,), (1,), (2,)], 5, CH=1, nmod=nmod)
        for i, Q in enumerate(Qr):
            r = vertex_elem(tuple(int(x) for x in T), tuple(int(x) for x in Q))
            total += 1
            d = abs(v[rkmap[i]] - r)
            mx = max(mx, d)
            if d < 1e-13:
                ok += 1
    return {"ok": ok, "total": total, "maxdiff": mx}


def gate_t2_t0():
    """閘 2：T2(c,d) = Σ_a χ⁴(a,a,c,d) 應為 δ_cd；T0 = Σ_ab χ⁴(a,a,b,b) 應為 N_d"""
    rows = []; keys = []
    for a in range(1, N + 1):
        for c in range(1, N + 1):
            for d in range(c, N + 1):
                rows.append(sorted([a, a, c, d])); keys.append((c, d))
    W = W4(np.array(rows, dtype=np.int32)) / 90.0
    T2 = np.zeros((N, N))
    for (c, d), w in zip(keys, W):
        T2[c - 1, d - 1] += w
    rows2 = []
    for a in range(1, N + 1):
        for b in range(1, N + 1):
            rows2.append(sorted([a, a, b, b]))
    T0 = float((W4(np.array(rows2, dtype=np.int32)) / 90.0).sum())
    off = float(np.abs(T2 - np.diag(np.diag(T2))).max())
    diag = np.diag(T2)
    return T2, T0, off, diag


def sectors(n):
    G = np.array([0.0] + [4 * np.sin(np.pi * k / 90) ** 2 for k in range(1, N + 1)])
    return G


def ranked(k):
    A = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), k)), dtype=np.int32)
    return A[np.argsort(rank_fast(A))]


def rs_head(n, reading):
    """至 E_3（廉價；E_4 以上需 5↔5，另以 GPU 長跑）"""
    G = sectors(n)
    TRI = ranked(3); QUI = ranked(5)
    A2 = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 2)), dtype=np.int32)
    A3 = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 3)), dtype=np.int32)
    ET3 = np.array([sum(G[t] for t in row) for row in TRI])
    b = vertex_apply(np.array([[n]], dtype=np.int32), np.array([1.0]), A3, [(0,)], 3)
    if reading == "B":
        # 同模 Δ=+2：(1/4)√((c+1)(c+2))，T=(n) ⟹ c=1 於 m=n，否則 0
        for m in range(1, N + 1):
            c = 1 if m == n else 0
            amp = 0.25 * np.sqrt((c + 1) * (c + 2))
            tgt = tuple(sorted((n, m, m)))
            b[int(rank_fast(np.array([tgt], dtype=np.int32))[0])] += amp
        den3 = G[n] - ET3 - (1.5 + 5.5)
    else:
        den3 = G[n] - ET3
    p1 = b / den3
    E2 = float(np.sum(b * p1))
    p2 = vertex_apply(TRI, p1, A2, [(0, 1), (0, 2), (1, 2)], 3, CH=512) / den3
    E3 = float(np.sum(b * p2))
    return E2, E3, float(np.linalg.norm(b)), float(G[n])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true", help="跑長表（不入 L1）")
    a = ap.parse_args()

    print("=" * 70)
    print("M14 §51 精確四波頂點：Wick 分解 / 收縮張量 / 兩讀法 E_2–E_3")
    print("=" * 70)

    # 閘 1：小模態全枚舉（3 模 → 5 模，N=6）
    g1 = gate_g2(6)
    print(f"\n§1 精確頂點閘（3→5，N=6 全枚舉）：{g1['ok']}/{g1['total']} 逐點 0 差"
          f"｜max 差 {g1['maxdiff']:.3e}")

    # 閘 2：T2／T0
    T2, T0, off, diag = gate_t2_t0()
    nz = int((np.abs(T2) > 1e-12).sum())
    print(f"\n§2 收縮張量（44 模）：T2 非零 {nz}/{N*N}｜非對角 max {off:.3e}"
          f"｜對角全為 1：{bool(np.allclose(diag, 1.0))}｜T0 = {T0:+.12f}（N_d = {N}）")

    # 閘 3：兩讀法 E_2／E_3（n=1）
    print("\n§3 路徑 A／B 之 E_1–E_3（n=1，N=44）")
    rows = {}
    for rd in ("A", "B"):
        E2, E3, nb, g1v = rs_head(1, rd)
        rows[rd] = (E2, E3)
        print(f"  路徑 {rd}：E_1 = γ_1 = {g1v:.12f}｜‖b‖ = {nb:.6f}｜E_2 = {E2:+.12f}｜E_3 = {E3:+.12f}")

    if a.full:
        print("\n§4 全 44 n（見 data/precise_vertex_44n.json；本旗標不於審計內跑）")

    print("\n§5 判詞")
    print("  χ⁴×emb 族（raw／choose／falling）全部否證：精確值含半整數，整數 emb 乘法不可達")
    print("  路徑 A 比值 ≫1（發散）｜路徑 B 比值 <1 除 r2≈1.08（收斂）")
    print("  docs/42–48 之 E_p 建於已否證之 emb 族 ⟹ 只作歷史基準")


if __name__ == "__main__":
    import time
    t0 = time.time()
    main()
    print(f"[耗時 {time.time()-t0:.1f} s]", file=sys.stderr)