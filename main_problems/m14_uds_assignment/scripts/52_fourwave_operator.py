#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
52 — 通用四波算子（CPU，float64）＋ 六通道閘
================================================================================
背景：docs/52 之 E_7 全 CPU 鏈所用算子；本腳本把算子固化並附**快速閘**，供 L1 重跑。

算子（amp(Q|T) 對單一 (removal, add) 組合）：
    amp = (W/90)·√(noA·noC)·φ(T) / (fpA·fpC·mulpos)
      A = 被湮滅多重集（na 個）      C = 被創造多重集（ncre 個）
      c_T(v) = T 中 v 之重數          c_S(v) = 倖存 S 中 v 之重數
      noA  = ∏_{v∈A}(c_T − t_v)      noC = ∏_{v∈C}∏_{u=1}^{q_v}(c_S + u)
      fpA  = ∏ c_A!   fpC = ∏ c_C!    mulpos = ∏_{v∈A distinct} C(c_T(v), t_v)
      W    = χ⁴(A ∪ C)（八個帶號和 ≡ 0 mod 90 之加權計數）

閘：
  G1 三通道（3→3、3→5、5→3）對 scripts/51 `vertex_apply` 逐點
  G2 同模 Δ=±2 通道之伴隨性 ⟨φ_Q|sm_cre|φ_T⟩ = ⟨φ_T|sm_ann|φ_Q⟩
  G3 χ⁴ 表非零率 ≈ 8.5%
本腳本只讀 scripts/51 之函數；不寫任何 docs。
"""
import itertools
import sys
from math import comb
from pathlib import Path

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util as _iu

_spec = _iu.spec_from_file_location("pv51", HERE / "51_precise_vertex.py")
PV = _iu.module_from_spec(_spec)
_spec.loader.exec_module(PV)

N = 44
NB = np.zeros((64, 12), dtype=np.int64)
for _n in range(64):
    for _r in range(12):
        NB[_n, _r] = comb(_n, _r) if 0 <= _r <= _n else 0
G = np.array([0.0] + [4 * np.sin(np.pi * k / 90) ** 2 for k in range(1, N + 1)])
SIG = np.array([(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
                (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)], dtype=np.int32)
CS = np.array([1, 1, -1, -1, -1, -1, 1, 1], dtype=np.int32)

MULTI = {k: np.array(list(itertools.combinations_with_replacement(range(1, N + 1), k)), dtype=np.int32)
         for k in (1, 2, 3, 5)}


def _rank(M):
    r = np.zeros(M.shape[0], dtype=np.int64)
    for j in range(M.shape[1]):
        r += NB[np.clip(M[:, j].astype(np.int64) + j - 1, 0, 63), j + 1]
    return r


RANK = {k: _rank(MULTI[k]) for k in MULTI}
COLEX = {k: MULTI[k][np.argsort(RANK[k])] for k in MULTI}


def chi4(rows):
    U = np.sort(rows, axis=1).astype(np.int32)
    W = np.zeros(U.shape[0], dtype=np.int32)
    for c, sg in zip(CS, SIG):
        W += np.where(np.abs(U @ sg) % 90 == 0, c, 0)
    return W


def build_WT13():
    tri = MULTI[3]
    WT = np.zeros((N, len(tri)), dtype=np.int8)
    rows = np.zeros((len(tri), 4), dtype=np.int32)
    rows[:, 1:] = tri
    for a in range(1, N + 1):
        rows[:, 0] = a
        WT[a - 1, RANK[3]] = chi4(rows).astype(np.int8)
    return WT


def build_WT31():
    tri = MULTI[3]
    WT = np.zeros((N, len(tri)), dtype=np.int8)
    rows = np.zeros((len(tri), 4), dtype=np.int32)
    rows[:, :3] = tri
    for x in range(1, N + 1):
        rows[:, 3] = x
        WT[x - 1, RANK[3]] = chi4(rows).astype(np.int8)
    return WT


def build_WT22():
    A2 = MULTI[2]
    WT = np.zeros((len(A2), N, N), dtype=np.int8)
    grid = np.arange(1, N + 1)
    VA, VB = np.meshgrid(grid, grid, indexing="ij")
    for t, (x, y) in enumerate(A2.tolist()):
        rows = np.stack([VA.ravel(), VB.ravel(), np.full(VA.size, x, np.int32),
                         np.full(VA.size, y, np.int32)], axis=1)
        WT[t] = chi4(rows).reshape(N, N).astype(np.int8)
    return WT


WT13, WT31, WT22 = build_WT13(), build_WT31(), build_WT22()


def op(src, phi, na, ncre, kt, out_dim, CH=1_500_000):
    ksrc = src.shape[1]
    rems = list(itertools.combinations(range(ksrc), na))
    add = MULTI[ncre]
    addrank = RANK[ncre]
    out = np.zeros(out_dim, dtype=np.float64)
    for rm in rems:
        keep = [j for j in range(ksrc) if j not in rm]
        for s0 in range(0, src.shape[0], CH):
            s1 = min(s0 + CH, src.shape[0])
            S = src[s0:s1]; ph = phi[s0:s1]
            A = S[:, list(rm)]
            A_sort = np.sort(A, axis=1)
            Ssur = S[:, keep]
            m = S.shape[0]
            if na == 1:
                fpA = np.ones(m)
            elif na == 2:
                fpA = np.where(A_sort[:, 0] == A_sort[:, 1], 2.0, 1.0)
            else:
                e01 = A_sort[:, 0] == A_sort[:, 1]; e12 = A_sort[:, 1] == A_sort[:, 2]
                al = e01 & e12
                fpA = np.where(al, 6.0, np.where((e01 | e12) & ~al, 2.0, 1.0))
            noA = np.ones(m); mulpos = np.ones(m)
            for j in range(na):
                v = A_sort[:, j]
                cT = (S == v[:, None]).sum(1).astype(np.float64)
                cb = np.zeros(m)
                for jj in range(j):
                    cb += (A_sort[:, jj] == v)
                noA *= (cT - cb)
                mulpos *= (cT - cb) / (cb + 1.0)
            rkA = _rank(A_sort) if na == 3 else None
            for ai in range(len(add)):
                if na == 1 and ncre == 3:
                    W = WT13[A_sort[:, 0].astype(np.int64) - 1, addrank[ai]]
                elif na == 3 and ncre == 1:
                    W = WT31[ai, rkA]
                elif na == 2 and ncre == 2:
                    W = WT22[ai][A_sort[:, 0].astype(np.int64) - 1, A_sort[:, 1].astype(np.int64) - 1]
                else:
                    raise ValueError("組合未支援")
                nz = W != 0
                k = int(nz.sum())
                if k == 0:
                    continue
                Sh = Ssur[nz]
                addv = np.sort(add[ai])
                noC = np.ones(k); fpC = np.ones(k)
                for j in range(ncre):
                    v = addv[j]
                    cS = (Sh == v).sum(1).astype(np.float64)
                    cb = 0.0
                    for jj in range(j):
                        cb += float(addv[jj] == v)
                    noC *= (cS + cb + 1.0)
                    fpC *= (cb + 1.0)
                T = np.sort(np.concatenate([Sh, np.broadcast_to(addv, (k, ncre))], axis=1), axis=1)
                amp = (W[nz] / 90.0) * np.sqrt(noA[nz] * noC) * ph[nz] / (fpA[nz] * fpC * mulpos[nz])
                np.add.at(out, _rank(T), amp)
    return out


def sm_op(src, phi, kt, create, out_dim):
    """同模 Δ=±2：create ⟹ (1/4)√((c+1)(c+2))φ；annihilate ⟹ (1/4)√(c(c−1))φ"""
    ksrc = src.shape[1]
    out = np.zeros(out_dim, dtype=np.float64)
    for m in range(1, N + 1):
        c = (src == m).sum(1).astype(np.float64)
        if create:
            amp = 0.25 * np.sqrt((c + 1.0) * (c + 2.0)) * phi
            T = np.sort(np.concatenate([src, np.full((src.shape[0], 2), m, np.int32)], axis=1), axis=1)
        else:
            mask = c >= 2
            if not mask.any():
                continue
            amp = 0.25 * np.sqrt(c[mask] * (c[mask] - 1.0)) * phi[mask]
            key = np.where(src[mask] == m, 99, src[mask])
            order = np.argsort(key, axis=1, kind="stable")
            T = np.sort(np.take_along_axis(src[mask], order, axis=1)[:, :ksrc - 2], axis=1)
        np.add.at(out, _rank(T), amp)
    return out


def main():
    TRI, QUI = COLEX[3], COLEX[5]
    rng = np.random.default_rng(20261109)
    nsub = 6000
    fail = 0

    print("=" * 78)
    print("G1 三通道 vs scripts/51 vertex_apply")
    print("=" * 78)
    for name, src, na, ncre, kt, od, rems, addt in (
            ("3→3 (2a2c)", TRI[:nsub], 2, 2, 3, len(TRI), [(0, 1), (0, 2), (1, 2)], MULTI[2]),
            ("3→5 (1a3c)", TRI[:nsub], 1, 3, 5, len(QUI), [(0,), (1,), (2,)], MULTI[3]),
            ("5→3 (3a1c)", QUI[:nsub], 3, 1, 3, len(TRI),
             [(a, b, c) for a in range(5) for b in range(a + 1, 5) for c in range(b + 1, 5)], MULTI[1])):
        phi = rng.standard_normal(len(src))
        got = op(src, phi, na, ncre, kt, od)
        ref = PV.vertex_apply(src, phi, addt, rems, kt, CH=2000)
        d = np.abs(got - ref).max() / max(1e-300, np.abs(ref).max())
        ok = d < 1e-12
        fail += 0 if ok else 1
        print(f"  {name}: 相對最大差 {d:.3e}  {'PASS' if ok else 'FAIL'}")

    print("=" * 78)
    print("G2 同模通道伴隨性")
    print("=" * 78)
    for k, kt, create in ((3, 5, True), (5, 3, False)):
        src = COLEX[k][:nsub]
        phi = rng.standard_normal(len(src))
        a = sm_op(src, phi, kt, create, len(MULTI[kt]))
        Q = COLEX[kt][:3 * nsub]
        phiq = rng.standard_normal(len(Q))
        b = sm_op(Q, phiq, k, not create, len(MULTI[k]))
        lhs = float(np.dot(phiq, a[:len(Q)]))
        rhs = float(np.dot(phi, b[:len(src)]))
        d = abs(lhs - rhs) / max(1e-300, abs(lhs))
        ok = d < 1e-10
        fail += 0 if ok else 1
        print(f"  {k}→{kt} create={create}: 伴隨相對差 {d:.3e}  {'PASS' if ok else 'FAIL'}")

    print("=" * 78)
    print("G3 χ⁴ 表統計")
    print("=" * 78)
    nz = (WT13 != 0).mean()
    print(f"  WT13 非零率 {nz:.6f}（期望 ≈ 0.085）{'PASS' if 0.07 < nz < 0.10 else 'FAIL'}")
    print("=" * 78)
    print(f"檢查項 = 6，FAIL = {fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())