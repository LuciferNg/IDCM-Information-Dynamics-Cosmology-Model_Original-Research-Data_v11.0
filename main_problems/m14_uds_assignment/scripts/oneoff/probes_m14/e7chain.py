#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_7 全 CPU 鏈（重建）：通用四波算子 + 黃金標準閘 + checkpoint
amp = (W/90) · √(noA·noC) · φ / (fpA·fpC·mulpos)
"""
import os, sys, json, time, itertools
from math import comb
import numpy as np

N = 44
CK = "/home/wsl/.hermes/cache/scratch/e7ckpt2"
os.makedirs(CK, exist_ok=True)
REPO = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/main_problems/m14_uds_assignment"

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
NRK = {k: len(MULTI[k]) for k in MULTI}


def _rank(M):
    r = np.zeros(M.shape[0], dtype=np.int64)
    k = M.shape[1]
    for j in range(k):
        r += NB[np.clip(M[:, j].astype(np.int64) + j - 1, 0, 63), j + 1]
    return r


RANK = {k: _rank(MULTI[k]) for k in MULTI}
# ★ 算子輸出為 colex（combinadic rank）序；源表與能量表必須同序
COLEX = {k: MULTI[k][np.argsort(RANK[k])] for k in MULTI}
COLEXRANK = {k: np.sort(RANK[k]) for k in MULTI}


def chi4(rows):
    U = np.sort(rows, axis=1).astype(np.int32)
    W = np.zeros(U.shape[0], dtype=np.int32)
    for c, sg in zip(CS, SIG):
        W += np.where(np.abs(U @ sg) % 90 == 0, c, 0)
    return W


def build_WT13():
    tri = MULTI[3]
    WT = np.zeros((N, NRK[3]), dtype=np.int8)
    rows = np.zeros((len(tri), 4), dtype=np.int32); rows[:, 1:] = tri
    rk = RANK[3]
    for a in range(1, N + 1):
        rows[:, 0] = a
        WT[a - 1, rk] = chi4(rows).astype(np.int8)
    return WT


def build_WT31():
    tri = MULTI[3]
    WT = np.zeros((N, NRK[3]), dtype=np.int8)
    rows = np.zeros((len(tri), 4), dtype=np.int32); rows[:, :3] = tri
    rk = RANK[3]
    for x in range(1, N + 1):
        rows[:, 3] = x
        WT[x - 1, rk] = chi4(rows).astype(np.int8)
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


_t0 = time.time()
WT13 = build_WT13(); WT31 = build_WT31(); WT22 = build_WT22()


def op(src, phi, na, ncre, kt, out_dim, CH=1_500_000):
    """src (M,ksrc)；回傳 dense out_dim"""
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
            # fpA
            if na == 1:
                fpA = np.ones(m)
            elif na == 2:
                fpA = np.where(A_sort[:, 0] == A_sort[:, 1], 2.0, 1.0)
            else:
                e01 = A_sort[:, 0] == A_sort[:, 1]; e12 = A_sort[:, 1] == A_sort[:, 2]
                al = e01 & e12
                fpA = np.where(al, 6.0, np.where((e01 | e12) & ~al, 2.0, 1.0))
            # noA、mulpos
            noA = np.ones(m); mulpos = np.ones(m)
            for j in range(na):
                v = A_sort[:, j]
                cT = (S == v[:, None]).sum(1).astype(np.float64)
                cb = np.zeros(m)
                for jj in range(j):
                    cb += (A_sort[:, jj] == v)
                noA *= (cT - cb)
                mulpos *= (cT - cb) / (cb + 1.0)
            if na == 3:
                rkt = _rank(A_sort)
            for ai in range(len(add)):
                # W 查表（每次只取一欄）
                if na == 1 and ncre == 3:
                    W = WT13[A_sort[:, 0].astype(np.int64) - 1, addrank[ai]]
                elif na == 3 and ncre == 1:
                    W = WT31[ai, rkt]
                elif na == 2 and ncre == 2:
                    W = WT22[ai][A_sort[:, 0].astype(np.int64) - 1, A_sort[:, 1].astype(np.int64) - 1]
                elif na == 1 and ncre == 1:
                    raise ValueError("na=1,ncre=1 未用")
                nz = W != 0
                k = int(nz.sum())
                if k == 0:
                    continue
                Sh = Ssur[nz]
                addv = add[ai]
                # noC、fpC
                noC = np.ones(k); fpC = np.ones(k)
                as_ = np.sort(addv)
                for j in range(ncre):
                    v = as_[j]
                    cS = (Sh == v).sum(1).astype(np.float64)
                    cb = 0.0
                    for jj in range(j):
                        cb += float(as_[jj] == v)
                    noC *= (cS + cb + 1.0)
                    fpC *= (cb + 1.0)
                T = np.sort(np.concatenate([Sh, np.broadcast_to(addv, (k, ncre))], axis=1), axis=1)
                amp = (W[nz] / 90.0) * np.sqrt(noA[nz] * noC) * ph[nz] / (fpA[nz] * fpC * mulpos[nz])
                r = _rank(T)
                np.add.at(out, r, amp)
    return out


def opener(tab_file, src_slice=None):
    tab = np.load(f"{CK}/{tab_file}", mmap_mode="r")
    return np.array(tab[src_slice]) if src_slice is not None else np.array(tab)


if __name__ == "__main__":
    print("模組就緒：", time.time() - _t0, "s", flush=True)

def sm_op(src, phi, kt, create, out_dim):
    """同模 Δ=±2 通道（(1/4)Σ_c :A_c²: 之創建／湮滅部分）
    create: amp = (1/4)·√((c+1)(c+2))·φ；annihilate: amp = (1/4)·√(c(c-1))·φ
    """
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
