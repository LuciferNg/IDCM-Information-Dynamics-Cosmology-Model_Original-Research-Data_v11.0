#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
53 — GPU 四波算子（CuPy RawModule）＋ 對 `scripts/52` CPU 算子之逐點閘
================================================================================
四個 kernel：k22（na=2,cre=2）、k13（na=1,cre=3）、k31（na=3,cre=1）、ksm（同模 Δ=±2）
公式與 CPU 版逐項相同：amp = (W/90)·√(noA·noC)·φ / (fpA·fpC·mulpos)
φ 以 float32 上傳（頻寬／記憶體權衡），累加用 float64 ⟹ 預期相對差 ~1e−8。

閘：五通道（3→3、3→5、5→3、同模創建、同模湮滅）對 CPU 之最大相對差 < 1e−6。
"""
import importlib.util as _iu
import itertools
import re
import sys
from math import comb
from pathlib import Path

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

import numpy as np
import cupy as cp

HERE = Path(__file__).resolve().parent
_spec = _iu.spec_from_file_location("fw52", HERE / "52_fourwave_operator.py")
FW = _iu.module_from_spec(_spec)
_spec.loader.exec_module(FW)

N = 44
NB = FW.NB
CUDA = r"""
__device__ inline long long rk(const long long* __restrict__ b, const int* T, const int k) {
    long long r = 0;
    for (int j = 0; j < k; ++j) {
        int nn = T[j] + j - 1;
        if (nn < 0) nn = 0;
        if (nn > 63) nn = 63;
        r += b[nn * 12 + (j + 1)];
    }
    return r;
}
__device__ inline void srt(int* T, const int k) {
    for (int i = 1; i < k; ++i) {
        int v = T[i], j = i - 1;
        while (j >= 0 && T[j] > v) { T[j + 1] = T[j]; --j; }
        T[j + 1] = v;
    }
}
extern "C" __global__ void k22(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt, const int i1, const int i2,
    const int x, const int y, const int addidx)
{
    long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= m) return;
    const int* T = S + g * ksrc;
    int va = T[i1], vb = T[i2];
    int cTa = 0, cTb = 0;
    for (int j = 0; j < ksrc; ++j) { if (T[j] == va) ++cTa; if (T[j] == vb) ++cTb; }
    int W = (int)WT[addidx * 1936 + (va - 1) * 44 + (vb - 1)];
    if (W == 0) return;
    int Q[9]; int n = 0;
    for (int j = 0; j < ksrc; ++j) { if (j != i1 && j != i2) Q[n++] = T[j]; }
    int cSx = 0, cSy = 0;
    for (int j = 0; j < n; ++j) { if (Q[j] == x) ++cSx; if (Q[j] == y) ++cSy; }
    Q[n++] = x; Q[n++] = y; srt(Q, kt);
    int same = (va == vb);
    double fpA = same ? 2.0 : 1.0;
    double mulpos = same ? (double)cTa * (cTa - 1) / 2.0 : (double)cTa * cTb;
    double noA = same ? (double)cTa * (cTa - 1) : (double)cTa * cTb;
    int sc = (x == y);
    double fpC = sc ? 2.0 : 1.0;
    double noC = sc ? (double)(cSx + 1) * (cSx + 2) : (double)(cSx + 1) * (cSy + 1);
    double amp = ((double)W / 90.0) * sqrt(noA * noC) * (double)PHI[g] / (fpA * fpC * mulpos);
    atomicAdd(&OUT[rk(b, Q, kt)], amp);
}
extern "C" __global__ void k13(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt, const int i1,
    const int x, const int y, const int z, const int addrank)
{
    long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= m) return;
    const int* T = S + g * ksrc;
    int va = T[i1];
    int cTa = 0;
    for (int j = 0; j < ksrc; ++j) if (T[j] == va) ++cTa;
    int W = (int)WT[(va - 1) * 15180 + addrank];
    if (W == 0) return;
    int Q[9]; int n = 0;
    for (int j = 0; j < ksrc; ++j) if (j != i1) Q[n++] = T[j];
    int V3[3] = {x, y, z}; srt(V3, 3);
    int nsurv = n;
    double noC = 1.0, fpCl = 1.0;
    int j0 = 0;
    while (j0 < 3) {
        int q = 1;
        while (j0 + q < 3 && V3[j0 + q] == V3[j0]) ++q;
        int cS = 0;
        for (int j = 0; j < nsurv; ++j) if (Q[j] == V3[j0]) ++cS;
        for (int u = 1; u <= q; ++u) { noC *= (double)(cS + u); fpCl *= (double)u; }
        j0 += q;
    }
    Q[n++] = x; Q[n++] = y; Q[n++] = z; srt(Q, kt);
    double amp = ((double)W / 90.0) * sqrt((double)cTa * noC) * (double)PHI[g] / ((double)cTa * fpCl);
    atomicAdd(&OUT[rk(b, Q, kt)], amp);
}
extern "C" __global__ void k31(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt, const int i1, const int i2, const int i3, const int x)
{
    long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= m) return;
    const int* T = S + g * ksrc;
    int A[3] = {T[i1], T[i2], T[i3]};
    srt(A, 3);
    int cT[3];
    for (int j = 0; j < 3; ++j) { int c = 0; for (int q = 0; q < ksrc; ++q) if (T[q] == A[j]) ++c; cT[j] = c; }
    long long rkA = 0;
    for (int j = 0; j < 3; ++j) { int nn = A[j] + j - 1; if (nn < 0) nn = 0; if (nn > 63) nn = 63; rkA += b[nn * 12 + (j + 1)]; }
    int W = (int)WT[(x - 1) * 15180 + (int)rkA];
    if (W == 0) return;
    double fpA = 1.0, mulpos = 1.0, noA = 1.0;
    for (int j = 0; j < 3; ++j) {
        int cb = 0;
        for (int q = 0; q < j; ++q) if (A[q] == A[j]) cb++;
        double c = (double)cT[j];
        noA *= (c - cb);
        mulpos *= (c - cb) / (cb + 1.0);
    }
    if (A[0] == A[2]) fpA = 6.0;
    else if (A[0] == A[1] || A[1] == A[2]) fpA = 2.0;
    int Q[9]; int n = 0;
    for (int j = 0; j < ksrc; ++j) if (j != i1 && j != i2 && j != i3) Q[n++] = T[j];
    int cSx = 0;
    for (int j = 0; j < n; ++j) if (Q[j] == x) ++cSx;
    Q[n++] = x; srt(Q, kt);
    double amp = ((double)W / 90.0) * sqrt(noA * (double)(cSx + 1)) * (double)PHI[g] / (fpA * mulpos);
    atomicAdd(&OUT[rk(b, Q, kt)], amp);
}
extern "C" __global__ void ksm(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const int m, const int ksrc, const int kt, const int mm, const int create)
{
    long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= m) return;
    const int* T = S + g * ksrc;
    int c = 0;
    for (int j = 0; j < ksrc; ++j) if (T[j] == mm) ++c;
    if (!create && c < 2) return;
    double a = create ? 0.25 * sqrt((double)(c + 1) * (c + 2)) : 0.25 * sqrt((double)c * (c - 1));
    int Q[9]; int n = 0;
    if (create) {
        for (int j = 0; j < ksrc; ++j) Q[n++] = T[j];
        Q[n++] = mm; Q[n++] = mm;
    } else {
        int rem = 2;
        for (int j = 0; j < ksrc; ++j) { if (T[j] == mm && rem > 0) { --rem; continue; } Q[n++] = T[j]; }
    }
    srt(Q, kt);
    atomicAdd(&OUT[rk(b, Q, kt)], a * (double)PHI[g]);
}
"""
_mod = cp.RawModule(code=re.sub(r"[^\x00-\x7F]", " ", CUDA))
K22, K13, K31, KSM = (_mod.get_function(n) for n in ("k22", "k13", "k31", "ksm"))
NB_G = cp.asarray(NB)


def gpu_op(src, phi, na, ncre, kt, out_dim, WT, CH=1_500_000, threads=256):
    out = cp.zeros(out_dim, dtype=cp.float64)
    WTg = cp.asarray(WT.astype(np.int8))
    ksrc, M = src.shape[1], src.shape[0]
    if na == 0:
        for s0 in range(0, M, CH):
            s1 = min(s0 + CH, M)
            Sg = cp.asarray(src[s0:s1]); Pg = cp.asarray(phi[s0:s1].astype(np.float32))
            for mm in range(1, N + 1):
                KSM(((s1 - s0 + threads - 1) // threads,), (threads,),
                    (Sg, Pg, out, NB_G, np.int32(s1 - s0), np.int32(ksrc), np.int32(kt),
                     np.int32(mm), np.int32(1 if ncre == 2 else 0)))
            cp.get_default_memory_pool().free_all_blocks()
        r = cp.asnumpy(out); del out; cp.get_default_memory_pool().free_all_blocks(); return r
    ADD = FW.MULTI[ncre]
    for rm in itertools.combinations(range(ksrc), na):
        for s0 in range(0, M, CH):
            s1 = min(s0 + CH, M)
            Sg = cp.asarray(src[s0:s1]); Pg = cp.asarray(phi[s0:s1].astype(np.float32))
            mm_ = s1 - s0
            grid = ((mm_ + threads - 1) // threads,)
            if na == 2 and ncre == 2:
                for ai, (x, y) in enumerate(ADD.tolist()):
                    K22(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc), np.int32(kt),
                                           np.int32(rm[0]), np.int32(rm[1]), np.int32(x), np.int32(y), np.int32(ai)))
            elif na == 1 and ncre == 3:
                for ai, t3 in enumerate(ADD.tolist()):
                    rk3 = int(sum(NB[min(max(int(t3[j]) + j - 1, 0), 63), j + 1] for j in range(3)))
                    K13(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc), np.int32(kt),
                                           np.int32(rm[0]), np.int32(t3[0]), np.int32(t3[1]), np.int32(t3[2]),
                                           np.int32(rk3)))
            elif na == 3 and ncre == 1:
                for x in range(1, N + 1):
                    K31(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc), np.int32(kt),
                                           np.int32(rm[0]), np.int32(rm[1]), np.int32(rm[2]), np.int32(x)))
            cp.get_default_memory_pool().free_all_blocks()
    r = cp.asnumpy(out); del out; cp.get_default_memory_pool().free_all_blocks(); return r


def main():
    rng = np.random.default_rng(20261109)
    nsub = 20000
    TRI, QUI = FW.COLEX[3], FW.COLEX[5]
    fail = 0
    print("=" * 78)
    print("G-A3 GPU 對 CPU 逐點閘（nsub = 20,000，容差 1e−6）")
    print("=" * 78)
    cases = [("3→3 (2a2c)", TRI[:nsub], 2, 2, 3, len(TRI), FW.WT22),
             ("3→5 (1a3c)", TRI[:nsub], 1, 3, 5, len(QUI), FW.WT13),
             ("5→3 (3a1c)", QUI[:nsub], 3, 1, 3, len(TRI), FW.WT31)]
    for name, src, na, ncre, kt, od, WT in cases:
        phi = rng.standard_normal(len(src))
        a = FW.op(src, phi, na, ncre, kt, od)
        b = gpu_op(src, phi, na, ncre, kt, od, WT)
        d = float(np.abs(a - b).max() / max(1e-300, np.abs(a).max()))
        ok = d < 1e-6
        fail += 0 if ok else 1
        print(f"  {name}: 相對最大差 {d:.3e}  {'PASS' if ok else 'FAIL'}")
    for name, k, kt, create, od in (("同模創建 3→5", 3, 5, True, len(QUI)),
                                    ("同模湮滅 5→3", 5, 3, False, len(TRI))):
        src = FW.COLEX[k][:nsub]
        phi = rng.standard_normal(len(src))
        a = FW.sm_op(src, phi, kt, create, od)
        b = gpu_op(src, phi, 0, 2 if create else -2, kt, od, FW.WT22)
        d = float(np.abs(a - b).max() / max(1e-300, np.abs(a).max()))
        ok = d < 1e-6
        fail += 0 if ok else 1
        print(f"  {name}: 相對最大差 {d:.3e}  {'PASS' if ok else 'FAIL'}")
    print("=" * 78)
    print(f"檢查項 = 5，FAIL = {fail}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())