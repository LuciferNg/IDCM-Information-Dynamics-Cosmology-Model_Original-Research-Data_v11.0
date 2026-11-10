#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GPU 四波算子（CuPy RawModule）：三通道 ＋ 同模通道
公式與 CPU e7chain.op 逐項相同：
  amp = (W/90)·√(noA·noC)·φ / (fpA·fpC·mulpos)
閘：對 CPU 版逐點比對（見 gate_gpu.py）
"""
import numpy as np
import cupy as cp
from math import comb

N = 44
NB = np.zeros((64, 12), dtype=np.int64)
for _n in range(64):
    for _r in range(12):
        NB[_n, _r] = comb(_n, _r) if 0 <= _r <= _n else 0

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

// ---- na=2, ncre=2（3→3、5→5、7→7）----
extern "C" __global__ void k22(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt,
    const int i1, const int i2, const int x, const int y, const int addidx)
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

// ---- na=1, ncre=3（3→5、5→7、1→3）----
extern "C" __global__ void k13(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt, const int i1,
    const int x, const int y, const int z, const int addrank, const double fpC, const int nsame)
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

// ---- na=3, ncre=1（5→3、7→5）----
extern "C" __global__ void k31(const int* __restrict__ S, const float* __restrict__ PHI, double* __restrict__ OUT,
    const long long* __restrict__ b, const signed char* __restrict__ WT,
    const int m, const int ksrc, const int kt,
    const int i1, const int i2, const int i3, const int x)
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
    int rem = 0;
    // 移除三個位置（以索引為準）
    for (int j = 0; j < ksrc; ++j) if (j != i1 && j != i2 && j != i3) Q[n++] = T[j];
    int cSx = 0;
    for (int j = 0; j < n; ++j) if (Q[j] == x) ++cSx;
    Q[n++] = x; srt(Q, kt);
    double amp = ((double)W / 90.0) * sqrt(noA * (double)(cSx + 1)) * (double)PHI[g] / (fpA * 1.0 * mulpos);
    atomicAdd(&OUT[rk(b, Q, kt)], amp);
}

// ---- 同模 Δ=±2 ----
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

_CUDA_ASCII = __import__("re").sub(r"[^\x00-\x7F]", " ", CUDA)
_mod = cp.RawModule(code=_CUDA_ASCII)
K22 = _mod.get_function("k22")
K13 = _mod.get_function("k13")
K31 = _mod.get_function("k31")
KSM = _mod.get_function("ksm")
NB_G = cp.asarray(NB)

_WT_CACHE = {}


def _tables(WT_np, key):
    if key not in _WT_CACHE:
        _WT_CACHE[key] = cp.asarray(WT_np.astype(np.int8))
    return _WT_CACHE[key]


def gpu_op(src, phi, na, ncre, kt, out_dim, WT_np, CH=1_500_000, threads=256):
    """src (M,ksrc) int32 numpy；回傳 float64 numpy（out_dim）"""
    import itertools
    out = cp.zeros(out_dim, dtype=cp.float64)
    WTg = _tables(WT_np, id(WT_np))
    ksrc = src.shape[1]
    M = src.shape[0]
    if na == 0 or ncre == 0:
        for s0 in range(0, M, CH):
            s1 = min(s0 + CH, M)
            Sg = cp.asarray(src[s0:s1]); Pg = cp.asarray(np.asarray(phi[s0:s1]).astype(np.float32))
            for mm in range(1, N + 1):
                KSM(((s1 - s0 + threads - 1) // threads,), (threads,),
                    (Sg, Pg, out, NB_G, np.int32(s1 - s0), np.int32(ksrc), np.int32(kt),
                     np.int32(mm), np.int32(1 if ncre == 2 else 0)))
            del Sg, Pg
            cp.get_default_memory_pool().free_all_blocks()
        r = cp.asnumpy(out); del out; cp.get_default_memory_pool().free_all_blocks()
        return r

    rems = list(itertools.combinations(range(ksrc), na))
    ADD = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), ncre)), dtype=np.int32)
    # 附加表次序
    if ncre == 2:
        ADD = np.array([(x, y) for x in range(1, N + 1) for y in range(x, N + 1)], dtype=np.int32)
    elif ncre == 3:
        ADD = np.array(list(itertools.combinations_with_replacement(range(1, N + 1), 3)), dtype=np.int32)
    for rm in rems:
        for s0 in range(0, M, CH):
            s1 = min(s0 + CH, M)
            Sg = cp.asarray(src[s0:s1]); Pg = cp.asarray(np.asarray(phi[s0:s1]).astype(np.float32))
            mm_ = s1 - s0
            grid = ((mm_ + threads - 1) // threads,)
            if na == 2 and ncre == 2:
                for ai, (x, y) in enumerate(ADD.tolist()):
                    K22(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc),
                                           np.int32(kt), np.int32(rm[0]), np.int32(rm[1]),
                                           np.int32(x), np.int32(y), np.int32(ai)))
            elif na == 1 and ncre == 3:
                for ai, t3 in enumerate(ADD.tolist()):
                    _rk3 = int(sum(NB[min(max(int(t3[j]) + j - 1, 0), 63), j + 1] for j in range(3)))
                    x, y, z = int(t3[0]), int(t3[1]), int(t3[2])
                    nsame = 2 if (x == y == z) else (1 if (x == y or y == z or x == z) else 0)
                    fpC = 6.0 if nsame == 2 else (2.0 if nsame == 1 else 1.0)
                    K13(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc),
                                           np.int32(kt), np.int32(rm[0]), np.int32(x), np.int32(y),
                                           np.int32(z), np.int32(_rk3), np.float64(fpC), np.int32(nsame)))
            elif na == 3 and ncre == 1:
                for x in range(1, N + 1):
                    K31(grid, (threads,), (Sg, Pg, out, NB_G, WTg, np.int32(mm_), np.int32(ksrc),
                                           np.int32(kt), np.int32(rm[0]), np.int32(rm[1]), np.int32(rm[2]),
                                           np.int32(x)))
            del Sg, Pg
            cp.get_default_memory_pool().free_all_blocks()
    r = cp.asnumpy(out); del out; cp.get_default_memory_pool().free_all_blocks()
    return r