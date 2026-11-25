#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E_8 第①步（並行）：raw = V^{7→9} φ_3^7，按「目標最大模態 M」分片

分片原理：9-多重集之 combinadic rank 對最大元素單調 ⟹ 最大元素 = M 之態佔
          rank ∈ [C(M+7,9), C(M+8,9)) ⟹ 每 worker 擁有自己 rank 區間，無寫衝突。
跳過成本：maxT = max(kept 最大元素, add 最大元素) ⟹ 兩次比較剔除 ~11/12 候選。

用法：python3 e8_79p.py <worker> <nworker>
產物：e8ckpt/p79_w<k>.npy ＋ p79_w<k>.json
"""
import os, sys, json, time
import numpy as np
from math import comb

sys.path.insert(0, "/home/wsl/.hermes/cache/scratch")
import e7chain as C

N = 44
NQ7 = comb(50, 7)
NQ9 = comb(52, 9)
CK = "/home/wsl/.hermes/cache/scratch/e8ckpt"
os.makedirs(CK, exist_ok=True)
TABLE7 = "/home/wsl/.hermes/cache/scratch/e7ckptB/table7.npy"
PHI37 = "/home/wsl/.hermes/cache/scratch/e7ckpt2/phi3_7.npy"

EDGE = [comb(M + 7, 9) for M in range(1, N + 1)] + [NQ9]      # EDGE[M-1] = 最大元素 ≤ M 之態數
BLOCK = [EDGE[M] - EDGE[M - 1] for M in range(1, N + 1)]

WT13 = C.build_WT13()
PRUNE = {a: np.nonzero(WT13[a - 1])[0] for a in range(1, N + 1)}
MULTI3 = C.MULTI[3]
ADDMAX = MULTI3.max(1).astype(np.int32)


def acc(S, ph, MLO, MHI, out, RLO):
    """把 V^{7→9} 作用於 (S 7-多重集, ph) 之貢獻，只寫 max∈[MLO,MHI] 之 slice"""
    # (0,2) 同模通道：先篩後排（只對落在本片之態做 sort/rank）
    Smax = S.max(1)
    smin = int(Smax.min()); smax = int(Smax.max())
    for v in range(1, MHI + 1):
        if v < MLO and smax < MLO:
            continue
        if smin > MHI:
            break
        mxt = np.maximum(Smax, v)
        keep = (mxt >= MLO) & (mxt <= MHI)
        if not keep.any():
            continue
        Sk = S[keep]; phk = ph[keep]
        cS = (Sk == v).sum(1).astype(np.float64)
        T2 = np.sort(np.concatenate([Sk, np.full((Sk.shape[0], 2), v, dtype=np.int32)], axis=1), axis=1)
        r = C._rank(T2) - RLO
        np.add.at(out, r, 0.25 * np.sqrt((cS + 1.0) * (cS + 2.0)) * phk)
    # (1,3) 四波通道
    for rm in range(7):
        keepc = [j for j in range(7) if j != rm]
        a = S[:, rm]
        Ssur = S[:, keepc]
        cT = (S == a[:, None]).sum(1).astype(np.float64)
        noA = cT.copy(); mp = cT.copy()
        kpm = Ssur.max(1)
        for av in range(1, N + 1):
            sel = a == av
            if not sel.any():
                continue
            Ss = Ssur[sel]; phs = ph[sel]; noa = noA[sel]; mps = mp[sel]; kpms = kpm[sel]
            kmin = int(kpms.min()); kmax = int(kpms.max())
            if kmin > MHI:                       # 全部 kept 已超上界
                continue
            pr = PRUNE[av]
            for ai in pr:
                am = int(ADDMAX[ai])
                if am > MHI:                     # 目標最大 ≥ am > MHI
                    continue
                if am < MLO and kmax < MLO:      # 兩者皆低於下界
                    continue
                mx = np.maximum(kpms, am)
                k = (mx >= MLO) & (mx <= MHI)
                if not k.any():
                    continue
                Sk = Ss[k]; phk = phs[k]; noak = noa[k]; mpsk = mps[k]
                addv = MULTI3[ai]
                W = float(WT13[av - 1, ai])
                kk = Sk.shape[0]
                noC = np.ones(kk); fpC = np.ones(kk)
                as_ = np.sort(addv)
                for j in range(3):
                    v = as_[j]
                    cS = (Sk == v).sum(1).astype(np.float64)
                    cb = 0.0
                    for jj in range(j):
                        cb += float(as_[jj] == v)
                    noC *= (cS + cb + 1.0)
                    fpC *= (cb + 1.0)
                T = np.sort(np.concatenate([Sk, np.broadcast_to(addv, (kk, 3))], axis=1), axis=1)
                r = C._rank(T) - RLO
                amp = (W / 90.0) * np.sqrt(noak * noC) * phk / (fpC * mpsk)
                np.add.at(out, r, amp)


def bounds_for(NW):
    """依「add 三元最大元素 = M」之數量 M(M+1)/2 均衡切連續 M 範圍
    （工作正比於 add 數，非狀態數；記憶體仍 ≤ 最大單塊 M=44 之 18.6 GiB）"""
    cost = {M: M * (M + 1) // 2 for M in range(1, N + 1)}
    tot = sum(cost.values())
    b = [1]; accc = 0
    for M in range(1, N + 1):
        accc += cost[M]
        rem_workers = NW - (len(b) - 1)
        if rem_workers > 1 and accc >= tot * (len(b)) / NW and (N + 1 - M) >= rem_workers - 1:
            b.append(M + 1)
    b.append(N + 1)
    while len(b) - 1 < NW:            # 補足空位（理論上不會）
        b.insert(-1, b[-2] + 1)
    return b


def main():
    WK = int(sys.argv[1]); NW = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    T0 = time.time()
    log = lambda m: print(f"[w{WK} {time.time()-T0:7.1f}s] {m}", flush=True)
    b = bounds_for(NW)
    MLO, MHI = b[WK], b[WK + 1] - 1
    RLO, RHI = EDGE[MLO - 1], EDGE[MHI]
    log(f"M∈[{MLO},{MHI}] rank∈[{RLO:,},{RHI:,}) 態數 {RHI-RLO:,}（{(RHI-RLO)*8/2**30:.2f} GiB）")
    P79 = f"{CK}/p79_w{WK}.npy"; J79 = f"{CK}/p79_w{WK}.json"
    if os.path.exists(P79) and os.path.exists(J79):
        led = json.load(open(J79)); done = set(led["done"])
        out = np.load(P79, mmap_mode="r+"); out = np.array(out)
        log(f"續跑：已完成 {len(done)} 塊")
    else:
        out = np.zeros(RHI - RLO, dtype=np.float64); done = set()
    tab = np.load(TABLE7, mmap_mode="r"); phi = np.load(PHI37)
    CHUNK = int(os.environ.get("CHUNK", "2000000"))
    NCH = (NQ7 + CHUNK - 1) // CHUNK
    t00 = time.time()
    for ci in range(NCH):
        if ci in done:
            continue
        c0, c1 = ci * CHUNK, min((ci + 1) * CHUNK, NQ7)
        S = np.array(tab[c0:c1], dtype=np.int32); ph = np.array(phi[c0:c1])
        m = ph != 0; S, ph = S[m], ph[m]
        acc(S, ph, MLO, MHI, out, RLO)
        done.add(ci)
        np.save(P79, out)
        json.dump({"wk": WK, "MLO": MLO, "MHI": MHI, "RLO": int(RLO), "RHI": int(RHI),
                   "done": sorted(done), "nchunk": NCH}, open(J79, "w"))
        dt = time.time() - t00; nd = len(done)
        log(f"塊 {ci+1}/{NCH}；{dt:.0f}s；ETA {dt/nd*(NCH-nd)/60:.1f} min")
    log("DONE")


if __name__ == "__main__":
    main()