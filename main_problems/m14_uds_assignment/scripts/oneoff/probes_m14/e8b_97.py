#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""路徑 B 之 V^{9→7}（未剪枝 (3,1) 四波 ＋ (2,0) 同模湮滅）

與 A 版 e8_97b.py 之差異：
  ① 讀 p79B_w<k>.npy（B 含同模，未剝）而非 p79f_w<k>.npy
  ② 加 (2,0) 同模湮滅通道：直接呼叫 e7chain 正典 sm_op（零自寫風險）
  ③ DEN₉ 位移：phi9 = raw / (γ₁ − E9e − 10)，10 = 9/2 + SH(5.5)
  ④ 輸出 p97B_<tag>.npy / .json，原子寫

閘（300 個隨機 9-態 vs C.op(3,1,7)+C.sm_op）：‖ref‖=‖mine‖=8.888447，最大差 0.000e+00

用法：e8b_97.py <tag> <planfile>
"""
import os, sys, json, time
import numpy as np
from math import comb
import itertools

for _p in ("/home/wsl/idcm_ckpt/e8", "/home/wsl/.hermes/cache/scratch"):
    if os.path.exists(_p + "/e7chain.py"):
        sys.path.insert(0, _p); break
import e7chain as C

N = 44
NQ7 = comb(50, 7)
NQ9 = comb(52, 9)
CK = os.environ.get("CK", "/home/wsl/idcm_ckpt/e8")
CKB = "/home/wsl/idcm_ckpt/e7/e7ckptB"
SH = 5.5
D9 = 9 / 2 + SH
G = C.G
OUTP = os.environ.get("OUT", "p97B")


def unrank9(r):
    r = np.asarray(r, dtype=np.int64).copy()
    out = np.empty((r.shape[0], 9), dtype=np.int32)
    for j in range(8, -1, -1):
        n = np.arange(0, N + 2, dtype=np.int64)
        cn = np.array([comb(int(x) + j - 1, j + 1) if int(x) + j - 1 >= 0 else 0 for x in n])
        idx = np.searchsorted(cn, r, side="right") - 1
        idx = np.clip(idx, 1, N)
        out[:, j] = idx
        r = r - cn[idx]
    return out


def chunk_apply(M9, phi9, out7):
    """(3,1) 四波未剪枝（每移除模式 × 44 add 全掃）＋ (2,0) 同模湮滅"""
    m = M9.shape[0]
    for rm in itertools.combinations(range(9), 3):
        keep = [j for j in range(9) if j not in rm]
        A = np.sort(M9[:, list(rm)], axis=1)
        Ssur = M9[:, keep]
        e01 = A[:, 0] == A[:, 1]; e12 = A[:, 1] == A[:, 2]
        al = e01 & e12
        fpA = np.where(al, 6.0, np.where((e01 | e12) & ~al, 2.0, 1.0))
        noA = np.ones(m); mulpos = np.ones(m)
        for j in range(3):
            v = A[:, j]
            cT = (M9 == v[:, None]).sum(1).astype(np.float64)
            cb = np.zeros(m)
            for jj in range(j):
                cb += (A[:, jj] == v)
            noA *= (cT - cb)
            mulpos *= (cT - cb) / (cb + 1.0)
        rkt = C._rank(A)
        for v in range(1, N + 1):
            W = C.WT31[v - 1, rkt].astype(np.float64)
            nz = W != 0
            if not nz.any():
                continue
            Sh = Ssur[nz]
            cS = (Sh == v).sum(1).astype(np.float64)
            T = np.sort(np.concatenate([Sh, np.full((Sh.shape[0], 1), v, dtype=np.int32)], axis=1), axis=1)
            r = C._rank(T)
            amp = (W[nz] / 90.0) * np.sqrt(noA[nz] * (cS + 1.0)) * phi9[nz] / (fpA[nz] * mulpos[nz])
            np.add.at(out7, r, amp)
    # (2,0) 同模湮滅：直接用 e7chain 正典 sm_op
    out7 += C.sm_op(M9, phi9, 7, False, NQ7)


def main():
    tag = int(sys.argv[1]); plan = json.load(open(sys.argv[2]))
    segs = plan[str(tag)] if isinstance(plan, dict) else plan[tag]
    T0 = time.time()
    log = lambda s: print(f"[B{tag} {time.time()-T0:7.1f}s] {s}", flush=True)
    log(f"區間 {segs}｜未剪枝 (3,1) ＋ (2,0)｜DEN9 位移 {D9}")
    out7 = np.zeros(NQ7)
    J = f"{CK}/{OUTP}_{tag}.json"
    done = {}
    if os.path.exists(J) and os.path.exists(f"{CK}/{OUTP}_{tag}.npy") \
            and os.path.getsize(f"{CK}/{OUTP}_{tag}.npy") == NQ7 * 8 + 128:
        d = json.load(open(J))
        done = {}
        for _k, _v in d["done"].items():
            _p = _k.split(",")
            if len(_p) != 2 or not all(_x.isdigit() for _x in _p):
                raise SystemExit(f"FATAL ledger 鍵格式錯: {_k!r}（degenerate key）")
            done[(int(_p[0]), int(_p[1]))] = int(_v)
        _segset = {(int(a), int(b)) for a, b in segs}
        _bad = [k for k in done if k not in _segset]
        if _bad:
            raise SystemExit(f"FATAL ledger 鍵不屬 plan 段: {_bad[:3]} ⟹ 禁止續跑（否則重複計數）")
        out7 = np.load(f"{CK}/{OUTP}_{tag}.npy")
        log(f"續跑已做 {len(done)} 段")
    CH = 2_000_000
    for (lo, hi) in segs:
        nch = (hi - lo + CH - 1) // CH
        dk = done.get((lo, hi), 0)
        for ci in range(dk, nch):
            r0 = lo + ci * CH; r1 = min(r0 + CH, hi)
            ranks = np.arange(r0, r1, dtype=np.int64)
            M9 = unrank9(ranks)
            rawv = load_raw(ranks)
            E9e = G[M9].sum(1)
            phi9 = rawv / (G[1] - E9e - D9)
            chunk_apply(M9, phi9, out7)
            done[(lo, hi)] = ci + 1
            _tmpf = f"{CK}/{OUTP}_{tag}.tmp.npy"
            np.save(_tmpf, out7)
            os.replace(_tmpf, f"{CK}/{OUTP}_{tag}.npy")
            json.dump({"tag": tag, "done": {f"{k[0]},{k[1]}": v for k, v in done.items()}},
                      open(J + ".tmp", "w"))
            os.replace(J + ".tmp", J)
            el = time.time() - T0
            tot_done = sum(done.values()); tot_all = sum((b - a + CH - 1) // CH for a, b in segs)
            log(f"段({lo:,},{hi:,}) 塊 {ci+1}/{nch}｜{tot_done}/{tot_all}｜{el:.0f}s")
    log("DONE ‖out7‖=%.9e" % np.linalg.norm(out7))


_SHARDS = None


def _shards():
    global _SHARDS
    if _SHARDS is None:
        out = []
        for k in range(12):
            jf = f"{CK}/p79B_w{k}.json"
            if not os.path.exists(jf):
                continue
            with open(jf) as fh:
                led = json.load(fh)
            out.append((k, led["RLO"], led["RHI"]))
        _tot = sum(b - a for _, a, b in out)
        if _tot != NQ9:
            raise SystemExit(f"FATAL 分片覆蓋 {_tot:,} ≠ NQ9 {NQ9:,} ⟹ load_raw 會靜靜回零")
        _SHARDS = out
    return _SHARDS


def load_raw(ranks):
    n = ranks.shape[0]
    out = np.zeros(n)
    if n == 0:
        return out
    contig = (n == 1) or bool(np.all(np.diff(ranks) == 1))
    if contig:
        r0 = int(ranks[0]); r1 = int(ranks[-1]) + 1
        for k, RLO, RHI in _shards():
            if RLO <= r0 and r1 <= RHI:
                a = np.load(f"{CK}/p79B_w{k}.npy", mmap_mode="r")
                v = np.array(a[r0 - RLO: r1 - RLO], dtype=np.float64)
                del a
                return v
    for k, RLO, RHI in _shards():
        sel = (ranks >= RLO) & (ranks < RHI)
        if not sel.any():
            continue
        a = np.load(f"{CK}/p79B_w{k}.npy", mmap_mode="r")
        out[sel] = a[ranks[sel] - RLO]
        del a
    return out


if __name__ == "__main__":
    main()
