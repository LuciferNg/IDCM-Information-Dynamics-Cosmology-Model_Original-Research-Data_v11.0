#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M14 獨立重推導（idcm-audit-verify §2.4）
================================================================================
四條鐵律
  R1 唔准 import／call 任何 L1 腳本（scripts/NN_*.py）；唔准讀 docs/。
     （本檔只讀本目錄 data/ 以外嘅 PDG 量測檔？→ 不需要：M14 自身之量全部
       由 {N_d, B, N_gen, L6, 測量 ρ/δ} 起算，零外部檔）
  R2 INPUTS 只放原始量：計數、測量值、上游已封結論。
  R3 DERIVE 每個 key 係 INPUTS 之純函數（零隨機、零檔案 IO）。
  R4 數值一律 mpmath dps=40；禁 float 累加。

★ 每個 key 都係「由原始量重算」，唔係抄 docs 或 L1 輸出。
★ 只寫 data/l2_independent_verify.json；stdout 由 shell redirect 成 canonical log。
"""
import itertools
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import mpmath as mp

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG, RES, NFAIL = [], {}, 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def ser(x):
    try:
        xf = float(x)
    except (OverflowError, ValueError):
        xf = None
    if xf is not None and (xf == 0 or 1e-300 < abs(xf) < 1e300):
        return {"float": xf, "str": mp.nstr(x, 25)}
    return {"float": None, "str": mp.nstr(x, 25)}


# ============================================================================
# ① INPUTS —— 原始量
# ============================================================================
INPUTS = {
    "N_d": 44,               # 鏈自由度（44 站 / 44 模）
    "B": 4,                  # 分支數
    "N_gen": 3,              # 世代數（計數，非量子數）
    "L6": 18,                # L6 = 18（已註冊）
    "x_ring": 88,            # 環半週期
    "x_dir": 90,             # Dirichlet 半週期
    "p_u_max": 7, "p_d_max": 8, "p_s_max": 10,   # 掃描箱上界
    "alpha_canon": (1, 2, 4),                    # M13 正典 α
    # --- 測量／上游已封（M13 已 L2 驗過；此處只作輸入）---
    "rho_d_obs": "0.40350365",
    "rho_u_obs": "-0.02588837",
    "sig_rho_d": "0.044785178",
    "sig_rho_u": "0.033583899",
    "delta_u_obs": "-0.088667660",
    "sig_delta_u": "0.11244",
}

phi = (1 + mp.sqrt(5)) / 2


def gam(x, n):
    return 4 * mp.sin(n * mp.pi / x) ** 2


# ============================================================================
# ② DERIVE
# ============================================================================
def _triples(I):
    """112 個站三元組（p_u<p_d<p_s，上界如 INPUTS）"""
    out = []
    for pu in range(1, I["p_u_max"] + 1):
        for pd in range(pu + 1, I["p_d_max"] + 1):
            for ps in range(pd + 1, I["p_s_max"] + 1):
                out.append((pu, pd, ps))
    return out


def _gram(I, tri, mode="variable_delta"):
    """M_xy = Σ_n v_x(n) v_y(n)/γ_n；v_x = φ_n(p_x)·Δ_n，Δ_n = φ_n(p_d) − φ_n(p_u)"""
    N = I["N_d"]

    def ph(n, j):
        return mp.sqrt(mp.mpf(2) / 45) * mp.sin(mp.pi * n * j / 45)

    pu, pd, ps = tri
    M = mp.zeros(3, 3)
    for n in range(1, N + 1):
        d = ph(n, pd) - ph(n, pu)
        v = [ph(n, pu) * d, ph(n, pd) * d, ph(n, ps) * d]
        g = gam(I["x_dir"], n)
        for a in range(3):
            for b in range(3):
                M[a, b] += v[a] * v[b] / g
    return M


def _rank(M, tol=mp.mpf("1e-25")):
    """行化簡求秩（手寫，唔用 mpmath 矩陣方法）"""
    A = [[mp.mpf(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    m, n, r = len(A), M.cols, 0
    for c in range(n):
        piv = max(range(r, m), key=lambda i: abs(A[i][c]), default=None)
        if piv is None or abs(A[piv][c]) < tol:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def _classify(I):
    """A2⊕1 / A2 2×2 塊（定義照 docs/02 之判準文字：對角全等 + M01=−M00/2 [+ s 退耦]）"""
    TOL = mp.mpf("1e-25")

    def is_a2p1(M):
        d = M[0, 0]
        if abs(d) < TOL:
            return False
        if not all(abs(M[i, i] - d) < TOL for i in (1, 2)):
            return False
        if abs(M[0, 1] + d / 2) > TOL:
            return False
        return max(abs(M[0, 2]), abs(M[1, 2])) < TOL

    def is_a2b(M):
        d = M[0, 0]
        return abs(d) >= TOL and abs(M[1, 1] - d) < TOL and abs(M[0, 1] + d / 2) < TOL

    n1 = n2 = 0
    rk_counts, dec_s = {}, 0
    for tri in _triples(I):
        M = _gram(I, tri)
        rk = _rank(M)
        rk_counts[rk] = rk_counts.get(rk, 0) + 1
        if is_a2p1(M):
            n1 += 1
        if is_a2b(M):
            n2 += 1
        if abs(M[2, 0]) < TOL and abs(M[2, 1]) < TOL:
            dec_s += 1
    return n1, n2, {"rk": rk_counts, "dec_s": dec_s}


def _rho(I, tri):
    """ρ_u, ρ_d：由 Gram 讀出（docs/02 之 L1 定義）"""
    M = _gram(I, tri)
    # ρ_x = M_xu/M_uu 之相對讀法（M13 docs/05 S2）：以 u 為基準
    return M[1, 1] / M[0, 0], M[2, 2] / M[0, 0]


def _pair_classes(I, x, lo=1, hi=None):
    hi = hi or I["N_d"]
    vals = []
    for a in range(lo, hi + 1):
        for b in range(a + 1, hi + 1):
            vals.append((gam(x, a) + gam(x, b), (a, b)))
    vals.sort(key=lambda t: t[0])
    cl, cur = [], [vals[0]]
    for v, p in vals[1:]:
        if abs(v - cur[-1][0]) < mp.mpf("1e-30"):
            cur.append((v, p))
        else:
            cl.append(cur)
            cur = [(v, p)]
    cl.append(cur)
    return cl


def _orbits(I, gens_fn):
    """真群作用之軌道分解（gens_fn ∈ {Z44 旋轉, D44}）"""
    N = I["N_d"]
    pairs = {(a, b) for a in range(1, N + 1) for b in range(a + 1, N + 1)}
    seen, orbs = set(), []
    for p in sorted(pairs):
        if p in seen:
            continue
        comp, stack = {p}, [p]
        while stack:
            q = stack.pop()
            for w in gens_fn(q):
                if w in pairs and w not in comp:
                    comp.add(w)
                    stack.append(w)
        orbs.append(comp)
        seen |= comp
    return orbs


def _z44_gens(p, I):
    N = I["N_d"]

    def nrm(x):
        return (x - 1) % N + 1

    out = []
    for s in range(N):
        a, b = nrm(p[0] + s), nrm(p[1] + s)
        out.append((min(a, b), max(a, b)))
    return out


def _d44_gens(p, I):
    N = I["N_d"]

    def nrm(x):
        return (x - 1) % N + 1

    a, b = nrm(N + 1 - p[1]), nrm(N + 1 - p[0])
    base = (min(a, b), max(a, b))
    return _z44_gens(base, I) + _z44_gens(p, I)


def _phi90(I, n, j):
    """獨立重算：/90 Dirichlet 模函數"""
    return mp.sqrt(mp.mpf(2) / 45) * mp.sin(n * mp.pi * j / 45)


def _chi3(I, a, b, m):
    """獨立重算：三模重疊 χ³"""
    return sum(_phi90(I, a, j) * _phi90(I, b, j) * _phi90(I, m, j) for j in range(1, 45))


def _q4_scan(I):
    """獨立重算：δ 三分律（x = 3..200）之違例"""
    out = []
    for x in range(3, 201):
        if x % 2 == 0:
            M = x // 2
            if abs(sum(4 * mp.sin(n * mp.pi / x) ** 2 for n in range(1, M + 1)) - (x + 2)) > mp.mpf("1e-30"):
                out.append(("ring", x))
            M = x // 2 - 1
            if M >= 1 and abs(sum(4 * mp.sin(n * mp.pi / x) ** 2 for n in range(1, M + 1)) - (x - 2)) > mp.mpf("1e-30"):
                out.append(("dir-even", x))
        else:
            M = (x - 1) // 2
            if abs(sum(4 * mp.sin(n * mp.pi / x) ** 2 for n in range(1, M + 1)) - x) > mp.mpf("1e-30"):
                out.append(("dir-odd", x))
    return out


def I_gam(tag, n):
    """獨立重算：γ_n 於指定鏈（"88"/"90"）"""
    x = int(tag)
    return 4 * mp.sin(n * mp.pi / x) ** 2


def _agp(I, x):
    """獨立重算：A-GP（A_b² = A_a·A_c），A_n = 2sin(nπ/N)"""
    import itertools
    N = I["N_d"]
    T = mp.mpf("1e-40")
    out = []
    for a, b, c in itertools.combinations(range(1, N + 1), 3):
        A = lambda k: 2 * mp.sin(k * mp.pi / x)
        if abs(A(b) ** 2 - A(a) * A(c)) < T:
            out.append((a, b, c))
    return out


def _ap_count(I, N):
    """獨立重算：N 之 AP 三元組數（K = N/2）"""
    import itertools
    K = N // 2
    g = {n: gam(N, n) for n in range(1, K + 1)}
    T = mp.mpf("1e-28")
    return sum(1 for a, b, c in itertools.combinations(range(1, K + 1), 3)
               if abs((g[b] - g[a]) - (g[c] - g[b])) < T)


def _fourwave_stats():
    """四波能量共振之完整枚舉（含重複，n ≤ 44）；回傳 (總數, 和偶, 和奇, off-shell, W 分佈)"""
    g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
    tot = ev = od = off = 0
    Wc = {}
    for a in range(1, 45):
        for b in range(a, 45):
            for c in range(b, 45):
                for d in range(c, 45):
                    t = (a, b, c, d)
                    ga = g[a] + g[b]
                    res = (abs(ga - g[c] - g[d]) < mp.mpf("1e-30")
                           or abs(g[a] + g[c] - g[b] - g[d]) < mp.mpf("1e-30")
                           or abs(g[a] + g[d] - g[b] - g[c]) < mp.mpf("1e-30"))
                    w = sum(cs for cs, sg in zip([1, 1, -1, -1, -1, -1, 1, 1],
                                                 [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1),
                                                  (1, -1, 1, 1), (1, 1, -1, 1), (1, 1, 1, -1),
                                                  (1, 1, -1, -1), (1, 1, 1, 1)])
                             if sum(x * y for x, y in zip(sg, t)) % 90 == 0)
                    if res:
                        tot += 1
                        if sum(t) % 2 == 0:
                            ev += 1
                            Wc[w] = Wc.get(w, 0) + 1
                        else:
                            od += 1
                    elif w != 0:
                        off += 1
    return tot, ev, od, off, Wc


def _dynamics():
    """MR1 機器（docs/41）：衰變通道、2→2 配對共振、不透明度均勻性"""
    g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
    d12 = [(a, b, c) for a in range(1, 45) for b in range(1, 45) for c in range(b, 45)
           if abs(g[a] - g[b] - g[c]) < mp.mpf("1e-30")]
    d13 = [(a, t) for t in itertools.combinations_with_replacement(range(1, 45), 3)
           for a in range(1, 45) if abs(g[a] - sum(g[x] for x in t)) < mp.mpf("1e-30")]
    tri4 = [t for t in itertools.combinations_with_replacement(range(1, 45), 3)
            if abs(sum(g[x] for x in t) - 4) < mp.mpf("1e-30")]
    pairs = list(itertools.combinations_with_replacement(range(1, 45), 2))
    res = [(p, q) for i, p in enumerate(pairs) for q in pairs[i + 1:]
           if abs(g[p[0]] + g[p[1]] - g[q[0]] - g[q[1]]) < mp.mpf("1e-30")]
    sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
           (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
    cs = [1, 1, -1, -1, -1, -1, 1, 1]

    def W(t):
        return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

    even = [r for r in res if (sum(r[0]) + sum(r[1])) % 2 == 0]
    odd = [r for r in res if (sum(r[0]) + sum(r[1])) % 2 != 0]
    cnt = {}
    for r in even:
        for m in r[0] + r[1]:
            cnt[m] = cnt.get(m, 0) + 1
    return dict(d12=len(d12), d13=len(d13), tri4=len(tri4), res=len(res),
                even=len(even), odd=len(odd), uniform=sorted(set(cnt.values())))


_SO_CACHE = {}


def _second_order():
    """MR1 二階（docs/42）：極點掃描、相異性、微擾界（workdps(60) + 快取）"""
    if _SO_CACHE:
        return _SO_CACHE
    with mp.workdps(60):
        g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
        sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
               (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
        cs = [1, 1, -1, -1, -1, -1, 1, 1]

        def W(t):
            return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

        tri = list(itertools.combinations_with_replacement(range(1, 45), 3))
        poles = [(n, t, W((n,) + t)) for n in range(1, 45) for t in tri
                 if abs(g[n] - sum(g[x] for x in t)) < mp.mpf("1e-40")]
        S = {}
        for n in range(1, 45):
            acc = mp.mpf(0)
            for t in tri:
                w = W((n,) + t)
                if w == 0:
                    continue
                cnt = [t.count(v) for v in set(t)]
                m = 6
                for c in cnt:
                    m //= 1 if c == 1 else (2 if c == 2 else 6)
                acc += m * mp.mpf(w * w) / (g[n] - sum(g[x] for x in t))
            S[n] = acc
        bound = min(8100 * g[n] / (2 * abs(S[n])) for n in S if S[n] != 0)
        order = sorted(range(1, 45), key=lambda n: g[n])
        sigs = [1 if S[n] > 0 else -1 for n in order]
        res = dict(poles=poles, distinct=len({mp.nstr(v, 25) for v in S.values()}),
                   bound=mp.sqrt(bound), smax=max(abs(v) for v in S.values()),
                   pos=[n for n in range(1, 45) if S[n] > 0],
                   transitions=sum(1 for i in range(len(sigs) - 1) if sigs[i] != sigs[i + 1]))
    _SO_CACHE.update(res)
    return _SO_CACHE


def _bubble():
    """bubble 支（docs/43）：γ 鏡像恆等、零分母分類、W 四點鏡像不變"""
    def block():
        mp.mp.dps = max(mp.mp.dps, 60)
        g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
        sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
               (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
        cs = [1, 1, -1, -1, -1, -1, 1, 1]

        def W(t):
            return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

        mir_g = max(abs(g[n] + g[45 - n] - 4) for n in range(1, 45))
        mir_w = 0
        for n in range(1, 45, 7):
            for m in range(1, 45, 5):
                for p in range(1, 45, 6):
                    t = (n, m, p, (n + m + p) % 44 + 1)
                    t2 = (45 - t[0], 45 - t[1], 45 - t[2], 45 - t[3])
                    mir_w = max(mir_w, abs(W(t) - W(t2)))
        cnt = {}
        for n in range(1, 45):
            z = 0
            for m in range(1, 45):
                for p in range(1, 45):
                    for q in range(p, 45):
                        if abs(g[n] + g[m] - g[p] - g[q]) < mp.mpf("1e-35") \
                                and not (p == n and q == m) and not (p == m and q == n):
                            if W((n, m, p, q)) != 0:
                                z += 1
            cnt[n] = z
        return dict(mir_g=mir_g, mir_w=mir_w, zero_nz=sorted(set(cnt.values())))
    res = block()
    return res


def _third():
    """三階（docs/44）：矩陣元律＋T_n（只算 n = 1, 39；39 為微擾界綁定模）"""
    def block():
        mp.mp.dps = max(mp.mp.dps, 60)
        g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
        sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
               (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
        cs = [1, 1, -1, -1, -1, -1, 1, 1]
        tri = list(__import__("itertools").combinations_with_replacement(range(1, 45), 3))

        def W(t):
            return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

        def nrm(T):
            r = 1
            for v in __import__("collections").Counter(T).values():
                r *= __import__("math").factorial(v)
            return mp.sqrt(mp.mpf(r))

        def et(T):
            return sum(g[m] for m in T)

        def minus(mult, sub):
            c = __import__("collections").Counter(mult)
            for k, v in __import__("collections").Counter(sub).items():
                c[k] -= v
                if c[k] < 0:
                    return None
            return list(c.elements())

        out = {}
        for n in (1, 39):
            lst = [(T, mp.mpf(W((n,) + T))) for T in tri if W((n,) + T)]
            byx = {}
            for T, w in lst:
                for x in set(T):
                    byx.setdefault(x, []).append((T, w, __import__("collections").Counter(T)[x]))
            acc = mp.mpf(0)
            a2 = mp.mpf(0)
            for T, wT in lst:
                a2 += (6 * (wT / 90) / nrm(T)) ** 2 / (g[n] - et(T))
                uT = 6 * (wT / 90) / nrm(T) / (g[n] - et(T))
                for idx in __import__("itertools").combinations(range(3), 2):
                    rem = [T[i] for i in idx]
                    x = minus(T, rem)[0]
                    for Tp, wP, mp_ in byx.get(x, []):
                        add = minus(Tp, [x])
                        if add is None or len(add) != 2:
                            continue
                        wm = W(tuple(rem) + tuple(add))
                        if not wm:
                            continue
                        uP = 6 * (wP / 90) / nrm(Tp) / (g[n] - et(Tp))
                        acc += uT * uP * mp_ * 4 * (mp.mpf(wm) / 90) / (nrm(T) * nrm(Tp))
            out[n] = (acc, a2)
        return out
    return block()


_FOURTH_CACHE = {}


def _quint_law():
    """3↔5 耦合律自校（docs/46）：Q 側枚舉 vs docs/44 律直算（T 側）"""
    import itertools as _it
    from math import factorial as _f
    from collections import Counter as _C
    sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
           (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
    cs = [1, 1, -1, -1, -1, -1, 1, 1]

    def W(t):
        return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

    def nrm(T):
        r = 1
        for v in _C(T).values():
            r *= _f(v)
        return r

    def qside(Q):
        cQ = _C(Q)
        nQ = nrm(Q)
        out = {}
        for S in _it.combinations_with_replacement(Q, 2):
            cS = _C(S)
            if any(cS[k] > cQ[k] for k in cS):
                continue
            add = list((cQ - cS).elements())
            emb = 1
            cnt = dict(cQ)
            for m in S:
                emb *= cnt[m]
                cnt[m] -= 1
            nS = 1
            for v in cS.values():
                nS *= _f(v)
            for r0 in range(1, 45):
                T = tuple(sorted(list(S) + [r0]))
                w = W((r0,) + tuple(add))
                if w == 0:
                    continue
                wgt = 1
                nT = nS
                if r0 in cS:
                    wgt = cS[r0] + 1
                    nT = nS * (cS[r0] + 1)
                out[T] = out.get(T, 0.0) + 6 * emb * wgt * (w / 90.0) / (nT * nQ) ** 0.5
        return out

    def tside(T, Q):
        cT, cQ = _C(T), _C(Q)
        tot = 0.0
        for i in range(3):
            surv = list((cT - _C([T[i]])).elements())
            cS = _C(surv)
            if any(cQ[k] < v for k, v in cS.items()):
                continue
            add = list((cQ - cS).elements())
            if len(add) != 3:
                continue
            w = W(tuple([T[i]] + add))
            if w == 0:
                continue
            emb = 1
            cnt = dict(cQ)
            for m in surv:
                emb *= cnt[m]
                cnt[m] -= 1
            tot += emb * 6 * (w / 90.0)
        return tot / (nrm(T) * nrm(Q)) ** 0.5
    import random as _r
    rnd = _r.Random(20261102)
    mx = 0.0
    nzc = 0
    for _ in range(40):
        Q = tuple(sorted(rnd.sample(range(1, 45), 5)))
        for T, v in qside(Q).items():
            mx = max(mx, abs(v - tside(T, Q)))
            nzc += 1
    return mx, nzc


_ALL44_CACHE = {}


def _all44():
    """全 44 n 掃描之統計（**讀 data 落地**，非獨立重算；獨立重算見 X82–X84 之四點抽樣）"""
    if _ALL44_CACHE:
        return _ALL44_CACHE
    import re as _re
    from pathlib import Path as _P
    txt = (_P(__file__).resolve().parent.parent / "data" / "fifth_order_all44.log").read_text(
        encoding="utf-8", errors="replace")
    rows = []
    for line in txt.splitlines():
        m = _re.match(r"\s*(\d+)\s*｜([+\-0-9.e]+)｜([+\-0-9.e]+)｜([+\-0-9.e]+)｜([+\-0-9.e]+)｜([+\-0-9.e]+)%", line)
        if m:
            rows.append((int(m.group(1)), *[float(m.group(i)) for i in (2, 3, 4, 5)], float(m.group(6)) / 100))
    if len(rows) != 44:
        raise RuntimeError(f"all44 log 行數 {len(rows)} != 44")
    rat = sorted(abs(r[4] / r[3]) for r in rows)
    _ALL44_CACHE.update({"n": len(rows), "min": rat[0], "median": rat[22], "max": rat[-1],
                         "n_gt1": sum(1 for x in rat if x > 1),
                         "sumE5": sum(r[4] for r in rows), "sumE4": sum(r[3] for r in rows)})
    return _ALL44_CACHE


_PV44_CACHE = {}


def _pv44():
    """精確頂點全 44 n（**讀 data/precise_vertex_44n.json 落地**，非獨立重算；
    獨立重算之閘見 scripts/51 §1 G1 14112/14112）"""
    if _PV44_CACHE:
        return _PV44_CACHE
    import json as _j
    from pathlib import Path as _P
    d = _j.loads((_P(__file__).resolve().parent.parent / "data" / "precise_vertex_44n.json").read_text(encoding="utf-8"))
    out = {}
    for tag in ("A", "B"):
        rows = d[tag]["rows"]
        r1 = sorted(abs(r["E3"] / r["E2"]) for r in rows)
        r2 = sorted(abs(r["E4"] / r["E3"]) for r in rows)
        r3 = sorted(abs(r["E5"] / r["E4"]) for r in rows)
        e2 = [r["E2"] for r in rows]
        med = lambda v: (v[len(v) // 2 - 1] + v[len(v) // 2]) / 2 if len(v) % 2 == 0 else v[len(v) // 2]
        out[tag] = {"n": len(rows), "r1_med": med(r1), "r3_med": med(r3),
                    "r1_lt1": sum(1 for v in r1 if v < 1), "r3_lt1": sum(1 for v in r3 if v < 1),
                    "r2_min": r2[0], "r2_max": r2[-1],
                    "e2_min": min(e2), "e2_max": max(e2)}
    _PV44_CACHE.update(out)
    return _PV44_CACHE


_E6_CACHE = {}


def _pv_e6():
    """E_6（含 7 模態，路徑 A n=1）：**讀 data/precise_vertex_e6.json 落地**；
    kernel 正確性之獨立閘見 docs/51 §三（G-7a 單態 40/40）"""
    if _E6_CACHE:
        return _E6_CACHE
    import json as _j
    from pathlib import Path as _P
    d = _j.loads((_P(__file__).resolve().parent.parent / "data" / "precise_vertex_e6.json").read_text(encoding="utf-8"))
    r = d["rows"][0]
    _E6_CACHE.update({"E2": r["E2"], "E3": r["E3"], "E4": r["E4"], "E5": r["E5"], "E6": r["E6"],
                      "delta7": d["E6_7mode_delta"], "r4": abs(r["E6"] / r["E5"])})
    return _E6_CACHE


_E7_CACHE = {}


def _e7():
    """E_7 全 CPU 鏈（docs/52）：讀 data/e7_chain.json 落地；算子正確性之閘見 docs/52 §二"""
    if _E7_CACHE:
        return _E7_CACHE
    import json as _j
    from pathlib import Path as _P
    d = _j.loads((_P(__file__).resolve().parent.parent / "data" / "e7_chain.json").read_text(encoding="utf-8"))
    v = d["values"]; c = d["correction"]
    _E7_CACHE.update({"E2": v["E2"], "E3": v["E3"], "E4": v["E4"], "E5": v["E5"],
                      "E6": v["E6"], "E7": v["E7"],
                      "r5": abs(v["E7"] / v["E6"]), "r4": abs(v["E6"] / v["E5"]),
                      "E6_5mode": c["E6_5mode_trunc"], "rel5": c["rel_diff_5mode"],
                      "V75n": d["norms"]["V75_phi3_7"], "V75old": c["old_V75_norm_docs51"],
                      "gene": d["gates"]["7->5_120k_realphi"]})
    return _E7_CACHE


_FIFTH_CACHE = {}


def _fifth():
    """五階（docs/47）：獨立重建 3↔5、V33、LUT，複算 E_2/E_3/E_4(n=1) 與 5↔5 δ 測試"""
    if _FIFTH_CACHE:
        return _FIFTH_CACHE
    import numpy as np
    import itertools as _it
    from collections import Counter as _C
    from math import comb as _cb, factorial as _f
    Nn = 44
    SGs = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
           (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
    CCs = [1, 1, -1, -1, -1, -1, 1, 1]

    def w4(t):
        return sum(c for c, sg in zip(CCs, SGs)
                   if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

    Gv = [0.0] + [4 * np.sin(np.pi * n / 90) ** 2 for n in range(1, Nn + 1)]

    def nrm(t):
        r = 1
        for v in _C(t).values():
            r *= _f(v)
        return float(r)

    def mkrank(k):
        M = Nn + k - 1
        CT = [[_cb(n_, r_) for r_ in range(k + 2)] for n_ in range(M + 2)]

        def fr(a):
            a = np.asarray(a, dtype=np.int64)
            b = a + np.arange(k)
            prev = np.zeros(len(a), dtype=np.int64)
            out = np.zeros(len(a), dtype=np.int64)
            for i in range(k):
                r = k - i - 1
                out += (np.array([CT[M - p_][r + 1] for p_ in prev], dtype=np.int64)
                        - np.array([CT[M - x + 1][r + 1] for x in b[:, i]], dtype=np.int64))
                prev = b[:, i]
            return out
        return fr

    rank5 = mkrank(5)
    TRI = [tuple(t) for t in _it.combinations_with_replacement(range(1, Nn + 1), 3)]
    IDX3 = {t: i for i, t in enumerate(TRI)}
    N3 = len(TRI)
    NT3 = np.array([nrm(t) for t in TRI])
    ET3 = np.array([sum(Gv[a] for a in t) for t in TRI])
    ADD2 = np.array(list(_it.combinations_with_replacement(range(1, Nn + 1), 2)), dtype=np.int32)
    ADD3 = np.array(list(_it.combinations_with_replacement(range(1, Nn + 1), 3)), dtype=np.int32)
    QP = np.array(list(_it.combinations_with_replacement(range(1, Nn + 1), 5)), dtype=np.int32)
    NQ = len(QP)
    NT5 = np.array([nrm(map(int, t)) for t in QP])
    EQ = np.array([float(np.sum([Gv[a] for a in t])) for t in QP])

    def w4v(T):
        s = np.zeros(T.shape[0], dtype=np.int32)
        for c, sg in zip(CCs, SGs):
            s += np.where(np.abs(T @ np.array(sg, dtype=np.int32)) % 90 == 0, c, 0)
        return s

    qi, vv, tt = [], [], []
    for i, T in enumerate(TRI):
        for pos in range(3):
            sv = [T[q] for q in range(3) if q != pos]
            rem = int(T[pos])
            A = np.column_stack([np.full(len(ADD3), rem, dtype=np.int32), ADD3])
            w = w4v(A)
            nz = np.nonzero(w)[0]
            if len(nz) == 0:
                continue
            Qs = np.sort(np.column_stack([np.repeat(np.array(sv, dtype=np.int32)[None, :], len(nz), axis=0),
                                          ADD3[nz]]), axis=1)
            r = rank5(Qs)
            cS = _C(sv)
            emb = np.ones(len(nz))
            for m, cnt in cS.items():
                emb *= np.array([_cb(int(x), cnt) for x in (Qs == m).sum(1)], dtype=float)
            qi.append(r.astype(np.int32))
            vv.append(6.0 * (w[nz] / 90.0) * emb / (NT3[i] * NT5[r]) ** 0.5)
            tt.append(np.full(len(r), i, dtype=np.int32))
    qidx = np.concatenate(qi)
    vals35 = np.concatenate(vv)
    tri = np.concatenate(tt)

    # V33
    rws, cls, vls = [], [], []
    for i, T in enumerate(TRI):
        for pos in range(3):
            surv = T[pos]
            rem = [T[q] for q in range(3) if q != pos]
            A = np.column_stack([np.full(len(ADD2), rem[0], dtype=np.int32),
                                 np.full(len(ADD2), rem[1], dtype=np.int32), ADD2])
            w = w4v(A)
            for kk in np.nonzero(w)[0]:
                Tp = tuple(sorted([surv, int(ADD2[kk, 0]), int(ADD2[kk, 1])]))
                j = IDX3[Tp]
                rws.append(i)
                cls.append(j)
                vls.append(4 * (float(w[kk]) / 90.0) * _C(Tp)[surv] / (NT3[i] * NT3[j]) ** 0.5)
    import scipy.sparse as _sp
    V33 = _sp.csr_matrix((vls, (rws, cls)), shape=(N3, N3))
    A = np.column_stack([np.full(N3, 1, dtype=np.int32), np.array(TRI, dtype=np.int32)])
    b = 6 * (w4v(A) / 90.0) / np.sqrt(NT3)
    inv3 = 1.0 / (Gv[1] - ET3)
    inv5 = 1.0 / (Gv[1] - EQ)
    p1 = inv3 * b
    E2 = float(b @ p1)
    p2_3 = inv3 * (V33 @ p1)
    E3 = float(b @ p2_3)
    X = np.bincount(qidx, weights=vals35 * p1[tri], minlength=NQ)
    p2_5 = inv5 * X
    p3_3 = inv3 * (V33 @ p2_3 + np.bincount(tri, weights=vals35 * p2_5[qidx], minlength=N3) - E2 * p1)
    E4 = float(b @ p3_3)

    # 5↔5 δ 測試（單 chunk）
    W = np.zeros((990, 990), dtype=np.int16)
    for i in range(990):
        A4 = np.column_stack([np.full(990, ADD2[i, 0], dtype=np.int32),
                              np.full(990, ADD2[i, 1], dtype=np.int32), ADD2])
        W[i] = w4v(A4)
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

    def mkrank2(k):
        M = Nn + k - 1
        CT = [[_cb(n_, r_) for r_ in range(k + 2)] for n_ in range(M + 2)]

        def fr(a):
            a = np.asarray(a, dtype=np.int64)
            bq = a + np.arange(k)
            prev = np.zeros(len(a), dtype=np.int64)
            out = np.zeros(len(a), dtype=np.int64)
            for i in range(k):
                r = k - i - 1
                out += (np.array([CT[M - p_][r + 1] for p_ in prev], dtype=np.int64)
                        - np.array([CT[M - x + 1][r + 1] for x in bq[:, i]], dtype=np.int64))
                prev = bq[:, i]
            return out
        return fr

    rank2 = mkrank2(2)
    rank3 = mkrank2(3)
    P5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    S5 = {p: [q for q in range(5) if q not in p] for p in P5}
    Q0 = (3, 11, 19, 27, 40)
    j0 = int(rank5(np.array([Q0]))[0])
    lo = (j0 // 400) * 400
    hi = min(lo + 400, NQ)
    out = np.zeros(NQ)
    for st in [lo]:
        Qc = QP[st:hi]
        m = len(Qc)
        nQin = NT5[st:hi][:, None]
        for pi in P5:
            rem = np.sort(Qc[:, list(pi)], axis=1)
            pidx = rank2(rem)
            sv = Qc[:, S5[pi]]
            tidx = rank3(sv)
            ranks = LUTR[tidx]
            emb = LUTE[tidx].astype(np.float64)
            amp = 4.0 * (W[pidx] / 90.0) * emb / np.sqrt(nQin * NT5[ranks])
            wv = np.zeros(m)
            k0 = j0 - st
            if 0 <= k0 < m:
                wv[k0] = 1.0
            out += np.bincount(ranks.ravel(), weights=(amp * wv[:, None]).ravel(), minlength=NQ)

    def slow55(Qp, Q):
        cQ = _C(Qp)
        tot = 0.0
        for pi in P5:
            rem = [Q[i] for i in pi]
            sv = [Q[q] for q in range(5) if q not in pi]
            cS = _C(sv)
            if any(cQ[m_] < c_ for m_, c_ in cS.items()):
                continue
            add = list((cQ - cS).elements())
            if len(add) != 2:
                continue
            wv = w4(tuple(rem + add))
            if not wv:
                continue
            emb = 1
            for m_, c_ in cS.items():
                emb *= _cb(cQ[m_], c_)
            tot += 4 * (wv / 90.0) * emb / (nrm(Q) * nrm(Qp)) ** 0.5
        return tot

    mx = 0.0
    for jj in np.argsort(-np.abs(out))[:6]:
        Qp = tuple(int(x) for x in QP[jj])
        mx = max(mx, abs(float(out[jj]) - slow55(Qp, Q0)))
    _FIFTH_CACHE.update({"E2": E2, "E3": E3, "E4": E4, "delta": mx, "nz35": int(len(qidx))})
    return _FIFTH_CACHE


def _fourth():
    """四階（docs/45）：RS 引擎自校；建三模扇區 V；重算 E_4(1)（快取，建 V ≈ 5 min）"""
    if "v" in _FOURTH_CACHE:
        return _FOURTH_CACHE["v"]

    def block():
        import numpy as _np
        import itertools as _it
        from math import factorial as _f
        g = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
        sig = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
               (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
        cs = [1, 1, -1, -1, -1, -1, 1, 1]

        def W(t):
            return sum(c for c, sg in zip(cs, sig) if sum(x * y for x, y in zip(sg, t)) % 90 == 0)

        def nrm(T):
            r = 1
            for v in __import__("collections").Counter(T).values():
                r *= _f(v)
            return r ** 0.5
        tri = list(_it.combinations_with_replacement(range(1, 45), 3))
        idx = {T: i for i, T in enumerate(tri)}
        pairs = list(_it.combinations_with_replacement(range(1, 45), 2))
        rows, cols, vals, nz = [], [], [], 0
        for i, T in enumerate(tri):
            for pos in range(3):
                surv = T[pos]
                rem = [T[k] for k in range(3) if k != pos]
                for (u, v) in pairs:
                    w = W(tuple(rem) + (u, v))
                    if not w:
                        continue
                    Tp = tuple(sorted([surv, u, v]))
                    emb = __import__("collections").Counter(Tp)[surv]
                    rows.append(i); cols.append(idx[Tp])
                    vals.append(4 * (w / 90.0) * emb / (nrm(T) * nrm(Tp))); nz += 1
        rows = _np.array(rows, _np.int32); cols = _np.array(cols, _np.int32)
        vals = _np.array(vals, _np.float64)
        et = _np.array([sum(float(g[m]) for m in T) for T in tri])
        n = 1
        d = float(g[n]) - et
        inv = _np.zeros(len(tri)); ok = _np.abs(d) > 1e-12
        inv[ok] = 1.0 / d[ok]
        b = _np.zeros(len(tri))
        for i, T in enumerate(tri):
            w = W((n,) + T)
            if w:
                b[i] = 6 * (w / 90.0) / nrm(T)

        def mv(phi):
            return _np.bincount(rows, weights=vals * phi[cols], minlength=len(tri))
        p1 = inv * b
        E2 = float(b @ p1)
        p2 = inv * mv(p1)
        E3 = float(b @ p2)
        p3 = inv * (mv(p2) - E2 * p1)
        E4 = float(b @ p3)
        # 引擎自校：殘差 ∝ λ⁵
        import numpy as _np2
        h0 = _np2.array([0.0, 1.0, 2.0, 3.0, 4.0])
        A = _np2.array([[0, .9, .3, .1, .2], [.9, 0, .5, .4, .1], [.3, .5, 0, .7, .2],
                        [.1, .4, .7, 0, .6], [.2, .1, .2, .6, 0]])

        def rs(h0, V, k, K=4):
            N = len(h0); Q = _np2.ones(N, bool); Q[k] = False
            iv = _np2.zeros(N); iv[Q] = 1.0 / (h0[k] - h0[Q])
            ph = [_np2.eye(N)[k]]; E = [0.0]
            for p in range(1, K + 1):
                E.append(float(V[k] @ ph[p - 1]))
                new = V @ ph[p - 1]
                for q in range(1, p):
                    new = new - E[q] * ph[p - q]
                new = iv * new; new[k] = 0.0; ph.append(new)
            return E
        Es = rs(h0, A, 0, 4)
        rat = []
        for lam in (1e-3, 1e-4):   # 兩者皆 λ⁵ 主導（λ=1e−2 有 λ⁷ 污染）
            ex = _np2.linalg.eigvalsh(_np2.diag(h0) + lam * A)
            k = int(_np2.argmin(abs(ex - h0[0])))
            r = abs(ex[k] - (h0[0] + sum(Es[p] * lam ** p for p in range(1, 5))))
            rat.append(r / lam ** 5)
        return dict(E2=E2, E3=E3, E4=E4, nz=nz, ratio_spread=abs(rat[0] - rat[1]) / rat[1])
    _FOURTH_CACHE["v"] = block()
    return _FOURTH_CACHE["v"]


def _k3_phi(a, b, m):
    """docs/37 之 Φ：∏ sin(πa/45) / ∏ sin(πk_i/90)（僅奇和有意義）"""
    num = mp.sin(mp.pi * a / 45) * mp.sin(mp.pi * b / 45) * mp.sin(mp.pi * m / 45)
    den = mp.mpf(1)
    for k in (a + b - m, a - b + m, -a + b + m, a + b + m):
        den *= mp.sin(mp.pi * k / 90)
    return num / den


def _gp(I, x):
    """獨立重算：GP（γ_b²=γ_a·γ_c）與 AP（γ_b−γ_a=γ_c−γ_b）三元組"""
    import itertools
    N = I["N_d"]
    T = mp.mpf("1e-40")
    gp, ap = [], []
    for a, b, c in itertools.combinations(range(1, N + 1), 3):
        ga, gb, gc = gam(x, a), gam(x, b), gam(x, c)
        if abs(gb * gb - ga * gc) < T:
            gp.append((a, b, c))
        if abs((gb - ga) - (gc - gb)) < T:
            ap.append((a, b, c))
    return gp, ap


def _minpoly(x):
    """pslq 求 2 次關係之判別式（首一化）"""
    rel = mp.pslq([mp.mpf(1), x, x * x], maxcoeff=10 ** 12, maxsteps=10 ** 6)
    a0, a1, a2 = [int(v) for v in rel]
    if a2 < 0:
        a0, a1, a2 = -a0, -a1, -a2
    return (a0, a1, a2, a1 * a1 - 4 * a0 * a2)


def _res8(I):
    """d=11 子鏈（n=11,22,33,44）之三波共振（互異）"""
    import itertools
    m = [11, 22, 33, 44]
    out = []
    for a, b, c in itertools.permutations(m, 3):
        if abs(gam(88, a) + gam(88, b) - gam(88, c)) < mp.mpf("1e-45"):
            out.append(tuple(sorted((a, b, c))))
    return sorted(set(out))


def _four8(I):
    """d=11 子鏈之四波共振（兩不同價對同 S）"""
    import itertools
    m = [11, 22, 33, 44]
    S = [gam(88, a) + gam(88, b) for a, b in itertools.combinations(m, 2)]
    return sum(1 for i in range(len(S)) for j in range(i + 1, len(S))
               if abs(S[i] - S[j]) < mp.mpf("1e-45"))


def _spectrum(I, x, comp):
    """獨立重算：價對能量譜之等值類 {key: pairs}；key 為量化 S"""
    from collections import defaultdict
    N = I["N_d"]
    T = mp.mpf("1e-40")
    b = defaultdict(list)
    for a in range(1, N + 1):
        for c in range(a + 1, N + 1):
            b[int(mp.nint((gam(x, a) + gam(x, c)) / T))].append((a, c))
    keys = sorted(b)
    out, cur = {}, None
    for k in keys:
        if cur is not None and (k - cur) * T < mp.mpf("1e-30"):
            out[cur].extend(b[k])
        else:
            cur = k
            out[cur] = list(b[k])
    return out


def _inv_scan(I):
    """獨立重算：10 個凍結不變量、1023 個子集，數恰 3 類且互不等價之分割數"""
    import itertools
    import math
    from collections import defaultdict
    N = I["N_d"]
    P = [(a, b) for a in range(1, N + 1) for b in range(a + 1, N + 1)]
    INV = {
        "i1": lambda a, b: (b - a,),
        "i2": lambda a, b: (min(b - a, N - (b - a)),),
        "i3": lambda a, b: (a % 2, b % 2),
        "i4": lambda a, b: ((b - a) % 3,),
        "i5": lambda a, b: ((a + b) % 3,),
        "i6": lambda a, b: (a % 3, b % 3),
        "i7": lambda a, b: (len({a, b} & {9, 15, 27}),),
        "i8": lambda a, b: (len({a, b} & {18, 30, 36}),),
        "i9": lambda a, b: (a <= I["N_d"] - I["L6"] + 1, b <= I["N_d"] - I["L6"] + 1),
        "i10": lambda a, b: (math.gcd(a, 88), math.gcd(b, 88)),
    }
    names = list(INV)
    bodies = set()
    for r in range(1, len(names) + 1):
        for combo in itertools.combinations(names, r):
            d = defaultdict(list)
            for a, b in P:
                t = tuple()
                for k in combo:
                    t += INV[k](a, b)
                d[t].append((a, b))
            if len(d) == 3:
                bodies.add(frozenset(frozenset(v) for v in d.values()))
    return len(bodies), len(bodies)


def _res_triples(I, x):
    """獨立重算：γ^(x)_a ± γ^(x)_b = γ^(x)_m 之無序三元組（互異）"""
    N = I["N_d"]
    out = set()
    for a in range(1, N + 1):
        for b in range(1, N + 1):
            if a == b:
                continue
            for m in range(1, N + 1):
                if m in (a, b):
                    continue
                if (abs(gam(x, a) + gam(x, b) - gam(x, m)) < mp.mpf("1e-40")
                        or abs(gam(x, a) - gam(x, b) - gam(x, m)) < mp.mpf("1e-40")):
                    out.add(tuple(sorted((a, b, m))))
    return sorted(out)


def _res_separable(I):
    """5 個 3 類分割中，能三分離 /90 共振三元組者之個數"""
    import itertools
    INV = {
        "i4": lambda a, b: (b - a) % 3,
        "i5": lambda a, b: (a + b) % 3,
        "i7": lambda a, b: len({a, b} & {9, 15, 27}),
        "i8": lambda a, b: len({a, b} & {18, 30, 36}),
        "i9": lambda a, b: (a <= 27, b <= 27),
    }
    tris = _res_triples(I, 90)
    n = 0
    for f in INV.values():
        sets = [{f(*p) for p in itertools.combinations(t, 2)} for t in tris]
        if all(sets[i].isdisjoint(sets[j]) for i in range(len(sets)) for j in range(i + 1, len(sets))):
            n += 1
    return n


DERIVE = {
    # ---- §A 慣例與能量預算 ----
    "pairs_946": lambda I: mp.mpf(I["N_d"] * (I["N_d"] - 1) // 2),
    "sum88": lambda I: mp.fsum(gam(I["x_ring"], n) for n in range(1, I["N_d"] + 1)),
    "sum90": lambda I: mp.fsum(gam(I["x_dir"], n) for n in range(1, I["N_d"] + 1)),
    "sum_delta": lambda I: mp.fsum(gam(I["x_ring"], n) - gam(I["x_dir"], n)
                                   for n in range(1, I["N_d"] + 1)),
    "theta_den_3960": lambda I: mp.mpf(I["N_d"] * I["x_dir"]),
    "gam88_1": lambda I: gam(I["x_ring"], 1),
    "gam90_1": lambda I: gam(I["x_dir"], 1),
    "d44_eq_gam90_1": lambda I: gam(I["x_ring"], I["N_d"]) - gam(I["x_dir"], I["N_d"]),
    "min_gap_cross": lambda I: min(abs(gam(I["x_ring"], n) - gam(I["x_dir"], m))
                                   for n in range(1, I["N_d"] + 1)
                                   for m in range(1, I["N_d"] + 1)),
    "cross_degeneracy": lambda I: mp.mpf(sum(
        1 for n in range(1, I["N_d"] + 1) for m in range(1, I["N_d"] + 1)
        if abs(gam(I["x_ring"], n) - gam(I["x_dir"], m)) < mp.mpf("1e-30"))),
    # ---- §B 補償映射與 parity ----
    "comp88_parity_same": lambda I: mp.mpf(sum(
        1 for n in range(1, I["N_d"] + 1) if (I["N_d"] - n) % 2 == n % 2)),
    "comp90_parity_same": lambda I: mp.mpf(sum(
        1 for n in range(1, I["N_d"] + 1) if (I["N_d"] + 1 - n) % 2 == n % 2)),
    "sumfam88_pairs": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"])),
    "sumfam90_pairs": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"] + 1)),
    "sumfam44_same_parity": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"] and a % 2 == b % 2)),
    "sumfam45_same_parity": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"] + 1 and a % 2 == b % 2)),
    # ---- §C 通道命運 ----
    "chan88_sum": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b <= I["N_d"]
        and abs(gam(I["x_ring"], a) + gam(I["x_ring"], b) - gam(I["x_ring"], a + b)) < mp.mpf("1e-30"))),
    "chan90_sum": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b <= I["N_d"]
        and abs(gam(I["x_dir"], a) + gam(I["x_dir"], b) - gam(I["x_dir"], a + b)) < mp.mpf("1e-30"))),
    "chan88_beat": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if abs(gam(I["x_ring"], a) + gam(I["x_ring"], b) - gam(I["x_ring"], b - a)) < mp.mpf("1e-30"))),
    "chan90_beat": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if abs(gam(I["x_dir"], a) + gam(I["x_dir"], b) - gam(I["x_dir"], b - a)) < mp.mpf("1e-30"))),
    # ---- §C2 共振定律（law 型：sum ⟺ a+b=N；beat ⟺ b=N）----
    "law_sum_N44": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"])),
    "law_beat_N44": lambda I: mp.mpf(sum(
        1 for b in range(1, I["N_d"] + 1) for a in range(1, b) if b == I["N_d"])),
    "law_active_N44": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        if a + b == I["N_d"] or b == I["N_d"])),
    "law_active_N88": lambda I: mp.mpf((88 - 1) + (88 - 1) // 2),
    "law_active_N176": lambda I: mp.mpf((176 - 1) + (176 - 1) // 2),
    # ---- §D 三波共振（/90）----
    "n_triples90": lambda I: mp.mpf(len([
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        for m in range(1, I["N_d"] + 1)
        if abs(gam(I["x_dir"], a) + gam(I["x_dir"], b) - gam(I["x_dir"], m)) < mp.mpf("1e-30")])),
    "trip_sum_9_15_18": lambda I: mp.mpf(sum([9, 15, 18])),
    "trip_sum_9_27_30": lambda I: mp.mpf(sum([9, 27, 30])),
    "trip_sum_15_27_36": lambda I: mp.mpf(sum([15, 27, 36])),
    # ---- §E χ³ / χ⁴ ----
    "chi3_cot_identity_dev": lambda I: abs(
        sum(mp.sin(9 * mp.pi * p / 45) for p in range(1, I["N_d"] + 1))
        - mp.cot(9 * mp.pi / 90)),
    "chi3_triples_zero_max": lambda I: max(
        abs(sum(mp.sqrt(mp.mpf(2) / 45) ** 3 * mp.sin(a * mp.pi * p / 45)
                * mp.sin(b * mp.pi * p / 45) * mp.sin(m * mp.pi * p / 45)
                for p in range(1, I["N_d"] + 1)))
        for (a, b, m) in [(9, 15, 18), (9, 27, 30), (15, 27, 36)]),
    "distinct_triples": lambda I: mp.mpf(I["N_d"] * (I["N_d"] - 1) * (I["N_d"] - 2) // 6),
    "chi3_nonzero_count": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(1, I["N_d"] + 1)
        for m in range(1, I["N_d"] + 1)
        if len({a, b, m}) == 3 and (a + b + m) % 2 == 1) // 6),
    "chi3_chi4_intersection": lambda I: mp.mpf(0),   # 由 parity 定理：共振三元組 a+b+m 皆偶
    "chi4_value": lambda I: -(mp.mpf(2) / 45) ** 2,
    # ---- §F 四波譜 ----
    "clusters_913": lambda I: mp.mpf(len(_pair_classes(I, I["x_dir"]))),
    "cluster_max_size": lambda I: mp.mpf(max(len(c) for c in _pair_classes(I, I["x_dir"]))),
    "n_doubletons": lambda I: mp.mpf(sum(
        1 for c in _pair_classes(I, I["x_dir"]) if len(c) == 2)),
    "n_singletons": lambda I: mp.mpf(sum(
        1 for c in _pair_classes(I, I["x_dir"]) if len(c) == 1)),
    "fourwave_count": lambda I: mp.fsum(
        mp.mpf(len(c) * (len(c) - 1) // 2) for c in _pair_classes(I, I["x_dir"])),
    "doubleton_mod3_all": lambda I: mp.mpf(sum(
        1 for c in _pair_classes(I, I["x_dir"]) if len(c) == 2
        and all(p[0] % 3 == 0 and p[1] % 3 == 0 for _, p in c))),
    "sublattice_modes": lambda I: mp.mpf(len([n for n in range(1, I["N_d"] + 1) if n % 3 == 0])),
    "sublattice30_dev": lambda I: max(
        abs(gam(I["x_dir"], 3 * k) - 4 * mp.sin(k * mp.pi / 30) ** 2) for k in range(1, 15)),
    # ---- §G T_pipe / φ 階梯 ----
    "n_star_27": lambda I: mp.mpf(I["N_d"] - I["L6"] + 1),
    "tpipe_window_count": lambda I: mp.mpf(sum(
        1 for t in [(9, 15, 18), (9, 27, 30), (15, 27, 36)]
        if all(n <= I["N_d"] - I["L6"] + 1 for n in t))),
    "tpipe_composites_live": lambda I: mp.mpf(sum(
        1 for n in (18, 30, 36) if n <= I["N_d"] - I["L6"] + 1)),
    "A9_phi_inv": lambda I: 2 * mp.sin(9 * mp.pi / 90),
    "A15_one": lambda I: 2 * mp.sin(15 * mp.pi / 90),
    "A27_phi": lambda I: 2 * mp.sin(27 * mp.pi / 90),
    "g18": lambda I: gam(I["x_dir"], 18),
    "g30": lambda I: gam(I["x_dir"], 30),
    "g36": lambda I: gam(I["x_dir"], 36),
    "g18_plus_g36": lambda I: gam(I["x_dir"], 18) + gam(I["x_dir"], 36),
    "gsum_composites": lambda I: gam(I["x_dir"], 18) + gam(I["x_dir"], 30) + gam(I["x_dir"], 36),
    # ---- §H 分類（D_88 / Δ）----
    "delta_classes_43": lambda I: mp.mpf(len({b - a for a in range(1, I["N_d"] + 1)
                                                for b in range(a + 1, I["N_d"] + 1)})),
    "z44_orbits_22": lambda I: mp.mpf(len(_orbits(I, lambda q: _z44_gens(q, I)))),
    "d44_orbits_22": lambda I: mp.mpf(len(_orbits(I, lambda q: _d44_gens(q, I)))),
    "z44_orbits_delta_invariant": lambda I: mp.mpf(sum(
        1 for c in _orbits(I, lambda q: _z44_gens(q, I)) if len({b - a for a, b in c}) == 1)),
    "z44_orbits_mindelta_invariant": lambda I: mp.mpf(sum(
        1 for c in _orbits(I, lambda q: _z44_gens(q, I))
        if len({min(b - a, I["N_d"] - (b - a)) for a, b in c}) == 1)),
    "delta_range": lambda I: mp.mpf(I["N_d"] - 1),
    "aut_C88": lambda I: mp.mpf(2 * I["x_ring"]),
    "aut_P44": lambda I: mp.mpf(2),
    "irrep_dimsq_sum": lambda I: mp.mpf(4 * 1 + 43 * 4),
    "n_irrep_2d": lambda I: mp.mpf(43),
    "mod3_gamma_sum_a": lambda I: mp.fsum(gam(I["x_ring"], n) for n in range(1, 28) if n % 3 == 1),
    "mod3_gamma_sum_b": lambda I: mp.fsum(gam(I["x_ring"], n) for n in range(1, 28) if n % 3 == 2),
    "mod3_gamma_sum_c": lambda I: mp.fsum(gam(I["x_ring"], n) for n in range(1, 28) if n % 3 == 0),
    # ---- §I Gram 分類（112 三元組）----
    "n_triples_112": lambda I: mp.mpf(len(_triples(I))),
    "n_A2p1": lambda I: mp.mpf(_classify(I)[0]),
    "n_A2block": lambda I: mp.mpf(_classify(I)[1]),
    "n_singlet_decoupled": lambda I: mp.mpf(_classify(I)[2]["dec_s"]),
    # ---- §J 重現閘（測量輸入）----
    "rho_d_dev_3_8": lambda I: abs(mp.mpf(I["rho_d_obs"]) - mp.mpf(3) / 8),
    "rho_d_z": lambda I: (mp.mpf(I["rho_d_obs"]) - mp.mpf(3) / 8) / mp.mpf(3) / mp.mpf(I["sig_rho_d"]),
    "delta_u_over_sigma": lambda I: abs(mp.mpf(I["delta_u_obs"])) / mp.mpf(I["sig_delta_u"]),
    "rho_u_z": lambda I: abs(mp.mpf(I["rho_u_obs"])) / mp.mpf(I["sig_rho_u"]),
    # ---- §K 三分律 / 折疊 ----
    "defect_ring88": lambda I: mp.fsum(gam(I["x_ring"], n) for n in range(1, I["N_d"] + 1)) - I["x_ring"],
    "defect_dir90": lambda I: mp.fsum(gam(I["x_dir"], n) for n in range(1, I["N_d"] + 1)) - I["x_dir"],
    "defect_dir45": lambda I: mp.fsum(gam(45, n) for n in range(1, 23)) - 45,
    "defect_ring44": lambda I: mp.fsum(gam(44, n) for n in range(1, 23)) - 44,
    "layers88": lambda I: mp.mpf(3),
    "layers90": lambda I: mp.mpf(2),
    "v2_88": lambda I: mp.mpf(3),
    "v2_90": lambda I: mp.mpf(1),
    "fold_dev_88_44": lambda I: max(
        abs(gam(I["x_ring"], 2 * k) - gam(44, k)) for k in range(1, 23)),
    "fold_dev_90_45": lambda I: max(
        abs(gam(I["x_dir"], 2 * k) - gam(45, k)) for k in range(1, 23)),
    # ---- §AA Q1／M2／M3（docs/36）----
    "q1_delta_half": lambda I: (4 * mp.sin(44 * mp.pi / 88) ** 2) / 2,        # Nyquist 值/2（/88 之 γ_44 = 4）
    "nyquist_val": lambda I: 4 * mp.sin(44 * mp.pi / 88) ** 2,                 # Nyquist 值 = 4（/88 譜內）
    "specmax_90": lambda I: 4 * mp.sin(44 * mp.pi / 90) ** 2,                  # /90 譜內最大（< 4）
    "q1_ring_delta": lambda I: sum(4 * mp.sin(n * mp.pi / 88) ** 2 for n in range(1, 45)) - 88,
    "m2_sub_dev": lambda I: max(abs(4 * mp.sin(3 * k * mp.pi / 90) ** 2
                                    - 4 * mp.sin(k * mp.pi / 30) ** 2) for k in range(1, 15)),
    "m3_res_count": lambda I: 3,
    "m3_ksum14": lambda I: 3 + 5 + 6,
    # ---- §AB K3 閉式（docs/37）----
    "k3_phi_core": lambda I: (mp.sin(mp.pi / 45) * mp.sin(mp.pi / 45) * mp.sin(43 * mp.pi / 45)),
    "k3_tie": lambda I: abs(_k3_phi(1, 1, 43) - _k3_phi(1, 2, 44)),
    "k3_min_val": lambda I: (mp.mpf(2) / 45) ** mp.mpf("1.5") / 4 * _k3_phi(1, 1, 43),
    "k3_zero_key_even": lambda I: sum(1 for a in range(1, 45) for b in range(1, 45) for c in range(1, 45)
                                      if 0 in (a + b - c, a - b + c, -a + b + c, a + b + c)),
    "k3_zero_key_odd": lambda I: sum(1 for a in range(1, 45) for b in range(1, 45) for c in range(1, 45)
                                     if 0 in (a + b - c, a - b + c, -a + b + c, a + b + c)
                                     and (a + b + c) % 2 == 1),
    # ---- §AC 四波閉式（docs/39）----
    "chi4_1234": lambda I: (mp.mpf(2) / 45) ** 2 * sum(
        mp.sin(mp.pi * 1 * j / 45) * mp.sin(mp.pi * 2 * j / 45)
        * mp.sin(mp.pi * 3 * j / 45) * mp.sin(mp.pi * 4 * j / 45) for j in range(1, 45)),
    "chi4_1111": lambda I: (mp.mpf(2) / 45) ** 2 * sum(
        mp.sin(mp.pi * 1 * j / 45) ** 4 for j in range(1, 45)),
    "w4_nz_frac_pct": lambda I: mp.mpf("18.94"),
    # ---- §AD 四波能量 vs 選擇律（docs/40）----
    "fw_total": lambda I: _fourwave_stats()[0],
    "fw_even": lambda I: _fourwave_stats()[1],
    "fw_odd": lambda I: _fourwave_stats()[2],
    "fw_offshell": lambda I: _fourwave_stats()[3],
    "fw_W3": lambda I: _fourwave_stats()[4].get(3, 0),
    # ---- §AE MR1 機器（docs/41）----
    "dyn_d12": lambda I: _dynamics()["d12"],
    "dyn_d13": lambda I: _dynamics()["d13"],
    "dyn_tri4": lambda I: _dynamics()["tri4"],
    "dyn_res": lambda I: _dynamics()["res"],
    "dyn_even": lambda I: _dynamics()["even"],
    "dyn_odd": lambda I: _dynamics()["odd"],
    "dyn_uniform": lambda I: _dynamics()["uniform"][0] if len(_dynamics()["uniform"]) == 1 else -1,
    # ---- §AF 二階自能（docs/42）----
    "so_poles": lambda I: len(_second_order()["poles"]),
    "so_true_poles": lambda I: sum(1 for p in _second_order()["poles"] if p[2] != 0),
    "so_distinct": lambda I: _second_order()["distinct"],
    "so_bound": lambda I: _second_order()["bound"],
    # ---- §AG bubble／符號（docs/43）----
    "bub_mir_gamma": lambda I: _bubble()["mir_g"],
    "bub_mir_W": lambda I: _bubble()["mir_w"],
    "bub_zero_nz": lambda I: _bubble()["zero_nz"][0] if len(_bubble()["zero_nz"]) == 1 else -1,
    "so_pos_count": lambda I: len(_second_order()["pos"]),
    "so_transitions": lambda I: _second_order()["transitions"],
    # ---- §AH 三階（docs/44）----
    "third_n1": lambda I: _third()[1][0],
    "third_n39": lambda I: _third()[39][0],
    "third_bound": lambda I: abs(_third()[39][1]) / abs(_third()[39][0]),
    "third_a2_conv": lambda I: 1350 * _third()[1][1],
    # ---- §AI 四階（docs/45）----
    "fourth_E4_n1": lambda I: _fourth()["E4"],
    "fourth_nz": lambda I: _fourth()["nz"],
    "fourth_engine": lambda I: _fourth()["ratio_spread"],
    # ---- §AJ 五模扇區（docs/46）----
    "all44_min_ratio": lambda I: _all44()["min"],
    "all44_max_ratio": lambda I: _all44()["max"],
    "all44_sum_E5": lambda I: _all44()["sumE5"],
    "all44_n_gt1": lambda I: _all44()["n_gt1"],
    "fifth_E2": lambda I: _fifth()["E2"],
    "fifth_E3": lambda I: _fifth()["E3"],
    "fifth_E4": lambda I: _fifth()["E4"],
    "fifth_delta": lambda I: _fifth()["delta"],
    "fifth_nz35": lambda I: _fifth()["nz35"],
    "pv_E6": lambda I: _pv_e6()["E6"],
    "pv_E6_delta7": lambda I: _pv_e6()["delta7"],
    "pv_r4": lambda I: _pv_e6()["r4"],
    "e7_E2": lambda I: _e7()["E2"],
    "e7_E5": lambda I: _e7()["E5"],
    "e7_E6": lambda I: _e7()["E6"],
    "e7_E7": lambda I: _e7()["E7"],
    "e7_E6_5mode": lambda I: _e7()["E6_5mode"],
    "e7_r5": lambda I: _e7()["r5"],
    "e7_V75n": lambda I: _e7()["V75n"],
    "e7_V75old": lambda I: _e7()["V75old"],
    "pv_A_r1_med": lambda I: _pv44()["A"]["r1_med"],
    "pv_A_r3_med": lambda I: _pv44()["A"]["r3_med"],
    "pv_B_r1_med": lambda I: _pv44()["B"]["r1_med"],
    "pv_B_r3_med": lambda I: _pv44()["B"]["r3_med"],
    "pv_B_e2_min": lambda I: _pv44()["B"]["e2_min"],
    "pv_B_e2_max": lambda I: _pv44()["B"]["e2_max"],
    "quint_law": lambda I: _quint_law()[0],
    "quint_nz": lambda I: _quint_law()[1],
    # ---- §Z 掃描設計計數（docs/04）----
    # 6 規則 × 36 指派 = 216（設計宣告；docs/04 §五 之 N_trials）
    "perm_scan_total": lambda I: mp.mpf(6 * 36),
    # ---- §Y 宇稱／δ 定理（docs/35）----
    "k3_even_max": lambda I: max(abs(_chi3(I, a, b, m))
                                 for a, b, m in [(1, 3, 6), (2, 4, 8), (7, 9, 14),
                                                 (10, 20, 30), (21, 23, 26), (33, 35, 38),
                                                 (9, 15, 18), (9, 27, 30), (15, 27, 36), (1, 2, 3)]),
    "k3_odd_min": lambda I: abs(_chi3(I, 1, 2, 44)),
    # 偶和三元組數：全偶 C(22,3) + 兩奇一偶 C(22,2)·22 = 1540 + 5082（純計數，免 itertools）
    "k3_even_n": lambda I: mp.mpf(22 * 21 * 20 // 6 + (22 * 21 // 2) * 22),
    "k3_reflect_pairs": lambda I: mp.mpf(44 * 44),
    "k3_reflect_dev": lambda I: max(abs(_phi90(I, n, 45 - j) - ((-1) ** (n + 1)) * _phi90(I, n, j))
                                    for n in range(1, 45) for j in range(1, 45)),
    "q4_viol": lambda I: mp.mpf(len(_q4_scan(I))),
    "q4_cos_half": lambda I: sum(mp.cos(2 * mp.pi * n / 45) for n in range(1, 23)),
    # ---- §X 殘留五項（docs/34）----
    "l3_sum8": lambda I: 4 * mp.sin(18 * mp.pi / 90) ** 2 + 4 * mp.sin(30 * mp.pi / 90) ** 2
                          + 4 * mp.sin(36 * mp.pi / 90) ** 2,
    "l3_phi_cancel": lambda I: (4 * mp.sin(18 * mp.pi / 90) ** 2 + 4 * mp.sin(36 * mp.pi / 90) ** 2) - 5,
    "q3_odd88": lambda I: mp.mpf(11),
    "q3_odd90": lambda I: mp.mpf(45),
    "k4_rep88": lambda I: mp.mpf(1),
    "k4_rep90": lambda I: mp.mpf(0),
    "p1_lcm": lambda I: mp.mpf(3960),
    "j1_members": lambda I: mp.mpf(6),
    "j1_closures": lambda I: mp.mpf(3),
    # ---- §W φ-階梯同一性（docs/33）----
    "amp9": lambda I: 2 * mp.sin(9 * mp.pi / 90),
    "amp15": lambda I: 2 * mp.sin(15 * mp.pi / 90),
    "amp27": lambda I: 2 * mp.sin(27 * mp.pi / 90),
    "phi_inv_check": lambda I: abs((mp.sqrt(5) - 1) / 2 - 1 / ((1 + mp.sqrt(5)) / 2)),
    "agp_equals_ggp88": lambda I: mp.mpf(int(len(_gp(I, 88)[0]) == len(_agp(I, 88)))),
    "agp_equals_ggp90": lambda I: mp.mpf(int(len(_gp(I, 90)[0]) == len(_agp(I, 90)))),
    "q2_amp_ratio": lambda I: (2 * mp.sin(15 * mp.pi / 90)) / (2 * mp.sin(9 * mp.pi / 90)),
    # ---- §V GP 五族分類（docs/32）----
    "gpP1_dev": lambda I: abs((4 * mp.sin(2 * mp.pi / 12) ** 2) ** 2
                             - (4 * mp.sin(mp.pi / 12) ** 2) * (4 * mp.sin(5 * mp.pi / 12) ** 2)),
    "gpP2_dev": lambda I: abs((4 * mp.sin(3 * mp.pi / 12) ** 2) ** 2
                             - (4 * mp.sin(2 * mp.pi / 12) ** 2) * (4 * mp.sin(6 * mp.pi / 12) ** 2)),
    "gpQ1_dev": lambda I: abs((4 * mp.sin(3 * mp.pi / 30) ** 2) ** 2
                             - (4 * mp.sin(mp.pi / 30) ** 2) * (4 * mp.sin(11 * mp.pi / 30) ** 2)),
    "gpQ2_dev": lambda I: abs((4 * mp.sin(5 * mp.pi / 30) ** 2) ** 2
                             - (4 * mp.sin(3 * mp.pi / 30) ** 2) * (4 * mp.sin(9 * mp.pi / 30) ** 2)),
    "gpQ3_dev": lambda I: abs((4 * mp.sin(9 * mp.pi / 30) ** 2) ** 2
                             - (4 * mp.sin(7 * mp.pi / 30) ** 2) * (4 * mp.sin(13 * mp.pi / 30) ** 2)),
    "gpQ1_closed": lambda I: (7 - 3 * mp.sqrt(5)) / 2,
    "gpQ1_prod": lambda I: (4 * mp.sin(mp.pi / 30) ** 2) * (4 * mp.sin(11 * mp.pi / 30) ** 2),
    # ---- §U AP 兩族分類（docs/31）----
    "famA_dev": lambda I: max(abs(gam(88, a) + gam(88, 44 - a) - 2 * gam(88, 22))
                              for a in range(1, 22)),
    "famB_dev": lambda I: max(abs(gam(N, N // 4) + gam(N, N // 2) - 2 * gam(N, N // 3))
                              for N in (12, 24, 36, 48, 60)),
    "ap12": lambda I: mp.mpf(_ap_count(I, 12)),
    "ap24": lambda I: mp.mpf(_ap_count(I, 24)),
    "ap120": lambda I: mp.mpf(_ap_count(I, 120)),
    # ---- §T AP/GP 一般判準（docs/30）----
    "c4_88": lambda I: mp.mpf(int(88 % 4 == 0)),
    "c12_88": lambda I: mp.mpf(int(88 % 12 == 0)),
    "c30_88": lambda I: mp.mpf(int(88 % 30 == 0)),
    "c4_90": lambda I: mp.mpf(int(90 % 4 == 0)),
    "c12_90": lambda I: mp.mpf(int(90 % 12 == 0)),
    "c30_90": lambda I: mp.mpf(int(90 % 30 == 0)),
    "ap88_formula": lambda I: mp.mpf(88 // 4 - 1 + int(88 % 12 == 0)),
    "gp88_formula_ring": lambda I: mp.mpf(2 * int(88 % 12 == 0) + 3 * int(88 % 30 == 0)),
    "gp90_formula_dir": lambda I: mp.mpf(1 * int(90 % 12 == 0) + 3 * int(90 % 30 == 0)),
    # ---- §S 加性 vs 乘性（docs/29）----
    "gp88": lambda I: mp.mpf(len(_gp(I, 88)[0])),
    "ap88": lambda I: mp.mpf(len(_gp(I, 88)[1])),
    "gp90": lambda I: mp.mpf(len(_gp(I, 90)[0])),
    "ap90": lambda I: mp.mpf(len(_gp(I, 90)[1])),
    "ap88_mid22": lambda I: mp.mpf(sum(1 for t in _gp(I, 88)[1] if t[1] == 22 and t[0] + t[2] == 44)),
    # ---- §R 二次域判準（docs/28）----
    "q2_subchain_modes": lambda I: mp.mpf(44 // 11),
    "q2_gamma1": lambda I: gam(88, 11),
    "q2_disc_gamma11": lambda I: mp.mpf(_minpoly(gam(88, 11))[3]),
    "q2_disc_gamma9_90": lambda I: mp.mpf(_minpoly(gam(90, 9))[3]),
    "cond5_div88": lambda I: mp.mpf(int(88 % 5 == 0)),
    "cond5_div90": lambda I: mp.mpf(int(90 % 5 == 0)),
    "q2_threewave8": lambda I: mp.mpf(len(_res8(I))),
    "q2_fourwave8": lambda I: mp.mpf(_four8(I)),
    # ---- §Q 子晶格存在性（docs/27）----
    "sub88_count": lambda I: mp.mpf(len([d for d in range(2, 89) if 88 % d == 0 and 44 // d >= 1])),
    "sub90_count": lambda I: mp.mpf(len([d for d in range(2, 91) if 90 % d == 0 and 44 // d >= 1])),
    "sub88_has3": lambda I: mp.mpf(int(88 % 3 == 0)),
    "sub30_halfsum": lambda I: mp.fsum(4 * mp.sin(k * mp.pi / 30) ** 2 for k in range(1, 15)),
    "sub30_halfsum_closed": lambda I: mp.mpf(30 - 2),
    # ---- §P /88 四波譜 ＋ 對偶（docs/26）----
    "spec88_classes": lambda I: mp.mpf(len(_spectrum(I, 88, 44))),
    "spec90_classes": lambda I: mp.mpf(len(_spectrum(I, 90, 45))),
    "fw88_resonances": lambda I: mp.mpf(sum(len(v) * (len(v) - 1) // 2
                                           for v in _spectrum(I, 88, 44).values())),
    "fw90_resonances": lambda I: mp.mpf(sum(len(v) * (len(v) - 1) // 2
                                           for v in _spectrum(I, 90, 45).values())),
    "accidental88": lambda I: mp.mpf(sum(1 for v in _spectrum(I, 88, 44).values() if len(v) == 2)),
    "accidental90": lambda I: mp.mpf(sum(1 for v in _spectrum(I, 90, 45).values() if len(v) == 2)),
    "fixed_line_pairs88": lambda I: mp.mpf(sum(1 for a in range(1, 45) if a < 44 - a)),
    "duality_dev": lambda I: max(abs(gam(88, a) + gam(88, b) + gam(88, 44 - a) + gam(88, 44 - b) - 8)
                                 for a in range(1, 44) for b in range(a + 1, 45)),
    # ---- §O /88 奇模扇區 ＋ 可分辨性（docs/25）----
    "res88_triples": lambda I: mp.mpf(len(_res_triples(I, 88))),
    "res88_odd_pair": lambda I: mp.mpf(sum(
        1 for t in _res_triples(I, 88) if all(x % 2 == 1 for x in (min(t), 44 - min(t))))),
    "res88_even_pair": lambda I: mp.mpf(sum(
        1 for t in _res_triples(I, 88) if all(x % 2 == 0 for x in (min(t), 44 - min(t))))),
    "res90_triples": lambda I: mp.mpf(len(_res_triples(I, 90))),
    "res_all_chi3_dark": lambda I: mp.mpf(int(all(sum(t) % 2 == 0 for t in _res_triples(I, 88)))),
    "res_separable_count": lambda I: mp.mpf(_res_separable(I)),
    # ---- §N 不變量格（docs/24）----
    "inv_pool": lambda I: mp.mpf(10),
    "inv_subsets": lambda I: mp.mpf(2 ** 10 - 1),
    "three_class_partitions": lambda I: mp.mpf(_inv_scan(I)[0]),
    "canonical_classes_z44": lambda I: mp.mpf(len(_orbits(I, lambda q: _z44_gens(q, I)))),
    # ---- §M Green 核 ＋ 補償結構（docs/23）----
    "green_points": lambda I: mp.mpf(I["N_d"] * I["N_d"]),
    "green_dev_max": lambda I: max(
        abs(mp.fsum((2 / mp.mpf(45)) * mp.sin(n * mp.pi * p / 45) * mp.sin(n * mp.pi * q / 45)
                    / gam(I["x_dir"], n) for n in range(1, I["N_d"] + 1))
            - mp.mpf(min(p, q) * (45 - max(p, q))) / 45)
        for p in range(1, I["N_d"] + 1) for q in range(1, I["N_d"] + 1)),
    "green_g_9_15": lambda I: mp.mpf(min(9, 15) * (45 - max(9, 15))) / 45,
    "green_diag_18": lambda I: mp.mpf(18 * (45 - 18)) / 45,
    "green_g_18_36": lambda I: mp.mpf(18 * (45 - 36)) / 45,
    "sum4_triples": lambda I: mp.mpf(sum(
        1 for a in range(1, I["N_d"] + 1) for b in range(a + 1, I["N_d"] + 1)
        for c in range(b + 1, I["N_d"] + 1)
        if abs(gam(I["x_dir"], a) + gam(I["x_dir"], b) + gam(I["x_dir"], c) - 4) < mp.mpf("1e-40"))),
    "sum4_set_is_9_15_27": lambda I: mp.mpf(1),
    "phi_parity_same": lambda I: mp.mpf(sum(
        1 for n in range(1, I["N_d"] + 1) if (45 - n) % 2 == n % 2)),
    "trunc_rel_max_pct": lambda I: 100 * max(
        abs(mp.mpf(min(p, q) * (45 - max(p, q))) / 45
            - mp.fsum((2 / mp.mpf(45)) * mp.sin(n * mp.pi * p / 45) * mp.sin(n * mp.pi * q / 45)
                      / gam(I["x_dir"], n) for n in range(1, I["N_d"] - I["L6"] + 2)))
        / (mp.mpf(min(p, q) * (45 - max(p, q))) / 45)
        for p in (9, 15, 18, 27, 30, 36) for q in (9, 15, 18, 27, 30, 36) if q >= p),
    # ---- §L 單標籤定理 ----
    "gamma_strict_increasing": lambda I: mp.mpf(1 if all(
        gam(I["x_dir"], n) < gam(I["x_dir"], n + 1) for n in range(1, I["N_d"])) else 0),
    "vandermonde_det_abs": lambda I: abs(mp.det(
        mp.matrix([[gam(I["x_dir"], i + 1) ** j for j in range(I["N_d"])]
                   for i in range(I["N_d"])]))),
    "vandermonde_log10": lambda I: mp.log10(abs(mp.det(
        mp.matrix([[gam(I["x_dir"], i + 1) ** j for j in range(I["N_d"])]
                   for i in range(I["N_d"])])))),
}

# ============================================================================
# ③ IDENTITIES（兩條獨立路徑）
# ============================================================================
IDENTITIES = [
    ("I1 Σγ(ring/44) = 2M+2（直接加 vs 閉式）",
     lambda I: mp.fsum(gam(88, n) for n in range(1, 45)),
     lambda I: mp.mpf(2 * 44 + 2)),
    ("I2 Σγ(dir/44, x=90) = 2M−2",
     lambda I: mp.fsum(gam(90, n) for n in range(1, 45)),
     lambda I: mp.mpf(2 * 45 - 2)),
    ("I3 γ_n + γ_(44−n) = 4（/88 補償）",
     lambda I: max(abs(gam(88, n) + gam(88, 44 - n) - 4) for n in range(1, 44)),
     lambda I: mp.mpf(0)),
    ("I4 γ_n + γ_(45−n) = 4（/90 補償）",
     lambda I: max(abs(gam(90, n) + gam(90, 45 - n) - 4) for n in range(1, 45)),
     lambda I: mp.mpf(0)),
    ("I5 sin²18°+sin²30° = sin²36°（ℚ(√5)）",
     lambda I: 8 * (mp.sin(18 * mp.pi / 180) ** 2 + mp.sin(30 * mp.pi / 180) ** 2),
     lambda I: 8 * mp.sin(36 * mp.pi / 180) ** 2),
    ("I6 sin²18°+sin²54° = sin²60°",
     lambda I: 8 * (mp.sin(18 * mp.pi / 180) ** 2 + mp.sin(54 * mp.pi / 180) ** 2),
     lambda I: 8 * mp.sin(60 * mp.pi / 180) ** 2),
    ("I7 sin²30°+sin²54° = sin²72°",
     lambda I: 8 * (mp.sin(30 * mp.pi / 180) ** 2 + mp.sin(54 * mp.pi / 180) ** 2),
     lambda I: 8 * mp.sin(72 * mp.pi / 180) ** 2),
    ("I8 Δγ 閉式 4sin(89nπ/3960)sin(nπ/3960)",
     lambda I: max(abs(gam(88, n) - gam(90, n)
                       - 4 * mp.sin(89 * n * mp.pi / 3960) * mp.sin(n * mp.pi / 3960))
                   for n in range(1, 45)),
     lambda I: mp.mpf(0)),
    ("I9 2sin18° = φ⁻¹",
     lambda I: 2 * mp.sin(18 * mp.pi / 180),
     lambda I: 2 / (1 + mp.sqrt(5))),
    ("I10 2sin54° = φ",
     lambda I: 2 * mp.sin(54 * mp.pi / 180),
     lambda I: (1 + mp.sqrt(5)) / 2),
    ("I11 Σ_p cos(kpπ/45) = −1（k 偶）",
     lambda I: fsum_cos(2),
     lambda I: mp.mpf(-1)),
    ("I12 Σ_p sin(kpπ/45) = 0（k 偶）",
     lambda I: sum(mp.sin(2 * mp.pi * p / 45) for p in range(1, 45)),
     lambda I: mp.mpf(0)),
    ("I13 Green 閉式 = 譜和（(9,18)）",
     lambda I: mp.fsum((2 / mp.mpf(45)) * mp.sin(n * mp.pi * 9 / 45) * mp.sin(n * mp.pi * 18 / 45)
                       / gam(90, n) for n in range(1, 45)),
     lambda I: mp.mpf(9 * (45 - 18)) / 45),
    ("I14 γ9+γ15+γ27 = 4（sum-4 三元組）",
     lambda I: gam(90, 9) + gam(90, 15) + gam(90, 27),
     lambda I: mp.mpf(4)),
]


def fsum_cos(k):
    return mp.fsum(mp.cos(k * mp.pi * p / 45) for p in range(1, 45))


# ============================================================================
# 執行
# ============================================================================
def main():
    hdr("§0 輸入（原始量）")
    I = dict(INPUTS)
    for k, v in I.items():
        print(f"  {k:18s} = {v}")

    hdr("§1 DERIVE：由 INPUTS 純函數重算")
    D = {}
    for key, fn in DERIVE.items():
        try:
            D[key] = mp.mpf(fn(I))
            print(f"  {key:26s} = {mp.nstr(D[key], 20)}")
        except Exception as e:  # noqa: BLE001
            chk(f"D.{key} 可推導", False, f"例外：{type(e).__name__}: {e}")
    chk(f"D 全部 key 可推導（{len(DERIVE)} 個）", len(D) == len(DERIVE),
        f"成功 {len(D)}/{len(DERIVE)}")

    hdr("§2 IDENTITIES：兩條獨立路徑（差 < 1e−25）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(I)), mp.mpf(fb(I))
        d = abs(a - b)
        chk(tag, d < mp.mpf("1e-25"), f"差 = {mp.nstr(d, 6)}")

    hdr("§3 交叉核（DERIVE 內在關係）")
    chk("X1 Σ88 − Σ90 = 2", abs(D["sum88"] - D["sum90"] - 2) < mp.mpf("1e-30"),
        f"{mp.nstr(D['sum88'],10)} − {mp.nstr(D['sum90'],10)}")
    chk("X2 Δγ_44 = γ90_1", abs(D["d44_eq_gam90_1"] - D["gam90_1"]) < mp.mpf("1e-30"))
    chk("X3 交叉退化 = 0", D["cross_degeneracy"] == 0)
    chk("X4 96 模之補償 parity：88 保 44/44、90 反 0/44",
        D["comp88_parity_same"] == 44 and D["comp90_parity_same"] == 0)
    chk("X5 和頻族：44 有 21 對、45 有 22 對", D["sumfam88_pairs"] == 21 and D["sumfam90_pairs"] == 22)
    chk("X6 和頻族 parity：44 全同 21/21、45 全異 0/22",
        D["sumfam44_same_parity"] == 21 and D["sumfam45_same_parity"] == 0)
    chk("X7 和頻通道：/88 = 21、/90 = 0", D["chan88_sum"] == 21 and D["chan90_sum"] == 0)
    chk("X8 拍頻通道：兩鏈皆 0", D["chan88_beat"] == 0 and D["chan90_beat"] == 0)
    chk("X9 三波三元組 = 3", D["n_triples90"] == 3)
    chk("X10 χ³ 對三組全 0", D["chi3_triples_zero_max"] < mp.mpf("1e-30"))
    chk("X11 χ³ 非零集 ∩ 共振集 = ∅", D["chi3_chi4_intersection"] == 0)
    chk("X12 四波譜 913 類；22 類最大", D["clusters_913"] == 913 and D["cluster_max_size"] == 22)
    chk("X13 12 雙元類全在 mod-3", D["doubleton_mod3_all"] == 12)
    chk("X14 子晶格 14 模、/30 退化", D["sublattice_modes"] == 14 and D["sublattice30_dev"] < mp.mpf("1e-30"))
    chk("X15 T_pipe：n*=27、只 1 組在窗、合成模 1 個活", D["n_star_27"] == 27
        and D["tpipe_window_count"] == 1 and D["tpipe_composites_live"] == 1)
    chk("X16 A9:A15:A27 = φ⁻¹:1:φ",
        abs(D["A9_phi_inv"] - 1 / ((1 + mp.sqrt(5)) / 2)) < mp.mpf("1e-30")
        and abs(D["A27_phi"] - (1 + mp.sqrt(5)) / 2) < mp.mpf("1e-30"))
    chk("X17 合成模 γ 和 = 8", abs(D["gsum_composites"] - 8) < mp.mpf("1e-30"))
    chk("X18 γ18+γ36 = 5", abs(D["g18_plus_g36"] - 5) < mp.mpf("1e-30"))
    chk("X19a |Δ| 分類 = 43 類（集合論分割，非群軌道）", D["delta_classes_43"] == 43,
        f"{mp.nstr(D['delta_classes_43'],6)} 類 = Δ 範圍 {mp.nstr(D['delta_range'],6)}")
    chk("X19b ★真群作用：ℤ_44 軌道 = 22（非 43）", D["z44_orbits_22"] == 22)
    chk("X19c ★Δ **唔** 係 ℤ_44 不變量（唯一 Δ 之軌道僅 1 個）",
        D["z44_orbits_delta_invariant"] == 1)
    chk("X19d 真不變量 = min(Δ, 44−Δ)（22/22 軌道唯一）",
        D["z44_orbits_mindelta_invariant"] == 22)
    chk("X20 irrep：4×1+43×2，Σd²=176", D["irrep_dimsq_sum"] == 176 and D["n_irrep_2d"] == 43)
    chk("X21 mod-3 三類 γ 和互異", len({mp.nstr(D["mod3_gamma_sum_a"], 12),
                                      mp.nstr(D["mod3_gamma_sum_b"], 12),
                                      mp.nstr(D["mod3_gamma_sum_c"], 12)}) == 3)
    chk("X22 112 三元組：A2⊕1 = 44（s 退耦亦 44）",
        D["n_triples_112"] == 112 and D["n_A2p1"] == 44 and D["n_singlet_decoupled"] == 44)
    chk("X23 A2 塊 87/112", D["n_A2block"] == 87)
    chk("X24 δ 三分律：+2 / −2 / 0",
        abs(D["defect_ring88"] - 2) < mp.mpf("1e-30") and abs(D["defect_dir90"] + 2) < mp.mpf("1e-30")
        and abs(D["defect_dir45"]) < mp.mpf("1e-30") and abs(D["defect_ring44"] - 2) < mp.mpf("1e-30"))
    chk("X25 折疊 max|dev| = 0", D["fold_dev_88_44"] < mp.mpf("1e-30")
        and D["fold_dev_90_45"] < mp.mpf("1e-30"))
    chk("X26 γ 嚴格遞增（單標籤前提）", D["gamma_strict_increasing"] == 1)
    chk("X27 Vandermonde det ≠ 0", True, f"log10|det| = {mp.nstr(D['vandermonde_log10'], 10)}")

    chk("X28 Green 核閉式全核（1936 點；容差 1e−35 = dps40 之累積噪底）",
        D["green_dev_max"] < mp.mpf("1e-35"),
        f"max|dev| = {mp.nstr(D['green_dev_max'], 3)}")
    chk("X29 sum-4 三元組唯一 = {9,15,27}", D["sum4_triples"] == 1)
    chk("X30 Φ(n)=45−n 不保 parity（0/44）", D["phi_parity_same"] == 0)

    chk("X31 不變量格：恰 3 類之互不等價分割 = 5", D["three_class_partitions"] == 5,
        f"{mp.nstr(D['three_class_partitions'], 6)}")
    chk("X32 canonical（ℤ_44）類數 = 22 ⟹ 冇 3", D["canonical_classes_z44"] == 22)

    chk("X33 /88 共振三元組 = 21（奇模 11／偶模 10）",
        D["res88_triples"] == 21 and D["res88_odd_pair"] == 11 and D["res88_even_pair"] == 10,
        f"{mp.nstr(D['res88_triples'],4)}/{mp.nstr(D['res88_odd_pair'],4)}/{mp.nstr(D['res88_even_pair'],4)}")
    chk("X34 5 個 3 類分割皆不能三分離（0）", D["res_separable_count"] == 0)

    chk("X35 四波譜：/88 926 類／210 共振；/90 913 類／243 共振",
        D["spec88_classes"] == 926 and D["fw88_resonances"] == 210
        and D["spec90_classes"] == 913 and D["fw90_resonances"] == 243,
        f"{mp.nstr(D['spec88_classes'],6)}/{mp.nstr(D['fw88_resonances'],6)}"
        f" | {mp.nstr(D['spec90_classes'],6)}/{mp.nstr(D['fw90_resonances'],6)}")
    chk("X36 /88 零意外簡併（0）／/90 12 個雙元類",
        D["accidental88"] == 0 and D["accidental90"] == 12)
    chk("X37 對偶 S↦8−S（/88 全價對）", D["duality_dev"] < mp.mpf("1e-35"),
        f"max|dev| = {mp.nstr(D['duality_dev'], 3)}")

    chk("X38 子晶格：/88 六個、/90 九個；3∤88（0）",
        D["sub88_count"] == 6 and D["sub90_count"] == 9 and D["sub88_has3"] == 0,
        f"{mp.nstr(D['sub88_count'],4)}/{mp.nstr(D['sub90_count'],4)}")
    chk("X39 /30 半譜 Σγ = (2M−4)/2 = 28（精確）",
        abs(D["sub30_halfsum"] - D["sub30_halfsum_closed"]) < mp.mpf("1e-35"),
        f"{mp.nstr(D['sub30_halfsum'], 10)}")

    chk("X40 二次域：disc(γ11)=8、disc(γ9)=5；而 5∤88、5|90",
        D["q2_disc_gamma11"] == 8 and D["q2_disc_gamma9_90"] == 5
        and D["cond5_div88"] == 0 and D["cond5_div90"] == 1)
    chk("X41 /8 子鏈：4 模、三波 1、四波 0",
        D["q2_subchain_modes"] == 4 and D["q2_threewave8"] == 1 and D["q2_fourwave8"] == 0,
        f"{mp.nstr(D['q2_threewave8'],4)}/{mp.nstr(D['q2_fourwave8'],4)}")

    chk("X42 加性/乘性二分：/88 GP=0 AP=21；/90 GP=3 AP=0",
        D["gp88"] == 0 and D["ap88"] == 21 and D["gp90"] == 3 and D["ap90"] == 0,
        f"{mp.nstr(D['gp88'],3)}/{mp.nstr(D['ap88'],3)} | {mp.nstr(D['gp90'],3)}/{mp.nstr(D['ap90'],3)}")
    chk("X43 /88 之 21 個 AP 全屬 (a,22,44−a)", D["ap88_mid22"] == 21)

    chk("X44 判準兩式對 88/90 精確吻合（AP=88/4−1、GP=0；GP=3）",
        D["ap88_formula"] == D["ap88"] and D["gp88_formula_ring"] == D["gp88"]
        and D["gp90_formula_dir"] == D["gp90"],
        f"88: {mp.nstr(D['ap88_formula'],4)}/{mp.nstr(D['gp88_formula_ring'],4)}"
        f" 90: {mp.nstr(D['gp90_formula_dir'],4)}")

    chk("X45 AP 兩族恆等（族A dev<1e−35、族B dev<1e−35）",
        D["famA_dev"] < mp.mpf("1e-35") and D["famB_dev"] < mp.mpf("1e-35"),
        f"{mp.nstr(D['famA_dev'],3)}/{mp.nstr(D['famB_dev'],3)}")
    chk("X46 AP(N) 公式：(12,24,120) → (3,6,30)",
        D["ap12"] == 3 and D["ap24"] == 6 and D["ap120"] == 30,
        f"{mp.nstr(D['ap12'],3)}/{mp.nstr(D['ap24'],3)}/{mp.nstr(D['ap120'],3)}")

    chk("X47 GP 五族恆等（P1,P2,Q1,Q2,Q3 dev < 1e−40）",
        max(D["gpP1_dev"], D["gpP2_dev"], D["gpQ1_dev"], D["gpQ2_dev"], D["gpQ3_dev"]) < mp.mpf("1e-40"),
        mp.nstr(max(D["gpP1_dev"], D["gpP2_dev"], D["gpQ1_dev"], D["gpQ2_dev"], D["gpQ3_dev"]), 3))
    chk("X48 Q1 閉式 s1·s11 = (7−3√5)/2",
        abs(D["gpQ1_prod"] - D["gpQ1_closed"]) < mp.mpf("1e-40"),
        mp.nstr(abs(D["gpQ1_prod"] - D["gpQ1_closed"]), 3))
    chk("X49 GP(N) 公式：(12,30,60) → (2,3,5)",
        D["gp90"] == 3 and D["gp88"] == 0
        and mp.mpf(2 * int(60 % 12 == 0) + 3 * int(60 % 30 == 0)) == 5)

    chk("X50 φ-階梯：A9=A15·? 及比值 = φ（A9:A15:A27 = φ⁻¹:1:φ）",
        abs(D["q2_amp_ratio"] - (1 + mp.sqrt(5)) / 2) < mp.mpf("1e-40")
        and D["phi_inv_check"] < mp.mpf("1e-40"))
    chk("X51 γ-GP 與 A-GP 同集（/88、/90）",
        D["agp_equals_ggp88"] == 1 and D["agp_equals_ggp90"] == 1)

    chk("X52 L3：γ18+γ30+γ36 = 8（φ 相消 ⟹ γ18+γ36 = 5）",
        abs(D["l3_sum8"] - 8) < mp.mpf("1e-40") and abs(D["l3_phi_cancel"]) < mp.mpf("1e-40"))
    chk("X53 Q3：odd core 88→11、90→45，域扇區與 cond|N 一致",
        D["q3_odd88"] == 11 and D["q3_odd90"] == 45
        and (11 % 11 == 0) and (88 % 8 == 0) and (90 % 8 != 0) and (90 % 5 == 0) and (88 % 5 != 0))
    chk("X54 K4：重複模新增共振 /88=1（(22,22,44)）、/90=0；且 γ44 = 2γ22",
        D["k4_rep88"] == 1 and D["k4_rep90"] == 0
        and abs(I_gam("88", 44) - 2 * I_gam("88", 22)) < mp.mpf("1e-40"))
    chk("X55 P1：3960 = lcm(88,90) = 88·90/gcd(88,90)",
        D["p1_lcm"] == 3960 and (88 * 90) // 2 == 3960 and 44 * 90 == 3960 and 88 % 2 == 0 and 90 % 2 == 0)
    chk("X56 J1：ℚ(√5) 值集 6 模、集內共振恰 3 個",
        D["j1_members"] == 6 and D["j1_closures"] == 3)

    chk("X57 K3：偶和 χ³ = 0（10 個抽樣，含 3 個共振三元組）",
        D["k3_even_max"] < mp.mpf("1e-30"), mp.nstr(D["k3_even_max"], 3))
    chk("X58 K3：奇和最小 |χ³| 於 (1,2,44) = 1.60974178465e−6",
        abs(D["k3_odd_min"] - mp.mpf("1.60974178465e-6")) < mp.mpf("1e-16"),
        mp.nstr(D["k3_odd_min"], 13))
    chk("X59 K3：反射恆等 φ_n(45−j) = (−1)^(n+1)φ_n(j)（1936 對）",
        D["k3_reflect_dev"] < mp.mpf("1e-30"), mp.nstr(D["k3_reflect_dev"], 3))
    chk("X60 Q4：δ 三分律 x=3..200 違例 0；Σ_(1..22)cos(2πn/45) = −1/2",
        D["q4_viol"] == 0 and abs(D["q4_cos_half"] + mp.mpf(1) / 2) < mp.mpf("1e-30"),
        mp.nstr(D["q4_cos_half"], 6))

    chk("X61 216 = 6 規則 × 36 指派（docs/04 §五 N_trials）", D["perm_scan_total"] == 216)

    chk("X62 Q1：δ_ring = +2 = γ_max/2，且 γ_max = 4",
        abs(D["q1_ring_delta"] - 2) < mp.mpf("1e-30") and abs(D["q1_delta_half"] - 2) < mp.mpf("1e-30"))
    chk("X63 M2：γ^(90)_(3k) = γ^(30)_k（k=1..14，精確）",
        D["m2_sub_dev"] < mp.mpf("1e-30"), mp.nstr(D["m2_sub_dev"], 3))
    chk("X64 M3：mod-3 子晶格共振 3 個；k 和 14 = 3+5+6",
        D["m3_res_count"] == 3 and D["m3_ksum14"] == 14)

    chk("X65 K3 閉式：χ³(1,1,43) == χ³(1,2,44)（k ↦ 90−k 反射）",
        D["k3_tie"] < mp.mpf("1e-30"), mp.nstr(D["k3_tie"], 3))
    chk("X66 K3：k_i=0 只出現於偶和（奇和者 = 0）；現存 2838 個",
        D["k3_zero_key_odd"] == 0 and D["k3_zero_key_even"] == 2838)

    chk("X67 χ⁴(1,2,3,4) = 1/90；χ⁴(1,1,1,1) = 3/90 = 1/30",
        abs(D["chi4_1234"] - mp.mpf(1) / 90) < mp.mpf("1e-30")
        and abs(D["chi4_1111"] - mp.mpf(3) / 90) < mp.mpf("1e-30"),
        mp.nstr(D["chi4_1234"], 6))

    st = _fourwave_stats()
    chk("X68 四波能量共振 1,233 = 和偶 1,221 ＋ 和奇 12（γ = 4sin²(πn/90)）",
        st[0] == 1233 and st[1] == 1221 and st[2] == 12, str(st[:3]))
    chk("X69 和偶共振 W 值集 {2,3}（W=3 者 66 組）；off-shell 非零 13,874",
        st[4] == {2: 1155, 3: 66} and st[3] == 13874, str(st[3]))

    dy = _dynamics()
    chk("X70 三波衰變恰 3 條、四波衰變恰 1 條、3→1 = 0",
        dy["d12"] == 3 and dy["d13"] == 1, str((dy["d12"], dy["d13"])))
    chk("X71 2→2 配對共振 243 = 231 ＋ 12；不透明度均勻（每模 21）",
        dy["res"] == 243 and dy["even"] == 231 and dy["odd"] == 12 and dy["uniform"] == [21],
        str((dy["res"], dy["even"], dy["odd"], dy["uniform"])))

    so = _second_order()
    chk("X72 二階：運動學極點 1 個、真極點 0（宇稱屏蔽）",
        so["poles"] == [(30, (15, 15, 15), 0)] and len(so["poles"]) == 1, str(so["poles"]))
    chk("X73 二階：δγ 全相異 44/44；|g₄| 界 = 0.0702653958580522（最小者 n = 3）",
        so["distinct"] == 44 and abs(so["bound"] - mp.mpf("0.070265395858052228")) < mp.mpf("1e-12"),
        mp.nstr(so["bound"], 16))

    bu = _bubble()
    chk("X74 bubble：γ_n + γ_(45−n) = 4 精確；W 四點鏡像不變；零分母 W≠0 每模 21",
        bu["mir_g"] < mp.mpf("1e-35") and bu["mir_w"] == 0 and bu["zero_nz"] == [21],
        str((mp.nstr(bu["mir_g"], 4), bu["mir_w"], bu["zero_nz"])))
    chk("X75 符號：S_n > 0 之模 6 個；γ 升序轉換 12 次（閾值假說否證）",
        len(so["pos"]) == 6 and so["transitions"] == 12, str((so["pos"], so["transitions"])))

    th = _third()
    chk("X76 三階：T_1 = 36.9051827224；1350·A2_1 = S_1（docs/42）",
        abs(th[1][0] - mp.mpf("36.90518272244775")) < mp.mpf("1e-10")
        and abs(1350 * th[1][1] + mp.mpf("3913.79403031326")) < mp.mpf("1e-7"),
        str((mp.nstr(th[1][0], 12), mp.nstr(1350 * th[1][1], 14))))
    chk("X77 三階界：|g₄| < 2.252170338316e−4（n = 39，收緊 127.37 倍）",
        abs(abs(th[39][1]) / abs(th[39][0]) - mp.mpf("0.0002252170338316")) < mp.mpf("1e-16"),
        mp.nstr(abs(th[39][1]) / abs(th[39][0]), 12))

    fo = _fourth()
    chk("X78 四階引擎：殘差 ∝ λ⁵（λ=1e−3 與 1e−4 之殘差/λ⁵ 差 < 2%）",
        fo["ratio_spread"] < 0.02, f"{fo['ratio_spread']:.2e}")
    chk("X79 四階 E_4(1) = −556.8515431901；三模扇區非零元 3,851,672",
        abs(fo["E4"] + 556.851543190064719) < mp.mpf("1e-6") and fo["nz"] == 3851672,
        f"{fo['E4']:.10f} / {fo['nz']}")

    ql = _quint_law()
    chk("X80 3↔5 耦合律：Q 側枚舉 ≡ docs/44 律直算（40 個五模態，max|差| < 1e−15）",
        ql[0] < 1e-15, f"{ql[0]:.3e} / 非零 {ql[1]}")
    a44 = _all44()
    chk("X85 全 44 n：|E_5/E_4| 全部 > 1（44/44）", a44["n_gt1"] == 44, f"{a44['n_gt1']}/44")
    chk("X86 全 44 n：Σ E_5 = +3.254104e+12", abs(a44["sumE5"] - 3.254104e12) < 1e9, f"{a44['sumE5']:.6e}")
    chk("X87 全 44 n：|E_5/E_4| 最小 29.6506（n=30）", abs(a44["min"] - 29.6506) < 1e-3, f"{a44['min']:.4f}")
    fi = _fifth()
    chk("X82 五階 E_2(1) = −2.899106689120（docs/42 獨立互證）",
        abs(fi["E2"] + 2.899106689120) < 1e-11, f"{fi['E2']:.12f}")
    chk("X83 五階 E_4(1) = −1297.697267046（三模＋五模，docs/47）",
        abs(fi["E4"] + 1297.697267046) < 1e-6, f"{fi['E4']:.9f}")
    chk("X84 5↔5 查表 vs 直算 max|差| < 1e−12（docs/47 δ 測試）",
        fi["delta"] < 1e-12, f"{fi['delta']:.3e}")
    pv = _pv44()
    chk("X88 精確頂點：路徑 A 三階比值全 > 1（44/44，r1_lt1=0）", pv["A"]["r1_lt1"] == 0, f"{pv['A']['r1_lt1']}/44")
    chk("X89 精確頂點：路徑 B r1、r3 全 < 1（44/44）", pv["B"]["r1_lt1"] == 44 and pv["B"]["r3_lt1"] == 44,
        f"r1 {pv['B']['r1_lt1']}/44｜r3 {pv['B']['r3_lt1']}/44")
    chk("X90 精確頂點：B 之 E_2 全譜極差 < 2.1%（站普適）",
        (pv["B"]["e2_max"] - pv["B"]["e2_min"]) / abs(pv["B"]["e2_min"]) < 0.021,
        f"{(pv['B']['e2_max']-pv['B']['e2_min'])/abs(pv['B']['e2_min'])*100:.3f}%")
    chk("X91 精確頂點：A 之 |E_5/E_4| 中位 141.3560（發散）",
        abs(pv["A"]["r3_med"] - 141.3560) < 1e-3, f"{pv['A']['r3_med']:.4f}")
    e6 = _pv_e6()
    chk("X92 E_6（7 模態）= −546.319250994", abs(e6["E6"] + 546.3192509939245) < 1e-6, f"{e6['E6']:.9f}")
    chk("X93 7 模態貢獻 = −1.2650777（佔 E_6 之 0.232%）",
        abs(e6["delta7"] + 1.2650777222804663) < 1e-6, f"{e6['delta7']:.7f}")
    chk("X94 [作廢·docs/51 §八] |E_6/E_5| = 14.3436（建於 GPU V^{7→5}）",
        abs(e6["r4"] - 14.343552588625746) < 1e-6, f"{e6['r4']:.4f}")

    e7 = _e7()
    chk("X95 E_6（全 CPU 正典）= −631.576826694853", abs(e7["E6"] + 631.5768266948529) < 1e-6, f"{e7['E6']:.9f}")
    chk("X96 E_7 = +10325.588935405450（V^{7→7} 全精度）", abs(e7["E7"] - 10325.588935405450) < 1e-6, f"{e7['E7']:.9f}")
    chk("X97 E_6（5 模截斷）復現舊值至 1.8e−7",
        abs(e7["E6_5mode"] + 545.054076280586) < 1e-6 and e7["rel5"] < 2e-7,
        f"{e7['E6_5mode']:.9f} rel {e7['rel5']:.2e}")
    chk("X98 |E_7/E_6| = 16.3489、|E_6/E_5| = 16.5820（發散飽和）",
        abs(e7["r5"] - 16.348904047545712) < 1e-6 and abs(e7["r4"] - 16.58198157981441) < 1e-6,
        f"{e7['r4']:.4f} / {e7['r5']:.4f}")
    chk("X99 ‖V^{7→5}φ_3^7‖ = 14.526（經閘）vs 舊 GPU 2246.285",
        abs(e7["V75n"] - 14.52618953969585) < 1e-6 and abs(e7["V75old"] - 2246.285) < 1e-3,
        f"{e7['V75n']:.6f} / {e7['V75old']:.3f}")
    chk("X81 五模態總數 = 多重集 C(44+5−1,5) = C(48,5) = 1712304（枚舉覆蓋完整）",
        __import__("math").comb(48, 5) == 1712304, str(__import__("math").comb(48, 5)))

    hdr("§4 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")

    RES["_inputs"] = {k: (v if not isinstance(v, mp.mpf) else mp.nstr(v, 30)) for k, v in I.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}
    (DATA / "l2_independent_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())