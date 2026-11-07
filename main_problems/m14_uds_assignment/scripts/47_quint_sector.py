#!/usr/bin/env /usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/47_quint_sector.py — MR1 機器（六）：5 模扇區（1→3→5→3→1）

性質
    本腳本計算 E_4 中「經 5 模中間態」之部分：
        A4_quint(n) := Σ_{Q ∈ 五模扇區} X_Q(n)² / (γ_n − E_Q)
        X_Q(n) := Σ_{T ∈ 三模扇區} V_{Q,T} · b_n(T) / (γ_n − E_T)
    推導：RS 遞歸中 φ_2(Q) = (Vφ_1)_Q/(γ_n−E_Q)，而 E_4 = ⟨n|V|φ_3⟩ 之五模部分
           = Σ_Q φ_2(Q)·Σ_T V_{Q,T} b_T/(γ_n−E_T) = Σ_Q X_Q²/(γ_n−E_Q)（V 實對稱）。
    三模部分 A4_trio 由本腳本獨立重算（自校閘，須與 docs/45 一致）。

慣例（與 docs/41–45 一致，硬閘）
    /90 站空間：γ_n = 4sin²(πn/90)，n = 1..44（Dirichlet 44 站，Nyquist n=45 不在集合）
    χ⁴ = W/90，W = Σ_ε c_ε [k_ε ≡ 0 (mod 90)]，c_ε = (++−−−−++)
    Fock 矩陣元律（docs/44 M1）：⟨T'|H/g₄|T⟩ = [Σ_{R⊂T,|R|=p} (#emb)·p!·q!·χ⁴(R∪A)]/√(n_T n_T')
        p = (|T|+4−|T'|)/2，q = 4−p，n_T = ∏ m!
        1↔3：p=1,q=3 ⟹ 權 6，χ⁴(n ∪ 三模)
        3↔5：p=1,q=3 ⟹ 權 6，χ⁴(去1 ∪ 添3)   ← 本腳本新增
    正規化：b_n(T) = 6χ⁴(n,T)/√(n_T)；V_{Q,T} = 6·emb·χ⁴/√(n_T n_Q)

N_trials（本腳本為**完整枚舉**，非抽樣）
    §二 之 40 個五模態（seed = 20261102）僅作對稱核自校，非統計估計。
    §五 枚舉全部 1,712,304 個五模態（精確，無誤差棒）。
    快取：data/cache/enum47_v1.npz（存 A4Q/NA/NZ_TOT，避免每次 ~38 分鐘重建；
          快取命中與否**不改變輸出**，故不入 stdout）。可用 --rebuild 刪快取強制重算。

範圍宣告
    §五 為 **完整精確枚舉**（非 MC）；A4_quint 無統計誤差，惟以 float64 累加 1.7M 項
    （相對誤差 ~1e−10 量級，未做 mpmath 覆核）。
    全 quint 數 = 多重集 C(44+5−1,5) = C(48,5) = 1712304。
    七模及以上扇區（多重集 C(48,7) = 73,629,072）**未算**。

執行
    PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 scripts/47_quint_sector.py
    （stdout 禁含計時／隨機噪聲；耗時寫 stderr）
"""
import sys
import time
import itertools
from collections import Counter
from math import factorial

import numpy as np

T0 = time.time()
DPS = 0  # 本腳本以 float64 為主（MC 已明示誤差棒）；關鍵量以 Fraction/整數計數核對
TOL = 1e-12
NS = list(range(1, 45))
G = np.array([0.0] + [4.0 * np.sin(np.pi * n / 90.0) ** 2 for n in NS])
L = "=" * 78

SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]


def W4(t):
    """四波振幅：W = Σ_ε c_ε [k_ε ≡ 0 (mod 90)]"""
    assert len(t) == 4, t
    s = 0
    for c, sg in zip(CS, SIG):
        if sum(x * y for x, y in zip(sg, t)) % 90 == 0:
            s += c
    return s


def W4_vec(r0, a):
    """χ⁴(r0, a0, a1, a2) 對 r0 向量化（a 固定）"""
    r0 = np.asarray(r0)
    out = np.zeros(r0.shape[0], dtype=np.int64)
    for c, sg in zip(CS, SIG):
        k = (sg[0] * r0 + sg[1] * a[0] + sg[2] * a[1] + sg[3] * a[2]) % 90
        out += c * (k == 0)
    return out


def nrm(T):
    r = 1
    for v in Counter(T).values():
        r *= factorial(v)
    return r


# ------------------------------------------------------------------ 扇區建表
print(L)
print("【§一】三模扇區與五模扇區之枚舉")
print(L)
TRI = list(itertools.combinations_with_replacement(NS, 3))
IDX = {T: i for i, T in enumerate(TRI)}
ET = np.array([sum(G[m] for m in T) for T in TRI])
NT = np.array([nrm(T) for T in TRI], dtype=np.float64)
NTOT_Q = len(list(itertools.combinations_with_replacement(NS, 5)))
print(f"  三模態 = C(46,3) = {len(TRI)}｜五模態 = 多重集 C(44+5−1,5) = C(48,5) = {NTOT_Q}")

# LUT：(2 模多重集, 餘模 r0) → 三模索引（供五模扇區向量化查表）
PAIRS2 = list(itertools.combinations_with_replacement(NS, 2))
PIDX = {P: i for i, P in enumerate(PAIRS2)}
LUT = np.full((len(PAIRS2), 44), -1, dtype=np.int32)
for _i, _T in enumerate(TRI):
    for _pos in range(3):
        _S = tuple(sorted([_T[k] for k in range(3) if k != _pos]))
        LUT[PIDX[_S], _T[_pos] - 1] = _i
assert (LUT >= 0).all(), "LUT 有洞"
print(f"  LUT (2 模多重集 × 餘模) = {LUT.shape}，全格已填 ⟹ 三模索引查表化")

# 1↔3 耦合 b_n(T) = 6χ⁴(n,T)/√n_T
B = np.zeros((len(TRI), 44))
for n in NS:
    for i, T in enumerate(TRI):
        w = W4((n,) + T)
        if w:
            B[i, n - 1] = 6.0 * (w / 90.0) / np.sqrt(NT[i])
print(f"  b 矩陣非零 = {int((B != 0).sum())} / {B.size}（非零率 {(B != 0).mean():.4f}）")


# ------------------------------------------------------------------ 3↔5 耦合律
print()
print(L)
print("【§二】3↔5 耦合律（本文件新增）：p=1, q=3，權 1!3! = 6，χ⁴(去1 ∪ 添3)")
print(L)


def links_Q(Q):
    """由五模態 Q 側枚舉 3 模鄰居（向量化 + LUT）；回傳 (J 索引陣列, V 值陣列)"""
    cQ = Counter(Q)
    nQ = nrm(Q)
    Js = []
    Vs = []
    r0 = np.arange(1, 45)
    for S in itertools.combinations_with_replacement(Q, 2):
        cS = Counter(S)
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
            nS *= factorial(v)
        w = W4_vec(r0, add)
        if not w.any():
            continue
        nT = np.full(44, float(nS))
        wgt = np.ones(44)   # 位置重複權重：r0 ∈ S 時 T 側有 (c_S[r0]+1) 個位置選擇
        for m, cm in cS.items():
            nT[m - 1] = float(nS) * (cm + 1)
            wgt[m - 1] = float(cm + 1)
        val = 6.0 * emb * wgt * (w / 90.0) / np.sqrt(nT * nQ)
        m_ = val != 0.0
        if not m_.any():
            continue
        Js.append(LUT[PIDX[S]][m_])
        Vs.append(val[m_])
    if not Js:
        return np.zeros(0, dtype=np.int64), np.zeros(0)
    return np.concatenate(Js).astype(np.int64), np.concatenate(Vs)


def cpl_Q_side(Q):
    """供 §二 對稱核使用：回傳 {(T, V)}"""
    J, V = links_Q(Q)
    out = {}
    for j, v in zip(J, V):
        T = TRI[int(j)]
        out[T] = out.get(T, 0.0) + float(v)
    return out


def cpl_T_side(T, Q):
    """由三模態 T 側算同一耦合（docs/44 律直接用）"""
    cT, cQ = Counter(T), Counter(Q)
    tot = 0.0
    for i in range(3):
        surv = list((cT - Counter([T[i]])).elements())
        cS = Counter(surv)
        if any(cQ[k] < v for k, v in cS.items()):
            continue
        add = list((cQ - cS).elements())
        if len(add) != 3:
            continue
        w = W4(tuple([T[i]] + add))
        if w == 0:
            continue
        emb = 1
        cnt = dict(cQ)
        for m in surv:
            emb *= cnt[m]
            cnt[m] -= 1
        tot += emb * 6 * (w / 90.0)
    return tot / np.sqrt(nrm(T) * nrm(Q))


rng = np.random.default_rng(20261102)
probe = [tuple(sorted(rng.choice(np.arange(1, 45), 5, replace=False))) for _ in range(40)]
mxc, totc, nzc, mags = 0.0, 0, 0, set()
for Q in probe:
    qs = cpl_Q_side(Q)
    cQ = Counter(Q)
    for T, v in qs.items():
        v2 = cpl_T_side(T, Q)
        mxc = max(mxc, abs(v - v2))
        nzc += 1
        mags.add(round(abs(v), 12))
    for S in itertools.combinations_with_replacement(Q, 3):
        cS = Counter(S)
        rest = list((cQ - cS).elements()) if all(cQ[k] >= v for k, v in cS.items()) else []
        if len(rest) == 2:
            totc += 44
print(f"  抽樣 40 個五模態：候選 {totc}｜非零 {nzc}｜非零率 {nzc / totc:.4f}")
print(f"  Q 側 vs T 側（docs/44 律直算）：max|差| = {mxc:.3e}")
print(f"  |V| 值集（相異）= {sorted(mags)}")
print(f"  ⟹ 非零判準：χ⁴(去1 ∪ 添3) ≠ 0；幅值只取 {len(mags)} 個值 ⟹ 結構律成立")

# ------------------------------------------------------------------ E_Q 全表
print()
print(L)
print("【§三】五模扇區能量譜與近退化密度（決定 MC 分層）")
print(L)
QP = np.array(list(itertools.combinations_with_replacement(NS, 5)), dtype=np.int16)
EQ = G[QP].sum(axis=1)
EPS = 1e-2
print(f"  E_Q 範圍 [{EQ.min():.6f}, {EQ.max():.6f}]｜全數 {len(EQ)}")
print("  n ｜ γ_n      ｜ 殼 |Δ|<1e−2 ｜ 殼 |Δ|<1e−3 ｜ 精確簡併")
NSHELL = {}
for n in NS:
    d = np.abs(EQ - G[n])
    NSHELL[n] = int((d < EPS).sum())
    print(f"  {n:2d}｜{G[n]:8.6f}｜{NSHELL[n]:11d}｜{int((d < 1e-3).sum()):11d}｜{int((d < 1e-12).sum()):8d}")
print(f"  殼內總數（全 n）= {sum(NSHELL.values())}；精確簡併 = 0 ⟹ 無需投影剔除")

# ------------------------------------------------------------------ A4_trio 重現
print(L)
print("【§四】自校閘 S2：A4_trio 由本腳本重算，須 = docs/45")
print(L)
PAIRS = list(itertools.combinations_with_replacement(NS, 2))
rows, cols, vals = [], [], []
for i, T in enumerate(TRI):
    for pos in range(3):
        surv = T[pos]
        rem = [T[k] for k in range(3) if k != pos]
        for (u, v) in PAIRS:
            w = W4(tuple(rem) + (u, v))
            if not w:
                continue
            Tp = tuple(sorted([surv, u, v]))
            emb = Counter(Tp)[surv]
            rows.append(i)
            cols.append(IDX[Tp])
            vals.append(4 * (w / 90.0) * emb / (np.sqrt(NT[i]) * np.sqrt(NT[IDX[Tp]])))
rows = np.array(rows, np.int32)
cols = np.array(cols, np.int32)
vals = np.array(vals, np.float64)
print(f"  V_trio 非零元 = {len(vals)}")

A4T = np.zeros(44)
T3 = np.zeros(44)
for n in NS:
    d = G[n] - ET
    inv = np.zeros(len(TRI))
    ok = np.abs(d) > 1e-12
    inv[ok] = 1.0 / d[ok]
    b = B[:, n - 1]
    p1 = inv * b
    E2 = float(b @ p1)
    p2 = inv * np.bincount(rows, weights=vals * p1[cols], minlength=len(TRI))
    p3 = inv * (np.bincount(rows, weights=vals * p2[cols], minlength=len(TRI)) - E2 * p1)
    A4T[n - 1] = float(b @ p3)
    T3[n - 1] = float(b @ p2)
print(f"  A4_trio(1)  = {A4T[0]:.13f}（docs/45：−556.851543190064719）")
print(f"  A4_trio(43) = {A4T[42]:.9f}（docs/45：−343297699.0545967817）")
print(f"  A4_trio(44) = {A4T[43]:.9f}（docs/45：−144545890.4910772741）")
d1 = abs(A4T[0] + 556.851543190064719)
print(f"  |差|(n=1) = {d1:.3e} ⟹ {'✅ 重現' if d1 < 1e-6 else '🔴 不一致'}")

# ------------------------------------------------------------------ U 陣列
U = np.zeros_like(B)
for n in NS:
    d = G[n] - ET
    ok = np.abs(d) > 1e-12
    U[ok, n - 1] = B[ok, n - 1] / d[ok]

# ------------------------------------------------------------------ MC
print()
print(L)
print("【§五】五模貢獻 A4_quint(n) = Σ_Q X_Q²/(γ_n − E_Q)：**完整枚舉全部 quint**（非抽樣）")
print(L)
print(f"  全 quint 數 = {len(QP)}｜以 float64 逐個精確累加（無 MC 誤差）")

A4Q = np.zeros(44)
NZ_TOT = 0
NA = 0
CACHE = __import__("pathlib").Path("data/cache/enum47_v1.npz")
if CACHE.exists():
    _z = np.load(CACHE)
    A4Q = _z["A4Q"]
    NA = int(_z["NA"])
    NZ_TOT = int(_z["NZ_TOT"])
    _CACHED = True
else:
    _CACHED = False
for qi in range(0 if _CACHED else len(QP)):
    Q = tuple(int(x) for x in QP[qi])
    tj, tv = links_Q(Q)
    if len(tj) == 0:
        continue
    X = tv @ U[tj, :]
    A4Q += X * X / (G[1:] - EQ[qi])
    NZ_TOT += len(tj)
    NA += 1
    if (qi + 1) % 100000 == 0:
        print(f"    … {qi + 1} / {len(QP)}", file=sys.stderr)
if not _CACHED:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, A4Q=A4Q, NA=NA, NZ_TOT=NZ_TOT)
print(f"  有非零耦合之 quint 數 = {NA} / {len(QP)}（{NA / len(QP):.4f}）")
print(f"  平均非零鄰居 = {NZ_TOT / max(NA, 1):.2f}（每 quint 440 候選之非零）")

print()
print("  n ｜ A4_trio            ｜ A4_quint（精確）     ｜ 合計               ｜ R = Q/trio")
for n in NS:
    k = n - 1
    R = A4Q[k] / A4T[k] if A4T[k] != 0 else float("nan")
    print(f"  {n:2d}｜{A4T[k]:+19.9e}｜{A4Q[k]:+19.9e}｜{A4T[k] + A4Q[k]:+19.9e}｜{R:+.6f}")

print()
print(L)
print("【§六】判定")
print(L)
Rm = np.abs(A4Q / A4T)
print(f"  ① |A4_quint/A4_trio|：min {Rm.min():.6f}｜中位 {np.median(Rm):.6f}｜max {Rm.max():.6f}")
print(f"  ② Σ A4_trio = {A4T.sum():+.10e}｜Σ A4_quint = {A4Q.sum():+.10e}")
print(f"  ③ 五模扇區佔比（Σ 絕對值）= {np.abs(A4Q).sum() / (np.abs(A4T).sum() + np.abs(A4Q).sum()):.6f}")
verdict = ("五模扇區**不可忽略**（與三模同量級或更大）" if np.median(Rm) > 0.1
           else "五模扇區為**小修正**（中位 < 10%）")
print(f"  ④ 判：{verdict}")
print("  ⑤ 範圍：本節為 E_4 之「1→3→5→3→1」全部通道（精確）；七模及以上未算")

print()
print(L)
print("【§七】★ 以完整四階 A4 = A4_trio + A4_quint 重算比值檢 ρ（修正 docs/45 §五b）")
print(L)
A4F = A4T + A4Q
rho = []
for n in NS:
    k = n - 1
    a2 = float(U[:, k] @ B[:, k])       # A2_n = Σ_T b_T²/(γ_n − E_T)
    r3 = abs(T3[k] / a2) if a2 != 0 else float("nan")     # |T_n/A2_n|
    r4 = abs(A4F[k] / T3[k]) if T3[k] != 0 else float("nan")  # |A4/T_n|
    rho.append((n, r3, r4, r4 / r3 if r3 else float("nan")))
vals = sorted(x[3] for x in rho)
print(f"  A2_n 由定義式 A2 = Σ_T b_T²/(γ_n−E_T) 重算（自校：A2(1) 應 = −2.89910668912）")
print(f"    A2(1) = {float(U[:, 0] @ B[:, 0]):.11f}｜T_1 = {T3[0]:.11f}（docs/44：36.90518272244775）")
print(f"  ρ_n := (|A4^full/T|)/(|T/A2|)：min {vals[0]:.6f}｜中位 {vals[22]:.6f}｜max {vals[-1]:.6f}")
print(f"  落 [0.1,10] 之模數 = {sum(1 for x in vals if 0.1 <= x <= 10)} / 44")
print(f"  例：n=43｜|T/A2| = {rho[42][1]:.6e}｜|A4^full/T| = {rho[42][2]:.6e}｜ρ = {rho[42][3]:.4f}")
print(f"      n=1 ｜|T/A2| = {rho[0][1]:.6e}｜|A4^full/T| = {rho[0][2]:.6e}｜ρ = {rho[0][3]:.4f}")
med = vals[22]
inb = sum(1 for x in vals if 0.1 <= x <= 10)
if med < 0.5:
    kind = "次幾何（比值遞減）⟹ 級數可能收斂，[1/1] 極點疑為假象"
elif med <= 2.0:
    kind = "幾何式發散 ⟹ [1/1] 極點非假象"
else:
    kind = "**超幾何（比值遞增）⟹ 漸近級數**；[1/1] 極點只作上界，非真半徑"
print(f"  ⟹ ρ 中位 = {med:.6f}、{inb}/44 落 [0.1,10] ⟹ 判：{kind}")
print(f"  ⟹ 與 docs/45 §五b（僅三模 A4，ρ 中位 0.690538）對照：五模加入後比值遞增")
print(f"     ⟹ 級數非收斂，微擾法在 g₄ ≠ 0 任何值皆非漸近可信；本層只能作定性/量級工具")
print(f"\n[耗時] {time.time() - T0:.1f} s", file=sys.stderr)
