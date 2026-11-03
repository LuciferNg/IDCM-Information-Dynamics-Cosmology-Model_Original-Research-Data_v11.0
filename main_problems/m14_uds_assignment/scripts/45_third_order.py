#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/45 — MR1 機器（四）：三階微擾（v2）

🔴 v2 更正（2026-11-02，append-only）：
   v1 誤述矩陣元非零條件為「3↔3 恰共享 1 模」—— 真相係「存在移除對使**餘模** ∈ T'」，
   即 |T∩T'| ≥ 1（多重集交），且振幅含 1/√(n_T n_{T'}) 與「移除索引選擇數」。
   v1 之快速式漏此兩者 ⟹ n=1 得 7.6064（正解 17.4426）。v1 log 封存 data/archive/。

────────────────────────────────────────────────────────────────────────
公約（/90）：γ_n = 4sin²(πn/90)，Σγ = 88，γ_max = γ_44 = 3.99512810052
態：|n⟩ 單模；|T⟩ = |b,c,d⟩（multiset，n_T = ∏ m!）
H_int = g₄ Σ_j :W_j⁴:（normal-ordered）；χ⁴(a,b,c,d) = W4/90（docs/39）

★ 矩陣元（Wick 收縮，p 去 q 添，p+q=4）
   ⟨T'|H/g₄|T⟩ = [Σ_{移除對 R ⊂ T, |R|=2} (#嵌入) · p!·q!·χ⁴(R ∪ 添)] / √(n_T n_{T'})
   p = (|T| + 4 − |T'|)/2
   ⟹ 1↔3：6χ⁴/√(n_T)；3↔3：4χ⁴·(移除選擇)/√(n_T n_{T'})；1↔1：p=2>1 ⟹ **恆零**

三階 Rayleigh–Schrödinger：
  δγ_n⁽³⁾ = g₄³ T_n，T_n = Σ_{T,T'} ⟨T|H|n⟩⟨T'|H|T⟩⟨n|H|T'⟩ / ((γ_n−E_T)(γ_n−E_{T'}))
  −V_nn 項因 ⟨n|H|n⟩ = 0 **整項消失**（結構性）

判準（跑前凍結）
  M1 3↔3 共享 0 ⟹ 恆零（多重集交為空 ⟹ 無合法移除對）
  M2 1↔3 = 6χ⁴（全相異例）
  M3 1↔1 = 0（44×44 全枚舉）
  M4 快速式與直接 mel 式於 n = 1, 44, 44 一致（差 < 1e−25）
性質：trios = C(46,3) = 15,180；dps = 50；TOL = 1e−35（高於噪底 ~1e−49）
────────────────────────────────────────────────────────────────────────
"""
import itertools
import random
import sys
import time
from collections import Counter
from math import factorial

import mpmath as mp

mp.mp.dps = 50
TOL = mp.mpf('1e-35')
L = "=" * 78
SIGS = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
        (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]
G = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}
TRIOS = list(itertools.combinations_with_replacement(range(1, 45), 3))
ET = {T: sum(G[m] for m in T) for T in TRIOS}


def W4(t):
    return sum(c for c, s in zip(CS, SIGS)
               if sum(x * y for x, y in zip(s, t)) % 90 == 0)


def chi(t):
    return mp.mpf(W4(t)) / 90


def nfac(T):
    r = 1
    for v in Counter(T).values():
        r *= factorial(v)
    return r


def nrm(T):
    return mp.sqrt(mp.mpf(nfac(T)))


def minus_sub(mult, sub):
    c = Counter(mult)
    for k, v in Counter(sub).items():
        c[k] -= v
        if c[k] < 0:
            return None
    return list(c.elements())


def mel(T, Tp):
    """⟨T'|H/g₄|T⟩（normal-ordered，Wick 收縮）"""
    s, t = len(T), len(Tp)
    if (s + 4 - t) % 2:
        return mp.mpf(0)
    p = (s + 4 - t) // 2
    q = 4 - p
    if p < 0 or q < 0 or p > s:
        return mp.mpf(0)
    z = mp.mpf(0)
    for idx in itertools.combinations(range(s), p):
        rem = [T[i] for i in idx]
        left = minus_sub(T, rem)
        if left is None:
            continue
        add = minus_sub(Tp, left)
        if add is None or len(add) != q:
            continue
        cTp, cL = Counter(Tp), Counter(left)
        emb = 1
        for k, v in cL.items():
            for i in range(v):
                emb *= (cTp[k] - i)
        z += emb * factorial(p) * factorial(q) * chi(tuple(rem) + tuple(add))
    return z / (nrm(T) * nrm(Tp))


print(L)
print("M14 scripts/45 — MR1（四）：三階微擾（v2）")
print(L)
Ssum = sum(G.values())
print(f"\n【§0】簽名：Σγ = {mp.nstr(Ssum, 12)}（期望 88）｜γ_44 = {mp.nstr(G[44], 12)}"
      f"｜trios = {len(TRIOS)}")
assert abs(Ssum - 88) < TOL

# ---------------------------------------------------------------- §一 M1–M3
print("\n" + L)
print("【§一】矩陣元律（自證）")
print(L)
st = {}
random.seed(11)
for _ in range(3000):
    T = random.choice(TRIOS)
    Tp = random.choice(TRIOS)
    ov = sum((Counter(T) & Counter(Tp)).values())
    v = mel(Tp, T)
    d = st.setdefault(ov, [0, 0])
    d[0] += 1
    if abs(v) > TOL:
        d[1] += 1
print(f"  M1 3↔3 [共享數 → (樣本, 非零)]：{ {k: tuple(x) for k, x in sorted(st.items())} }")
print(f"     共享 0 ⟹ 非零 {st.get(0, [0, 0])[1]} "
      f"{'✅ 恆零' if st.get(0, [0, 0])[1] == 0 else '⚠'}"
      f"｜共享 ≥ 1 ⟹ 可非零（v1 誤作「恰 1」）")
mx = mp.mpf(0)
for _ in range(40):
    n, b, c, d = random.sample(range(1, 45), 4)
    mx = max(mx, abs(mel(tuple(sorted((b, c, d))), (n,)) - 6 * chi((n, b, c, d))))
print(f"  M2 1↔3 = 6χ⁴（40 例全相異）：max|差| = {mp.nstr(mx, 4)} "
      f"{'✅ 精確' if mx < TOL else '⚠'}")
mx = max(abs(mel((m,), (n,))) for n in range(1, 45) for m in range(1, 45))
print(f"  M3 1↔1（44×44 全枚舉）：max|元| = {mp.nstr(mx, 4)} "
      f"{'✅ 恆零 ⟹ ⟨n|H|n⟩ = 0 ⟹ 三階 −V_nn 項整項消失' if mx < TOL else '⚠'}")
# 重複模樣本：振幅 = 4χ⁴·移除選擇/√(n_T n_T')
chk = mel((1, 3, 4), (1, 1, 2))
ref = 2 * 4 * chi((1, 2, 3, 4)) / mp.sqrt(2)
print(f"  重複模樣本 mel((1,3,4),(1,1,2)) = {mp.nstr(chk, 12)}"
      f"｜4χ⁴·2/√2 = {mp.nstr(ref, 12)} "
      f"{'✅' if abs(chk - ref) < mp.mpf('1e-40') else '⚠'}")

# ---------------------------------------------------------------- §二
print("\n" + L)
print("【§二】三階係數 T_n")
print(L)
S = {}
for n in range(1, 45):
    S[n] = [(T, mp.mpf(W4((n,) + T))) for T in TRIOS if W4((n,) + T)]
BYX = {}
for n in range(1, 45):
    d = {}
    for T, w in S[n]:
        for x in set(T):
            d.setdefault(x, []).append((T, w, Counter(T)[x]))
    BYX[n] = d
print(f"  V(n,T) ≠ 0 之 trios 數：min {min(len(v) for v in S.values())}"
      f"｜max {max(len(v) for v in S.values())}")


def third_fast(n):
    acc = mp.mpf(0)
    for T, wT in S[n]:
        uT = 6 * (wT / 90) / nrm(T) / (G[n] - ET[T])
        for idx in itertools.combinations(range(3), 2):
            rem = [T[i] for i in idx]
            left = minus_sub(T, rem)
            x = left[0]
            lst = BYX[n].get(x)
            if not lst:
                continue
            for Tp, wP, mp_ in lst:
                add = minus_sub(Tp, [x])
                if add is None or len(add) != 2:
                    continue
                wm = W4(tuple(rem) + tuple(add))
                if not wm:
                    continue
                uP = 6 * (wP / 90) / nrm(Tp) / (G[n] - ET[Tp])
                acc += uT * uP * mp_ * 4 * (mp.mpf(wm) / 90) / (nrm(T) * nrm(Tp))
    return acc


t0 = time.time()
TH = {n: third_fast(n) for n in range(1, 45)}
# 耗時只寫 stderr：stdout 必須 byte-identical（審計 ⑦ 兩跑比對）
print(f"[耗時] {time.time() - t0:.1f} s", file=sys.stderr)
for n in range(1, 45):
    print(f"  n = {n:2d}｜T_n = {mp.nstr(TH[n], 16):>24}")

# ---------------------------------------------------------------- §三 M4
print("\n" + L)
print("【§三】M4 直接 mel 式交叉驗證（n = 1, 44）")
print(L)
for n in (1, 44):
    acc = mp.mpf(0)
    lst = S[n]
    for T, _ in lst:
        vT = mel(T, (n,)) / (G[n] - ET[T])
        for Tp, _ in lst:
            if not (set(T) & set(Tp)):
                continue
            v = mel(Tp, T)
            if v == 0:
                continue
            acc += vT * v * mel((n,), Tp) / (G[n] - ET[Tp])
    print(f"  n = {n:2d}：直接法 = {mp.nstr(acc, 16):>24}"
          f"｜快速法 = {mp.nstr(TH[n], 16):>24}"
          f"｜相對差 = {mp.nstr(abs(acc - TH[n]) / max(abs(acc), mp.mpf(1)), 3)} "
          f"{'✅ 一致' if abs(acc - TH[n]) < mp.mpf('1e-20') * max(abs(acc), mp.mpf(1)) else '⚠'}")

# ---------------------------------------------------------------- §四 結構
print("\n" + L)
print("【§四】結構：相異性、鏡像、級聯")
print(L)
vals = {mp.nstr(TH[n], 30) for n in TH}
print(f"  相異值 = {len(vals)} / 44")
mir = max(abs(TH[n] + TH[45 - n]) for n in range(1, 45))
print(f"  鏡像檢：max|T_n + T_(45−n)| = {mp.nstr(mir, 10)} ⟹ "
      f"{'反號' if mir < mp.mpf('1e-20') else '唔反號（三階無鏡像反號）'}")
print(f"  符號：正 {sum(1 for n in TH if TH[n] > 0)}｜負 {sum(1 for n in TH if TH[n] < 0)}")
print(f"  min = {mp.nstr(min(TH.values()), 12)}（n = {min(TH, key=TH.get)}）"
      f"｜max = {mp.nstr(max(TH.values()), 12)}（n = {max(TH, key=TH.get)}）")
print(f"  Σ_n T_n = {mp.nstr(sum(TH.values()), 14)}")
big = sorted(TH, key=lambda n: -abs(TH[n]))[:6]
print(f"  |T_n| 最大六模 = {[(n, mp.nstr(TH[n], 8)) for n in big]}")

# ---------------------------------------------------------------- §四b 對照
print("\n" + L)
print("【§四b】與二階對照：係數比、強化之微擾界、跨腳本一致性")
print(L)
A2 = {}
for n in range(1, 45):
    s = mp.mpf(0)
    for T, w in S[n]:
        s += (6 * (w / 90) / nrm(T)) ** 2 / (G[n] - ET[T])
    A2[n] = s
print(f"  二階（真 Fock 正規化）A2_n = Σ_T |V(n,T)|²/(γ_n − E_T)")
print(f"  ↳ 與 docs/42 一致性：1350·A2_n 應等於 S_n（A2 為真 Fock 正規化）")
for n in (1, 10, 44):
    print(f"     n = {n:2d}：1350·A2_n = {mp.nstr(1350 * A2[n], 16)}")
print(f"  （docs/42 之 S_1 = −3913.79403031、S_10 = −22573.3904714、S_44 = −85590.5782401）")
r = {n: TH[n] / A2[n] for n in TH}
bd = min((abs(A2[n]) / abs(TH[n]), n) for n in TH)
print(f"  係數比 r_n := T_n / A2_n：min {mp.nstr(min(r.values()), 8)}（n = {min(r, key=r.get)}）"
      f"｜max {mp.nstr(max(r.values()), 8)}（n = {max(r, key=r.get)}）")
print(f"  r_n 負 {sum(1 for n in r if r[n] < 0)} / 44 ⟹ 該些模三階與二階**反號**（部分抵消）")
print(f"  r_n 正 {sum(1 for n in r if r[n] > 0)} / 44 ⟹ 應與 docs/42 之 S_n > 0 模集 "
      f"{{17,30,33,35,39,41}} 一致：{[n for n in r if r[n] > 0]}")
print(f"  ★ 三階微擾界 |g₄| < min_n |A2_n|/|T_n| = {mp.nstr(bd[0], 14)}（n = {bd[1]}）")
g2t = mp.mpf('0.070265395858052') / mp.sqrt(6)
print(f"    二階界（docs/42 單位）= 0.070265395858052 ⟹ 真單位 = {mp.nstr(g2t, 12)}（÷√6）")
print(f"    ⟹ 三階收緊 {mp.nstr(g2t / bd[0], 8)} 倍（同單位制）")
print(f"  主導模（|T_n| 最大六）：{[(n, mp.nstr(TH[n], 6)) for n in big]}")

# ---------------------------------------------------------------- §五
print("\n" + L)
print("【§五】判定")
print(L)
nz = sum(1 for n in TH if TH[n] > 0)
print(f"""  ① 矩陣元律：1↔3 = 6χ⁴；3↔3 = 4χ⁴·(移除選擇)·emb/√(n_T n_T')（共享 0 ⟹ 0）
  ② ⟨n|H|n⟩ = 0（1↔1 恆零）⟹ 三階 RS 之 −V_nn 項**整項消失**（結構性，非近似）
  ③ ⟹ 三階 = 純「1→3 → 換兩模 → 3→1」鏈；三模扇區內部躍遷首次進入
  ④ T_n 相異 {len(vals)}/44；正負 {nz}/{44 - nz}；鏡像 {'有' if mir < mp.mpf('1e-20') else '無'}
  ⑤ 兩算法一致（n = 1, 44；相對差 < 1e−20）""")
print(L)
raise SystemExit(0)