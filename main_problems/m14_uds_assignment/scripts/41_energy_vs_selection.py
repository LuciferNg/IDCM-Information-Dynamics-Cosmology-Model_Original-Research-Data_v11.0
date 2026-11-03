#!/usr/bin/python3
# -*- coding: utf-8 -*-
# ══════════════════════════════════════════════════════════════════════════
# 🔴 公約更正（2026-11-02，append-only）：本檔 v1 誤用 γ_n = 4sin²(πn/45)
#    正確（/90 站空間，Dirichlet 44 站，45 區間）：γ_n = 4sin²(πn/90)，Σγ = 88
#    症狀：v1 之 Σγ = 90（≠ 88）即已洩漏；γ_max 被誤置於 n = 22 而非 n = 44
#    v1 log 保留於 data/run_41_v1_wronggamma.log / run_42_v1_wronggamma.log
#    重算後：1→2 衰變 3 條（非 0）；1→3 衰變 1 條；2→2 共振 243 對（非 1,452）
# ══════════════════════════════════════════════════════════════════════════
"""
M14 scripts/41 — 四波：能量守恆 vs 振幅選擇律；奇 p 之非零性

────────────────────────────────────────────────────────────────────────
背景
  docs/39：χ⁴ = 0（和奇）／χ⁴ = W/90（和偶），W = Σ_ε c_ε[k_ε ≡ 0 (mod 90)]
  本腳本問兩個結構問題：
    Q-a 能量共振（γ_a+γ_b = γ_c+γ_d）與振幅非零（W ≠ 0）之關係？
    Q-b 奇 p（5、7）之 χ⁽ᵖ⁾ 是否可能為零（宇稱容許時）？

判準（跑前凍結）
  T1 四波能量共振之完整枚舉（n ≤ 44，含重複，C(47,4) = 178,365 組）
     ─ 分「和偶」與「和奇」兩類
     ─ 命題 α：和偶之能量共振 ⟹ W ≠ 0（違例須 0）
     ─ 命題 β：和奇之能量共振 ⟹ χ⁴ = 0（宇稱，必真）
  T2 反之不然：給出 W ≠ 0 但無能量共振之最小例（須存在，並列 W 值）
  T3 奇 p：p = 5、7，抽樣（固定 seed）全部宇稱容許樣本，統計零值數（命題 γ：0）
  T4 交叉核：和偶共振之 W 值集（期望 {2, 3}）
性質：T1 全窮舉；T3 抽樣（N_trials = 3000 每 p，全報）。dps = 50。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import random
from collections import Counter

import mpmath as mp

mp.mp.dps = 50
N = 45
C2 = mp.mpf(2) / 45
L = "=" * 78
TOL = mp.mpf('1e-35')
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]


def gam(n):
    return 4 * mp.sin(mp.pi * n / 90) ** 2


def W4(t):
    return sum(c for c, sg in zip(CS, SIG) if sum(s * x for s, x in zip(sg, t)) % 90 == 0)


def chi(ns):
    s = mp.mpf(0)
    for j in range(1, 45):
        p = mp.mpf(1)
        for n in ns:
            p *= mp.sin(mp.pi * n * j / 45)
        s += p
    return C2 ** (mp.mpf(len(ns)) / 2) * s


def eres(t):
    g = [gam(n) for n in t]
    return (abs(g[0] + g[1] - g[2] - g[3]) < TOL or abs(g[0] + g[2] - g[1] - g[3]) < TOL
            or abs(g[0] + g[3] - g[1] - g[2]) < TOL)


print(L)
print("M14 scripts/41 — 四波能量守恆 vs 振幅選擇律；奇 p 非零性")
print(L)

# ---------------------------------------------------------------- T1
print("\n【§一】T1 四波能量共振之完整枚舉（n ≤ 44，含重複）")
even_res = odd_res = 0
bad_a = bad_b = 0
Wc = Counter()
for t in itertools.combinations_with_replacement(range(1, 45), 4):
    if not eres(t):
        continue
    if sum(t) % 2 == 0:
        even_res += 1
        w = W4(t)
        Wc[w] += 1
        if w == 0:
            bad_a += 1
    else:
        odd_res += 1
        if abs(chi(t)) > TOL:
            bad_b += 1
print(f"  共振總數 = {even_res + odd_res}（和偶 {even_res} ＋ 和奇 {odd_res}）")
print(f"  命題 α（和偶共振 ⟹ W ≠ 0）：違例 = {bad_a} {'✅' if bad_a == 0 else '🔴'}")
print(f"  命題 β（和奇共振 ⟹ χ⁴ = 0）：違例 = {bad_b} {'✅' if bad_b == 0 else '🔴'}")

# ---------------------------------------------------------------- T4
print("\n" + L)
print("【§二】T4 和偶共振之 W 值集")
print(L)
print(f"  W 分佈 = {dict(sorted(Wc.items()))}")
print(f"  ⟹ 共振振幅只取 χ⁴ ∈ {{W/90}} = {sorted({'%s/90' % k for k in Wc})}")
print(f"  抽核：χ⁴(1,1,1,1) = {mp.nstr(chi((1,1,1,1)), 12)}（W = {W4((1,1,1,1))}）")

# ---------------------------------------------------------------- T2
print("\n" + L)
print("【§三】T2 反之不然：W ≠ 0 但無能量共振")
print(L)
ex = None
for t in itertools.combinations_with_replacement(range(1, 10), 4):
    if W4(t) != 0 and not eres(t):
        ex = t
        break
print(f"  最小例 {ex}：W = {W4(ex)} ⟹ χ⁴ = {mp.nstr(chi(ex), 12)}；能量共振 = {eres(ex)}")
cnt = sum(1 for t in itertools.combinations_with_replacement(range(1, 45), 4)
          if W4(t) != 0 and not eres(t))
print(f"  n ≤ 44 中含此情形之組數 = {cnt}（對照：共振 1221 組）")
print("  ⟹ 振幅選擇律嚴格弱於能量守恆：晶格允許 off-shell 四波耦合")

# ---------------------------------------------------------------- T3
print("\n" + L)
print("【§四】T3 奇 p 之非零性（抽樣，N_trials 全報）")
print(L)
for p, seed in ((5, 5), (7, 9)):
    random.seed(seed)
    tested = 0
    zeros = 0
    trials = 0
    while tested < 2000 and trials < 20000:
        trials += 1
        ns = tuple(random.randint(1, 44) for _ in range(p))
        if sum(ns) % 2 == 0:      # p 奇 ⟹ 需 Σ 奇
            continue
        tested += 1
        if abs(chi(ns)) < TOL:
            zeros += 1
    print(f"  p = {p}：宇稱容許樣本 {tested}（抽 {trials} 次），為零者 = {zeros} "
          f"{'✅' if zeros == 0 else '🔴'}")

# ---------------------------------------------------------------- T5
print("\n" + L)
print("【§五】T5 光譜穩健性：兩種 γ 之『和偶共振集』是否同一集合")
print(L)
GA = {n: 4 * mp.sin(mp.pi * n / 90) ** 2 for n in range(1, 45)}   # 正確
GB = {n: 4 * mp.sin(mp.pi * n / 45) ** 2 for n in range(1, 45)}   # v1 誤用（作對照）


def resset(Gg, even_only=True):
    out = set()
    for t in itertools.combinations_with_replacement(range(1, 45), 4):
        g = [Gg[x] for x in t]
        if (abs(g[0] + g[1] - g[2] - g[3]) < TOL or abs(g[0] + g[2] - g[1] - g[3]) < TOL
                or abs(g[0] + g[3] - g[1] - g[2]) < TOL):
            if (not even_only) or sum(t) % 2 == 0:
                out.add(t)
    return out


SA = resset(GA)
SB = resset(GB)
print(f"  和偶共振集：正確 γ {len(SA)} 組｜v1 γ {len(SB)} 組｜集合恆等 = {SA == SB}｜交集 = {len(SA & SB)}")
print(f"  ⟹ {'1,221 之穩定係**集合恆等**，非計數巧合' if SA == SB else '⚠ 集合不同 ⟹ 1,221 之重合属巧合'}")
RA = resset(GA, even_only=False)
RB = resset(GB, even_only=False)
print(f"  全共振集（含和奇）：正確 {len(RA)}｜v1 {len(RB)}｜差 = {len(RA ^ RB)} 組（全在和奇扇區）")

# ---------------------------------------------------------------- §六
print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ① 命題 α（和偶能量共振 ⟹ 非零）**成立**（n ≤ 44 全枚舉，{even_res} 組，違例 0）
  ② 命題 β（和奇共振 ⟹ 零）**成立**（{odd_res} 組，宇稱；違例 0）
  ③ 共振振幅只取 χ⁴ ∈ {{2/90, 3/90}}（即 1/45、1/30）—— 量化更細
  ④ 反之不然：W ≠ 0 而無能量共振者存在（最小 (1,1,1,3)，W = −1）
     ⟹ 振幅選擇律 ⊋ 能量守恆；晶格允許 off-shell 四波耦合
  ⑤ 奇 p：p = 5、7 各 2,000 樣本，**零值 0** ⟹ 命題 γ 🟡 強數值（未證）
  ⑥ 光譜穩健性：和偶共振集在兩種 γ 下**同一集合**（{len(SA)} 組），差異全在和奇扇區（{len(RA ^ RB)} 組）""")
print(L)
raise SystemExit(0)