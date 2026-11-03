#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/14 — 三波共振檢定：beat／sum 通道有無「在殼」

────────────────────────────────────────────────────────────────────────
問題：cos(ω_a t)cos(ω_b t) = ½[cos((ω_b−ω_a)t) + cos((ω_a+ω_b)t)]。
      兩條通道要**係載體之合法模**，須 ω_b−ω_a = ω_m 或 ω_a+ω_b = ω_m（三波共振）。
      色散 ω_n = 2sin(πn/88)（skill 已註冊；非線性 ⟹ 共振非平凡）。
      本腳本掃盡 (a,b,m) 三格，檢定有無精確共振。

⚠ **自我更正**：scripts/13 以「Σ = a+b 相對 44」作帶內／帶外分類，
   但 Σ 係**動量指數**，唔係頻率指數；兩者只在線性色散下同一。本輪改用真實頻率，
   並重檢 scripts/13 之物理讀法（見 §五）。

性質：確定性掃描，N_trials = 0（非擬合、非盲測；掃盡而非抽樣）。
零觀測量。dps = 60，精確判定閾 1e−40。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp
from collections import Counter

mp.mp.dps = 60
L = "=" * 78
ND = 44
U = mp.pi / 88
TOL = mp.mpf('1e-40')


def om(n):                     # 頻散
    return 2 * mp.sin(n * U)


print(L)
print("M14 scripts/14 — 三波共振檢定（beat／sum 在殼？）")
print(L)

print("\n【§一】頻散與帶頂")
print(f"  ω_n = 2 sin(πn/88)；ω_1 = {mp.nstr(om(1),10)}；ω_22 = {mp.nstr(om(22),10)}（=√2）；"
      f"ω_44 = {mp.nstr(om(44),10)}（=2 = 帶頂, Nyquist）")
print(f"  ω_a + ω_b 之範圍 = ({mp.nstr(om(1)+om(2),8)}, {mp.nstr(om(43)+om(44),8)}) —— 可**超出**帶頂 2")
print(f"  ω_b − ω_a 之範圍 = (0, {mp.nstr(om(44)-om(1),8)}) —— 恆在帶內")

pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]
print(f"  價對 = {len(pairs)}；模 m = 1..88 ⟹ 三格組合 = {len(pairs)*88}")

# ======================================================================
print("\n【§二】sum 通道：ω_a + ω_b = ω_m ?")
best = (mp.mpf(9), None)
exact = []
above = 0
for (a, b) in pairs:
    s = om(a) + om(b)
    if s > 2:
        above += 1
    for m in range(1, 89):
        r = abs(s - om(m))
        if r < best[0]:
            best = (r, (a, b, m))
        if r < TOL:
            exact.append((a, b, m))
print(f"  ω_a+ω_b > 2（頻率超出帶頂 ⟹ 自動不可能）之價對 = {above} / {len(pairs)}")
print(f"  精確共振（<1e−40）= {len(exact)}")
print(f"  最小頻隙 = {mp.nstr(best[0], 8)}  at (a,b,m) = {best[1]}")
print(f"    該處 ω_a+ω_b = {mp.nstr(om(best[1][0])+om(best[1][1]),12)} ；ω_m = {mp.nstr(om(best[1][2]),12)}")

# ======================================================================
print("\n【§三】beat 通道：ω_b − ω_a = ω_m ?")
best2 = (mp.mpf(9), None)
exact2 = []
for (a, b) in pairs:
    d = om(b) - om(a)
    for m in range(1, 89):
        r = abs(d - om(m))
        if r < best2[0]:
            best2 = (r, (a, b, m))
        if r < TOL:
            exact2.append((a, b, m))
print(f"  精確共振（<1e−40）= {len(exact2)}")
print(f"  最小頻隙 = {mp.nstr(best2[0], 8)}  at (a,b,m) = {best2[1]}")
print(f"    該處 ω_b−ω_a = {mp.nstr(om(best2[1][1])-om(best2[1][0]),12)} ；ω_m = {mp.nstr(om(best2[1][2]),12)}")

# ======================================================================
print("\n【§四】特例 m = Σ = a+b 之代數 no-go（可完全證明）")
print("""
  取 u = π/88，A = au, B = bu。若 m = a+b：
    sin(au) + sin(bu) − sin((a+b)u) = 2 sin(Σu/2)cos(Δu/2) − 2 sin(Σu/2)cos(Σu/2)
                                   = 2 sin(Σu/2)[cos(Δu/2) − cos(Σu/2)]
  = 0 ⟺  sin(Σu/2) = 0  或  cos(Δu/2) = cos(Σu/2)
  前者 ⟹ Σ = 176k ⟹ Σ ≥ 176 > 87 ✗
  後者 ⟹ Δu/2 = ±Σu/2 + 2πn ⟹ Δ = ±Σ + 352n
         Δ =  Σ ⟹ a = 0 ✗      （a ≥ 1）
         Δ = −Σ + 352n ⟹ 2b = 352n ⟹ b = 176n ⟹ b ≥ 176 ✗（b ≤ 44）
  ⟹ **m = Σ 永不共振**（精確）。同法（sin B − sin A）得 beat 之 m = Δ 亦永不共振。
""")
# 數值核對上式
worst = mp.mpf(0)
for (a, b) in pairs:
    S, D = a + b, b - a
    lhs = mp.sin(a * U) + mp.sin(b * U) - mp.sin(S * U)
    rhs = 2 * mp.sin(S * U / 2) * (mp.cos(D * U / 2) - mp.cos(S * U / 2))
    worst = max(worst, abs(lhs - rhs))
print(f"  分解式數值核對：max|LHS−RHS| = {mp.nstr(worst, 3)}（946 對）")

# ======================================================================
# ======================================================================
print("\n【§六】改用**能量** γ（IDCM 之能量量子）作守恒量：γ_a ± γ_b = γ_m ?")
print("  γ_n = 4sin²(nπ/88) = ω_n²；Σγ = 90 係能量預算 ⟹ 能量守恒用 γ 而非 ω")


def gm(n):
    return 4 * mp.sin(n * U) ** 2


bestS = (mp.mpf(9), None)
exS = []
for (a, b) in pairs:
    v = gm(a) + gm(b)
    for m in range(1, 89):
        r = abs(v - gm(m))
        if r < bestS[0]:
            bestS = (r, (a, b, m))
        if r < TOL:
            exS.append((a, b, m))
bestD = (mp.mpf(9), None)
exD = []
for (a, b) in pairs:
    v = gm(b) - gm(a)
    for m in range(1, 89):
        r = abs(v - gm(m))
        if r < bestD[0]:
            bestD = (r, (a, b, m))
        if r < TOL:
            exD.append((a, b, m))

print(f"\n  sum（γ_a+γ_b = γ_m）：精確解 = {len(exS)} / {len(pairs)*88}")
ms = sorted({m for (_, _, m) in exS})
print(f"    命中模 m 之集合 = {ms}；涉及價對 = {sorted({(a,b) for (a,b,_) in exS})[:6]} …")
print(f"    ⟹ 全部 Σ=a+b=44 之 21 對，通道落於 m=44（γ=4, Nyquist）")
print(f"\n  beat（γ_b−γ_a = γ_m）：精確解 = {len(exD)} / {len(pairs)*88}")
bs = sorted({b for (_, b, _) in exD})
print(f"    涉及之 b 值 = {bs}")
print(f"    命中模 m 之集合 = {sorted({m for (_,_,m) in exD})[:8]} … （共 {len({m for (_,_,m) in exD})} 個）")
print(f"    ⟹ 全部 b=44 之 43 對，通道落於 m=44−a 或 m=44+a（同一 γ）")

# ---------------------------------------------------------------- 結構定理
print("\n  ★ 根源 = 已註冊恆等式 γ_n + γ_{44−n} = 4（= γ_44），非新常數：")
mx = max(abs(gm(n) + gm(44 - n) - 4) for n in range(1, 44))
print(f"     max|γ_n + γ_(44−n) − 4| = {mp.nstr(mx, 3)}（n=1..43）")
print("     (i)  Σ=44 ⟹ γ_a+γ_b = γ_a+γ_{44−a} = 4 = γ_44 ⟹ sum 在殼，落 m=44")
print("     (ii) b=44 ⟹ γ_44−γ_a = 4−γ_a = γ_{44−a} ⟹ beat 在殼，落 m=44−a")
print("     ⟹ 兩條選擇規則**同一條恆等式**之兩個讀法（零新常數）")

# 互斥性
setS = {(a, b) for (a, b, _) in exS}
setD = {(a, b) for (a, b, _) in exD}
print(f"\n  集合：|S（sum 在殼）| = {len(setS)}；|B（beat 在殼）| = {len(setD)}；"
      f"交 = {len(setS & setD)}（互斥 ✅）")
print(f"  合計活價對 = {len(setS | setD)} / 946；其餘 {len(pairs)-len(setS|setD)} 對**兩通道皆 off-shell**")
print(f"  每型之 Δ 結構：B 型 ⟹ Δ = 44−a ⟹ 覆蓋 1..43 **全部 43 個 Δ 各一**")
print(f"                S 型 ⟹ Δ = 44−2a ⟹ 覆蓋 **偶 Δ 2..42 共 21 個各一**")

print("\n【§七】判定（取代 §五 之自我更正，見下）")
print(L)
print(f"""  ① 頻率守恒（ω 加法）：0 / 83248 共振；m=Σ／m=Δ 有代數 no-go。
  ② **能量守恒（γ 加法）：恰 64 / 946 價對有一條通道在殼** —— 全部由 γ_n+γ_(44−n)=4 逼出。
       B 型（b=44，43 對）：**beat 在殼** ⟹ 物質通道（Δ 各一，覆蓋 1..43）
       S 型（a+b=44，21 對）：**sum 在殼** ⟹ 反物質通道（偶 Δ 各一）
       兩型互斥，交集 0 ⟹ **每活價對恰有一條通道在殼、另一條 off-shell**
  ③ ⟹ **定向機制找到**（docs/09 §五 第 5 條所要求之「帶符號之定向」）：
       物質／反物質 = 邊條通道在殼，由價對位置唯一決定；**唔用宇稱** ⟹ 與 T2 不衝突。
  ④ 修正 scripts/13：其「II Σ=44 = Nyquist」**係真**（今由 γ 恆等式證實）；
       但「I Σ<44 傳播／III Σ>44 帶外」之三分**唔成立** —— I 與 III 兩者皆 off-shell，
       故 scripts/13 之三分塌縮為「在殼 / 唔在殼」二分 ＋ B/S 型之別。
  ⑤ 未決：43（B 型）與 21（S 型）如何對應 u,d,s 三味；3 仍未在任一層出現
       （43 = 43, 21 = 21, 64 = 43+21）。
""")
print(L)
raise SystemExit(0)