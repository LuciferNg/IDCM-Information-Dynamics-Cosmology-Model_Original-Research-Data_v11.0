#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/13 — 扇區之**通道存續**分類：beat / sum 之帶內、Nyquist、帶外

⚠ **讀法已更正（2026-10-31，見 docs/13）**：
   本腳本以 Σ=a+b 相對 N_d=44 作三分，但 Σ 係**動量指數**，唔係頻率指數。
   scripts/14 之真實頻率掃描：ω 加法 0/83,248 共振 ⟹ I 與 III 兩類**皆 off-shell**；
   而 II 類（Σ=44，21 對）經 γ 恆等式確認**真在殼**（γ_a+γ_(44−a)=4=γ_44）。
   ⟹ 三分**塌縮**為「在殼／唔在殼」＋「beat 型／sum 型」。
   程式與 log 保留不動（append-only）；正確結構見 docs/12、docs/14。

────────────────────────────────────────────────────────────────────────
起點（T2 之代數核心，docs/10 §四）：
  cos(A)cos(B) = ½[cos(A−B) + cos(A+B)] ⟹ beat 與 sum 之係數**恆等**（+½/+½）。
  ⟹ 物質／反物質之差**不可能**來自係數大小或符號，只可能來自
     「該通道係唔係載體之合法模」＝**帶內／帶外**。
  而 docs/09 §五 第 5 條正要求「必須帶符號之定向」——本腳本檢定此門。

對象：44-站 Dirichlet 鏈，模 n = 1..44，γ_n = 4sin²(nπ/88) 嚴格遞增。
  價對 (a,b)：a<b。兩條派生通道
    beat  Δ = b − a ∈ 1..43   （恆在帶內；Δ = 44 不可達）
    sum   Σ = a + b ∈ 3..87   （可達 44 = 44-鏈 Nyquist；可 > 44 = 帶外）
  結構特殊模（skill 已註冊）：n=22 (γ=2, v_g=c)  ·  n=44 (γ=4, Nyquist, v_g=0)
  ⚠ 故 Σ = 44 之 sum 通道係**駐波（非傳播）**，與 Σ<44（傳播）／Σ>44（不存在）三分。

性質：確定性計數；零觀測量；N_trials = 0。
⚠ **誠實標記**：本分類**並非** docs/10 §一 預先凍結之清單（P1–P8）成員 ——
   佢係在「尋找 3 類」之動機下**後驗**發現 ⟹ 依禁區規則，只能作假設，
   必須經預註冊盲測方可升格。本腳本只報結構，不作「命中 3」之宣稱。
────────────────────────────────────────────────────────────────────────
"""
import math
from collections import Counter
import mpmath as mp

mp.mp.dps = 30
L = "=" * 78
ND = 44


def gam(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print(L)
print("M14 scripts/13 — 通道存續分類（帶內／Nyquist／帶外）")
print(L)

# ======================================================================
print("\n【§一】地基：γ 嚴格遞增 ＋ 特殊模")
gs = [gam(n) for n in range(1, ND + 1)]
inc = all(gs[i] < gs[i + 1] for i in range(ND - 1))
print(f"  γ_n 於 1..44 嚴格遞增？{inc}")
print(f"  γ_22 = {mp.nstr(gam(22), 10)}（=2，γ_22 = 4−γ_22 ⟹ 自對偶, v_g = c）")
print(f"  γ_44 = {mp.nstr(gam(44), 10)}（=4，Nyquist, v_g = 0）")
print(f"  Δ 之可達範圍 = 1..{ND-1}（Δ=44 需 a=0，不可達）")
print(f"  Σ 之可達範圍 = {1+2}..{43+44}（從未到 88）")

# ======================================================================
print("\n【§二】價對之分類：sum 通道之命運")
pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]
cls = {"I  Σ<44 傳播": [], "II Σ=44 Nyquist(駐波)": [], "III Σ>44 帶外": []}
for (a, b) in pairs:
    s = a + b
    if s < ND:
        cls["I  Σ<44 傳播"].append((a, b))
    elif s == ND:
        cls["II Σ=44 Nyquist(駐波)"].append((a, b))
    else:
        cls["III Σ>44 帶外"].append((a, b))
tot = 0
for k, v in cls.items():
    tot += len(v)
    ex = v[:5]
    print(f"  {k:<24s} 價對 = {len(v):>4d}   例 {ex}")
print(f"  合計 = {tot}（應 = C(44,2) = {len(pairs)}）⟹ {'✅' if tot == len(pairs) else '🔴'}")
print(f"\n  算術副產物（只報，不解釋）：")
print(f"    |I| = {len(cls['I  Σ<44 傳播'])} = {round(len(cls['I  Σ<44 傳播'])**0.5)}² ?"
      f"  |II| = {len(cls['II Σ=44 Nyquist(駐波)'])}  |III| = {len(cls['III Σ>44 帶外'])} = "
      f"{round(len(cls['III Σ>44 帶外'])**0.5)}² ?")
print(f"    Σ=44 之價對：a+b=44, a<b ⟹ a=1..21 ⟹ 恰 {len(cls['II Σ=44 Nyquist(駐波)'])} 對"
      f"（21 = (N_d−2)/2 = F₈）")

# ======================================================================
print("\n【§三】對載體唯一對稱之承載性（reflection n ↦ 標籤不變）")
# φ_n(45−p) = sin(nπ − nπp/45) = (−1)^(n+1) sin(nπp/45)
# ⚠ 容差：dps=30 ⟹ sin(nπ) 之浮點殘餘 ~1e−28 ⟹ 容差取 1e−20（首跑誤用 1e−30 致 1170 假違例）
bad = 0
worst = mp.mpf(0)
for n in range(1, ND + 1):
    for p in range(1, ND + 1):
        lhs = mp.sin(n * mp.pi * (ND + 1 - p) / (ND + 1))
        rhs = (-1) ** (n + 1) * mp.sin(n * mp.pi * p / (ND + 1))
        dev = abs(lhs - rhs)
        worst = max(worst, dev)
        if dev > mp.mpf('1e-20'):
            bad += 1
print(f"  φ_n(45−p) = (−1)^(n+1) φ_n(p)：違例 {bad} / {ND*ND}；max|dev| = {mp.nstr(worst, 3)}")
print("  ⟹ 反射作用於模**標籤恆等**（只改符號）⟹ (a,b)、Δ、Σ 全部 reflection-invariant")
print("  ⟹ 本分類唔依賴座標frame ⟹ canonical（對 44-鏈之唯一對稱而言）")
# D_88 之反射/旋轉亦保 Δ：τ(a,b) = (88−b, 88−a) ⟹ Δ' = (88−a)−(88−b) = b−a
print("  另一檢：D_88 之 τ(a,b)=(88−b,88−a) ⟹ Δ' = b−a = Δ（精確恆等，無需計算）")

# ======================================================================
print("\n【§四】每 Δ 之三類分佈（檢查有無 Δ 之通道結構特殊）")
tab = {}
for (a, b) in pairs:
    d, s = b - a, a + b
    key = "I" if s < ND else ("II" if s == ND else "III")
    tab.setdefault(d, Counter())[key] += 1
full3 = [d for d, c in tab.items() if len(c) == 3]
only1 = [d for d, c in tab.items() if len(c) == 1]
print(f"  Δ 值數 = {len(tab)}；三類皆出現之 Δ = {len(full3)}；只出現一類之 Δ = {len(only1)}")
print(f"    Δ 只一類者：{only1}")
print(f"  {'Δ':>3s} {'I':>4s} {'II':>4s} {'III':>4s}")
for d in range(1, ND):
    c = tab[d]
    print(f"  {d:>3d} {c.get('I',0):>4d} {c.get('II',0):>4d} {c.get('III',0):>4d}")

# ---------------------------------------------------------------- §四b
print("\n【§四b】Nyquist 代表元：每個**偶** Δ 恰有一個 Σ=44 之價對？")
even_ok = all(tab[d].get("II", 0) == 1 for d in range(2, ND, 2))
odd_ok = all(tab[d].get("II", 0) == 0 for d in range(1, ND, 2))
print(f"  偶 Δ (2..42, 共 21 個)：II 數目全 = 1 ？{even_ok}")
print(f"  奇 Δ (1..43, 共 22 個)：II 數目全 = 0 ？{odd_ok}")
print("  公式：Σ=44 ⟹ a = (44−Δ)/2, b = (44+Δ)/2 ⟹ 偶 Δ 唯一解，奇 Δ 無整數解")
print(f"  逐例：Δ=2 ⟹ {((44-2)//2, (44+2)//2)} ；Δ=42 ⟹ {((44-42)//2, (44+42)//2)}")
print("  ⟹ **每個偶 Δ 之軌道內恰有一個駐波（Nyquist）成員** ⟹ 係 canonical representative")
print("  ⟹ 奇 Δ 之軌道**冇**駐波成員 ⟹ 兩族之通道結構不同（但只 2 族，非 3）")

# ======================================================================
print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ① 三類之價對數 = {len(cls['I  Σ<44 傳播'])} / {len(cls['II Σ=44 Nyquist(駐波)'])} / {len(cls['III Σ>44 帶外'])}
     ⟹ **存在一個 canonical（reflection-invariant）、零參數之 3 類分割**
  ② 但三類**大小極不等**（{len(cls['I  Σ<44 傳播'])} : {len(cls['II Σ=44 Nyquist(駐波)'])} : {len(cls['III Σ>44 帶外'])}）
     ⟹ 與「3 個世代等價」之讀法**不相容**；此 3 類係「通道命運」，非世代。
  ③ 本分類**不在** docs/10 §一 之凍結清單 ⟹ **後驗**，只可作假設，須預註冊盲測。
  ④ 對 docs/09 §五 第 5 條（必須帶符號之定向）之意義：
     cos·cos 之係數恆等 ⟹ 定向唔可能來自係數；本分類指出**唯一尚存之機制**
     = 「該通道係唔係合法模」（帶內傳播／Nyquist 駐波／帶外）。
     此機制**唔用宇稱**，故與 T2 並不衝突（T2 只封宇稱路線）。
  ⑤ 未決：三類之中，「III 帶外」意味 sum 通道不存在 ⟹ 該價對只有 beat（物質）通道；
     「I」兩通道俱在；「II」sum 為駐波（不傳播，等效半通道）。
     若要由三類生成 ℤ₂ 電荷，須額外一個「選邊」規則 —— 未找到。
""")
print(L)
raise SystemExit(0)