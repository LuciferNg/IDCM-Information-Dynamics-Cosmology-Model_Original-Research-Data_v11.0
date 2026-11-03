#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/22 — 相位差之二階結構：能量預算「缺陷」δ = Σγ − x 之三分律 ＋ 2-adic 折疊深度

────────────────────────────────────────────────────────────────────────
承接 scripts/21（docs/20）：/88 與 /90 之差 = 一位點平移（parity flip）＋ 88/90 互換。
本輪追問：該相位差在**每一層**（dyadic fold）上留下什麼？

三條閉式（待驗）：
    ring（x 偶，M = x/2）           ⟹ Σγ = x + 2   ⟹ δ = **+2**
    Dirichlet（x 偶，M = x/2 − 1）  ⟹ Σγ = x − 2   ⟹ δ = **−2**
    Dirichlet（x 奇，M = (x−1)/2）  ⟹ Σγ = x       ⟹ δ = **0**
折疊：γ^(x)_(2k) = γ^(x/2)_k（x、M 皆偶時可摺）⟹ 2-adic 深度 = v₂(x)。

性質：確定性掃盡；dps = 50；判定閾 1e−30。N_trials = 0。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp

mp.mp.dps = 50
L = "=" * 78
tol = mp.mpf('1e-30')


def spec(x, M):
    return [4 * mp.sin(n * mp.pi / x) ** 2 for n in range(1, M + 1)]


print(L)
print("M14 scripts/22 — 相位差之二階結構：δ = Σγ − x 之三分律 ＋ 2-adic 折疊深度")
print(L)

# ======================================================================
print("\n【§一】三條閉式：掃 x = 3..200")
bad = []
for x in range(3, 201):
    if x % 2 == 0:
        # ring：M = x/2
        s = sum(spec(x, x // 2))
        if abs(s - (x + 2)) > tol:
            bad.append(("ring", x, mp.nstr(s, 10)))
        # Dirichlet：M = x/2 − 1
        s = sum(spec(x, x // 2 - 1))
        if abs(s - (x - 2)) > tol:
            bad.append(("dir-even", x, mp.nstr(s, 10)))
    else:
        # Dirichlet：M = (x−1)/2
        s = sum(spec(x, (x - 1) // 2))
        if abs(s - x) > tol:
            bad.append(("dir-odd", x, mp.nstr(s, 10)))
print(f"  違例 {len(bad)} 個（x = 3..200，三類全掃）")
for t in bad[:5]:
    print("   ", t)
print("""  閉式確立：
    ring（x 偶，M = x/2）           Σγ = x + 2   δ = +2
    Dirichlet（x 偶，M = x/2 − 1）  Σγ = x − 2   δ = −2
    Dirichlet（x 奇，M = (x−1)/2）  Σγ = x       δ =  0""")

print("\n  ★ 三分律之意義：δ **由 (邊界類型, x 之奇偶) 完全決定**")
print("     ring ⟹ 恆 +2（與 x 無關）")
print("     Dirichlet ⟹ x 偶則 −2、x 奇則 0 ⟹ **過 odd/even 界線時 δ 由 −2 跳 0**")

# ======================================================================
print("\n【§二】本題兩載體：δ 實測")
for lab, x, M in (("/88 ring", 88, 44), ("/90 Dirichlet", 90, 44)):
    s = sum(spec(x, M))
    print(f"  {lab:<15s} x={x:3d} M={M:2d}  Σγ={mp.nstr(s, 12):>8s}  δ = Σγ − x = {mp.nstr(s - x, 6)}")
print("  ⟹ /88 之 δ = +2；/90 之 δ = −2 ⟹ 兩者相隔 **4**（= 2B）")
print("  ⟹ 已註冊之 Σγ=90 正是 /88 之 x+2；/90 之 88 正是 x−2")

# ======================================================================
print("\n【§三】★ 2-adic 折疊鏈（odd/even 之迭代結構）")


def chain(x, M, lab):
    print(f"  {lab}：")
    lv = [(x, M, sum(spec(x, M)))]
    while x % 2 == 0 and M % 2 == 0:
        x2, M2 = x // 2, M // 2
        dev = max(abs(spec(x, M)[2 * k - 1] - spec(x2, M2)[k - 1]) for k in range(1, M2 + 1))
        lv.append((x2, M2, sum(spec(x2, M2))))
        print(f"    fold {x:>3d}(M={M:>2d}) → {x2:>3d}(M={M2:>2d})   max|dev| = {mp.nstr(dev, 3)}"
              f"   Σγ = {mp.nstr(sum(spec(x2, M2)), 10)}   δ = {mp.nstr(sum(spec(x2, M2)) - x2, 4)}")
        x, M = x2, M2
    print(f"    終止：x = {x}（{'奇' if x % 2 else 'M 奇'}）⟹ 層數 = {len(lv)}，折疊數 = {len(lv) - 1}")
    return lv


lv8 = chain(88, 44, "/88 ring")
print()
lv9 = chain(90, 44, "/90 Dirichlet")
print(f"\n  v₂(88) = {(88 & -88).bit_length() - 1} = 3 ； v₂(90) = {(90 & -90).bit_length() - 1} = 1")
print("  ⟹ **/88 有 3 層（44→22→11 模），/90 只有 2 層（44→22）**")
print("  ⟹ odd/even 相位差之迭代身份 = **2-adic 深度之差（3 層 vs 2 層）**")
print(f"  ⟹ odd core：/88 → {88 // 8}（質數）；/90 → {90 // 2} = 4·11 + 1 = {4 * 11 + 1}")

print("\n【§三b】折疊條件 ≡ x mod 4（掃 x = 6..200）")
bad4 = 0
for x in range(6, 201, 2):
    Mring, Mdir = x // 2, x // 2 - 1
    f_ring = (x % 2 == 0 and Mring % 2 == 0)
    f_dir = (x % 2 == 0 and Mdir % 2 == 0)
    if f_ring != (x % 4 == 0):
        bad4 += 1
    if f_dir != (x % 4 == 2):
        bad4 += 1
print(f"  違例 {bad4} 個 ⟹ ring 可摺 ⟺ x ≡ 0 (mod 4)；Dirichlet 可摺 ⟺ x ≡ 2 (mod 4)")
print("  ⟹ 因 ring 之 M = x/2、Dirichlet 之 M = x/2 − 1 差 1，故**兩者之可摺性恰好相反**")
print("     ⟹ /88 ≡ 0 (mod 4) 走 ring（可摺 2 次）；/90 ≡ 2 (mod 4) 走 Dirichlet（可摺 1 次）")
print("  ⟹ **相位差之最深不可約形式：88 ≡ 0 (mod 4) 對 90 ≡ 2 (mod 4)**")

# ======================================================================
print("\n【§四】各層之能量預算序列")
print("  /88：", " → ".join(f"{mp.nstr(s, 8)}(x={x})" for x, M, s in lv8))
print("  /90：", " → ".join(f"{mp.nstr(s, 8)}(x={x})" for x, M, s in lv9))
print("  ⟹ /88 之層預算 ≡ x+2（+2 缺陷恆定）；/90 第 2 層 x=45（奇）⟹ δ = 0")
print("  ⟹ 兩載體在共同層（M=22）之預算：46 vs 45，差 **1**")

# ======================================================================
print("\n" + L)
print("【§五】判定")
print(L)
print("""  ① **三分律**（掃 x=3..200 違例 0）：δ = Σγ − x ∈ {+2, −2, 0}
       ring ⟹ +2（恆）；Dirichlet x 偶 ⟹ −2；Dirichlet x 奇 ⟹ 0
     ⟹ **odd/even 相位差之最簡不可約內容 = 呢個 ±2／0 之缺陷分類。**
  ② 本題：/88 之 δ = +2、/90 之 δ = −2 ⟹ 相隔 4 = 2B。
  ③ 折疊：/88 三層（88→44→22，odd core 11）；/90 兩層（90→45，odd core 45 = 4·11+1）。
     ⟹ 相位差之迭代身份 = **2-adic 深度 3 vs 2**。
  ④ 共同層（M=22）預算 46 vs 45 ⟹ 差 1（每上一層差翻倍：2 → 1）。
  ⑤ ⚠ 「3 層」是否即 `N_gen = 3`：**不作主張**。依 docs/19 單標籤定理，3 只能係**計數**；
     本輪給出一個新的計數身份 `v₂(88) = 3`（層數 3），但與世代之對應須**先凍結再盲測**。
  ⑥ 新計數身份（記錄，不裁決）：88 = 2³·11 ⟹ 折疊深度 3；90 = 2·3²·5 ⟹ 深度 2。
""")
print(L)
raise SystemExit(0)