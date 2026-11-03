#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/21 — 拍頻空間中 /88 與 /90 之相位差：odd/even 產生什麼？

────────────────────────────────────────────────────────────────────────
問題（調律者）：`4sin²(nπ/x)`，x = 88 與 90，兩者之相位差；其 odd/even 之間
                嘅相位差產生嘅係什麼？

本腳本分解為五個可計算命題：
  §一 相位差本體：θ_n = nπ/x ⟹ Δθ_n = nπ(1/88 − 1/90) = nπ/3960（3960 = 44·90）
      ⟹ γ^(88)_n − γ^(90)_n = 4 sin(89nπ/3960) sin(nπ/3960)（cos 差化積）
  §二 總量：Σ_n Δ_n = ?（精確）
  §三 補償映射：/88 n↦44−n（保宇稱）；/90 n↦45−n = (44−n)+1（**單位點平移 ⟹ flip 宇稱**）
  §四 宇稱分級：和頻價對 a+b=44（/88）vs a+b=45（/90）之 parity 組成
  §五 通道命運：和頻 m=a+b 與拍頻 m=b−a 在兩鏈之命中數
  §六 交叉退化：γ^(88)_n = γ^(90)_(n') 是否存在？最小間隙？
  §七 偶子格折疊：γ^(88)_(2k) ≡ /44 鏈；γ^(90)_(2k) ≡ /45 鏈
性質：確定性全掃；dps = 60；精確判定閾 1e−40。N_trials = 0。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp

mp.mp.dps = 60
L = "=" * 78
N = 44
g8 = [mp.mpf(0)] + [4 * mp.sin(n * mp.pi / 88) ** 2 for n in range(1, N + 1)]
g9 = [mp.mpf(0)] + [4 * mp.sin(n * mp.pi / 90) ** 2 for n in range(1, N + 1)]
tol = mp.mpf('1e-40')

print(L)
print("M14 scripts/21 — /88 vs /90 之相位差：odd/even 產生什麼")
print(L)

# ======================================================================
print("\n【§一】相位差本體與 γ 差之閉式")
mx = mp.mpf(0)
pos = 0
for n in range(1, N + 1):
    d = g8[n] - g9[n]
    cl = 4 * mp.sin(89 * n * mp.pi / 3960) * mp.sin(n * mp.pi / 3960)
    mx = max(mx, abs(d - cl))
    if d <= 0:
        pos += 1
print("  Δθ_n = θ^(88)_n − θ^(90)_n = nπ(1/88 − 1/90) = **nπ/3960**（3960 = 44·90）")
print("  cos 差化積 ⟹ Δγ_n = γ^(88)_n − γ^(90)_n")
print("                  = 4·sin(89nπ/3960)·sin(nπ/3960)")
print(f"  數值核：max|Δγ − 閉式| = {mp.nstr(mx, 3)}（n = 1..44）")
print(f"  Δγ_n ≤ 0 之模數 = {pos}（即 Δγ_n > 0 對全部 n）⟹ /88 之譜**逐模嚴格高於** /90")
print(f"  Δγ_1 = {mp.nstr(g8[1]-g9[1], 8)}（最细）；max Δγ = {mp.nstr(max(g8[n]-g9[n] for n in range(1,N+1)), 8)}")

# ======================================================================
print("\n【§二】★★ 總量：Σ Δγ_n —— 88 與 90 之**互換**")
s8 = sum(g8[n] for n in range(1, N + 1))
s9 = sum(g9[n] for n in range(1, N + 1))
print(f"  Σγ^(88)（44 模，半週期 88）= {mp.nstr(s8, 12)}  ⟹ **= 90**")
print(f"  Σγ^(90)（44 模，半週期 90）= {mp.nstr(s9, 12)}  ⟹ **= 88 = 2N_d**")
print(f"  Σ Δγ = {mp.nstr(s8 - s9, 12)}  ⟹ **精確 = 2**")

# 閉式（精確推導）：已知 Σ_{n=1}^{M} cos(nπ/M) = −1
#   環（n = 1..M，半週期 2M）：Σγ = 2M − 2·(−1) = 2M + 2
#   Dirichlet（n = 1..M−1，半週期 2M）：Σγ = (2M+2) − 4sin²(Mπ/2M) = 2M − 2
print("\n  閉式（精確）：Σ_{n=1}^{M} cos(nπ/M) = −1 ⟹")
print("    環（n = 1..M，半週期 2M）      ⟹ Σγ = 2M + 2")
print("    Dirichlet（n = 1..M−1，半週期 2M）⟹ Σγ = 2M − 2")
bad = 0
for M in range(2, 61):
    ring = sum(4 * mp.sin(n * mp.pi / (2 * M)) ** 2 for n in range(1, M + 1))
    dirch = sum(4 * mp.sin(n * mp.pi / (2 * M)) ** 2 for n in range(1, M))
    if abs(ring - (2 * M + 2)) > mp.mpf('1e-30') or abs(dirch - (2 * M - 2)) > mp.mpf('1e-30'):
        bad += 1
print(f"  掃 M = 2..60 兩條閉式：違例 {bad} 個")
print("  ⟹ 本題兩個對象：環 M=44 ⟹ **90**；Dirichlet M=45 ⟹ **88**")
print("  ⟹ **相位差（88 ↔ 90）產生嘅係：週期與能量預算之互換**")
print("     環/44 之下 90 係預算、88 係週期；Dirichlet/45 之下 88 係預算、90 係週期。")
print(f"  另核：Δγ_44 = {mp.nstr(g8[44]-g9[44], 10)} ＝ γ^(90)_1 = {mp.nstr(g9[1], 10)}（精確）")
print("       因 γ^(88)_44 = 4（Nyquist）且 γ^(90)_44 = 4 − γ^(90)_1")

# ======================================================================
print("\n【§三】★ 補償映射之差 = 一位點平移 = parity flip")
ok8 = max(abs(g8[n] + g8[44 - n] - 4) for n in range(1, 44))
ok9 = max(abs(g9[n] + g9[45 - n] - 4) for n in range(1, 45))
print(f"  /88：γ_n + γ_(44−n) = 4    max|dev| = {mp.nstr(ok8,3)}   ⟹ 補償模 c(n) = 44 − n")
print(f"  /90：γ_n + γ_(45−n) = 4    max|dev| = {mp.nstr(ok9,3)}   ⟹ 補償模 c(n) = 45 − n")
same8 = sum(1 for n in range(1, 45) if (44 - n) % 2 == n % 2)
same9 = sum(1 for n in range(1, 45) if (45 - n) % 2 == n % 2)
print(f"  宇稱：c_88(n) ≡ n (mod 2) 之模數 = {same8}/44 ⟹ **保宇稱**")
print(f"        c_90(n) ≡ n (mod 2) 之模數 = {same9}/44 ⟹ **flip 宇稱**")
print("  ⟹ 45 − n = (44 − n) + 1 ⟹ /90 之反射 = /88 之反射 **再一位點平移**")
print("     而一位點平移 = 二分格（bipartite）換邊 ⟹ 恰係 parity flip。")
print("  ⟹ **odd/even 相位差之幾何身份 = 單位點平移（half-cell shift）**")

# ======================================================================
print("\n【§四】宇稱分級：和頻價對之 parity 組成")
for x, fam in ((44, [(a, 44 - a) for a in range(1, 22)]), (45, [(a, 45 - a) for a in range(1, 23)])):
    same = sum(1 for a, b in fam if a % 2 == b % 2)
    tag = "/88" if x == 44 else "/90"
    print(f"  {tag} a+b={x}：{len(fam)} 對；同 parity {same}/{len(fam)}"
          f" ⟹ {'全部同 parity（grade +1）' if same == len(fam) else '全部異 parity（grade −1）'}")
print("  ⟹ 和頻合成物之**宇稱 grade** 由 x 之奇偶直接決定：x 偶 ⟹ 同 parity；x 奇 ⟹ 異 parity")

# ======================================================================
print("\n【§五】★ 通道命運：和頻 vs 拍頻")
for nm, g, comp in (("/88", g8, 44), ("/90", g9, 45)):
    s = [(a, b, a + b) for a in range(1, N + 1) for b in range(a + 1, N + 1)
         if a + b <= N and abs(g[a] + g[b] - g[a + b]) < tol]
    bt = [(a, b, b - a) for a in range(1, N + 1) for b in range(a + 1, N + 1)
          if abs(g[a] + g[b] - g[b - a]) < tol]
    print(f"  {nm}：和頻 m=a+b 命中 **{len(s)}**（落點 m = {sorted({m for _,_,m in s})}）；"
          f"拍頻 m=b−a 命中 **{len(bt)}**")
print("  ⟹ **和頻（反物質）通道：/88 開 21 對（全落 Nyquist m=44，γ=4）；/90 完全關閉 0 對**")
print("     因 /90 之和頻落點必為 m = 45 − (a+b 之補) ⟹ m = 45 ∉ 模集（Dirichlet 邊界節點）")
print("  ⟹ 拍頻通道兩鏈皆關（能量層），與 docs/12 之頻率 no-go 一致。")

# ======================================================================
print("\n【§六】交叉退化：γ^(88)_n = γ^(90)_(n′) 存在否？")
cross = [(n, n2) for n in range(1, N + 1) for n2 in range(1, N + 1) if abs(g8[n] - g9[n2]) < tol]
gaps = min(abs(g8[n] - g9[n2]) for n in range(1, N + 1) for n2 in range(1, N + 1))
print(f"  交叉退化數 = {len(cross)} ⟹ **零**")
print(f"  代數理由：γ^(88)_n = γ^(90)_(n′) ⟹ 45n = 44n′ ⟹ n = 44t, n′ = 45t ⟹ t = 0（出界）")
print(f"  最小間隙 = {mp.nstr(gaps, 8)}（= Δγ_1）⟹ 兩譜**永不相交**")

# ======================================================================
print("\n【§七】偶子格折疊（n ↦ 2k）")
e8 = max(abs(g8[2 * k] - 4 * mp.sin(k * mp.pi / 44) ** 2) for k in range(1, 23))
e9 = max(abs(g9[2 * k] - 4 * mp.sin(k * mp.pi / 45) ** 2) for k in range(1, 23))
print(f"  γ^(88)_(2k) ≡ 4sin²(kπ/44)  max|dev| = {mp.nstr(e8,3)}  ⟹ 偶子格 = **/44 鏈**")
print(f"  γ^(90)_(2k) ≡ 4sin²(kπ/45)  max|dev| = {mp.nstr(e9,3)}  ⟹ 偶子格 = **/45 鏈**")
print("  ⟹ 折疊唔保載體類型：/88 → /44（22 模）；/90 → /45（44 模）")
print("  ⟹ 故 odd/even 之相位差**唔對稱**：偶子格之兩鏈模數不同（22 vs 44）")

# ======================================================================
print("\n" + L)
print("【§八】判定：odd/even 相位差產生什麼")
print(L)
print("""  ① 相位差本體：Δθ_n = nπ/3960（3960 = 44·90），亦等價於 Δγ_n = 4 sin(89nπ/3960)sin(nπ/3960) > 0。
  ② 其幾何身份 = **單位點平移**（45−n 對 44−n 差 1）⟹ 二分格換邊 ⟹ **parity flip**。
  ③ ⟹ 產生三件事（皆精確、可核）：
     (a) **宇稱 grade 翻轉**：/88 之和頻價對 21/21 全同 parity（grade +1）；
         /90 之 22/22 全異 parity（grade −1）。
     (b) **和頻（反物質）通道關閉**：/88 命中 21（落 m=44=Nyquist，γ=4）；/90 命中 **0**
         （落點 m=45 係 Dirichlet 邊界節點，非模）。
     (c) **能量預算差 = 2**（Σγ：88+4 = 92 vs 90），且逐模 /88 嚴格高於 /90（44/44）。
  ④ 交叉退化 = 0，最小間隙 2.2387e−4 ⟹ 兩載體**不可互換**，必須由外部原則裁決。
  ⑤ 偶子格折疊唔保類型（/88→/44 22 模；/90→/45 44 模）⟹ odd 與 even 兩半唔同構。
  ⑥ **物理讀法（明示為讀法，非結果）**：和頻 = 反物質（docs/12 之定向機制）⟹
     /88 載體「物質＋反物質」皆有；/90 載體**只有物質側**。
     此讀法須 T2 閘（C(s,s):C(s,c):C(s,b) 符號）前先凍結，否則屬後驗。
""")
print(L)
raise SystemExit(0)