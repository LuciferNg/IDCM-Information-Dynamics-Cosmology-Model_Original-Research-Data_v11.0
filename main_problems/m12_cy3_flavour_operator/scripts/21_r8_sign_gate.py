#!/usr/bin/env python3
"""R8 — 重-重相關路線之可行性閘：核之符號結構 ＋ 站點指派前置審計

背景（R7 §九）：R7 證明加性味荷（任何標籤）全死，唯一有結構支撐之殘路係
   Δm = （2 體對荷）⊕（兩個重旁觀者之成對相關能 C(f1,f2)）
且 C 必須給出 (s,s) < 0、(s,c) > 0、(s,b) ≈ 0（R6/R7 之偏移符號）。

預註冊判準（跑前凍結）
--------------------------------------------------
G1 正性檢定：44-site 鏈 Dirichlet Green 函數 G(x,y) = Σ_n φ_n(x)φ_n(y)/γ_n。
   若 min G > 0（全正）⟹ G 為 sign-definite ⟹ 任何以 G 為核之相關能
   C(f1,f2) = G(x_{f1}, x_{f2}) 對所有味對**同號** ⟹ **範疇判死**（不能生 (s,s)<0 而 (s,c)>0）。
   若出現負值 ⟹ 該核可生異號，路線在原理上存活。
G2 正性之來源：γ_n > 0 且 φ_n(x)φ_n(y) 之符號由 sin 決定；
   對 x,y 同側（同半鏈）全正；跨中點可負。故須**分區**報 sign census。
G3 站點指派前置審計：v11 已註冊之味→模式指數指派有幾多？c、b 有無？
   若無 ⟹ 任何 C(s,c)、C(s,b) 之評估都係後驗指數選擇 ⟹ 依禁區規則不可作判據。
禁用：不得為 c/b 臨場挑指數；不得調整核以遷就符號。
"""
import pathlib
import mpmath as mp

mp.mp.dps = 40
N = 44


def gamma(n):
    """IDCM 鏈譜 γ_n = 4 sin²(πn/88)（44-site 雙鏈之半譜）"""
    return 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2


def phi(n, j):
    """44-site 鏈 Dirichlet 標準正交基（DST-I）：√(2/45)·sin(πnj/45)"""
    return mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))


print("=" * 92)
print("R8 — 重-重相關路線之可行性閘（核符號結構 ＋ 指派前置審計）")
print("=" * 92)

print("\n[G1] 44-site 鏈 Dirichlet Green 函數之符號普查")
G = [[mp.mpf(0)] * (N + 1) for _ in range(N + 1)]
for i in range(1, N + 1):
    for j in range(1, N + 1):
        s = mp.mpf(0)
        for n in range(1, N + 1):
            s += phi(n, i) * phi(n, j) / gamma(n)
        G[i][j] = s
pos = neg = zero = 0
mn = mp.mpf(10) ** 9
mnpair = None
for i in range(1, N + 1):
    for j in range(1, N + 1):
        v = G[i][j]
        if v > 1e-30:
            pos += 1
        elif v < -1e-30:
            neg += 1
        else:
            zero += 1
        if v < mn:
            mn, mnpair = v, (i, j)
print("   總元素 %d：正 %d、負 %d、零 %d" % (N * N, pos, neg, zero))
print("   最小值 = %s（於 i,j = %s）；最大值 = %s"
      % (mp.nstr(mn, 8), mnpair, mp.nstr(max(max(r[1:]) for r in G[1:]), 8)))
print("   對角 G(1,1) = %s ；G(1,44) = %s ；G(22,23) = %s"
      % (mp.nstr(G[1][1], 8), mp.nstr(G[1][44], 8), mp.nstr(G[22][23], 8)))
print("   ⟹ %s" % ("**全正 ⟹ sign-definite ⟹ 範疇判死** 🔴" if neg == 0
                  else "含負值 ⟹ 原理上可生異號 🟡"))

print("\n[G2] 分區 sign census（同半鏈 vs 跨中點）")
same = sum(1 for i in range(1, 23) for j in range(1, 23) if G[i][j] < 0)
cross = sum(1 for i in range(1, 23) for j in range(23, 45) if G[i][j] < 0)
print("   上半鏈×上半鏈 22×22：負值 %d" % same)
print("   上半×下半 22×22：負值 %d" % cross)
print("   核構造：G = Σ_n φ_n(i)φ_n(j)/γ_n，γ_n > 0 且 φ 之符號由 sin(πnj/45) 決定")
print("   ⟹ 對稱、正定、以 1/γ_n 加權 ⟹ 屬 Stieltjes 型；符號完全由基底之交疊決定")

print("\n[G3] 站點指派前置審計（v11 已註冊之味 → 模式指數）")
REG = [("u", "2", "O209 d14 / O233（已裁：觀測輸入角色，譜指派 n=2）"),
       ("d", "3", "同上（n=3）"),
       ("s", "4 / 5 / 7（三讀未定）", "R4-O4 顯示三讀全死 ⟹ 非指派選擇問題"),
       ("c", "**無**", "v11 無註冊指派"),
       ("b", "**無**", "v11 無註冊指派")]
for f, k, note in REG:
    print("   %-2s → %-22s %s" % (f, k, note))
print("   ⟹ 要評估 C(s,c)、C(s,b) 之數值，須先有 c、b 之站點指數；")
print("     現時**不存在** ⟹ 任何此類評估即屬後驗指數選擇，依禁區規則不可作判據。")

print("\n" + "=" * 92)
print("[匯總 — 只報事實]")
print("   G1 Green 核：負元素 = %d ⟹ %s" % (neg, "sign-definite（路線死）" if neg == 0 else "可生異號"))
print("   G3 指派：c、b 均無註冊指數 ⟹ 相關能數值評估**前置受阻**")
print("   ⟹ 相關能路線之狀態：%s" % ("🔴 範疇判死" if neg == 0 else "🟡 前置受阻（缺註冊指派）"))
print("=" * 92)