#!/usr/bin/env python3
"""30 — 外部模型主張裁決：Gemini「δ = 拓撲代數空間之非對角交叉干涉本徵投影」

被裁決之命題（用戶轉貼，屬**未驗證散文**，非工具鑄造之 claim）
--------------------------------------------------
G1（消極）：「R13 終結了『用線性／非線性關聯核之有限差分生成 δ_x』這整個範式」。
G2（積極）：「δ 不是空間／網路站點上的位移差分，而是拓撲代數空間裏非對角交叉干涉（Non-diagonal
            Phase Interference）之直接本徵投影」。

本腳本把 G2 化為**可計算**之三種讀法，逐一核（不引入任何新參數）：
  R0  a(k,x) = ⟨k,x|A|k,x⟩，A = T_c ⊗ I                     （T 之對角元 ⟹ 恆零）
  R1  a(k,x) = T_c(k,x)（以 hopping 矩陣元作對荷）          （非對角但非同 x）
  R2  δ_x = ⟨3,x|A|2,x⟩（真非對角），A = T_c ⊗ I            （同 x 之非對角矩陣元）
  R3  δ_x = ⟨3,x|A|2,x⟩，A = T_c ⊗ D_t（D_t 對角）          （張量分解之一般形）
  R4  δ_x = ⟨3,x|A|2,x⟩，A = T_c ⊗ T_t（雙 hopping）        （真「交叉」干涉）
另核 G1 之範圍：確認 R13 實際殺咗邊個核族（正核／Green／截斷），**未**殺邊個（帶號、非線性乘積）。

判準（先寫死）：任何讀法若 |z(ρ_d)| ≥ 3 或 |z(ρ_u)| ≥ 3 ⟹ 該讀法判死。
恆等零亦判死（δ ≡ 0 ⟹ 無分裂）。
"""
import mpmath as mp

mp.mp.dps = 40
N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

print("=" * 96)
print("30 — 外部模型主張裁決（Gemini：非對角交叉干涉）")
print("=" * 96)
print("  γ_2 = %s  γ_3 = %s  γ_4 = %s" % (mp.nstr(gam(2), 8), mp.nstr(gam(3), 8), mp.nstr(gam(4), 8)))

# ---- 44-site Dirichlet 鏈之 hopping 算子 T（adjacency，無自環）
T = {}
for i in range(1, N + 1):
    for j in range(1, N + 1):
        T[(i, j)] = mp.mpf(1) if abs(i - j) == 1 else mp.mpf(0)

print("\n  鏈 hopping 核實：T(3,2) = %s（非對角、常數）；T(x,x) = %s（無自環）"
      % (T[(3, 2)], T[(22, 22)]))


def judge(tag, du, dd, ds, note=""):
    if dd == 0 and ds == 0 and du == 0:
        print("   %-46s δ ≡ 0 ⟹ **恆等零，判死** %s" % (tag, note))
        return False
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-46s ρ_d = %-12s（z %+8.2f） ρ_u = %-12s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return ok


print("\n[讀法 R0] a(k,x) = ⟨k,x|T_c⊗I|k,x⟩ = T_c(k,k) = 0 ∀k ⟹ δ ≡ 0")
judge("R0 對角元（T 無自環）", mp.mpf(0), mp.mpf(0), mp.mpf(0), "（恆等零）")

print("\n[讀法 R1] a(k,x) = T_c(k,x)：以 hopping 矩陣元作對荷，δ_x = T(3,x) − T(2,x)")
judge("R1 hopping 矩陣元作對荷",
      T[(3, KU)] - T[(KU, KU)], T[(3, KD)] - T[(KU, KD)], T[(3, KS)] - T[(KU, KS)],
      "（僅 x=2,3,4 非零 ⟹ 三點定三值）")

print("\n[讀法 R2] δ_x = ⟨3,x|T_c⊗I|2,x⟩ = T_c(3,2) = 常數 （真非對角，同 x）")
c = T[(3, 2)]
print("   δ_u = δ_d = δ_s = %s ⟹ δ 與 x 無關" % c)
print("   ⟹ 比值未定義（分母同值、分子同值 ⟹ ρ ≡ 1），且**與實測 δ 之 2.5× 變動矛盾** ⟹ 判死")

print("\n[讀法 R3] δ_x = ⟨3,x|T_c⊗D_t|2,x⟩ = T_c(3,2)·d(x) ⟹ ρ 只由 d(x) 決定")
print("   即：δ_x ∝ d(x)，與 x 之關係**等於任意譜函數 f(γ_x)** ⟹ 屬 R10 已測之 18 形式族")
judge("R3 ∝ d(x)（取 d = γ 為代表）",
      gam(KU), gam(KD), gam(KS), "（= R10 之 f=γ、k_s=4：z +0.41 但 6/18、P=0.115 零證據）")
print("   註：d 可換任何函數 ⟹ 本讀法之全部自由度 = R10 之形狀族；R10 已判零證據、R12b 已判 TI 全死。")

print("\n[讀法 R4] δ_x = ⟨3,x|T_c⊗T_t|2,x⟩ = T_c(3,2)·T_t(x,x)")
print("   任何圖之 adjacency 皆無自環 ⟹ T_t(x,x) = 0 ∀x ⟹ **δ ≡ 0**")
print("   核實（44 鏈）：T(5,5) = %s ；T(22,22) = %s ；T(44,44) = %s"
      % (T[(5, 5)], T[(22, 22)], T[(44, 44)]))
judge("R4 雙 hopping（交叉干涉之字面讀法）", mp.mpf(0), mp.mpf(0), mp.mpf(0), "（恆等零）")

print("\n" + "=" * 96)
print("[G1 範圍核實] R13 實際殺咗邊個核族")
print("   已殺：正定／PSD 核（R7-A1 範疇）、平移不變核（R12b-S1，15.58×）、")
print("         單調 Green 族 44/88（R11/R12）、T_pipe 截斷全 n*（R13）")
print("   未殺：①帶號／staggered 核（(−1)^{i+j}G，即 checkerboard 對偶）—— 即 S3")
print("         ②核之非線性乘積（G·G、G²、off-diagonal 乘積）—— 未形式化")
print("         ③非張量分解之 A（真交叉項）—— 無已註冊載體，即 S5")
print("   ⟹ G1 之「終結整個範式」**超出實際範圍**：終結嘅係『正核 × 模式指標差分』，非『所有核』。")

print("\n" + "=" * 96)
print("[裁決]")
print("   G1（消極主張）：🟡 **部分成立**。確實終結「正核／Green／截斷 × 模式差分」，")
print("                    但**未**終結帶號核（S3）與非線性乘積 ⟹ 「整個範式」係過度概括。")
print("   G2（積極主張）：🔴 **按可計算讀法全部判死**：")
print("                    R0 恆等零；R1 ρ_d=−1、ρ_u=+1（符號全反）；R2 與 x 無關（與 2.5× 變動矛盾）；")
print("                    R3 退化為 R10 形狀族（零證據）；R4 恆等零。")
print("   ⟹ G2 之唯一存活讀法 = **非張量分解之真交叉項** = 本 O 之 S5，而 S5 目前**無可算載體**。")
print("   ⟹ 依 SOUL：外部模型散文只作假設源；此處四個可算讀法皆已用真實計算判死。")
print("=" * 96)