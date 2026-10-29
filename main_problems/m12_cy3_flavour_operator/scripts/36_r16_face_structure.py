#!/usr/bin/env python3
"""36 — R16：面型（2-cell）路線之結構檢定

未做之推：docs/27 §六 明示「面（2-cell）尚未真正測試」。本輪做三項結構檢定。

預註冊
--------------------------------------------------
Q1 承載體有無面？𝒯 = P₄₄ ⊗ P₆₈（O210 charter C9：K_tree 係**路徑 P₆₈**之計數 ⟹ 兩因子皆路徑）。
   Cartesian 積之方格數 = (44−1)(68−1)；代數 2-形式層 dim Λ²(I₄₄) = C(44,2)。
Q2 通量是否均勻？已註冊 link phase（charter C6）：鏈 2θ_k = kπ/88、樹 2θ'_ℓ = ℓπ/68。
   plaquette flux F(a,b) = [2θ_{a+1} − 2θ_a] + [2θ'_{b+1} − 2θ'_b]。
   若 F 與 (a,b) 無關 ⟹ **均勻曲率 ⟹ 面無結構 ⟹ 味盲**（範疇性，非數值不中）。
Q3 面荷之前提：阻挫（frustration）需要**非二部**（奇環）。鏈係二部 ⟹ 阻挫邊 = 0。
   對照 docs/10（R3-P2）：RH 60v/120e 之阻挫邊 = 32 ⟹ 面之所以有內容，係因 J_ij 有 ±1 非均勻。

判準（先寫死）
--------------------------------------------------
F1  F(a,b) 之 max−min = 0（精確）⟹ 通量均勻 ⟹ 面型在鏈載體上判死（結構性）
F2  鏈阻挫邊 = 0（二部）⟹ 無阻挫 ⟹ 無面荷來源
F3  若 F 非均勻 ⟹ 面型仍活，須另立協議
"""
import mpmath as mp

mp.mp.dps = 40
NC, NT = 44, 68          # 鏈站數、樹站數（K_tree^(node) = 2(B²+1) = 68）

print("=" * 96)
print("36 — R16：面型（2-cell）結構檢定")
print("=" * 96)

print("\n[Q1] 承載體之面計數")
sq = (NC - 1) * (NT - 1)
lam2 = NC * (NC - 1) // 2
print("   𝒯 = P₄₄ ⊗ P₆₈ ⟹ Cartesian 方格數 = (44−1)×(68−1) = %d" % sq)
print("   代數 2-形式層 dim Λ²(I₄₄) = C(44,2) = %d（O210 Ix: (a,b) ∈ Λ²(I₄₄)）" % lam2)
print("   ⟹ 面**存在**（撤回 docs/27 §六「未註冊」之說：面存在，問題係下面 Q2）")

print("\n[Q2] plaquette 通量之均勻性（charter C6 之已註冊 link phase）")
th_c = lambda k: mp.mpf(k) * mp.pi / 88          # 鏈：2θ_k = kπ/88
th_t = lambda l: mp.mpf(l) * mp.pi / 68          # 樹：同構約定 2θ'_l = lπ/68
F = lambda a, b: (th_c(a + 1) - th_c(a)) + (th_t(b + 1) - th_t(b))
vals = [F(a, b) for a in range(1, NC) for b in range(1, NT)]
mn, mx = min(vals), max(vals)
print("   F(a,b) 之值域：min = %s  max = %s  差 = %s"
      % (mp.nstr(mn, 12), mp.nstr(mx, 12), mp.nstr(mx - mn, 12)))
print("   解析值 F = π/88 + π/68 = %s" % mp.nstr(mp.pi / 88 + mp.pi / 68, 12))
print("   ⟹ **F 與 (a,b) 無關（均勻曲率）** ⟹ 任何面函數皆為常數 ⟹ 面內差 δ ≡ 0")
print("   ⟹ F1 成立：面型在**已註冊 link phase** 下判死（結構性，非數值不中）")

print("\n[Q3] 面荷之前提：阻挫 vs 二部性")
chain_edges = NC - 1
print("   鏈 P₄₄：邊 = %d，頂點 = %d，環秩 = E − V + 1 = %d ⟹ 無環 ⟹ 二部 ⟹ **阻挫邊 = 0**"
      % (chain_edges, NC, chain_edges - NC + 1))
print("   對照 docs/10（R3-P2，RH 60v/120e/4-regular）：阻挫邊 = **32** ⟹ J_ij 有 ±1 非均勻")
print("   ⟹ 面之所以在 CMB 層有內容，係因 RH 有奇環（三角/五邊）＋ 阻挫；鏈載體兩者皆無")
print("   ⟹ F2 成立：鏈載體上『面荷』無來源（非「未做」，係前提不成立）")

print("\n" + "=" * 96)
print("[判定]")
print("   面型（2-cell）在 𝒯 = P₄₄ ⊗ P₆₈ 上：面存在（%d 方格、%d 個 Λ² 指標）" % (sq, lam2))
print("   但 ①通量均勻（F ≡ π/88+π/68）②鏈二部 ⟹ 阻挫 = 0 ⟹ **面無結構**")
print("   ⟹ 味分裂不可能由面荷生成（結構性判死，強於數值不中）")
print("   ⟹ 若要面型可行，必須先註冊一個**非二部、帶 ±1 阻挫耦合**之因子 —— 即 S5 之載體缺口，")
print("      且該因子不得係 P₄₄ 或 P₆₈（兩者皆二部）。")
print("=" * 96)