#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/19 — 四波（W⁴）共振結構 ＋ mod-3 子晶格約化

────────────────────────────────────────────────────────────────────────
動機：IDCM 之自耦合係 V(W) = ½φW² + (ε/2)W⁴（已註冊）⟹ **物理頂點係四波，唔係三波**。
      而 scripts/17 已證三波 χ³ 對三個共振三元組全為 0（宇稱禁戒）。
      ⟹ 本輪檢四波：能量守恒 γ_a+γ_b = γ_c+γ_d，以及 χ⁴ 選擇規則。

性質：確定性掃盡；dps = 80，完全相同判定閾 1e−50。N_trials = 0。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp
from collections import Counter

mp.mp.dps = 80
L = "=" * 78
ND = 44
g = [mp.mpf(0)] + [4 * mp.sin(n * mp.pi / 90) ** 2 for n in range(1, ND + 1)]
pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]
TRIP = [(9, 15, 18), (9, 27, 30), (15, 27, 36)]
TM = {m for t in TRIP for m in t}

print(L)
print("M14 scripts/19 — 四波共振（pair-sum 譜）＋ mod-3 子晶格")
print(L)

# ======================================================================
print("\n【§一】pair-sum 譜：γ_a + γ_b 之等值類")
vals = sorted([(g[a] + g[b], (a, b)) for (a, b) in pairs], key=lambda t: t[0])
cl, cur = [], [vals[0]]
for v, p in vals[1:]:
    if abs(v - cur[-1][0]) < mp.mpf('1e-50'):
        cur.append((v, p))
    else:
        cl.append(cur)
        cur = [(v, p)]
cl.append(cur)
hist = Counter(len(c) for c in cl)
print(f"  價對 {len(pairs)} ⟹ 等值類 {len(cl)}；類大小分佈 = {dict(sorted(hist.items()))}")
print(f"  四波共振（兩對同和、互不相等）= 每類內任取兩對 ⟹ Σ_k C(size,2)"
      f" = {sum(len(c)*(len(c)-1)//2 for c in cl)}")

# ======================================================================
print("\n【§二】22 類：γ_a + γ_b = 4 ⟺ a + b = 45（精確）")
big = [c for c in cl if len(c) == 22][0]
ps = sorted(p for _, p in big)
print(f"  該類 22 對：{ps[:10]} … {ps[-3:]}")
ok45 = all(a + b == 45 for (a, b) in ps)
sum4 = all(abs(g[a] + g[b] - 4) < mp.mpf('1e-60') for (a, b) in ps)
print(f"  全部 a+b = 45 ？{ok45}；全部 γ_a+γ_b = 4 ？{sum4}")
print("  證明：γ_n = 2−2cos(4n°) ⟹ γ_a+γ_b = 4 ⟺ cos4a° + cos4b° = 0")
print("        ⟺ 4(a+b) ≡ 180 (mod 360) ⟺ a+b = 45（因 a+b ≤ 87，唯一解）")
print("        ⟹ 22 對（a=1..22, b=45−a）；**對所有 a 成立** ⟹ 係**對稱恆等式**，非場巧合")

# ======================================================================
print("\n【§三】12 個雙元類：全部落於 mod-3 子晶格？")
dbl = [c for c in cl if len(c) == 2]
allmod3 = True
for c in dbl:
    p1, p2 = c[0][1], c[1][1]
    if (p1[0] % 3 or p1[1] % 3 or p2[0] % 3 or p2[1] % 3):
        allmod3 = False
print(f"  雙元類 = {len(dbl)}；全部四模皆 ≡ 0 (mod 3) ？{allmod3}")
print(f"  3ℤ ∩ [1,{ND}] 之模數 = {len([n for n in range(1, ND+1) if n % 3 == 0])}"
      f"（= {sorted(n for n in range(1, ND+1) if n % 3 == 0)}）")
sub = [p for p in pairs if p[0] % 3 == 0 and p[1] % 3 == 0]
print(f"  子晶格價對 = {len(sub)} = C(14,2)；其中構成雙元類者 = {len(dbl)*2}")
print("  （22 類之 a+b=45 家族**不**受 mod-3 限制 —— 因佢係對稱恆等式）")

# ======================================================================
print("\n【§四】子晶格約化：γ_(3k) = 4sin²(kπ/30) ⟹ 15-模 /30 鏈")
mx = max(abs(g[3 * k] - 4 * mp.sin(k * mp.pi / 30) ** 2) for k in range(1, 15))
print(f"  max|γ_(3k) − 4sin²(kπ/30)| = {mp.nstr(mx, 3)}（k=1..14）")
print("  ⟹ 共振活性扇區 = k = 1..14 之 /30 鏈（半週期 30 = N_gen × 10）")
print(f"  該鏈之互補恆等式：γ_k + γ_(15−k) = 4 ⟹ 在 k=1..14 內可實現（k=1..7 配對）")

# ======================================================================
print("\n【§五】場論理由：為何只有 mod-3 會有非平凡重合")
print("""  γ_n = 2 − 2cos(4n°)。γ_a+γ_b 之重合 = 4 個 cos 值間之加減關係。
  cos(4n°) 之最小場：4n° ≡ 0 (mod 12°) ⟺ n ≡ 0 (mod 3) ⟹ cos ∈ ℚ(ζ₃₀)⁺ = ℚ(√3,√5)（次數 4）
  反之 n ≢ 0 (mod 3) 者，4n° 不在 30 分圓內 ⟹ 其 cos 與 mod-3 模之 cos 無低次關係。
  ⟹ **非平凡重合只可能發生於 n ≡ 0 (mod 3)**（本輪 12/12 全符）。
  三個三波三元組亦全部在 3ℤ，且其角度 {18°,36°,54°,72°} ⊂ ℚ(√5)∪ℚ，故進一步收窄至 √5 子域。
  ⚠ 此為結構解釋，非嚴格證明（嚴格證明需 ℚ(√3,√5) 之加性關係窮舉）。
""")

# ======================================================================
print("【§六】χ³ 與 χ⁴ 選擇規則（互補）")
print("""  χ³：C = Σ_p φ_aφ_bφ_m，四項 k = a+b−m, a−b+m, −a+b+m, a+b+m（同 parity）
      S(k) = cot(kπ/90)（k 奇）/ 0（k 偶）⟹ **C ≠ 0 ⟺ a+b+m 奇**
  χ⁴：D = Σ_p φ_aφ_bφ_cφ_d，八項 cos 之 k 亦同 parity
      Σ_{p=1}^{44} cos(kπp/45) = 0（k 奇）/ −1（k 偶）⟹ **D ≠ 0 ⟺ a+b+c+d 偶**
  ⟹ 兩規則**互補**：三波所禁者，四波所允。
""")
print("  三個共振三元組（a+b+m 皆偶）⟹ χ³ = 0（禁）、χ⁴ ≠ 0（允）")
for t in TRIP:
    print(f"    {t}: a+b+m = {sum(t)}（偶）⟹ χ³ 禁 / χ⁴ 允")

# ======================================================================
print("\n" + L)
print("【§七】判定")
print(L)
print(f"""  ① pair-sum 譜：946 對 ⟹ 913 類（900 單、12 雙、1 個 22 類）。
  ② 22 類 = **對稱恆等式** a+b=45 ⟹ γ_a+γ_b = 4（22 對，全 a）。
  ③ 12 個雙元類 = **非平凡場重合**，全部落於 mod-3 子晶格 n=3k（12/12）。
  ④ 子晶格約化：γ_(3k) = 4sin²(kπ/30) ⟹ 共振活性扇區 = 14 模 /30 鏈。
  ⑤ 選擇規則互補：χ³ 需奇和、χ⁴ 需偶和 ⟹ 三個三波三元組「χ³ 暗、χ⁴ 亮」。
  ⑥ 物理解讀（謹慎）：IDCM 之頂點係 W⁴ ⟹ 三波暗正係**自然**（三波本非物理頂點）；
     四波亮 ⟹ 三個三元組可經四波過程參與。
  ⑦ 未決：「14 模 /30 鏈」是否即「味活性扇區」；仍無 3 之直接來源
     （14 = N_d−30? 非；14 = L₆−4? 非 —— 不强行解釋）。
""")
print(L)
raise SystemExit(0)