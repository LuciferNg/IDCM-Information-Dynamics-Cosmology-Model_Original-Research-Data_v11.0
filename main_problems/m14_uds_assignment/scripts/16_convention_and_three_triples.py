#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/16 — ★ 慣例決定結果：/90 Dirichlet 鏈恰有 **3** 個三波共振三元組

────────────────────────────────────────────────────────────────────────
背景：scripts/14-15 之共振定律係在 **/88 週期** 慣例（γ_n = 4sin²(nπ/88)）下導出。
      但 M14 之站空間（φ_n(j) = √(2/45) sin(πnj/45)）之特徵值係
      **/90 Dirichlet**：γ_n = 4sin²(nπ/90)，n = 1..44（M14 docs/01 §3 已分立兩慣例）。
      ⚠ 兩者係唔同對象，**唔可以混用**。本腳本把兩者分開掃盡。

三條鏈（全掃 a<b、m ∈ 1..M，dps = 60，精確閾 1e−45）：
  A. /88 週期，44 模式  γ_n = 4sin²(nπ/88)    （2N = 88 ⟹ N = 44）
  B. /88 雙鏈，88 模式  γ_k = 4sin²(kπ/176)
  C. /90 Dirichlet，44 站 γ_n = 4sin²(nπ/90)  ← M14 站空間所用者

性質：確定性掃盡，N_trials = 0。零觀測量。精確部分用 ℚ(√5) 之 Fraction 算術。
────────────────────────────────────────────────────────────────────────
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 60
TOL = mp.mpf('1e-45')
L = "=" * 78


class Q5:
    """a + b√5，a,b ∈ ℚ —— 精確算術。"""
    def __init__(self, a, b=0):
        self.a, self.b = F(a), F(b)

    def __add__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    def __sub__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return Q5(self.a - o.a, self.b - o.b)

    def __eq__(self, o):
        o = o if isinstance(o, Q5) else Q5(o)
        return self.a == o.a and self.b == o.b

    def __repr__(self):
        return f"({self.a}{'+' if self.b >= 0 else '−'}{abs(self.b)}√5)"


def scan(den, M, label):
    g = [mp.mpf(0)] * (M + 1)
    for n in range(1, M + 1):
        g[n] = 4 * mp.sin(mp.pi * n / den) ** 2
    exS, exD = [], []
    for a in range(1, M + 1):
        for b in range(a + 1, M + 1):
            v, d = g[a] + g[b], g[b] - g[a]
            for m in range(1, M + 1):
                if abs(v - g[m]) < TOL:
                    exS.append((a, b, m))
                if abs(d - g[m]) < TOL:
                    exD.append((a, b, m))
    print(f"\n  {label}")
    print(f"    價對 = {M*(M-1)//2}；三格 = {M*M*(M-1)//2}")
    print(f"    sum  精確解 = {len(exS):>3d} {exS if len(exS) <= 8 else exS[:6] + ['…']}")
    print(f"    beat 精確解 = {len(exD):>3d} {exD if len(exD) <= 8 else exD[:6] + ['…']}")
    return exS, exD


print(L)
print("M14 scripts/16 — 慣例隔離：/88 週期 vs /90 Dirichlet")
print(L)

print("\n【§一】三條鏈全掃")
sA, dA = scan(88, 44, "A. /88 週期，44 模式（= 88-鏈偶扇區）")
sB, dB = scan(176, 88, "B. /88 雙鏈，88 模式")
sC, dC = scan(90, 44, "C. /90 Dirichlet，44 站  ← M14 站空間")
print("\n  ⟹ **慣例決定結果**：A/B 給大量解；C 只給 3（sum）＋ 6（beat）")

# ======================================================================
print("\n" + L)
print("【§二】C（/90）之三元組：完整枚舉與精確 ℚ(√5) 驗證")
print(L)
trip = {}
for (a, b, m) in sC:
    trip.setdefault(tuple(sorted((a, b, m))), []).append(("sum", (a, b, m)))
for (a, b, m) in dC:
    trip.setdefault(tuple(sorted((a, b, m))), []).append(("beat", (a, b, m)))
print(f"  三元組數 = {len(trip)}")
for t, how in sorted(trip.items()):
    print(f"    {t}  ↔  {how}")

print("\n  精確值（ℚ(√5)/8 表示）：")
S2 = {18: Q5(3, -1), 30: Q5(2, 0), 36: Q5(5, -1),
      54: Q5(3, 1), 60: Q5(6, 0), 72: Q5(5, 1)}
for k, v in S2.items():
    print(f"    sin²({k}°) = (1/8){v}")
print("\n  逐三元組驗證（γ = 4sin²，故等價於 sin² 相加）：")
for t in sorted(trip):
    i, j, k = t
    Ai, Aj, Ak = (2 * x for x in t)          # γ_n = 4sin²(nπ/90) ⟹ 角度 = 2n 度
    lhs = S2[Ai] + S2[Aj]
    rhs = S2[Ak]
    ok = "✅" if lhs == rhs else "🔴"
    print(f"    {{{i},{j},{k}}} ⟹ sin²{Ai}° + sin²{Aj}° = {lhs}/8 ＝ sin²{Ak}° = {rhs}/8  {ok}")

print("""
  三條恆等式（皆 ℚ(√5) 精確）：
    (i)   sin²18° + sin²30° = sin²36°   ⟺ cos36° − cos72° = cos60° = ½
    (ii)  sin²18° + sin²54° = sin²60°   ⟺ cos36° + cos108° = ½
    (iii) sin²30° + sin²54° = sin²72°   ⟺ cos60° + cos108° = 1 + cos144°
""")

# ======================================================================
print("【§三】為何只有 3 個？—— 角度落於 ℚ(√5 ∪ √3) 之交集")
print("""
  index n ⟹ 角度 2n°（因 γ_n = 4sin²(2n·π/180)）。參與者之角度：
    9→18°, 15→30°, 18→36°, 27→54°, 30→60°, 36→72°
  ⟹ 全部落於 {18°,30°,36°,54°,60°,72°} = 18° 之倍數 ∪ {30°}
     = ℚ(ζ₅)⁺（√5）與 ℚ(√3) 之交角集；而 sin² 在其中取值只有 6 個
     {(3−√5)/8, 1/4, (5−√5)/8, (3+√5)/8, 3/4, (5+√5)/8}
  ⟹ 三波共振 = 呢 6 個值之間之加法關係，窮盡得 3 條（上列）。
  ⟹ 非擬合、非選取：由 ℚ(√5) 之加法表**逼出**，且掃盡確認無第四條。
""")

# ======================================================================
print(L)
print("【§四】判定")
print(L)
print(f"""  ① 慣例隔離為硬閘：/88 週期 ⟹ {len(sA)}+{len(dA)}；/88 雙鏈 ⟹ {len(sB)}+{len(dB)}；
     **/90 Dirichlet ⟹ {len(sC)}+{len(dC)}，恰 3 個三元組**。
  ② M14 站空間（φ_n(j)=√(2/45)sin(πnj/45)，特徵值 4sin²(nπ/90)）用嘅係 **/90**
     ⟹ 對本 M 之味對問題，正確載體係 C，**唔係** A/B。
  ③ 故 scripts/14-15 之 64 活價對（/88）**唔適用於站空間**；本輪結果取代其角色。
  ④ ★ 首次出現 canonical 3：{{9,15,18}}、{{9,27,30}}、{{15,27,36}}，
     全部為 ℚ(√5) 精確恆等式，掃盡確認無第四條。
  ⑤ ⚠ 未宣稱：「3 個三元組 = 3 個世代／3 個味」。此為**結構事實**，
     映射到 u,d,s 須另立規則（見 docs/15 §五）。
""")
print(L)
raise SystemExit(0)