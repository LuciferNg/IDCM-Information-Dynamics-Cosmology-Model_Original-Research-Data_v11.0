#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/09 — 價對層之分類：Galois 群之**對角作用**（canonical，零參數）

────────────────────────────────────────────────────────────────────────
性質：把 P1（Galois gcd 類，即 Gal ≅ (ℤ/88)^× 之作用）**對角**作用於價對空間
      {(a,b) : 1 ≤ a < b ≤ 44}（946 對），求其軌道分解。
理據：M12 docs/16 §九 規格第 1/3/7 條 ⟹ 對象係「階 2 價對」；P1 係唯一已註冊、
      作用於模集之**群作用**（其餘 P2–P8 係分割，非作用）。群作用才能誘導價對層之分類。
      本腳本**不引用任何觀測量**（無 ρ、無 δ、無質量表、無 14.86%）。
N_trials：0（軌道分解為確定性計算；無掃描、無擬合）

⚠ 誠實標記：我**已知** M12 docs/13 O3 謂「自由味對荷 = 13 型（6 qq ＋ 7 q q̄）」。
   故本輪若得出 13，其對照屬**後驗**，**不得**作證據。本輪只報軌道基數與結構，
   不作「命中」宣稱。
────────────────────────────────────────────────────────────────────────
"""
import math
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
L = "=" * 78

N = 88
ND = 44


def fold(r):
    """γ_n = γ_{88−n} ⟹ 把 ℤ/88 之代表元摺回 1..44。"""
    r %= N
    if r == 0:
        return 0
    return r if r <= ND else N - r


def gam(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print(L)
print("M14 scripts/09 — 價對層分類：Gal ≅ (ℤ/88)^× 之對角作用")
print(L)

# ======================================================================
print("\n【§一】群元素與作用之良定義（先核）")
units = [u for u in range(1, N) if math.gcd(u, N) == 1]
print(f"  |(ℤ/88)^×| = {len(units)}（φ(88) = {len(units)}）")
# ± 同一作用（因 fold）⟹ 有效群階
pairs_u = {min(u, N - u) for u in units}
print(f"  fold 後 u ~ 88−u ⟹ 有效群階 = {len(pairs_u)}（skill 記 20 ✅）")
act = None  # 不預先建表：作用須逐 u 計算（舊版 bug：dict 被覆蓋 ⟹ 946 個假單點軌道）
# 作用保持 γ？
maxdev = max(abs(gam(fold(n * u)) - gam(n)) for u in units for n in range(1, ND + 1))
print(f"  作用保持 γ ？max|γ_σ(n) − γ_n| = {mp.nstr(maxdev, 3)}  ⟹ {'✅' if maxdev < mp.mpf('1e-30') else '🔴'}")
# 作用在模集上之軌道數（應 = 7 = P1）
seen, orbits = set(), []
for n in range(1, ND + 1):
    if n in seen:
        continue
    o = sorted({fold(n * u) for u in units})
    seen |= set(o)
    orbits.append(o)
print(f"  模集軌道數 = {len(orbits)}（P1 之 gcd 類 = 7）⟹ {'✅ 一致' if len(orbits) == 7 else '🔴'}")
print(f"  軌道大小 = {sorted(len(o) for o in orbits)}")

# ======================================================================
print("\n【§二】對角作用於價對空間（946 對）")
pairs = [(a, b) for a in range(1, ND + 1) for b in range(a + 1, ND + 1)]
print(f"  價對總數 = {len(pairs)}")


def diag(p, u):
    a, b = fold(p[0] * u), fold(p[1] * u)
    return (a, b) if a < b else (b, a)


seen_p, porbits = set(), []
for p in pairs:
    if p in seen_p:
        continue
    o = {diag(p, u) for u in units}
    seen_p |= o
    porbits.append(sorted(o))
print(f"  對角軌道數 = {len(porbits)}")
sizes = sorted(len(o) for o in porbits)
from collections import Counter
print(f"  軌道大小分佈 = {dict(sorted(Counter(sizes).items()))}")
print(f"  最大軌道 = {max(sizes)}；最小 = {min(sizes)}")
print(f"  Σ 軌道大小 = {sum(sizes)}（應 = 946）⟹ {'✅' if sum(sizes) == 946 else '🔴'}")

# ======================================================================
print("\n【§三】軌道之結構清單（前 40 軌；全表見 log 尾部摘要）")
print(f"  {'#':>3s} {'代表 (a,b)':>12s} {'size':>5s} {'Δ':>3s} {'Σ':>3s}  {'gcd(Δ,88)':>9s} {'gcd(Σ,88)':>9s}")
rows = []
for i, o in enumerate(sorted(porbits, key=lambda o: (min(o), len(o))), 1):
    a, b = min(o)
    d, s = b - a, a + b
    rows.append((i, (a, b), len(o), d, s, math.gcd(d, N), math.gcd(s, N)))
    if i <= 40:
        print(f"  {i:>3d} {str((a,b)):>12s} {len(o):>5d} {d:>3d} {s:>3d}  {math.gcd(d,N):>9d} {math.gcd(s,N):>9d}")
if len(rows) > 40:
    print(f"  … 其餘 {len(rows)-40} 軌略（本 log 只列前 40；全部軌之基數已於 §二 統計）")

# ======================================================================
print("\n【§四】不變量檢定：軌道是否被 (gcd(Δ,88), gcd(Σ,88)) 唯一決定？")
key = {}
for i, p, sz, d, s, gd, gs in rows:
    key.setdefault((gd, gs), []).append(i)
multi = {k: v for k, v in key.items() if len(v) > 1}
print(f"  不同 (gcd(Δ),gcd(Σ)) 鍵數 = {len(key)}；對應多於一軌之鍵數 = {len(multi)}")
print(f"  最多軌之鍵 = {max((len(v) for v in key.values()), default=0)} 軌")
print("  ⟹ 若「多於一軌」為 0，則 (gcd(Δ),gcd(Σ)) 係價對層之**完備不變量**（分類閉合）")

# ======================================================================
print("\n【§五】對照：M12 之「13 型」（後驗，不作證據）")
print(f"""  本輪軌道數 = {len(porbits)}。M12 docs/13 O3 謂自由味對荷 13 型。
  ⚠ 此對照**後驗**（計算前我已知 13）⟹ 依禁區規則，不得引為命中或證據。
  可說者僅為：價對層之 canonical 群作用軌道數 = {len(porbits)}（確定性、零參數）。
  是否「等於味對型之數」需另立預註冊盲測（目標值由外部（PDG/群論）於跑前封存）。""")

# ======================================================================
print("\n【§六】beat/sum 之同軌檢定（定向可否由本作用給出？）")
same = 0
for i, p, sz, d, s, gd, gs in rows:
    if (d, s) in [(x[3], x[4]) for x in rows] and math.gcd(d, N) == math.gcd(s, N):
        same += 1
print(f"  軌道中 gcd(Δ) == gcd(Σ) 者 = {same} / {len(rows)}")
print("  ⟹ 若全部相等，則本作用**不區分** beat 與 sum ⟹ 與 P9 之宇稱定理同結論（定向需另一機制）")

print("\n" + L)
print(f"結論：價對層 canonical 分類軌道數 = {len(porbits)}；"
      f"完備不變量 = {'(gcd(Δ),gcd(Σ)) ✅' if not multi else '否（需更多不變量）'}")
print(L)
raise SystemExit(0)