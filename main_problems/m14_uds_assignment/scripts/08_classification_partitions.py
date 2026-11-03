#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/08 — 味分類：**不變量先驗清單（凍結）** ＋ 各分割之精確計算

────────────────────────────────────────────────────────────────────────
性質：分類不變量之**先驗凍結**與其誘導分割之精確計算。
      本腳本**不引用任何觀測量**（無 ρ、無 δ、無 14.86%、無質量表）。
      呢點係刻意的：分類必須在見到目標之前定死，否則即重犯「盲目地找」。

凍結聲明（跑前寫定，不得事後增刪）：
  清單 P1..P10 之每一項，其出處皆為 corpus 已註冊之結構陳述（逐項標明）。
  任何**未列於此**之分割，不得用於本輪之分類判定。

慣例（依 M12 docs/40 §六 約束 3，不得混用）：
  44-site 鏈：γ_n = 4sin²(nπ/88)，n = 1..44，Σ = 90
  88-site 雙鏈：γ'_k = 4sin²(kπ/176)，k = 1..88，Σ = 178 = 2F₁₁
  ⚠ 1-indexed parity（0-indexed slicing 會把 even/odd 互換 — 已知陷阱）

N_trials：0（無掃描、無擬合；本腳本只做分割與恆等式核對）
────────────────────────────────────────────────────────────────────────
"""
import math
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
L = "=" * 78


def sec(x):
    return 1 / mp.cos(x)


def gam44(n):                      # 44-site K_chain
    return 4 * mp.sin(mp.pi * n / 88) ** 2


def gam88(k):                      # 88-site dual chain
    return 4 * mp.sin(mp.pi * k / 176) ** 2


def nstr(x, d=12):
    return mp.nstr(x, d)


print(L)
print("M14 scripts/08 — 味分類：不變量先驗清單（凍結）＋ 分割精確計算")
print(L)

# ======================================================================
print("\n【§〇】凍結清單（跑前寫定）")
print("""
  P1  Galois 軌道（gcd(n,88) 等價類）           出處: skill idcm-88site-dual-chain
  P2  88-chain mod-4 扇區（4 類，精確閉式）      出處: O148 S5 (43/43 PASS)
  P3  88-chain 奇偶（2 類：90 / 88 精確）        出處: O218 / skill
  P4  44-chain 奇偶（2 類：46 / 44 精確，1-idx） 出處: O210 C2 裁決
  P5  T_pipe 支撐（n≤27 投影 / n≥28 核，2 類）   出處: O121 (n*=27, F_proj=17)
  P6  block(n)=⌈rank/3⌉（15 類）                 出處: O137 docs/35 §一 ③
  P7  k mod 11（11 類，88=8×11）                 出處: 88 之因數分解
  P8  k mod 8（8 類，88=8×11）                   出處: 同上
  P9  價對層 (a,b) → (Δ,Σ) 之宇稱定理            出處: 本輪（代數，見 §九）
  P10 計數恆等式（3 之來源枚舉，非分割）         出處: 逐項標明
""")

# ======================================================================
print("【§一】地基恆等式（先核，後方可用）")
s44 = sum(gam44(n) for n in range(1, 45))
s88e = sum(gam88(2 * m) for m in range(1, 45))
s88o = sum(gam88(2 * m - 1) for m in range(1, 45))
print(f"  Σγ^(44)  n=1..44        = {nstr(s44, 20)}   dev(90) = {nstr(abs(s44-90), 3)}")
print(f"  88-chain even k=2m      = {nstr(s88e, 20)}   dev(90) = {nstr(abs(s88e-90), 3)}")
print(f"  88-chain odd  k=2m-1    = {nstr(s88o, 20)}   dev(88) = {nstr(abs(s88o-88), 3)}")
print(f"  Σγ^(88) total           = {nstr(s88e + s88o, 20)}   (2F₁₁ = 178)")
carry_max = max(abs(gam88(2 * n) - gam44(n)) for n in range(1, 45))
print(f"  O218 carry γ'_(2n)=γ_n  : max|Δ| = {nstr(carry_max, 3)}")
pair_max = max(abs((gam44(n) + gam44(44 - n)) - 4) for n in range(1, 45))
print(f"  γ_n + γ_(44−n) = 4      : max|Δ| = {nstr(pair_max, 3)}")

# ======================================================================
print("\n【§二】P1 — Galois 軌道（gcd(n,88)）")
orb = {}
for n in range(1, 45):
    orb.setdefault(math.gcd(n, 88), []).append(n)
sizes = {}
print(f"  {'gcd':>4s} {'size':>5s}  成員")
for g in sorted(orb, key=lambda g: (len(orb[g]), -g)):
    sizes[g] = len(orb[g])
    mem = orb[g]
    show = "{" + ",".join(map(str, mem[:6])) + ("…" if len(mem) > 6 else "") + "}"
    print(f"  {g:>4d} {len(mem):>5d}  {show}")
print(f"  ⟹ 類數 = {len(orb)}；size 分佈 = {sorted(sizes.values())}；總 = {sum(sizes.values())}")
print("  ⟹ 對照 skill 記載 [1,1,2,5,5,10,20]：", "✅ 吻合" if sorted(sizes.values()) == [1, 1, 2, 5, 5, 10, 20] else "🔴 不符")

# ======================================================================
print("\n【§三】P2 — 88-chain mod-4 扇區（精確閉式核對）")
sec88 = sec(mp.pi / 88)
cls = {0: [], 1: [], 2: [], 3: []}
for k in range(1, 89):
    cls[k % 4].append(k)
want = {0: mp.mpf(46), 1: 44 - sec88, 2: mp.mpf(44), 3: 44 + sec88}
print(f"  sec(π/88) = {nstr(sec88, 15)}")
for r in (0, 1, 2, 3):
    got = sum(gam88(k) for k in cls[r])
    print(f"  s{r} (k≡{r}, {len(cls[r])} 模)  Σ = {nstr(got, 14):>18s}   期望 {nstr(want[r], 14):>18s}"
          f"   dev = {nstr(abs(got - want[r]), 3)}")
print(f"  ⟹ 世界側 s0+s2 = {nstr(want[0] + want[2], 12)}（精確整數）；"
      f"觀察者側 s1+s3 = {nstr(want[1] + want[3], 12)}（±sec 對消）")

# ======================================================================
print("\n【§四】P3／P4 — 奇偶分割")
se = sum(gam44(n) for n in range(2, 45, 2))
so = sum(gam44(n) for n in range(1, 45, 2))
print(f"  P4 44-chain（1-indexed）：even n = {nstr(se, 12)}，odd n = {nstr(so, 12)}（22/22 模）")
print(f"  P3 88-chain：even k = {nstr(s88e, 12)}，odd k = {nstr(s88o, 12)}（44/44 模）")
print("  ⚠ 兩者係唔同對象（46/44 vs 90/88）；0-index slicing 會互換 — 本輪用 range(2,45,2) 顯式 1-index")

# ======================================================================
print("\n【§五】P5 — T_pipe 支撐")
nd, l6 = 44, 18
nstar = nd - l6 + 1
s_le = sum(gam44(n) for n in range(1, nstar + 1))
s_gt = sum(gam44(n) for n in range(nstar + 1, nd + 1))
fproj = (mp.mpf(57) / 2 - sum(gam44(n) for n in range(1, nstar))) / gam44(nstar)
print(f"  n* = N_d − L₆ + 1 = {nstar}；|ker T_pipe| = N_d − n* = {nd - nstar} = F_proj")
print(f"  Σγ_(1..27) = {nstr(s_le, 12)}；Σγ_(28..44) = {nstr(s_gt, 12)}")
print(f"  f_proj = (57/2 − Σγ_(1..26))/γ_27 = {nstr(fproj, 12)}（skill 記 0.771960）")
print(f"  ⟹ 類大小 = {nstar} / {nd - nstar}（{nstar/nstar:.0f} 與 {1}）")

# ======================================================================
print("\n【§六】P6 — block(n) = ⌈rank/3⌉（O137 docs/35 §一 ③）")
nblk = (44 + 2) // 3
blk_sizes = {}
for r in range(1, 45):
    b = (r + 2) // 3
    blk_sizes[b] = blk_sizes.get(b, 0) + 1
print(f"  rank = 1..44 ⟹ 類數 = ⌈44/3⌉ = {nblk}；類大小 = {sorted(blk_sizes.values(), reverse=True)[:3]}…"
      f"（3 個一組，末組 {blk_sizes[nblk]}）")
print("  ⚠ 此分割依 **rank 次序**（排序後位置），非站/模之座標 ⟹ 對「邊個味」無幾何內容")

# ======================================================================
print("\n【§七】P7／P8 — 因數分割（88 = 8 × 11）")
for d in (11, 8):
    c = {}
    for k in range(1, 89):
        c.setdefault(k % d, 0)
        c[k % d] += 1
    print(f"  k mod {d}: {len(c)} 類，類大小 {sorted(set(c.values()))}")
print("  ⟹ 88 之因數 = {1,2,4,8,11,22,44,88} ⟹ 週期性分割可得之類數只有呢批（3 ∤ 88）")

# ======================================================================
print("\n【§九】P9 — 價對層宇稱定理（代數）")
n_pair = 0
bad = 0
for a in range(1, 45):
    for b in range(a + 1, 45):
        n_pair += 1
        if ((a + b) - (a - b)) % 2 != 0:
            bad += 1
print(f"  (a,b) 對數 = C(44,2) = {n_pair}")
print(f"  (Σ − Δ) = 2b ⟹ 恆為偶數 ⟹ parity(Σ) ≡ parity(Δ) ；違例 {bad}/{n_pair}")
print("  ⟹ **定理**：beat（Δ）與 sum（Σ）必同宇稱 ⟹ 宇稱**不可能**分辨物質／反物質扇區")
print("     — 與 M12 docs/16 §九 第 6 條「定向不可為單一 ℤ₂」同向（獨立推得）")
inj = len({((a + b), (b - a)) for a in range(1, 45) for b in range(a + 1, 45)})
print(f"  (a,b) ↦ (Δ=b−a, Σ=a+b) 單射？{'是' if inj == n_pair else '否'}（image = {inj}）")
print("     ⟹ 價對空間 = {(Δ,Σ) : Δ<Σ, Δ≡Σ (mod 2), 1≤Δ≤43, 3≤Σ≤87}")

# ======================================================================
print("\n【§十】P10 — 「3」之來源枚舉（計數恆等式，非分割）")
items = [
    ("B − 1", F(3), "B=4 分枝數（已註冊）"),
    ("N_gen", F(3), "世代數（輸入，見 M14 docs/06 Q4）"),
    ("√(9) = √(N_gen²)", F(3), "O137: 9 模式 / I-60 軌道 = 9"),
    ("3³", F(27), "= n* = N_d − L₆ + 1（T_pipe 截斷點）"),
    ("L₆ = φ⁶+φ⁻⁶", F(18), "= 18（O121 F_proj = L₆−1）"),
    ("⌈N_d/3⌉", F(15), "= 15（O137 block 數）"),
]
for name, val, src in items:
    print(f"  {name:>18s} = {str(val):>4s}   {src}")
print("  ⚠ 本節**只列**，不作任何「因此 3 係咩」之宣稱（避免凑數）")

# ======================================================================
print("\n" + L)
print("【§十一】判定：有無「canonical 且 3 類」之分割？")
print(L)
tbl = [
    ("P1 Galois gcd", 7, "不屬於 3"),
    ("P2 mod-4", 4, "不屬於 3"),
    ("P3 88 奇偶", 2, "不屬於 3"),
    ("P4 44 奇偶", 2, "不屬於 3"),
    ("P5 T_pipe 支撐", 2, "不屬於 3"),
    ("P6 block(rank/3)", 15, "以 rank 分組，非幾何"),
    ("P7 k mod 11", 11, "不屬於 3"),
    ("P8 k mod 8", 8, "不屬於 3"),
]
for name, c, note in tbl:
    print(f"  {name:<20s} 類數 = {c:>2d}   {note}   {'★ 3 類' if c == 3 else ''}")
print(f"""
  ⟹ 凍結清單 P1–P8 之中，類數為 3 者：**0 個**（N_trials = 8，全報）。
  ⟹ **定理（結構，非枚舉）**：3 ∤ 88 ⟹ 任何「站/模座標之週期性」分割，其類數必係 88 之因數
     （或該因數之倍數），而 88 = 2³·11 ⟹ {3} ∉ divisors(88)。
     故「3 個味」**不可能**由座標週期性直接讀出。同源 no-go：O216 之 Δn_color（3 ∤ 88）。
  ⟹ 剩餘可行來源收窄為三類（皆需前置工作，非掃描）：
     ① **非阿貝爾群之 3 維不可約表示**（A₅ 有兩個 3-dim irrep）—— 但 O137 docs/35 §一
        之歸約定理謂「A₅ 幾何對每模式供給之數據完全由 n 決定」⟹ 此路受該定理約束，
        須先用 O137 之軌道定義（achiral 24 / chiral 9）核實 A₅ 作用於 **88 站**（非模式）。
     ② **計數恆等式**（§十）—— 只能給「3」之出現位置，不能給「邊三類」（無分割內容）。
     ③ **價對層組合**（P9：946 對，單射至 (Δ,Σ) 格）—— 分類對象在此層，但格數 = 946，
        需一個 group action 把它摺疊；候選 = P1／P2 在 Δ、Σ 上之同時作用。
""")
print(L)
print("結論：味分類**未完成**。已排除「座標週期性 ⟹ 3 類」。本輪零掃描、零擬合。")
print(L)
raise SystemExit(0)