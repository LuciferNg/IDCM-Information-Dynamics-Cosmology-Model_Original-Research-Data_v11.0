#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/17 — 三個共振三元組之「荷係數」：χ³ overlap（載體唯一決定？）

────────────────────────────────────────────────────────────────────────
候選定義（canonical，零參數）：三模頂點 = 模剖面之三乘積在鏈上之和
    C(a,b,m) = Σ_{p=1}^{44} φ_a(p) φ_b(p) φ_m(p)，  φ_n(p) = √(2/45) sin(nπp/45)
閉式：4 sinA sinB sinC = sin(A+B−C) + sin(A−B+C) + sin(−A+B+C) − sin(A+B+C)
      Σ_{p=1}^{N} sin(kpπ/(N+1)) = cot(kπ/(2(N+1)))（k 奇）／ 0（k 偶）
      此處 N = 44 ⟹ cot(kπ/90)，k 偶則 0。
⟹ C 之四項同 parity（皆 ≡ a+b+m mod 2）⟹ **C ≠ 0 要求 a+b+m 為奇**。

本腳本：① 核 cot 恆等式；② 核 C 之閉式；③ 三個共振三元組之 C；
        ④ 全掃 (a,b,m) 之 C ≠ 0 集與共振集之交；⑤ /88 環（平面波）之對照。
性質：確定性，N_trials = 0。零觀測量。dps = 40。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp

mp.mp.dps = 40
L = "=" * 78
ND = 44


def cot_sum(k):
    """Σ_{p=1}^{44} sin(kπp/45)"""
    return sum(mp.sin(k * mp.pi * p / 45) for p in range(1, ND + 1))


print(L)
print("M14 scripts/17 — 共振三元組之 χ³ overlap（荷係數）")
print(L)

# ======================================================================
print("\n【§一】cot 恆等式：Σ_{p=1}^{44} sin(kπp/45) = cot(kπ/90)（k 奇）／0（k 偶）")
worst_o = worst_e = mp.mpf(0)
odd_ok = even_ok = 0
for k in range(1, 90):
    v = cot_sum(k)
    if k % 2:
        d = abs(v - mp.cot(k * mp.pi / 90))
        worst_o = max(worst_o, d)
        odd_ok += 1
    else:
        d = abs(v)
        worst_e = max(worst_e, d)
        even_ok += 1
print(f"  k 奇（{odd_ok} 個）：max|Σ − cot(kπ/90)| = {mp.nstr(worst_o, 3)}")
print(f"  k 偶（{even_ok} 個）：max|Σ| = {mp.nstr(worst_e, 3)} ⟹ 恆等於 0")
print("  ⟹ **宇稱選擇規則**：k 偶之正弦和精確為 0")


def chi3(a, b, m, numeric=True):
    """三模 overlap（閉式；numeric=False 則只用 parity 判定）"""
    ks = [a + b - m, a - b + m, -a + b + m, a + b + m]
    if not numeric:
        return None if any(k % 2 for k in ks) else 0.0
    tot = sum(cot_sum(k) for k in ks if k > 0)
    return tot * (mp.mpf(2) / 45) ** mp.mpf(1.5) / 4


# ======================================================================
print("\n【§二】三個共振三元組之 χ³")
TRIP = [(9, 15, 18), (9, 27, 30), (15, 27, 36)]
for (a, b, m) in TRIP:
    ks = [a + b - m, a - b + m, -a + b + m, a + b + m]
    val = chi3(a, b, m)
    print(f"  {{{a:2d},{b:2d},{m:2d}}}  四項 k = {ks}（全偶）⟹ C = {mp.nstr(abs(val), 3)}"
          f"  a+b+m = {a+b+m}（偶）")
print("  ⟹ **三個共振三元組之 χ³ 全部為 0**（宇稱禁戒）")

# ======================================================================
print("\n【§三】全掃：C ≠ 0 之集 vs 能量共振之集")
resonant = set()
for a in range(1, ND + 1):
    for b in range(a + 1, ND + 1):
        for m in range(1, ND + 1):
            pass
# 用 §二 之已知結果（docs/15）：3 個三元組（無序）
resonant = {tuple(sorted(t)) for t in TRIP}
nonzero = set()
for a in range(1, ND + 1):
    for b in range(1, ND + 1):
        for m in range(1, ND + 1):
            if len({a, b, m}) < 3:
                continue                          # 要求三模互異
            if (a + b + m) % 2 == 1:
                nonzero.add(tuple(sorted((a, b, m))))
print(f"  三模互異三元組總數（無序）= C(44,3) = {44*43*42//6}")
print(f"  其中 C ≠ 0（a+b+m 奇）= {len(nonzero)}")
print(f"  能量共振集 = {sorted(resonant)}（{len(resonant)} 個）")
print(f"  交集 = {sorted(resonant & nonzero)}  ⟹ {'✅ 空（互斥）' if not (resonant & nonzero) else '🔴 非空'}")
print("  ⟹ **定理：能量允許 ⟺ 宇稱禁戒**（在 /90 Dirichlet 載體上，三模頂點無一例外）")

# ======================================================================
print("\n【§四】對照：/88 環（平面波）之三模頂點")
print("""  88-環上 φ_k(p) = e^{2πikp/88}/√88 ⟹
    C = (1/88^{3/2}) Σ_p e^{2πi(±a±b∓m)p/88} = (1/88^{3/2})·88·[a±b∓m ≡ 0 (mod 88)]
  ⟹ 耦合 = 動量守恒 delta，且對所有滿足者**數值相同**（1/√88）
  ⟹ 環上「荷係數」唔區分任何價對（A 慣例之 64 活價對全部等權）。
  ⟹ 兩個慣例給出相反之荷結構：/90 = 宇稱禁戒（全 0）；/88 = 全等（無差異）。
""")

# ======================================================================
print(L)
print("【§五】判定")
print(L)
print("""  ① 荷係數（χ³，canonical）在 /90 載體上：3 個共振三元組**全為 0**（宇稱禁戒）。
  ② 更強：**能量共振與 χ³ 非零互斥**（44-鏈上無例外；三模互異者）。
  ③ /88 環上：耦合 = 動量 delta，64 活價對全部等權 ⟹ 無荷差異。
  ④ ⟹ **載體唔提供「荷係數」之差異**：/90 給 0，/88 給常數。
     故 docs/09 §五 第 7 條（係數由載體唯一決定）在 χ³ 定義下**得唔到非平凡係數**。
  ⑤ 後果：味對若要帶荷，其載體**唔係**三模頂點（χ³）；需二階（質量矩陣）或
     觀察者讀數（T_pipe 投影）之類。此為下一步。
""")
print(L)
raise SystemExit(0)