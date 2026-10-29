#!/usr/bin/env python3
"""R12 — S2 三攻：88-site 雙鏈載體之有限差分（先驗形式；`N_trials = 2`）

載體（**已註冊定義**，見 skill `idcm-88site-dual-chain`）
--------------------------------------------------
88-site 雙鏈：γ'_k = 4sin²(πk/176)，k = 1..88；Σγ'_even = 90（世界側，k=2n
重構 44-site 譜 γ_n = 4sin²(πn/88)）、Σγ'_odd = 88（觀察者側）。
Even-mode carry（O218，G3 裁決之唯一遷移）：γ'_{2n} = γ_n^{(44)}。

★ 為何張量載體 `K_chain ⊗ K_tree` **不**構成獨立測試（代數，非假設）
  K 可分解 ⟹ G_tensor(i,j) = Σ_{n,ℓ} φ_n(i)φ_n(j)/(γ_n μ_ℓ)
                            = (Σ_ℓ 1/μ_ℓ) · Σ_n φ_n(i)φ_n(j)/γ_n = const × G_chain(i,j)
  樹因子係**共同常數**，在比值中精確相消 ⟹ 與 R11 之 P1 逐位相同 ⟹ 退化，`N_trials` 不增。
  故本輪只剩兩條**真正獨立**之協議。

預註冊協議（跑前凍結）
--------------------------------------------------
A1  carry-by-2：把已註冊指數（u→2, d→3, s→4）按 O218 偶模 carry 遷移至 88-site
    ⟹ u→4, d→6, s→8；δ_x = C[G₈₈(6, k_x) − G₈₈(4, k_x)]。
A2  direct：直接把已註冊指數當 88-site 指數（u→2, d→3, s→4）；
    δ_x = C[G₈₈(3, k_x) − G₈₈(2, k_x)]。
判準 K1：主協議 A1 之 |z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中。
      K2：A2 全報（不計入 K1；若 A1 不中而 A2 中，須標明 `N_trials = 2` 之選擇效應）。
      K3：報 C 與由 δ_b、δ_c 反解之 k（後驗，不作證據）。
負控制：同 harness 跑 R11 之已死形式 `G₄₄` 有限差分，必須判死（驗 pipeline）。
禁：不得改 BC、不得改 carry 規則、不得新增形式、不得揀 k_s。
"""
import mpmath as mp

mp.mp.dps = 40

# ---- 44-site（R11 載體，作負控制）
N44 = 44
g44 = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi44 = lambda n, j: mp.sqrt(mp.mpf(2) / (N44 + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N44 + 1))

# ---- 88-site 雙鏈
N88 = 88
g88 = lambda k: 4 * mp.sin(mp.mpf(k) * mp.pi / 176) ** 2
phi88 = lambda k, j: mp.sqrt(mp.mpf(2) / (N88 + 1)) * mp.sin(mp.mpf(k) * mp.pi * j / (N88 + 1))

print("=" * 96)
print("R12 — S2 三攻：88-site 雙鏈之有限差分（先驗形式；N_trials = 2）")
print("=" * 96)
print("  88-site：γ'_k = 4sin²(πk/176), k=1..88；Σ even/odd = %s / %s"
      % (mp.nstr(sum(g88(2 * n) for n in range(1, 45)), 10),
         mp.nstr(sum(g88(2 * n - 1) for n in range(1, 45)), 10)))
print("  even-mode carry 檢查：γ'_{2n} = γ_n⁽⁴⁴⁾ ⟹ max|Δ| = %s"
      % mp.nstr(max(abs(g88(2 * n) - g44(n)) for n in range(1, 45)), 5))

IDX = {}
for i in set([2, 3, 4, 6, 8]):
    for j in range(1, N88 + 1):
        if i <= 8 or j <= 8:
            IDX[(i, j)] = sum(phi88(k, i) * phi88(k, j) / g88(k) for k in range(1, N88 + 1))
G88 = lambda i, j: IDX[(i, j)]

G44 = {}
for i in (2, 3):
    for j in range(1, N44 + 1):
        G44[(i, j)] = sum(phi44(n, i) * phi44(n, j) / g44(n) for n in range(1, N44 + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
DS, DB, DC = mp.mpf("3.4250"), mp.mpf("1.8933"), mp.mpf("-0.2820")

print("\n  G₈₈ 取樣：G(4,4)=%s  G(6,6)=%s  G(6,4)=%s  G(6,8)=%s  G(4,8)=%s"
      % tuple(mp.nstr(G88(*p), 8) for p in [(4, 4), (6, 6), (6, 4), (6, 8), (4, 8)]))
print("  目標：ρ_d = %s ± %s ；ρ_u = %s ± %s" % (mp.nstr(RD, 8), mp.nstr(RD_S, 5),
                                               mp.nstr(RU, 8), mp.nstr(RU_S, 5)))


def run(tag, ku, kd, ks, Gf, main=False):
    du = Gf(kd, ku) - Gf(ku, ku)
    dd = Gf(kd, kd) - Gf(ku, kd)
    ds = Gf(kd, ks) - Gf(ku, ks)
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("\n  [%s] k_u=%d k_d=%d k_s=%d" % (tag, ku, kd, ks))
    print("     δ_u=%s  δ_d=%s  δ_s=%s" % (mp.nstr(du, 8), mp.nstr(dd, 8), mp.nstr(ds, 8)))
    print("     ρ_d = %s（z %+6.2f）   ρ_u = %s（z %+6.2f） ⟹ %s"
          % (mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "不中"))
    return du, dd, ds, pd, pu, zd, zu, ok


print("\n" + "=" * 96)
print("[A1] carry-by-2（O218 偶模遷移：u→4, d→6, s→8）—— 主協議")
a1 = run("A1", 4, 6, 8, G88, main=True)
print("\n[A2] direct（已註冊指數直接當 88-site 指數：u→2, d→3, s→4）")
a2 = run("A2", 2, 3, 4, G88)
print("\n[負控制] G₄₄ 有限差分（R11 已判死之形式）")
c0 = run("CTRL", 2, 3, 4, lambda i, j: G44[(i, j)])

print("\n" + "=" * 96)
print("[判定]")
print("   K1（A1，主協議）：%s" % ("命中 🟢" if a1[7] else "不中 🔴（z_d = %+0.2f）" % a1[5]))
print("   K2（A2，須標 N_trials=2）：%s" % ("命中（但屬第二協議，須報選擇效應）" if a2[7] else "不中"))
print("   負控制：%s" % ("判死 ✓（pipeline 有效）" if ((abs(c0[5]) >= 3) or (abs(c0[6]) >= 3)) else "**反而通 ⟹ pipeline 失效**"))

print("\n" + "=" * 96)
print("[K3] 標度與反解（A1；後驗，不作證據）")
C = DS / a1[2]
print("   C = δ_s/[G₈₈(6,8) − G₈₈(4,8)] = %s MeV" % mp.nstr(C, 8))
for nm, dv in (("δ_b", DB), ("δ_c", DC)):
    tgt = dv / C
    best = min(range(1, 89), key=lambda k: abs((G88(6, k) - G88(4, k)) - tgt) / abs(tgt))
    rel = (G88(6, best) - G88(4, best)) / tgt
    print("   由 %s 反解：最近整數 k = %d（值/目標 = %s）⟹ %s"
          % (nm, best, mp.nstr(rel, 6), "整數解存在" if abs(rel - 1) < mp.mpf("0.05") else "**無整數解**"))

print("\n" + "=" * 96)
print("[匯總 — 只報事實]")
print("   A1（carry-by-2）：ρ_d=%s（z %+0.2f）、ρ_u=%s（z %+0.2f）" % (mp.nstr(a1[3], 8), a1[5], mp.nstr(a1[4], 8), a1[6]))
print("   A2（direct）    ：ρ_d=%s（z %+0.2f）、ρ_u=%s（z %+0.2f）" % (mp.nstr(a2[3], 8), a2[5], mp.nstr(a2[4], 8), a2[6]))
print("   C（A1）= %s MeV" % mp.nstr(C, 8))
print("   ⟹ 判詞：%s" % ("A1 命中 ⟹ 零參數形狀成立（仍需標度）" if a1[7] else
                        "A1 不中 ⟹ 「88-site 雙鏈 + 有限差分」class 判死"))
print("=" * 96)