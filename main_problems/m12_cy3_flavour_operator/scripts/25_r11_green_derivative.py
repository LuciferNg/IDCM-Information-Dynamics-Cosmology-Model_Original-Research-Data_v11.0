#!/usr/bin/env python3
"""R11 — S2 第二攻：δ 作為「模式指標之有限差分」（N_trials = 2，先驗形式）

結構動機（先於計算，非後驗）
--------------------------------------------------
1. 可分離性排除：若對荷 a(k,x) 可分離（f(k)+g(x) 或 f(k)g(x)），則
   δ_x = a(3,x) − a(2,x) 與 x 無關；實測 δ_d:δ_s:δ_b = 1.382 : 3.425 : 1.893
   （變動 2.5×）⟹ **可分離全死**（代數陳述，不依賴任何形式）。
2. 故 a 必須非可分離 ⟹ δ_x = a(3,x) − a(2,x) 係 a 對**模式指標之有限差分**，
   即 δ = ⟨∂_k a⟩（離散導數）。
3. 載體唯一提供之兩點關聯子 = K_chain⁻¹（Dirichlet Green 函數）。故最簡先驗形式：
        δ_x = C·[ G(3, k_x) − G(2, k_x) ] ,     G = K_chain⁻¹
   注意此形式**唔係** R4-P2 所測（該處測 G 之絕對值對分割比；此處測 G 對 k 之導數）。
4. 定標 C 在比值中相消 ⟹ ρ_x 係**零參數**預言。

預註冊判準（跑前凍結，不得事後改）
--------------------------------------------------
K1  主協議 P1（G 之有限差分）：|z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中。
K2  負控制 P2（已知死形式 F1：g(k₁,k₂)=γ_{|k₁−k₂|} 之有限差分，scripts/23 已判不合）：
    必須**判死**（至少一項 |z| ≥ 3）⟹ 證明本 harness 有分辨力；若 P2 反而通，則
    本 harness 失效，P1 之任何結果不得引用。
K3  標度：報 C = δ_s/[G(3,4) − G(2,4)]；並由 δ_b、δ_c 反解其等效 k（非擬合後果）。
K4  k_s 敏感度：{4,5,7} 全報，但**只有 k_s=4（已註冊主方案）計入 K1**；其餘作敏感度。
禁用：不得加入任何自由函數；不得改 BC；不得揀 k_s；不得因結果而換協議。
"""
import mpmath as mp

mp.mp.dps = 40

N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))

# --- 目標（R7 錨解；σ 線性傳播）
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
DS = mp.mpf("3.4250")
DB = mp.mpf("1.8933")
DC = mp.mpf("-0.2820")
KU, KD = 2, 3
KS_MAIN = 4
KS_SENS = (4, 5, 7)

# --- K_chain⁻¹（Dirichlet 44-site；譜 γ_n = 4sin²(πn/88)）
G = {}
for i in range(1, N + 1):
    for j in range(1, N + 1):
        if i <= 4 or j <= 4:
            G[(i, j)] = sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, N + 1))

print("=" * 96)
print("R11 — S2 第二攻：δ = 模式指標之有限差分（N_trials = 2，先驗形式）")
print("=" * 96)
print("  K_chain：44-site Dirichlet，γ_n = 4sin²(πn/88)，φ_n(j)=√(2/45)sin(πnj/45)")
print("  G(i,j) 取樣：G(2,2)=%s  G(3,3)=%s  G(3,2)=%s  G(3,4)=%s"
      % (mp.nstr(G[(2, 2)], 8), mp.nstr(G[(3, 3)], 8), mp.nstr(G[(3, 2)], 8), mp.nstr(G[(3, 4)], 8)))
print("  目標：ρ_d = %s ± %s ；ρ_u = %s ± %s" % (mp.nstr(RD, 8), mp.nstr(RD_S, 5),
                                               mp.nstr(RU, 8), mp.nstr(RU_S, 5)))

def report(tag, d_of, ks_list, main_ks):
    """d_of(kx) -> δ_x（未定標）；回傳主協議判定所需之 (ρ_d, ρ_u, z_d, z_u)"""
    out = {}
    den = d_of(main_ks) if main_ks is not None else None
    print("\n" + "-" * 96)
    print("[%s]" % tag)
    for ks in ks_list:
        den = d_of(ks)
        pd = d_of(KD) / den
        pu = d_of(KU) / den
        zd = (pd - RD) / RD_S
        zu = (pu - RU) / RU_S
        mark = ""
        if ks == main_ks:
            ok = (abs(zd) < 3) and (abs(zu) < 3)
            mark = "  ★ K1 命中" if ok else "  ✗ K1 不中"
            out = (pd, pu, zd, zu, ok)
        print("   k_s=%d: δ_u=%s δ_d=%s δ_s=%s ⟹ ρ_d = %s（z %+6.2f）  ρ_u = %s（z %+6.2f）%s"
              % (ks, mp.nstr(d_of(KU), 8), mp.nstr(d_of(KD), 8), mp.nstr(den, 8),
                 mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, mark))
    return out

# ---- P1：G 之有限差分（主協議）
dG = lambda kx: G[(KD, kx)] - G[(KU, kx)]
p1 = report("P1 主協議：δ_x ∝ G(3,k_x) − G(2,k_x)", dG, KS_SENS, KS_MAIN)

# ---- P2：負控制（已知死形式 F1 之有限差分）
g1 = lambda n: gam(n)
dF1 = lambda kx: g1(abs(KD - kx)) - g1(abs(KU - kx))
p2 = report("P2 負控制：δ_x ∝ γ_|3−k_x| − γ_|2−k_x|（scripts/23 之 F1，已知不合）", dF1, KS_SENS, KS_MAIN)

print("\n" + "=" * 96)
print("[K2] 負控制判定")
dead2 = (abs(p2[2]) >= 3) or (abs(p2[3]) >= 3)
print("   P2：z(ρ_d) = %+0.2f、z(ρ_u) = %+0.2f ⟹ %s"
      % (p2[2], p2[3], "判死 ✓（harness 有分辨力）" if dead2 else "**反而通 ⟹ harness 失效，P1 不得引用**"))
print("   ⟹ K2：%s" % ("通過" if dead2 else "不通過"))

print("\n" + "=" * 96)
print("[K3] 標度與非擬合後果（P1）")
C = DS / dG(KS_MAIN)
print("   C = δ_s / [G(3,4) − G(2,4)] = %s MeV" % mp.nstr(C, 8))
print("   標度結構性：C 須由載體標度給出；M₀ 係 O38 開放項 ⟹ 仍屬同一個未解標度桶")

def implied_k(dv):
    """由 δ = C[G(3,k) − G(2,k)] 反解 k（在 k∈[1,44] 上單調搜尋）"""
    tgt = dv / C
    lo, hi = 1, 44
    def f(k):
        # 連續化：以 sin² 形式之連續 G 不方便；改用最近整數 k 之插值
        k = max(mp.mpf(1), min(mp.mpf(44), k))
        k0 = int(mp.floor(k)); k1 = min(k0 + 1, 44)
        w = k - k0
        v0 = G[(KD, k0)] - G[(KU, k0)]
        v1 = G[(KD, k1)] - G[(KU, k1)]
        return (1 - w) * v0 + w * v1
    # 掃描
    best = min(range(1, 45), key=lambda k: abs(f(mp.mpf(k)) - tgt))
    return best, f(mp.mpf(best)) / tgt

for nm, dv in (("δ_b", DB), ("δ_c", DC)):
    k, ratio = implied_k(dv)
    print("   由 %s = %s 反解：最近整數 k = %d（該點值/目標 = %s）⟹ %s"
          % (nm, mp.nstr(dv, 6), k, mp.nstr(ratio, 6),
             "整數解存在" if abs(ratio - 1) < mp.mpf("0.05") else "**無整數解**（偏離 >5%）"))

print("\n" + "=" * 96)
print("[匯總 — 只報事實]")
print("   K1（P1, k_s=4）：ρ_d = %s（z %+0.2f）、ρ_u = %s（z %+0.2f）⟹ %s"
      % (mp.nstr(p1[0], 8), p1[2], mp.nstr(p1[1], 8), p1[3], "命中" if p1[4] else "不中"))
print("   K2（負控制 P2）：%s" % ("通過（已知死形式確實判死）" if dead2 else "不通過"))
print("   K3：C = %s MeV（無結構表達）" % mp.nstr(C, 8))
print("   ⟹ 判詞：%s" % ("P1 命中且 harness 有效 ⟹ 零參數形狀成立（仍需 K3 標度）"
                        if (p1[4] and dead2) else
                        "P1 不中 ⟹ 「G 之有限差分」class 死；S2 收窄"))
print("=" * 96)