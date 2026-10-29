#!/usr/bin/env python3
"""R10 — S2：δ 之 isovector 算子身份（形狀族掃描 ＋ 選擇性控制 ＋ 可實現性閘）

背景
----
R4–R9 後唯一未被盲測之靶（`final_result.md §殘餘 S2`）：
  δ_x ≡ a(d,x) − a(u,x) 係**兩構型之差**，等價於「翻轉 k_u↔k_d 之矩陣元」。
已註冊指數：u→2、d→3；s 三讀 {4,5,7}（未定）。c、b 未註冊。

本輪形式族（★ 後驗，見 §H 誠實標記）
  H0 形狀假設：δ_x = C·[ f(γ_{k_x}) − f(γ_{k_u}) ]
  ⟹ ρ_x ≡ δ_x/δ_s = [f(γ_{k_x}) − f(γ_{k_u})] / [f(γ_{k_s}) − f(γ_{k_u})]
  f ∈ {γ, √γ, γ², 1/γ, 4−γ, ln γ}（6 個）；k_s ∈ {4,5,7}；共 18 組合。
  注意：此形狀令 δ_u ≡ 0（因 k_x = k_u），故 ρ_u 之預言恆為 0。

預註冊判準（跑前凍結）
--------------------------------------------------
P1 命中：|z(ρ_d)| < 3 **且** |z(ρ_u)| < 3。
P2 選擇性（**證據之必要條件**）：命中組合數 ≤ 2（18 中），且空檢定（10⁴ 亂數 target）
   P(命中數 ≥ 觀測值) < 0.05；否則「命中」不構成證據。
P3 可實現性：須存在對稱對荷 a(k₁,k₂)=a(k₂,k₁) 使 a(3,x)−a(2,x) 等於該式。若不存在 ⟹
   該形狀唔係「對荷之差」之陳述，只係 δ 之 fit。
P4 標度：C ≡ δ_s / [f(γ_{k_s}) − f(γ_{k_u})] 必須有結構表達；否則只得「形狀」無「算子」。
判死門檻 |z| ≥ 3。禁：不得由本輪回擬任何 δ；不得挑 k_s。
"""
import mpmath as mp

mp.mp.dps = 40

# λ 譜
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
KU, KD = 2, 3
KS_LIST = (4, 5, 7)

# 目標（R7 錨解；σ 線性傳播）
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")

FS = {
    "γ":       (lambda g: g,           "1"),
    "√γ":      (lambda g: mp.sqrt(g),  "2"),
    "γ²":      (lambda g: g ** 2,      "3"),
    "1/γ":     (lambda g: 1 / g,       "4"),
    "4−γ":     (lambda g: 4 - g,       "5"),
    "ln γ":    (lambda g: mp.log(g),   "6"),
}

print("=" * 96)
print("R10 — δ 之 isovector 算子身份：形狀族掃描 ＋ 選擇性控制 ＋ 可實現性閘")
print("=" * 96)
print("  γ_n = 4sin²(πn/88) ；γ_1..γ_7 = %s" % [mp.nstr(gam(n), 7) for n in range(1, 8)])
print("  目標：ρ_d = %s ± %s ；ρ_u = %s ± %s" % (mp.nstr(RD, 8), mp.nstr(RD_S, 5),
                                               mp.nstr(RU, 8), mp.nstr(RU_S, 5)))

print("\n" + "=" * 96)
print("[T1] 形狀族掃描（ρ_u 之預言恆為 0，構造使然）")
print("  %-6s %-4s %10s %9s %8s %10s %9s %8s %s" %
      ("f", "k_s", "預言 ρ_d", "z", "z(ρ_u)", "命中", "", "", ""))
hits = []
for fname, (f, _) in FS.items():
    for ks in KS_LIST:
        den = f(gam(ks)) - f(gam(KU))
        pd = (f(gam(KD)) - f(gam(KU))) / den
        zd = (pd - RD) / RD_S
        zu = (mp.mpf(0) - RU) / RU_S
        hit = (abs(zd) < 3) and (abs(zu) < 3)
        if hit:
            hits.append((fname, ks, pd, zd))
        print("  %-6s %-4d %10s %+9.2f %8.2f %10s"
              % (fname, ks, mp.nstr(pd, 8), zd, zu, "★ 命中" if hit else "—"))

print("\n  命中組合 = %d / 18" % len(hits))
for fname, ks, pd, zd in hits:
    print("     %s（k_s=%d）：ρ_d = %s（目標 %s，z = %+.2f）" % (fname, ks, mp.nstr(pd, 8), mp.nstr(RD, 8), zd))

print("\n" + "=" * 96)
print("[T2] 空檢定：target 換成亂數（ρ_d ∈ [0.10,0.90]，ρ_u ∈ [−0.30,0.30]），10⁴ 次")
import random
rng = random.Random(20261029)
grid = []
for fname, (f, _) in FS.items():
    for ks in KS_LIST:
        den = f(gam(ks)) - f(gam(KU))
        grid.append(float((f(gam(KD)) - f(gam(KU))) / den))
cnt = []
for _ in range(10000):
    td = rng.uniform(0.10, 0.90)
    tu = rng.uniform(-0.30, 0.30)
    c = 0
    for v in grid:
        # 用固定 σ 作尺度（與 T1 相同門檻 3σ）
        if abs(v - td) / float(RD_S) < 3 and abs(0 - tu) / float(RU_S) < 3:
            c += 1
    cnt.append(c)
p_ge = sum(1 for c in cnt if c >= max(len(hits), 1)) / len(cnt)
print("   每 target 命中數：均值 %.2f ；P(≥%d) = %.3f ；P(≥1) = %.3f"
      % (sum(cnt) / len(cnt), max(len(hits), 1), p_ge, sum(1 for c in cnt if c >= 1) / len(cnt)))
print("   ⟹ P2 選擇性：%s" % ("通過（P < 0.05）" if (p_ge < 0.05 and len(hits) <= 2) else "**不通過**（命中可由偶然產生）"))

print("\n" + "=" * 96)
print("[T3] 可實現性（P3）：是否存在對稱對荷 a 使 a(3,x)−a(2,x) = C[f(γ_x)−f(γ_u)] ?")
print("   代數：令 a(3,x) := a(2,x) + C[f(γ_x)−f(γ_u)]。對稱性只需 a(3,2) = a(2,3)，")
print("   而 x=2 給 a(3,2) = a(2,2) + C[f(γ_2)−f(γ_2)] = a(2,2)；a(2,3) 由 x=3 之定義給")
print("   a(3,3) = a(2,3) + C[f(γ_3)−f(γ_2)] ⟹ 一致，**故存在**。")
print("   但 a(2,x) 對每個 x 皆自由（44 個自由度）⟹ 該形狀**不決定** a，只約束其差。")
print("   ⟹ P3：形式可實現但**非強制**；唔可以由『對荷結構』推出此形狀。")

print("\n" + "=" * 96)
print("[T4] 標度（P4）：C = δ_s / [f(γ_{k_s}) − f(γ_{k_u})]")
DS = mp.mpf("3.4250")
for fname, ks, pd, zd in hits:
    f = FS[fname][0]
    C = DS / (f(gam(ks)) - f(gam(KU)))
    print("   %s（k_s=%d）：C = %s MeV" % (fname, ks, mp.nstr(C, 8)))
print("   結構檢查：C 之絕對值須由質量標度給出；v11 之質量標度 M₀ 本身係 O38 開放項")
print("   （`M₀ 絕對數值 = 單位制 artifact 偽命題`）⟹ **P4 不通過**：只得形狀，無算子標度。")

print("\n" + "=" * 96)
print("[T5] 非擬合後果（只就 f = γ 之形狀，且去重：γ 與 4−γ 給同一比值）")
uniq = {}
for fname, ks, pd, zd in hits:
    key = mp.nstr(pd, 10)
    uniq.setdefault(key, []).append((fname, ks))
print("   命中組合去重後 = %d 個相異比值（γ 與 4−γ 之比值恆同 ⟹ 非獨立形式）" % len(uniq))
for key, lst in uniq.items():
    print("     比值 %s ：%s" % (key, "、".join("%s(k_s=%d)" % (a, b) for a, b in lst)))

f = FS["γ"][0]
C = DS / (f(gam(4)) - f(gam(KU)))
print("\n   就 f = γ、k_s = 4（C = %s MeV）：" % mp.nstr(C, 8))
for nm, dv in (("δ_b", mp.mpf("1.8933")), ("δ_c", mp.mpf("-0.2820"))):
    tgt = f(gam(KU)) + dv / C
    k = 88 / mp.pi * mp.asin(mp.sqrt(tgt / 4))
    print("     由 %s 反解：γ_k = %s ⟹ k = %s（**非整數**）" % (nm, mp.nstr(tgt, 8), mp.nstr(k, 6)))
print("     k_s 三讀：%s"
      % "；".join("k_s=%d → %s（%+0.1f%%）"
                  % (k2, mp.nstr((f(gam(KD)) - f(gam(KU))) / (f(gam(k2)) - f(gam(KU))), 8),
                     float((((f(gam(KD)) - f(gam(KU))) / (f(gam(k2)) - f(gam(KU)))) - RD) / RD * 100))
                  for k2 in KS_LIST))

print("\n" + "=" * 96)
print("★ H 誠實標記：H0 之『參考點 = γ_{k_u}』係**後驗**發現")
print("   （原為解 ρ_d = [γ_x−c]/[γ_s−c] 得 c = 0.021838，注意到 c ≈ γ_2 = 0.020357）。")
print("   故本輪全部結果屬**後驗形式選擇**；P2 選擇性檢定即為此而設。")
print("=" * 96)
print("\n[匯總 — 只報事實]")
print("   命中：%d/18（清單見 T1）" % len(hits))
print("   P1 命中：%s" % ("通過" if hits else "無"))
print("   P2 選擇性：%s" % ("通過" if (p_ge < 0.05 and len(hits) <= 2) else "不通過"))
print("   P3 可實現但非強制：不通過（形狀不決定 a）")
print("   P4 標度：不通過（C 無結構表達，繫於 O38）")
print("   ⟹ 本輪只得『形狀候選 🟡』，未達算子身份。")
print("=" * 96)