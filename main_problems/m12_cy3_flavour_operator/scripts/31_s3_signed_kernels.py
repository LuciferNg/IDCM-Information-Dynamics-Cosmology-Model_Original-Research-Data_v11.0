#!/usr/bin/env python3
"""31 — S 路線（最舊之殘餘）：S3 帶號核相關能（可算部分）

背景：S1–S5 之中，S3「帶號核」在 R7-A1 判死「sign-definite 族」之後**從未測過**。
      S1（二體對荷）已達可測上限（10 通道，R9）；S4（非 ℤ₂ 定向）與 S5（非加性）
      分別阻塞於 c／b 模式指數未註冊、與無可算載體 ⟹ 本輪**只執行 S3**。

預註冊判準（跑前凍死，不得事後改）
--------------------------------------------------
F1：|z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中（符號須自然吻合，不另設門檻）
F2：`N_trials = 3`（三個帶號核變體，全報；只中一個須報選擇效應）
F3：負控制 = R13-A1（`n≤27` 截斷，已判死 `z=+10.35`）**必須仍判死**
禁：不得加入新變體、不得改 k_s、不得用 c／b（未註冊）。
"""
import mpmath as mp

mp.mp.dps = 40
N = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (N + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (N + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

print("=" * 96)
print("31 — S 路線：S3 帶號核相關能（N_trials = 3，含負控制）")
print("=" * 96)


def mkG(w):
    """G_w(i,j) = Σ_n w(n) φ_n(i)φ_n(j) ；w 決定核之帶號結構"""
    def G(i, j):
        return sum(w(n) * phi(n, i) * phi(n, j) for n in range(1, N + 1))
    return G


def run(tag, G, note=""):
    du = G(KD, KU) - G(KU, KU)
    dd = G(KD, KD) - G(KU, KD)
    ds = G(KD, KS) - G(KU, KS)
    if ds == 0:
        print("   %-40s δ_s = 0 ⟹ 比值未定義 ⟹ 判死 %s" % (tag, note))
        return None
    pd, pu = dd / ds, du / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-40s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


print("\n[A] 帶號核變體（S3）")
r1 = run("S3a 交錯譜權 (−1)^n/γ_n", mkG(lambda n: (mp.mpf(-1) ** n) / gam(n)),
         "（譜權帶號 ⟹ G 非 sign-definite）")
Gplain = mkG(lambda n: 1 / gam(n))
Gs = lambda i, j: (mp.mpf(-1) ** (i + j)) * Gplain(i, j)      # ← 真 staggered：site 空間加號
r2 = run("S3b staggered (−1)^{i+j}G", Gs,
         "（δ = (−1)^x[G(3,x)+G(2,x)] ⟹ 和型而非差型）")
try:
    r3 = run("S3c 對偶核 1/(4−γ_n)", mkG(lambda n: 1 / (4 - gam(n))),
             "（γ↔4−γ 對偶，k↔88−k）")
except ZeroDivisionError:
    r3 = None
    print("   S3c 對偶核 1/(4−γ_n)               **結構奇異：γ_44 = 4 ⟹ 1/(4−γ) 於 Nyquist 有極點**")
    print("                                     ⟹ 該核非良定義（4−K 奇異）⟹ **判死**（不引入任何正則化）")

print("\n[B] 負控制（R13-A1：n≤27 截斷，已判死 z=+10.35）")
Gc = mkG(lambda n: 1 / gam(n))
Gtrunc = lambda i, j: sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, 28))
c = run("CTRL n≤27（頂點差分）", Gtrunc)

print("\n" + "=" * 96)
print("[判定]")
for tag, r in (("S3a", r1), ("S3b", r2), ("S3c", r3)):
    print("   %s：%s" % (tag, "— 未定義" if r is None else ("★ 命中" if r[4] else "判死（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]))))
print("   負控制：%s" % ("判死 ✓ pipeline 有效" if (c is None or not c[4]) else "**反而通 ⟹ pipeline 失效**"))
print("\n[未執行之 S 項（誠實標記）]")
print("   S1 二體對荷：10 通道已通（R9），僅描述性、非 F̂ ⟹ 無新可測通道（b 側需 k_b，未註冊）")
print("   S4 非 ℤ₂ 定向：阻塞於 c／b 模式指數未註冊 ⟹ 不可算")
print("   S5 非加性關係型：無可算載體 ⟹ 不可算（Gemini 讀法見 docs/24、E 路線見 docs/26）")
print("=" * 96)