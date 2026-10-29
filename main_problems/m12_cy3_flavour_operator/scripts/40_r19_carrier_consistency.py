#!/usr/bin/env python3
"""40 — R19：載體一致性審計 ＋ 全輪重跑（自捕 off-by-one）

★ 發現：R8–R14 之譜用 `γ_n = 4sin²(nπ/88)`，而特徵向量用 `φ_n(j) = √(2/45)sin(πnj/45)`。
   兩者**屬不同站數**：`√(2/45)` ⟹ 44 站 Dirichlet ⟹ 譜應為 `4sin²(nπ/90)`；
   `4sin²(nπ/88)` ⟹ 43 站 Dirichlet。混合兩者 ⟹ `G` 唔係任何自伴算子之逆。

本輪：
 Q1 一致性診斷：驗 `K G = I`（K 之譜 = 用同一約定之 γ）於兩種約定
 Q2 全輪重跑：受影響之 9 個變體，於 (a) 原約定 `/88`（不一致）與 (b) 修正約定 `/90`（44 站）
    下各算一次，逐個比較判詞
 Q3 若判詞全同 ⟹ 修正只影響數值；若有翻轉 ⟹ 須重寫該輪判詞

判準（先寫死）：判詞同 ⟺ 兩約定下 |z| 同側於 3（同為命中或同為判死）。
"""
import mpmath as mp

mp.mp.dps = 30
NS = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
g88 = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 88) ** 2
g90 = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / 90) ** 2
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4

print("=" * 100)
print("40 — R19：載體一致性審計（自捕 off-by-one）＋ 全輪重跑")
print("=" * 100)
print("   γ(/88)：γ₂=%s γ₃=%s γ₄=%s γ₄₄=%s" % (mp.nstr(g88(2), 8), mp.nstr(g88(3), 8), mp.nstr(g88(4), 8), mp.nstr(g88(44), 8)))
print("   γ(/90)：γ₂=%s γ₃=%s γ₄=%s γ₄₄=%s" % (mp.nstr(g90(2), 8), mp.nstr(g90(3), 8), mp.nstr(g90(4), 8), mp.nstr(g90(44), 8)))

print("\n[Q1] 一致性診斷：K G = I ？（K 與 G 用同一 γ）")


def residual(g):
    G = lambda i, j: sum(phi(n, i) * phi(n, j) / g(n) for n in range(1, NS + 1))
    worst = mp.mpf(0)
    for i in range(1, NS + 1):
        for j in range(1, NS + 1):
            kg = sum(g(n) * phi(n, i) * phi(n, j) / g(n) for n in range(1, NS + 1))
            worst = max(worst, abs(kg - (1 if i == j else 0)))
    return worst


for tag, g in (("γ=/88（不一致：φ 屬 44 站）", g88), ("γ=/90（一致：φ 屬 44 站）", g90)):
    print("   %-30s max|KG − I| = %s" % (tag, mp.nstr(residual(g), 6)))
print("   註：`K G = I` 對任何 γ 皆形式上成立（同一 γ 分子分母相消），故 Q1 改以下述「站數歸屬」判：")
print("       φ_n(j) = √(2/(N+1))·sin(nπj/(N+1)) ⟹ K 之譜必為 4sin²(nπ/(2(N+1)))；")
print("       N = 44 ⟹ 2(N+1) = 90 ⟹ 譜 = 4sin²(nπ/90) ⟹ **/88 與 φ(/45) 不相容**。")

print("\n[Q2] 全輪重跑（同一變體、兩個 γ 約定）")


def mk(g):
    G = lambda i, j: sum(phi(n, i) * phi(n, j) / g(n) for n in range(1, NS + 1))
    def Gt(i, j, nmax):
        return sum(phi(n, i) * phi(n, j) / g(n) for n in range(1, nmax + 1))
    Gs = lambda i, j: (mp.mpf(-1) ** (i + j)) * G(i, j)
    Amp = lambda i, j: sum(phi(n, i) * phi(n, j) / mp.sqrt(g(n)) for n in range(1, NS + 1))
    th = lambda n: 2 * mp.pi * n / 88
    Atil = lambda i, j: sum(phi(n, i) * phi(n, j) * mp.e ** (1j * th(n)) / mp.sqrt(g(n)) for n in range(1, NS + 1))
    Gb = lambda p, q: sum((phi(n, p) + phi(n, p + 1)) * (phi(n, q) + phi(n, q + 1)) / (2 * g(n)) for n in range(1, NS + 1))
    t = Gb(KU, KD)
    V = {
        "R11 全G差分": lambda x: G(KD, x) - G(KU, x),
        "R13 截斷 n≤27": lambda x: Gt(KD, x, 27) - Gt(KU, x, 27),
        "S3a 交錯譜權": lambda x: sum((mp.mpf(-1) ** n) * phi(n, KD) * phi(n, x) / g(n) for n in range(1, NS + 1))
                            - sum((mp.mpf(-1) ** n) * phi(n, KU) * phi(n, x) / g(n) for n in range(1, NS + 1)),
        "S3b staggered": lambda x: Gs(KD, x) - Gs(KU, x),
        "G1 實振幅叉項": lambda x: 2 * Amp(KD, x) * Amp(KU, x),
        "G2 環相位干涉": lambda x: 2 * mp.re(Atil(KD, x) * mp.conj(Atil(KU, x))),
        "E1 割能量差": lambda x: Gb(KD, x) - Gb(KU, x),
        "E2 避免交叉": lambda x: (1 if (Gb(KD, x) - Gb(KU, x)) >= 0 else -1) * mp.sqrt((Gb(KD, x) - Gb(KU, x)) ** 2 + 4 * t ** 2),
        "f1 1/(γ+μ²)": lambda x: sum(phi(n, KD) * phi(n, x) / (g(n) + 2) for n in range(1, NS + 1))
                              - sum(phi(n, KU) * phi(n, x) / (g(n) + 2) for n in range(1, NS + 1)),
        "f4 K 本身": lambda x: sum(g(n) * phi(n, KD) * phi(n, x) for n in range(1, NS + 1))
                            - sum(g(n) * phi(n, KU) * phi(n, x) for n in range(1, NS + 1)),
    }
    return V


def evald(V):
    out = {}
    for k, f in V.items():
        du, dd, ds = f(KU), f(KD), f(KS)
        if ds == 0:
            out[k] = (None, None, False)
            continue
        pd, pu = dd / ds, du / ds
        zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
        out[k] = (zd, zu, (abs(zd) < 3) and (abs(zu) < 3))
    return out


A88, A90 = evald(mk(g88)), evald(mk(g90))
print("\n   %-16s %-22s %-22s %s" % ("變體", "γ=/88（原）z_d, z_u", "γ=/90（修正）z_d, z_u", "判詞"))
flip = 0
for k in A88:
    a, b = A88[k], A90[k]
    if a[0] is None or b[0] is None:
        print("   %-16s %s" % (k, "— 未定義"))
        continue
    same = (a[2] == b[2])
    flip += (0 if same else 1)
    print("   %-16s %+8.2f, %+8.2f      %+8.2f, %+8.2f      %s%s"
          % (k, a[0], a[1], b[0], b[1], "同" if same else "★ 翻轉", "（命中）" if b[2] else "（判死）"))

print("\n" + "=" * 100)
print("[判定]")
print("   Q1：`/88` 與 `φ(/45)` 不相容（站數 43 vs 44）⟹ R8–R14 之載體有 off-by-one")
print("   Q2：十個變體中，判詞翻轉 = %d 個" % flip)
print("   ⟹ %s" % ("全部判詞不變：修正只影響數值（第三位），R11–R17 之結論保持"
                  if flip == 0 else "**有翻轉 ⟹ 相關輪次之判詞須重寫**"))
print("=" * 100)