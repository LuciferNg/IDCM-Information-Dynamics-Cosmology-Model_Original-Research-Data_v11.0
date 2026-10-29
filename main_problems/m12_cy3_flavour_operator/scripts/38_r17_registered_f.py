#!/usr/bin/env python3
"""38 — R17：方向①「註冊新 f」—— 已註冊譜權函數之窮舉

docs/28 之結論：S2 之 closure 條件 = 出現**新結構 f**（T0：w = f(A) 須有來歷），非新核族。
本輪把 v11 charter（O210 docs/03）內**已註冊**且可作譜權之 f 全部列出並逐一測（**非擬合**：清單跑前凍死）。

預註冊清單（N_trials = 5；全部註冊於 charter，無一為臨場造）
--------------------------------------------------
 f1  w = 1/(γ + μ²)      μ² = 2 = γ₄₄（C3）              有質量傳播子
 f2  w = 1/(γ + m_ρ²)    m_ρ² = 4μ² = 8（C3）             ρ 質量傳播子
 f3  w = 1/γ²                                             兩點傳播子（resolvent²）
 f4  w = γ                                                K 本身（Laplacian）
 f5  w = γ(4−γ) = γ_n·γ_{44−n}                            對偶乘積（γ↔4−γ 恆等）
 對照 f0 = 1/γ（R11 已判死 z=+17.96）—— 不計入 N_trials，只作 pipeline 核對

判準（先寫死）
--------------------------------------------------
F1  |z(ρ_d)| < 3 **且** |z(ρ_u)| < 3 ⟹ 命中
F2  全報（不得只報命中者）；命中者須另立選擇性檢定
F3  f 之清單不得事後增刪；不得改 k_s=4
"""
import mpmath as mp

mp.mp.dps = 30
NS = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))

RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
KU, KD, KS = 2, 3, 4
MU2 = mp.mpf(2)          # C3: μ² = 2 = γ₄₄
MRHO2 = mp.mpf(8)        # C3: m_ρ² = 4μ² = 8

print("=" * 96)
print("38 — R17：已註冊譜權函數 f 之窮舉（N_trials = 5，清單跑前凍死）")
print("=" * 96)


def run(tag, w, note=""):
    G = lambda i, j: sum(w(n) * phi(n, i) * phi(n, j) for n in range(1, NS + 1))
    ds = G(KD, KS) - G(KU, KS)
    if ds == 0:
        print("   %-30s δ_s = 0 ⟹ 未定義 ⟹ 判死 %s" % (tag, note))
        return None
    pd, pu = (G(KD, KD) - G(KU, KD)) / ds, (G(KD, KU) - G(KU, KU)) / ds
    zd, zu = (pd - RD) / RD_S, (pu - RU) / RU_S
    ok = (abs(zd) < 3) and (abs(zu) < 3)
    print("   %-30s ρ_d = %-13s（z %+8.2f）  ρ_u = %-13s（z %+8.2f）⟹ %s %s"
          % (tag, mp.nstr(pd, 8), zd, mp.nstr(pu, 8), zu, "★ 命中" if ok else "判死", note))
    return (pd, pu, zd, zu, ok)


res = {}
print("\n[對照] f0（R11 已知判死）")
res["f0"] = run("f0 = 1/γ", lambda n: 1 / gam(n), "（pipeline 核對）")

print("\n[五個已註冊 f]")
res["f1"] = run("f1 = 1/(γ+μ²), μ²=2", lambda n: 1 / (gam(n) + MU2), "（C3）")
res["f2"] = run("f2 = 1/(γ+m_ρ²), m_ρ²=8", lambda n: 1 / (gam(n) + MRHO2), "（C3）")
res["f3"] = run("f3 = 1/γ²", lambda n: 1 / gam(n) ** 2, "（resolvent²）")
res["f4"] = run("f4 = γ", lambda n: gam(n), "（K 本身）")
res["f5"] = run("f5 = γ(4−γ)", lambda n: gam(n) * (4 - gam(n)), "（對偶乘積）")

print("\n" + "=" * 96)
print("[判定]")
hit = [k for k in ("f1", "f2", "f3", "f4", "f5") if res[k] and res[k][4]]
for k in ("f1", "f2", "f3", "f4", "f5"):
    r = res[k]
    print("   %s：%s" % (k, "— 未定義" if r is None else ("★ 命中（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]) if r[4] else "判死（z_d %+0.2f、z_u %+0.2f）" % (r[2], r[3]))))
print("   命中數 = %d ／ 5" % len(hit))
print("   對照 f0：%s" % ("判死 ✓ pipeline 有效" if not (res["f0"] and res["f0"][4]) else "★ 異常"))
print("\n[結論] 五個已註冊 f 之窮舉%s" % ("全數判死 ⟹ 方向① 在此清單內無解" if not hit else "有命中，須另立選擇性檢定"))
print("=" * 96)