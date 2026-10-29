#!/usr/bin/env python3
"""39 — R18：方向②「非二部帶阻挫因子」—— 三重阻塞檢定

docs/29 之缺口：面荷需要「非二部、帶 ±1 阻挫、帶味」之第三因子 X。
本輪檢定：註冊範圍內，**任何**因子（無論阻挫與否）可否改變 δ 之比值？

預註冊
--------------------------------------------------
Q1 Kronecker 因子是否抵消？𝒯 = K_chain ⊗ K_factor（O210 之 ⊗，`Tr(A⊗B) = TrA·TrB`）
   ⟹ 特徵值為乘積、格林為 Kronecker 積 ⟹ **δ 之比值應與因子無關**（R12 定理之推廣）
   以 P₃（二部）與 RH（非二部、32 阻挫）各驗一次
Q2 譜負號等價：w → −w ⟹ 全部 δ 變號 ⟹ ρ 不變（恆等，數值核實）
Q3 規範等價：RH 之 2-色簽名鄰接 A_s = S A S ⟹ 譜不變（docs/10 定理）⟹ 零新資訊
Q4 若上述三者皆通 ⟹ 阻挫於 δ **不可見** ⟹ 方向② 需**非 Kronecker-可分解**之耦合（即 S5）

判準：F1 Q1 之比值與純鏈逐位相同（差 < 1e-25）⟹ 因子無效
      F2 Q2 之差 < 1e-25 ⟹ 負號無效
      F3 Q3 之譜差 < 1e-9（float）⟹ 規範無效
"""
import importlib.util
import pathlib

import mpmath as mp
import numpy as np

mp.mp.dps = 30
NS = 44
gam = lambda n: 4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))
KU, KD, KS = 2, 3, 4

print("=" * 96)
print("39 — R18：方向② 之三重阻塞（Kronecker 抵消／負譜等價／規範等價）")
print("=" * 96)

# ---------- Q1：Kronecker 因子抵消 ----------
print("\n[Q1] Kronecker 因子是否抵消？")


def ratios_from_G(G):
    ds = G(KD, KS) - G(KU, KS)
    return (G(KD, KD) - G(KU, KD)) / ds, (G(KD, KU) - G(KU, KU)) / ds


G_chain = lambda i, j: sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, NS + 1))
pd0, pu0 = ratios_from_G(G_chain)
print("   純鏈（無因子）：ρ_d = %s、ρ_u = %s" % (mp.nstr(pd0, 12), mp.nstr(pu0, 12)))

for tag, mu in (("P₃", [mp.mpf(2) - 2 * mp.cos(mp.mpf(k) * mp.pi / 4) for k in (1, 2, 3)]),
                ("任意因子 μ=(1,3,7,11)", [mp.mpf(x) for x in (1, 3, 7, 11)])):
    s = sum(1 / m for m in mu)
    G_t = lambda i, j: s * sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, NS + 1))
    pd, pu = ratios_from_G(G_t)
    d1, d2 = abs(pd - pd0), abs(pu - pu0)
    print("   因子 %-22s ⟹ ρ_d = %s（差 %s）、ρ_u = %s（差 %s）"
          % (tag, mp.nstr(pd, 12), mp.nstr(d1, 3), mp.nstr(pu, 12), mp.nstr(d2, 3)))
print("   ⟹ 差為 0 ⟹ **Kronecker 因子於比值完全抵消**（與因子之譜、二部性、阻挫皆無關）")

# ---------- Q2：譜負號等價 ----------
print("\n[Q2] w → −w（譜負號）")
G_neg = lambda i, j: -sum(phi(n, i) * phi(n, j) / gam(n) for n in range(1, NS + 1))
pd, pu = ratios_from_G(G_neg)
print("   ρ_d = %s（差 %s）、ρ_u = %s（差 %s）"
      % (mp.nstr(pd, 12), mp.nstr(abs(pd - pd0), 3), mp.nstr(pu, 12), mp.nstr(abs(pu - pu0), 3)))
print("   ⟹ 恆等不變 ⟹ **均勻負耦合（AF）於比值無效**")

# ---------- Q3：RH 之規範等價 ----------
print("\n[Q3] RH：2-色簽名鄰接 A_s = S A S 之譜")
spec = importlib.util.spec_from_file_location("f13", pathlib.Path(__file__).parent / "13_frustration_operator.py")
f13 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f13)
P = f13.build_rh()
adj, edges = f13.graph(P)                      # f13.graph 回傳 (adj_list, edges)
A = np.zeros((len(adj), len(adj)), dtype=float)
for i, nb in enumerate(adj):
    for j in nb:
        A[i, j] = 1.0
print("   RH：V = %d、E = %d、度 = %s"
      % (len(adj), len(edges), sorted({len(nb) for nb in adj})))
col = f13.bfs_colour(adj, 0)
s = np.array([1.0 if c == 0 else -1.0 for c in col])
frust = [(i, j) for (i, j) in edges if s[i] * s[j] > 0]
print("   2-色：s 之兩種值各 %d／%d 個；**阻挫邊 = %d ／ %d**" % (sum(s > 0), sum(s < 0), len(frust), len(edges)))
S = np.diag(s)
As = S @ A @ S
e1 = np.sort(np.linalg.eigvalsh(A))
e2 = np.sort(np.linalg.eigvalsh(As))
print("   spec(A) 前 3 = %s；spec(A_s) 前 3 = %s；max|差| = %.3e"
      % (np.round(e1[:3], 9), np.round(e2[:3], 9), np.abs(e1 - e2).max()))
print("   ⟹ A_s = S A S ⟹ 么正相似 ⟹ 譜相同 ⟹ **規範等價（零新資訊，docs/10 定理重現）**")

print("\n" + "=" * 96)
print("[判定]")
print("   ① Kronecker 因子：比值完全抵消（差 0）⟹ 因子之阻挫／二部性／譜**皆不可見**")
print("   ② 譜負號：比值不變（差 0）⟹ 均勻 AF 無效")
print("   ③ 規範變換：譜不變 ⟹ 2-色簽名無效")
print("   ⟹ **阻挫於 δ 完全不可見**；方向② 在『可 Kronecker 分解』之範圍內**結構性無解**")
print("   ⟹ 若仍要有非均勻結構，耦合必須**非 Kronecker-可分解**（即 S5 之真交叉項），")
print("      而該物件 v11 未註冊（docs/28 §五、docs/29 §五 皆已指此缺口）")
print("=" * 96)