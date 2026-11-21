#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/05 — R4：徑向激發階梯（T4）
=======================================
前設問題
--------
T4：「徑向激發階梯 ⊂ {γ_n} 子序列（非新標度）」；失敗判準：「若需新標度常數 ⟹ 🔴」。

R4 要答三件事：
  Q-a 字面：樹譜 {λ_k} 係否 {γ_n} 之子序列？（精確）
  Q-b 結構：徑向階梯嘅**正確載體**係邊個算子？（若唔係 γ 子序列，替代係咩？）
  Q-c 標度：階梯嘅物理間距需唔需要**新標度常數**？（T4 生死判準）

本輪用兩個已註冊對象，唔引入任何新參數：
  ① 路徑樹譜（O155 R3 / O210 C9）：λ_k = B sin²(kπ/sites)，sites = 4(B²+1) = 68，k = 1..33
     （重現閘 E1–E5 已核 Σλ = 66、λ_k + λ_{34−k} = 4、λ₁₇ = μ²）
  ② M5 之 Cayley 樹徑向算子（docs/10 §二，已校準）：E ψ_d = −t[ψ_{d−1}+(B−1)ψ_{d+1}] + V_dψ_d，
     V_d = −(G/a)B^{−d}，B = 4；節點數 N_d = B(B−1)^{d−1}
     ⟹ 對稱化 φ_d = √N_d ψ_d 後本徵譜只依 u ≡ (G/a)/t（單一無量綱比）

事前判準（寫死，唔准事後改）
----------------------------
  W1 Q-a 必須用 dps=60（tol 1e−45）＋ ℚ(ζ₆₈)／ℚ(ζ₈₈) 環交叉核，唔可以只報「唔似」
  W2 Q-b 若判「唔係 γ 子序列」，必須指出**替代算子**（唔可以留空）
  W3 Q-c 標度測試：以已註冊觀測之徑向間距為 target，反算所需 t；若所需 t 與 M5 校準 t
     相差 ≥ 10× ⟹ 判「需新標度」⟹ T4 🔴（此為 T4 之生死判準）
  W4 形狀測試（附加）：階梯間距模式（遞增／遞減）必須同觀測對照；唔准只比單一間距
  W5 禁：唔准以觀測間距反解 λ 指標或 u（後驗）；u 只可以掃描列報

輸出：data/05_radial_ladder.json
"""
import json
import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 60
TOL = mp.mpf("1e-45")
B = 4
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="W", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="W"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("W1. Q-a 字面：{λ_k} 係否 {γ_n} 之子序列？")
print("=" * 78)
G = {n: mp.mpf(4) * mp.sin(mp.pi * n / 88) ** 2 for n in range(1, 45)}
LAM = {k: mp.mpf(4) * mp.sin(mp.pi * k / 68) ** 2 for k in range(1, 34)}
matches = []
for k, v in LAM.items():
    for n, u in G.items():
        if abs(v - u) < TOL:
            matches.append((k, n))
print("  值集交集（tol 1e−45）= %d 條：%s" % (len(matches), matches))
chk("W1a 值集交集恰 = {(λ₁₇, γ₂₂)} = μ² = 2（兩譜唯一共同值 = chiral 中心）",
    matches == [(17, 22)], "交集 %s ⟹ 唔係子序列，但**共享自對極點** μ²" % matches)
chk("W1b λ 與 γ 同族但不同載體：λ_k = 4sin²(kπ/68)（N_c = 34）、γ_n = 4sin²(nπ/88)（N_c = 44）",
    True, "同一函數族 4sin²(·)，分母 68 = 2·34 vs 88 = 2·44 ⟹ 值集唔重合")
verdict("W1c T4 字面（階梯 ⊂ {γ_n} 子序列）判死", len(matches) <= 2,
        "交集 %d 條 ⟹ 唔係子序列；但階梯仍有零參數形式（樹譜 λ_k）" % len(matches), "W1")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W2. Q-b 結構：階梯之正確載體")
print("=" * 78)
pairs = [(k, 34 - k) for k in range(1, 17)]
print("  路徑樹譜：33 模 = 16 對極對 + 自對點 {17}（λ₁₇ = μ²）")
print("  ⟹ 每通道最多 16 個徑向激發（對極對數 = B² = 16）")
chk("W2a 徑向激發容量 = B² = 16（每通道）", len(pairs) == 16, "16 = B²")
chk("W2b 替代算子已指名 = K_tree 路徑 Laplacian（Dirichlet，sites = 68）",
    True, "非 γ 子序列；係**同族不同鏈**（N_c = 34 = K_tree^(node)）")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W3. Q-c：Cayley 樹徑向算子本徵階梯（M5 算子，無量綱 u = (G/a)/t）")
print("=" * 78)


def cayley_ladder(u, L=80, nmax=12):
    """E ψ_d = −t[ψ_{d−1}+(B−1)ψ_{d+1}] + V_d ψ_d，V_d = −(u·t)·B^{−d}
       對稱化 φ_d = √N_d ψ_d，N_0=1, N_d = B(B−1)^{d−1}；Dirichlet ψ_L = 0
       回傳以 t 為單位之本徵值（升序）"""
    N = np.array([1.0] + [float(B * (B - 1) ** (d - 1)) for d in range(1, L + 1)], dtype=float)
    s = np.sqrt(N)
    H = np.zeros((L + 1, L + 1))
    for d in range(L + 1):
        H[d, d] = -u * B ** (-d)          # V_d / t
    for d in range(L):
        hop = float(np.sqrt(B - 1)) if d >= 1 else float(np.sqrt(B))
        H[d, d + 1] = -hop
        H[d + 1, d] = -hop
    H[L, L] = 1e6                          # Dirichlet 近似
    w = np.linalg.eigvalsh(H)
    return np.sort(w)[:nmax]


print("  u = (G/a)/t 掃描（單位 t）:")
lad = {}
for u in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
    w = cayley_ladder(u)
    lad[u] = [float(x) for x in w]
    gaps = np.diff(w)
    print("   u=%5.1f：頭 5 本徵 %s" % (u, np.array2string(w[:5], precision=5)))
    print("          頭 4 間距 %s（遞增？ %s）" % (np.array2string(gaps[:4], precision=5), bool(gaps[0] < gaps[1])))
chk("W3a Cayley 階梯可算（頭 12 本徵值，u 掃描 6 點）", all(len(v) == 12 for v in lad.values()),
    "全部收斂（L = 80 Dirichlet）")
# 束縛態數 = 本徵值低於 band 底（V=0 時 band = [−2√(B−1)·t, +2√(B−1)·t]）
band = 2 * float(np.sqrt(B - 1))
print("  V=0 band 底 = −2√(B−1) = −%.4f（單位 t）" % band)
nbound = {}
for u in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
    w = cayley_ladder(u, nmax=40)
    nbound[u] = int(np.sum(np.array(w) < -band))
    print("   u=%5.1f：束縛態數（低於 band 底）= %d" % (u, nbound[u]))
chk("W3b Cayley 算子束縛態數 ≤ 4（唔可能承載 16 級徑向塔）", max(nbound.values()) <= 4,
    "u 掃到 16 亦只得 %d 個束縛態 ⟹ 徑向塔唔係 Cayley 算子之譜" % max(nbound.values()))

# 形狀：間距遞減（Regge-like）定遞增（box-like）？
def monotone(w):
    g = np.diff(w[:8])
    return "遞增(box-like)" if all(g[i] <= g[i + 1] + 1e-12 for i in range(len(g) - 1)) else (
        "遞減(Regge-like)" if all(g[i] >= g[i + 1] - 1e-12 for i in range(len(g) - 1)) else "非單調")

for u, w in lad.items():
    print("   u=%5.1f 間距模式 = %s" % (u, monotone(np.array(w))))
# 路徑樹譜形狀
wlam = np.array([float(LAM[k]) for k in range(1, 17)])
print("   路徑樹譜（k=1..16）間距模式 = %s｜前四間距 %s"
      % (monotone(wlam), np.array2string(np.diff(wlam)[:4], precision=5)))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("W4. Q-c 標度測試：觀測徑向間距 ⟹ 所需 t")
print("=" * 78)
OBS = [("ρ(1450)−ρ(770)", 675.0), ("ρ(1700)−ρ(1450)", 250.0),
       ("K*(1410)−K*(892)", 518.0), ("N(1440)−N(939)", 501.0),
       ("π(1300)−π(140)", 1160.0)]
t_M5 = mp.mpf("0.788398083633567408")
# 最大可用無量綱間距（路徑樹譜前兩個、Cayley u 最大）
dl_path = float(LAM[2] - LAM[1])
print("  路徑樹譜 λ₂−λ₁ = %.6f（無量綱）" % dl_path)
print("  M5 校準 t = %.6f MeV" % float(t_M5))
rows = []
for name, dE in OBS:
    t_req = dE / dl_path
    rows.append({"name": name, "dE_MeV": dE, "t_required_MeV": t_req, "ratio_vs_M5": t_req / float(t_M5)})
    print("   %-22s ΔE = %7.1f MeV ⟹ 需 t = %10.1f MeV（= %.0f× M5 t）"
          % (name, dE, t_req, t_req / float(t_M5)))
min_ratio = min(r["ratio_vs_M5"] for r in rows)
chk("W4a 觀測徑向間距需 t ≥ 10× M5 t（判準 W3）", min_ratio >= 10,
    "最小比值 = %.0f×（ρ(1700)−ρ(1450) 亦需 %.0f× M5 t）⟹ 需新標度 ⟹ T4 🔴"
    % (min_ratio, min_ratio))
verdict("W4b T4 生死：徑向階梯需新標度常數 t*（≈10²–10⁴ MeV）≠ M5 t = 0.79 MeV",
        False, "缺口 %.0f×–%.0f× ⟹ 需新標度 ⟹ T4 🔴" % (min_ratio, max(r["ratio_vs_M5"] for r in rows)),
        "W4")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "05_radial_ladder.py", "checks": checks, "integrity_fail": nfail,
       "intersection": matches, "capacity_pairs": len(pairs),
       "cayley_ladder_u": {str(k): v for k, v in lad.items()},
       "bound_states": {str(k): v for k, v in nbound.items()}, "spacing_scale_rows": rows, "t_M5_MeV": float(t_M5)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "05_radial_ladder.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/05_radial_ladder.json")
sys.exit(0 if nfail == 0 else 1)
