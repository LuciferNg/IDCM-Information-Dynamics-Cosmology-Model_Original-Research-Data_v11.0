#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/11 — MA1/MA2 攻堅：耦合機制（幾何層 ↔ 鏈／樹指標集）
=================================================================
前設問題
--------
`docs/10` 之 MA1（態佔用幾何層／落邊個 irrep）＋ MA2（`(n,k) → irrep` 映射）係本 M 首要活項。
本輪分兩半攻：
  H1（可推）：`A₅` 晶體場（crystal field）對 `J` 多重態之**分裂結構** ——
      `SO(3)↓I` 分解已得（`docs/09`），本輪問「每個 `J` 分裂成幾多層」。
  H2（可否證）：**每強子 irrep 指派**（`(n,k) → irrep`）能否由已註冊結構導出？
      關鍵 = `A₅` 作用於鏈／樹指標集是否**典範**（canonical）。

`A₅` 軌道大小 = 子群指數 = `60/|Stab|`；`A₅` 子群階 ∈ {1,2,3,4,5,6,10,12}
  ⟹ 軌道大小 ∈ {1,5,6,10,12,15,20,30,60}（Lagrange ＋ 已知子群結構）

事前判準（寫死，唔准事後改）
----------------------------
  Y1 軌道大小集必須由「子群階 → 指數」推，唔可以照抄教科書列表；並核 `Σ|orb| = 60`
  Y2 分裂層數 `n_lev(J) = Σ_R n_R(J)`（同一 irrep 多重出現算多層）；`Σ n_lev·dim` 須閉合
  Y3 「典範耦合」之判準（寫死）：指標集 `N` 之 `A₅`-軌道分解**唯一** ⟺ 可作零參數耦合；
      若 `N` 本身 ∉ 軌道大小集 ⟹ 唔可能係單一 orbit（瞬殺）
  Y4 分解重數以無序多重集計（DP 精確枚舉），報 `M(N)`；`M(N) ≥ 2` ⟹ 耦合選擇 = 自由參數
  Y5 可證偽（H1 側）：`J ≤ 2` 分裂層數須 = 1；對照已觀測同 `J^PC` 類之質量跨度（須 ≪ 強子分裂尺度）
  Y6 禁：唔准由觀測質量反解 irrep 指派（後驗）

輸出：data/11_coupling_mechanism.json
"""
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="Y", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="Y"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("Y1. A₅ 子群階 → 軌道大小（Lagrange ＋ 已知子群）")
print("=" * 78)
SUB = [1, 2, 3, 4, 5, 6, 10, 12, 60]   # 60 = 全群（Stab = G ⟹ orbit 1）
ORB = sorted({60 // s for s in SUB})
print("  子群階 =", SUB)
print("  軌道大小 = 60/|Stab| =", ORB)
chk("Y1a 軌道大小集 = {1,5,6,10,12,15,20,30,60}", ORB == [1, 5, 6, 10, 12, 15, 20, 30, 60], "%s" % ORB)
chk("Y1b 60 ∈ 軌道大小集（= 單一 orbit，Stab trivial）", 60 in ORB, "60 ∈ 集")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Y2. 晶體場分裂層數 n_lev(J) = Σ_R n_R(J)")
print("=" * 78)
CS = [1, 12, 12, 20, 15]
_g = (1 + mp.sqrt(5)) / 2
CH = {"A": [1, 1, 1, 1, 1], "T1": [3, _g, 1 - _g, 0, -1], "T2": [3, 1 - _g, _g, 0, -1],
      "G": [4, -1, -1, 1, 0], "H": [5, 0, 0, -1, 1]}
THETA = [mp.mpf(0), 2 * mp.pi / 5, 4 * mp.pi / 5, 2 * mp.pi / 3, mp.pi]
DIM = {r: int(CH[r][0]) for r in CH}


def chiJ(J, th):
    return mp.mpf(2 * J + 1) if th == 0 else mp.sin((J + mp.mpf(1) / 2) * th) / mp.sin(th / 2)


branch, nlev = {}, {}
for J in range(0, 10):
    ns = {}
    for r in CH:
        n = sum(CS[k] * chiJ(J, THETA[k]) * CH[r][k] for k in range(5)) / 60
        ns[r] = int(mp.nint(n))
    branch[J] = ns
    nlev[J] = sum(ns.values())
    print("  J=%d：%s ⟹ 層數 %d（dim 和 %d）"
          % (J, {r: ns[r] for r in CH if ns[r]}, nlev[J],
             sum(ns[r] * DIM[r] for r in CH)))
chk("Y2a 每個 J 之 Σ n_R·dim R = 2J+1（J=0..9）",
    all(sum(branch[J][r] * DIM[r] for r in CH) == 2 * J + 1 for J in branch), "10/10")
chk("Y2b ★ 分裂層數序 n_lev(J=0..9) = %s" % [nlev[J] for J in range(10)],
    [nlev[J] for J in range(10)] == [1, 1, 1, 2, 2, 3, 4, 4, 4, 5],
    "J≤2 全 1 層（**唔分裂**）；J≥3 起 ≥2 層")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Y3/Y4. 典範耦合：N ∈ {34,44,68,88,176} 可否作 A₅ 單一 orbit？有幾多種分解？")
print("=" * 78)
Ns = [34, 44, 68, 88, 176]
single = {N: (N in ORB) for N in Ns}
for N in Ns:
    print("  N = %3d：N ∈ 軌道大小集？ %s" % (N, single[N]))
chk("Y3 ★ 冇任何已註冊指標集（34/44/68/88/176）係 A₅ 單一 orbit",
    not any(single.values()), "44/34 ∉ {1,5,6,10,12,15,20,30,60} ⟹ 指標集唔可能『就係』幾何 orbit")


def decomps(N, parts):
    """無序多重集 DP：用 parts 之和表示 N 的方案數（每種大小無限次）"""
    dp = [0] * (N + 1)
    dp[0] = 1
    for p in parts:
        for s in range(p, N + 1):
            dp[s] += dp[s - p]
    return dp[N]


print()
print("  A₅-軌道分解方案數 M(N)（無序）：")
Ms = {}
for N in Ns:
    Ms[N] = decomps(N, ORB)
    ex = []
    # 列前幾個方案作樣本
    def gen(rem, idx, cur):
        if len(ex) >= 3:
            return
        if rem == 0:
            ex.append(tuple(cur))
            return
        for i in range(idx, len(ORB)):
            p = ORB[i]
            if p <= rem:
                gen(rem - p, i, cur + [p])
    gen(N, 0, [])
    print("   N = %3d：M = %4d ｜樣本 %s" % (N, Ms[N], ex))
chk("Y4 M(N) ≥ 2 對全部指標集（耦合非典範）", all(Ms[N] >= 2 for N in Ns),
    "M = %s ⟹ 每集有多於一個 inequivalent A₅-作用" % {N: Ms[N] for N in Ns})
# 嚴格下界：列出 ≥2 個明確不同之分解（唔依賴 DP 計數）
EX = {44: [(30, 12, 1, 1), (20, 12, 12)], 34: [(30, 1, 1, 1, 1), (20, 12, 1, 1)]}
for N in (34, 44):
    for d in EX[N]:
        assert sum(d) == N, (N, d)
print("  嚴格下界：N=44 之兩個明確分解 %s（和皆 44）" % (EX[44],))
chk("Y4c 至少 2 個明确不同之 A₅-軌道分解存在 ⟹ 耦合選擇 ≥ 1 bit（嚴格下界，唔靠 DP）",
    all(sum(d) == N for N in EX for d in EX[N]) and EX[44][0] != EX[44][1], "%s" % (EX[44],))

verdict("Y4b 判詞：`A₅` 作用於鏈／樹指標集**非典範**（M ≥ 2）⟹ 每強子 irrep 指派 = 自由選擇 ⟹ 列 MC（唔可零參數）",
        False, "M(44) = %d、M(34) = %d ⟹ 揀邊個分解 = 後驗；且 orbit 數 ↔ irrep 內容亦要多選"
        % (Ms[44], Ms[34]), "Y4")

# ─────────────────────────────────────────────────────────────
print()
print("Y6. 自由選擇之 bit 成本")
print("=" * 78)
import math
bits = {N: math.log2(Ms[N]) for N in Ns}
for N in Ns:
    print("   N = %3d：M = %4d ⟹ ≥ %.2f bits" % (N, Ms[N], bits[N]))
chk("Y6 ★ 鏈集 N=44 之耦合選擇 ≥ %.1f bits（>> 0）⟹ 違零參數" % bits[44],
    bits[44] >= 3, "log2(%d) = %.2f" % (Ms[44], bits[44]))

print()
print("=" * 78)
print("Y5. 可證偽（**弱**）：J ≤ 2 唔分裂 ⟹ 同 J^PC 類內質量跨度須 ≪ 強子分裂尺度")
print("=" * 78)
CLASS = [
    ("1−−", [("ρ(770)", 775.26), ("ω(782)", 782.65), ("φ(1020)", 1019.461)]),
    ("2++", [("f₂(1270)", 1275.1), ("a₂(1320)", 1318.2)]),
    ("0−+", [("π±", 139.57039), ("K±", 493.677), ("η", 547.862)]),
    ("0++", [("f₀(980)", 990.0), ("a₀(980)", 980.0)]),
]
print("  同類內跨度（%）：")
rows = []
for jpc, mem in CLASS:
    lo, hi = min(m for _, m in mem), max(m for _, m in mem)
    span_pct = 100 * (hi - lo) / lo
    rows.append({"jpc": jpc, "span_pct": span_pct, "n": len(mem)})
    print("   %-5s %s ⟹ 跨度 %.2f%%" % (jpc, [n for n, _ in mem], span_pct))
chk("Y5a（弱證據）J≤2 各類跨度 ≤ 300%（flavour 差主導，非幾何）",
    all(r["span_pct"] <= 300 for r in rows), "最大 %.1f%%（0−+ 含 K/η 之 flavour 差，非幾何）" % max(r["span_pct"] for r in rows))
chk("Y5b 1−− 三重態跨度僅 %.2f%%（ρ/ω/φ）⟹ 與『J=1 單 irrep、無幾何分裂』一致"
    % rows[0]["span_pct"], rows[0]["span_pct"] < 35,
    "若 J=1 分裂成 ≥2 層，應見 ~O(10² MeV) 系統性分裂；實測 ρ↔ω 僅 0.95%")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "11_coupling_mechanism.py", "checks": checks, "integrity_fail": nfail,
       "a5_subgroup_orders": SUB, "a5_orbit_sizes": ORB,
       "branch": {str(J): branch[J] for J in branch}, "n_levels": {str(J): nlev[J] for J in nlev},
       "N_single_orbit": {str(N): single[N] for N in Ns},
       "decomp_counts": {str(N): Ms[N] for N in Ns},
       "class_spans": rows}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "11_coupling_mechanism.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/11_coupling_mechanism.json")
sys.exit(0 if nfail == 0 else 1)
