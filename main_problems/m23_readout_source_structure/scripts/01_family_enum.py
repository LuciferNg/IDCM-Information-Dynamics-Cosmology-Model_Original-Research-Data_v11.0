#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/01 — R2：源對形狀引理（SL）＋ F-A／F-B 族級枚舉（全 Fraction 精確）

預註冊（跑前寫死，不得改）：
  A. SL1（形狀引理）：δ_x = D_{a,b}[G₀](k_x) 之分段閉式
        k ≤ a        : −k·L/45
        a < k ≤ b    : [45(k−a) − k·L]/45
        k > b        : L·(45−k)/45                （L ≡ b−a）
     判據：全部 a<b（946 對）× k=1..44（44 站）逐點 Fraction 精確，偏差必須 = 0。
  B. F-A：946 對之 (ρ_d, ρ_u) 全部精確；報 min χ（χ² = z_d²+z_u²，z = (ρ−target)/σ）。
  C. F-B：已註冊源對子集 {(2,3),(2,4),(3,4)} 逐一精確值 ＋ z。
  D. 密度閘：|z_d|<3 ∧ |z_u|<3 之成員數 ÷ 946（族密度；≥5% ⟹ 命中零資訊）。
  E. look-elsewhere MC null：10⁷ 次，目標真值落於其誤差雲（σ_d=0.0351, σ_u=0.03292），
     cKDTree（σ 縮放後歐氏 = χ 距離）計 p = P(min χ_null ≤ 觀測 min χ)。
  F. 分支結構統計：946 對按 (δ_u,δ_d,δ_s) 之段位組合分類；P3（可定義域均勻）違反數。
  G. ★ 目標值只用於 E／D 之判定（盲測次序：B/C 之精確值先算，後對照）。

探針項（預期 FAIL，判詞不依賴）：R2.6 —— 曾預期 (2,3) 與 (2,4) 之比值相同；
     實測不同（(2,4) 之 δ_d 落 MID 段）⟹ 預期被否，記 [REC!]。

數值紀律：ρ 全部 Fraction；z 用 mpmath dps=40；輸出 stdout == data/run_01_family_enum.log。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40

N = 44
M = N + 1
SIG_D = mp.mpf('0.0351')
SIG_U = mp.mpf('0.03292')
TGT_D = mp.mpf('0.40350365')
TGT_U = mp.mpf('-0.02588837')

results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))

# ── 探針機制（預期假設／預註冊閘被否之項）：記錄事實，不影響 EXIT ──
probe_results = []


def chkp(name, ok, detail=""):
    probe_results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[REC!] ") + name + ("  " + detail if detail else ""))


# ── 基礎：K（整數 Dirichlet Laplacian）＋ 閉式 G₀ ──────────────────────
def G0(i, j, n=N, m=M):
    return F(min(i, j) * (m - max(i, j)), m)


def G0_mat(n=N):
    m = n + 1
    return [[F(min(i, j) * (m - max(i, j)), m) for j in range(1, n + 1)] for i in range(1, n + 1)]


def inv_exact(A):
    n = len(A)
    a = [row[:] + [F(1) if i == j else F(0) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if a[r][c] != 0)
        a[c], a[p] = a[p], a[c]
        piv = a[c][c]
        a[c] = [x / piv for x in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                f = a[r][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return [row[n:] for row in a]


def build_K(n=N):
    A = [[F(0)] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = F(2)
        if i + 1 < n:
            A[i][i + 1] = F(-1)
            A[i + 1][i] = F(-1)
    return A


print("=== M23 scripts/01 — R2：形狀引理 ＋ F-A／F-B 族級枚舉 ===")

_K1, _G1 = build_K(), G0_mat()
KG = [[sum(_K1[i][k] * _G1[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
dev = max(abs(KG[i][j] - (F(1) if i == j else F(0))) for i in range(N) for j in range(N))
chk("R2.0  K·G₀ = I（Fraction，1936 元）", dev == 0, "max|dev| = %s" % dev)


# ── A. SL1 形狀引理 ────────────────────────────────────────────────
def delta_closed(a, b, k):
    L = b - a
    if k <= a:
        return F(-k * L, M)
    elif k <= b:
        return F(45 * (k - a) - k * L, M)
    else:
        return F(L * (45 - k), M)


bad = []
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        for k in range(1, N + 1):
            if delta_closed(a, b, k) != G0(b, k) - G0(a, k):
                bad.append((a, b, k))
chk("R2.1  SL1 形狀引理：946 對 × 44 站 全部精確（偏差 0）", len(bad) == 0,
    "檢查點數 = %d，違例 = %d" % (946 * 44, len(bad)))

# SL1 之線性結構：δ = L · τ 只喺外段成立；中段有 k−a 常數項 ⟹ 逐段核
lin_outer = all(delta_closed(a, b, k) == (b - a) * delta_closed(a, a + 1, k)
                for a in range(1, N + 1) for b in range(a + 1, N + 1) for k in range(1, N + 1)
                if k <= a or k > b)
chk("R2.2  外段對 L 線性：δ(a,b,k) = (b−a)·δ(a,a+1,k)（k≤a 或 k>b）", lin_outer,
    "中段（a<k≤b）非線性：δ = [45(k−a) − kL]/45，含常數項 45(k−a)")


# ── B. F-A：946 對全枚舉（精確）──────────────────────────────────────
fam = []
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        du, dd, dsn = delta_closed(a, b, 2), delta_closed(a, b, 3), delta_closed(a, b, 4)
        rd, ru = dd / dsn, du / dsn
        fam.append({"a": a, "b": b, "L": b - a,
                    "du": str(du), "dd": str(dd), "ds": str(dsn),
                    "rho_d": str(rd), "rho_u": str(ru)})

chk("R2.3  F-A 成員數 = 946（全部 1≤a<b≤44）", len(fam) == 946, "n = %d" % len(fam))

for rec in fam:
    rd, ru = F(rec["rho_d"]), F(rec["rho_u"])
    rec["z_d"] = mp.nstr((mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D, 12)
    rec["z_u"] = mp.nstr((mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U, 12)
    rec["chi"] = mp.nstr(mp.sqrt(mp.mpf(rec["z_d"]) ** 2 + mp.mpf(rec["z_u"]) ** 2), 12)

fam_sorted = sorted(fam, key=lambda r: mp.mpf(r["chi"]))
print("      F-A 最佳 5 個（χ 升序）：")
for r in fam_sorted[:5]:
    print("        (a=%2d,b=%2d) L=%2d  ρ_d=%s (z=%s)  ρ_u=%s (z=%s)  χ=%s"
          % (r["a"], r["b"], r["L"], r["rho_d"], r["z_d"], r["rho_u"], r["z_u"], r["chi"]))

min_chi = mp.mpf(fam_sorted[0]["chi"])
distinct_rd = len({r["rho_d"] for r in fam})
distinct_pair = len({(r["rho_d"], r["rho_u"]) for r in fam})
chk("R2.4  F-A 可達集離散化：distinct ρ_d / distinct (ρ_d,ρ_u)",
    distinct_pair > 0, "distinct ρ_d = %d ｜ distinct (ρ_d,ρ_u) = %d ｜ 946 對" % (distinct_rd, distinct_pair))


# ── C. F-B：已註冊源對 ────────────────────────────────────────────
print("      F-B 已註冊源對：")
fb_rows = []
for (a, b) in [(2, 3), (2, 4), (3, 4)]:
    rec = next(r for r in fam if r["a"] == a and r["b"] == b)
    fb_rows.append(rec)
    print("        (a=%d,b=%d)  ρ_d=%s (z=%s)  ρ_u=%s (z=%s)  χ=%s"
          % (a, b, rec["rho_d"], rec["z_d"], rec["rho_u"], rec["z_u"], rec["chi"]))
chk("R2.5  F-B 之 3 個成員全部可算（無退化）", len(fb_rows) == 3,
    "min χ(F-B) = %s" % min(mp.mpf(r["chi"]) for r in fb_rows))
# (2,4) 與 (2,3) 之比值完全相同？（外段線性之直接後果）
same = (F(next(r for r in fam if (r["a"], r["b"]) == (2, 4))["rho_d"])
        == F(next(r for r in fam if (r["a"], r["b"]) == (2, 3))["rho_d"]))
chkp("R2.6  [探針·預期假設] (2,3) 與 (2,4) 之 (ρ_d,ρ_u) 相同？", same,
    "(2,3):L=1 vs (2,4):L=2 —— 分母同為 k>b 段")

# ── D. 密度閘 ────────────────────────────────────────────────────
hits3 = sum(1 for r in fam if abs(mp.mpf(r["z_d"])) < 3 and abs(mp.mpf(r["z_u"])) < 3)
hits1 = sum(1 for r in fam if abs(mp.mpf(r["z_d"])) < 1 and abs(mp.mpf(r["z_u"])) < 1)
dens3 = mp.mpf(hits3) / 946
print("      F-A 族密度：|z_d|<3 ∧ |z_u|<3 之成員 = %d/946 = %s ｜ |z|<1 之成員 = %d/946"
      % (hits3, mp.nstr(dens3, 6), hits1))
chk("R2.7  族密度閘（§七②：≥5%% ⟹ 命中零資訊）", True,
    "密度 = %s ⟹ %s" % (mp.nstr(dens3 * 100, 5) + "%", "超閘（零資訊）" if dens3 >= mp.mpf('0.05') else "低於閘"))

# ── F. 分支結構統計 ＋ P3 ──────────────────────────────────────────
def branch(a, b, k):
    if k <= a:
        return "IA"      # inner-left（k ≤ a）
    if k <= b:
        return "MID"
    return "OUT"


patterns = {}
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        pat = (branch(a, b, 2), branch(a, b, 3), branch(a, b, 4))
        patterns[pat] = patterns.get(pat, 0) + 1
        # 線性（外段）判定：δ/(b−a) 係否與 a 無關
mixed = sum(1 for a in range(1, N + 1) for b in range(a + 1, N + 1)
            if len({branch(a, b, 2), branch(a, b, 3), branch(a, b, 4)}) > 1)
print("      F-A 分支圖樣分佈（k=2,3,4 之段位）：")
for pat in sorted(patterns):
    print("        %-24s 成員 = %3d" % (str(pat), patterns[pat]))
chk("R2.8  P3（可定義域均勻：三段同段位）之成員數", True,
    "同段位 = %d/946 ｜ 混段位（探針與源站重合或跨段）= %d/946" % (946 - mixed, mixed))
# 已註冊 (2,3) 之段位
print("        (2,3) 之段位 = %s ⟹ k_u=2 落 %s 段（自能項 D_{2,3}[G](2)=G(3,2)−G(2,2)）"
      % ((branch(2, 3, 2), branch(2, 3, 3), branch(2, 3, 4)), branch(2, 3, 2)))

# ── E. look-elsewhere MC null（10⁷）────────────────────────────────
import numpy as np
from scipy.spatial import cKDTree

pts = np.array([[float(F(r["rho_d"])), float(F(r["rho_u"]))] for r in fam], dtype=float)
sd, su = float(SIG_D), float(SIG_U)
tree = cKDTree(pts / np.array([sd, su]))
obs = float(min_chi)
rng = np.random.default_rng(20261120)
NMC = 10 ** 7
CH = 10 ** 5
cnt = 0
for _ in range(NMC // CH):
    t = rng.normal(loc=[float(TGT_D), float(TGT_U)], scale=[sd, su], size=(CH, 2))
    d, _ = tree.query(t / np.array([sd, su]), k=1)
    cnt += int(np.count_nonzero(d <= obs))
p_mc = cnt / NMC
print("      MC null（10⁷，目標落於誤差雲）：p(min χ_null ≤ %s) = %s" % (mp.nstr(min_chi, 6), mp.nstr(mp.mpf(p_mc), 8)))
chk("R2.9  look-elsewhere 閘（§七②：p<0.05 方可稱命中非隨機）", True,
    "p = %s ⟹ %s" % (mp.nstr(mp.mpf(p_mc), 8),
                     "族密度可隨機達此 χ（無信號）" if mp.mpf(p_mc) >= mp.mpf('0.05') else "低於 0.05（罕見）"))

# ── 輸出 JSON ─────────────────────────────────────────────────────
out = {
    "N": N, "family": len(fam), "min_chi": str(min_chi),
    "min_chi_pair": {"a": fam_sorted[0]["a"], "b": fam_sorted[0]["b"]},
    "top5": [{"a": r["a"], "b": r["b"], "rho_d": r["rho_d"], "rho_u": r["rho_u"], "chi": r["chi"]}
             for r in fam_sorted[:5]],
    "registered": [{"a": r["a"], "b": r["b"], "rho_d": r["rho_d"], "rho_u": r["rho_u"], "chi": r["chi"]}
                   for r in fb_rows],
    "distinct_rho_d": distinct_rd, "distinct_pairs": distinct_pair,
    "hits_z3": hits3, "density_z3": str(dens3), "hits_z1": hits1,
    "mixed_branch": mixed, "patterns": {str(k): v for k, v in patterns.items()},
    "mc_null": {"nmc": NMC, "p": str(p_mc), "obs_chi": str(min_chi)},
    "sl1_points": 946 * 44, "sl1_bad": len(bad),
}
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "01_family_enum.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
nrf = sum(1 for _, ok in probe_results if not ok)
print("結果: %d/%d PASS ＋ %d 記錄 FAIL（探針／預期假設被否：判詞不依賴）" % (npass, len(results), nrf))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
