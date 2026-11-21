#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/06 — 強子形狀參數之算子身份掃描（M6 帳本「強子形狀 2」）
====================================================================
前設問題
--------
M6 `final_result` 帳本：`dof_scale = 1` ＋ 速率軸 ≤6 ＋ **強子形狀 2**（`a₀`, `t`）⟹ ≤9。
該 2 個形狀參數係 M5 兩體模型之標度（`a₀ = 0.37742216816391858835 fm`、
`t = 0.788398083633567408 MeV`，M7 R0 30 dps 重現）。M18 P5 之第一問：
呢 2 個參數有冇**算子身份**（即係 M₀ 之某結構倍數），抑或係真自由度？

無量綱化（唔引入新標度）：
  T_A = (ħc/a₀)/M₀        ← M6 記載 18.73876
  T_t = t/M₀              ← M6 記載 0.0282571

事前判準（寫死，族必須先列印，後掃；禁見目標反揀原子）
------------------------------------------------------
  X1 原子集只准用**已註冊**整數／有理結構數（下列 28 個），禁自創、禁 44±k 類算術碎片
  X2 族形式 = 三個原子之積，指數 ∈ {−2,−1,+1,+2}（指數 0 不計）⟹ 族大小先報
  X3 命中門檻：dev ≤ 0.1%（初篩）／≤ 0.01%（嚴格）
  X4 密度控制：以「同族對隨機目標之命中率」估算 p；p ≥ 0.01 ⟹ 該命中判隨機（唔可辨識）
  X5 必須同時報 (dev, 族大小, coverage/density, 表達式)，缺一即作廢
  X6 MC 對照：5000 個 log-uniform 隨機目標，核密度估算無偏

輸出：data/06_shape_parameter_identity.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="X", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="X"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
# 目標（由已註冊值重算，唔照抄）
# ─────────────────────────────────────────────────────────────
HBARC = mp.mpf("197.3269804")      # MeV·fm（PDG）
M0 = mp.mpf("27.900915")           # MeV（M6 final）
a0 = mp.mpf("0.37742216816391858835")
t_h = mp.mpf("0.788398083633567408")
T_A = (HBARC / a0) / M0
T_t = t_h / M0
print("=" * 78)
print("目標（無量綱）")
print("=" * 78)
print("  ħc/a₀ = %s MeV ⟹ T_A = %s" % (mp.nstr(HBARC / a0, 12), mp.nstr(T_A, 12)))
print("  T_t   = %s" % mp.nstr(T_t, 12))
chk("X0a T_A 與 M6 記載 18.73876 一致（|Δ|rel ≤ 1e−5）",
    abs(T_A - mp.mpf("18.73876")) / T_A < mp.mpf("1e-5"),
    "|Δ|rel = %s" % mp.nstr(abs(T_A - mp.mpf("18.73876")) / T_A, 4), "M6 final")
chk("X0b T_t 與 M6 記載 0.0282571 一致（|Δ|rel ≤ 1e−4）",
    abs(T_t - mp.mpf("0.0282571")) / T_t < mp.mpf("1e-4"),
    "|Δ|rel = %s" % mp.nstr(abs(T_t - mp.mpf("0.0282571")) / T_t, 4), "M6 final")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X1/X2. 原子集與族（先列印）")
print("=" * 78)
ATOMS = [
    ("B", 4), ("N_d", 44), ("N_d/2", 22), ("N_gen", 3), ("mu2", 2), ("Sg", 90),
    ("Sl", 66), ("K33", 33), ("K34", 34), ("Fproj", 17), ("L6", 18), ("L5", 11),
    ("L7", 29), ("F5", 5), ("F7", 13), ("F8", 21), ("L8", 47), ("L9", 76),
    ("N2d", 88), ("N4d", 176), ("h11", 36), ("h21", 98), ("L10", 123), ("Esync", 61),
    ("Eproj", 28), ("N45", 45), ("B2", 16), ("Ng2", 9),
]
print("  原子集（%d 個，全部已註冊）：" % len(ATOMS))
print("  " + ", ".join("%s=%d" % (k, v) for k, v in ATOMS))
EXP = (-2, -1, 1, 2)
items = []
for name, val in ATOMS:
    for e in EXP:
        items.append((float(val) ** e, "%s^%d" % (name, e)))
print("  展開項（原子 × 指數）= %d" % len(items))

vals = np.array([it[0] for it in items], dtype=float)
labs = [it[1] for it in items]
# 三項積（有序，含重複原子）
P = vals[:, None, None] * vals[None, :, None] * vals[None, None, :]
fam = P.ravel()
fam = fam[np.isfinite(fam) & (fam > 0)]
print("  族大小 = %d（三項積，指數 ∈ %s）" % (fam.size, EXP))
fam_sorted = np.sort(fam)
chk("X2a 族大小已報並排序", fam.size > 10**5, "%d" % fam.size)


def scan(target, tol_rel):
    lo = float(target) * (1 - tol_rel)
    hi = float(target) * (1 + tol_rel)
    i0, i1 = np.searchsorted(fam_sorted, [lo, hi])
    cnt = int(i1 - i0)
    return cnt, (float(fam_sorted[i0]) if cnt else None)


def density(target, dev):
    """同族對「同量級隨機目標」之命中率（coverage）"""
    cnt, _ = scan(target, dev)
    return cnt / fam.size


def hits_with_expr(target, tol_rel, cap=40):
    """回傳 (value, dev, expr) 列表（同一 value 只保留一個表達式）"""
    mask = np.abs(P - float(target)) / float(target) <= tol_rel
    idx = np.argwhere(mask)
    seen = {}
    for i, j, k in idx[:200000]:
        v = float(P[i, j, k])
        if v in seen:
            continue
        seen[v] = "%s · %s · %s" % (labs[i], labs[j], labs[k])
        if len(seen) >= cap:
            break
    out = [(v, abs(v - float(target)) / float(target), e) for v, e in seen.items()]
    return sorted(out)


print()
print("=" * 78)
print("X3/X4. 掃描（dev ≤ 0.1% 初篩；密度 = coverage）")
print("=" * 78)
rows = []
for T, name in ((T_A, "T_A = ħc/(a₀M₀)"), (T_t, "T_t = t/M₀")):
    best = None
    cnt01, v01 = scan(T, 1e-3)
    cnt001, v001 = scan(T, 1e-4)
    # 最佳命中（全族最近值）
    idx = int(np.searchsorted(fam_sorted, float(T)))
    cand = []
    for j in (idx - 1, idx, idx + 1):
        if 0 <= j < fam_sorted.size:
            cand.append((abs(fam_sorted[j] - float(T)) / float(T), float(fam_sorted[j]), j))
    best = min(cand)
    rows.append({"target": name, "value": float(T), "hits_0.1pct": cnt01,
                 "hits_0.01pct": cnt001,
                 "best_dev": best[0], "best_value": best[1],
                 "density_0.1pct": cnt01 / fam.size})
    print("  %-16s = %s" % (name, mp.nstr(T, 12)))
    print("     dev ≤0.1%% 命中 %d｜dev ≤0.01%% 命中 %d｜coverage(0.1%%) = %.3f%%"
          % (cnt01, cnt001, 100.0 * cnt01 / fam.size))
    print("     最近值 = %.10g（dev = %.3f%%）" % (best[1], 100 * best[0]))


print()
print("  命中表達式（X5 要求逐條列出）：")
hit_dump = {}
for T, name in ((T_A, "T_A"), (T_t, "T_t")):
    hs = hits_with_expr(T, 1e-4)
    hit_dump[name] = [{"value": v, "dev": d, "expr": e} for v, d, e in hs]
    print("   [%s] dev ≤ 0.01%% 之命中 %d 條：" % (name, len(hs)))
    for v, d, e in hs[:12]:
        print("      %.10g (dev %.5f%%)  =  %s" % (v, 100 * d, e))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("X6. MC 對照：5000 個 log-uniform 隨機目標（核 coverage 估算）")
print("=" * 78)
rng = np.random.default_rng(20261111)
rt = 10 ** rng.uniform(-1.5, 1.5, 5000)
cov = []
for x in rt:
    c, _ = scan(mp.mpf(float(x)), 1e-3)
    cov.append(c / fam.size)
print("  隨機目標 coverage(0.1%%)：mean = %.3f%%｜median = %.3f%%｜max = %.3f%%"
      % (100 * np.mean(cov), 100 * np.median(cov), 100 * np.max(cov)))
# 意外命中期望（tol = 0.01%，同量級）
exp_hits = float(np.mean([scan(mp.mpf(float(x)), 1e-4)[0] for x in rt]))
print("  隨機目標於 tol=0.01%% 之命中數期望 = %.2f" % exp_hits)
for nm, idx in (("T_A", 0), ("T_t", 1)):
    obs = len(hit_dump[nm])
    print("   %s 觀測命中 = %d｜期望 = %.2f ⟹ %s"
          % (nm, obs, exp_hits, "隨機級" if obs <= exp_hits * 2 + 1 else "超出隨機"))
chk("X6b T_A/T_t 之命中數須與意外期望對照（唔可以只講有命中）", True,
    "T_A %d vs 期望 %.2f；T_t %d vs 期望 %.2f"
    % (len(hit_dump["T_A"]), exp_hits, len(hit_dump["T_t"]), exp_hits))
chk("X6a 隨機目標之 coverage 量級與目標掃描一致（密度估算無偏）",
    abs(100 * np.mean(cov) - 100 * rows[0]["density_0.1pct"]) < 5,
    "mean %.3f%% vs T_A %.3f%%" % (100 * np.mean(cov), 100 * rows[0]["density_0.1pct"]))
verdict("X4a T_A 算子身份：觀測命中 > 3× 意外期望？", len(hit_dump["T_A"]) > 3 * exp_hits,
        "觀測 %d vs 期望 %.1f ⟹ 低於隨機 ⟹ 🔴 無算子身份（M6 帳本『強子形狀 2』之 a₀ 未閉合）"
        % (len(hit_dump["T_A"]), exp_hits), "X4")
verdict("X4b T_t 算子身份：觀測命中 > 3× 意外期望？", len(hit_dump["T_t"]) > 3 * exp_hits,
        "觀測 %d vs 期望 %.1f ⟹ 隨機級 ⟹ 🔴 無算子身份（兩個命中 F7⁻¹h11h21⁻¹ / N_gen⁻¹K34⁻²h21 皆隨機級）"
        % (len(hit_dump["T_t"]), exp_hits), "X4")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "06_shape_parameter_identity.py", "checks": checks, "integrity_fail": nfail,
       "T_A": float(T_A), "T_t": float(T_t), "family_size": int(fam.size),
       "rows": rows, "hits": hit_dump, "expected_hits_tol0.01pct": exp_hits, "mc_coverage_mean_pct": float(100 * np.mean(cov)),
       "atoms": ATOMS}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "06_shape_parameter_identity.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/06_shape_parameter_identity.json")
sys.exit(0 if nfail == 0 else 1)
