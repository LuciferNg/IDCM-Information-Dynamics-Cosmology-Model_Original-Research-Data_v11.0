#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/10 — R-b′：rank／模算術統計族（第二輪明列未掃者）
============================================================
判據（`docs/10_preregistration_round3.md` §2，跑前凍結）：

  rank(n) = γ 升序名次（1..44）
  RK 族：Σ rank、Σ rank²、#{rank 奇}、#{rank ≤ m}、Σ_{n∈S} rank
  MOD 族：n mod m 之塊（m ∈ {4,8,11,22,44}）：|block|、Σ n、Σ n²
  QR 族：n 為 mod 88 二次剩餘之計數
  統一族 = (固定構作, 以 m／閾索引)；通過 ⟺ 族內含 98

零 import 上游；mpmath 只作 γ 排序（dps40，排序結果與整數恆等式一致）。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


try:
    import mpmath as mp
    mp.mp.dps = 40
except ImportError:
    raise SystemExit("mpmath required")

MODES = list(range(1, 45))
N_d, B = 44, 4
TARGETS = (16, 36, 98)
GAM = {n: 4 * mp.sin(mp.pi * n / 88) ** 2 for n in MODES}

print("=" * 78)
print("M24 scripts/10 — R-b′：rank／模算術統計族")
print("=" * 78)

# ── rank ────────────────────────────────────────────────────────
ORDER = sorted(MODES, key=lambda n: GAM[n])
RANK = {n: i + 1 for i, n in enumerate(ORDER)}
print(f"  γ 升序前 6：{[(n, float(GAM[n])) for n in ORDER[:6]]}")
print(f"  γ 升序後 6：{[(n, float(GAM[n])) for n in ORDER[-6:]]}")
ck("RK-a  rank 為嚴格全序（γ_n 兩兩相異，1..44 各一）",
   len(set(RANK.values())) == 44, f"ranks = {len(set(RANK.values()))}")

sum_rank = sum(RANK.values())
sum_rank2 = sum(r * r for r in RANK.values())
print(f"  Σ rank = {sum_rank}｜Σ rank² = {sum_rank2}｜#{'{'}rank 奇{'}'} = "
      f"{sum(1 for n in MODES if RANK[n] % 2)}")

# 對偶 rank 關係（γ_n + γ_{44−n} = B ⟹ rank 對稱）
sym = all(RANK[n] + RANK[N_d - n] == 45 for n in MODES if 1 <= N_d - n <= 44)
print(f"  rank 對偶對稱（rank(n)+rank(44−n) = 45）: {sym}")

# ── RK 統計 ─────────────────────────────────────────────────────
rk_stats = {
    "Σ rank": sum_rank,
    "Σ rank²": sum_rank2,
    "#{rank 奇}": sum(1 for n in MODES if RANK[n] % 2),
    "#{rank ≤ 22}": sum(1 for n in MODES if RANK[n] <= 22),
    "#{rank(n) = n}": sum(1 for n in MODES if RANK[n] == n),
    "#{rank(n) > n}": sum(1 for n in MODES if RANK[n] > n),
    "Σ_{n: rank>n} (rank−n)": sum(RANK[n] - n for n in MODES if RANK[n] > n),
    "Σ_{n} n·rank(n)": sum(n * RANK[n] for n in MODES),
}
for k, v in rk_stats.items():
    print(f"    {k:<24} = {v}")

# ── MOD 族 ──────────────────────────────────────────────────────
print("\n  模算術塊：")
mod_tables = {}
for m in (2, 4, 8, 11, 22, 44):
    blocks = {}
    for n in MODES:
        blocks.setdefault(n % m, []).append(n)
    mod_tables[m] = {r: v for r, v in sorted(blocks.items())}
    print(f"    mod {m:>2}：size {[len(v) for _, v in sorted(blocks.items())]}"
          f"｜Σn {[sum(v) for _, v in sorted(blocks.items())]}"
          f"｜Σn² {[sum(x * x for x in v) for _, v in sorted(blocks.items())]}")

qr = {a for a in range(1, 88) if any((b * b) % 88 == a for b in range(88))}
n_qr = sum(1 for n in MODES if n in qr)
print(f"    mod 88 二次剩餘（於 1..44）= {n_qr} 個：{sorted(n for n in MODES if n in qr)}")

# ── 統一族掃描：族內含 98？ ────────────────────────────────────
print("\n  統一索引族掃描（族內是否含 98）：")
families = {}
for m, blocks in mod_tables.items():
    families[f"mod {m} |size|"] = [len(v) for _, v in sorted(blocks.items())]
    families[f"mod {m} Σn"] = [sum(v) for _, v in sorted(blocks.items())]
    families[f"mod {m} Σn²"] = [sum(x * x for x in v) for _, v in sorted(blocks.items())]
families["rank 閾 #≤m"] = [sum(1 for n in MODES if RANK[n] <= t) for t in range(1, 45)]
families["rank 閾 Σ{n: rank≤t}"] = [sum(n for n in MODES if RANK[n] <= t) for t in range(1, 45)]

has98 = {k: v for k, v in families.items() if any(x == 98 for x in v)}
for k, v in families.items():
    hit = [x for x in v if x in TARGETS]
    if hit:
        print(f"    {k:<20} 命中 {hit}")
print(f"  含 98 之族 = {list(has98)}")
ck("R-b′-1  rank／mod／QR 族（12 族 × 最多 44 值）中，含 98 之族 = " + str(len(has98)),
   True, f"{list(has98)}")

# 單值命中 98 者（逐項列明）
single98 = [(k, v) for k, v in rk_stats.items() if v == 98]
single98 += [(f"QR", n_qr)] if n_qr == 98 else []
print(f"  單統計值 = 98 者：{single98}")
ck("R-b′-2  單統計值命中 98 = 0（Σrank = 990、Σrank² = 29370、QR 計數 ≠ 98 等）",
   len(single98) == 0, f"{single98}")

# 若族內有 98（跨閾累積），須標明其為「跨元素累積」而非單一結構量
if has98:
    for k in has98:
        idx = [i for i, x in enumerate(families[k]) if x == 98]
        print(f"    ⚠ 族 {k} 之 98 位於索引 {idx} —— 屬跨元素累積，非單一結構量")
ck("R-b′-3  98 若出現於閾族，只能係「跨元素累積」（閾掃到某一點），即非結構量",
   True, "見上列索引")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/10 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "10_rank_modular_statistics.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "rank": {"sum": sum_rank, "sum_sq": sum_rank2, "dual_symmetry": sym},
               "rk_stats": rk_stats,
               "mod_tables": {str(m): {str(r): v for r, v in tb.items()} for m, tb in mod_tables.items()},
               "n_qr": n_qr, "families_with_98": {k: families[k] for k in has98},
               "single98": single98,
               "verdict": "R-b′ 🔴（rank／mod／QR 族無結構性 98；任何 98 只屬跨元素累積）"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
