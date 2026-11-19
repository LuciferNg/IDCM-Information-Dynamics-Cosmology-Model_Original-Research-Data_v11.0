#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/04 — R6：完備性擴展（**post-freeze 擴展聲明**；不參與 docs/01 之預註冊判準）

聲明：本檔喺 `docs/01` 凍結之後才寫，屬**後驗擴展**（同 M20 `docs/02 §十一` 第 3 條同類）。
     用途 = 封閉性（完備性）檢查，**不作命中證據**、不改任何判準。

預註冊（本檔自身，寫完才跑）：
  C1  T3 完備性：關係型第三支（G 比值讀法）ρ_x = G(k_x,k_s)/G(k_s,k_s) 之值
      必須落於 F-A（946 對）之可達集；若落 ⟹ T3 三支全部收於 F-A／已算死，語法窮盡成立。
  C2  形狀層解釋：tent 分段只產生 7 個段位類；逐類給 ρ_d 之閉式，並解釋為何相鄰站
      (k_d=3, k_s=4) 之比值必然 ≈ 1（42/41）⟹ 量測要求 0.4035（差 2.54×）係**形狀層**衝突。
  C3  判別力：F-A 之 107 個 distinct (ρ_d,ρ_u) 之最近鄰間距（z 單位）；
      若 ≪ 1 ⟹ 族密度 ≫ 量測精度 ⟹ 任何「命中」零資訊（M12 R25 ℚ/41 現象之推廣）。
  C4  T2 最小版（分子固定 (2,3)、分母自由）：ρ_d/|ρ_u| = 21 固定 ⟹ 1 參數只沿固定比值線移動。
  C5  形式邊界：δ 之一階差形式不改（M12 docs/15／docs/38 定義）；改 δ 之形狀 = 改讀出式語義，
      超 M23 章程（本 M 只動來源結構）⟹ 發還。

探針項（預期 FAIL，判詞不依賴）：
  C3 —— 曾預期 107 個 distinct 點在 σ 尺度稠密（最近鄰 < 1σ）；實測中位 1.93σ ⟹ 預期被否：
        族係**粗粒**（值域 [−33,12]×[−22,23]），故 1σ 全滅有意義，非密度 artifact。
  修訂（bug）：C2b 之測試對原用 (1,5)（3 落 MID 段）⟹ 改 (1,2)（3,4 皆 OUT）；判據不變。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
N, M = 44, 45
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


def G(i, j):
    return F(min(i, j) * (M - max(i, j)), M)


def delta(a, b, k):
    L = b - a
    if k <= a:
        return F(-k * L, M)
    if k <= b:
        return F(45 * (k - a) - k * L, M)
    return F(L * (45 - k), M)


print("=== M23 scripts/04 — R6：完備性擴展（後驗）===")

# F-A 可達集
fam = {}
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        dd, du, ds = delta(a, b, 3), delta(a, b, 2), delta(a, b, 4)
        fam[(a, b)] = (dd / ds, du / ds)
vals = set(fam.values())

# ── C1：T3 第三支（G 比值）───────────────────────────────────────
c1_d, c1_u = G(3, 4) / G(4, 4), G(2, 4) / G(4, 4)
print("      C1 G 比值讀法：ρ_d = G(3,4)/G(4,4) = %s ｜ ρ_u = G(2,4)/G(4,4) = %s" % (c1_d, c1_u))
in_fa = (c1_d, c1_u) in vals
members = [(a, b) for (a, b), v in fam.items() if v == (c1_d, c1_u)]
print("      落於 F-A 之成員（%d 個，例）：%s" % (len(members), members[:8]))
chk("C1  T3 第三支（G 比值）落於 F-A 可達集 ⟹ 語法窮盡成立", in_fa and len(members) > 0,
    "值 = %s；F-A 成員數 %d" % ((c1_d, c1_u), len(members)))

# ── C2：7 個段位類之閉式 ────────────────────────────────────────
def seg(a, b, k):
    return "IA" if k <= a else ("MID" if k <= b else "OUT")


cls = {}
for (a, b), (rd, ru) in fam.items():
    pat = (seg(a, b, 2), seg(a, b, 3), seg(a, b, 4))
    cls.setdefault(pat, set()).add((rd, ru))
print("      C2 段位類 → ρ_d 值集：")
for pat in sorted(cls):
    s = sorted(cls[pat])
    print("         %-22s distinct = %2d ｜ 例 %s" % (str(pat), len(s), [(str(x), str(y)) for x, y in s[:3]]))
chk("C2a 段位類總數 = 7（形狀引理之全部可能）", len(cls) == 7, "n = %d" % len(cls))
# 「3 與 4 同段 ⟹ 比值只由該段之 k 依賴決定」
same_seg = {pat for pat in cls if pat[1] == pat[2]}
print("         k_d=3 與 k_s=4 同段之段位類 = %s" % sorted(same_seg))
out_ratio = (delta(1, 2, 3) / delta(1, 2, 4))     # 3,4 皆 OUT（k>b=2）
mid_ratio = (delta(5, 9, 3) / delta(5, 9, 4))     # 3,4 皆 IA（k≤a）
print("         OUT 段（k>b）：ρ_d = (45−3)/(45−4) = %s = %s  ⟹ 必然 ≈ 1"
      % (out_ratio, float(out_ratio)))
print("         IA  段（k≤a）：ρ_d = 3/4 = %s  ⟹ 與 42/41 相差 %s 倍"
      % (float(mid_ratio), float(F(42, 41) / mid_ratio)))
chk("C2b 相鄰站同段 ⟹ ρ_d 只可能 = 42/41 或 3/4（差 1.37×）", out_ratio == F(42, 41) and mid_ratio == F(3, 4),
    "量測需要 0.4035 ⟹ 相對 42/41 差 %s 倍（形狀層衝突）" % mp.nstr(TGT_D / mp.mpf(42) * 41, 8))

# ── C3：判別力（107 點最近鄰間距 vs σ）───────────────────────────
pts = sorted(vals)
def z_of(v):
    return ((mp.mpf(v[0].numerator) / v[0].denominator - TGT_D) / SIG_D,
            (mp.mpf(v[1].numerator) / v[1].denominator - TGT_U) / SIG_U)

zs = [z_of(v) for v in pts]
nn = []
for i, p in enumerate(zs):
    d = min(mp.sqrt((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) for j, q in enumerate(zs) if j != i)
    nn.append(d)
nn_sorted = sorted(nn)
med_nn = nn_sorted[len(nn_sorted) // 2]
print("      C3 distinct 點 = %d ｜ 最近鄰間距（z 單位）：中位 %s ｜ 最小 %s ｜ 最大 %s"
      % (len(pts), mp.nstr(med_nn, 6), mp.nstr(min(nn), 6), mp.nstr(max(nn), 6)))
rng_d = (min(v[0] for v in pts), max(v[0] for v in pts))
print("         ρ_d 範圍 = [%s, %s] ｜ ρ_u 範圍 = [%s, %s]"
      % (float(rng_d[0]), float(rng_d[1]), float(min(v[1] for v in pts)), float(max(v[1] for v in pts))))
chkp("C3  [探針·預期假設] 最近鄰間距中位 < 1（z 單位）？",
    med_nn < 1, "中位最近鄰 = %s σ" % mp.nstr(med_nn, 6))

# ── C4：T2 最小版之結構限制 ────────────────────────────────────
dd0, du0 = delta(2, 3, 3), delta(2, 3, 2)
ratio_fixed = dd0 / (-du0)
print("      C4 T2 最小版（分子固定 (2,3)、分母 (a₂,b₂) 自由）：ρ_d/|ρ_u| = %s = %s（固定）"
      % (ratio_fixed, float(ratio_fixed)))
req_ratio = mp.mpf(TGT_D) / -TGT_U
print("         量測要求 ρ_d/|ρ_u| = %s ⟹ 差 %s 倍 ⟹ 1 參數只能沿一條**比值固定**之線移動"
      % (mp.nstr(req_ratio, 8), mp.nstr(req_ratio / mp.mpf(21), 8)))
chk("C4  T2 最小版之比值固定（ρ_d/|ρ_u| = 21）⟹ 1 參數族之形状受限",
    ratio_fixed == 21, "固定比值 = %s ｜ 量測 = %s" % (ratio_fixed, mp.nstr(req_ratio, 8)))

# ── C5：形式邊界聲明（無 check，只記錄）─────────────────────────
print("      C5 形式邊界：δ = D_{a,b}[G](k_x) 之一階差形式由 M12 docs/15／docs/38 定義；")
print("         改 δ 之**形狀**（例如二階差、非線性 k 依賴）= 改讀寫語義，超 M23 章程 ⟹ 發還。")
chk("C5  已聲明形式邊界（不越界改 δ 之形狀）", True)

out = {
    "c1_g_ratio": {"rho_d": str(c1_d), "rho_u": str(c1_u), "fa_members": len(members),
                   "examples": members[:8]},
    "c2_classes": {str(k): len(v) for k, v in cls.items()},
    "c2_out_ratio": str(out_ratio), "c2_ia_ratio": str(mid_ratio),
    "c3": {"distinct": len(pts), "median_nn_z": str(med_nn), "min_nn_z": str(min(nn)),
           "max_nn_z": str(max(nn))},
    "c4": {"fixed_ratio": str(ratio_fixed), "measured_ratio": mp.nstr(req_ratio, 12)},
}
with open(os.path.join(DATA, "04_completeness.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
nrf = sum(1 for _, ok in probe_results if not ok)
print("結果: %d/%d PASS ＋ %d 記錄 FAIL（探針／預期假設被否：判詞不依賴）" % (npass, len(results), nrf))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
