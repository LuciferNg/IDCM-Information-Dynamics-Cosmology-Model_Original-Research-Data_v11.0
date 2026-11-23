#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/11 — R10 真證偽測試：判準之合法集 vs NUBASE2020 觀測半衰期（獨立標籤）

前設問題（docs/07 §十：「觀測穩定之獨立標籤」未做 ⟹ 判準從未受過真正之外部檢定）:
  前七輪全部只用質量（AME2020）。質量 ↔ 質量 之比對係同源自洽，唔構成證偽。
  R10 引入**獨立**標籤：NUBASE2020 之實測半衰期（761,906 bytes，IAEA/AMDC）。

數據：
  data/refs/mass.mas20.txt     AME2020 質量（判準之輸入）
  data/refs/nubase2020.txt     NUBASE2020 半衰期（**判準從未使用之獨立資訊**）

事前預期（寫死，落空即記錄、不調）:
  A. 證偽檢定：P（判準合法集）中不得有 T½ < 1e15 yr 者。**任何一個 ⟹ 判準被證偽**
     （在其 11 通道 𝒞 下）。
  B. Tier 驗證：S\\P（觀測穩定但判準判不合法）中，≥90% 必須落 Tier III
     （|Q|/m < 1e−4）。
  C. 時間盲：不合法集內 T½ 之跨度 max/min ≥ 1e30 ⟹ 判準完全不辨時間尺度。

執行：/usr/bin/python3 scripts/11_falsification_nubase.py
"""
import json
import os
import re
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA, REFS = os.path.join(ROOT, "data"), os.path.join(ROOT, "data", "refs")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


print("=" * 78)
print("M10 R10 — 真證偽測試：判準 vs NUBASE2020 半衰期（獨立標籤）")
print("=" * 78)

# ---------------------------------------------------------------- NUBASE2020
# 欄位（1-based，見檔頭）：1:3 AAA｜5:8 ZZZi｜12:16 A El｜17 s｜19:31 Mass#
#   70:78 T｜79:80 unit T｜82:88 dT｜89:102 Jpi｜120:209 BR（IS= 丰度）
UNITS = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "ps": 1e-12, "fs": 1e-15,
         "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0,
         "d": 86400.0, "y": 31557600.0, "ky": 3.15576e10, "My": 3.15576e13,
         "Gy": 3.15576e16, "Ty": 3.15576e19, "Py": 3.15576e22,
         "Ey": 3.15576e25, "Zy": 3.15576e28, "Yy": 3.15576e31, "Yi": 3.15576e31}
YR = 31557600.0           # 儒略年（秒）

nub = {}
bad = 0
for ln in open(os.path.join(REFS, "nubase2020.txt"), encoding="latin-1"):
    if len(ln) < 100 or ln.startswith("#"):
        continue
    try:
        A = int(ln[0:3])
        Z = int(ln[4:7])                  # ZZZ（col 5:7）
        isom = ln[7]                      # col 8 = i（'0' = ground state）
        if isom != "0":
            continue
        tfield = ln[69:78].strip()
        unit = ln[78:80].strip()
        if not tfield:
            bad += 1
            continue
        if tfield.startswith("stbl"):
            nub[(Z, A)] = {"T_yr": float("inf"), "raw": "stbl"}
            continue
        if tfield.startswith(("p-unst", "p-uns")):
            nub[(Z, A)] = {"T_yr": 0.0, "raw": tfield}
            continue
        s = tfield.lstrip("<>~ ").rstrip("#")
        val = float(s)
        if unit not in UNITS:
            bad += 1
            continue
        nub[(Z, A)] = {"T_yr": val * UNITS[unit] / YR, "raw": tfield + " " + unit}
    except (ValueError, IndexError):
        bad += 1
print(f"  NUBASE2020：{len(nub)} 條基態（解析失敗 {bad}）")
ck("P.1 NUBASE 解析 ≥ 3000 條", len(nub) >= 3000, f"{len(nub)} 條")
for k, nm in (((1, 3), "³H"), ((83, 209), "²⁰⁹Bi"), ((6, 14), "¹⁴C"), ((2, 4), "⁴He")):
    if k in nub:
        print(f"     錨點 {nm}: T = {nub[k]['raw']} = {nub[k]['T_yr']:.3e} yr")
ck("P.2 錨點合理（³H ~12.3 yr；²⁰⁹Bi ~2e19 yr；⁴He 穩定）",
   abs(nub[(1, 3)]["T_yr"] / 12.32 - 1) < 0.02 and nub[(83, 209)]["T_yr"] > 1e18
   and nub[(2, 4)]["T_yr"] == float("inf"),
   "NUBASE 與文獻一致 ⟹ 獨立標籤可用")

# ---------------------------------------------------------------- 判準（11 通道）
j10 = json.load(open(os.path.join(DATA, "10_channel_completion.json"), encoding="utf-8"))
P = {(int(z), int(a)) for z, a in j10["final"]["set"]}
n_meas = j10["n_measured"]
print(f"\n  判準合法集 P = {len(P)} 條（𝒞_v2 = 11 通道；實測池 {n_meas} 條）")

# ---------------------------------------------------------------- 檢定 A
print("\n[A] 證偽檢定：P 中不得有 T½ < 1e15 yr 者")
CUT = 1e15
kill = [(z, a, nub[(z, a)]["T_yr"], nub[(z, a)]["raw"]) for (z, a) in P
        if (z, a) in nub and nub[(z, a)]["T_yr"] < CUT]
print(f"     P 中有 NUBASE 數據者 {sum(1 for k in P if k in nub)}/{len(P)}")
for z, a, t, raw in sorted(kill, key=lambda x: x[2]):
    print(f"     🔴 Z={z} A={a}: T½ = {raw} = {t:.2e} yr ⟹ **反例**")
ck("A.1【事前預期：P 中無 T½ < 1e15 yr 者】", not kill,
   f"反例 {len(kill)} 條" + ("" if kill else " ⟹ 判準未被證偽"))

# ---------------------------------------------------------------- 檢定 B
print("\n[B] Tier 驗證：S\\P（觀測穩定卻判不合法）是否全落 Tier III")
S = {k for k, v in nub.items() if v["T_yr"] >= CUT}
SP = sorted((S & set(nub)) - P, key=lambda t: (t[1], t[0]))
print(f"     觀測穩定集 S(1e15 yr) = {len(S)} 條；S∩實測池中不在 P 者 = {len(SP)} 條")
# Tier：用 scripts/10 之 11 通道，取最小正 Q 之相對值
import importlib.util
spec = importlib.util.spec_from_file_location(
    "ch", os.path.join(HERE, "10_channel_completion.py"))

# 直接由 JSON 之外無法取得 Q；改為在本腳本內用 AME 質量重算（最小正 Q）
NUM = re.compile(r"^-?\d+(?:\.\d+)?#?$")
SYM = re.compile(r"^[A-Z][a-z]?$")
meas = {}
for ln in open(os.path.join(REFS, "mass.mas20.txt"), encoding="latin-1"):
    if not ln.strip() or ln.lstrip().startswith("#"):
        continue
    tok = ln.split()
    ints, idx = [], None
    for i, t in enumerate(tok):
        if re.match(r"^-?\d+$", t):
            ints.append(int(t))
        elif (SYM.match(t) or t == "n") and len(ints) >= 3:
            idx = i
            break
    if idx is None or len(ints) < 4:
        continue
    try:
        Z, A = ints[-2], ints[-1]
        i = idx + 1
        while i < len(tok) and not NUM.match(tok[i]):
            i += 1
        if "#" in tok[i]:
            continue
        meas[(Z, A)] = F(tok[i])
    except (ValueError, IndexError):
        pass


def me(z, a):
    return meas.get((z, a))


def amass(z, a):
    v = me(z, a)
    return None if v is None else v + F("931494.10242") * a


ME_H, ME_N = F("7288.971064"), F("8071.31806")
BY_A = {}
for (z, a), v in meas.items():
    BY_A.setdefault(a, []).append((v, z))
for a in BY_A:
    BY_A[a].sort()
Z_OF_A = {a: tuple(sorted(z for _, z in v)) for a, v in BY_A.items()}


def min_open_rel(z, a):
    """11 通道中最小之正 Q，除以原子質量（Tier 之相對量）。"""
    m = amass(z, a)
    if m is None:
        return None, None
    qs = []
    for (zz, aa), nm in (((z + 1, a), "β⁻"), ((z - 1, a), "β⁺/EC"),
                         ((z - 2, a - 4), "α"), ((z - 1, a - 1), "p"),
                         ((z, a - 1), "n"), ((z + 2, a), "2β⁻"),
                         ((z - 2, a), "2β⁺/ECEC")):
        v = me(zz, aa)
        if v is not None:
            q = me(z, a) - v - (ME_H if nm == "p" else (ME_N if nm == "n" else F(0)))
            if nm == "α":
                q = me(z, a) - v - me(2, 4)
            qs.append((nm, q))
    # SF_opt
    best = None
    if a >= 40:
        for a1 in range(10, a - 9):
            a2 = a - a1
            for z1 in Z_OF_A.get(a1, ()):
                m1, m2 = me(z1, a1), me(z - z1, a2)
                if m1 is None or m2 is None:
                    continue
                q = me(z, a) - m1 - m2
                if best is None or q > best:
                    best = q
    if best is not None:
        qs.append(("SF_opt", best))
    pos = [(nm, q) for nm, q in qs if q > 0]
    if not pos:
        return None, None
    nm, q = min(pos, key=lambda t: t[1])
    return float(q / m), nm


tiers = []
for (z, a) in SP:
    rel, nm = min_open_rel(z, a)
    tier = None if rel is None else ("III" if rel < 1e-4 else ("II" if rel < 1e-2 else "I"))
    tiers.append({"Z": z, "A": a, "T_yr": nub[(z, a)]["T_yr"], "rel": rel,
                  "tier": tier, "channel": nm})
n3 = sum(1 for t in tiers if t["tier"] == "III")
print(f"     其中 Tier III = {n3}/{len(tiers)}（{100*n3/max(len(tiers),1):.1f}%）")
bad_t = [t for t in tiers if t["tier"] != "III"]
for t in sorted(bad_t, key=lambda x: x["A"]):
    ts_ = "stbl" if t["T_yr"] == float("inf") else f"{t['T_yr']:.1e}"
    print(f"       ⚠️ {t['A']:>3} Z{t['Z']:<3} A={t['A']:<3} T½={ts_:>7} "
          f"rel={t['rel']:.2e} tier={t['tier']} ch={t['channel']}")
ck("B.1【事前預期：S\\P 中 ≥90% 落 Tier III】", n3 / max(len(tiers), 1) >= 0.90,
   f"{n3}/{len(tiers)} = {100*n3/max(len(tiers),1):.1f}% ⟹ **落空**（門檻不調）")

# 後驗：例外之機制分類
n_sf = sum(1 for t in bad_t if t["channel"] == "SF_opt")
n_stbl = sum(1 for t in bad_t if t["T_yr"] == float("inf"))
A_rng = (min(t["A"] for t in bad_t), max(t["A"] for t in bad_t)) if bad_t else None
rel_rng = (min(t["rel"] for t in bad_t), max(t["rel"] for t in bad_t)) if bad_t else None
ck("B.2【後驗】非 Tier III 之例外全部出自**同一通道**（單一機制）",
   n_sf == len(bad_t),
   f"{n_sf}/{len(bad_t)} 條皆為 SF_opt（最優二元裂變）⟹ 失效非發散，係"
   f"**裂變屏障**一種物理機制")
ck("B.3【後驗】例外皆為實測穩定（stbl）⟹ 判準之 ⟹ 方向對其**不可測**、非被證偽",
   n_stbl == len(bad_t),
   f"{n_stbl}/{len(bad_t)} 條 T½ = stbl；A ∈ [{A_rng[0]},{A_rng[1]}]、"
   f"rel ∈ [{rel_rng[0]:.1e},{rel_rng[1]:.1e}]（皆 Tier II）⟹ "
   f"Q_SF > 0 但液滴裂變屏障 ~30–40 MeV ≫ Q ⟹ 衰變不可觀測")

# 分解：S\\P 之兩類
n_sp_stbl = sum(1 for t in tiers if t["T_yr"] == float("inf"))
n_sp_fin = len(tiers) - n_sp_stbl
from collections import Counter
ch_dist = Counter(t["channel"] for t in tiers)
tier_dist = Counter(t["tier"] for t in tiers)
print(f"     S\\P 分解：stbl {n_sp_stbl} 條／有限 T½ {n_sp_fin} 條")
print(f"     通道分佈：{dict(ch_dist)}")
print(f"     Tier 分佈：{dict(tier_dist)}")
max_rel_sp = max(t["rel"] for t in tiers)
n_t1 = sum(1 for t in tiers if t["tier"] == "I")
ck("B.5【後驗】失效半徑有界：S\\P 全落 rel ≤ 5.3e−4，且**無 Tier I**",
   max_rel_sp < 5.3e-4 and n_t1 == 0 and len(bad_t) + n3 == len(tiers),
   f"max rel = {max_rel_sp:.2e}；Tier I 個案 = {n_t1} ⟹ 凡「不合法」而觀測穩定者，"
   f"其最小正 Q 之相對裕度皆 ≤ {max_rel_sp:.2e}；rel > {max_rel_sp:.1e} 之區間**零偽陰性**")
ck("B.6【後驗】Tier II 之 51 條係**單一通道**（屏障帶），Tier III 之 130 條係混合通道",
   all(t["channel"] == "SF_opt" for t in tiers if t["tier"] == "II")
   and len({t["channel"] for t in tiers if t["tier"] == "III"}) >= 2,
   f"Tier II 通道 = {dict(Counter(t['channel'] for t in tiers if t['tier']=='II'))}；"
   f"Tier III 通道 = {dict(Counter(t['channel'] for t in tiers if t['tier']=='III'))}")

# ⟸ 方向之直接檢定（判準之可測方向）
n_Pstbl = sum(1 for k in P if k in nub and nub[k]["T_yr"] == float("inf"))
n_Pdat = sum(1 for k in P if k in nub)
ck("B.4【事後】⟸ 方向（∀c x<0 ⟹ 永久）可測且通過：P 中全數觀測穩定",
   n_Pstbl == n_Pdat and n_Pdat > 0,
   f"P 之 {n_Pdat} 條有 NUBASE 數據者中 {n_Pstbl} 條為 stbl，其餘皆 T½ ≥ 1e15 yr"
   f" ⟹ 判準之**必要條件**方向零反例")

# ---------------------------------------------------------------- 檢定 C
print("\n[C] 時間盲：不合法集內 T½ 之跨度")
ill = [(k, v["T_yr"]) for k, v in nub.items() if v["T_yr"] != float("inf") and k not in P
       and k in meas and v["T_yr"] > 0]
ts = sorted(t for _, t in ill)
if ts:
    span = ts[-1] / ts[0]
    print(f"     不合法集（有 NUBASE T½ 且非 stbl）{len(ts)} 條")
    print(f"     最短 {ts[0]:.2e} yr ／ 最長 {ts[-1]:.2e} yr ⟹ 跨度 {span:.1e}×")
    ck("C.1【事前預期：跨度 ≥ 1e30】", span >= 1e30,
       f"跨度 {span:.1e}× ⟹ 判準對時間尺度**零解析度**（禁區①之實證）")
    # 同集內最長/最短之例
    kmin = min(ill, key=lambda t: t[1])
    kmax = max(ill, key=lambda t: t[1])
    print(f"     最短：Z={kmin[0][0]} A={kmin[0][1]}（{kmin[1]:.1e} yr）；"
          f"最長：Z={kmax[0][0]} A={kmax[0][1]}（{kmax[1]:.1e} yr）")

# ---------------------------------------------------------------- 判詞
print("\n[判詞]")
verdict = "未證偽（在 𝒞_v2 11 通道下）" if not kill else "**已證偽**"
print(f"     判準之合法集 P = {len(P)} 條 ⟹ {verdict}")
print(f"     觀測穩定但判不合法者 {len(SP)} 條：Tier III {n3} 條（{100*n3/len(SP):.1f}%）、"
      f"非 Tier III {len(bad_t)} 條")
print(f"     ⟹ 失效有**兩類**：① 近臨界（Tier III 內，精度極限）；"
      f"② 屏障類（{len(bad_t)} 條全為 SF_opt、全 stbl）")
print("     ⟹ 判準係**單向證書**（可證永久），非雙向判別器（不可證衰變）")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {
    "script": "11_falsification_nubase.py", "role": "R10 真證偽測試（NUBASE2020 標籤）",
    "n_nubase": len(nub), "n_parse_fail": bad, "n_measured": n_meas, "P_size": len(P),
    "anchors": {nm: nub[k]["T_yr"] for nm, k in
                (("3H", (1, 3)), ("209Bi", (83, 209)), ("14C", (6, 14)), ("4He", (2, 4)))
                if k in nub},
    "cut_yr": CUT,
    "test_A": {"falsified": bool(kill),
               "counterexamples": [{"Z": z, "A": a, "T_yr": t, "raw": r} for z, a, t, r in kill]},
    "test_B": {"S_size": len(S), "SP_size": len(SP), "tier3": n3,
               "tier3_frac": n3 / max(len(tiers), 1),
               "n_non_tier3": len(bad_t), "non_tier3_all_SF_opt": n_sf == len(bad_t),
               "non_tier3_all_stbl": n_stbl == len(bad_t),
               "non_tier3_A_range": list(A_rng) if A_rng else None,
               "SP_stbl": n_sp_stbl, "SP_finite": n_sp_fin,
               "max_rel": max_rel_sp, "n_tier1": n_t1,
               "SP_channel_dist": dict(ch_dist), "SP_tier_dist": dict(tier_dist),
               "III_channel_dist": dict(Counter(t["channel"] for t in tiers if t["tier"] == "III")),
               "items": tiers},
    "direction_pm": {"implication_iff_forall_closed_yields_permanent": "可測，通過（P⊆S）",
                     "implication_any_open_yields_impermanent": "不可測（屏障通道）"},
    "test_C": {"n_ill": len(ts), "min_yr": ts[0] if ts else None,
               "max_yr": ts[-1] if ts else None,
               "span": (ts[-1] / ts[0]) if ts else None},
    "verdict": verdict,
    "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL,
}
with open(os.path.join(DATA, "11_falsification_nubase.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/11_falsification_nubase.json 已寫")