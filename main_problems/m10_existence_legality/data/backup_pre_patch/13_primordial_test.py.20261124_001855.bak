#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/13 — R12 第二獨立標籤：天然存在（primordial 丰度 IS>0）

前設問題:
  R10 之標籤係 NUBASE 半衰期（T½ ≥ 1e15 yr）。但判準宣稱之「永久」係**無限時間**，
  而 T½ 標籤本身亦有門檻任意性。故本輪用**更硬**之標籤：
      「此刻仍存在於自然界」= NUBASE 之 `IS=`（同位素丰度）欄 > 0。
  此標籤不需門檻：丰度 > 0 ⟹ 該核素此刻在地球上存在（來源可為 primordial 或宇宙成因，
  但皆意味 T½ 或生成率足以支撐其存在）。

事前預期（寫死）:
  E1 判準合法集 P ⊆ 天然存在集 Nat（每個永久核素都應天然存在）
  E2 Nat\\P 中 **Tier I 個案 = 0**（凡判不合法而天然存在者，必落近臨界或屏障帶）
  E3 Nat\\P 之 max rel ≤ 5.3e−4（與 R10 之失效半徑一致 ⟹ 兩個獨立標籤給同一半徑）

執行：/usr/bin/python3 scripts/13_primordial_test.py
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
print("M10 R12 — 第二獨立標籤：天然存在（IS>0）")
print("=" * 78)

# ---------------------------------------------------------------- 抽 IS=
IS = {}
for ln in open(os.path.join(REFS, "nubase2020.txt"), encoding="latin-1"):
    if len(ln) < 100 or ln.startswith("#"):
        continue
    if "IS=" not in ln:
        continue
    try:
        A = int(ln[0:3])
        Z = int(ln[4:7])
        if ln[7] != "0":
            continue
        m = re.search(r"IS=([0-9.]+)", ln)
        if m:
            IS[(Z, A)] = float(m.group(1))
    except ValueError:
        pass
print(f"  天然存在（IS>0）核素 = {len(IS)} 條")
ck("P.1 天然存在集規模 250–320", 250 <= len(IS) <= 320, f"{len(IS)} 條 ⟹ 與「~288 primordial」相符")

# ---------------------------------------------------------------- 質量與通道
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

ME_H, ME_N = F("7288.971064"), F("8071.31806")
BY_A = {}
for (z, a), v in meas.items():
    BY_A.setdefault(a, []).append((v, z))
for a in BY_A:
    BY_A[a].sort()
Z_OF_A = {a: tuple(sorted(z for _, z in v)) for a, v in BY_A.items()}


def me(z, a):
    return meas.get((z, a))


def amass(z, a):
    v = me(z, a)
    return None if v is None else v + F("931494.10242") * a


def min_open_rel(z, a):
    m = amass(z, a)
    if m is None:
        return None, None
    qs = []
    for (zz, aa), nm, sub in (((z + 1, a), "β⁻", F(0)), ((z - 1, a), "β⁺/EC", F(0)),
                              ((z - 2, a - 4), "α", me(2, 4)), ((z - 1, a - 1), "p", ME_H),
                              ((z, a - 1), "n", ME_N), ((z + 2, a), "2β⁻", F(0)),
                              ((z - 2, a), "2β⁺/ECEC", F(0))):
        v = me(zz, aa)
        if v is not None and sub is not None:
            qs.append((nm, me(z, a) - v - sub))
    if a >= 40:
        best = None
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


j9 = json.load(open(os.path.join(DATA, "10_channel_completion.json"), encoding="utf-8"))
P = {(int(z), int(a)) for z, a in j9["final"]["set"]}
Nat = set(IS)

# ---------------------------------------------------------------- E1
print("\n[E1] 合法集 ⊆ 天然存在集？")
inP_not_nat = sorted(P - Nat)
print(f"     |P| = {len(P)}；|P ∩ Nat| = {len(P & Nat)}；P \\ Nat = {len(inP_not_nat)}")
if inP_not_nat:
    print(f"     P 而不天然存在：{inP_not_nat}")
ck("E1【事前預期：P ⊆ Nat】", not inP_not_nat,
   f"{len(P & Nat)}/{len(P)} 條合法項天然存在 ⟹ " +
   ("判準之「永久」與「此刻仍在」一致" if not inP_not_nat
    else f"有 {len(inP_not_nat)} 條永久卻不天然存在（可能係純人工核素，須逐個查）"))

# ---------------------------------------------------------------- E2/E3
print("\n[E2] Nat\\P（天然存在卻判不合法）之 Tier 分佈")
NP = sorted(Nat - P)
rows = []
for (z, a) in NP:
    rel, nm = min_open_rel(z, a)
    tier = None if rel is None else ("III" if rel < 1e-4 else ("II" if rel < 1e-2 else "I"))
    rows.append({"Z": z, "A": a, "IS": IS[(z, a)], "rel": rel, "tier": tier, "channel": nm})
cnt_t = {}
for r in rows:
    cnt_t[r["tier"]] = cnt_t.get(r["tier"], 0) + 1
print(f"     |Nat\\P| = {len(NP)}；Tier 分佈 = {cnt_t}")
n_t1 = cnt_t.get("I", 0)
ck("E2【事前預期：Nat\\P 中 Tier I 個案 = 0】", n_t1 == 0,
   f"Tier I = {n_t1} 條 ⟹ 凡天然存在而判不合法者，皆落近臨界或屏障帶")
max_rel = max((r["rel"] for r in rows if r["rel"] is not None), default=None)
ck("E3【事前預期：max rel ≤ 5.3e−4】", max_rel is not None and max_rel <= 5.3e-4,
   f"max rel = {max_rel:.2e} ⟹ 與 R10（NUBASE T½ 標籤）之 5.21e−4 **一致**"
   if max_rel is not None and max_rel <= 5.3e-4 else f"max rel = {max_rel}")

# ---------------------------------------------------------------- 最強實證
print("\n[最強實證] 判「不合法」而此刻仍存在於自然界者")
top = sorted([r for r in rows if r["rel"] is not None], key=lambda r: -r["IS"])[:8]
for r in top:
    print(f"     Z={r['Z']:>3} A={r['A']:>3} 丰度 {r['IS']:>9.5f}%  "
          f"rel={r['rel']:.2e} tier={r['tier']} ch={r['channel']}")
print(f"     ⟹ {len(NP)} 條天然存在核素被判「不合法」；判準之 ⟹ 方向對其**不可操作**")
print("        （「不永久」≠「已衰變」—— 此即 D5.1(b) 之實證）")

# 反向：天然存在者中最短命者（若有 IS>0 而 T½ 極短 ⟹ 宇宙成因，非 primordial）
print("\n[註] 天然存在之最短命者（由 NUBASE T½ 排）")
nubT = {}
for ln in open(os.path.join(REFS, "nubase2020.txt"), encoding="latin-1"):
    if len(ln) < 100 or ln.startswith("#"):
        continue
    try:
        A = int(ln[0:3])
        Z = int(ln[4:7])
        if ln[7] != "0":
            continue
        t = ln[69:78].strip()
        u = ln[78:80].strip()
        if t and not t.startswith("stbl"):
            UN = {"ys": 1e-24 / 3.15576e7, "zs": 1e-21 / 3.15576e7,
                  "as": 1e-18 / 3.15576e7, "ps": 1e-12 / 3.15576e7,
                  "ns": 1e-9 / 3.15576e7, "us": 1e-6 / 3.15576e7,
                  "fs": 1e-15 / 3.15576e7,
                  "ms": 1e-3 / 3.15576e7, "s": 1 / 3.15576e7, "m": 1 / 525960.0,
                  "h": 1 / 8766.0, "d": 1 / 365.25, "y": 1.0, "ky": 1e3, "My": 1e6,
                  "Gy": 1e9, "Ty": 1e12, "Py": 1e15, "Ey": 1e18, "Zy": 1e21, "Yy": 1e24, "Yi": 1e24}
            nubT[(Z, A)] = float(t.lstrip("<>~ ").rstrip("#")) * UN.get(u, 1.0)
        else:
            nubT[(Z, A)] = float("inf")
    except (ValueError, IndexError):
        pass
short = sorted(((nubT[k], k) for k in Nat if k in nubT), key=lambda t: t[0])[:5]
for t, k in short:
    print(f"     Z={k[0]:>3} A={k[1]:>3} T½={t:.2e} yr（宇宙成因／短命天然）")
ck("P.2 天然存在集中最短命者 T½ ≥ 1e3 yr（宇宙成因之下限）", short and short[0][0] >= 1e3,
   f"最短 {short[0][0]:.2e} yr（Z={short[0][1][0]} A={short[0][1][1]}）⟹ "
   f"IS 標籤之最短命者仍達 10³ yr 級（宇宙成因），非瞬時核素")

# 單一同位素元素（IS = 100%）而被判不合法者
mono = [r for r in rows if abs(r["IS"] - 100.0) < 1e-9]
print(f"\n[單一同位素元素] IS = 100% 而被判「不合法」者 = {len(mono)} 條")
for r in sorted(mono, key=lambda r: r["A"]):
    print(f"     Z={r['Z']:>3} A={r['A']:>3} rel={r['rel']:.2e} tier={r['tier']} ch={r['channel']}")
ck("P.3 單一同位素元素而被判不合法者皆落 Tier II/III（無 Tier I）",
   all(r["tier"] in ("II", "III") for r in mono),
   f"{len(mono)} 條；其元素在自然界 100% 由此核素構成 ⟹ 判準之 ⟹ 方向最強反例（但皆屏障類）")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "13_primordial_test.py", "role": "R12 第二獨立標籤（IS>0）",
       "n_nat": len(Nat), "n_P": len(P), "n_P_and_nat": len(P & Nat),
       "P_not_nat": [[z, a] for z, a in inP_not_nat],
       "n_Nat_not_P": len(NP), "tier_dist": cnt_t, "n_tier1": n_t1,
       "max_rel": max_rel,
       "n_mono": len(mono), "mono": [{"Z": r["Z"], "A": r["A"], "rel": r["rel"],
                                      "tier": r["tier"], "channel": r["channel"]} for r in mono],
       "items": rows,
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "13_primordial_test.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/13_primordial_test.json 已寫")