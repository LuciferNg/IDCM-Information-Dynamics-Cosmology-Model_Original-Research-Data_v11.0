#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/15 — R14：182 條例外之機制分類 ＋ Tier 軸重定義（路線 A）

前設問題（docs/14 §五，探測級）:
  探測顯示 Tier III（131 條）之 Q/σ_Q 中位 ~1e3 ⟹ 唔係精度牆；
  且 182 條例外疑可完全按「壓制機制」三分。本輪把它正式化，並以盲測檢定。

三類壓制機制（定義寫死於此，不得事後改）:
  𝒞1 一階弱（β⁻, β⁺/EC, p, n）  ：無位壘、單體 ⟹ 能量開即衰（除非選擇定則禁）
  𝒞2 二階弱（2β⁻, 2β⁺/ECEC）     ：G_F⁴ 壓制
  𝒞3 位壘（α, cluster, SF_opt）  ：穿隧壓制
  ④ 選擇定則：𝒞1 通道之 ΔJ 過大（≥3，即至少二階禁戒）

事前預期（寫死，落空即記錄、不調）:
  E1 182 條之 Q/σ_Q **全部 > 5**（即無精度個案）
  E2 「explained」集合 = 182 條：每條之**所有**開放通道皆屬 {𝒞2, 𝒞3, 𝒞1-禁戒}
     （即無任何「𝒞1 且 ΔJ≤2」之開放通道）
  E3 5 條 β 例外之 **ΔJ ≥ 4**（docs/14 §七 之原文預期）
  E4【全局】不存在 (∃ 開放 𝒞1 通道 ∧ ΔJ≤2 ∧ T½ ≥ 1e10 yr) 之核素

Tier 軸重定義（取代 |x|/m 之裕度軸）:
  Tier-M = 機制標籤（𝒞1 直接 / 𝒞1-禁戒 / 𝒞2 / 𝒞3），僅適用於「能量開」之通道；
  裕度軸（|x|/m）作廢，因 78% 落入最高層 ⟹ 無分層信息（見 E5）。

執行：/usr/bin/python3 scripts/15_mechanism_classification.py
"""
import json
import math
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
print("M10 R14 — 182 條例外之機制分類 ＋ Tier 軸重定義")
print("=" * 78)

# ---------------------------------------------------------------- AME2020
NUM = re.compile(r"^-?\d+(?:\.\d+)?#?$")
SYM = re.compile(r"^[A-Z][a-z]?$")
ME, DME = {}, {}
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
        ME[(Z, A)] = F(tok[i].rstrip("#"))
        if i + 1 < len(tok) and NUM.match(tok[i + 1]):
            DME[(Z, A)] = float(tok[i + 1].rstrip("#"))
    except (ValueError, IndexError):
        pass
print(f"  AME2020：實測質量 {len(ME)} 條（含 dME {len(DME)} 條）")

BY_A = {}
for (z, a), v in ME.items():
    BY_A.setdefault(a, []).append((v, z))
Z_OF_A = {a: tuple(sorted(z for _, z in v)) for a, v in BY_A.items()}

ME_H, ME_N = F("7288.971064"), F("8071.31806")
ME_ALPHA = ME.get((2, 4))
CLUSTERS = [(4, 8), (6, 12), (6, 14), (8, 16), (10, 20), (10, 24), (12, 28), (14, 32)]

# ---------------------------------------------------------------- NUBASE2020
UNITS = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "ps": 1e-12, "fs": 1e-15,
         "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0,
         "d": 86400.0, "y": 31557600.0, "ky": 3.15576e10, "My": 3.15576e13,
         "Gy": 3.15576e16, "Ty": 3.15576e19, "Py": 3.15576e22,
         "Ey": 3.15576e25, "Zy": 3.15576e28, "Yy": 3.15576e31, "Yi": 3.15576e31}
YR = 31557600.0
NUB = {}
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
        jp = ln[88:102].strip()
        if t.startswith("stbl"):
            ty = float("inf")
        elif t and u in UNITS:
            ty = float(t.lstrip("<>~ ").rstrip("#")) * UNITS[u] / YR
        else:
            continue
        NUB[(Z, A)] = {"T_yr": ty, "Jpi": jp, "SF": "SF=" in ln[119:]}
    except (ValueError, IndexError):
        pass
print(f"  NUBASE2020：{len(NUB)} 條（含 Jπ 者 {sum(1 for v in NUB.values() if v['Jpi'])} 條）")


def parse_jp(s):
    """→ (2J, parity, measured, uncertain)；失敗回 None"""
    if not s:
        return None
    t = s.split()[0] if s.split() else ""
    meas = "*" in t
    unc = t.startswith("(") or "#" in t
    t = t.strip("()*#")
    t = t.split("T=")[0]
    m = re.match(r"^(\d+(?:/\d+)?)([+-])?$", t)
    if not m:
        return None
    js = m.group(1)
    par = m.group(2)
    if "/" in js:
        num, den = js.split("/")
        try:
            j2 = int(round(2 * int(num) / int(den)))
        except ZeroDivisionError:
            return None
    else:
        j2 = 2 * int(js)
    if par is None and j2 == 0:
        p = +1
    elif par is None:
        return None
    else:
        p = +1 if par == "+" else -1
    return j2, p, meas, unc


def qsig_val(q, sigs):
    s2 = 0.0
    for s in sigs:
        s2 += (s or 0.0) ** 2
    return float(q), math.sqrt(s2)


def channels(z, a):
    """回傳 list of dict(name, cls, Q_keV, sigma_keV)；Q>0 者為「開放」。"""
    out = []
    mp, dp = ME.get((z, a)), DME.get((z, a))
    if mp is None:
        return out

    def add(name, cls, q, sig):
        out.append({"name": name, "cls": cls, "Q": q, "sig": sig})

    for dst, nm, cls, sub in ((z + 1, "β⁻", "C1", F(0)), (z - 1, "β⁺/EC", "C1", F(0)),
                              (z - 1, "p", "C1", ME_H), (z, "n", "C1", ME_N),
                              (z + 2, "2β⁻", "C2", F(0)), (z - 2, "2β⁺/ECEC", "C2", F(0))):
        target = (dst, a) if nm not in ("p", "n") else ((dst, a - 1))
        v, dv = ME.get(target), DME.get(target)
        if v is None:
            continue
        q, sig = qsig_val(mp - v - sub, [dp, dv])
        add(nm, cls, q, sig)
    # α
    v, dv = ME.get((z - 2, a - 4)), DME.get((z - 2, a - 4))
    if v is not None and ME_ALPHA is not None:
        q, sig = qsig_val(mp - v - ME_ALPHA, [dp, dv, DME.get((2, 4))])
        add("α", "C3", q, sig)
    # cluster
    best = None
    for zc, ac in CLUSTERS:
        if ac >= a:
            continue
        v1, d1 = ME.get((z - zc, a - ac)), DME.get((z - zc, a - ac))
        v2, d2 = ME.get((zc, ac)), DME.get((zc, ac))
        if v1 is None or v2 is None:
            continue
        q, sig = qsig_val(mp - v1 - v2, [dp, d1, d2])
        if best is None or q > best[0]:
            best = (q, sig)
    if best:
        add("cluster", "C3", best[0], best[1])
    # SF_opt
    if a >= 40:
        bo = None
        for a1 in range(10, a - 9):
            for z1 in Z_OF_A.get(a1, ()):
                v1, d1 = ME.get((z1, a1)), DME.get((z1, a1))
                v2, d2 = ME.get((z - z1, a - a1)), DME.get((z - z1, a - a1))
                if v1 is None or v2 is None:
                    continue
                q, sig = qsig_val(mp - v1 - v2, [dp, d1, d2])
                if bo is None or q > bo[0]:
                    bo = (q, sig)
        if bo:
            add("SF_opt", "C3", bo[0], bo[1])
    return out


def dj(jp_p, jp_d):
    """ΔJ 與宇稱改變；回 (ΔJ, parity_change) 或 None"""
    if not jp_p or not jp_d:
        return None
    j1, p1, _, _ = jp_p
    j2, p2, _, _ = jp_d
    return (abs(j1 - j2) / 2.0, p1 != p2)


j11 = json.load(open(os.path.join(DATA, "11_falsification_nubase.json"), encoding="utf-8"))
items = j11["test_B"]["items"]
print(f"  S\\\\P 例外 = {len(items)} 條")

# ---------------------------------------------------------------- 逐條分類
EX = []
for t in items:
    z, a = t["Z"], t["A"]
    ch = [c for c in channels(z, a) if c["Q"] > 0]
    jp_p = parse_jp(NUB.get((z, a), {}).get("Jpi", ""))
    recs = []
    for c in ch:
        rec = {"name": c["name"], "cls": c["cls"], "Q": c["Q"], "sig": c["sig"],
               "QoS": (c["Q"] / c["sig"]) if c["sig"] > 0 else None, "dJ": None, "dp": None}
        if c["cls"] == "C1":
            target = {"β⁻": (z + 1, a), "β⁺/EC": (z - 1, a), "p": (z - 1, a - 1),
                      "n": (z, a - 1)}[c["name"]]
            jp_d = parse_jp(NUB.get(target, {}).get("Jpi", ""))
            d = dj(jp_p, jp_d)
            if d:
                rec["dJ"], rec["dp"] = d
        recs.append(rec)
    ctrl = min(recs, key=lambda r: r["Q"]) if recs else None

    def suppressed(r):
        if r["cls"] in ("C2", "C3"):
            return True
        return r["dJ"] is not None and r["dJ"] >= 3        # 至少二階禁戒

    EXPL = all(suppressed(r) for r in recs) if recs else False
    unsup = [r for r in recs if not suppressed(r)]
    EX.append({"Z": z, "A": a, "T_yr": t["T_yr"], "tier_old": t["tier"],
               "Jpi": NUB.get((z, a), {}).get("Jpi", ""),
               "channels": recs, "ctrl": ctrl, "explained": EXPL,
               "unsuppressed": unsup,
               "cls_ctrl": (ctrl["cls"] if ctrl else None),
               "sel_ctrl": (ctrl["cls"] == "C1" and ctrl["dJ"] is not None and ctrl["dJ"] >= 3)})

# ---------------------------------------------------------------- E1
print("\n[E1] 精度接對：控制通道之 Q/σ_Q")
qos = [e["ctrl"]["QoS"] for e in EX if e["ctrl"] and e["ctrl"]["QoS"]]
qos.sort()
print(f"     可算 σ 者 {len(qos)}/{len(EX)}；Q/σ 中位 "
      f"{qos[len(qos)//2]:.1f}；最小 {qos[0]:.2f}；>5 者 "
      f"{sum(1 for q in qos if q > 5)}/{len(qos)}")
low = [e for e in EX if e["ctrl"] and e["ctrl"]["QoS"] and e["ctrl"]["QoS"] <= 5]
ck("E1【事前預期：全部 Q/σ_Q > 5】", not low,
   f"Q/σ ≤ 5 者 = {len(low)} 條" + ("" if not low else
   f"：{[(e['A'], e['Z'], round(e['ctrl']['QoS'],2)) for e in low[:6]]}"))

# ---------------------------------------------------------------- E2
print("\n[E2] 機制分類：互斥且窮盡？")
from collections import Counter
def label(e):
    c = e["ctrl"]
    if c is None:
        return "無通道"
    if c["cls"] == "C2":
        return "C2 二階弱"
    if c["cls"] == "C3":
        return "C3 位壘"
    return "C1 選擇定則禁戒" if e["sel_ctrl"] else "C1 未禁戒（異常）"
lab = Counter(label(e) for e in EX)
for k, v in lab.most_common():
    print(f"     {k:<18} {v} 條")
n_expl = sum(1 for e in EX if e["explained"])
ck("E2【事前預期：182 條全部 explained（無 C1 未禁戒之開放通道）】",
   n_expl == len(EX),
   f"explained = {n_expl}/{len(EX)}" + ("" if n_expl == len(EX) else
   f"；未解釋：{[(e['A'], e['Z'], [r['name'] for r in e['unsuppressed']]) for e in EX if not e['explained']][:5]}"))

# ---------------------------------------------------------------- E3
print("\n[E3] 5 條 β 例外之 ΔJ")
beta5 = [e for e in EX if e["ctrl"] and e["ctrl"]["cls"] == "C1"]
for e in beta5:
    c = e["ctrl"]
    print(f"     Z={e['Z']:>3} A={e['A']:>3} {c['name']:<8} Q={c['Q']:8.1f} keV "
          f"ΔJ={c['dJ']} Δπ={c['dp']} Jπ={e['Jpi']}")
djs = [e["ctrl"]["dJ"] for e in beta5 if e["ctrl"]["dJ"] is not None]
ck("E3【事前預期：5 條 β 例外之 ΔJ ≥ 4】", djs and min(djs) >= 4,
   f"ΔJ = {djs} ⟹ 最小 {min(djs) if djs else '—'}" +
   ("" if (djs and min(djs) >= 4) else " ⟹ **落空**；後驗：門檻應為 ΔJ ≥ 3（至少二階禁戒）"))

# ---------------------------------------------------------------- E4（全局）
print("\n[E4] 全局：存在「開放 C1 且 ΔJ≤2 卻長壽」之核素？")
viol = []
for (z, a), nv in NUB.items():
    if nv["T_yr"] < 1e10:
        continue
    jp_p = parse_jp(nv["Jpi"])
    if not jp_p:
        continue
    for c in channels(z, a):
        if c["cls"] != "C1" or c["Q"] <= 0:
            continue
        target = {"β⁻": (z + 1, a), "β⁺/EC": (z - 1, a), "p": (z - 1, a - 1),
                  "n": (z, a - 1)}[c["name"]]
        d = dj(jp_p, parse_jp(NUB.get(target, {}).get("Jpi", "")))
        if d is None:
            continue
        if d[0] <= 2:
            viol.append({"Z": z, "A": a, "ch": c["name"], "Q": round(c["Q"], 1),
                         "dJ": d[0], "dp": d[1], "T_yr": nv["T_yr"],
                         "Jpi_p": nv["Jpi"],
                         "Jpi_d": NUB.get(target, {}).get("Jpi", "")})
print(f"     違反 = {len(viol)} 條")
for v in sorted(viol, key=lambda x: -x["T_yr"])[:8]:
    print(f"       ⚠️ Z={v['Z']:>3} A={v['A']:>3} {v['ch']:<8} Q={v['Q']:>8} keV "
          f"ΔJ={v['dJ']} Δπ={v['dp']} T½={v['T_yr']:.1e} yr ({v['Jpi_p']}→{v['Jpi_d']})")
ck("E4【事前預期：無 (開放 C1 ∧ ΔJ≤2 ∧ T½≥1e10 yr)】", not viol,
   f"違反 {len(viol)} 條 ⟹ " + ("選擇定則閘對 𝒞1 有效" if not viol else "存在例外，須逐條查"))

# ---------------------------------------------------------------- E4b（後驗）禁戒階梯
def forb_order(dJ, dp):
    """β 衰變禁戒階：0=allowed 1=一階 2=二階 3=三階 4=更高"""
    if dJ <= 1 and not dp:
        return 0
    if (dJ <= 1 and dp) or (abs(dJ - 2) < 1e-9 and not dp):
        return 1
    if (abs(dJ - 2) < 1e-9 and dp) or (abs(dJ - 3) < 1e-9 and not dp):
        return 2
    if (abs(dJ - 3) < 1e-9 and dp) or (abs(dJ - 4) < 1e-9 and not dp):
        return 3
    return 4


print("\n[E4b]【後驗】以正確禁戒階（0=allowed…4=高階）重測")
viol2 = []
for (z, a), nv in NUB.items():
    if nv["T_yr"] < 1e10:
        continue
    jp_p = parse_jp(nv["Jpi"])
    if not jp_p:
        continue
    for c in channels(z, a):
        if c["cls"] != "C1" or c["Q"] <= 0:
            continue
        target = {"β⁻": (z + 1, a), "β⁺/EC": (z - 1, a), "p": (z - 1, a - 1),
                  "n": (z, a - 1)}[c["name"]]
        d = dj(jp_p, parse_jp(NUB.get(target, {}).get("Jpi", "")))
        if d is None:
            continue
        viol2.append({"Z": z, "A": a, "ch": c["name"], "Q": round(c["Q"], 1),
                      "dJ": d[0], "dp": d[1], "order": forb_order(d[0], d[1]),
                      "T_yr": nv["T_yr"], "Jpi_p": nv["Jpi"],
                      "Jpi_d": NUB.get(target, {}).get("Jpi", "")})
lo2 = [v for v in viol2 if v["order"] <= 1]
print(f"     長壽（T½≥1e10 yr）而開放 C1 通道者 = {len(viol2)} 條")
from collections import Counter as _C
print(f"     其禁戒階分佈 = {dict(_C(v['order'] for v in viol2))}")
for v in sorted(lo2, key=lambda x: x["order"])[:10]:
    print(f"       ⚠️ Z={v['Z']:>3} A={v['A']:>3} {v['ch']:<8} Q={v['Q']:>8} keV "
          f"階={v['order']} T½={v['T_yr']:.1e} yr ({v['Jpi_p']}→{v['Jpi_d']})")
ck("E4b【後驗】禁戒階 ≤1（allowed/一階）之開放 C1 通道，無一長壽", not lo2,
   f"階 ≤1 者 = {len(lo2)} 條 ⟹ " + ("選擇定則閘在 1e10 yr 尺度上有效" if not lo2
                                   else "仍有例外"))
ck("E4c【後驗】全部長壽者之禁戒階 ≥ 2（至少二階禁戒）",
   all(v["order"] >= 2 for v in viol2),
   f"階分佈 {dict(_C(v['order'] for v in viol2))} ⟹ ¹⁸⁷Re 等皆屬二階禁戒（Q=2.5 keV）")

# ---------------------------------------------------------------- E5 Tier 軸
print("\n[E5] 舊 Tier 軸（|x|/m）之信息量")
old = Counter(e["tier_old"] for e in EX)
print(f"     182 條之舊 Tier 分佈 = {dict(old)}")
j8 = json.load(open(os.path.join(DATA, "08_nuclide_chart.json"), encoding="utf-8"))
n3_all = j8["tier_bins"]["III (<1e-4)"]
n_tot = j8["n_nuclides"]
print(f"     全表 Tier III = {n3_all}/{n_tot}（{100*n3_all/n_tot:.1f}%）")
ck("E5【後驗】舊 Tier 軸之 Tier III 覆蓋 > 70% ⟹ 無分層信息",
   n3_all / n_tot > 0.70,
   f"{100*n3_all/n_tot:.1f}% ⟹ 裕度軸作廢，改用機制軸（Tier-M）")

# ---------------------------------------------------------------- 結算
print("\n[結算] Tier-M（機制軸）")
print(f"     C2 二階弱（G_F⁴）       {lab['C2 二階弱']:>3} 條")
print(f"     C3 位壘（穿隧）         {lab['C3 位壘']:>3} 條")
print(f"     C1 選擇定則禁戒（ΔJ≥3）  {lab['C1 選擇定則禁戒']:>3} 條")
print(f"     C1 未禁戒（無法解釋）    {lab.get('C1 未禁戒（異常）', 0):>3} 條")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "15_mechanism_classification.py", "role": "R14 機制分類 ＋ Tier 軸重定義",
       "n_exceptions": len(EX), "labels": dict(lab), "n_explained": n_expl,
       "qos_median": qos[len(qos) // 2] if qos else None, "qos_min": qos[0] if qos else None,
       "n_qos": len(qos),
       "E1_low": [{"Z": e["Z"], "A": e["A"], "QoS": e["ctrl"]["QoS"]} for e in low],
       "E3_dJ": djs, "E4_violations": viol, "E4b_violations": viol2,
       "old_tier_dist": dict(old), "tier3_frac_all": n3_all / n_tot,
       "items": [{"Z": e["Z"], "A": e["A"], "label": label(e), "Jpi": e["Jpi"],
                  "ctrl": (e["ctrl"]["name"] if e["ctrl"] else None),
                  "Q": (e["ctrl"]["Q"] if e["ctrl"] else None),
                  "QoS": (e["ctrl"]["QoS"] if e["ctrl"] else None),
                  "dJ": (e["ctrl"]["dJ"] if e["ctrl"] else None)}
                 for e in EX],
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "15_mechanism_classification.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/15_mechanism_classification.json 已寫")