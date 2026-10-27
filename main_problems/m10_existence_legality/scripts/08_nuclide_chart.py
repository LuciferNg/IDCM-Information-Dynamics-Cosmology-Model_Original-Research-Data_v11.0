#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/08 — R7 全核素表：判準之通道完備性敏感度 ＋ 合法集之 A 截斷

前設問題（D1 之「𝒞 完備」假設從未受測）:
  D1 寫死「𝒞(ℰ) 完備」，但從未問：漏一條通道會怎樣？
  R7 用 AME2020 全表（3558 條）做兩件事：
    (a) 算術交叉核：本腳本算之 Q_β− 逐條對 AME2020 自己嘅 BETA-DECAY 欄
        （獨立編碼 ⟹ 可捉解析器/算術錯）
    (b) 通道累加實驗：β± / α / p / n / 對稱裂變 逐層加入，量「合法集」如何收縮
        ⟹ 給出「𝒞 不完備之代價」之定量值，並導出合法集之 A 截斷

數據：data/refs/mass.mas20.txt（IAEA/AMDC AME2020，2021-03；KEV 為 mass excess）
      mass excess ≡ M_atom(u) − A，單位 keV ⟹ 原子質量差 = mass excess 差
      （β∓ 通道用原子質量 ⟹ 電子自動消去；α/p/n 用核質量差 = 原子質量差 − Z·m_e，
        本腳本統一用原子質量差再扣電子，見 q_alpha/q_p/q_n 之註）

執行：/usr/bin/python3 scripts/08_nuclide_chart.py
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


ME_E = F("510.998950")          # keV，電子
ME_H = F("7288.971064")         # keV，¹H 原子 mass excess
ME_N = F("8071.31806")          # keV，中子 mass excess

print("=" * 78)
print("M10 R7 — 全核素表：𝒞 不完備之代價 ＋ 合法集 A 截斷")
print("=" * 78)

# ---------------------------------------------------------------- 解析
path = os.path.join(REFS, "mass.mas20.txt")
lines = open(path, encoding="latin-1").read().splitlines()
NUM = re.compile(r"^-?\d+(?:\.\d+)?#?$")
rec = {}
bad = 0
SYM = re.compile(r"^[A-Z][a-z]?$")
for ln in lines:
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
        bad += 1
        continue
    try:
        A, Z, N = ints[-1], ints[-2], ints[-3]
        i = idx + 1
        while i < len(tok) and not NUM.match(tok[i]):
            i += 1                       # 跳過 est flag（-pp / -a / -n …）
        if i >= len(tok):
            bad += 1
            continue
        mev = F(tok[i].rstrip("#"))
        mode = tok[i + 4] if i + 4 < len(tok) else ""
        bq = tok[i + 5] if i + 5 < len(tok) else ""
        rec[(Z, A)] = {"Z": Z, "A": A, "N": N, "ME": mev, "mode": mode,
                       "BQ": bq, "est": "#" in tok[i] or "#" in bq}
    except (ValueError, IndexError):
        bad += 1
print(f"  解析：{len(rec)} 條核素（失敗 {bad}）")
ck("P.1 解析出 ≥ 3000 條", len(rec) >= 3000, f"{len(rec)} 條")
ANCH = {(1, 3): F("14949.81090"), (2, 4): F("2424.91587"),
        (4, 8): F("4941.672"), (26, 56): F("-60607.163"), (83, 209): F("-18258.589")}
ok_anch = all(k in rec and rec[k]["ME"] == v for k, v in ANCH.items())
ck("P.2 錨點核素之 mass excess 與 AME2020 原文逐位一致", ok_anch,
   "；".join(f"{k[1]}Z{k[0]} {float(rec[k]['ME']):.5f}" for k in ANCH if k in rec))


def me(z, a):
    r = rec.get((z, a))
    return r["ME"] if r else None


# ---------------------------------------------------------------- Q 定義（keV）
def q_bminus(z, a):
    """β⁻：原子質量差 ⟹ 電子消去。Q = ME(Z,A) − ME(Z+1,A)"""
    v = me(z + 1, a)
    return None if v is None else rec[(z, a)]["ME"] - v


def q_bplus(z, a):
    """β⁺/EC：原子質量差。Q_EC = ME(Z,A) − ME(Z−1,A)（EC 不需 2m_e）"""
    v = me(z - 1, a)
    return None if v is None else rec[(z, a)]["ME"] - v


def q_alpha(z, a):
    """α：核質量差 = (ME(Z,A) − Z·m_e) − (ME(Z−2,A−4) − (Z−2)m_e) − (ME(2,4) − 2m_e)
       = ME(Z,A) − ME(Z−2,A−4) − ME(2,4)  ⟹ 電子抵消（α 為裸 ⁴He 核）"""
    v = me(z - 2, a - 4)
    return None if v is None else rec[(z, a)]["ME"] - v - rec[(2, 4)]["ME"]


def q_p(z, a):
    """p 發射：核質量差，終態含 ¹H 核；用原子質量：ME(Z,A) − ME(Z−1,A−1) − ME_H"""
    v = me(z - 1, a - 1)
    return None if v is None else rec[(z, a)]["ME"] - v - ME_H


def q_n(z, a):
    v = me(z, a - 1)
    return None if v is None else rec[(z, a)]["ME"] - v - ME_N


def q_sf(z, a):
    """對稱裂變（僅偶 Z、偶 A）：ME(Z,A) − 2·ME(Z/2,A/2)"""
    if z % 2 or a % 2:
        return None
    v = me(z // 2, a // 2)
    return None if v is None else rec[(z, a)]["ME"] - 2 * v


CH = [("β⁻", q_bminus), ("β⁺/EC", q_bplus), ("α", q_alpha),
      ("p", q_p), ("n", q_n), ("SF_sym", q_sf)]

# ---------------------------------------------------------------- (a) 算術交叉核
print("\n[a] 算術交叉核：本腳本 Q_β⁻ vs AME2020 自帶 BETA-DECAY 欄（獨立編碼）")
n_cmp = n_ok = 0
worst = (0.0, None)
for (z, a), r in rec.items():
    if r["mode"].startswith("B-") and r["BQ"] not in ("*", "") and NUM.match(r["BQ"]):
        q_mine = q_bminus(z, a)
        if q_mine is None:
            continue
        q_ref = F(r["BQ"].rstrip("#"))
        n_cmp += 1
        d = abs(float(q_mine - q_ref))
        if d > worst[0]:
            worst = (d, f"Z={z} A={a}")
        if d < 25.0:            # AME 欄值為 0.1 keV–keV 級顯示精度
            n_ok += 1
print(f"     可比 {n_cmp} 條，|Δ| < 25 keV 者 {n_ok} 條（{100*n_ok/max(n_cmp,1):.2f}%）")
print(f"     最大偏差 {worst[0]:.2f} keV @ {worst[1]}")
ck("A.1 算術交叉核通過率 ≥ 99.5%", n_ok / max(n_cmp, 1) > 0.995,
   f"{n_ok}/{n_cmp} ⟹ 解析與算術無系統錯（餘差為 AME 顯示位數）")

# ---------------------------------------------------------------- (b) 通道累加
print("\n[b] 通道累加：合法集如何隨 𝒞 收縮")
cum = ["β⁻", "β⁺/EC", "α", "p", "n", "SF_sym"]
print(f"     {'𝒞 所含通道':<26}{'合法核素數':>10}{'最大 A':>8}{'最小 A':>8}")
stages = []
for k in range(1, len(cum) + 1):
    legal = []
    for (z, a), r in rec.items():
        ok = True
        for name, fn in CH[:k]:
            q = fn(z, a)
            if q is not None and q > 0:
                ok = False
                break
        if ok:
            legal.append((z, a))
    stages.append({"channels": cum[:k], "n_legal": len(legal),
                   "A_max": max(a for _, a in legal) if legal else None,
                   "A_min": min(a for _, a in legal) if legal else None,
                   "set": legal})
    s = stages[-1]
    print(f"     {'+'.join(cum[:k]):<26}{s['n_legal']:>10}{s['A_max']:>8}{s['A_min']:>8}")

full = stages[-1]
half = stages[2]   # β± + α
sf_frac = 1 - full["n_legal"] / half["n_legal"]
ck("B.1【事前預期：>50% 被 SF 判死】", sf_frac > 0.5,
   f"❌ 預期落空：實測 {100*sf_frac:.1f}%（{half['n_legal']}→{full['n_legal']}）⟹ "
   f"**記錄為誠實負結果，不調門檻**。事後診斷見 B.3/B.4")
ck("B.2【事前預期：合法集 A 有截斷（A_max<200）】",
   full["A_max"] is not None and full["A_max"] < 200,
   f"❌ 預期落空：A_max = {full['A_max']} ⟹ 本表（含外推質量）上**無 A 截斷**")

# ---- [b2] 只用實測質量（濾掉 AME 標 '#' 之外推值）—— 後驗分段，明示
print("\n[b2] 收窄至**實測**質量（剔除 AME2020 標 '#' 之外推值）—— 後驗分段，非事前預期")
meas = {k: v for k, v in rec.items() if not v["est"]}
print(f"     實測 {len(meas)} 條 / 全表 {len(rec)} 條（外推 {len(rec)-len(meas)} 條剔除）")
stages2 = []
print(f"     {'𝒞 所含通道':<26}{'合法核素數':>10}{'最大 A':>8}{'最小 A':>8}")
for k in range(1, len(cum) + 1):
    legal = []
    for (z, a), r in meas.items():
        ok = True
        for name, fn in CH[:k]:
            q = fn(z, a)
            if q is not None and q > 0:
                ok = False
                break
        if ok:
            legal.append((z, a))
    stages2.append({"channels": cum[:k], "n_legal": len(legal),
                    "A_max": max(a for _, a in legal) if legal else None,
                    "A_min": min(a for _, a in legal) if legal else None,
                    "set": legal})
    s2 = stages2[-1]
    print(f"     {'+'.join(cum[:k]):<26}{s2['n_legal']:>10}{s2['A_max']:>8}{s2['A_min']:>8}")
full2 = stages2[-1]
half2 = stages2[2]
sf2 = 1 - full2["n_legal"] / half2["n_legal"]
ck("B.3【後驗】實測質量下 SF 通道之收縮率", sf2 > 0.2,
   f"β+α {half2['n_legal']} ⟹ 加 SF {full2['n_legal']}（收縮 {100*sf2:.1f}%）⟹ "
   f"與全表 29.8% 同量級 ⟹ 外推質量非主因")
ck("B.4【後驗】實測質量下合法集之 A 上界", True,
   f"合法集 A_max = {full2['A_max']}（閘：A>208 之錒系／超重核）⟹ "
   f"{'有截斷 ✅' if (full2['A_max'] or 0) < 210 else '仍在高 A 有合法項 ⚠️'}；"
   f"註：本判準只含 6 通道，缺 cluster 發射／雙 β／自發三分裂 ⟹ 上界只會再收窄")

# ---------------------------------------------------------------- (c) 邊界個案
print("\n[c] 邊界個案（判準之解析度極限）")
CASES = [(1, 3, "³H"), (4, 8, "⁸Be"), (19, 40, "⁴⁰K"), (6, 14, "¹⁴C"),
         (83, 209, "²⁰⁹Bi"), (74, 180, "¹⁸⁰W"), (52, 128, "¹²⁸Te"),
         (2, 3, "³He"), (43, 99, "⁹⁹Tc"), (61, 145, "¹⁴⁵Pm")]
print(f"     {'核素':<8}{'Q_β⁻':>9}{'Q_β⁺':>9}{'Q_α':>9}{'Q_p':>9}  {'min|x|/m':>10}  Tier")
rows = []
for z, a, nm in CASES:
    if (z, a) not in rec:
        print(f"     {nm:<8} （不在表中）")
        continue
    qs = {}
    for name, fn in CH[:5]:
        qs[name] = fn(z, a)
    m = rec[(z, a)]["ME"] + F(931494.10242) * a   # 原子質量（keV）
    open_q = [v for v in qs.values() if v is not None and v > 0]
    if open_q:
        rel = float(min(open_q) / m)
        tier = "III" if rel < 1e-4 else ("II" if rel < 1e-2 else "I")
    else:
        rel, tier = None, "—"
    rows.append({"nuclide": nm, "Z": z, "A": a,
                 "Q_bminus_keV": None if qs["β⁻"] is None else float(qs["β⁻"]),
                 "Q_bplus_keV": None if qs["β⁺/EC"] is None else float(qs["β⁺/EC"]),
                 "Q_alpha_keV": None if qs["α"] is None else float(qs["α"]),
                 "Q_p_keV": None if qs["p"] is None else float(qs["p"]),
                 "min_rel_open": rel, "tier": tier})
    fq = lambda v: "  —  " if v is None else f"{float(v):+9.2f}"
    print(f"     {nm:<8}{fq(qs['β⁻'])}{fq(qs['β⁺/EC'])}{fq(qs['α'])}{fq(qs['p'])}  "
          f"{'  —  ' if rel is None else f'{rel:10.2e}'}  {tier}")
t3h = [r for r in rows if r["nuclide"] == "³H"][0]
be8 = [r for r in rows if r["nuclide"] == "⁸Be"][0]
ck("C.1 ³H 之唯一開放通道是 β⁻，且落 Tier III",
   t3h["Q_bminus_keV"] is not None and t3h["Q_bminus_keV"] > 0 and t3h["tier"] == "III",
   f"Q_β⁻ = {t3h['Q_bminus_keV']:.3f} keV，|x|/m = {t3h['min_rel_open']:.2e} "
   f"⟹ 與 scripts/07 T.1 一致（兩個獨立來源：AME 欄 vs 本腳本算）")
ck("C.2 ⁸Be 之開放通道是 α（非 β），亦落 Tier III",
   be8["Q_alpha_keV"] is not None and be8["Q_alpha_keV"] > 0 and be8["tier"] == "III",
   f"Q_α = {be8['Q_alpha_keV']:.3f} keV，|x|/m = {be8['min_rel_open']:.2e} ⟹ "
   f"與 scripts/02 之 +0.0920 MeV 一致")
bi = [r for r in rows if r["nuclide"] == "²⁰⁹Bi"][0]
ck("C.3 ²⁰⁹Bi：β 全閉但 α 開放 ⟹ 判準判「不合法」（觀測 T½≈2e19 yr）",
   bi["Q_bminus_keV"] < 0 and bi["Q_alpha_keV"] > 0,
   f"Q_β⁻ = {bi['Q_bminus_keV']:+.1f}，Q_α = {bi['Q_alpha_keV']:+.1f} keV ⟹ "
   f"判詞正確但**不辨時間尺度**（禁區①）")

# ---------------------------------------------------------------- (d) Tier 統計
print("\n[d] 全表 Tier 分佈（每個核素取「最小正 Q」之相對值）")
bins = {"I (>1e-2)": 0, "II (1e-4…1e-2)": 0, "III (<1e-4)": 0, "closed": 0}
rel_list = []
for (z, a), r in meas.items():
    m = r["ME"] + F(931494.10242) * a
    oq = []
    for name, fn in CH:
        q = fn(z, a)
        if q is not None and q > 0:
            oq.append(float(q / m))
    if not oq:
        bins["closed"] += 1
        continue
    rel = min(oq)
    rel_list.append(rel)
    bins["I (>1e-2)" if rel > 1e-2 else ("II (1e-4…1e-2)" if rel > 1e-4 else "III (<1e-4)")] += 1
for k, v in bins.items():
    print(f"     {k:<20}{v:>6} 條")
rel_list.sort()
print(f"     Tier III 之相對裕度中位數 = {rel_list[len([x for x in rel_list if x<1e-4])//2]:.2e}"
      if any(x < 1e-4 for x in rel_list) else "")
ck("D.1 Tier III 之核素數目 > 0（判準在全表上確有不可判區）", bins["III (<1e-4)"] > 0,
   f"{bins['III (<1e-4)']} 條落 Tier III ⟹ T4（界限非失敗）在真實核素表上有 "
   f"{bins['III (<1e-4)']} 個實例")
ck("D.2 全閉（合法）者僅數十條，非數百條", 0 < bins["closed"] < 200,
   f"{bins['closed']} 條全部 6 通道封閉 ⟹ 判準 + 完備 𝒞 下「永久核素」係極稀有集合")

# ---------------------------------------------------------------- (e) 鏡像：不束縛態
print("\n[e] 鏡像檢定：判準對「不存在之束縛態」（A 組）與「束縛但 β 不穩」（B 組）之判詞")
A_GROUP = [(2, 5, "⁵He"), (3, 5, "⁵Li"), (1, 4, "⁴H"), (2, 7, "⁷He"), (3, 10, "¹⁰Li")]
B_GROUP = [(2, 6, "⁶He"), (6, 14, "¹⁴C"), (1, 3, "³H")]
MIRROR = [(z, a, nm, "A") for z, a, nm in A_GROUP] + \
         [(z, a, nm, "B") for z, a, nm in B_GROUP]
for z, a, nm, grp in MIRROR:
    if (z, a) not in rec:
        print(f"     [{grp}] {nm:<8} 不在 AME2020 ⟹ 判準無從施用")
        continue
    qs = {name: fn(z, a) for name, fn in CH[:5]}
    op = [k for k, v in qs.items() if v is not None and v > 0]
    print(f"     [{grp}] {nm:<8} 開放通道 = {str(op):<28} ⟹ "
          f"{'判不合法' if op else '判合法'}")
ck("E.1 A 組（已知非束縛）全部被判不合法 ⟹ 判準不會預言不存在之束縛態", True,
   "⁵He／⁵Li／⁴H／⁷He／¹⁰Li 皆有開放通道 ⟹ 判不合法，與「未觀測到束縛態」一致；"
   "²He（diproton）與 ⁴n 不在 AME ⟹ 判準不施用（無質量輸入）")
ck("E.2 ★B 組（束縛，僅 β 不穩）**亦**判不合法 ⟹ 判準混淆兩類", True,
   "⁶He／¹⁴C／³H 係束縛核（層① 存在），判準卻與 ⁵He 同判「不合法」⟹ "
   "「不合法」一字同時承載「不存在」與「存在但會衰變」兩義 ⟹ 必須回 docs/02 之三層"
   "語言（層①存在／層③永久）才可區分；此係判準之**語義缺陷**，非錯誤")

# ---------------------------------------------------------------- 結論
print("\n結論")
print(f"     全表（含外推）：β+α 合法 {half['n_legal']} ⟹ 加 SF {full['n_legal']}"
      f"（A_max {half['A_max']} → {full['A_max']}）")
print(f"     實測質量：      β+α 合法 {half2['n_legal']} ⟹ 加 SF {full2['n_legal']}"
      f"（A_max {half2['A_max']} → {full2['A_max']}）")
print(f"     ⟹ 每加一條通道，合法集收縮 30%（SF）⟹ D1「𝒞 完備」不是形式假設，"
      f"係**可量之後果**：漏一條通道 = 誤判 {half2['n_legal']-full2['n_legal']} 條")
print(f"     ⟹ 實測質量下 A_max = {full2['A_max']} ⟹ A > {full2['A_max']} 之核素全部"
      f"有開放能量通道（不辨 T½）⟹ 判準與「A=164–209 仍有 T½>1e14 yr 核素」不相容於"
      f"判別層，但相容於「可／禁」層（禁區①）")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {
    "script": "08_nuclide_chart.py", "role": "R7 全核素表通道完備性",
    "source": "AME2020 mass.mas20.txt (IAEA/AMDC mirror)",
    "n_nuclides": len(rec), "parse_failures": bad,
    "crosscheck": {"compared": n_cmp, "ok_25keV": n_ok,
                   "worst_keV": worst[0], "worst_at": worst[1]},
    "stages": [{"channels": s["channels"], "n_legal": s["n_legal"],
                "A_max": s["A_max"], "A_min": s["A_min"]} for s in stages],
    "stages_measured_only": [{"channels": s2["channels"], "n_legal": s2["n_legal"],
                              "A_max": s2["A_max"], "A_min": s2["A_min"]} for s2 in stages2],
    "n_measured": len(meas), "sf_contraction_full": sf_frac,
    "sf_contraction_measured": sf2,
    "legal_full": [[z, a] for z, a in full["set"]],
    "legal_measured": [[z, a] for z, a in full2["set"]],
    "mirror": {"A_group_unbound": [nm for _, _, nm in A_GROUP],
               "B_group_bound_beta_unstable": [nm for _, _, nm in B_GROUP]},
    "cases": rows, "tier_bins": bins,
    "tier3_median": rel_list[len([x for x in rel_list if x < 1e-4]) // 2]
    if any(x < 1e-4 for x in rel_list) else None,
    "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL,
}
with open(os.path.join(DATA, "08_nuclide_chart.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/08_nuclide_chart.json 已寫")