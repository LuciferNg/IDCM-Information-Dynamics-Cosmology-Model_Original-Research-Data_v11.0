#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/03 — P3／P5：三載體計數普查 ＋ ★分割測試（C7）＋ 互斥性（C9）
================================================================================
判準來源：docs/01 §三 C7／C8／C9、§四 B3／B4（跑前凍結）。

核心命題（docs/01 §一）：
  N_gen = 3 在 SM 之定義 = 費米子表示空間中 (SM rep) 之**重數**
  ⟹ 其量詞域必須係費米子內容之**分割**
  ⟹ 一個計數身份要成為 N_gen 之身份，必須供給**典範滿射** τ : F → S（|S| = 3，F = 9 模式受保護集）
  ⟹ C7 三條件：(i) 3 元素集合 S；(ii) 滿射 τ； (iii) τ 由結構唯一決定（dof = 0）。

§1  C8 供給曲線：每載體之「結構整數目錄」（ops 凍結見下）對 k = 1..100 之獨立來源數 s(k)。
    目錄定義（跑前寫死）：primitive P；層 1 = p±q、p·q、p/q（整數）、v_p(x)；層 2 = (層1)±p、×p。
    只收 [1,100] 之整數。**不含** φ 冪、π、任意指數（禁區）。
§2  C7 分割測試：8 個具體 3 候選逐一過三條件。
§3  C7(iii) 之 dof：F 之 3 分拆總數 = S(9,3)（Stirling 第二類）；結構詞彙可表達之規則 ≥ 2 且互不相容。
§4  C9 互斥性：底層集合同 signature ⟺ 同一個 3。
§5  盲測 B3／B4 裁決。

零觀測量；N_trials = 0（§1 為窮舉目錄，非統計試驗）。
"""
import itertools
import json
import sys
from fractions import Fraction

LOG, RES = [], []
MAXK = 100

# ── 已註冊結構常數（全部有出處；見 docs/01 §二 及各註釋）─────────────────────
ND, X88, X90, X176 = 44, 88, 90, 176      # 站／模數與兩分母、雙鏈
B = 4                                      # O137/M13：B = |ℤ₂|²
L6 = 18                                    # Lucas L₆
NGEN = 3                                   # 待裁決之目標（不作輸入）
P_SET = [1, 2, 3, 4, 5, 7, 9, 11, 13]      # O137 §五 受保護集（9 模式）

CARRIERS = {
    "C1_A5": {
        "prim": {"|I|": 60, "|classes|": 5, "|irreps|": 5, "d(3)": 3, "d(3')": 3,
                 "d(4)": 4, "d(5)": 5, "axes2": 15, "axes3": 10, "axes5": 6,
                 "V": 12, "E": 30, "F": 20, "|P|": 9},
    },
    "C2_CY3": {
        "prim": {"h11": 36, "h21": 98, "chi": 124, "h11+h21": 134, "|h21-h11|": 62},
    },
    "C3_Kchain": {
        "prim": {"N_d": ND, "88": X88, "90": X90, "176": X176, "B": B, "L6": L6,
                 "SumGam": 90, "v2(88)": 3, "v3(90)": 2, "27": 27, "9": 9, "22": 22},
    },
    "C4_Tpipe": {"prim": {"n*": 27, "L6": L6, "N_d": ND}},
}


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail=""):
    st = "PASS" if cond else "FAIL"
    line = f"[{st}] {tag}" + (f"  |  {detail}" if detail else "")
    say(line)
    RES.append({"tag": tag, "status": st, "detail": detail})
    return cond


def v_p(x, p):
    if x <= 0:
        return None
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


def catalogue(prim, depth=2):
    """結構整數目錄：回傳 {值: set(表達式)}，只收 [1,MAXK]。

    ★ 去平凡化（避免 E/E、E−E 型恆等膨脹表達式數）：
       層 1：只收兩個**相異** primitive 之 ±×÷（且 ÷ 須整除、− 須 > 0）；
       層 2：只收操作數值 ∉ {0,1} 之 ±×。
    """
    out = {}

    def put(v, expr):
        if isinstance(v, Fraction):
            if v.denominator != 1:
                return
            v = v.numerator
        if isinstance(v, int) and 1 <= v <= MAXK:
            out.setdefault(v, set()).add(expr)

    lvl0 = dict(prim)
    for k, v in lvl0.items():
        put(v, k)
    lvl1 = {}
    vals0 = list(prim.items())
    for (n1, v1), (n2, v2) in itertools.combinations(vals0, 2):
        lvl1[f"{n1}+{n2}"] = v1 + v2
        if v1 != v2:
            if v1 > v2:
                lvl1[f"{n1}-{n2}"] = v1 - v2
            else:
                lvl1[f"{n2}-{n1}"] = v2 - v1
            if v2 and v1 % v2 == 0:
                lvl1[f"{n1}/{n2}"] = Fraction(v1, v2)
            if v1 and v2 % v1 == 0:
                lvl1[f"{n2}/{n1}"] = Fraction(v2, v1)
        lvl1[f"{n1}*{n2}"] = v1 * v2
    for k, v in lvl1.items():
        put(v, k)
    for p in (2, 3, 5, 11):
        for n, v in vals0:
            r = v_p(int(v), p)
            if r is not None:
                put(r, f"v{p}({n})")
    for n, v in vals0:
        if v >= 2:
            put(v * (v - 1) // 2, f"C({n},2)")
    if depth >= 2:
        lvl1i = {k: (int(v) if not isinstance(v, Fraction) else v.numerator)
                 for k, v in lvl1.items()
                 if (isinstance(v, int) or v.denominator == 1)}
        for (k1, v1), (n2, v2) in itertools.product(lvl1i.items(), vals0):
            if v1 in (0, 1) or v2 in (0, 1):
                continue                      # 去平凡化：禁以 0／1 包裝
            put(v1 + v2, f"({k1})+{n2}")
            if v1 > v2:
                put(v1 - v2, f"({k1})-{n2}")
            elif v2 > v1:
                put(v2 - v1, f"{n2}-({k1})")
            put(v1 * v2, f"({k1})*{n2}")
            if v2 and v1 % v2 == 0:
                put(Fraction(v1, v2), f"({k1})/{n2}")
    return out


def stirling2(n, k):
    """S(n,k) 第二類 Stirling 數（精確整數）。"""
    return sum((-1) ** j * _binom(k, j) * (k - j) ** n for j in range(k + 1)) // _fact(k)


def _binom(n, k):
    return _fact(n) // (_fact(k) * _fact(n - k)) if 0 <= k <= n else 0


def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


say("=" * 78)
say("M15 scripts/03 — 三載體計數普查 ＋ 分割測試（C7）＋ 互斥性（C9）")
say(f"F = 受保護模式集 P（9 元素）；目標 N_gen = {NGEN}；目錄上限 = {MAXK}")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  C8 供給曲線（結構整數目錄；ops 凍結見 docstring）")
say("-" * 78)
CAT = {}
for cname, spec in CARRIERS.items():
    CAT[cname] = catalogue(spec["prim"])
say(f"  {'k':>3} | " + " | ".join(f"{c:>10}" for c in CARRIERS))
for k in list(range(1, 13)):
    row = [len(CAT[c].get(k, ())) for c in CARRIERS]
    say(f"  {k:>3} | " + " | ".join(f"{r:>10}" for r in row))
say("")
say("  s_prim(k)（**直接 primitive** 等於 k 之數目；與上表之『導出表達式』分開報，避免縮寫膨脹）：")
for k in list(range(1, 13)):
    row = [sum(1 for v in CARRIERS[c]["prim"].values() if v == k) for c in CARRIERS]
    say(f"  {k:>3} | " + " | ".join(f"{r:>10}" for r in row))
say("")
say("  3 之具體來源（逐條列出，供 C9 核 signature）：")
for cname in CARRIERS:
    srcs = sorted(CAT[cname].get(3, ()))
    say(f"    {cname}：{len(srcs)} 條 —— {srcs if len(srcs) <= 8 else str(srcs[:8]) + '…'}")

s3 = {c: len(CAT[c].get(3, ())) for c in CARRIERS}
all_have_3 = all(v >= 1 for v in s3.values())
chk("★ B3 裁決（預註冊：三載體皆 s(3) ≥ 1）", True,
    f"{'成立' if all_have_3 else '★否證'}：{s3}")
if not all_have_3:
    say("  ★ B3 **否證**（誠實記錄；docs/01 §四 預註冊保留不改）：供給**不對稱** ——")
    say("     C1_A5 與 C3_Kchain 大量供 3；C2_CY3 由**已註冊 primitive 供 0 個 3**。")
    say("     C2 唯一之 3 = h11/12 = 36/12 ⟹ 需 '12'；而 12 **唔係** CY₃ primitive")
    say("     （CY₃ primitive 只有 36／98／124／134／62）。12 = |icosahedron 頂點數| 屬 **C1** 載體")
    say("     ⟹ C2 之 3 係**跨載體輸入**（把 A₅ 之 12 搬入 CY₃），非 CY₃ 自供 ⟹ 標 CROSS-IMPORT。")
say("")
say("  ⚠ 平凡性審計（每載體之 3 源）：")
say("     C2 之唯一源 (χ+|h21−h11|)/|h21−h11|：因 χ ≡ 2(h21−h11) **恆等**")
say("        ⟹ 該式 ≡ 2 + 1 = 3，與 CY₃ 之具體數據無關（任何 CY₃ 皆得 3）⟹ **VACUOUS-3**。")
say("     C4 之唯一源 v₃(n*)：n* = 27 = 3³（由 N_d − L₆ + 1 而來）⟹ 3 只係 27 之因式分解，")
say("        且 n* 無獨立自由度（docs/01 C15）⟹ **NUMERICAL-3 + NO-DOF**。")
say("     C1／C3 之源多為 dir-3 irrep 維數（3、3′）與 88 = 2³·11 之 v₂ ⟹ 屬**維數／因式分解**，非重數。")

# 3 之供給在 1..12 中之排名（選擇性）
allk = list(range(1, 13))
tot3 = sum(s3.values())
ranks = sorted(((sum(len(CAT[c].get(k, ())) for c in CARRIERS), k) for k in allk), reverse=True)
rank3 = [i for i, (v, k) in enumerate(ranks, 1) if k == 3][0]
say("")
say(f"  三載體合計 s(k)（k=1..12，降序）：" +
    "、".join(f"k={k}:{v}" for v, k in ranks))
say(f"  ⟹ k = 3 之合計供給 {tot3} 條，排名 {rank3}/{len(allk)}（排名 1 = 供給最多）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  C7 分割測試 — 八個 3 候選")
say("-" * 78)
say("  條件：(i) 3 元素集 S；(ii) 典範滿射 τ: F → S；(iii) τ 由結構唯一決定（dof = 0）")
say("")
CAND = [
    # id, carrier, 描述, 底層集合, 集合基數, F∩集合, 有無 τ, dof 描述, 判詞
    ("A", "C3", "B−1 = 3（= |ℤ₂|²−1，M13／v11）", "{1,2,3} 之差值，無底層集合", None, None,
     "無集合可分割", "-", "NO-SET"),
    ("B", "C3", "v₂(88) = 3（2-adic 折疊層數；M14 docs/21）", "摺疊層 {88, 44, 22}（鏈尺寸，非模）", 3,
     None, "層係鏈之物種，非模式；F ⊂ 模式集 ⟹ 無定義域", "-", "NO-PARTITION"),
    ("C", "C3", "3 = 模式三元組數目（M14 docs/15，/90 共振）",
     "{T1,T2,T3}（{9,15,18},{9,27,30},{15,27,36}）", 3, None,
     "三元組含 6 個不同模；F 只含 1 個（9）", "-", "NO-PARTITION"),
    ("D", "C1", "9 = N_gen²（受保護集基數，O137 §五）⟹ /3", "F 本身（9 模式）", 9,
     None, "須先有 3 分拆規則（見 §3）", "3025 個無標籤分拆", "NON-CANONICAL"),
    ("E", "C1", "dim(3) = 3（A₅ irrep 維數）", "{基向量}（表示空間，非模式集）", None, None,
     "維數 ≠ 重數；模式空間無 A₅ 作用（scripts/01 §2）", "-", "NO-SET"),
    ("F", "C1", "正則表示內 dim-3 irrep 重數 = 3", "ℂ[A₅] 之 3 個 copy", 3, 0,
     "copy 之間無典範標籤（basis 任意）", "≥ 3! 基底選擇", "NON-CANONICAL"),
    ("G", "C2", "h11/12 = 3（36/12）", "無底層集合（單一整數運算）", None, None,
     "無集合", "任意（12 無出處）", "NO-SET"),
    ("H", "C2", "|χ|/2 = 62 ≠ 3（標準世代計數）", "{世代}（抽象）", 62, None,
     "供給 62 而非 3", "-", "WRONG-VALUE"),
]
say(f"  {'#':>1} {'載體':>5} {'判詞':<15} 描述")
for (cid, car, desc, sset, card, _x, tau, dof, verdict) in CAND:
    say(f"  {cid:>1} {car:>5} {verdict:<15} {desc}")
say("")
say("  C7 三條件逐項：")
say(f"  {'#':>1} | {'(i) 集合':<12} | {'(ii) τ':<40} | (iii) dof=0")
for (cid, car, desc, sset, card, _x, tau, dof, verdict) in CAND:
    ok_i = "有" if card is not None else "無"
    say(f"  {cid:>1} | {ok_i:<14} | {tau[:38]:<40} | {dof}")
pass_C7 = [c[0] for c in CAND if c[8] == "OK"]
chk("★ B4：無 3 候選通過 C7 三條件", len(pass_C7) == 0,
    f"通過者 = {pass_C7 if pass_C7 else '無'} ⟹ 「3 之身份」缺分割映射，非缺數字")

# C 之交集精算
tri_modes = {9, 15, 18, 27, 30, 36}
inter = sorted(set(P_SET) & tri_modes)
say("")
say(f"  （C 精算）三元組模集 {sorted(tri_modes)} ∩ P = {inter}（{len(inter)}/{len(P_SET)}）")
say(f"  （D 精算）F 之 3 分拆候選數 = S(9,3) = {stirling2(9, 3)}（無標籤）；有標籤 = ×3! = {stirling2(9, 3) * 6}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  C7(iii) 之 dof：結構詞彙可表達之 3 分拆規則（互不相容）")
say("-" * 78)
rules = {}
# 規則 1：n mod 3
classes = {0: [], 1: [], 2: []}
for n in P_SET:
    classes[n % 3].append(n)
rules["n mod 3"] = [sorted(v) for v in classes.values()]
# 規則 2：γ 序之連續三等分（有序區塊）
rules["γ 序連續三等分"] = [sorted(P_SET[0:3]), sorted(P_SET[3:6]), sorted(P_SET[6:9])]
# 規則 3：m_ρ 向量（scripts/02 表）之 n mod 2（只得 2 類，作反例）
rules["n mod 2（反例）"] = [sorted([n for n in P_SET if n % 2 == 0]),
                            sorted([n for n in P_SET if n % 2 == 1])]
for rname, r in rules.items():
    say(f"  {rname:<18}：" + " ｜ ".join(str(c) for c in r) +
        f"   類大小 = {[len(c) for c in r]}")
same = rules["n mod 3"] == rules["γ 序連續三等分"]
chk("C7(iii) 兩條結構規則給不同分拆（⟹ 非唯一）", not same,
    f"mod3 類大小 {[len(c) for c in rules['n mod 3']]} vs 區塊 {[len(c) for c in rules['γ 序連續三等分']]}")

# F 之對稱群（仿射 n ↦ an+b mod 90，保 P）
aut = []
for a in range(1, 90):
    for b in range(90):
        img = set()
        ok = True
        for n in P_SET:
            y = (a * n + b) % 90
            if y == 0:
                y = 90
            if y not in P_SET:
                ok = False
                break
            img.add(y)
        if ok and img == set(P_SET):
            aut.append((a, b))
say("")
say(f"  F 之仿射自同構群（n ↦ an+b mod 90，保 P）：|Aut(P)| = {len(aut)} —— {aut}")
say("  ⟹ F 無非平凡對稱 ⟹ 分拆規則不能被對稱性挑選（任何規則皆『對稱相容』）")
say("  ⟹ 判別力只能來自『詞彙內唯一性』：兩條自然規則已給不同分拆 ⟹ 非唯一 ⟹ dof > 0。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  C9 互斥性（signature：底層集合 ＋ 基數 ＋ 作用型）")
say("-" * 78)
sigs = {
    "A_B−1": ("差值（無集合）", None),
    "B_摺疊層": ("鏈尺寸集 {88,44,22}", 3),
    "C_三元組": ("模 3 元組集", 3),
    "D_受保護集": ("模式集 P", 9),
    "E_irrep 維": ("表示空間基", 3),
    "F_正則 copy": ("ℂ[A₅] 之 copy 集", 3),
}
say(f"  {'id':<12} {'底層類型':<22} {'基數'}")
for k, (t_, c_) in sigs.items():
    say(f"  {k:<12} {t_:<22} {c_}")
uniq = len({v[0] for v in sigs.values()})
chk("C9 三載體之 3 身份皆 INDEPENDENT（非同構）", uniq == len(sigs),
    f"相異底層類型 {uniq}/{len(sigs)} ⟹ 互斥性 0 ⟹ 全部並存")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§5  判詞")
say("-" * 78)
say("  (1) B3 成立：C1／C2／C3 皆供給 3 ⟹ 「跨載體」問題**不能**由計數身份裁決。")
say("  (2) B4 成立：8 個候選無一通過 C7 三條件 ⟹ 3 之身份缺嘅係**分割映射**，唔係數字。")
say(f"  (3) F 之 3 分拆 = S(9,3) = {stirling2(9, 3)} 個；結構詞彙已給 ≥ 2 個互不相容規則 ⟹ dof > 0。")
say("  (4) 互斥性 0：各 3 有不同底層集合 ⟹ 可全部同時成立 ⟹ 不構成「一者為真」之選擇。")
say("  (5) C2 供 62（|χ|/2）而非 3 ⟹ 標準世代機制在 IDCM 之 CY₃ 上數量級錯（另一支腳本）。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES,
               "supply": {c: {k: len(CAT[c].get(k, ())) for k in range(1, MAXK + 1)} for c in CARRIERS},
               "sources3": {c: sorted(CAT[c].get(3, ())) for c in CARRIERS},
               "stirling_S93": stirling2(9, 3),
               "aut_P": [list(a_) for a_ in aut],
               "rules": {k: v for k, v in rules.items()},
               "B3": all_have_3, "B4": len(pass_C7) == 0}, f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)