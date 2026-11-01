#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/05 — C2（CY₃）：標準世代計數之判死 ＋ 跨載體輸入審計
================================================================================
判準來源：docs/01 §三 C12／C13（跑前凍結）。

§1  C12：CY₃ 之唯一已註冊世代計數機制 = N_gen = |χ|/2（χ = 2(h11−h21)；heterotic／F-theory 標準）。
     IDCM 註冊值 h11 = 36、h21 = 98（M12）⟹ χ = −124 ⟹ |χ|/2 = 62。
     判準：|62 − 3|/3 > 0.1 ⟹ 數量級錯 ⟹ 機制判死。
§2  C13 密度閘：由 {h11,h21,χ} 經 {+,−,×,÷, C(·,2), v_p} 得 3 之表達式；並報目錄對 [1,100] 覆蓋率。
§3  跨載體輸入審計：C2 之 3 只在 h11/12 出現，需 '12'；12 ∉ CY₃ primitive ⟹ 需由 C1（icosahedron
     頂點數 = 12）輸入 ⟹ 標 CROSS-IMPORT（違 A1 出處閘）。
§4  非對稱核：3 | h11（36 = 4·9）但 3 ∤ h21（98 = 2·7²）⟹ 「3」只係單一 Hodge 數之因式，
     非任何集合之基數。
§5  判詞。

零觀測量；N_trials = 0（§2 為窮舉目錄）。
"""
import itertools
import json
import sys
from fractions import Fraction

LOG, RES = [], []
H11, H21 = 36, 98
CHI = 2 * (H11 - H21)          # = −124
MAXK = 100


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


say("=" * 78)
say("M15 scripts/05 — C2（CY₃）：標準世代計數判死 ＋ 跨載體輸入審計")
say("CY₃ 註冊數據（M12）：h11 = 36、h21 = 98、χ = 2(h11−h21)")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  C12 標準世代計數")
say("-" * 78)
say(f"  χ = 2({H11} − {H21}) = {CHI}")
say(f"  |χ|/2 = {abs(CHI) // 2}")
ng = abs(CHI) // 2
rel = abs(ng - 3) / 3
say(f"  與 N_gen = 3 之相對偏差 = |62 − 3|/3 = {rel:.4f}（{(ng / 3):.3f} 倍）")
say(f"  判準門檻 = 0.1 ⟹ {'超標 ⟹ 判死' if rel > 0.1 else '未超標'}")
chk("C12 標準機制 |χ|/2 = 62 ≠ 3（數量級錯）", rel > 0.1 and ng != 3,
    f"ng = {ng}；相對偏差 {rel:.3f}（{ng / 3:.2f}×）")
say("  ⚠ 範圍：此判死針對『CY₃ 之 χ/2 世代計數』此一**機制**；不宣稱 CY₃ 整體不可能供 3")
say("     （後者見 §3：唯一之 3 需跨載體輸入）。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  C13 密度閘：由 {h11, h21, |χ|, |h21−h11|} 得 3 之表達式")
say("-" * 78)
PRIM = {"h11": H11, "h21": H21, "abs_chi": abs(CHI), "abs_diff": abs(H21 - H11)}
cat = set()
expr3 = []


def note(v, expr):
    if isinstance(v, Fraction):
        if v.denominator != 1:
            return
        v = v.numerator
    if isinstance(v, int) and 1 <= v <= MAXK:
        cat.add(v)
        if v == 3:
            expr3.append(expr)


for k, v in PRIM.items():
    note(v, k)
for (n1, v1), (n2, v2) in itertools.combinations(list(PRIM.items()), 2):
    note(v1 + v2, f"{n1}+{n2}")
    if v1 != v2:
        note(abs(v1 - v2), f"|{n1}−{n2}|")
        if v1 % v2 == 0:
            note(Fraction(v1, v2), f"{n1}/{n2}")
        if v2 % v1 == 0:
            note(Fraction(v2, v1), f"{n2}/{n1}")
    note(v1 * v2, f"{n1}*{n2}")
for k, v in PRIM.items():
    note(v * (v - 1) // 2, f"C({k},2)")
    for p in (2, 3, 5, 7):
        r = v_p(v, p)
        if r is not None:
            note(r, f"v{p}({k})")
# 層 2（與 scripts/03 之去平凡化目錄一致：操作數值 ∉ {0,1}）
lvl1 = {}
for (n1, v1), (n2, v2) in itertools.combinations(list(PRIM.items()), 2):
    lvl1[f"{n1}+{n2}"] = v1 + v2
    if v1 != v2:
        lvl1[f"|{n1}−{n2}|"] = abs(v1 - v2)
    if v2 and v1 % v2 == 0:
        lvl1[f"{n1}/{n2}"] = Fraction(v1, v2)
    if v1 and v2 % v1 == 0:
        lvl1[f"{n2}/{n1}"] = Fraction(v2, v1)
    lvl1[f"{n1}*{n2}"] = v1 * v2
for (k1, v1), (n2, v2) in itertools.product(lvl1.items(), list(PRIM.items())):
    v1i = v1.numerator if isinstance(v1, Fraction) else v1
    if v1i in (0, 1) or v2 in (0, 1):
        continue
    note(v1i + v2, f"({k1})+{n2}")
    if v1i > v2:
        note(v1i - v2, f"({k1})-{n2}")
    elif v2 > v1i:
        note(v2 - v1i, f"{n2}-({k1})")
    note(v1i * v2, f"({k1})*{n2}")
    if v2 and v1i % v2 == 0:
        note(Fraction(v1i, v2), f"({k1})/{n2}")
    if v1i and v2 % v1i == 0:
        note(Fraction(v2, v1i), f"{n2}/({k1})")
say(f"  目錄（primitive {len(PRIM)} 個，ops 見 docstring）：覆蓋 [1,{MAXK}] 中 {len(cat)} 個整數"
    f"（覆蓋率 {len(cat) / MAXK:.1%}）")
say(f"  得 3 之表達式：{expr3 if expr3 else '（無）'}")
say(f"  其中含 '12' 者：{[e for e in expr3 if '12' in e]}")
if expr3:
    say("  ⚠ 平凡性核：abs_chi = 2·abs_diff **恆等**（χ = 2(h11−h21)）⟹ 該式 ≡ 2+1 = 3,")
    say("     與 CY₃ 之具體數據無關（任何 CY₃ 皆給 3）⟹ **VACUOUS-3**，不作證據。")
chk("C13 覆蓋率記錄（密度閘）", True,
    f"覆蓋率 {len(cat) / MAXK:.1%}；得 3 之式 = {len(expr3)}")
if expr3:
    say("  ⚠ 每一式皆須通過 A1 出處閘：catalogue 內之常數必須係 CY₃ 之**已註冊量**。")
    say("     → 下一節核 12 之出處。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  跨載體輸入審計：C2 之 3 = h11/12 需要 '12'")
say("-" * 78)
say(f"  h11/12 = {H11}/12 = {H11 // 12} ⟹ 需常數 12。")
say("  CY₃ 之已註冊 primitive：h11 = 36、h21 = 98、χ = −124（衍生：|χ| = 124、|h21−h11| = 62）")
say("  ⟹ **12 唔在列**；且 12 由 {36, 98, 124, 62} 之 ±×÷ 亦不可得（見 §2 目錄：12 ∈ catalogue？）")
say(f"  核查：12 ∈ 目錄？{'是' if 12 in cat else '否'}；若『是』，列出其表達式：")
if 12 in cat:
    # 反向列舉 12 之來源
    srcs = []
    for k, v in PRIM.items():
        if v == 12:
            srcs.append(k)
    for (n1, v1), (n2, v2) in itertools.combinations(list(PRIM.items()), 2):
        if v1 + v2 == 12:
            srcs.append(f"{n1}+{n2}")
        if abs(v1 - v2) == 12:
            srcs.append(f"|{n1}−{n2}|")
        if v2 and v1 % v2 == 0 and v1 // v2 == 12:
            srcs.append(f"{n1}/{n2}")
        if v1 and v2 % v1 == 0 and v2 // v1 == 12:
            srcs.append(f"{n2}/{n1}")
        if v1 * v2 == 12:
            srcs.append(f"{n1}*{n2}")
    say(f"    {srcs if srcs else '（無）'}")
say("  C1（A₅／icosahedron）之 primitive 含 V = 12（頂點數）⟹ 12 屬 **C1** 載體。")
say("  ⟹ 3 = h11 / |V_ico| = （CY₃ 量）／（A₅ 量）⟹ 兩者無已註冊映射 ⟹ **CROSS-IMPORT**（違 A1）。")
chk("C3 之 3 須跨載體輸入（12 ∉ CY₃ primitive 且不可由 CY₃ 量合成）",
    12 not in cat,
    "12 ∉ catalogue ⟹ 3 須由外部輸入 ⟹ CROSS-IMPORT")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  非對稱核：3 只係單一 Hodge 數之因式")
say("-" * 78)
say(f"  h11 = {H11} = 4·9 ⟹ 3 | h11，v₃(h11) = {v_p(H11, 3)}")
say(f"  h21 = {H21} = 2·7² ⟹ 3 ∤ h21，v₃(h21) = {v_p(H21, 3)}")
say(f"  χ = {CHI} ⟹ 3 ∤ |χ|（v₃ = {v_p(abs(CHI), 3)}）")
chk("§4 3 只整除 h11（不對稱）",
    v_p(H11, 3) >= 1 and (v_p(H21, 3) == 0) and (v_p(abs(CHI), 3) == 0),
    f"v₃(h11)={v_p(H11, 3)}, v₃(h21)={v_p(H21, 3)}, v₃(|χ|)={v_p(abs(CHI), 3)}")
say("  ⟹ 「3」若取 h11 之因式，則係**單一 Hodge 數之算術因式**，唔係任何集合之基數；")
say("     且 3 ∤ h21、3 ∤ χ ⟹ 冇任何 CY₃ 整體量攜帶 3。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§5  判詞")
say("-" * 78)
say("  (1) C12：|χ|/2 = 62，與 3 差 20.67 倍 ⟹ CY₃ 之標準世代計數機制**判死**（範圍：該機制）。")
say("  (2) C13：CY₃ 目錄覆蓋率低（見 §2），但唯一『3』須 12 ⟹ 12 ∉ CY₃ ⟹ CROSS-IMPORT（違 A1）。")
say("  (3) §4：3 只在 h11 之因式出現，3 ∤ h21、3 ∤ χ ⟹ 無整體量攜 3。")
say("  (4) ⟹ **C2（CY₃）不能供給 N_gen = 3**（範圍：已註冊 primitive 36／98；已註冊機制 χ/2）。")
say("  (5) 與 M12 S5（『無可算載體』）一致：本輪把 S5 由『阻塞』升為『scoped no-go』。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "h11": H11, "h21": H21, "chi": CHI, "ng_chi2": ng,
               "catalogue_coverage": len(cat), "expr3": expr3, "twelve_in_cat": 12 in cat},
              f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)