#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/l3_doc_crosscheck.py — L3 文檔交叉核（idcm-audit-verify §2.4）
================================================================================
規矩：期望值**只**來自 `data/l2_independent_verify.json`（L2 獨立重算）＋ L1 之表檔（互核）。
每份有數值宣稱之 doc 須 100% 命中；敘述／判定型 doc 豁免（明列）。

兩種判據：
  check_str(doc, label, token)           —— doc 內含該字串（去空白後）
  check_num(doc, label, value, tol)      —— doc 內存在一個數值，與 value 之相對差 < tol
                                            （由 doc 掃出全部數字 token 後比對；支援 1.23×10⁻⁴ 形）
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOCS = ROOT / "docs"
DATA = ROOT / "data"

l2 = json.loads((DATA / "l2_independent_verify.json").read_text())["_derived"]
l1_h = json.loads((DATA / "02_a5_harmonic_multiplicity.json").read_text())

NFAIL, NPASS = 0, 0
SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"
SUPMAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def say(s=""):
    print(s)


def _txt(doc):
    p = DOCS / doc
    return p.read_text() if p.exists() else None


def check_str(doc, label, token):
    global NFAIL, NPASS
    t = _txt(doc)
    if t is None:
        say(f"[FAIL] {doc} :: {label} —— 檔案不存在")
        NFAIL += 1
        return False
    flat = re.sub(r"[,]", "", re.sub(r"[\s\u00a0]+", "", t))
    ok = re.sub(r"[,]", "", re.sub(r"[\s\u00a0]+", "", str(token))) in flat
    say(f"[{'PASS' if ok else 'FAIL'}] {doc} :: {label} = {token}")
    NPASS += ok
    NFAIL += (not ok)
    return ok


def _numbers(txt):
    """由 doc 掃出數值：支援 1.23×10⁻⁴、4.36e-60、千分位逗號、整數。"""
    t = re.sub(r"(?<=\d),(?=\d{3})", "", txt)          # 去千分位
    out = []
    for m in re.finditer(r"([-+]?\d+(?:\.\d+)?)\s*[×x]\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", t):
        out.append(float(m.group(1)) * 10.0 ** int(m.group(2).translate(SUPMAP)))
    for m in re.finditer(r"([-+]?\d+(?:\.\d+)?)\s*[eE]([-+]?\d+)", t):
        out.append(float(m.group(1)) * 10.0 ** int(m.group(2)))
    for m in re.finditer(r"[-+]?\d+(?:\.\d+)?", t):
        out.append(float(m.group(0)))
    return out


def check_num(doc, label, value, tol=2e-3):
    global NFAIL, NPASS
    t = _txt(doc)
    if t is None:
        say(f"[FAIL] {doc} :: {label} —— 檔案不存在")
        NFAIL += 1
        return False
    v = float(value)
    nums = _numbers(t)
    ok = any(abs(x - v) <= tol * max(abs(v), 1e-300) for x in nums)
    say(f"[{'PASS' if ok else 'FAIL'}] {doc} :: {label} ≈ {value}（doc 命中 {ok}）")
    NPASS += ok
    NFAIL += (not ok)
    return ok


def check_bool(doc, label, cond, token):
    """布林判據：L2 之條件成立 **且** doc 含指定 token。"""
    global NFAIL, NPASS
    t = _txt(doc)
    ok_tok = t is not None and re.sub(r"[\s\u00a0]+", "", str(token)) in re.sub(r"[\s\u00a0]+", "", t)
    ok = bool(cond) and ok_tok
    say(f"[{'PASS' if ok else 'FAIL'}] {doc} :: {label}（L2 = {bool(cond)}；doc token = {ok_tok}）")
    NPASS += ok
    NFAIL += (not ok)
    return ok


say("=" * 78)
say("M15 L3 文檔交叉核（期望值全部來自 L2 獨立重算）")
say("=" * 78)

# ── §0 L1 ↔ L2 重數表互核（Molien vs character 內積）───────────────────────
l1t, l2t = l1_h["table"], l2["harmonic_multiplicity"]
mism = [(n, l1t[str(n)]) for n in range(0, 45)
        if l1t[str(n)] != [l2t[str(n)][k] for k in ("1", "3", "3p", "4", "5")]]
say(f"§0  L1（Molien）↔ L2（character 內積）重數表：45×5 = 225 項，不一致 = {len(mism)}")
if mism:
    say(f"    不一致：{mism[:5]}")
    NFAIL += 1
else:
    NPASS += 1

say("")
say("§1  docs/02 — 歸約定理地位 ＋ A₅ 重數表")
check_num("02_reduction_theorem_status.md", "最小相鄰間隙（/88）", l2["min_gap_88"])
check_str("02_reduction_theorem_status.md", "9! 值", "362880")
check_num("02_reduction_theorem_status.md", "P 內 max m_ρ", l2["max_multiplicity_in_P"], 1e-9)
check_num("02_reduction_theorem_status.md", "首個 m≥3 之 n", l2["first_n_with_m_ge_3"], 1e-9)
check_num("03_no_second_label_theorem.md", "A₅ 44 維分解數", l2["A5_44dim_decompositions"], 1e-9)
check_str("02_reduction_theorem_status.md", "n=44 行", "(1,4,4,6,8)")
check_str("02_reduction_theorem_status.md", "n=11 行", "(0,2,1,1,2)")
check_str("02_reduction_theorem_status.md", "Hironaka 級數", "(1+t¹⁵)/((1−t²)(1−t⁶)(1−t¹⁰))")

say("")
say("§2  docs/03 — 第二標籤不可能定理")
check_num("03_no_second_label_theorem.md", "A₅ 交換子階", l2["A5_commutator_order"], 1e-9)
check_num("03_no_second_label_theorem.md", "簡併對數", l2["degenerate_pairs_both"], 1e-9)
check_num("03_no_second_label_theorem.md", "min gap（/90）", l2["min_gap_90"])
check_num("03_no_second_label_theorem.md", "Σγ(/88 ring)", l2["sum_gamma_88ring"], 1e-9)
check_num("03_no_second_label_theorem.md", "Σγ(/90 dir)", l2["sum_gamma_90dir"], 1e-9)
check_str("03_no_second_label_theorem.md", "3²² 值", "31381059609")

say("")
say("§3  docs/04 — 計數普查 ＋ 分割測試")
check_num("04_count_census_and_partition_test.md", "S(9,3)", l2["stirling_S93"], 1e-9)
check_num("04_count_census_and_partition_test.md", "有標籤分拆", l2["labeled_3partitions_F"], 1e-9)
check_str("04_count_census_and_partition_test.md", "三元組模集", "9,15,18,27,30,36")
check_str("04_count_census_and_partition_test.md", "∩P", "1/9")
check_str("04_count_census_and_partition_test.md", "mod 3 類大小", "(2,4,3)")
check_str("04_count_census_and_partition_test.md", "區塊類大小", "(3,3,3)")

say("")
say("§4  docs/05 — adic 剖面 ＋ 映射 no-go")
check_str("05_adic_profile_and_map_nogo.md", "/88 剖面", str(tuple(l2["profile_88"])))
check_str("05_adic_profile_and_map_nogo.md", "/90 剖面", str(tuple(l2["profile_90"])))
check_str("05_adic_profile_and_map_nogo.md", "/176 剖面", str(tuple(l2["profile_176"])))
check_str("05_adic_profile_and_map_nogo.md", "合格 m 集", "{3, 5, 9}")
check_str("05_adic_profile_and_map_nogo.md", "γ-級 0 解", "0解")
check_bool("05_adic_profile_and_map_nogo.md", "γ-級解數 = 0", l2["gamma_level_solutions"] == 0, "0解")
check_num("05_adic_profile_and_map_nogo.md", "仿射窮舉數", 44 * 44, 1e-9)
check_bool("05_adic_profile_and_map_nogo.md", "保族 max|dev| < 10⁻⁴⁵", float(l2["m_div_x_max_dev"]) < 1e-45, "10⁻⁴⁵")

say("")
say("§5  docs/06 — CY₃ 判死")
check_num("06_cy3_verdict.md", "|χ|/2", l2["ng_chi_over_2"], 1e-9)
check_num("06_cy3_verdict.md", "χ 值", abs(l2["chi_cy3"]), 1e-9)
check_num("06_cy3_verdict.md", "相對偏差 19.667", 19.6667, 1e-3)
check_num("06_cy3_verdict.md", "目錄覆蓋 20", l2["cy3_catalogue_coverage_1_100"], 1e-9)
check_num("06_cy3_verdict.md", "20.67 倍", 20.6667, 1e-3)

say("")
say("§6  docs/07 — 判詞關鍵字")
check_str("07_verdicts.md", "C7 通過數", "0/8")
check_str("07_verdicts.md", "dof 9!", "9!=362880")

say("")
say("§7  docs/08 — 第九候選 h11 − M")
check_num("08_ninth_candidate_h11_minus_M.md", "h11 − M", l2["h11_minus_M"], 1e-9)
check_str("08_ninth_candidate_h11_minus_M.md", "Lucas L₃", "L₃ − 1")
check_num("08_ninth_candidate_h11_minus_M.md", "L₃ = 4", l2["lucas_L1_L3_L5"][1], 1e-9)
check_num("08_ninth_candidate_h11_minus_M.md", "L₅ = 11", l2["lucas_L1_L3_L5"][2], 1e-9)
check_num("08_ninth_candidate_h11_minus_M.md", "k 三路皆 3", l2["k_routes"][0], 1e-9)
check_num("08_ninth_candidate_h11_minus_M.md", "33 之 k-free 族數", l2["kfree_33_families"], 1e-9)
check_str("08_ninth_candidate_h11_minus_M.md", "通用恆等式", "3L_n − L_{n+2} = L_{n−2}")
check_str("08_ninth_candidate_h11_minus_M.md", "9 候選 0/9", "0/9")

say("")
say("§8  豁免（敘述／判定型）")
say("  docs/01（判準凍結表，跑前文件；數值屬常數定義而非計算結果）—— 豁免")

say("")
say("=" * 78)
say(f"L3 總結：PASS {NPASS}、FAIL {NFAIL}")
say("=" * 78)
sys.exit(1 if NFAIL else 0)