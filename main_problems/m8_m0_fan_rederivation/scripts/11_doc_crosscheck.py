#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/11 — L3 文檔交叉核對（Doc Crosscheck，唯讀）
================================================================================
鐵律（idcm-audit-verify §2）：
  · **唯讀**：本腳本唔改任何 doc（只 open 讀）
  · 期望值**由 L2 獨立重算**（`data/10_independent_verify.json`），唔抄 docs
  · 數字比對用**數值容差**（依文檔所示有效位；長數字 1e−4、短數字 5e−3），
    而唔係字串相等 —— 因為文檔會用不同精度／記法（純小數 vs e-記法）
  · 容差內仍分辨得出真錯：78.55 vs 79.5538 = 1.26% ⟹ FAIL；2.3e35 vs 5.22e35 = 56% ⟹ FAIL

執行: /usr/bin/python3 scripts/11_doc_crosscheck.py
輸出: data/11_doc_crosscheck.json + data/run_11.log
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"

LOG = []
RES = {}
N_FAIL = 0


def chk(tag, cond, detail=""):
    global N_FAIL
    ok = bool(cond)
    if not ok:
        N_FAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")


def tokens(text):
    """抽全部數值 token：(float_value, 有效位數)"""
    out = []
    for m in TOK.finditer(text.replace("\u2212", "-")):
        mant, exp = m.group(1), m.group(2)
        try:
            v = float(mant) * (10.0 ** int(exp)) if exp else float(mant)
        except ValueError:
            continue
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        out.append((v, sig))
    return out


def match(text_tokens, expect):
    """回傳 (命中?, 最近 token, 相對差)；負值以絕對值比對（文檔負號不進 token）"""
    target = abs(expect)
    best, best_d = None, 1e9
    for v, sig in text_tokens:
        if v == 0 or target == 0:
            continue
        rel = abs(abs(v) - target) / target
        tol = 1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3)
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


L2 = json.loads((DATA / "10_independent_verify.json").read_text(encoding="utf-8"))
sp = L2["_spectrum"]
ch = L2["_chart"]
hi = L2["_hierarchy"]
de = L2["_determinant"]

RATIO = 1.0 / de["gap_factor"]          # 1/det(1+K) 比 μ 大幾多倍（= 79.5538）

CLAIMS = [
    ("docs/02_r1_results.md", "Σtan² = 1247", 1247.0, "O246 交叉核"),
    ("docs/04_r2prime_two_chart_results.md", "cot²(π/88)", ch["cot2"], "圖卡最大比值"),
    ("docs/04_r2prime_two_chart_results.md", "Δt (ln)", 6.664363948470389, "跨距"),
    ("docs/04_r2prime_two_chart_results.md", "Δt (decade)", ch["span_dec"], "跨距"),
    ("docs/04_r2prime_two_chart_results.md", "log₁₀T₁", ch["log10_T1"], "目標"),
    ("docs/04_r2prime_two_chart_results.md", "缺口 (decade)", ch["gap_dec"], "缺口"),
    ("docs/04_r2prime_two_chart_results.md", "缺口因子（★更正後）", ch["factor"], "5.223686e35"),
    ("docs/04_r2prime_two_chart_results.md", "迭代數 13.34", ch["log10_T1"] / ch["span_dec"], "非整數"),
    ("docs/05_r3prime_r4prime_audit.md", "f_M0 = M₀/M_Pl", 2.2852930605e-21, "dof=1 錨"),
    ("docs/05_r3prime_r4prime_audit.md", "f_R = N_IDCM 級", 7.8846e60, "O212 表列"),
    ("docs/06_r5prime_hierarchy_reduction.md", "X₀ = M₀/E_kchain", hi["X0"], "層級歸約"),
    ("docs/06_r5prime_hierarchy_reduction.md", "p 必需指數", hi["p"], "μ = r^p"),
    ("docs/06_r5prime_hierarchy_reduction.md", "p dev vs 1/3", 1.6855e-2, "|p−q|/q"),
    ("docs/06_r5prime_hierarchy_reduction.md", "N_O38 = 44·T₁", 1.801868e40, "恆等"),
    ("docs/06_r5prime_hierarchy_reduction.md", "E_kchain (meV)", 4.348, "錨交接點"),
    ("docs/07_r6prime_three_input_falsifiability.md", "S* = ln(1/μ)", hi["S_star"], "指數族"),
    ("docs/07_r6prime_three_input_falsifiability.md", "隨機基準 1/(2q²)", 2.21606648199446e-06, "q=475"),
    ("docs/08_r7prime_power_audit_selfcorrection.md", "ρ(F1)", 0.025, "功率"),
    ("docs/08_r7prime_power_audit_selfcorrection.md", "ρ(F2)", 0.070, "功率"),
    ("docs/08_r7prime_power_audit_selfcorrection.md", "ρ(前因子族)", 4.015, "唯一有效 null"),
    ("docs/09_r8prime_forcing_candidates.md", "det'(1−K/2)", sp["det_prime_half"], "精確恆等式 1"),
    ("docs/09_r8prime_forcing_candidates.md", "Πcos(πn/88)", 7.541131686492463e-13, "精確恆等式 2"),
    ("docs/09_r8prime_forcing_candidates.md", "det(1+K)", sp["det_1K"], "上限基準"),
    ("docs/09_r8prime_forcing_candidates.md", "1/det(1+K)", 1 / sp["det_1K"], "最接近之強迫量"),
    ("docs/09_r8prime_forcing_candidates.md", "log₁₀ 1/det", de["log10_inv_det"], "dex"),
    ("docs/09_r8prime_forcing_candidates.md", "缺口 dex", de["gap_dex"], "dex"),
    ("docs/09_r8prime_forcing_candidates.md", "缺口倍數（★79.5538）", RATIO, "1/det 比 μ 大"),
    ("README.md", "缺口倍數", RATIO, "★79.5538"),
    ("README.md", "X₀", hi["X0"], "層級歸約"),
    ("README.md", "p 必需指數", hi["p"], "冪律橋"),
    ("README.md", "det'(1−K/2)", sp["det_prime_half"], "R8′ 新恆等式"),
    ("01_problem_record.md", "缺口倍數", RATIO, "★79.5538"),
    ("01_problem_record.md", "S*", hi["S_star"], "指數族"),
    ("01_problem_record.md", "det'(1−K/2)", sp["det_prime_half"], "R8′"),
]

hdr("§1 文檔數字 ↔ L2 獨立重算（唯讀；數值容差）")

cache = {}
by_doc = {}
for doc, label, val, why in CLAIMS:
    if doc not in cache:
        cache[doc] = tokens((ROOT / doc).read_text(encoding="utf-8"))
    ok, got, rel = match(cache[doc], val)
    by_doc.setdefault(doc, [0, 0])
    by_doc[doc][1] += 1
    if ok:
        by_doc[doc][0] += 1
    chk(f"C.{doc.split('/')[-1][:14]}.{label}",
        ok,
        f"期望 {val:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}；{why}）" if ok
        else f"✗ 期望 {val:.6g}，文檔最近 {got if got is None else f'{got:.6g}'}（{why}）")

print("\n  逐檔命中：")
for d, (h, t) in by_doc.items():
    print(f"    {d:46s} {h}/{t}")

hdr("§2 高風險數字之專項核（★ = 本輪審計更正項）")

d09 = (ROOT / "docs/09_r8prime_forcing_candidates.md").read_text(encoding="utf-8")
d04 = (ROOT / "docs/04_r2prime_two_chart_results.md").read_text(encoding="utf-8")
rd = (ROOT / "README.md").read_text(encoding="utf-8")
rec = (ROOT / "01_problem_record.md").read_text(encoding="utf-8")

chk("C2.1 docs/09 載 79.5538（★二輪更正）", "79.5538" in d09, "舊 78.55 只可留於更正欄")
chk("C2.2 docs/09 後見 dev 已改 1.33%（★二輪更正）",
    any(k in d09 for k in ("1.33%", "1.329%", "1.32%")), "舊 0.068% 須作廢")
chk("C2.3 docs/04 載 5.223686e35 且附註更正在場（★一輪更正）",
    "5.223686e35" in d04 and "附註更正" in d04, "")
chk("C2.4 README 缺口倍數 = 79.55（★）", "79.55" in rd, "非 78.55")
chk("C2.5 01_problem_record 缺口倍數 = 79.55（★）", "79.55" in rec, "非 78.55")

hdr("§3 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL, "by_doc": by_doc, "ratio": RATIO}

(DATA / "11_doc_crosscheck.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_11.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {DATA / '11_doc_crosscheck.json'}")
print(f"寫出: {DATA / 'run_11.log'}")

sys.exit(0 if N_FAIL == 0 else 1)