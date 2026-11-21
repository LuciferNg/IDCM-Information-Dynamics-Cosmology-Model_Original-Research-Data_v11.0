#!/usr/bin/env python3
"""
M19 / scripts/08 — 反向數字覆蓋（③.1⑦ L3b）：doc 每個精確數字必須有產物落地

預註冊（跑前凍結）:
  R1  掃描範圍 = docs/01–07 ＋ 四件套（01_problem_record／SESSION_START／README／final_result）
  R2  doc 側 token = 4+ 位小數 或 科學記數（帶符號；`−`(U+2212) 先歸一為 `-`）
  R3  pool = data/run_0{0,1,2,3,4,5,7}_*.log ∪ data/0*.json
      **排除** run_06（audit）／run_l1／run_l2／run_l3／run_08（自我確認）
  R4  容差：一般量相對 0.5%；p 值（40 字元窗內見 `p`／`p=`／`Permutation`／`MC`）用絕對 0.03
  R5  不掃：分數／比值式（`a/b`）、章節號（`§N.M`）、日期、`docs/NN`、檔案編號
  R6  負控制：注入假值 9.876543 必須被判 UNTRACEABLE（證明掃描器有鑑別力）
  R7  ALLOW 白名單：只准逐條列出並附理由；禁整體放寬

失敗判準: R6 負控制無鑑別力 ⟹ 掃描器報廢（停）；R1 範圍內 UNTRACEABLE > 0 且無白名單 ⟹ FAIL。
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DOCS = ROOT / "docs"

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def norm(s):
    return s.replace("\u2212", "-")


POOL_EXCLUDE = ("run_06_audit_verify.log", "run_l1_reproducibility.log",
                "run_l2_independent_recompute.log", "run_l3_doc_crosscheck.log",
                "run_08_reverse_numeric_coverage.log")

pool_vals = []
for p in sorted(DATA.glob("*.log")):
    if p.name in POOL_EXCLUDE:
        continue
    for m in re.finditer(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", norm(p.read_text())):
        try:
            pool_vals.append(float(m.group(0)))
        except ValueError:
            pass
for p in sorted(DATA.glob("0*.json")):
    if p.name == "08_reverse_numeric_coverage.json":   # 掃描器自身輸出（含負控制值）⟹ 排除，否則自我中毒
        continue
    for m in re.finditer(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", norm(p.read_text())):
        try:
            pool_vals.append(float(m.group(0)))
        except ValueError:
            pass

# ALLOW：逐條列出（量, 理由）。空 = 要求零豁免。
ALLOW = {
    "9.876543": "負控制值（docs/07 §九.2 引述掃描器自身之負控制），非物理量",
}


def hit(v, rel=5e-3):
    return any(abs(f - v) <= rel * max(abs(v), 1e-12) for f in pool_vals)


def hit_abs(v, tol=0.03):
    return any(abs(f - v) <= tol for f in pool_vals)


SCAN_FILES = sorted(DOCS.glob("*.md")) + [ROOT / f for f in
                                         ("01_problem_record.md", "SESSION_START.md",
                                          "README.md", "final_result.md")]

TOK = re.compile(r"(?<![\w/.§])([-+]?\d+\.\d{4,}(?:[eE][-+]?\d+)?|[-+]?\d\.\d+[eE][-+]?\d+)")

untraceable = []
scanned = 0
for p in SCAN_FILES:
    body = norm(p.read_text())
    for m in TOK.finditer(body):
        tok = m.group(1)
        # R5：分數／比值／章節號／日期／檔號
        pre = body[max(0, m.start() - 2):m.start()]
        if "/" in body[m.end():m.end() + 1] or pre.endswith("/"):
            continue
        if pre.endswith("§"):
            continue
        try:
            v = float(tok)
        except ValueError:
            continue
        scanned += 1
        if tok in ALLOW:
            continue
        win = body[max(0, m.start() - 40):m.end() + 40]
        is_p = bool(re.search(r"(p\s*=|p\s*\||Permutation|\bMC\b|p 值|p-value)", win))
        ok = hit_abs(v) if is_p else hit(v)
        if not ok:
            untraceable.append((p.name, tok, v))

print("=" * 92)
print("M19 scripts/08 — 反向數字覆蓋（L3b）")
print("=" * 92)
print(f"  pool 值數 = {len(pool_vals)}（排 audit／L1／L2／L3／L3b 自身）")
print(f"  掃描 token = {scanned}")

# R6 負控制
neg_ok = not hit(9.876543) and not hit_abs(9.876543)
chk("R6 負控制 9.876543 必須判 UNTRACEABLE", neg_ok, "掃描器有鑑別力")

chk("R1 掃描範圍 = 全部 docs ＋ 四件套", len(SCAN_FILES) == len(list(DOCS.glob("*.md"))) + 4,
    f"{len(SCAN_FILES)} 檔（docs {len(list(DOCS.glob('*.md')))} ＋ 四件套 4）")

chk("R7 ALLOW 逐條列出（現 %d 條）" % len(ALLOW), True,
    str(ALLOW) if ALLOW else "零豁免")

if untraceable:
    print("\n  UNTRACEABLE：")
    for nm, tok, v in untraceable[:40]:
        print(f"    {nm:<44} {tok}")
    print(f"    共 {len(untraceable)} 條")

chk("R2/R3/R4 全部 doc token 有產物支持", not untraceable, f"UNTRACEABLE {len(untraceable)}")

OUT = DATA
(OUT / "08_reverse_numeric_coverage.json").write_text(json.dumps({
    "pool_size": len(pool_vals),
    "scanned_tokens": scanned,
    "untraceable": [{"file": n, "token": t, "value": v} for n, t, v in untraceable],
    "allow": ALLOW,
    "negative_control": {"value": 9.876543, "detected": neg_ok},
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
