#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M6 scripts/16_claim_coverage.py — 文檔宣稱覆蓋率台賬（唯讀）
================================================================================
第 2／3 項交付之依據：列出每份 docs／四件套被「數據層（`04_doc_consistency`）」與
「獨立層（`l3_doc_crosscheck`）」覆蓋嘅 claim 數，並標出**零覆蓋**文檔。

判據（事前寫死）
  P1 覆蓋來源只認兩個已註冊驗證器：`scripts/04_doc_consistency.py`（data↔doc）與
     `scripts/l3_doc_crosscheck.py`（L2↔doc）。唔自行新增第三套判準。
  P2 「零覆蓋」文檔必須逐份交代**為何可以零覆蓋**（純定義／無數值宣稱／由他記錄持有）。
  P3 輸出 `data/16_claim_coverage.json` ＋ `data/run_16.log`（log 印 `n/n PASS` 供審計 A6）。
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
DATA = ROOT / "data"

DOCS = sorted((ROOT / "docs").glob("*.md")) + [ROOT / n for n in
                                               ("01_problem_record.md", "README.md",
                                                "final_result.md", "SESSION_START.md")]

# 零覆蓋之豁免（必須有理由）
EXEMPT_ZERO = {
    "SESSION_START.md": "第 1 輪不可變快照（原則不改）；其數字全部係 docs/01 之前身，已由 04 覆蓋對應 docs",
}

NAME_RE = re.compile(r'["\']?(?:docs/[0-9a-z_]+\.md|(?:01_problem_record|final_result|README|SESSION_START)\.md)["\']?')


def claims_of(script):
    p = ROOT / "scripts" / script
    if not p.exists():
        return {}
    txt = p.read_text(encoding="utf-8")
    # 只數「claim 區塊」內嘅檔名（CLAIMS = [ ... ] / by_doc 統計用）
    m = re.search(r"CLAIMS\s*=\s*\[(.*?)\n\]", txt, re.S)
    body = m.group(1) if m else txt
    out = {}
    for d in NAME_RE.findall(body):
        d = d.strip("\"'")
        out[d] = out.get(d, 0) + 1
    return out


def main():
    print("=" * 96)
    print("M6 scripts/16 — 文檔宣稱覆蓋率台賬（唯讀）")
    print("=" * 96)
    cov04 = claims_of("04_doc_consistency.py")
    covl3 = claims_of("l3_doc_crosscheck.py")

    rows, zero = [], []
    for doc in DOCS:
        if not doc.exists():
            continue
        rel = ("docs/" + doc.name) if doc.parent.name == "docs" else doc.name
        n04, n3 = cov04.get(rel, 0), covl3.get(rel, 0)
        rows.append({"doc": rel, "claims_04": n04, "claims_l3": n3, "total": n04 + n3})
        mark = "✅" if (n04 + n3) else "🔴"
        print(f"  {mark} {rel:44s} 04={n04:2d}  L3={n3:2d}  合計={n04+n3:2d}")
        if not (n04 + n3):
            zero.append(rel)

    print(f"\n零覆蓋文檔 = {len(zero)}：{zero}")
    for z in zero:
        why = EXEMPT_ZERO.get(z)
        if why is None:
            print(f"  🔴 {z} 未豁免 —— 須加 claim 或寫明理由")

    unexempt = [z for z in zero if z not in EXEMPT_ZERO]
    npairs = len(rows)
    npass = npairs - len(unexempt)
    print(f"\n{npass}/{npairs} PASS（覆蓋台賬：零覆蓋且未豁免 = {len(unexempt)}）")

    out = {"rows": rows, "zero_coverage": zero, "exempt": EXEMPT_ZERO,
           "unexempt": unexempt, "total_claims_04": sum(cov04.values()),
           "total_claims_l3": sum(covl3.values())}
    (DATA / "16_claim_coverage.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("→ data/16_claim_coverage.json")
    return 0 if not unexempt else 1


if __name__ == "__main__":
    raise SystemExit(main())