#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / 全鏈審計（idcm-audit-verify §3.1 六項 + §7 陷阱）
2026-10-29

只讀 docs/scripts/data（除自身輸出），另以 shell redirect 重跑 L1 腳本至 /tmp 作 byte 比對。

檢查：
  A1 檔案清單（docs/scripts/data/root）
  A2 腳本集對稱：磁碟 scripts/NN_*.py 是否連續、無缺號
  A3 可重現性：每支 NN 腳本跑 RUNS 次 → 與 data/run_NN_*.log byte-identical
  A4 產物綁定：data/*.json ↔ 產者腳本；腳本宣告之輸出存在（無孤兒）
  A5 路徑衛生：scripts/*.py 出現 /tmp/ 寫入模式 = 違規
  A6 四件套狀態行（01_problem_record / SESSION_START / README / final_result）
  A7 文檔衛生：``` 圍欄偶數、U+FFFD 為 0
  A8 README 宣稱檔數 vs 實際

跑法：/usr/bin/python3 scripts/15_full_audit.py [--runs N]
"""
from __future__ import annotations
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOCS = ROOT / "docs"
SCRIPTS = HERE
DATA = ROOT / "data"
SELF = "15_full_audit.py"

LOG, RES = [], {}
NFAIL, NWARN = 0, 0


def say(s=""):
    LOG.append(s)
    print(s)


def chk(tag, ok, detail="", warn=False):
    global NFAIL, NWARN
    if not ok:
        if warn:
            NWARN += 1
            status = "WARN"
        else:
            NFAIL += 1
            status = "FAIL"
    else:
        status = "PASS"
    line = "[%s] %s%s" % (status, tag, ("  |  " + detail) if detail else "")
    say(line)
    RES[tag] = {"status": status, "detail": detail}
    return ok


def hdr(s):
    say("\n" + "=" * 82)
    say(s)
    say("=" * 82)


def main():
    runs = 2
    if "--runs" in sys.argv:
        runs = int(sys.argv[sys.argv.index("--runs") + 1])

    hdr("A1 檔案清單")
    docs = sorted(p.name for p in DOCS.glob("*.md"))
    scripts = sorted(p.name for p in SCRIPTS.glob("*.py"))
    nn = sorted(s for s in scripts if re.match(r"^\d\d_", s) and s != SELF)
    dataj = sorted(p.name for p in DATA.glob("*.json"))
    datal = sorted(p.name for p in DATA.glob("*.log"))
    rootf = sorted(p.name for p in ROOT.iterdir() if p.is_file())
    say("  docs    = %d 份" % len(docs))
    for d in docs:
        say("     %s" % d)
    say("  scripts = %d 支（NN 前綴 %d）" % (len(scripts), len(nn)))
    for s in scripts:
        say("     %s" % s)
    say("  data    = %d json / %d log" % (len(dataj), len(datal)))
    say("  root    = %s" % rootf)
    RES["A1"] = {"docs": docs, "scripts": scripts, "data_json": dataj,
                 "data_log": datal, "root": rootf}

    hdr("A2 腳本集對稱（NN 連續性）")
    # 連續性係「磁碟真相」：須含 SELF（15_full_audit.py 本身亦是 NN 腳本），
    # 否則加入 16/17 後會誤報缺 15。SELF 只在 A3 重跑時排除（避免遞歸）。
    nums = sorted(int(s[:2]) for s in scripts if re.match(r"^\d\d_", s))
    missing = [n for n in range(1, max(nums) + 1) if n not in nums] if nums else []
    chk("A2 NN 編號連續（1..%d）" % (max(nums) if nums else 0), not missing,
        "缺 %s" % missing if missing else "共 %d 支" % len(nums))

    hdr("A3 可重現性（每支跑 %d 次 → byte-identical vs canonical log）" % runs)
    for s in nn:
        stem = s[:-3]
        canon = DATA / ("run_%s.log" % stem)
        if not canon.exists():
            chk("A3 %s 有 canonical log" % s, False, "缺 %s" % canon.name)
            continue
        ok = True
        detail = ""
        first = None
        snap = canon.read_bytes()          # ★ 跑之前快照：腳本會自行覆寫 canonical，避免自證
        for r in range(runs):
            tmp = Path("/tmp") / ("audit_%s_%d.log" % (stem, r))
            with open(tmp, "w") as fh:
                p = subprocess.run(["/usr/bin/python3", str(SCRIPTS / s)],
                                   stdout=fh, stderr=subprocess.STDOUT,
                                   cwd=str(ROOT), timeout=1200)
            if p.returncode != 0:
                ok = False
                detail = "EXIT=%d" % p.returncode
                break
            b = tmp.read_bytes()
            if first is None:
                first = b
            elif b != first:
                ok = False
                detail = "run0 vs run1 不同"
                break
            if b != snap:
                ok = False
                detail = "與跑前 canonical 快照不同（內容漂移）"
                break
        if ok and canon.read_bytes() != snap:
            ok = False
            detail = "canonical 被重寫成不同內容"
        chk("A3 %s 重現 byte-identical" % s, ok, detail or "×%d" % runs)

    hdr("A4 產物綁定")
    orphans = []
    for j in dataj:
        if j == "15_full_audit.json":   # 審計自身產物，非 L1 腳本輸出
            continue
        m = re.match(r"^(\d\d)_", j)
        if not m:
            if j not in ("l2_independent_verify.json", "l3_doc_crosscheck.json",
                             "15_full_audit.json"):
                orphans.append(j)
            continue
        src = [s for s in nn if s.startswith(m.group(1) + "_")]
        if not src:
            orphans.append(j)
    chk("A4a 每個 data/*.json 有產者腳本", not orphans, "孤兒：%s" % orphans if orphans else "")
    # 腳本宣告之輸出：抓 data/ '檔名' 字串
    declared_missing = []
    for s in nn:
        txt = (SCRIPTS / s).read_text(encoding="utf-8")
        for m in re.finditer(r'"(1[0-9]_\w+\.json|[0-9][0-9]_\w+\.json)"', txt):
            name = m.group(1)
            if not (DATA / name).exists():
                declared_missing.append((s, name))
    chk("A4b 腳本宣告之輸出全部存在", not declared_missing,
        "缺：%s" % declared_missing if declared_missing else "")

    hdr("A5 路徑衛生（/tmp 寫入）")
    bad = []
    for s in scripts:
        txt = (SCRIPTS / s).read_text(encoding="utf-8")
        for i, line in enumerate(txt.splitlines(), 1):
            if "/tmp/" not in line:
                continue
            if re.search(r"open\([^)]*/tmp/|write_text\([^)]*/tmp/|Path\(['\"]/tmp/", line):
                bad.append("%s:%d" % (s, i))
    chk("A5 scripts 無 /tmp 寫入", not bad, ", ".join(bad) if bad else "")

    hdr("A6 四件套狀態行")
    for f in ("01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"):
        p = ROOT / f
        if not p.exists():
            chk("A6 %s 存在" % f, False, "檔案不存在", warn=(f == "final_result.md"))
            continue
        txt = p.read_text(encoding="utf-8")
        line = ""
        for ln in txt.splitlines():
            if "狀態" in ln and ("partial" in ln or "closed" in ln or "registered" in ln or "🟡" in ln or "⚪" in ln):
                line = ln.strip()
                break
        say("  %-24s → %s" % (f, line[:110] if line else "（找不到狀態行）"))
        RES.setdefault("A6", {})[f] = line

    hdr("A7 文檔衛生（圍欄偶數 / U+FFFD）")
    for d in docs:
        txt = (DOCS / d).read_text(encoding="utf-8")
        fences = len(re.findall(r"^```", txt, flags=re.M))
        fffd = txt.count("\ufffd")
        chk("A7 %s 圍欄偶數 & 無 U+FFFD" % d,
            fences % 2 == 0 and fffd == 0,
            "fences=%d fffd=%d" % (fences, fffd))

    hdr("A8 README 宣稱 vs 實際")
    rd = (ROOT / "README.md").read_text(encoding="utf-8")
    n_doc_claim = len(re.findall(r"docs/\d\d_", rd))
    chk("A8 README 檔表覆蓋全部 docs", n_doc_claim >= len(docs),
        "README 提及 %d 個 docs/ 條目 vs 實際 %d 份" % (n_doc_claim, len(docs)),
        warn=True)

    hdr("匯總")
    say("  FAIL = %d   WARN = %d" % (NFAIL, NWARN))
    RES["_summary"] = {"fail": NFAIL, "warn": NWARN,
                       "docs": len(docs), "scripts_nn": len(nn), "runs": runs}
    say("")
    say("[written] data/15_full_audit.json , data/run_15_full_audit.log")
    (DATA / "15_full_audit.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_15_full_audit.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())