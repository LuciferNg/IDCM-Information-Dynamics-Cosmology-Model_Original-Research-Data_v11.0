#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/10_audit_verify.py — M16 全鏈審計（idcm-audit-verify §3.1 六項 ＋ §4.8 機器閘）
================================================================
檢查
  ① 腳本集對稱      磁碟 scripts/ ↔ 文檔／README 宣告集
  ② 逆向綁定        headline 數字 → 有 producer script 且該 script 有落地 log
  ③ 路徑衛生        scripts/ 內零 /tmp、.hermes/cache、/var/tmp、~scratch；GAP 走 data/gap/
  ④ 產物綁定        每個 data/run_*.log 有 producer；每個腳本宣告輸出存在
  ⑤ 幽靈推導／⑦ 反向覆蓋  docs ≥6 位有效數字 token 必須可在 data 池找到（rel ≤ 1e-3）
  ⑥ driver          scripts/run_chain.sh 存在且覆蓋全部 .py/.g（除基礎設施）
  A8a 舊值掃描      已更正之數值每出現點須帶標註詞
  A10 時序閘        data 產物 mtime ≥ producer script mtime
  A12 commit 完整性  git status 對本目錄乾淨（排除審計自身輸出）＋ 逐檔在 git 內
  G1  GAP 部署      data/gap/ 必須被 git ignore（二進位不入庫）；setup 腳本存在
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DOCS = ROOT / "docs"

LOG, RES, NFAIL = [], {}, 0
NPASS_BY_SEC = {}


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append(s)


def sh(cmd):
    return subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True).stdout.strip()


# ── 宣告集（腳本）────────────────────────────────────────────
CORE_PY = ["01_arithmetic_reduction.py", "03_variants_and_lookelsewhere.py",
           "05_candidate_principles.py", "06_minimum_cost_principle.py",
           "07_variational_axiom_attack.py", "08_mf1prime_mf2_attack.py",
           "09_mf1prime_dissolution.py"]
CORE_G = ["00_gap_named_degrees.g", "02_gap_full_scan.g", "02c_order_counts.g",
          "04_gap_alternative_principles.g"]
INFRA = ["02b_chunk_driver.py", "run_02_parallel.sh", "gap_run.sh", "gap_setup.sh",
         "run_chain.sh", "l2_independent_verify.py", "l3_doc_crosscheck.py",
         "10_audit_verify.py", "11_l1_reproducibility.py"]


def sec(name):
    hdr(name)
    NPASS_BY_SEC.setdefault(name, [0, 0])


def main():
    # ───────── ① 腳本集對稱 ─────────
    sec("① 腳本集對稱")
    disk = sorted(p.name for p in HERE.iterdir() if p.suffix in (".py", ".g", ".sh"))
    declared = sorted(CORE_PY + CORE_G + INFRA)
    missing_disk = [x for x in declared if x not in disk]
    unlisted = [x for x in disk if x not in declared]
    chk("①a 宣告集全部在磁碟", not missing_disk, f"缺 {missing_disk}")
    chk("①b 磁碟無未宣告腳本", not unlisted, f"未列 {unlisted}")

    # 文檔／README 內引用之 scripts/NN 是否全部存在
    refs = set()
    for f in list(DOCS.glob("*.md")) + [ROOT / "README.md", ROOT / "SESSION_START.md",
                                        ROOT / "01_problem_record.md", ROOT / "final_result.md"]:
        if not f.exists():
            continue
        txt = f.read_text(encoding="utf-8")
        for mm in re.finditer(r"scripts/(\d{2}[a-z]?_[a-z0-9_]+\.(?:py|g))", txt):
            pre = txt[max(0, mm.start() - 90): mm.start()] + mm.group(0)
            # 跨記錄引用（o171/scripts/…、m15_…/scripts/…、`o238_…/`）不屬本 O，不計
            if re.search(r"[a-z0-9_]+/", pre) or re.search(r"\bo\d{3}_|\bm\d{2}_", pre):
                continue
            refs.add(mm.group(1))
    chk("①c 文檔引用之腳本全部存在", all((HERE / r).exists() for r in refs),
        f"缺 {[r for r in refs if not (HERE / r).exists()]}")

    # ───────── ② 逆向綁定 ─────────
    sec("② 逆向綁定（headline 數字 → producer ＋ log）")
    PAIRS = [
        ("92570", "03_variants_and_lookelsewhere.py"),
        ("92570", "03_variants_and_lookelsewhere.py"),
        ("17", "03_variants_and_lookelsewhere.py"),
        ("267", "03_variants_and_lookelsewhere.py"),
        ("586", "03_variants_and_lookelsewhere.py"),
        ("37.77", "03_variants_and_lookelsewhere.py"),
        ("11", "05_candidate_principles.py"),
        ("60", "06_minimum_cost_principle.py"),
        ("0.281304", "07_variational_axiom_attack.py"),
        ("1.0142", "08_mf1prime_mf2_attack.py"),
        ("55.93", "08_mf1prime_mf2_attack.py"),
        ("1.79", "09_mf1prime_dissolution.py"),
    ]
    LOGG = {"01_arithmetic_reduction.py": "run_01.log",
            "03_variants_and_lookelsewhere.py": "run_03.log",
            "05_candidate_principles.py": "run_05.log",
            "06_minimum_cost_principle.py": "run_06.log",
            "07_variational_axiom_attack.py": "run_07.log",
            "08_mf1prime_mf2_attack.py": "run_08.log",
            "09_mf1prime_dissolution.py": "run_09.log"}
    for num, scr in PAIRS:
        lg = DATA / LOGG[scr]
        ok = lg.exists() and num in lg.read_text(encoding="utf-8", errors="ignore")
        chk(f"② {num} ← {scr}", ok, f"log={lg.name}")

    # ───────── ③ 路徑衛生 ─────────
    sec("③ 路徑衛生")
    BAD = re.compile(r"/tmp/|\.hermes/cache|/var/tmp|~/scratch|/home/wsl/\.cache")
    hits = []
    for f in HERE.iterdir():
        if f.suffix not in (".py", ".g", ".sh"):
            continue
        if f.name == "10_audit_verify.py":
            continue
        for i, ln in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if BAD.search(ln) and not ln.lstrip().startswith("#"):
                hits.append(f"{f.name}:{i}:{ln.strip()[:60]}")
    chk("③a scripts 內零 prune 路徑（非註釋）", not hits, f"{hits}")
    chk("③b gap_run.sh 指向記錄內 data/gap/root",
        "data/gap/root" in (HERE / "gap_run.sh").read_text(encoding="utf-8"))
    chk("③c gap_setup.sh 存在（可重建部署）", (HERE / "gap_setup.sh").exists())

    # ───────── ④ 產物綁定 ─────────
    sec("④ 產物綁定")
    PRODUCER = {"run_00_named_degrees.log": "00_gap_named_degrees.g",
                "run_01.log": "01_arithmetic_reduction.py",
                "run_02_scan_all_500.txt": "02_gap_full_scan.g",
                "order_counts.txt": "02c_order_counts.g",
                "run_03.log": "03_variants_and_lookelsewhere.py",
                "run_04_gap_alt.log": "04_gap_alternative_principles.g",
                "run_05.log": "05_candidate_principles.py",
                "run_06.log": "06_minimum_cost_principle.py",
                "run_07.log": "07_variational_axiom_attack.py",
                "run_08.log": "08_mf1prime_mf2_attack.py",
                "run_09.log": "09_mf1prime_dissolution.py",
                "l2_independent_verify.json": "l2_independent_verify.py",
                "l3_doc_crosscheck.json": "l3_doc_crosscheck.py"}
    orphan = [k for k in PRODUCER if not (DATA / k).exists()]
    chk("④a 宣告產物全部存在", not orphan, f"缺 {orphan}")
    HARNESS = {"run_l2.log", "run_l3.log", "run_audit.log", "run_l1.log",
               "run_l1_boot.log", "run_l2_boot.log", "run_l3_boot.log", "run_audit_boot.log",
               "run_gap_setup.log", "run_l2_raw.log", "run_02.log", "run_02_driver.log",
               "run_02_parallel.log"}
    logs = [p.name for p in DATA.glob("run_*.log") if p.name not in HARNESS]
    unknown = [x for x in logs if x not in PRODUCER]
    chk("④b 無孤兒 data/run_*.log", not unknown, f"孤兒 {unknown}")

    # ───────── ⑤／⑦ 反向數值覆蓋 ─────────
    sec("⑤⑦ 反向數值覆蓋（docs ≥6 位有效數字 → data 池）")
    pool_txt = []
    for p in DATA.glob("**/*"):
        if p.is_file() and p.suffix in (".log", ".txt", ".json") and "backup" not in str(p):
            try:
                pool_txt.append(p.read_text(encoding="utf-8", errors="ignore"))
            except Exception:
                pass
    pool = "\n".join(pool_txt)
    pool_vals = []
    for m in re.finditer(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?", pool.replace("\u2212", "-")):
        try:
            pool_vals.append(float(m.group(1) + (("e" + m.group(2)) if m.group(2) else "")))
        except ValueError:
            pass
    TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")
    orphan_tokens = {}
    for f in sorted(DOCS.glob("*.md")):
        txt = f.read_text(encoding="utf-8").replace("\u2212", "-")
        txt = re.sub(r"(?<=\d),(?=\d{3})", "", txt)
        bad = []
        for m in TOK.finditer(txt):
            mant, exp = m.group(1), m.group(2)
            sig = len(mant.replace(".", "").lstrip("0"))
            if sig < 6:
                continue
            if re.match(r"\s*ms", txt[m.end():m.end() + 4]) or " ms" in txt[max(0, m.start() - 40):m.start()]:
                continue  # 執行時間屬環境量，非推導值
            try:
                v = float(mant + (("e" + exp) if exp else ""))
            except ValueError:
                continue
            if v == 0:
                continue
            if not any(abs(abs(pv) - abs(v)) / abs(v) <= 1e-3 for pv in pool_vals if pv != 0):
                bad.append(m.group(0))
        if bad:
            orphan_tokens[f.name] = sorted(set(bad))[:8]
    chk("⑤⑦ docs 高精度數值全部可追溯", not orphan_tokens, f"{orphan_tokens}")

    # ───────── ⑥ driver ─────────
    sec("⑥ driver 覆蓋")
    drv = HERE / "run_chain.sh"
    chk("⑥a run_chain.sh 存在", drv.exists())
    if drv.exists():
        txt = drv.read_text(encoding="utf-8")
        cover = set(re.findall(r"(?:(?:\d{2}[a-z]?|l[23])_[a-z0-9_]+\.(?:py|g))", txt))
        # 傳遞覆蓋：L1 哈尼斯所跑之腳本集
        l1 = HERE / "11_l1_reproducibility.py"
        if l1.exists() and "11_l1_reproducibility.py" in txt:
            cover |= set(re.findall(r"\"((?:\d{2}[a-z]?|l[23])_[a-z0-9_]+\.(?:py|g))\"",
                                    l1.read_text(encoding="utf-8")))
        # 全量掃描 02 之入口：run_02_parallel.sh（由 02b 生成分片）
        for helper in ["run_02_parallel.sh", "02b_chunk_driver.py"]:
            hp = HERE / helper
            if hp.exists() and helper in txt:
                cover |= set(re.findall(r"((?:\d{2}[a-z]?)_[a-z0-9_]+\.(?:py|g))",
                                        hp.read_text(encoding="utf-8")))
        need = set(CORE_PY + CORE_G + ["l2_independent_verify.py", "l3_doc_crosscheck.py",
                                       "10_audit_verify.py"])
        chk("⑥b driver 覆蓋全部核心＋L2/L3/audit", need <= cover, f"缺 {sorted(need - cover)}")

    # ───────── A8a 舊值掃描 ─────────
    sec("A8a 舊值掃描")
    STALE = {
        "1/6 = 16.7%": ("MF1′ 代價", ["docs/06", "docs/07", "docs/08"]),
        "γ-等差單獨唯一": None,
    }
    ann = re.compile(r"更正|舊|已廢|正解|審計|撤回|不適用|失效|降級|superseded")
    bad_lines = []
    for f in sorted(DOCS.glob("*.md")):
        lines = f.read_text(encoding="utf-8").splitlines()
        for i, ln in enumerate(lines, 1):
            if "16.7%" not in ln:
                continue
            ctx = "\n".join(lines[max(0, i - 12): i + 1])
            if ann.search(ln) or ann.search(ctx):
                continue
            bad_lines.append(f"{f.name}:{i}")
    chk("A8a-1 「16.7%」每出現點帶標註", not bad_lines, f"{bad_lines}")

    # ───────── A10 時序閘 ─────────
    sec("A10 時序閘（log mtime ≥ producer mtime）")
    stale = []
    for prod, scr in PRODUCER.items():
        src = HERE / scr
        dst = DATA / prod
        if src.exists() and dst.exists() and dst.stat().st_mtime + 1 < src.stat().st_mtime:
            stale.append(f"{prod} < {scr}")
    chk("A10 全部產物新於 producer", not stale, f"{stale}")

    # ───────── A12 commit 完整性 ─────────
    sec("A12 commit 完整性")
    AUDIT_OUT = ("data/run_audit.log data/audit_verify.json data/run_l1.json data/run_l1.log "
                 "data/run_l1_boot.log data/run_l2.log data/run_l3.log data/l2_independent_verify.json "
                 "data/l3_doc_crosscheck.json data/run_l2_boot.log data/run_l3_boot.log "
                 "data/run_audit_boot.log data/run_gap_setup.log data/run_l2_raw.log")
    excl = " ".join(f'| grep -v "{p}"' for p in AUDIT_OUT.split())
    dirty = sh(f'git status --porcelain . | grep -v "data/backup/" | grep -v "data/gap/" '
               f'| grep -v "^?? data/chunk_" {excl} | wc -l')
    chk("A12a 工作樹乾淨（排除審計輸出）", dirty.strip() in ("0", ""), f"dirty={dirty}")
    tracked_missing = []
    for p in ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"] + \
             [f"docs/{p.name}" for p in sorted(DOCS.glob("*.md"))] + \
             [f"scripts/{p.name}" for p in sorted(HERE.iterdir()) if p.suffix in (".py", ".g", ".sh")]:
        if not (ROOT / p).exists():
            tracked_missing.append(f"NODISC {p}")
        elif not sh(f'git ls-files --error-unmatch "{p}" 2>/dev/null'):
            tracked_missing.append(f"UNTRACKED {p}")
    chk("A12b 全部文檔／腳本已入 git", not tracked_missing, f"{tracked_missing[:6]}")

    # ───────── G1 GAP 部署 ─────────
    sec("G1 GAP 部署")
    ig = sh('git check-ignore -v data/gap/root/usr/share/gap 2>/dev/null || true')
    chk("G1a data/gap/ 被 git 忽略", "data/gap" in ig, ig[:80])
    chk("G1b GAP 可執行存在", (DATA / "gap/root/usr/lib/x86_64-linux-gnu/gap/gap").exists())

    # ───────── L1/L2/L3 閘 ─────────
    sec("L1／L2／L3 閘")
    l1 = DATA / "run_l1.json"
    if l1.exists():
        j = json.loads(l1.read_text(encoding="utf-8"))
        chk("L1 雙跑 byte-identical 全過", j.get("fail") == 0,
            f"{j.get('pass')}/{j.get('total')} fail={j.get('fail')}")
    else:
        chk("L1 結果存在", False, "缺 data/run_l1.json")
    l2 = json.loads((DATA / "l2_independent_verify.json").read_text(encoding="utf-8"))
    chk("L2 FAIL=0", l2["_summary"]["fail"] == 0, f"{l2['_summary']}")
    l3 = json.loads((DATA / "l3_doc_crosscheck.json").read_text(encoding="utf-8"))
    chk("L3 FAIL=0", l3["_summary"]["fail"] == 0, f"{l3['_summary']['fail']}")

    hdr("匯總")
    print(f"audit 檢查項 = {len(RES)}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": len(RES), "fail": NFAIL}
    (DATA / "run_audit.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    (DATA / "audit_verify.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2),
                                            encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())