#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/11_l1_reproducibility.py — L1 重現性閘（idcm-audit-verify §2.4）
每支核心腳本跑兩次，stdout 必須逐位元組相同且 rc=0。
02_gap_full_scan.g（全量掃描，~27 min）例外：以「重生成檔 md5 == 入庫 md5」證明（A4 等價閘）。
輸出：data/run_l1.json ＋ data/run_l1.log
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"

PY = ["01_arithmetic_reduction.py", "03_variants_and_lookelsewhere.py",
      "05_candidate_principles.py", "06_minimum_cost_principle.py",
      "07_variational_axiom_attack.py", "08_mf1prime_mf2_attack.py",
      "09_mf1prime_dissolution.py"]
G = ["00_gap_named_degrees.g", "02c_order_counts.g", "04_gap_alternative_principles.g"]

LOG, RES, NFAIL, NPASS = [], {}, 0, 0


def chk(tag, ok, detail=""):
    global NFAIL, NPASS
    if ok:
        NPASS += 1
    else:
        NFAIL += 1
    ln = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(ln)
    LOG.append(ln)
    RES[tag] = {"pass": bool(ok), "detail": detail}


def run(cmd):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return subprocess.run(cmd, shell=True, cwd=ROOT, capture_output=True, text=True, env=env)


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def main():
    print("=" * 78 + "\nL1 雙跑重現性\n" + "=" * 78)
    for s in PY:
        a = run(f"/usr/bin/python3 scripts/{s}")
        b = run(f"/usr/bin/python3 scripts/{s}")
        ok = a.returncode == 0 and b.returncode == 0 and a.stdout == b.stdout
        chk(f"L1 py {s}", ok,
            f"rc={a.returncode}/{b.returncode} bytes={len(a.stdout)} identical={a.stdout == b.stdout}")
    for s in G:
        a = run(f"bash scripts/gap_run.sh scripts/{s}")
        b = run(f"bash scripts/gap_run.sh scripts/{s}")
        ok = a.returncode == 0 and b.returncode == 0 and a.stdout == b.stdout and a.stdout.strip() != ""
        chk(f"L1 gap {s}", ok,
            f"rc={a.returncode}/{b.returncode} bytes={len(a.stdout)} identical={a.stdout == b.stdout}")

    # 02：全量掃描以 md5 對入庫值（重生成後比對）
    ref = DATA / "backup/scan_all_500.committed.md5"
    cur = DATA / "run_02_scan_all_500.txt"
    if ref.exists() and cur.exists():
        want = ref.read_text().split()[0]
        got = md5(cur)
        chk("L1 gap 02_gap_full_scan.g（md5 ≡ 入庫）", want == got, f"{got} vs {want}")
    else:
        chk("L1 gap 02_gap_full_scan.g（md5 ≡ 入庫）", False, "缺備份或掃描檔")

    print(f"\nL1 PASS={NPASS} FAIL={NFAIL}")
    (DATA / "run_l1.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    (DATA / "run_l1.json").write_text(
        json.dumps({"pass": NPASS, "fail": NFAIL, "total": NPASS + NFAIL, "detail": RES},
                   ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())