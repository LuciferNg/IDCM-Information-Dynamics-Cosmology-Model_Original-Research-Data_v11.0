#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L1 原樣可重現性 runner（模板；每輪 audit-verify 直接用）
================================================================================
用法
    /usr/bin/python3 scripts/l1_reproducibility.py [--dir <記錄目錄>]
    → data/l1_reproducibility.json ＋ 打印 N/N

做四件事（缺一即假綠）
  1. glob 只用 `[0-9][0-9]_*.py`（L2/L3/audit 用 l2_／l3_／audit_ 前綴，自動排除，唔會自我遞歸）
  2. 每支跑**兩次**，兩次都要求 rc == 0（TIMEOUT／非零＝FAIL，唔可以靠 hash 相同過關）
  3. 比較前剝除**不穩定行**（耗時、時間戳、pid、路徑），否則永遠 FAIL
  4. 兩跑 sha256 相同才算 byte-identical（★ 用檔案內 stable() 邏輯，唔好直接比 raw stdout）

判別訊號
  · 只有一支 FAIL 且其輸出含時間戳／隨機 id ⟹ 非確定性工具，唔係內容漂移（見 SKILL §7 陷阱 82）
  · FAIL 但兩次 rc 都 0 ⟹ 內容漂移：查該腳本是否讀 docs（改完 doc 要重跑）
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

PY = "/usr/bin/python3"
VOLATILE = re.compile(
    r"(耗時|elapsed|耗时|time[:：]\s*[\d.]+|took\s+[\d.]+|[\d.]+\s*s\b|pid[=: ]\d+|0x[0-9a-f]{6,})",
    re.I,
)


def stable(text):
    """剝除不穩定行後嘅 stdout"""
    return "\n".join(ln for ln in text.splitlines() if not VOLATILE.search(ln))


def run_once(path, cwd):
    env = {"PATH": "/usr/bin:/bin", "PYTHONDONTWRITEBYTECODE": "1",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    try:
        p = subprocess.run([PY, str(path)], capture_output=True, text=True,
                           cwd=str(cwd), timeout=1800, env=env)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return None, "", "TIMEOUT"          # ★ None ≠ 0：唔可以當成功


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=".", help="記錄目錄（含 scripts/ 與 data/）")
    a = ap.parse_args()
    root = Path(a.dir).resolve()
    scripts = sorted((root / "scripts").glob("[0-9][0-9]_*.py"))
    data = root / "data"
    data.mkdir(exist_ok=True)
    print(f"L1 腳本數 = {len(scripts)}  （dir={root}）")

    rows, nfail = [], 0
    for i, sc in enumerate(scripts, 1):
        rc1, so1, se1 = run_once(sc, root)
        rc2, so2, se2 = run_once(sc, root)
        s1, s2 = stable(so1), stable(so2)
        h1 = hashlib.sha256(s1.encode()).hexdigest()[:16]
        h2 = hashlib.sha256(s2.encode()).hexdigest()[:16]
        ident = h1 == h2
        ok = (rc1 == 0 and rc2 == 0 and ident)
        if not ok:
            nfail += 1
        rows.append({"script": sc.name, "rc1": rc1, "rc2": rc2,
                     "md5_1": h1, "md5_2": h2, "identical": ident,
                     "stdout_lines": len(s1.splitlines()),
                     "stderr_tail": (se1 or se2 or "")[-200:]})
        print(f"[{'PASS' if ok else 'FAIL'}] {sc.name:44s} rc={rc1}/{rc2} "
              f"identical={ident} ({i}/{len(scripts)})")
        if not ok:
            print(f"        stderr: {(se1 or se2)[:300]}")

    out = {"n_scripts": len(scripts), "n_fail": nfail, "rows": rows}
    (data / "l1_reproducibility.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nL1 RESULT: {len(scripts) - nfail}/{len(scripts)} PASS")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
