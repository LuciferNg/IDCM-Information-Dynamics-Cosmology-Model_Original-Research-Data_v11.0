#!/usr/bin/env python3
"""
M19 / l1 — 可重現性：每支腳本兩跑，stdout byte-identical，且與入庫 log 相同

判據（跑前凍結）:
  L1.1  每支 scripts/0N 兩跑 md5 相同
  L1.2  兩跑結果與 data/run_0N_<rest>.log 之 md5 相同（入庫即實跑）
  L1.3  rc = 0
"""
import hashlib
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
DATA = ROOT / "data"

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 l1 — 可重現性（兩跑 byte-identical）")
print("=" * 92)

files = sorted(p for p in SCRIPTS.glob("0*.py") if ".pre_patch_" not in p.name)
# 排除自審腳本（其自身 log 會因 A8 之 docs 狀態而變；仍可兩跑比較）
for p in files:
    nn, rest = p.name[:2], p.name[3:-3]
    outs = []
    rcs = []
    for _ in range(2):
        r = subprocess.run(["/usr/bin/python3", str(p)], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"))
        outs.append(r.stdout)
        rcs.append(r.returncode)
    h = [hashlib.md5(o.encode()).hexdigest() for o in outs]
    lg = DATA / f"run_{nn}_{rest}.log"
    log_txt = lg.read_text() if lg.exists() else ""
    h_log = hashlib.md5(log_txt.encode()).hexdigest()
    ok = h[0] == h[1] and all(rc == 0 for rc in rcs)
    chk(f"L1 {p.name} 兩跑一致 ＋ rc=0", ok, f"md5 {h[0][:8]} ×2")
    if not p.name.startswith("06_"):
        chk(f"L1 {p.name} 與入庫 log 一致", h[0] == h_log, f"log md5 {h_log[:8]}")

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
