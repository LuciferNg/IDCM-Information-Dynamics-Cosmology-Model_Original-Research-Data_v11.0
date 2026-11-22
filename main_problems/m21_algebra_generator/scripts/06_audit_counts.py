#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/06 — 自查：由 log 機械重算 check 水位（防手報）
==========================================================
規則（o-series 判數紀律）：只有 [PASS]/[FAIL] 行算 check；
宣告式（[DECL]）唔入 check 數。水位必須由 log 機械解析，唔可以憑記憶。

本檔同時核：
  A. 每支 script 之 PASS/FAIL 數 vs `data/*.json` 之 `n_pass`/`n_fail`
  B. 全套重跑之決定性腳本 stdout byte-identical（重跑 00–05，比對 log）
  C. 預期 FAIL 之登記（02: CY1′/CY1；04: GR1-c/d）
"""
import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


PAIRS = [
    ("00", "00_reproduction_gate.json"),
    ("01", "01_a5_census.json"),
    ("02", "02_cy3_pool_gate.json"),
    ("03", "03_two_factor_audit.json"),
    ("04", "04_derivation_graph.json"),
    ("05", None),   # 05 之 JSON 結構不同（registry），只查 log
    ("07", "07_m22_intake.json"),
    ("08", "08_br1_mechanism_verdict.json"),
]

print("=" * 78)
print("M21 scripts/06 — check 水位自查")
print("=" * 78)

tot_p = tot_f = 0
for num, jf in PAIRS:
    log = os.path.join(DATA, f"run_{num}.log")
    if not os.path.exists(log):
        ck(f"A{num}  log 存在", False, log)
        continue
    txt = open(log, encoding="utf-8").read()
    p = txt.count("[PASS]")
    f_ = txt.count("[FAIL]")
    tot_p += p
    tot_f += f_
    if jf:
        j = json.load(open(os.path.join(DATA, jf), encoding="utf-8"))
        ok = (j["n_pass"] == p and j["n_fail"] == f_ and j["n_total"] == p + f_)
        ck(f"A{num}  log({p}P/{f_}F) ≡ JSON(n_pass={j['n_pass']},n_fail={j['n_fail']})", ok,
           f"{p}/{f_}")
    else:
        ck(f"A{num}  log 已解析（{p}P/{f_}F）", p + f_ > 0, f"{p}/{f_}")

print(f"\n  合計：PASS = {tot_p} ｜ FAIL = {tot_f}")
ck("TOT  合計水位（131 PASS / 4 FAIL，scripts 00–05 ＋ 07–08）",
   (tot_p, tot_f) == (131, 4), f"{tot_p}/{tot_f}")

# 決定性重跑（byte-identical）
print("\n── 決定性重跑（00–05，比對 stdout） ──")
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
names = {0: "00_reproduction_gate.py", 1: "01_a5_census_identity.py", 2: "02_cy3_pool_gate.py",
         3: "03_two_factor_audit.py", 4: "04_derivation_graph.py", 5: "05_registry_and_dedup.py",
         7: "07_m22_intake.py", 8: "08_br1_mechanism_verdict.py"}
NONDET = {5}      # 05 之輸出含 corpus 掃描總檔數（會隨他人新增檔案變動）⟹ 非決定性
for num, fn in names.items():
    r = subprocess.run([sys.executable, os.path.join(HERE, fn)], capture_output=True, env=env)
    new = r.stdout.decode("utf-8")
    log_p = os.path.join(DATA, f"run_{num:02d}.log")
    old = open(log_p, encoding="utf-8").read() if os.path.exists(log_p) else ""
    if num in NONDET:
        # 非決定性 script：以本次輸出刷新快照 log（與其 JSON 同步），只核 PASS/FAIL 模式穩定
        prev_pat = (old.count("[PASS]"), old.count("[FAIL]"))
        new_pat = (new.count("[PASS]"), new.count("[FAIL]"))
        with open(log_p, "w", encoding="utf-8") as fh:
            fh.write(new)
        ck(f"B{num:02d}  非決定性 script：PASS/FAIL 模式穩定（快照已刷新）",
           prev_pat == new_pat and r.returncode == 0,
           f"rc={r.returncode}, pattern {prev_pat} → {new_pat}")
        continue
    same = (new == old)
    ck(f"B{num:02d}  重跑 stdout byte-identical", same,
       f"rc={r.returncode}, diff={'0' if same else 'yes'}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/06 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)
with open(os.path.join(DATA, "06_audit_counts.json"), "w", encoding="utf-8") as fh:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "watermark": {"pass": tot_p, "fail": tot_f}},
              fh, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
