#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M10 audit_verify — 全鏈審計器（idcm-audit-verify §3／§8）
============================================================
A1 檔案齊全（scripts／docs／data／refs）
A2 斷鏈掃描（docs／四件套／README 內引用的 docs|scripts|data 路徑全部存在）
A3 四件套存在 ＋ 水位節存在
A4 L2／L3／L1 三份結果存在且零 FAIL
A5 文檔數字交叉核（scripts/06）C1 100%
A6 register：M10 相關事件之 ref 全部可解析；反向推送 active = 0
A7 README 檔案表涵蓋全部 script／doc（無孤兒）

★ SKIP 清單（審計器不得遞歸）：本檔自身、l1/l2/l3 前綴、append 工具。
"""
import json
import glob
import os
import re
import subprocess
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = "/usr/bin/python3"
PASS, FAIL = [], []


def ck(tag, cond, msg=""):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


print("=" * 78)
print("M10 audit_verify — 全鏈審計")
print("=" * 78)

# ---------------------------------------------------------------- A1
scripts = sorted(os.path.basename(p) for p in glob.glob(os.path.join(BASE, "scripts", "[0-9][0-9]_*.py")))
docs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(BASE, "docs", "*.md")))
jsons = sorted(os.path.basename(p) for p in glob.glob(os.path.join(BASE, "data", "*.json")))
refs = sorted(os.listdir(os.path.join(BASE, "data", "refs")))
print(f"\n[A1] scripts {len(scripts)}｜docs {len(docs)}｜json {len(jsons)}｜refs {len(refs)}")
ck("A1.1 腳本數 = 20", len(scripts) == 20, f"{len(scripts)}")
def has_json(sc):
    pre = sc[:2]
    return any(j.startswith(pre) for j in jsons)
ck("A1.2 每支腳本有對應 data JSON（前綴配對；06 為交叉核不算）",
   all(has_json(s_) for s_ in scripts if not s_.startswith("06_")),
   f"缺：{[s_ for s_ in scripts if not s_.startswith('06_') and not has_json(s_)]}")
ck("A1.3 外部數據 2 份", len(refs) == 2, str(refs))

# ---------------------------------------------------------------- A2 斷鏈
print("\n[A2] 引用路徑斷鏈掃描")
texts = {}
for fn in ["README.md", "SESSION_START.md", "01_problem_record.md", "final_result.md"] + \
          ["docs/" + d for d in docs]:
    p = os.path.join(BASE, fn)
    if os.path.exists(p):
        texts[fn] = open(p, encoding="utf-8").read()
pat_doc = re.compile(r"`?(docs/[0-9]{2}_[A-Za-z0-9_]+\.md)")
pat_scr = re.compile(r"`?(scripts/[0-9]{2}_[A-Za-z0-9_]+\.py)")
pat_dat = re.compile(r"`?(data/[0-9]{2}_[A-Za-z0-9_]+\.json)")
broken = []
for fn, t in texts.items():
    for pat, kind in ((pat_doc, "docs"), (pat_scr, "scripts"), (pat_dat, "data")):
        for m in pat.finditer(t):
            rel = m.group(1)
            if not os.path.exists(os.path.join(BASE, rel)):
                broken.append((fn, rel))
print(f"     引用總數：docs {sum(len(pat_doc.findall(t)) for t in texts.values())}／"
      f"scripts {sum(len(pat_scr.findall(t)) for t in texts.values())}／"
      f"data {sum(len(pat_dat.findall(t)) for t in texts.values())}")
ck("A2.1 零斷鏈", not broken, f"斷鏈 {len(broken)}：" + str(broken[:6]))

# ---------------------------------------------------------------- A3 四件套 ＋ 水位
print("\n[A3] 四件套與水位")
four = ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"]
missing = [f for f in four if not os.path.exists(os.path.join(BASE, f)) or os.path.getsize(os.path.join(BASE, f)) < 1000]
ck("A3.1 四件套齊且非空", not missing, str(missing))
rec = open(os.path.join(BASE, "01_problem_record.md"), encoding="utf-8").read()
ck("A3.2 水位節存在（§N 檔案閉合水位）", "檔案閉合水位" in rec, "有" if "檔案閉合水位" in rec else "缺")

# ---------------------------------------------------------------- A4 L1/L2/L3
print("\n[A4] L1／L2／L3 結果")
for nm, f in (("L1", "l1_reproducibility.json"), ("L2", "l2_independent_verify.json"), ("L3", "l3_doc_crosscheck.json")):
    p = os.path.join(BASE, "data", f)
    if not os.path.exists(p):
        ck(f"A4.{nm} 結果存在", False, f"缺 {f}（先跑 scripts/{nm.lower()}_*.py）")
        continue
    j = json.load(open(p, encoding="utf-8"))
    fail = j.get("fail", j.get("n_fail"))
    ck(f"A4.{nm} 存在且零 FAIL", fail == 0, f"fail = {fail}")

# ---------------------------------------------------------------- A5 06 交叉核
print("\n[A5] 文檔數字交叉核（06）")
r = subprocess.run([PY, os.path.join(BASE, "scripts", "06_doc_crosscheck.py")],
                   cwd=BASE, capture_output=True, text=True)
m = re.search(r"C1 全部 claim 字串在 doc 出現\s+(\d+)/(\d+)", r.stdout)
ck("A5.1 06 C1 全命中", bool(m) and m.group(1) == m.group(2), m.group(0) if m else "未取得")
ck("A5.2 06 全 PASS", "CHECK 合計: PASS" in r.stdout and "FAIL 0" in r.stdout.split("CHECK 合計")[1][:40],
   r.stdout.split("CHECK 合計")[-1].strip()[:60] if "CHECK 合計" in r.stdout else "")

# ---------------------------------------------------------------- A6 register
print("\n[A6] register（M10）")
REG = "/home/wsl/Projects/IDCM-Original/problem_register"
V10 = "/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat"
V11 = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture"
rows = []
for f in glob.glob(os.path.join(REG, "*.jsonl")):
    for ln in open(f, encoding="utf-8"):
        ln = ln.strip()
        if ln:
            try:
                rows.append(json.loads(ln))
            except Exception:
                pass
m10 = [x for x in rows if x.get("pid") == "M10" or x.get("target_pid") == "M10"]


def ok(ref):
    if not ref:
        return True
    cands = [os.path.join(REG, ref), os.path.join(V10, ref), os.path.join(V11, ref)]
    if ref.startswith("/"):
        cands.append(ref)
    return any(os.path.exists(c) for c in cands)


bad = [x for x in m10 if not ok(x.get("ref"))]
ck("A6.1 M10 事件 ref 全部可解析", not bad, f"{len(m10)} 事件，缺 {len(bad)}")
active_back = [x for x in m10 if x.get("pid") == "M10" and x.get("kind") == "handed_off"
               and str(x.get("target_pid") or "").startswith("O")]
withdrawn = {(x.get("target_pid")) for x in rows if x.get("kind") == "handoff_withdrawn"}
ck("A6.2 反向推送（v11→v10）已全撤回", all(str(x.get("target_pid")) in withdrawn for x in active_back),
   f"{len(active_back)} 條，已撤回目標 {sorted(withdrawn & set(str(x.get('target_pid')) for x in active_back))}")

# ---------------------------------------------------------------- A7 README 覆蓋
print("\n[A7] README 檔案表覆蓋")
readme = open(os.path.join(BASE, "README.md"), encoding="utf-8").read()
miss_doc = [d for d in docs if d not in readme]
ck("A7.1 每份 doc 見於 README", not miss_doc, f"未列：{miss_doc[:5]}")
ck("A7.2 腳本範圍已列（01–20）", "`scripts/01`–`20`" in readme or "scripts/01`–`20" in readme,
   "README 之腳本行")

print("\n" + "=" * 78)
print(f"AUDIT 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)
out = {"script": "audit_verify.py", "checks": PASS + FAIL, "pass": len(PASS), "fail": len(FAIL),
       "failures": FAIL, "n_scripts": len(scripts), "n_docs": len(docs), "n_json": len(jsons)}
json.dump(out, open(os.path.join(BASE, "data", "audit_verify.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=2)
sys.exit(0 if not FAIL else 1)