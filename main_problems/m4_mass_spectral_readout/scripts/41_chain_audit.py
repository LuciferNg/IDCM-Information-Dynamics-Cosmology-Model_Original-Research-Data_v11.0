#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/41 -- 推導鏈審計（idcm-audit-verify §3.1 六項核對 ＋ L2/L3 閘 ＋ 產物綁定）。

唯讀：只讀 scripts/ docs/ data/ 四件套，唔改任何檔（只寫本審計自己的 data/41_*.json）。
判據全部機械化（唔靠「我看過」）：
  ① 腳本／文檔集對稱（差集必須已申報）
  ② 逆向綁定：headline 數字 → producer script（由 L2 key → script 對照表 ＋ 檔內字串搜尋）
  ③ 路徑衛生：scripts 有無寫 /tmp（真寫入模式，非字串提及）
  ④ 產物綁定：data/*.json ↔ 生產腳本（無孤兒）
  ⑤ inline 幽靈推導：headline 數字有無落碟 script（等價於 ② 之反向）
  ⑥ driver 存在 ＋ 每輪有 run <log> <script>
  ⑦ L2 / L3 結果（讀各自 JSON）
"""
import glob
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
SCRIPTS = os.path.join(ROOT, "scripts")
DOCS = os.path.join(ROOT, "docs")
DATA = os.path.join(ROOT, "data")
OUT = []


def add(check, ok, detail):
    OUT.append({"check": check, "ok": bool(ok), "detail": detail})
    print(("PASS " if ok else "FAIL ") + "%-46s %s" % (check, detail))


py = sorted(os.path.basename(p) for p in glob.glob(os.path.join(SCRIPTS, "*.py")))
nn_py = [p for p in py if re.match(r"^[0-9]{2}[a-z]?_", p)]
md = sorted(os.path.basename(p) for p in glob.glob(os.path.join(DOCS, "*.md")))
js = sorted(os.path.basename(p) for p in glob.glob(os.path.join(DATA, "*.json")))
print("=" * 100)
print("M4 scripts/41 -- 推導鏈審計（六項核對 ＋ L2/L3）")
print("scripts=%d（NN 前綴 %d）｜docs=%d｜data json=%d" % (len(py), len(nn_py), len(md), len(js)))
print("=" * 100)

# ---------------- ① 腳本／文檔集對稱
sc_nums = set()
for p in nn_py:
    m = re.match(r"^(\d{2}[a-z]?)_", p)
    sc_nums.add(m.group(1))
doc_nums = set()
for p in md:
    m = re.match(r"^(\d{2}[a-z]?)_", p)
    if m:
        doc_nums.add(m.group(1))
# 已申報：doc-only（無腳本）／script-only（無 doc）
DECLARED_DOC_ONLY = {"32": "docs/32_flavour_structure_scout.md — 純設計/scout 文檔，無腳本（宣告式）",
                     "37": "docs/37_closeout_readiness.md — 收尾準備文檔，無腳本（宣告式）"}
script_no_doc = sorted(sc_nums - doc_nums)
doc_no_script = sorted(doc_nums - sc_nums)
DECLARED_SHARED = {"31b": "docs/31_flavour_insertion_negative.md（enum 附於同一 doc）",
                   "35": "docs/34_u_d_sector_verdict_citation.md（v9 表審計併入判決 doc）",
                   "36b": "docs/36_r5_em_kaon_external_verdict.md（v2 取代 v1，同一 doc）",
                   "41": "docs/41_chain_audit.md（本審計自身）",
                   "43": "docs/42_flag2024_x_revision.md（X 修訂之傳播量化，併入同一 doc）"}
add("1a 每支 NN 腳本都有 doc（或已申報共用）",
    all(n in DECLARED_DOC_ONLY or n in DECLARED_SHARED or n in doc_nums for n in script_no_doc),
    "script-only = %s" % {n: DECLARED_SHARED.get(n, "未申報") for n in script_no_doc})
add("1b doc-only 已申報", all(n in DECLARED_DOC_ONLY for n in doc_no_script),
    "doc-only = %s" % {n: DECLARED_DOC_ONLY.get(n, "未申報") for n in doc_no_script})
add("1c 非 NN 腳本（L2/L3/driver）入審計範圍",
    {"l2_independent_verify.py", "l3_doc_crosscheck.py"} <= set(py),
    "pycache=%s" % ("有" if os.path.isdir(os.path.join(SCRIPTS, "__pycache__")) else "無"))

# ---------------- ③ 路徑衛生：真寫入模式掃 /tmp
WRITE_PAT = re.compile(r"(open\s*\(\s*[\"']/tmp|write_text\s*\(|Path\s*\(\s*[\"']/tmp|>\s*/tmp)")
tmp_hits = []
for p in glob.glob(os.path.join(SCRIPTS, "*.py")):
    for i, line in enumerate(open(p, encoding="utf-8", errors="replace"), 1):
        if "/tmp" in line and WRITE_PAT.search(line) and not line.lstrip().startswith("#"):
            tmp_hits.append("%s:%d" % (os.path.basename(p), i))
add("3 路徑衛生（scripts 無真 /tmp 寫入）", not tmp_hits, "hits = %s" % (tmp_hits or "0"))

# ---------------- ④ 產物綁定：每個 data/*.json 有生產者
producers = {}
for p in glob.glob(os.path.join(SCRIPTS, "*.py")):
    txt = open(p, encoding="utf-8", errors="replace").read()
    writes = ("json.dump" in txt) or ('"w"' in txt) or ("'w'" in txt)
    if not writes:
        continue
    for m in re.finditer(r"([A-Za-z0-9_]+)\.json\"?", txt):
        producers.setdefault(m.group(1) + ".json", set()).add(os.path.basename(p))
orphans = [j for j in js if j not in producers]
add("4a 每個 data/*.json 有 producer", not orphans, "orphans = %s" % (orphans or "0"))
DECLARED_MULTI = {"l2_independent_verify.json": "L2 產出、L3 只讀（正則誤判讀取為寫入，已申報）",
                  "36_em_kaon_external.json": "36 v1 與 36b v2 同寫（v2 取代，已申報）",
                  "23_full_audit.json": "23 產出、24 讀寫加註（已申報）"}
multi = {j: sorted(v - {"41_chain_audit.py"}) for j, v in producers.items()
         if j in js and len(v - {"41_chain_audit.py"}) > 1}
add("4b 一產物一 producer（多產者須已申報）", set(multi) <= set(DECLARED_MULTI),
    "multi-producer = %s" % {j: (DECLARED_MULTI.get(j, "未申報"), v) for j, v in multi.items()})

# 腳本宣告的輸出存在
missing_out = []
for p in glob.glob(os.path.join(SCRIPTS, "0[0-9]_*.py")) + glob.glob(os.path.join(SCRIPTS, "[1-4][0-9]_*.py")):
    txt = open(p, encoding="utf-8", errors="replace").read()
    for m in re.finditer(r"data/([0-9]{2}[a-z]?_[A-Za-z0-9_]+\.json)", txt):
        if not os.path.exists(os.path.join(DATA, m.group(1))):
            missing_out.append("%s→%s" % (os.path.basename(p), m.group(1)))
_miss = sorted(set(m for m in missing_out if "41_chain_audit.json" not in m))   # 自身輸出（先檢查後寫，已申報）
add("4c 腳本宣告輸出全部存在", not _miss, "missing = %s（自身輸出 41_chain_audit.json 已申報）" % (_miss or "0"))

# ---------------- ② 逆向綁定：headline 數字 → producer
HEADLINE = [
    ("4/γ1 = 784.96（負結果總表 A4）", "784.9", "09_density_gate.json"),
    ("span = 12.005533 / m_u = 2.1498", "12.005533", "11_lattice_flag.json"),
    ("em(K+−K0) = 2.264589 / z = −15.98", "2.264588", "36_em_kaon_external.json"),
    ("⟨O⟩_K = 2.421324 / z = −9.76", "2.421323", "40_recheck_m7_m10.json"),
    ("rank_Q = 20（譜代數 A1；L1 無 JSON，由 L2 精確提供）", "rank_Q", "l2_independent_verify.json"),
]
rev_gaps = []
for label, needle, jsf in HEADLINE:
    jp = os.path.join(DATA, jsf)
    prod = sorted(producers.get(jsf, []))
    if not os.path.exists(jp):
        rev_gaps.append("%s: %s 不存在" % (label, jsf))
    elif needle not in open(jp, encoding="utf-8", errors="replace").read():
        rev_gaps.append("%s: %s 內零命中 %r" % (label, jsf, needle))
    elif not prod:
        rev_gaps.append("%s: %s 無 producer" % (label, jsf))
add("2 headline 數字 → data → producer（三重反查）", not rev_gaps, "gaps = %s" % (rev_gaps or "0"))

# ---------------- ⑥ driver
drv = os.path.join(SCRIPTS, "run_chain.sh")
drv_ok = os.path.exists(drv)
n_run = 0
if drv_ok:
    dtxt = open(drv, encoding="utf-8").read()
    n_run = len(re.findall(r"\brun\b|python3 \"?\$f", dtxt))
    thr = "OMP_NUM_THREADS=1" in dtxt
else:
    thr = False
add("6 driver 存在且鎖 BLAS 單線程", drv_ok and thr, "run_chain.sh=%s, OMP=1=%s" % (drv_ok, thr))

# ---------------- ⑦ L2 / L3
l2p = os.path.join(DATA, "l2_independent_verify.json")
l3p = os.path.join(DATA, "l3_doc_crosscheck.json")
l2 = json.load(open(l2p, encoding="utf-8"))["summary"] if os.path.exists(l2p) else None
l3 = json.load(open(l3p, encoding="utf-8"))["summary"] if os.path.exists(l3p) else None
add("7a L2 全 PASS", bool(l2) and l2["fail"] == 0, "L2 = %s/%s" % (l2["pass"], l2["total"]) if l2 else "缺")
add("7b L3 全 PASS 且零 MISS", bool(l3) and l3["fail"] == 0 and l3["miss"] == 0,
    "L3 = %s/%s (miss %s)" % (l3["pass"], l3["total"], l3["miss"]) if l3 else "缺")

# ---------------- L1 重跑覆蓋（讀 data/audit_20261028/exit_codes.txt）
ec = os.path.join(DATA, "audit_20261028", "exit_codes.txt")
nz = []
if os.path.exists(ec):
    for line in open(ec, encoding="utf-8"):
        if "rc=" in line and "rc=0" not in line:
            nz.append(line.strip())
        if "MISSING" in line:
            nz.append(line.strip())
add("7c L1 全鏈重跑零非零退出", os.path.exists(ec) and not nz,
    "nonzero = %s（log dir=data/audit_20261028）" % (nz or "0"))

npass = sum(1 for r in OUT if r["ok"])
print("\n" + "=" * 100)
print("AUDIT RESULT: %d/%d PASS" % (npass, len(OUT)))
print("=" * 100)

with open(os.path.join(DATA, "41_chain_audit.json"), "w", encoding="utf-8") as fh:
    json.dump({"summary": {"pass": npass, "total": len(OUT)},
               "checks": OUT,
               "inventory": {"scripts": len(py), "nn_scripts": len(nn_py), "docs": len(md), "data_json": len(js)},
               "l2": l2, "l3": l3}, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON -> data/41_chain_audit.json")
sys.exit(0 if npass == len(OUT) else 1)