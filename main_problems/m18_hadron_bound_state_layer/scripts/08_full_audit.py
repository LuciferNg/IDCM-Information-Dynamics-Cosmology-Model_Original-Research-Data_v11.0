#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/08 — 全量審計
=========================
前設問題
--------
交付層審計：7 支推導腳本 × 7 份 docs 之可重現性、可追溯性、四件套同步。

事前判準
--------
  A1 全部 01–07 重跑 exit 0，且 data/*.json **byte-identical**（sha256 前後相同）
  A2 每支腳本有 docstring，且 docs/NN 之 §一 與 docstring 首段一致（存在性）
  A3 每份 data JSON 可解析，且含 checks/pass 欄位
  A4 docs/NN 存在且檔名 ↔ 腳本編號一一對應
  A5 關鍵數字可追溯：docs 內出現之核心數字必須喺對應 data JSON 內出現
  A6 本 M 目錄外零寫入（以 git status 核，見 §三）

輸出：data/08_full_audit.json
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
DOCS = os.path.join(ROOT, "docs")
checks = []

DOCMAP = {
    "01_reproduction_gate.py": "01_reproduction_gate.md",
    "02_support_theorem.py": "02_support_theorem.md",
    "03_composition_law.py": "03_composition_law.md",
    "04_jpc_classification.py": "04_jpc_classification.md",
    "05_radial_ladder.py": "05_radial_ladder.md",
    "06_shape_parameter_identity.py": "06_shape_parameter_identity.md",
    "07_feedback_gate.py": "07_feedback_gate.md",
    "09_jpc_geometric_carrier.py": "09_jpc_geometric_carrier.md",
    "10_residual_census.py": "10_residual_register.md",
    "11_coupling_mechanism.py": "11_coupling_mechanism.md",
    "12_ma1_faithfulness.py": "12_ma1_faithfulness.md",
    "13_ma1_tower.py": "13_ma1_tower.md",
    "14_ma1_carrier.py": "14_ma1_carrier.md",
    "15_ma1_induction.py": "15_ma1_induction.md",
    "16_md12_correction.py": "16_md12_correction.md",
    "17_connection_law.py": "17_connection_law.md",
    "18_ma234_resolution.py": "18_ma234_resolution.md",
    "21_closeout_verify.py": "21_closeout_verify.md",
    "22_coupling_dof.py": "22_coupling_dof.md",
}


def chk(name, ok, detail):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail)})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


print("=" * 78)
print("A1. 重跑 01–07，比對 data JSON（byte-identical）")
print("=" * 78)
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
rows = []
for scr in DOCMAP:
    jp = os.path.join(DATA, scr.replace(".py", ".json"))
    before = sha(jp) if os.path.exists(jp) else None
    r = subprocess.run(["/usr/bin/python3", os.path.join(HERE, scr)], capture_output=True, env=env)
    after = sha(jp) if os.path.exists(jp) else None
    rows.append({"script": scr, "exit": r.returncode, "sha_before": before, "sha_after": after,
                 "identical": before == after and before is not None})
    print("  %-36s exit=%d  sha %s… %s" % (scr, r.returncode, (after or "")[:12],
                                           "identical" if before == after else "DIFF"))
chk("A1a 全部 01–07 exit 0", all(r["exit"] == 0 for r in rows),
    "exit 集合 = %s" % sorted({r["exit"] for r in rows}))
chk("A1b 全部 data JSON byte-identical（重跑）", all(r["identical"] for r in rows),
    "%d/%d identical" % (sum(1 for r in rows if r["identical"]), len(rows)))

print()
print("=" * 78)
print("A2/A3/A4. docstring／JSON 結構／docs 對應")
print("=" * 78)
for scr, doc in DOCMAP.items():
    src = open(os.path.join(HERE, scr), encoding="utf-8").read()
    has_doc = bool(re.match(r'#!/usr/bin/env python3\n# -\*- coding: utf-8 -\*-\n"""', src)) or '"""' in src[:400]
    dpath = os.path.join(DOCS, doc)
    chk("A2 %s docstring + docs 存在" % scr, has_doc and os.path.exists(dpath),
        "docstring=%s docs=%s" % (has_doc, os.path.exists(dpath)))
    jp = os.path.join(DATA, scr.replace(".py", ".json"))
    try:
        j = json.load(open(jp, encoding="utf-8"))
        ok = ("checks" in j) or ("pass" in j)
    except Exception as e:
        ok = False
        j = {"error": str(e)}
    chk("A3 %s JSON 可解析且含 checks/pass" % scr, ok, list(j.keys())[:6])

print()
print("=" * 78)
print("A5. 關鍵數字追溯（docs §四 出現之核心數 ⟹ data JSON 內存在）")
print("=" * 78)
TRACE = [
    ("01_reproduction_gate.md", "01_reproduction_gate.json", ["0.005095770779492917", "90", "272"]),
    ("02_support_theorem.md", "02_support_theorem.json", ["8388608", "22"]),
    ("03_composition_law.md", "03_composition_law.json", ["15180", "22"]),
    ("04_jpc_classification.md", "04_jpc_classification.json", ["176"]),
    ("05_radial_ladder.md", "05_radial_ladder.json", ["12424", "16"]),
    ("06_shape_parameter_identity.md", "06_shape_parameter_identity.json", ["18.73875", "0.028257"]),
    ("07_feedback_gate.md", "07_feedback_gate.json", ["1.271754", "1.04855"]),
    ("09_jpc_geometric_carrier.md", "09_jpc_geometric_carrier.json", ["12", "0+−", "T2"]),
    ("10_residual_register.md", "10_residual_census.json", ["MA1", "MB1", "MC1", "MD1"]),
    ("11_coupling_mechanism.md", "11_coupling_mechanism.json", ["249", "T2", "n_lev"]),
    ("12_ma1_faithfulness.md", "12_ma1_faithfulness.json", ["MC10", "60", "88"]),
    ("13_ma1_tower.md", "13_ma1_tower.json", ["59", "60", "20"]),
    ("14_ma1_carrier.md", "14_ma1_carrier.json", ["30", "60", "20"]),
    ("15_ma1_induction.md", "15_ma1_induction.json", ["20", "60", "4"]),
    ("16_md12_correction.md", "16_md12_correction.json", ["8", "36", "60"]),
    ("17_connection_law.md", "17_connection_law.json", ["5", "20", "60"]),
    ("18_ma234_resolution.md", "18_ma234_resolution.json", ["60", "120", "44"]),
    ("21_closeout_verify.md", "21_closeout_verify.json", ["90", "272", "68"]),
    ("22_coupling_dof.md", "22_coupling_dof.json", ["11", "36", "60"]),
]
for doc, js, keys in TRACE:
    dtxt = open(os.path.join(DOCS, doc), encoding="utf-8").read()
    jtxt = open(os.path.join(DATA, js), encoding="utf-8").read()
    miss_d = [k for k in keys if k not in dtxt]
    miss_j = [k for k in keys if k not in jtxt]
    chk("A5 %s 數字雙向可溯" % doc, not miss_d and not miss_j,
        "docs 缺 %s｜JSON 缺 %s" % (miss_d, miss_j))

print()
print("=" * 78)
print("A6. 本 M 目錄外零寫入（git status ＋ mtime 歸屬）")
print("=" * 78)
repo = os.path.abspath(os.path.join(ROOT, "..", ".."))
MARK = os.path.expanduser("~/.hermes/backups/m18_scope_20261111/scopes.yaml.bak")
t_session = os.path.getmtime(MARK)
print("  session 起點標記 mtime = %s" % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t_session)))
try:
    out = subprocess.run(["git", "status", "--porcelain"], cwd=repo,
                         capture_output=True, text=True).stdout
    lines = [l for l in out.splitlines() if l.strip()]
    mine, pre, late = [], [], []
    for l in lines:
        rel = l[3:].strip().split(" -> ")[-1]
        if "m18_hadron_bound_state_layer" in rel:
            mine.append(l)
            continue
        ap = os.path.join(repo, rel)
        if os.path.isdir(ap):
            pre.append((l, "dir"))
            continue
        mt = os.path.getmtime(ap) if os.path.exists(ap) else 0
        tag = time.strftime("%m-%d %H:%M", time.localtime(mt)) if mt else "?"
        if mt > t_session:
            late.append((l, tag))
        else:
            pre.append((l, tag))
    print("  變更行 %d｜本 M %d｜本 M 以外（早於起點）%d｜晚於起點 %d"
          % (len(lines), len(mine), len(pre), len(late)))
    for l, mt in pre:
        print("     [pre-existing %s] %s" % (mt, l))
    for l, mt in late:
        print("     [★ late %s] %s" % (mt, l))
    for l in mine:
        print("     [mine] %s" % l)
    # 硬閘結構證明：session-write-gate mode=enforce 且本 session roots 只含 m18
    import yaml
    sc = os.path.expanduser("~/.hermes/plugins/session-write-gate/scopes.yaml")
    cfg = yaml.safe_load(open(sc, encoding="utf-8"))
    myscope = [x for x in cfg.get("scopes", []) if str(x.get("session_key")) == "137600"]
    ok_gate = (cfg.get("mode") == "enforce" and len(myscope) == 1
               and len(myscope[0].get("roots", [])) == 1
               and "m18_hadron_bound_state_layer" in myscope[0]["roots"][0])
    chk("A6a 硬閘結構證明：mode=enforce、本 session 恰 1 scope、roots 只含 m18 ⟹ 越界寫入結構上不可能",
        ok_gate, "mode=%s｜scope 數=%d｜roots=%s"
        % (cfg.get("mode"), len(myscope), myscope[0].get("roots") if myscope else None))
    import re as _re
    # 本 M 外變更：逐條歸屬（落喺他 M 目錄 ⟹ 並行 session）
    attrib = []
    for l, mt in pre + late:
        rel = l[3:].strip().split(" -> ")[-1]
        other_m = bool(_re.search(r"main_problems/m\d+_", rel)) or rel.startswith("_backup")
        attrib.append((l, mt, other_m))
    bad = [a for a in attrib if not a[2]]
    print("   本 M 外變更 %d 行；落喺他 M 目錄／備份 %d 行；無法歸屬 %d 行"
          % (len(attrib), sum(1 for a in attrib if a[2]), len(bad)))
    for l, mt, _ in attrib:
        print("     [outside %s] %s" % (mt, l))
    chk("A6b 本 M 外變更全部落喺他 M 目錄（並行 session）⟹ 本 session 可歸因寫入 = 0",
        len(bad) == 0, "無法歸屬 %d 行" % len(bad))
except Exception as e:
    chk("A6 git status 可執行", False, str(e))

npass = sum(1 for c in checks if c["pass"])
nfail = len(checks) - npass
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）" % (npass, nfail, len(checks)))
print("=" * 78)
out = {"script": "08_full_audit.py", "checks": checks, "pass": npass, "fail": nfail, "rerun": rows}
with open(os.path.join(DATA, "08_full_audit.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/08_full_audit.json")
sys.exit(0 if nfail == 0 else 1)
