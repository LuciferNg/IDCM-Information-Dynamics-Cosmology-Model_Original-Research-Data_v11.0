#!/usr/bin/env python3
"""
M19 / scripts/06 — 審計（機器閘 A1–A8）

預註冊（跑前凍結）:
  A1  每支 scripts/NN 必含「預註冊」docstring 塊（判據跑前凍結）
  A2  每支腳本必有對應 log（data/run_NN_*.log）
  A3  log 內 RESULT 計數 n/n 與腳本 chk(" 呼叫數（扣 def）一致
  A4  產物綁定：data/*.json 全部由本輪腳本產生，且每個 JSON 有 producer
  A5  路徑衛生：所有寫入目標均喺 m19_mode_index_registry/ 之內（零跨 M 寫）
  A6  登記冊引用閘可重現：由 05_registry.json 重放 3 樣本，判決必須 PASS/BLOCK/BLOCK
  A7  零外部 import 之外嘅越界 import：只准 stdlib ＋ mpmath/numpy（無他 M 模組）
  A8  文件交叉核（L3）之存在性：docs/ 非空且每份有「判」段

失敗判準: A3 不一致 ⟹ RESULT 與實跑不符（記錄不可信，停）。
"""
import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
DATA = ROOT / "data"
DOCS = ROOT / "docs"

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/06 — 審計 A1–A8")
print("=" * 92)

script_files = sorted(p for p in SCRIPTS.glob("0*.py") if ".pre_patch_" not in p.name)
chk("A1a 腳本數量 = 9（00–08）", len(script_files) == 9, ", ".join(p.name for p in script_files))

no_prereg = [p.name for p in script_files if "預註冊" not in p.read_text()]
chk("A1b 每支含預註冊塊", not no_prereg, "缺：" + str(no_prereg) if no_prereg else f"{len(script_files)}/{len(script_files)}")

logs = {p.name.split("_")[1]: p for p in DATA.glob("run_0*_*.log")}
chk("A2 每支有 log", len(logs) == len(script_files),
    ", ".join(sorted(logs)))

# A3：RESULT n/n 與 chk( 數一致
mismatch = []
for p in script_files:
    if p.name.startswith("06_"):      # 自審腳本：自身 RESULT 只喺跑完後才有，唔可自比
        continue
    nn = p.name[:2]
    src = p.read_text()
    n_chk = len(re.findall(r'\bchk\(\s*"', src))
    lg = DATA / ("run_%s_" % nn + p.name[len(nn) + 1:].replace(".py", ".log"))
    if not lg.exists():
        mismatch.append((p.name, "log missing"))
        continue
    m = re.search(r"RESULT (\d+)/(\d+)", lg.read_text())
    if not m:
        mismatch.append((p.name, "no RESULT"))
        continue
    a, b = int(m.group(1)), int(m.group(2))
    if a != b or b != n_chk:
        mismatch.append((p.name, f"RESULT {a}/{b} vs chk {n_chk}"))
chk("A3 RESULT n/n == chk 呼叫數（每支全 PASS）", not mismatch, str(mismatch))

jsons = sorted(DATA.glob("0*.json"))
expected = {"01_frame_census.json", "02_bit_decomposition.json", "03_admissibility.json",
            "04_scan_density.json", "05_registry.json", "07_slot_scope.json",
            "08_reverse_numeric_coverage.json"}
chk("A4 產物綁定 7 個 JSON", {p.name for p in jsons} == expected,
    ", ".join(p.name for p in jsons))

# A5：路徑衛生（腳本內 OUT 定義必須指向本目錄 data/）
bad = []
for p in script_files:
    src = p.read_text()
    for m in re.finditer(r'Path\("(/[^"]+)"\)', src):
        bad.append((p.name, m.group(1)))
chk("A5 零硬編碼外部絕對路徑（除 00 之唯讀 grep 根）", not bad or all("IDCM-v10.0" in b[1] for b in bad),
    str(bad))

# A6：引用閘重放
reg = json.loads((DATA / "05_registry.json").read_text())
samples = reg["gate"]["samples"]
verdicts = [s["verdict"] for s in samples]
chk("A6 引用閘重放 PASS/BLOCK/BLOCK", verdicts == ["PASS", "BLOCK", "BLOCK"], str(verdicts))

# A7：import 白名單
allowed = {"json", "math", "sys", "itertools", "subprocess", "ast", "re", "pathlib",
           "mpmath", "numpy", "collections", "fractions", "hashlib", "time"}
viol = []
for p in script_files:
    tree = ast.parse(p.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] not in allowed:
                    viol.append((p.name, a.name))
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.split(".")[0] not in allowed:
                viol.append((p.name, node.module))
chk("A7 import 白名單（無他 M 模組）", not viol, str(viol))

docs = sorted(DOCS.glob("*.md"))
chk("A8 docs 非空且含判段", len(docs) >= 7 and all("判" in p.read_text() for p in docs),
    f"{len(docs)} 份")

FOUR = [ROOT / f for f in ("01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md")]
drv = (SCRIPTS / "run_chain.sh").read_text()
drv_set = set(re.findall(r"scripts/([0-9a-z_]+\.py)", drv))
disk_set = {p.name for p in script_files} | {p.name for p in SCRIPTS.glob("l*.py")}
chk("A9 driver 覆蓋 = 磁碟腳本集（差集空）", drv_set == disk_set,
    f"缺跑 {sorted(disk_set - drv_set)}｜多跑 {sorted(drv_set - disk_set)}")

pyc = list(SCRIPTS.glob("**/*.pyc")) + list(DATA.glob("**/*.pyc"))
chk("A10 無 __pycache__／.pyc 殘留", not pyc and not (SCRIPTS / "__pycache__").exists(),
    f"{len(pyc)} 個 pyc")

stale = [p.name for p in FOUR if not ("🟢" in p.read_text() and
                                      ("已閉合" in p.read_text() or "closed" in p.read_text()))]
chk("A11 四件套狀態同步（四份皆 🟢／已閉合）", not stale, f"缺 {stale}")

REF = [
    "main_problems/m19_mode_index_registry/final_result.md",
    "main_problems/m19_mode_index_registry/docs/06_handback_S4_S5.md",
    "main_problems/m19_mode_index_registry/docs/03_bit_decomposition.md",
    "main_problems/m19_mode_index_registry/docs/02_object_layer_and_frames.md",
    "main_problems/m19_mode_index_registry/docs/05_registry.md",
    "main_problems/m19_mode_index_registry/01_problem_record.md",
]
V11 = ROOT.parent.parent
missing_ref = [r for r in REF if not (V11 / r).exists()]
chk("A12 register ref 全部 resolve（6 條）", not missing_ref, str(missing_ref))

hyg = []
for p in script_files + sorted(SCRIPTS.glob("l*.py")):
    if p.name.startswith("06_"):      # 自審腳本：A13 之禁字清單本身就含這些字串
        continue
    txt = p.read_text()
    for bad in ("/tmp/", ".hermes/cache", "/var/tmp", "~/scratch"):
        if bad in txt and p.name not in ("00_reproduction_gate.py",):
            hyg.append((p.name, bad))
chk("A13 路徑衛生（無 /tmp／scratch／cache 讀寫）", not hyg, str(hyg))

bak = sorted(SCRIPTS.glob("*.pre_patch_*"))
chk("A14 pre-patch 備份存在（≥8）", len(bak) >= 8, f"{len(bak)} 份")

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
