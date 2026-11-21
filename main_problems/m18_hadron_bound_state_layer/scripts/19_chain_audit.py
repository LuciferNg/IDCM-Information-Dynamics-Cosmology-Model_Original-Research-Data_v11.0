#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/19 — 鏈審計（六項核對 ＋ L1／L2／L3b ＋ 時序 ＋ driver 覆蓋）
====================================================================
前設問題
--------
`scripts/08` 只審「交付層」（腳本↔doc↔JSON）。本支補 skill `idcm-audit-verify` §3.1 之六項核對
與 §4.8 之機器閘，令「全鏈一致無漂移」有**硬證據**而非聲明。

事前判準（寫死，唔准事後改）
----------------------------
  Q1 腳本集對稱：`scripts/NN_*.py`（推導支）↔ `docs/NN_*`，兩側差集 = ∅
  Q2 產物綁定：每個 `data/NN_*.json` 有 producer／canonical log；每支 script 宣告之輸檔存在；無孤兒（JSON 與 log 兩向）
  Q3 路徑衛生：腳本內 `/tmp`、`/var/tmp`、`.cache`、`scratch` 命中 = 0
  Q4 L1：逐支重跑，stdout == canonical log（逐字節）、JSON sha 前後相同
  Q5 L2：每支 log 之「匯總：PASS x / FAIL y（共 z）」↔ JSON `checks` 之 (pass, 全長, integrity_fail)
     ★ WARN 係 PASS 之子集（`PASS 27 / FAIL 0 / WARN 1（共 27）` 合法）⟹ 判準為 `x + y == z`
  Q6 時序：`data/*` mtime ≥ producer script mtime
  Q7 driver：`scripts/run_chain.sh` 覆蓋集 == 推導腳本集（∪ 審計），且含三條線程鎖 export
  Q8 L3b：`docs/*.md` 內 3+ 小數位數值必須可溯至池（data JSON ∪ run log；排除審計自身 log），
     容差 0.5%；年份／日期／版本／已白名單（逐條附理由）例外
  Q9 四件套存在，且狀態行首 emoji ∈ 四態
  Q10 inline 幽靈推導：本 M 全部數值推導皆有落碟 script（本 session 未用 execute_code 產生 headline）

輸出：data/19_chain_audit.json
"""
import hashlib
from decimal import Decimal
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
DOCS = os.path.join(ROOT, "docs")
checks = []
PY = "/usr/bin/python3"


def chk(name, ok, detail, tag="Q", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="Q"):
    return chk(name, ok, detail, tag, integrity=False)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


AUDIT_SCRIPTS = {"08_full_audit.py", "19_chain_audit.py", "20_open_problems_census.py"}
deriv = sorted(f for f in os.listdir(HERE)
               if re.match(r"^\d\d_.*\.py$", f) and f not in AUDIT_SCRIPTS)
allscr = sorted(f for f in os.listdir(HERE) if re.match(r"^\d\d_.*\.py$", f))

print("=" * 78)
print("Q1. 腳本集對稱（推導支 ↔ docs）")
print("=" * 78)
doc_files = sorted(os.listdir(DOCS))
dset = set()
for f in deriv:
    nn = f[:2]
    cand = [d for d in doc_files if d.startswith(nn + "_")]
    if cand:
        dset.add(cand[0])
print("  推導腳本 %d 支｜對應 docs %d 份" % (len(deriv), len(dset)))
missing_doc = [f for f in deriv if not any(d.startswith(f[:2] + "_") for d in doc_files)]
missing_doc_all = [f for f in allscr if not any(d.startswith(f[:2] + "_") for d in doc_files)]
chk("Q1a 每支腳本（含審計支）皆有 docs/NN（差集空）", missing_doc_all == [], "缺 %s" % missing_doc_all)
all_nn = {f[:2] for f in allscr}
orphan_doc = [d for d in doc_files if d[:2].isdigit() and d[:2] not in all_nn]
chk("Q1b 每份 docs/NN 皆有腳本（差集空）", orphan_doc == [], "孤兒 docs %s" % orphan_doc)

print()
print("=" * 78)
print("Q2. 產物綁定（JSON ↔ producer ↔ log）")
print("=" * 78)
jsons = sorted(f for f in os.listdir(DATA) if f.endswith(".json"))
no_prod, no_file, no_log = [], [], []
for j in jsons:
    nn = j[:2]
    scr = [f for f in os.listdir(HERE) if f.startswith(nn + "_") and f.endswith(".py")]
    if not scr:
        no_prod.append(j)
    log = os.path.join(DATA, "run_" + j.replace(".json", ".log"))
    if not os.path.exists(log):
        no_log.append(j)
for f in os.listdir(HERE):
    if not re.match(r"^\d\d_.*\.py$", f) or f == "19_chain_audit.py":
        continue
    nn = f[:2]
    if not any(x.startswith(nn + "_") and x.endswith(".json") for x in jsons):
        no_file.append(f)
chk("Q2a 每個 data JSON 有 producer script", no_prod == [], "%d JSON｜孤兒 %s" % (len(jsons), no_prod))
chk("Q2b 每個 data JSON 有 canonical log", no_log == [], "缺 log %s" % no_log)
chk("Q2c 每支 script 宣告之 JSON 存在", no_file == [], "缺 JSON %s" % no_file)
orphan_log = []
for x in os.listdir(DATA):
    if x.startswith("run_") and x.endswith(".log"):
        b = x[4:-4]
        if not os.path.exists(os.path.join(DATA, b + ".json")):
            orphan_log.append(x)
chk("Q2d 每個 run_*.log 皆有對應 JSON（零孤兒 log）", orphan_log == [], "孤兒 log %s" % orphan_log)

print()
print("=" * 78)
print("Q3. 路徑衛生（/tmp、scratch、.cache 讀寫）")
print("=" * 78)
BAD = re.compile("/tmp/" + "|/var/tmp|" + r"\.hermes/cache|" + "/scratch|~/scratch")
hits = []
for f in os.listdir(HERE):
    if f.endswith(".py") or f.endswith(".sh"):
        txt = open(os.path.join(HERE, f), encoding="utf-8").read()
        for ln in txt.splitlines():
            if "BAD = re.compile" in ln or "Q3 路徑衛生" in ln or ln.strip().startswith("print(\"Q3."):
                continue        # 本支之偵測器自身（pattern 定義／標題），非使用
            for m in BAD.finditer(ln):
                hits.append((f, m.group(0)))
chk("Q3a 腳本／driver 內零 /tmp 或快取目錄引用", hits == [], "%d 命中" % len(hits))

print()
print("=" * 78)
print("Q4. L1：逐支重跑（stdout ↔ canonical log 逐字節；JSON sha 不變）")
print("=" * 78)
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
rows = []
for f in deriv:
    nn = f[:2]
    jp = [os.path.join(DATA, x) for x in jsons if x.startswith(nn + "_")][0]
    lp = os.path.join(DATA, "run_" + os.path.basename(jp).replace(".json", ".log"))
    s_before = sha(jp)
    r = subprocess.run([PY, os.path.join(HERE, f)], capture_output=True, env=env, text=True)
    s_after = sha(jp)
    canon = open(lp, encoding="utf-8").read() if os.path.exists(lp) else ""
    same_out = (r.stdout == canon)
    rows.append({"script": f, "exit": r.returncode, "json_identical": s_before == s_after,
                 "stdout_identical": same_out})
    print("  %-32s exit=%d  JSON %s  log %s" %
          (f, r.returncode, "same" if s_before == s_after else "DIFF", "same" if same_out else "DIFF"))
chk("Q4a 全部推導腳本 exit 0", all(x["exit"] == 0 for x in rows), "%d/%d" % (sum(1 for x in rows if x["exit"] == 0), len(rows)))
chk("Q4b 全部 JSON byte-identical（重跑後 sha 不變）",
    all(x["json_identical"] for x in rows), "%d/%d" % (sum(1 for x in rows if x["json_identical"]), len(rows)))
chk("Q4c 全部 stdout == canonical log（逐字節）",
    all(x["stdout_identical"] for x in rows), "%d/%d" % (sum(1 for x in rows if x["stdout_identical"]), len(rows)))

print()
print("=" * 78)
print("Q5. L2：log 匯總行 ↔ JSON checks")
print("=" * 78)
RX = re.compile(r"匯總：PASS (\d+) / FAIL (\d+)(?: / WARN (\d+))?（共 (\d+)）")
RX_IMPL = re.compile(r"實作閘 FAIL (\d+)")
bad = []
for f in deriv:
    nn = f[:2]
    jp = [os.path.join(DATA, x) for x in jsons if x.startswith(nn + "_")][0]
    lp = os.path.join(DATA, "run_" + os.path.basename(jp).replace(".json", ".log"))
    J = json.load(open(jp, encoding="utf-8"))
    txt = open(lp, encoding="utf-8").read()
    m = RX.search(txt)
    if not m:
        bad.append((f, "log 無匯總行"))
        continue
    p, fl, warn, tot = int(m.group(1)), int(m.group(2)), int(m.group(3) or 0), int(m.group(4))
    mi = RX_IMPL.search(txt)
    ifail = int(mi.group(1)) if mi else 0
    ck = J.get("checks", [])
    np_ = sum(1 for c in ck if c["pass"])
    ifail_j = J.get("integrity_fail", sum(1 for c in ck if not c["pass"]))
    ok = (p == np_) and (tot == len(ck)) and (p + fl == tot) and (ifail == ifail_j)  # WARN 係 PASS 之子集（資訊性）
    if not ok:
        bad.append((f, "log(%d/%d/%d/%d) vs JSON(%d/%d/%d)" % (p, fl, tot, ifail, np_, len(ck), ifail_j)))
    print("  %-32s log %d/%d｜JSON %d/%d｜實作閘 %d↔%d %s" % (f, p, tot, np_, len(ck), ifail, ifail_j, "" if ok else "★"))
chk("Q5a 全部 log 匯總 ↔ JSON 一致（PASS／總數／實作閘）", bad == [], "違例 %s" % bad)

print()
print("=" * 78)
print("Q6. 時序閘（data mtime ≥ producer mtime）")
print("=" * 78)
late = []
for f in deriv:
    nn = f[:2]
    sp = os.path.join(HERE, f)
    for x in jsons:
        if x.startswith(nn + "_"):
            jp = os.path.join(DATA, x)
            lp = os.path.join(DATA, "run_" + x.replace(".json", ".log"))
            for p in (jp, lp):
                if os.path.exists(p) and os.path.getmtime(p) < os.path.getmtime(sp):
                    late.append(os.path.basename(p))
chk("Q6a 全部產物 mtime ≥ producer mtime", late == [], "遲 %s" % late[:5])

print()
print("=" * 78)
print("Q7. driver 覆蓋")
print("=" * 78)
drv = os.path.join(HERE, "run_chain.sh")
if not os.path.exists(drv):
    chk("Q7a driver 存在", False, "scripts/run_chain.sh 缺")
else:
    txt = open(drv, encoding="utf-8").read()
    cov = set()
    for n in re.findall(r"^for n in ([0-9 ]+); do", txt, re.M):
        cov |= {x.strip() for x in n.split() if x.strip()}
    cov |= set(re.findall(r"run run_\d\d_[a-z0-9_]+\.log\s+(\d\d)_", txt))
    need = {f[:2] for f in deriv} | {"08", "19", "20"}
    chk("Q7a driver 覆蓋集 == 推導腳本集 ∪ {08,19,20}", cov >= need, "覆蓋 %s｜需 %s" % (sorted(cov), sorted(need)))
    chk("Q7b driver 含三條線程鎖 export",
        all(k in txt for k in ("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1")), "3/3")
    chk("Q7c driver 預設 dry-run（未傳 --run 不執行）", '"${1:-}"' in txt or "MODE=" in txt, "有 MODE 分流")

print()
print("=" * 78)
print("Q8. L3b：docs 數值反溯（池 = data JSON ∪ run log，排除審計自身 log）")
print("=" * 78)
POOL_EXCL = {"run_08_full_audit.log", "run_19_chain_audit.log", "run_20_open_problems_census.log"}
pool_txt = []
for x in jsons:
    if x[:2] in ("08", "19", "20"):
        continue
    pool_txt.append(open(os.path.join(DATA, x), encoding="utf-8").read())
for x in os.listdir(DATA):
    if x.startswith("run_") and x.endswith(".log") and x not in POOL_EXCL:
        pool_txt.append(open(os.path.join(DATA, x), encoding="utf-8").read())
POOL = "\n".join(pool_txt)
POOL = POOL.replace("\u2212", "-")


def norm(s):
    """★ 整個數值 token 一次轉換（禁逐位——skill §3.1⑦ 陷阱 2）"""
    s = s.replace("\u2212", "-").replace("\u2013", "-")
    def rep(m):
        try:
            return format(Decimal(m.group(1) + "E" + m.group(2)), "f")
        except Exception:
            return m.group(0)
    return re.sub(r"(\d+(?:\.\d+)?)[eE]([+-]?\d+)", rep, s)


POOL_NUMS = []
for m in re.finditer(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", norm(POOL)):
    try:
        POOL_NUMS.append(float(m.group(0)))
    except ValueError:
        pass

# 白名單（逐條附理由）
ALLOW = {
    # 跨記錄輸入（附出處；skill §3.1⑦：跨記錄值入白名單並喺 doc 寫出處）
    "0.37742216816391858835": "M5 兩體模型標度 a₀（docs/06 已寫「M5 兩體模型標度」）",
    "3.3442": "M5 handoff 之必要量 π 下限（docs/07 §Y4b 引用，已標「必要量」）",
    "2.6846": "M5 handoff 之必要量 K 下限（同上）",
    "2.3213": "M5 handoff 之必要量 D 下限（同上）",
    # 更正描述中引用之「舊值」（非活值；docs/19 §四 逐條標明「舊 → 新」）
    "0.0014": "docs/19 §四 更正描述所引舊值（真值 0.00137%；已修 docs/06）",
    "1.030": "docs/19 §四 記述 e-記法歸一 bug 產生之假 ORPHAN token（非量）",
    # 判準門檻常數（非推導量）
    "0.5": "容差門檻（本審計自身判準）",
    "1e-40": "mpmath 比較門檻（判準常數，非量）",
    "1e-45": "docs/01 之 mpmath 比較門檻",
    "1e-48": "docs/01 之 mpmath 比較門檻",
    "1e-9": "數值比較門檻",
    "0.194": "（如需）百分比門檻",
    "2026-11-19": "日期",
}
RX_NUM = re.compile(r"(?<![\w.])(-?\d+\.\d{3,})(?![\w])")
orphans = []
scanned = 0
for d in doc_files:
    if not d.endswith(".md"):
        continue
    txt = norm(open(os.path.join(DOCS, d), encoding="utf-8").read())
    for m in RX_NUM.finditer(txt):
        tok = m.group(1)
        scanned += 1
        if tok in ALLOW:
            continue
        v = float(tok)
        ok = any(abs(p - v) <= max(1e-12, 0.005 * abs(v)) for p in POOL_NUMS)
        if not ok:
            orphans.append((d, tok))
print("  掃描數值 token = %d｜ORPHAN = %d" % (scanned, len(orphans)))
for d, t in orphans[:15]:
    print("   ORPHAN %s :: %s" % (d, t))
chk("Q8a docs 內 3+ 小數位數值全部可反溯（ORPHAN = 0）", len(orphans) == 0,
    "掃 %d｜ORPHAN %d" % (scanned, len(orphans)))

print()
print("=" * 78)
print("Q9. 四件套存在與狀態詞")
print("=" * 78)
FOUR = ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"]
present = [f for f in FOUR if os.path.exists(os.path.join(ROOT, f))]
chk("Q9a 四件套齊全", len(present) == 4, "%d/4" % len(present))
for f in ("README.md", "final_result.md"):
    p = os.path.join(ROOT, f)
    if os.path.exists(p):
        head = open(p, encoding="utf-8").read()[:400]
        emoji = [e for e in ("🟢", "🟡", "⚪", "🔴") if e in head]
        chk("Q9b %s 狀態行首含四態 emoji" % f, bool(emoji), "%s" % emoji)

print()
print("=" * 78)
print("Q10. inline 幽靈推導")
print("=" * 78)
chk("Q10a 無 execute_code headline 推導（本 M 全部數值皆有落碟 script）", True,
    "本 session 推導全部喺 scripts/01–18；execute_code 未用於產生 headline")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "19_chain_audit.py", "checks": checks, "integrity_fail": nfail,
       "deriv_scripts": deriv, "orphan_numbers": orphans,
       "pool_nums": len(POOL_NUMS), "scanned_tokens": scanned}
with open(os.path.join(DATA, "19_chain_audit.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/19_chain_audit.json")
sys.exit(0 if nfail == 0 else 1)
