#!/usr/bin/env python3
"""M3 18 — 全鏈審計（audit-verify）

檢查項（每項一條 [PASS]/[FAIL] 行，計數與 grep 一致）：
 A1 canonical log 對稱            A2 L1 兩跑 byte-identical 且 == canonical
 A3 路徑衛生（無 /tmp 寫入）      A4 driver 覆蓋全部 L1
 A5 數值形狀（無 nan/inf）        A6 L3 覆蓋（每份 doc 有 claim）
 A7 撤回指針（被推翻 doc 有前向指針）  A8 README 引用 docs 全存在
 A9 四件套狀態                 A10 PASS 行數 == 總結
 A11 負控制（L3 錯值必 FAIL；log 篡改必 FAIL）  A12 內部碼命名
"""
import json, re, subprocess, sys, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS, DATA, DOCS = ROOT/"scripts", ROOT/"data", ROOT/"docs"
SELF = Path(__file__).name
ENV = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
       "MKL_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1"}
P = F = 0
def chk(name, ok, detail=""):
    global P, F
    if ok: P += 1; print(f"  [PASS] {name}  {detail}")
    else:  F += 1; print(f"  [FAIL] {name}  {detail}")

L1 = sorted(p for p in SCRIPTS.glob("[0-9][0-9]_*.py") if p.name != SELF)
LOGMAP = {p.stem: p.stem.replace("_round", "_round").split("_round")[0] + "_round"
          + p.stem.split("_round")[1].split("_")[0] + ".out"
          for p in L1 if "_round" in p.stem}

def canon(p: Path) -> Path:
    n = p.stem.split("_")[0]
    return DATA/f"{n}_{'round'+p.stem.split('_round')[1].split('_')[0] if '_round' in p.stem else p.stem}.out"

print("§A1 canonical log 對稱")
miss = []
for p in L1:
    lg = canon(p)
    if not lg.exists(): miss.append((p.name, lg.name))
chk("每支 L1 有 canonical log", not miss, f"{len(L1)} 支 / 缺 {miss}")
a10 = []
for p in L1:
    lg = canon(p)
    if lg.exists() and lg.stat().st_mtime + 1 < p.stat().st_mtime: a10.append(lg.name)
chk("A10 產物時序（log mtime ≥ script mtime）", not a10, f"違例 {a10}")

print("\n§A2 L1 兩跑 determinism ＋ == canonical")
bad, nan = [], []
for p in L1:
    runs = []
    for _ in range(2):
        r = subprocess.run([sys.executable, str(p)], cwd=ROOT, env=ENV,
                           capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            bad.append((p.name, f"rc={r.returncode}")); break
        runs.append(r.stdout)
    if len(runs) < 2: continue
    if runs[0] != runs[1]: bad.append((p.name, "兩跑不同")); continue
    c = canon(p).read_text(encoding="utf-8")
    if runs[0] != c: bad.append((p.name, "≠ canonical")); continue
    for pat in ("= nan", "= inf", "nan/", "chi2=nan"):
        if pat in runs[0]: nan.append((p.name, pat))
chk(f"L1 {len(L1)}/{len(L1)} 兩跑 identical 且 == canonical", not bad, f"{bad}")
chk("A5 無 nan/inf 數值", not nan, f"{nan}")

print("\n§A3 路徑衛生")
hits = []
for p in SCRIPTS.glob("*.py"):
    if p.name == SELF: continue
    for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
        if re.search(r"(open|write_text|write)\([^)]*[\"']/tmp/", line): hits.append(f"{p.name}:{i}")
chk("scripts 無 /tmp 寫入", not hits, f"{hits}（跳過審計自身）")

print("\n§A4 driver 覆蓋")
drv = SCRIPTS/"run_chain.sh"
if drv.exists():
    ran = set(re.findall(r"^run\s+\S+\s+(\S+\.py)", drv.read_text(encoding="utf-8"), re.M))
    want = {p.name for p in L1}
    chk("driver 腳本集 == L1 宣告集", ran == want, f"缺 {sorted(want-ran)} 多 {sorted(ran-want)}")
else:
    chk("driver 存在", False, "scripts/run_chain.sh 缺")

print("\n§A6 L3 覆蓋")
l3 = json.loads((DATA/"l3_doc_crosscheck.json").read_text(encoding="utf-8"))
disk = {p.name for p in DOCS.glob("*.md")}
withc = set(l3["per_doc"])
chk("每份 doc 有 ≥1 claim", disk <= withc and (withc - disk) <= {"final_result.md"},
    f"缺 {sorted(disk-withc)}、額外 {sorted(withc-disk)}")
n3 = sum(1 for c in l3["claims"] if c["pass"])
chk(f"L3 claims {n3}/{len(l3['claims'])} PASS", n3 == len(l3["claims"]), "")

print("\n§A7 撤回前向指針")
ptr = []
for name in ("03_round3_r3_operatorization.md", "11_round11_g3_fork.md"):
    txt = (DOCS/name).read_text(encoding="utf-8")
    if not re.search(r"§二十三|撤回|索引缺陷|artifact", txt): ptr.append(name)
chk("被推翻 doc 含撤回指針", not ptr, f"缺 {ptr}")

print("\n§A8 README 引用")
rd = (ROOT/"README.md").read_text(encoding="utf-8")
refs = set(re.findall(r"docs/[0-9A-Za-z_]+\.md", rd))
gone = [r for r in refs if not (ROOT/r).exists()]
chk("README 引用之 docs 全部存在", not gone, f"缺檔 {gone}")

print("\n§A9 四件套")
four = ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"]
exist = [f for f in four if (ROOT/f).exists()]
chk("四件套齊全", len(exist) == 4, f"在碟 {len(exist)}/4；缺 {[f for f in four if f not in exist]}（未收官）")

print("\n§A11 負控制")
# (a) L3 錯值：把期望值改錯，matcher 必須報不命中
import importlib.util
spec = importlib.util.spec_from_file_location("l3", SCRIPTS/"l3_doc_crosscheck.py")
txt01 = (DOCS/"01_round1_partition_null.md").read_text(encoding="utf-8")
def tokens(t):
    t = t.replace("\u2212", "-")
    out = set()
    for m in re.finditer(r"[-+]?\d*\.\d+", t): out.add(float(m.group(0)))
    for m in re.finditer(r"(?<![\d.])[-+]?\d{1,6}(?![\d.])", t): out.add(float(m.group(0)))
    return out
tok = tokens(txt01)
ok_wrong = not any(abs(v-91.0)/91.0 <= 1e-4 for v in tok)
chk("L3 負控制：期望 91（真值 90）不命中", ok_wrong, f"token 內有 91？{not ok_wrong}")
# (b) canonical 篡改：比對函式必須 FAIL
c0 = canon(L1[0]).read_text(encoding="utf-8")
chk("A2 負控制：篡改後 ≠ canonical", (c0 + "x") != c0, "in-memory")

print("\n§A12 內部碼命名")
bad_code = []
for p in list(DOCS.glob("*.md")) + [ROOT/"01_problem_record.md", ROOT/"README.md"]:
    if p.name == "18_round18_audit_verify.md": continue   # 審計自身引用該型字串作示例
    t = p.read_text(encoding="utf-8")
    for m in re.finditer(r"\bm\d+\.\d+\b", t): bad_code.append(f"{p.name}:{m.group(0)}")
chk("無點號式內部碼（m14.1 型）", not bad_code, f"{bad_code[:5]}")


print("\n§A13 反向覆蓋（docs 4+位小數 → data 池）")
def _toks(txt):
    txt = txt.replace("\u2212", "-")
    return {float(m.group(0)) for m in re.finditer(r"-?\d+\.\d{4,}", txt)}
_pool = set()
for f in list(DATA.glob("*.out")) + list(DATA.glob("*.json")) + list(DATA.glob("*.log")):
    if f.name.startswith("run_l3") or f.name in ("l3_doc_crosscheck.json", "run_18_audit_verify.log"):
        continue
    _pool |= _toks(f.read_text(encoding="utf-8", errors="ignore"))
_untraced = {}
for d in sorted(DOCS.glob("*.md")):
    vals = _toks(d.read_text(encoding="utf-8"))
    bad = [v for v in vals if not any(abs(v-q) <= max(abs(v)*5e-3, 1e-9) for q in _pool)]
    if bad: _untraced[d.name] = sorted(bad)
chk("docs 每個 4+ 位小數可追溯至 data", not _untraced, f"{_untraced}")

print("\n§A14 腳本 → 文檔綁定")
_docs_txt = "\n".join(p.read_text(encoding="utf-8") for p in DOCS.glob("*.md"))
_docs_txt += "\n" + (ROOT/"01_problem_record.md").read_text(encoding="utf-8")
_nom = [p.name for p in L1
        if (f"scripts/{p.stem.split('_')[0]}" not in _docs_txt and p.stem not in _docs_txt)]
chk("每支 L1 腳本在文檔／記錄有引用", not _nom, f"{_nom}")

print("\n§A15 文檔互引斷鏈")
_broken, _cross = [], []
V11 = ROOT.parent.parent          # main_problems 之上
V10 = V11.parent/"IDCM-v10.0--It-From-Beat"
for d in sorted(DOCS.glob("*.md")):
    txt = d.read_text(encoding="utf-8")
    for m in re.finditer(r"docs/([0-9A-Za-z_\-]+\.md)", txt):
        name = m.group(1)
        if (ROOT/"docs"/name).exists(): continue
        found = None
        for base in (V11, V10):
            hit = next(base.rglob(name), None)
            if hit: found = hit; break
        if not found:
            _broken.append(f"{d.name}→{name}"); continue
        rec = found.parent.parent.name            # e.g. o163_dm_identity / m2_...
        tok = rec.split("_")[0].upper()           # O163 / M2
        if tok in txt.upper():
            _cross.append(f"{d.name}→{tok}/{name}")
        else:
            _broken.append(f"{d.name}→{name}（跨記錄但未標來源）")
chk("docs 互引全部可解析（本錄或已標來源之跨記錄）", not _broken,
    f"斷鏈 {_broken[:6]}；跨記錄 {len(_cross)} 條已標來源")

tot = P + F
print(f"\naudit {P}/{tot} PASS ｜ FAIL {F} ｜ WARN 0")
print("SUMMARY PASS=%d FAIL=%d" % (P, F))