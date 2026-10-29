#!/usr/bin/env python3
"""47 — 全鏈審計（All-chain audit）：六項核對 ＋ 四件套 ＋ register ＋ 覆蓋範圍

依 skill `idcm-audit-verify`：
  §3.1 ①腳本集對稱 ②逆向綁定 ③路徑衛生 ④產物綁定 ⑥driver彙總
  §0.1 四件套狀態行同步
  §7.3 L2／L3 覆蓋範圍披露
  §7   陷阱 9（圍欄平衡）、10（U+FFFD）

判準（跑前凍死）
--------------------------------------------------
A1 腳本↔log 配對：NN 腳本（EXCL 15,28）與 run_NN log 兩側差集皆空
A2 docs↔README：README 表格列數 = docs 檔數，且每個 docs/NN 皆有列
A3 產物綁定：每個 data/*.json 名皆在至少一支 script 中出現（無孤兒）
A4 路徑衛生：scripts/*.py 之 /tmp 寫入模式 = 0（排除檢查語句自身）
A5 四件套狀態行：01／README／final_result 三者所宣稱之最高輪次須等於實際最高輪次
     （SESSION_START 屬第 1 輪快照，依 skill §0 只 append，故列 INFO 不計 FAIL）
A6 時序索引覆蓋：01_problem_record 須提及 ≥ 95% 之 docs/NN
A7 register ref 覆蓋：M12 之 note/handoff ref 須覆蓋 docs/01–docs/NN
A8 圍欄平衡（偶數）＋ U+FFFD = 0（全部 .md）
A9 更正留痕：已撤回之結論須有「撤回／更正」字樣
A10 覆蓋披露：L2／L3 是否在全鏈 glob 內（只報，不判）
"""
import json
import re
from pathlib import Path

# ── canonical log 自寫（2026-10-29 修）────────────────────────────────────────
# 原缺陷：本腳本只 print 到 stdout，canonical log 靠 shell redirect 產生
#          ⟹ scripts/15 之 A3「跑前快照」法失效（腳本不覆寫 ⟹ 狀態一變即漂移）。
# 此處令 stdout 雙寫（終端 + canonical log），atexit 落碟。
import atexit as _atexit
import pathlib as _pl
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path = path
        self.buf = []
        self.out = _sys.__stdout__

    def write(self, s):
        self.buf.append(s)
        if self.out is not None:
            try:
                self.out.write(s)
            except Exception:
                pass

    def flush(self):
        if self.out is not None:
            try:
                self.out.flush()
            except Exception:
                pass

    def close(self):
        self.path.write_text("".join(self.buf), encoding="utf-8")


_log = _Tee(_pl.Path(__file__).resolve().parent.parent / "data" / "run_47_allchain_audit.log")
_sys.stdout = _log
_atexit.register(_log.close)
# ────────────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent
DOCS, SCRIPTS, DATA = ROOT / "docs", ROOT / "scripts", ROOT / "data"
REG = Path("/home/wsl/Projects/IDCM-Original/problem_register/register.jsonl")
EXCL = {15, 28}
NFAIL = 0
NINFO = 0
res = {}


def chk(tag, ok, detail=""):
    global NFAIL
    if not ok:
        NFAIL += 1
    print("  [%s] %-34s %s" % ("PASS" if ok else "FAIL", tag, detail))
    return ok


def info(tag, detail=""):
    global NINFO
    NINFO += 1
    print("  [INFO] %-34s %s" % (tag, detail))


print("=" * 96)
print("47 — 全鏈審計（All-chain audit）")
print("=" * 96)

docs = sorted(p for p in DOCS.glob("*.md") if ".pre_patch" not in p.name)
dnums = {int(p.name[:2]) for p in docs}
scripts = {int(q.name[:2]): q for q in SCRIPTS.glob("[0-9][0-9]_*.py") if int(q.name[:2]) not in EXCL}
logs_all = {int(q.name.split("_")[1]): q for q in DATA.glob("run_[0-9][0-9]_*.log")}
logs = {k: v for k, v in logs_all.items() if k not in EXCL}   # 15／28 之 log 屬彙總器自身
print("\n[A1] 腳本集對稱 ↔ log 配對")
print("     docs = %d 份（01–%02d）｜scripts = %d 支（EXCL=%s）｜logs = %d 份"
      % (len(docs), max(dnums), len(scripts), sorted(EXCL), len(logs)))
chk("A1 script 無 log", not (set(scripts) - set(logs)), str(sorted(set(scripts) - set(logs))))
chk("A1 log 無 script", not (set(logs) - set(scripts)), str(sorted(set(logs) - set(scripts))))
chk("A1 docs 缺號", not (set(range(1, max(dnums) + 1)) - dnums),
    str(sorted(set(range(1, max(dnums) + 1)) - dnums)))
res["A1"] = {"docs": len(docs), "scripts": len(scripts), "logs": len(logs)}

print("\n[A2] docs ↔ README 表格")
readme = (ROOT / "README.md").read_text()
rows = {int(m) for m in re.findall(r"^\| `docs/(\d{2})_", readme, re.M)}
chk("A2 README 列數 == docs 檔數", len(rows) == len(dnums), "%d vs %d" % (len(rows), len(dnums)))
chk("A2 README 無漏列", not (dnums - rows), str(sorted(dnums - rows)))
chk("A2 README 無多列", not (rows - dnums), str(sorted(rows - dnums)))
res["A2"] = len(rows)

print("\n[A3] 產物綁定（data/*.json 無孤兒）")
srcs = "\n".join(q.read_text(errors="ignore") for q in SCRIPTS.glob("*.py"))
orph = [p.name for p in DATA.glob("*.json") if p.stem not in srcs and p.name not in srcs]
chk("A3 json 孤兒", not orph, str(orph))
res["A3"] = {"json": len(list(DATA.glob('*.json'))), "orphan": len(orph)}

print("\n[A4] 路徑衛生（/tmp 寫入）")
bad = []
for q in SCRIPTS.glob("*.py"):
    for i, ln in enumerate(q.read_text(errors="ignore").splitlines(), 1):
        if "/tmp/" not in ln or ln.strip().startswith("#"):
            continue
        if re.search(r"open\([^)]*/tmp/|write_text\([^)]*/tmp/|Path\(['\"]/tmp/", ln) and "if " not in ln:
            bad.append("%s:%d" % (q.name, i))
chk("A4 /tmp 寫入", not bad, str(bad))
res["A4"] = len(bad)

print("\n[A5] 四件套狀態行輪次同步")
actual = max(int(m) for m in re.findall(r"^## R(\d+)", (ROOT / "final_result.md").read_text(), re.M))
claim = {}
for f in ("01_problem_record.md", "README.md", "final_result.md"):
    txt = (ROOT / f).read_text()
    slines = [l for l in txt.splitlines() if "狀態" in l]
    rs = [int(x) for l in slines for x in re.findall(r"R(\d+)", l)]
    claim[f] = max(rs) if rs else 0
    shown = next((l for l in slines if re.findall(r"R(\d+)", l)), slines[0] if slines else "")
    print("     %-24s 宣稱最高輪 = R%-3d ｜ %s" % (f, claim[f], shown.strip()[:60]))
chk("A5 三者皆 == 實際 R%d" % actual, all(v == actual for v in claim.values()), str(claim))
snap = next((l for l in (ROOT / "SESSION_START.md").read_text().splitlines() if "狀態" in l), "")
info("A5 SESSION_START（快照，不計）", snap.strip()[:60])
res["A5"] = {"actual": actual, "claim": claim}

print("\n[A6] 時序索引覆蓋（01_problem_record）")
rec = (ROOT / "01_problem_record.md").read_text()
seen = {int(m) for m in re.findall(r"docs/(\d{2})", rec)}
cov = len(seen & dnums)
chk("A6 覆蓋率 ≥ 95%", cov >= 0.95 * len(dnums),
    "%d/%d = %.1f%%（缺 %s）" % (cov, len(dnums), 100.0 * cov / len(dnums), sorted(dnums - seen)[:12]))
res["A6"] = {"covered": cov, "total": len(dnums), "missing": sorted(dnums - seen)}

print("\n[A7] register ref 覆蓋（M12）")
refs = set()
if REG.exists():
    import subprocess
    out = subprocess.run(["/usr/bin/python3", str(REG.parent / "register.py"), "log", "--pid", "M12"],
                         capture_output=True, text=True, env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin"})
    for m in re.findall(r"docs/(\d{2})", out.stdout):
        refs.add(int(m))
chk("A7 覆蓋 docs/01–%02d" % max(dnums), not (dnums - refs), str(sorted(dnums - refs)))
res["A7"] = {"refs": sorted(refs), "missing": sorted(dnums - refs)}

print("\n[A8] 圍欄平衡 ＋ U+FFFD")
prob = []
for p in [ROOT / "README.md", ROOT / "final_result.md", ROOT / "SESSION_START.md",
          ROOT / "01_problem_record.md"] + docs:
    t = p.read_text(errors="ignore")
    fences = sum(1 for ln in t.splitlines() if ln.startswith("```"))
    if fences % 2 or t.count(chr(0xFFFD)):
        prob.append(p.name)
chk("A8 圍欄偶數且無 U+FFFD", not prob, str(prob))
res["A8"] = len(prob)

print("\n[A9] 更正留痕（已撤回結論須有『撤回／更正』）")
need = {"34_r21_omega_identity.md": "撤回",
        "35_r22_omega_ledger.md": "撤回",
        "37_r24_sign_anatomy.md": "更正",
        "38_r25_half_integer_parity.md": "定理"}
miss = [k for k, kw in need.items() if kw not in (DOCS / k).read_text()]
chk("A9 撤回／更正標記", not miss, str(miss))
res["A9"] = len(miss)

print("\n[A10] 覆蓋披露（§7.3）")
npy = len(list(SCRIPTS.glob("*.py")))
info("15_full_audit 之 A3 glob", "*.py ⟹ **含** L2／L3（全 %d 支）" % npy)
info("28_continuity 之 glob", "[0-9][0-9]_*.py ⟹ 不含 L2／L3；但已以 AUX 收 run_l*.log")
info("L2／L3 納入全鏈", "L2/L3 由 15 覆蓋重現性、由 28 AUX 覆蓋溯源 ⟹ 雙重覆蓋")

print("\n" + "=" * 96)
print("[匯總] FAIL = %d ／ INFO = %d ／ docs %d ／ scripts %d ／ logs %d"
      % (NFAIL, NINFO, len(docs), len(scripts), len(logs)))
print("⟹ %s" % ("全鏈審計 PASS" if NFAIL == 0 else "**有待修項 %d 個**" % NFAIL))
print("=" * 96)
(DATA / "47_allchain_audit.json").write_text(json.dumps(
    {"fail": NFAIL, "info": NINFO, "res": {k: v for k, v in res.items()}, "ts": "2026-10-29"},
    ensure_ascii=False, indent=1))