"""M5 docs/11 配套 — 全量審計：11 支腳本 × 10 份文檔（重跑 byte-identical／docstring／數字回溯／完整性六掃）
=================================================================================================================
審計範圍（2026-10-24，調律者指示「先驗證和審計 M4 及 M5 所有之前未經驗證的文檔和腳本」）：
  - `scripts/01..11_*.py`（11 支；`scripts/deprecated/` 之被否證版本按設計排除，另列 INFO）
  - `docs/01..10_*.md`（10 份；`11` 同 `10` 共用 docs/10）
  - `data/run_NN_*.log`（11 支正典 log ＋ 1 支 v1 buggy log 存證）
  - `data/NN_*.json`（11 個）＋ `data/pdg_mass_width_2024.mcd`（外部輸入，PDG 官方下載）

四類主檢（T1–T4，機械化）＋ 六項完整性掃描（I1–I6）：
  v2（2026-10-24 增補）：PAIR_OVERRIDE 加 "13"→"12"，覆蓋 A3 腳本↔文檔；標題改動態計數。
  T1 重跑一致：臨時目錄（預置 PDG 輸入）重跑 ⟹ stdout 與正典 log byte 相同、JSON 語義相同
  T2 docstring 對齊：docs/NN「## 一、前設問題」fenced block 逐字等於腳本 docstring
  T3 數字回溯：docs 內 ≥3 位有效數字須能回溯（字串／0.5% 相對／三桶分流＋運算式復算 0.6%）
  T4 結構：docs H1 前綴、四個必要節、README 宣稱檔數 = 實際
  I1 幽靈路徑（舊 v10 p16 路徑／舊 repo root）／I2 `/tmp` 寫入／I3 空洞 check（`check(..., True)`）
  I4 data JSON 孤兒／I5 腳本↔文檔覆蓋／I6 交付檔案齊全（README/final_result/01_problem_record/SESSION_START）
輸出：`data/12_full_audit.json` ＋ stdout 表
"""
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS, SCRIPTS, DATA = ROOT / "docs", ROOT / "scripts", ROOT / "data"
PY = "/usr/bin/python3"
NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
PAIR_OVERRIDE = {"11": "10", "13": "12", "14": "13", "15": "14", "16": "15", "17": "16", "18": "17", "19": "18", "20": "19", "21": "19", "22": "20", "23": "21", "24": "22", "25": "23", "26": "24", "27": "25", "28": "26", "29": "27"}
#   11→10（同 doc 兩腳本）；13→12；14→13；15→14；16→15（docs 比 scripts 少 1）


def nums(text):
    out = []
    for tok in NUM.findall(text):
        try:
            v = float(tok)
        except ValueError:
            continue
        if len(re.sub(r"[-.eE+]", "", tok).lstrip("0")) >= 3 and abs(v) > 0:
            out.append((tok, v))
    return out


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def doc_block(txt, header):
    """回傳該節內**全部** fenced block（docs/10 收錄 10、11 兩支腳本之 docstring）"""
    i = txt.find(header)
    if i < 0:
        return None
    j = txt.find("\n## ", i + 1)
    seg = txt[i: j if j > 0 else len(txt)]
    return [m.group(1) for m in re.finditer(r"```[a-z]*\n(.*?)```", seg, re.S)]


report = {"meta": {"script": "12_full_audit", "date": "2026-10-24", "python": PY, "root": str(ROOT)},
          "pairs": {}, "counts": {}, "violations": [], "info": []}
print("=" * 96)
print(f"M5 docs/11 — 全量審計（{len(list(SCRIPTS.glob('[0-9][0-9]_*.py')))-1} 支腳本（不含本審計器） × {len(list(DOCS.glob('[0-9][0-9]_*.md')))-1} 份推導文檔（不含 docs/11 本身））")
print("=" * 96)

SELF = Path(__file__).name
pairs = []
for sc in sorted(SCRIPTS.glob("[0-9][0-9]_*.py")):
    if sc.name == SELF:
        continue
    nn = sc.name[:2]
    target = PAIR_OVERRIDE.get(nn, nn)
    doc = next(iter(sorted(DOCS.glob(f"{target}_*.md"))), None)
    pairs.append((nn, sc, doc, DATA / f"run_{sc.stem}.log"))

dep = sorted((SCRIPTS / "deprecated").glob("*.py")) if (SCRIPTS / "deprecated").is_dir() else []
print(f"\n配對 {len(pairs)} 組；deprecated（按設計排除，只列 INFO）={[p.name for p in dep]}")
print(f"docs={len(list(DOCS.glob('[0-9][0-9]_*.md')))} scripts={len(list(SCRIPTS.glob('[0-9][0-9]_*.py')))} "
      f"logs={len(list(DATA.glob('run_*.log')))} json={len(list(DATA.glob('*.json')))}")

tmp = Path(tempfile.mkdtemp(prefix="m5audit_", dir=os.environ.get("BH_AGENT_WORKSPACE", "/tmp")))
print(f"臨時目錄：{tmp}")

for nn, sc, doc, log in pairs:
    rec = {"script": sc.name, "doc": doc.name if doc else None, "log": log.name, "log_exists": log.exists()}
    work = tmp / nn
    (work / "data").mkdir(parents=True, exist_ok=True)
    shutil.copy(sc, work / sc.name)
    mcd = DATA / "pdg_mass_width_2024.mcd"
    if mcd.exists():
        shutil.copy(mcd, work / "data" / mcd.name)
    # 全項目驗證器（scripts/16）需要 root docs、docs/、data/*.json ⟹ 一併複製（隔離性保持：
    # 只讀複本，寫入仍落 work 內）
    for f in ROOT.glob("*.md"):
        shutil.copy(f, work / f.name)
    (work / "docs").mkdir(exist_ok=True)
    for f in DOCS.glob("*.md"):
        shutil.copy(f, work / "docs" / f.name)
    for f in DATA.glob("*.json"):
        shutil.copy(f, work / "data" / f.name)
    for f in DATA.glob("run_*.log"):
        shutil.copy(f, work / "data" / f.name)
    (work / "scripts").mkdir(exist_ok=True)
    for f in SCRIPTS.glob("*.py"):
        shutil.copy(f, work / "scripts" / f.name)
    try:
        r = subprocess.run([PY, sc.name], cwd=work, capture_output=True, text=True, timeout=1800)
        out, err, rc = r.stdout, r.stderr, r.returncode
    except subprocess.TimeoutExpired:
        out, err, rc = "", "TIMEOUT>1800s", 124
    rec["exit"] = rc
    rec["stderr_tail"] = err.strip().splitlines()[-2:] if err.strip() else []
    rec["rerun_identical"] = None
    if log.exists():
        canon = log.read_text()
        rec["rerun_identical"] = (canon == out)
        if canon != out:
            cl, ol = canon.splitlines(), out.splitlines()
            rec["diff"] = {"canon_lines": len(cl), "rerun_lines": len(ol),
                           "first_diff": next((i for i, (a, b) in enumerate(zip(cl, ol)) if a != b), None)}
    rec["json"] = []
    for jf in sorted((work / "data").glob("*.json")):
        cj = DATA / jf.name
        if not cj.exists():
            rec["json"].append((jf.name, "canonical-missing"))
            continue
        try:
            same = json.loads(cj.read_text()) == json.loads(jf.read_text())
        except Exception as e:
            same = f"parse-error {e}"
        rec["json"].append((jf.name, same))
    src = sc.read_text()
    m = re.search(r'"""(.*?)"""', src, re.S)
    docstr = m.group(1) if m else None
    rec["docstring_match"] = None
    if doc and docstr:
        blks = doc_block(doc.read_text(), "## 一、前設問題")
        if blks:
            rec["docstring_match"] = any(norm(b) == norm(docstr) for b in blks)
        elif blks is not None:
            rec["docstring_match"] = False
    if doc:
        dtxt = doc.read_text()
        dtxt = re.sub(r"(?<=\d) (?=\d{3}(?!\d))", "", dtxt)      # 「39 304」→「39304」（千位空隔）
        LOGS = {"11": ["10_tree_two_body", "11_scale_fixed_test"]}   # docs/10 由 10、11 兩腳本共用
        stems = LOGS.get(nn, [sc.stem])
        ltxt = "".join((DATA / f"run_{s}.log").read_text() for s in stems if (DATA / f"run_{s}.log").exists())
        jsontxt = "".join((DATA / f"{s}.json").read_text() for s in stems if (DATA / f"{s}.json").exists())
        corpus = ltxt + jsontxt + sc.read_text()                     # 程式常數亦屬可查語料
        log_nums = [v for _, v in nums(ltxt + jsontxt)][:400]
        HC = 197.3269804                                             # ħc [MeV·fm]
        buckets = {"meta": 0, "external": 0, "derived": 0, "pending": 0}
        unexplained = []

        def near(v, c):
            return abs(v - c) <= max(0.006 * abs(v), 1e-12)

        for line in dtxt.splitlines():
            line_nums = [v for _, v in nums(line)][:12]
            for tok in NUM.findall(line):
                try:
                    v = float(tok)
                except ValueError:
                    continue
                if len(re.sub(r"[-.eE+]", "", tok).lstrip("0")) < 3 or v == 0:
                    continue
                if tok in corpus:
                    continue
                if re.fullmatch(r"20\d\d", tok) or re.fullmatch(r"[0-9a-f]{6,}", tok):
                    buckets["meta"] += 1
                    continue
                cands = []
                for a in line_nums:
                    cands += [a, 100 * a, a / 100, a * HC, a / HC, 1 / a if a else 0]
                    cands += [a ** 2, a ** 3]
                for a in line_nums:
                    for b in log_nums:
                        if b:
                            cands += [a + b, a - b, a * b, a / b, b / a]
                for a in log_nums:
                    cands += [a, 100 * a, a / 100, a * HC, a / HC, 1 / a if a else 0]
                if any(near(v, c) for c in cands if c):
                    buckets["derived"] += 1
                    continue
                if any(h in line for h in ("PDG", "FLAG", "CODATA", "Planck", "DESI", "lattice", "docs/", "O1", "O2",
                                           "O8", "O9", "O4", "register", "sha", "MeV @", "@ 2 GeV", "首版", "v1",
                                           "deprecated", "外部")):
                    buckets["external"] += 1
                    continue
                for a in log_nums:
                    if a and near(v, a):
                        buckets["pending"] += 1
                        break
                else:
                    unexplained.append((tok, line.strip()[:70]))
        rec["t3"] = {**buckets, "unexplained": [u[0] for u in unexplained[:12]],
                     "unexplained_n": len(unexplained), "ctx": [u[1] for u in unexplained[:4]]}
    if doc:
        dt = doc.read_text()
        rec["sections_ok"] = ("前設問題" in dt) and ("結論" in dt)
        rec["has_criteria_section"] = "事前判準" in dt
        rec["h1_ok"] = dt.splitlines()[0].startswith("# M5 docs/")
    report["pairs"][nn] = rec
    t3 = rec.get("t3", {})
    print(f"  [{nn}] exit={rec['exit']} 重跑={rec['rerun_identical']} docstring={rec['docstring_match']} "
          f"T3 {t3.get('meta')}/{t3.get('external')}/{t3.get('pending')} 未解釋={t3.get('unexplained_n')} "
          f"結構={rec.get('sections_ok')} json={[c[1] for c in rec['json']] or '—'}")

# ---------- I1–I6 完整性掃描 ----------
print("\n[完整性掃描]")
scan = {}
srcs = {p.name: p.read_text() for p in SCRIPTS.glob("*.py")}
mdtext = {p.name: p.read_text() for p in list(DOCS.glob("*.md")) + [ROOT / "README.md", ROOT / "01_problem_record.md", ROOT / "final_result.md", ROOT / "SESSION_START.md"] if p.exists()}
I1 = [(n, pat) for n, t in srcs.items() if n != SELF for pat in ("/tmp/", "/var/tmp") if f"open(" in t and pat in t]
scan["I1_tmp_write"] = I1
I2 = [(n, pat) for n, t in srcs.items() if n != SELF for pat in ("IDCM-v10.0--It-From-Beat", "p16_hadron_bound_state_spectrum") if pat in t]
scan["I2_ghost_path"] = I2
I3 = []
for n, t in srcs.items():
    for m in re.finditer(r"check\(\s*[\"'][^\"']*[\"']\s*,\s*True\b", t):
        I3.append((n, m.group(0)[:60]))
scan["I3_vacuous_check"] = I3
I4 = []
for jf in sorted(DATA.glob("*.json")):
    producer = [n for n, t in srcs.items() if jf.name in t]
    if not producer:
        I4.append(jf.name)
scan["I4_orphan_json"] = I4
I5 = []
for sc in sorted(SCRIPTS.glob("[0-9][0-9]_*.py")):
    if sc.name == SELF:
        continue
    stem = sc.stem
    if not any(stem in t for t in mdtext.values()):
        I5.append(sc.name)
scan["I5_nondocumented_script"] = I5
need = ["README.md", "final_result.md", "01_problem_record.md", "SESSION_START.md"]
scan["I6_missing_root_files"] = [f for f in need if not (ROOT / f).exists()]
for k, v in scan.items():
    print(f"  {k}: {len(v)} {v[:6] if v else ''}")
report["integrity"] = {k: (v if len(v) < 50 else v[:50]) for k, v in scan.items()}

n_doc = len(list(DOCS.glob("[0-9][0-9]_*.md")))
n_sc = len(list(SCRIPTS.glob("[0-9][0-9]_*.py")))
rd = (ROOT / "README.md").read_text()
claims = {"docs": re.findall(r"docs/` \| 推導記錄（(\d+) 檔", rd),
          "scripts": re.findall(r"計算腳本（核心 (\d+) 檔", rd),
          "log_json": re.findall(r"log ＋ JSON（(\d+)＋(\d+)）", rd)}
report["counts"] = {"docs": n_doc, "scripts": n_sc,
                    "logs": len(list(DATA.glob("run_*.log"))), "json": len(list(DATA.glob("*.json"))),
                    "claims": claims}
print(f"\n[計數] actual docs={n_doc} scripts={n_sc} logs={len(list(DATA.glob('run_*.log')))} "
      f"json={len(list(DATA.glob('*.json')))}；README 宣稱 {claims}")

viol = []
for nn, r in report["pairs"].items():
    if r["exit"] != 0:
        viol.append(f"{nn}: exit {r['exit']}")
    if r["rerun_identical"] is False:
        viol.append(f"{nn}: 重跑非 byte-identical (diff={r.get('diff')})")
    if r["docstring_match"] is False:
        viol.append(f"{nn}: docstring ≠ docs §一")
    if r.get("sections_ok") is False:
        viol.append(f"{nn}: 缺必要節")
    if r.get("h1_ok") is False:
        viol.append(f"{nn}: doc H1 非「# M5 docs/」")
    for jn, ok in r["json"]:
        if ok is not True:
            viol.append(f"{nn}: JSON {jn} 不一致（{ok}）")
    if r.get("t3", {}).get("unexplained_n", 0) > 0:
        viol.append(f"{nn}: T3 未解釋數字 {r['t3']['unexplained']}")
for k in ("I2_ghost_path", "I3_vacuous_check", "I4_orphan_json", "I6_missing_root_files"):
    if scan[k]:
        viol.append(f"{k}: {scan[k]}")
report["violations"] = viol
print(f"\n[違規] {len(viol)}")
for v in viol[:40]:
    print("  -", v)
with open(DATA / "12_full_audit.json", "w") as f:
    json.dump(report, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/12_full_audit.json")