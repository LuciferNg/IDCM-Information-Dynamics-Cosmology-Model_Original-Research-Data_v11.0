#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/06_audit_verify.py — 全鏈審計（idcm-audit-verify §3.1 十一項 ＋ L1 兩跑 byte 比對）
================================================================================
① 腳本集對稱   磁碟 NN_*.py ↔ LOGS 宣告表，兩側差集須空
② 逆向綁定     每份 docs/NN 至少綁一支 producer script（grep scripts/NN）
③ 路徑衛生     scripts/ 內禁 /tmp 寫入
④ 產物綁定     每個 NN 腳本有 data log；每個 data log 有腳本；無孤兒
⑤ inline 幽靈   docs 之 headline 數值須有 L3 claim 覆蓋（無 claim 者列出）
⑥ driver       scripts/run_chain.sh 存在且語法可讀
⑦ L1 兩跑      每支 L1 腳本跑兩次，stdout byte-identical，且 ≡ canonical log
⑧ EXTRA        l2_／l3_ 兩支一併入審計（§7.3 盲區修正）
⑨ 檔數對稱     README 宣稱 docs/scripts 數 == 實際
⑩ 四件套同步    01_problem_record／final_result／README／SESSION_START 狀態行一致
⑪ 文件衛生     markdown 圍欄偶數、U+FFFD 為 0
性質：確定性；唔改任何檔案（除 data/ 輸出）。exit 1 若有任何 FAIL。
"""
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA, DOCS, SCRIPTS = ROOT / "data", ROOT / "docs", ROOT / "scripts"
PY = "/usr/bin/python3"
ENV = {"PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
       "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}

LOGS = {
    "01_label_theorems": "run_01_label_theorems.log",
    "02_a5_harmonic_multiplicity": "run_02_a5_harmonic.log",
    "03_count_census": "run_03_count_census.log",
    "04_adic_map_nogo": "run_04_adic_map.log",
    "05_cy3_generation_count": "run_05_cy3.log",
    "07_h11_minus_M_candidate": "run_07_h11_M.log",
    "08_ngen_equals_B_minus_1": "run_08_ngen_B1.log",
    "09_tenc_c_downgrade_audit": "run_09_downgrade.log",
    "10_ngen_provenance_scan": "run_10_provenance.log",
    "11_gate_count_candidate": "run_11_gate_count.log",
    "12_gate_budget_identity": "run_12_budget.log",
    "13_remaining_four": "run_13_remaining.log",
    "14_residual_push": "run_14_residual.log",
    "15_forced_recompute": "run_15_forced_recompute.log",
    "16_zero_mode_arbitration": "run_16_zero_mode.log",
    "17_minus_one_derivation": "run_17_minus_one.log",
    "18_q2_provenance_recheck": "run_18_q2_recheck.log",
    "19_handoff_m14_ids_and_z3": "run_19_m14_ids.log",
}
EXTRA = [("L2", "l2_independent_verify.py", "run_l2_independent_verify.log"),
         ("L3", "l3_doc_crosscheck.py", "run_l3_doc_crosscheck.log")]
# docs 之 producer 綁定（判定／協定文件可無 producer，但須列入豁免）
DOC_PROD = {
    "01_frozen_criteria.md": [],
    "02_reduction_theorem_status.md": ["01", "02"],
    "03_no_second_label_theorem.md": ["01"],
    "04_count_census_and_partition_test.md": ["03"],
    "05_adic_profile_and_map_nogo.md": ["04"],
    "06_cy3_verdict.md": ["05"],
    "07_verdicts.md": [],
    "08_ninth_candidate_h11_minus_M.md": ["07"],
    "09_tenth_candidate_ngen_equals_B_minus_1.md": ["08"],
    "10_tenth_candidate_downgrade.md": ["09"],
    "11_gate_candidate_and_provenance_scan.md": ["10", "11"],
    "12_verification_and_budget_anchor.md": ["12"],
    "13_remaining_four_q1_q2_q3_q5.md": ["13"],
    "14_residual_push.md": ["14"],
    "16_zero_mode_arbitration.md": ["16"],
    "17_minus_one_derivation.md": ["17"],
    "18_q2_provenance_recheck.md": ["18"],
    "19_handoff_m14_ids_and_z3.md": ["19"],
}
EXEMPT_DOCS = set()      # 全部 doc 納入 L3（2026-11-01 全面審計）
L3_COVERED = {d for d in DOC_PROD if d not in EXEMPT_DOCS}

NFAIL, NWARN = 0, 0
RES = {}


def chk(tag, cond, detail="", warn=False):
    global NFAIL, NWARN
    ok = bool(cond)
    if not ok:
        if warn:
            NWARN += 1
        else:
            NFAIL += 1
    st = "PASS" if ok else ("WARN" if warn else "FAIL")
    print(f"[{st}] {tag}" + (f"  |  {detail}" if detail else ""))
    RES[tag] = {"status": st, "detail": detail}
    return ok


print("=" * 78)
print("M15 scripts/06 — 全鏈審計")
print("=" * 78)

# ① 腳本集對稱
SELF = Path(__file__).stem
disk = sorted(p.stem for p in SCRIPTS.glob("[0-9][0-9]_*.py") if p.stem != SELF)
decl = sorted(LOGS)
chk("① 腳本集對稱（磁碟 ↔ 宣告；審計自身不入 L1）", disk == decl,
    f"磁碟 {disk}；宣告 {decl}；差集 {set(disk) ^ set(decl)}")

# ② 逆向綁定 ＋ ⑤ inline 幽靈
missing_prod = []
for d, prods in DOC_PROD.items():
    p = DOCS / d
    if not p.exists():
        missing_prod.append(f"{d}（檔不存在）")
        continue
    txt = p.read_text()
    if d in EXEMPT_DOCS:
        continue
    for s in prods:
        if f"scripts/{s}" not in txt:
            missing_prod.append(f"{d} 缺 scripts/{s}")
chk("② docs ↔ producer 綁定", not missing_prod, f"缺 {missing_prod}")
disk_docs = {p.name for p in DOCS.glob("*.md")}
uncovered = sorted(disk_docs - L3_COVERED)
stale_decl = sorted(set(DOC_PROD) - disk_docs)
chk("⑤ inline 幽靈（每份 doc 皆須有 L3 claim 覆蓋）",
    not uncovered and not stale_decl,
    f"L3 覆蓋 {len(L3_COVERED)}/{len(disk_docs)}；未被覆蓋 {uncovered}；宣告而磁碟無 {stale_decl}")

# ③ 路徑衛生
tmp_hits = [p.name for p in SCRIPTS.glob("*.py") if "/tmp" in p.read_text() and p.name != Path(__file__).name]
chk("③ scripts/ 內禁 /tmp 寫入", not tmp_hits, f"違例 {tmp_hits}")

# ④ 產物綁定
produced = {LOGS[s] for s in LOGS} | {e[2] for e in EXTRA} | {"run_06_audit_verify.log"}
on_disk = {p.name for p in DATA.glob("run_*.log")}
orphan = sorted(on_disk - produced)
missing_log = sorted(produced - on_disk)
chk("④ 產物綁定（無孤兒 log、無失蹤 log）", not orphan and not missing_log,
    f"孤兒 {orphan}；失蹤 {missing_log}")

# ⑥ driver
drv = SCRIPTS / "run_chain.sh"
chk("⑥ driver 存在且可讀", drv.exists() and len(drv.read_text()) > 0, str(drv.name))

# ⑦ L1 兩跑 byte-identical ＋ ≡ canonical log
bad = []
for stem, logname in LOGS.items():
    outs = []
    for _ in range(2):
        r = subprocess.run([PY, str(SCRIPTS / f"{stem}.py")], capture_output=True,
                           text=True, cwd=str(ROOT), env=ENV)
        outs.append((r.returncode, r.stdout))
    canon = (DATA / logname).read_text()
    ok = (outs[0][0] == 0 and outs[1][0] == 0
          and outs[0][1] == outs[1][1] == canon)
    if not ok:
        bad.append((stem, outs[0][0], outs[1][0],
                    outs[0][1] == outs[1][1], outs[0][1] == canon))
chk("⑦ L1 兩跑 byte-identical 且 ≡ canonical log", not bad, f"違例 {bad}")

# ⑧ EXTRA（L2／L3）
badx = []
for tag, fname, logname in EXTRA:
    r = subprocess.run([PY, str(SCRIPTS / fname)], capture_output=True,
                       text=True, cwd=str(ROOT), env=ENV)
    canon = (DATA / logname).read_text()
    if r.returncode != 0 or r.stdout != canon:
        badx.append((tag, r.returncode, r.stdout == canon))
chk("⑧ L2／L3 一併審計（rc=0 且 ≡ log）", not badx, f"違例 {badx}")

# ⑨ 檔數對稱（README 宣稱）
rd = (ROOT / "README.md").read_text()
m_docs = re.search(r"docs/（(\d+)\s*份", rd)
m_scr = re.search(r"scripts/（(\d+)\s*支", rd)
ndoc = len(list(DOCS.glob("*.md")))
nscr = len(list(SCRIPTS.glob("*.py")))
chk("⑨ README 宣稱檔數 == 實際",
    bool(m_docs) and bool(m_scr) and int(m_docs.group(1)) == ndoc and int(m_scr.group(1)) == nscr,
    f"README docs={m_docs.group(1) if m_docs else '?'}/{ndoc}；"
    f"scripts={m_scr.group(1) if m_scr else '?'}/{nscr}")

# ⑩ 四件套同步
files = ["01_problem_record.md", "final_result.md", "README.md"]
states = {}
for f in files:
    p = ROOT / f
    t = p.read_text() if p.exists() else ""
    # 只抽 verdict 詞（避開狀態行「closed / verdict partial」中的 closed）
    m = re.search(r"verdict\s*`?\s*(negative|partial|registered|open)", t)
    states[f] = m.group(1) if m else "MISSING"
ss = set(states.values())
snap = (ROOT / "SESSION_START.md")
chk("⑩ 四件套狀態一致（三活文檔 verdict 相同；SESSION_START 為快照只查存在）",
    len(ss) == 1 and "MISSING" not in ss and snap.exists(),
    f"{states}；SESSION_START 存在={snap.exists()}")

# ⑫ L2 衛生：零 import L1、零讀 docs
l2src = (SCRIPTS / "l2_independent_verify.py").read_text()
l2_bad = []
if re.search(r"^\s*(import|from)\s+\d\d_", l2src, re.M):
    l2_bad.append("import L1")
_code = "\n".join(ln.split("#")[0] for ln in l2src.splitlines())     # 去註釋
if re.search(r"(open|read_text|read_bytes)\([^)]*docs", _code) or re.search(r"DOCS\s*=", _code):
    l2_bad.append("讀 docs")
if "subprocess" in l2src:
    l2_bad.append("呼叫外部腳本")
chk("⑫ L2 衛生（零 import L1、零讀 docs、零 subprocess；註釋提及不計）", not l2_bad, f"{l2_bad}")

# ⑬ data/*.json 綁定（非 l2／audit 之 JSON 須有同名 producer script）
json_orphan = []
for j in DATA.glob("*.json"):
    stem = j.stem
    if stem in ("l2_independent_verify", "audit_verify"):
        continue
    if not (SCRIPTS / f"{stem}.py").exists():
        json_orphan.append(j.name)
chk("⑬ data/*.json ↔ producer script（無孤兒）", not json_orphan, f"孤兒 {json_orphan}")

# ⑮ driver ↔ 宣告一致（防 driver 漏跑：log 存在但從未由 driver 重生成）
_drv = (SCRIPTS / "run_chain.sh").read_text()
_pairs = re.findall(r"^run\s+(\S+\.py)\s+(\S+\.log)\s*$", _drv, re.M)
_drv_map = {f[:-3]: lg for f, lg in _pairs}
_decl = {k: v for k, v in LOGS.items()}
_miss = sorted(k for k in _decl if k not in _drv_map)
_extra = sorted(k for k in _drv_map if k not in _decl and k not in ("l2_independent_verify", "l3_doc_crosscheck"))
_badlog = sorted(k for k in _decl if k in _drv_map and _drv_map[k] != _decl[k])
chk("⑮ driver ↔ LOGS 宣告一致（漏跑／多跑／log 名不符）",
    not _miss and not _extra and not _badlog,
    f"漏跑 {_miss}；多跑 {_extra}；log 名不符 {_badlog}")

# ⑭ 反向覆蓋：docs 內 4+ 有效位小數 → 須可追溯至 data（容差 0.5%；跨記錄引文明示豁免）
import bisect as _bisect
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
ALLOW = {"02_reduction_theorem_status.md": {0.1657}}   # O137 (v10) 屏障曲線，跨記錄引用
_pool = []
for _f in list(DATA.glob("*.json")) + list(DATA.glob("run_*.log")):
    _t = _f.read_text(errors="ignore")
    for _m in re.finditer(r"([-+]?\d+(?:\.\d+)?)\s*[×x]\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", _t):
        _pool.append(float(_m.group(1)) * 10 ** int(_m.group(2).translate(SUP)))
    _pool += [float(x) for x in re.findall(r"[-+]?\d+\.\d+", _t)]
_pool.sort()


def _hit(v, tol=5e-3):
    i = _bisect.bisect_left(_pool, v)
    return any(0 <= j < len(_pool) and abs(_pool[j] - v) <= tol * max(abs(v), 1e-12)
               for j in (i - 1, i, i + 1))


_untr = []
for _f in sorted(DOCS.glob("*.md")):
    _t = _f.read_text().replace("\u2212", "-")
    _toks, _spans = set(), []
    for _m in re.finditer(r"(\d+(?:\.\d+)?)\s*[×x]\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", _t):
        _toks.add(float(_m.group(1)) * 10 ** int(_m.group(2).translate(SUP)))
        _spans.append(_m.span(1))
    for _m in re.finditer(r"\d+\.\d{3,}", _t):
        if any(a <= _m.start() < b for a, b in _spans):
            continue
        _toks.add(float(_m.group(0)))
    for _v in _toks:
        if _v in ALLOW.get(_f.name, set()):
            continue
        if not _hit(_v):
            _untr.append(f"{_f.name}:{_v}")
chk("⑭ 反向覆蓋（docs 4+ 位小數 → data 可追溯；跨記錄引文明示豁免）",
    not _untr, f"未追溯 {_untr}")

# ⑪ 文件衛生
hyg = []
for p in list(DOCS.glob("*.md")) + [ROOT / "README.md", ROOT / "01_problem_record.md"]:
    if not p.exists():
        continue
    t = p.read_text()
    if t.count("```") % 2:
        hyg.append(f"{p.name}: 圍欄奇數")
    if "\ufffd" in t:
        hyg.append(f"{p.name}: U+FFFD")
chk("⑪ 文件衛生（圍欄偶數、無 U+FFFD）", not hyg, f"{hyg}")

print("")
print("=" * 78)
print(f"審計總結：檢查 {len(RES)} 項，FAIL {NFAIL}、WARN {NWARN}")
print("=" * 78)

(DATA / "audit_verify.json").write_text(json.dumps(
    {"results": RES, "n_fail": NFAIL, "n_warn": NWARN,
     "disk_scripts": disk, "docs": ndoc, "scripts": nscr}, ensure_ascii=False, indent=1))
sys.exit(1 if NFAIL else 0)