#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/18_audit_verify.py — 全鏈審計（idcm-audit-verify §3.1 六項 ＋ L2/L3 納入）

A1 腳本集對稱      ：scripts/NN_*.py ↔ data/run_NN*.log（兩側差集須空；孤兒 log = 0）
A2 逆向綁定        ：每個 headline 數字 → 磁碟上 producer log 真的含該值
A3 路徑衛生        ：L1 腳本零 /tmp 寫入
A4 重跑 byte-identical：全部 L1 重跑，stdout 必須等於 canonical log（trap §7.2）
A5 產物綁定        ：L2／L3 JSON 存在且由各自腳本產出（EXPLICIT EXTRA，trap §7.3）
A6 文檔核          ：引用 data/l3_doc_crosscheck.json 之逐檔命中率
A7 四件套          ：01_problem_record／SESSION_START／final_result（階段性收尾已建）＋ README 狀態行
"""
import json
import re
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
SCRIPTS = ROOT / "scripts"

LOG, RES, NFAIL, NWARN = [], {}, 0, 0


def chk(tag, cond, detail="", warn=False):
    global NFAIL, NWARN
    ok = bool(cond)
    if not ok:
        if warn:
            NWARN += 1
        else:
            NFAIL += 1
    kind = "PASS" if ok else ("WARN" if warn else "FAIL")
    line = f"[{kind}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "warn": (not ok and warn), "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


L1 = sorted(p for p in SCRIPTS.glob("[0-2][0-9]_*.py"))
EXTRA = [SCRIPTS / "l2_independent_verify.py", SCRIPTS / "l3_doc_crosscheck.py"]

# 正典映射（2026-10-30 審計建立）：腳本 stem → log stem
#   歷史上部分 log 用短名（run_02_criteria.log），docs 已引用該名 ⟹ 保留原名，改為顯式映射。
LOGMAP = {
    "01_reproduction_gate": "run_01_reproduction_gate",
    "02_criteria_T1_T5": "run_02_criteria",
    "03_null_capacity": "run_03_null_capacity",
    "04_closure_scan": "run_04_closure_scan",
    "05_carry_by_2_candidate": "run_05_carry_by_2",
    "06_universality": "run_06_universality",
    "07_ktree_alpha_probe": "run_07_ktree",
    "08_alpha_source": "run_08_alpha_source",
    "09_precision_budget": "run_09_precision",
    "10_family_x_scan": "run_10_family_x",
    "11_sector_extension": "run_11_sector_extension",
    "12_t1_label_theorem": "run_12_t1_label",
    "13_rival_lattice_audit": "run_13_rival",
    "14_symmetric_trace_emergence": "run_14_symmetric",
    "15_trace_statistics": "run_15_trace_stats",
    "16_trace_target_precision": "run_16_trace_target",
    "17_hypothesis_verification": "run_17_verification",
    "18_audit_verify": "run_18_audit_verify",
    "19_handoff_replies": "run_19_handoff_replies",
    "20_channel_optimization": "run_20_channel",
    "21_r1_r3_push": "run_21_r1_r3",
}
CANON = {p: DATA / f"{LOGMAP[p.stem]}.log" for p in L1 if p.stem in LOGMAP}
# A4 排除本檔自身（否則 18→18 無限遞歸）；本檔之可重現性由「跑兩次 cmp」外部驗證
L1_RUN = [p for p in L1 if p.stem != "18_audit_verify"]


def _norm(s):
    """千分位逗號歸一，令 34,752,256 與 34752256 可比"""
    return re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", s)

# A2：headline 數字 → producer（值字串必須出現在該 producer 的 canonical log）
BIND = [
    ("rho_d 0.403503649635", SCRIPTS / "01_reproduction_gate.py", "0.403503649635"),
    ("sigma_rho_d 0.044785178", SCRIPTS / "09_precision_budget.py", "0.044785178"),
    ("delta_s 3.4250", SCRIPTS / "01_reproduction_gate.py", "3.425"),
    ("trace 0.37761528", SCRIPTS / "14_symmetric_trace_emergence.py", "0.377615281752"),
    ("sigma_trace 0.0116811", SCRIPTS / "15_trace_statistics.py", "0.011681"),
    ("corr −0.99635", SCRIPTS / "15_trace_statistics.py", "-0.9963504"),
    ("z_trace 0.22389", SCRIPTS / "17_hypothesis_verification.py", "0.2239"),
    ("chi2_x2 0.999", SCRIPTS / "10_family_x_scan.py", "0.999"),
    ("chi2_x4 18.207", SCRIPTS / "10_family_x_scan.py", "18.207"),
    ("z_r11_3 0.1248", SCRIPTS / "13_rival_lattice_audit.py", "0.1247"),
    ("N_trials 216", SCRIPTS / "08_alpha_source.py", "216"),
    ("K=5 類大小 34752256", SCRIPTS / "03_null_capacity.py", "34752256"),
]


def main():
    hdr("§A1 腳本集對稱（磁碟 ↔ log）")
    on_disk = {p.stem for p in L1}
    declared = {LOGMAP[s] for s in on_disk if s in LOGMAP}
    actual = {p.stem for p in DATA.glob("run_*.log")}
    actual_l1 = {s for s in actual if re.match(r"^run_[0-2][0-9]_", s)}
    missing_log = sorted(s for s in declared if s not in actual_l1)
    orphan_log = sorted(actual_l1 - declared)
    undeclared = sorted(s for s in on_disk if s not in LOGMAP)
    chk("A1a 每個 L1 腳本都有 canonical log（顯式 LOGMAP）", not missing_log, f"缺 {missing_log}")
    chk("A1b 無孤兒 log", not orphan_log, f"孤兒 {orphan_log}")
    chk("A1c 所有腳本已入 LOGMAP", not undeclared, f"未宣告 {undeclared}")
    chk(f"A1d 腳本數 = log 數（{len(on_disk)}）", len(declared) == len(actual_l1),
        f"{len(declared)} vs {len(actual_l1)}")

    hdr("§A2 逆向綁定（headline → producer log）")
    for label, prod, needle in BIND:
        lg = CANON.get(prod) or (DATA / f"run_{prod.stem}.log")
        if not lg.exists():
            chk(f"A2 {label}", False, f"缺 log {lg.name}")
            continue
        hit = needle in _norm(lg.read_text(encoding="utf-8", errors="replace"))
        chk(f"A2 {label} ← {prod.name}", hit, "" if hit else f"log 內無 {needle!r}")

    hdr("§A3 路徑衛生（/tmp 寫入）")
    bad = []
    for p in L1 + EXTRA:
        t = p.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"(open|write_text|write)\([^)\n]*[\"']/tmp/", t):
            bad.append((p.name, m.group(0)[:40]))
    chk("A3 L1/L2/L3 零 /tmp 寫入", not bad, f"{bad}")

    hdr("§A4 重跑 byte-identical（canonical log 比對）")
    env = {"PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    import os
    e = dict(os.environ); e.update(env)
    drift = []
    for p in L1_RUN:
        r = subprocess.run([sys.executable, str(p)], capture_output=True, env=e, timeout=280)
        canon = CANON[p].read_bytes()
        if r.returncode != 0:
            drift.append((p.name, f"rc={r.returncode}"))
        elif r.stdout != canon:
            drift.append((p.name, "stdout ≠ canonical log"))
    chk(f"A4 全部 L1（{len(L1_RUN)} 支，排除本檔）重跑 byte-identical",
        not drift, f"{drift}")

    hdr("§A5 產物綁定（L2／L3 納入審計；trap §7.3）")
    r2 = subprocess.run([sys.executable, str(EXTRA[0])], capture_output=True, env=e, timeout=280)
    r3 = subprocess.run([sys.executable, str(EXTRA[1])], capture_output=True, env=e, timeout=280)
    chk("A5a L2 重跑 rc=0", r2.returncode == 0, f"rc={r2.returncode}")
    chk("A5b L3 重跑 rc=0", r3.returncode == 0, f"rc={r3.returncode}")
    for f in ("l2_independent_verify.json", "l3_doc_crosscheck.json"):
        chk(f"A5c data/{f} 存在", (DATA / f).exists())

    hdr("§A6 文檔核（引用 L3 結果）")
    p3 = DATA / "l3_doc_crosscheck.json"
    if p3.exists():
        j3 = json.loads(p3.read_text(encoding="utf-8"))
        bd = j3.get("_summary", {}).get("by_doc", {})
        allok = all(v.split("/")[0] == v.split("/")[1] for v in bd.values())
        chk(f"A6 全部 {len(bd)} 份 doc 命中率 100%", allok, f"{bd}")
        chk("A6b L3 FAIL = 0", j3.get("_summary", {}).get("fail") == 0)
    else:
        chk("A6 L3 JSON 存在", False)

    hdr("§A7 四件套")
    pr = ROOT / "01_problem_record.md"
    ss = ROOT / "SESSION_START.md"
    chk("A7a 01_problem_record.md 存在", pr.exists())
    chk("A7b SESSION_START.md 存在", ss.exists())
    fr = ROOT / "final_result.md"
    chk("A7c final_result.md 存在（階段性收尾已建）", fr.exists())
    if fr.exists():
        ftxt = fr.read_text(encoding="utf-8")
        chk("A7c2 final_result 載 verdict `partial`（可 grep）", "verdict `partial`" in ftxt)
        chk("A7c3 final_result 保留判死／禁區段",
            all(k in ftxt for k in ("判死", "禁區", "已否證")) or ("判死／禁區" in ftxt))
        chk("A7c4 final_result 四數與 audit 一致",
            all(k in ftxt for k in ("20/20", "109/109", "35/35")))
    txt = pr.read_text(encoding="utf-8")
    chk("A7d 問題記錄載明 🟡 狀態", "🟡" in txt)

    hdr("§A8 禁用 token 掃描（修正後不得再有未標註之舊值）")
    STALE = [
        ("1.339", ("更正", "舊", "審計", "歷史", "為準", "supersede", "錯誤", "約定", "混尺", "正解")),
        ("0.00443", ("更正", "舊", "審計", "歷史", "誤", "已廢", "正解")),
        ("1.239146", ("更正", "舊", "已廢", "初版", "錯", "1.324")),
        # 註：938.2720882 係 log 中 10 位有效數字之正常顯示，非舊值 ⟹ 不作 stale token
    ]
    bad_stale = []
    for f in sorted(list((ROOT / "docs").glob("*.md")) + [ROOT / "01_problem_record.md"]
                    + sorted(DATA.glob("run_*.log"))):
        for i, line in enumerate(f.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            for tok, allow in STALE:
                if tok in line and not any(a in line for a in allow):
                    bad_stale.append(f"{f.name}:{i}「{tok}」")
    chk("A8a 無未標註之舊值（1.339／0.00443）", not bad_stale, f"{bad_stale[:6]}")

    # A8b：L2 之硬編碼質量必須等於 vendored PDG 原文（避免「少 6 位」重演）
    PDGF = (ROOT.parent / "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt")
    IDS = {"2112": "m_n", "2212": "m_p", "3112": "m_Sm", "3212": "m_S0",
           "3222": "m_Sp", "3312": "m_Xm", "3322": "m_X0"}
    raw = {}
    if PDGF.exists():
        for line in PDGF.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("*") or len(line) < 70:
                continue
            cols = [line[0:8], line[8:16], line[16:24], line[24:32]]
            for pid, key in IDS.items():
                if key not in raw and any(c.strip() == pid for c in cols):
                    raw[key] = Decimal(line[33:51].strip())*1000    # GeV → MeV（十進位精確）
    src = (SCRIPTS / "l2_independent_verify.py").read_text(encoding="utf-8")
    mism = []
    for key, v in raw.items():
        m = re.search(r'"%s": "([0-9.eE+-]+)"' % key, src)
        if not m:
            mism.append((key, "L2 無此輸入"))
        elif Decimal(m.group(1)) != v:
            mism.append((key, f"L2={m.group(1)} vs PDG={v}"))
    chk("A8b L2 硬編碼質量 == vendored PDG 原文（7 個）", not mism, f"{mism}")

    hdr("§A9 跨記錄引用掃描（v11 main_problems，排除 venv／其他版本）")
    V11 = ROOT.parent.parent
    pat = re.compile(r"0\.403503649635|0\.025888367883|0\.377615281?7?5?2?|0\.044785178")
    refs = []
    for f in V11.rglob("*"):
        if not f.is_file() or f.suffix not in (".md", ".py", ".json", ".log"):
            continue
        if "m13_flavour_orthogonality" in str(f) or ".venv" in str(f):
            continue
        try:
            t = f.read_text(encoding="utf-8", errors="replace")
        except Exception:                                   # noqa: BLE001
            continue
        for i, line in enumerate(t.splitlines(), 1):
            if pat.search(line):
                refs.append(f"{f.relative_to(V11)}:{i}")
    chk("A9 跨記錄引用已列舉（供傳播表）", True, f"{len(refs)} 處: {refs[:6]}")

    hdr("§A10 時序閘（修正後下游必須重算）")
    late = []
    for p in L1_RUN:
        lg = CANON[p]
        if lg.exists() and lg.stat().st_mtime < p.stat().st_mtime - 1e-6:
            late.append((p.name, "log 早於腳本"))
    for js, sc in ((DATA / "l2_independent_verify.json", SCRIPTS / "l2_independent_verify.py"),
                   (DATA / "l3_doc_crosscheck.json", SCRIPTS / "l3_doc_crosscheck.py")):
        if js.exists() and js.stat().st_mtime < sc.stat().st_mtime - 1e-6:
            late.append((js.name, "json 早於腳本"))
    chk("A10 全部 log／json 喺其 producer 修改之後重算", not late, f"{late}")

    hdr("§匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}，WARN = {NWARN}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "warn": NWARN,
                       "n_scripts": len(L1),
                       "n_logs": len({LOGMAP[s] for s in on_disk if s in LOGMAP})}
    (DATA / "audit_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())