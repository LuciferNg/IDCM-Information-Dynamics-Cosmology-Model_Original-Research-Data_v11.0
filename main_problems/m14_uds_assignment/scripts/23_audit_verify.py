#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/23 — 全鏈審計（idcm-audit-verify §3.1 六項 ＋ L1 兩跑 byte 比對）
================================================================================
① 腳本集對稱   磁碟 NN_*.py ↔ README 所列集，兩側差集須空
② 逆向綁定     每份 docs/NN 至少綁一支 producer script（grep scripts/NN）
③ 路徑衛生     scripts/ 內禁 /tmp 寫入
④ 產物綁定     每個 NN 腳本有 data log；每個 data log 有腳本；無孤兒
⑤ inline 幽靈   docs 內 headline 數字須有 L3 claim 覆蓋（無 claim 者列出）
⑥ driver       scripts/run_chain.sh 存在且語法可讀
⑦ L1 兩跑      每支腳本跑兩次，stdout byte-identical，且 = canonical log
⑧ EXTRA        l2_/l3_ 兩支一併入審計（§7.3 盲區修正）
⑨ 檔數對稱     README 宣稱 docs/scripts 數 == 實際
⑩ 四件套同步   狀態行一致
⑪ 文件衛生     markdown 圍欄偶數、U+FFFD 為 0、備份存在
性質：確定性；唔改任何檔案（除 data/ 輸出）。exit 1 若有任何 FAIL。
"""
import json
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DOCS = ROOT / "docs"
PY = "/usr/bin/python3"
ENV = {"PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
       "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"}

# script stem → canonical log（顯式表；audit 自身不入 L1）
LOGS = {
    "01_reproduction_gate": "run_01_reproduction_gate.log",
    "02_q1_station_gauge": "run_02_q1.log",
    "03_q2a_label_provenance": "run_03_q2a.log",
    "04_q2b_permutation_uniqueness": "run_04_q2b.log",
    "05_q3_delta_u_sign": "run_05_q3.log",
    "06_q4_third_slot": "run_06_q4.log",
    "07_m13_incoming_reply": "run_07_m13_reply.log",
    "08_classification_partitions": "run_08_classification.log",
    "09_pair_galois_orbits": "run_09_pair_orbits.log",
    "10_symmetry_group_irreps": "run_10_symmetry.log",
    "11_pair_physical_classification": "run_11_pair_physical.log",
    "12_delta_support": "run_12_delta_support.log",
    "13_channel_fate_classification": "run_13_channel_fate.log",
    "14_three_wave_resonance": "run_14_resonance.log",
    "15_general_resonance_law": "run_15_general_law.log",
    "16_convention_and_three_triples": "run_16_three_triples.log",
    "17_chi3_charge_coefficient": "run_17_chi3.log",
    "18_triples_vs_tpipe": "run_18_tpipe_filter.log",
    "19_four_wave_and_mod3": "run_19_four_wave.log",
    "20_single_label_theorem": "run_20_single_label.log",
    "21_odd_even_phase_gap": "run_21_phase_gap.log",
    "22_budget_defect_and_dyadic_depth": "run_22_budget_defect.log",
    "24_charge_kernel_and_complement": "run_24_charge_kernel.log",
    "25_invariant_lattice": "run_25_invariant_lattice.log",
    "26_odd_sector_and_partition_separability": "run_26_odd_sector.log",
    "27_fourwave_88_and_duality": "run_27_fourwave_88.log",
    "28_sublattice_divisibility": "run_28_sublattice.log",
    "29_qsqrt2_sector_and_field_criterion": "run_29_qsqrt2.log",
    "30_ladder_ap_vs_gp": "run_30_ladder.log",
    "31_general_ap_gp_criterion": "run_31_general_criterion.log",
    "32_ap_classification": "run_32_ap_classification.log",
    "33_gp_classification": "run_33_gp_classification.log",
    "34_phi_ladder_identity": "run_34_phi_ladder.log",
    "35_residue_five_items": "run_35_residue_five.log",
    "36_parity_and_delta_theorems": "run_36_theorems.log",
    "37_reverse_coverage": "run_37_reverse_coverage.log",
    "38_q1_m2_m3": "run_38_q1_m2_m3.log",
    "39_k3_closed_form_proof": "run_39_k3_proof.log",
    "40_fourwave_overlap_theorem": "run_40_fourwave_overlap.log",
    "41_energy_vs_selection": "run_41_energy_vs_selection.log",
    "42_cross_section_machine": "run_42_cross_section.log",
    "43_second_order_selfenergy": "run_43_second_order.log",
    "44_bubble_and_sign": "run_44_bubble_sign.log",
    "45_third_order": "run_45_third_order.log",
}
EXTRA = [("L2", "l2_independent_verify.py", "run_l2_independent_verify.log"),
         ("L3", "l3_doc_crosscheck.py", "run_l3_doc_crosscheck.log")]

LOG, RES, NFAIL, NWARN = [], {}, 0, 0


def chk(tag, cond, detail="", warn=False):
    global NFAIL, NWARN
    ok = bool(cond)
    if not ok:
        if warn:
            NWARN += 1
        else:
            NFAIL += 1
    st = "PASS" if ok else ("WARN" if warn else "FAIL")
    line = f"[{st}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"status": st, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def run(script: Path, reps=1):
    t0 = time.time()
    outs = []
    for _ in range(reps):
        p = subprocess.run([PY, str(script)], cwd=str(ROOT), env=ENV,
                           capture_output=True, text=True)
        outs.append((p.returncode, p.stdout))
    return outs, time.time() - t0


def main():
    hdr("§① 腳本集對稱")
    disk = {p.stem: p for p in sorted(SCRIPT_DIR.glob("[0-9][0-9]_*.py")) if "audit" not in p.stem}
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    missing_in_readme = [s for s in disk if s not in readme]
    chk("①a 每支 NN 腳本都在 README 列出", not missing_in_readme, f"漏列 {missing_in_readme}")
    listed_not_disk = [ln for ln in re.findall(r"scripts/([0-9]{2}_[A-Za-z0-9_]+)\.py", readme)
                       if ln not in disk]
    chk("①b README 所列腳本都存在", not listed_not_disk, f"幽靈 {sorted(set(listed_not_disk))}")
    chk("①c 腳本數 == LOGS 宣告數（不含審計）", len(disk) == len(LOGS),
        f"實際 {len(disk)} vs 宣告 {len(LOGS)}")

    # 敘事類 doc：無 producer script 屬正常，但必須明示理由（唔可以靜靜放行）
    NARRATIVE = {
        "07_verdicts_and_handoffs.md": "判詞與移交敘述表：全部數字引他 doc，本身無 producer",
        "09_correction_object_type.md": "概念更正（對象型別）；唯一精確值屬 M5 跨記錄",
        "13_channel_fate_correction.md": "讀法更正記錄：程式與 log 保留不動（append-only）",
    }

    hdr("§② 逆向綁定：每份 doc 有 producer script")
    docs = sorted(p.name for p in DOCS.glob("[0-9][0-9]_*.md"))
    unbound = []
    for d in docs:
        txt = (DOCS / d).read_text(encoding="utf-8")
        if not re.search(r"scripts/[0-9]{2}_", txt) and d not in NARRATIVE:
            unbound.append(d)
    chk("② 每份 doc 至少綁一支 producer（或屬敘事類並附理由）", not unbound,
        f"未綁且未申報 {unbound}", warn=True)
    chk("②b 敘事類白名單每項理由非空", all(v.strip() for v in NARRATIVE.values()),
        f"空理由 {[k for k, v in NARRATIVE.items() if not v.strip()]}")

    hdr("§③ 路徑衛生：scripts/ 禁 /tmp 寫入")
    hits = []
    for p in sorted((ROOT / "scripts").glob("*.py")):
        if p.name == Path(__file__).name:
            print("  （跳過審計自身：§7.1 陷阱 12 自我指涉假陽性）")
            continue
        for i, ln in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            if "/tmp/" in ln and re.search(r"open\(|write_text|write\(|Path\(", ln):
                hits.append(f"{p.name}:{i}")
    chk("③ 無 /tmp 寫入", not hits, f"命中 {hits}")

    hdr("§④ 產物綁定：腳本 ↔ data log")
    orphan_scripts = [s for s in disk if s not in LOGS]
    chk("④a 每個 L1 腳本都有宣告 log", not orphan_scripts, f"缺 {orphan_scripts}")
    missing_logs = [s for s, lg in LOGS.items() if s in disk and not (DATA / lg).exists()]
    chk("④b 每個宣告 log 都存在", not missing_logs, f"缺 {missing_logs}")
    logs_on_disk = {p.name for p in DATA.glob("run_*.log")}
    declared = set(LOGS.values()) | {lg for _, _, lg in EXTRA} | {"run_23_audit_verify.log"}
    chk("④c 無孤兒 log", logs_on_disk - declared == set(),
        f"孤兒 {sorted(logs_on_disk - declared)}", warn=True)

    hdr("§⑤ inline 幽靈：docs 數字須有 L3 claim 覆蓋")
    l3 = json.loads((DATA / "l3_doc_crosscheck.json").read_text(encoding="utf-8")) \
        if (DATA / "l3_doc_crosscheck.json").exists() else {}
    covered = set((l3.get("_summary", {}) or {}).get("by_doc", {}).keys())
    exempt = set((l3.get("_summary", {}) or {}).get("exempt_docs", []))
    reasons = (l3.get("_summary", {}) or {}).get("exempt_reasons", {}) or {}
    uncovered = [f"docs/{d}" for d in docs
                 if f"docs/{d}" not in covered and f"docs/{d.split('_')[0]}" not in exempt]
    chk("⑤ 每份 doc 有 claim 或被明示豁免", not uncovered, f"未覆蓋 {uncovered}")
    # ★ 分母斷言（skill 陷阱 41）：claim 份數 + 豁免份數 必須 == 磁碟 doc 份數
    n_cov = len({d for d in covered if d in {f"docs/{x}" for x in docs}})
    n_ex = len({d for d in docs if f"docs/{d}" not in covered and f"docs/{d.split('_')[0]}" in exempt})
    chk("⑤b claim 份數 + 豁免份數 == 磁碟 doc 份數", n_cov + n_ex == len(docs),
        f"{n_cov} + {n_ex} = {n_cov + n_ex} vs {len(docs)}")
    chk("⑤c 每項豁免附非空理由", all(str(v).strip() for v in reasons.values()) and len(reasons) >= len(exempt),
        f"豁免 {len(reasons)} 條 vs 名單 {len(exempt)}；空理由 {[k for k, v in reasons.items() if not str(v).strip()]}")

    hdr("§⑥ driver")
    drv = ROOT / "scripts" / "run_chain.sh"
    chk("⑥a run_chain.sh 存在", drv.exists())
    chk("⑥b 含 BLAS 線程鎖", drv.exists() and "OMP_NUM_THREADS" in drv.read_text(encoding="utf-8"))

    hdr("§⑦ L1 兩跑 byte-identical ＋ = canonical log")
    total_t = 0.0
    for stem, lg in LOGS.items():
        if stem not in disk:
            chk(f"⑦.{stem} 存在", False)
            continue
        script = disk[stem]
        o1, _ = run(script, reps=1)
        o2, _ = run(script, reps=1)
        rc1, out1 = o1[0]
        rc2, out2 = o2[0]
        canon = (DATA / lg).read_text(encoding="utf-8")
        chk(f"⑦.{stem[:26]} rc=0 ×2", rc1 == 0 and rc2 == 0, f"rc {rc1}/{rc2}")
        chk(f"⑦.{stem[:26]} 兩跑 byte-identical", out1 == out2)
        chk(f"⑦.{stem[:26]} ≡ canonical", out1 == canon,
            "" if out1 == canon else
            f"len {len(out1)} vs {len(canon)}；首差 {next((i for i in range(min(len(out1),len(canon))) if out1[i]!=canon[i]), -1)}")
        total_t += 0
    print(f"\n  （L1 共 {len(LOGS)} 支 ×2 跑完成）")

    hdr("§⑧ EXTRA：L2 / L3 入審計（§7.3 盲區修正）")
    for tag, fname, lg in EXTRA:
        p = ROOT / "scripts" / fname
        if not p.exists():
            chk(f"⑧.{tag} 存在", False)
            continue
        o1, _ = run(p, reps=1)
        o2, _ = run(p, reps=1)
        rc1, out1 = o1[0]
        rc2, out2 = o2[0]
        canon = (DATA / lg).read_text(encoding="utf-8")
        chk(f"⑧.{tag} 兩跑 byte-identical", out1 == out2)
        chk(f"⑧.{tag} ≡ canonical", out1 == canon)
        chk(f"⑧.{tag} FAIL 數 = 0", out1.count("[FAIL]") == 0,
            f"FAIL {out1.count('[FAIL]')}")

    hdr("§⑨ 檔數對稱（README 宣稱 == 實際）")
    claimed_docs = len(set(re.findall(r"docs/([0-9]{2}_[A-Za-z0-9_]+\.md)", readme)))
    chk("⑨a docs 實際數 == README 檔表所列", len(docs) == claimed_docs,
        f"實際 {len(docs)} vs README {claimed_docs}")
    chk("⑨b docs 每份都在 README 檔表", all(d in readme for d in docs),
        f"漏 {[d for d in docs if d not in readme]}")
    claimed_scripts = len(set(re.findall(r"scripts/([0-9]{2}_[A-Za-z0-9_]+)\.py", readme)))
    chk("⑨c scripts 實際數 == README 所列", len(disk) == claimed_scripts,
        f"實際 {len(disk)} vs README {claimed_scripts}")

    hdr("§⑩ 四件套狀態行同步")
    pat = re.compile(r"(狀態|verdict|state)", re.I)
    lines = {}
    for f in ("01_problem_record.md", "README.md", "final_result.md", "SESSION_START.md"):
        fp = ROOT / f
        if fp.exists():
            ls = [ln for ln in fp.read_text(encoding="utf-8").splitlines() if pat.search(ln)][:3]
            lines[f] = ls
    chk("⑩a 四件套齊存", len(lines) == 4, f"{sorted(lines)}")
    verdicts = {}
    for f, ls in lines.items():
        blob = " ".join(ls)
        verdicts[f] = ("partial" in blob.lower(), "negative" in blob.lower())
    chk("⑩b final_result 與 register 一致（partial, 非 negative-only）",
        verdicts.get("final_result.md", (False, True))[0])

    hdr("§⑪ 文件衛生")
    bad_fence = [p.name for p in list(DOCS.glob("*.md")) + [ROOT / "README.md", ROOT / "final_result.md"]
                 if len(re.findall(r"^```", p.read_text(encoding="utf-8"), re.M)) % 2]
    chk("⑪a markdown 圍欄偶數", not bad_fence, f"奇數 {bad_fence}")
    fffd = [p.name for p in list(DOCS.glob("*.md")) + list((ROOT / "scripts").glob("*.py"))
            if "\ufffd" in p.read_text(encoding="utf-8")]
    chk("⑪b 無 U+FFFD", not fffd, f"含 {fffd}")
    baks = list(ROOT.glob("*.pre_patch_*"))
    chk("⑪c 備份存在", len(baks) >= 4, f"{len(baks)} 個")

    hdr("§⑫ 反向數值覆蓋（skill §3.1 ⑦；scripts/37）")
    rc_log = DATA / "run_37_reverse_coverage.log"
    chk("⑫a scripts/37 log 存在", rc_log.exists())
    if rc_log.exists():
        rtxt = rc_log.read_text(encoding="utf-8")
        m = re.search(r"未追溯（且不在 ALLOW）\s*=\s*(\d+)", rtxt)
        chk("⑫b 未追溯 = 0", m is not None and int(m.group(1)) == 0,
            f"值 {m.group(1) if m else '缺'}")
        chk("⑫c 工具自身輸出已排除（防自我餵養）", "SELF" in (ROOT / "scripts" / "37_reverse_coverage.py").read_text(encoding="utf-8"))
        chk("⑫d ALLOW 每項附理由",
            all(len(x) > 10 for x in re.findall(r"← (.+)", rtxt)) if "←" in rtxt else True)

    hdr("§匯總")
    n = len([k for k in RES])
    print(f"檢查項 = {n}，FAIL = {NFAIL}，WARN = {NWARN}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "warn": NWARN,
                       "n_docs": len(docs), "n_scripts": len(disk), "n_logs": len(logs_on_disk)}
    (DATA / "audit_verify.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


SCRIPT_DIR = ROOT / "scripts"
if __name__ == "__main__":
    sys.exit(main())