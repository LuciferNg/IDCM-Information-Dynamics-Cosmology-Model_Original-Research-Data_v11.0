#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — 文檔宣稱交叉核（idcm-audit-verify §2.4 / §3.4）
2026-10-29

三層：
  §1 CLAIMS   文檔字面值 → L2 獨立重推值（kind: num / last / multi / set）
  §2 REFUTED  已否證值須仍存在（append-only）且文檔須有更正節
  §3 溯源掃描  文檔內 ≥3 位小數之數值，須能在 machine artifact
              （data/*.json + data/*.log + scripts/*.py）中找到文檔精度下四捨五入相等之數
              ⟹ 否則 = 真孤兒（文檔有數、機器無源）

跑法：/usr/bin/python3 scripts/l3_doc_crosscheck.py
"""
import bisect
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOCS = ROOT / "docs"
DATA = ROOT / "data"

L2 = json.loads((DATA / "l2_independent_verify.json").read_text(encoding="utf-8"))
DER = {k: v["float"] for k, v in L2["_derived"].items()}

LOG, RES = [], {}
NFAIL, NMISS = 0, 0


def say(s=""):
    LOG.append(s)
    print(s)


def hdr(s):
    say("\n" + "=" * 80)
    say(s)
    say("=" * 80)


def norm(s):
    return s.replace("−", "-").replace("–", "-")


def nums_in(lit):
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", norm(lit))]


CLAIMS = [
    ("01_problem_statement.md", "1.1263", "num", ["O_np"], 2e-4),
    ("01_problem_statement.md", "1.9090", "num", ["O_S"], 2e-4),
    ("01_problem_statement.md", "2.1719", "num", ["O_X"], 2e-4),
    ("01_problem_statement.md", "2.4213", "num", ["O_K"], 2e-4),
    ("02_r1_domain_additivity.md", "1.1263", "num", ["O_np"], 2e-4),
    ("02_r1_domain_additivity.md", "1.9090", "num", ["O_S"], 2e-4),
    ("02_r1_domain_additivity.md", "2.1719", "num", ["O_X"], 2e-4),
    ("02_r1_domain_additivity.md", "2.4213", "num", ["O_K"], 2e-4),
    ("02_r1_domain_additivity.md", "14.11", "num", ["z_K"], 1e-3),
    ("02_r1_domain_additivity.md", "11.29", "num", ["z_X"], 1e-3),
    ("02_r1_domain_additivity.md", "5.77", "num", ["z_S"], 1e-3),
    ("02_r1_domain_additivity.md", "0.68", "num", ["z_np"], 5e-3),
    ("02_r1_domain_additivity.md", "0.1848", "num", ["sigO_np"], 5e-4),
    ("02_r1_domain_additivity.md", "0.1575", "num", ["sigO_S"], 5e-4),
    ("02_r1_domain_additivity.md", "0.1038", "num", ["sigO_X"], 5e-4),
    ("02_r1_domain_additivity.md", "0.1007", "num", ["sigO_K"], 5e-4),
    ("03_r2_method_gate.md", "3.223", "num", ["z_decay"], 1e-3),
    ("03_r2_method_gate.md", "0.683", "num", ["seg_rel"], 1e-3),
    ("03_r2_method_gate.md", "0.1848", "num", ["sigO_np"], 5e-4),
    ("04_r2a_two_body_exhaustion.md", "3.223", "num", ["z_decay"], 1e-3),
    ("05_r2b_resolution_density.md", "0.2640", "num", ["sigT_3ch"], 1e-3),
    ("05_r2b_resolution_density.md", "1.9090", "num", ["O_S"], 2e-4),
    ("06_r2c_readout_sigma_budget.md", "0.184766", "num", ["sigO_np"], 1e-5),
    ("06_r2c_readout_sigma_budget.md", "0.180078", "num", ["sigEm_np"], 1e-5),
    ("06_r2c_readout_sigma_budget.md", "0.041356", "num", ["floorX_np"], 1e-4),
    ("06_r2c_readout_sigma_budget.md", "0.1038", "num", ["sigO_X"], 5e-4),
    ("06_r2c_readout_sigma_budget.md", "2.1719", "num", ["O_X"], 2e-4),
    ("07_r2d_pairtype_counting.md", "14.11", "num", ["z_K"], 1e-3),
    ("07_r2d_pairtype_counting.md", "11.29", "num", ["z_X"], 1e-3),
    ("07_r2d_pairtype_counting.md", "5.77", "num", ["z_S"], 1e-3),
    ("08_r2e_splitting_ratios.md", "3.041756", "num", ["r_K"], 1e-6),
    ("08_r2e_splitting_ratios.md", "3.716756", "num", ["r_S"], 1e-6),
    ("08_r2e_splitting_ratios.md", "5.296397", "num", ["r_X"], 1e-6),
    ("08_r2e_splitting_ratios.md", "6.246656", "num", ["r_SS"], 1e-6),
    ("08_r2e_splitting_ratios.md", "0.015940", "num", ["sig_r_K"], 1e-5),
    ("08_r2e_splitting_ratios.md", "0.029705", "num", ["sig_r_S"], 1e-5),
    ("08_r2e_splitting_ratios.md", "0.163837", "num", ["sig_r_X"], 1e-5),
    ("08_r2e_splitting_ratios.md", "0.058885", "num", ["sig_r_SS"], 1e-5),
    ("08_r2e_splitting_ratios.md", "1.5796", "num", ["r_X_minus_r_S"], 1e-4),
    ("08_r2e_splitting_ratios.md", "0.6750", "num", ["r_S_plus_r_K"], 1e-4),
    ("09_r3p1_kappa_blindness.md", "−8", "num", ["kappa_T3"], 1e-9),
    ("09_r3p1_kappa_blindness.md", "−380", "num", ["kappa_T1"], 1e-9),
    ("09_r3p1_kappa_blindness.md", "| 1341 |", "num", ["kappa_nz_full"], 1e-9),
    ("09_r3p1_kappa_blindness.md", "10 / 1296", "num", ["kappa_K1_nonzero"], 1e-9),
    ("09_r3p1_kappa_blindness.md", "{1, 3, 21, 34}", "set", ["kappa_support"], 1e-9),
    ("10_r3p2_frustration_gauge.md", "−88 + 32 = −56", "last", ["rh_ising_H"], 1e-9),
    ("10_r3p2_frustration_gauge.md", "54 / 60", "num", ["rh_rank"], 1e-9),
    ("10_r3p2_frustration_gauge.md", "8 / 40 / 12", "multi",
     ["rh_hist_0", "rh_hist_1", "rh_hist_2"], 1e-9),
    ("11_r3p3_canonical_incidence.md", "Σn = 55", "num", ["sum_matter_n"], 1e-9),
    ("11_r3p3_canonical_incidence.md", "2.30906440147", "num", ["sum_gamma_9"], 1e-11),
    ("11_r3p3_canonical_incidence.md", "(14, 18, 23)", "multi",
     ["slot_col_1", "slot_col_2", "slot_col_3"], 1e-9),
    ("11_r3p3_canonical_incidence.md", "−81.8", "num", ["z_K1_slot2"], 1e-3),
    ("11_r3p3_canonical_incidence.md", "−22.3", "num", ["z_K1_slot3"], 1e-3),
    ("11_r3p3_canonical_incidence.md", "−74.3", "num", ["z_K2_slot2"], 1e-3),
    ("11_r3p3_canonical_incidence.md", "−18.6", "num", ["z_K2_slot3"], 1e-3),
]

REFUTED = [
    ("06_r2c_readout_sigma_budget.md", "0.2076",
     "M11 §14.6 判為導式錯（多乘一個 ⟨O⟩；正解 0.103777）", "六、更正"),
]


def main():
    global NFAIL, NMISS
    # docs/ ＋ 根目錄四件套（含 final_result.md）一併納入溯源
    texts = {p.name: p.read_text(encoding="utf-8") for p in DOCS.glob("*.md")}
    for p in sorted(ROOT.glob("*.md")):
        texts["(root)" + p.name] = p.read_text(encoding="utf-8")

    hdr("§1 CLAIMS：文檔字面值 → L2 獨立重推")
    npass = 0
    for doc, lit, kind, keys, tol in CLAIMS:
        txt = texts.get(doc)
        if txt is None:
            say("[FAIL] %s :: 檔案不存在" % doc)
            NFAIL += 1
            continue
        if lit not in txt:
            say("[MISS] %-34s %-22s 文檔無此字面值" % (doc, lit))
            NMISS += 1
            RES["%s|%s" % (doc, lit)] = {"status": "MISSING", "keys": keys}
            continue
        vs = nums_in(lit)
        bad = []
        if kind == "set":
            prod = 1.0
            for v in vs:
                prod *= v
            ref = abs(DER[keys[0]])
            if abs(prod - ref) > 1e-6 * max(ref, 1):
                bad.append("prod=%.1f ref=%.1f" % (prod, ref))
        elif kind == "last":
            ref = abs(DER[keys[0]])
            if not vs or abs(abs(vs[-1]) - ref) > tol * max(ref, 1):
                bad.append("doc=%s ref=%.10g" % (vs[-1] if vs else "?", ref))
        elif kind == "multi":
            if len(vs) < len(keys):
                bad.append("數字不足")
            for v, k in zip(vs, keys):
                ref = abs(DER[k])
                if abs(abs(v) - ref) > tol * max(ref, 1):
                    bad.append("%s: doc=%.6g ref=%.10g" % (k, v, ref))
        else:
            v = vs[0] if vs else None
            ref = abs(DER[keys[0]])
            if v is None:
                bad.append("無數字")
            elif abs(abs(v) - ref) > tol * max(ref, 1):
                bad.append("doc=%.6g ref=%.10g" % (v, ref))
        if bad:
            NFAIL += 1
            say("[FAIL] %-34s %-22s %s" % (doc, lit, "; ".join(bad)))
            RES["%s|%s" % (doc, lit)] = {"status": "FAIL", "keys": keys, "bad": bad}
        else:
            npass += 1
            RES["%s|%s" % (doc, lit)] = {"status": "PASS", "keys": keys}
    say("  CLAIMS %d：PASS %d / FAIL %d / MISSING %d" % (len(CLAIMS), npass, NFAIL, NMISS))

    hdr("§2 REFUTED：已否證值仍存在 + 文檔有更正節")
    for doc, lit, why, sect in REFUTED:
        txt = texts.get(doc, "")
        present = lit in txt
        fixed = sect in txt
        ok = present and fixed
        if not ok:
            NFAIL += 1
        say("[%s] %s :: %s 仍存=%s 更正節「%s」=%s ｜ %s"
            % ("PASS" if ok else "FAIL", doc, lit, present, sect, fixed, why))
        RES["REFUTED|%s|%s" % (doc, lit)] = {"status": "PASS" if ok else "FAIL"}

    hdr("§3 溯源掃描：文檔數值須可回溯至機器產物")
    parts = []
    skip = ("l3_", "run_l3", "15_full")
    for p in list(DATA.glob("*.json")) + list(DATA.glob("*.log")) + list(HERE.glob("*.py")):
        if p.name.startswith(skip):   # 排除審計自身產物（否則循環自證）
            continue
        parts.append(norm(p.read_text(encoding="utf-8", errors="ignore")))
    blob = "\n".join(parts)
    bset = set()
    for m in re.finditer(r"-?\d+\.\d+", blob):
        try:
            bset.add(abs(float(m.group(0))))
        except ValueError:
            pass
    blobs = sorted(bset)

    def traced(v, dec):
        w = 10 ** (-dec) * 0.51 + 1e-12
        lo = bisect.bisect_left(blobs, v - w)
        hi = bisect.bisect_right(blobs, v + w)
        for x in blobs[lo:hi]:
            if round(abs(x), dec) == round(v, dec):
                return True
        return False

    ACCEPTED = {
        "02_r1_domain_additivity.md": {
            "-2.0695": "跨記錄引用 M2 docs/02 §三（λ 下界）",
            "4.9538": "跨記錄引用 M2 docs/02 §三（兩模值）",
            "4.2257": "★ R1 內 inline 枚舉（34³）結果，未落地腳本 ⟹ 僅旁證，非判據",
        },
        "05_r2b_resolution_density.md": {
            "0.0026": "文內算術 0.1263/48（由已落地數推得）",
            "0.026": "文內算術 ~0.185/7（同上）",
            "0.079750": "M11 §14.6 之 |⟨O⟩|·σ_X/X（上游值，6dp）",
            "0.132744": "M11 §14.6 含 σ_Δm 之 σ_O(Ξ)（上游值）",
            "0.0798": "M11 §14.6 之 0.079750 取 4dp 寫法",
        },
        "06_r2c_readout_sigma_budget.md": {
            "0.079750": "M11 §14.6 之 |⟨O⟩|·σ_X/X（上游值，6dp）",
            "0.132744": "M11 §14.6 含 σ_Δm 之 σ_O(Ξ)（上游值）",
            "0.0798": "M11 §14.6 之 0.079750 取 4dp 寫法",
        },
        "(root)SESSION_START.md": {
            "0.0026": "轉引 docs/05（文內算術）",
            "0.0798": "轉引 M11 §14.6 之 4dp 寫法",
            "4.2257": "轉引 docs/02（inline 旁證，不作判據）",
        },
        "22_r1_r12_continuity_audit.md": {
            "4.2257": "轉引 docs/02 之 F6 已降級值（明文「旁證，不作判據」）⟹ 與 docs/02／SESSION_START 同類",
        },
        "(root)final_result.md": {
            "0.0798": "M11 §14.6 之 0.079750 取 4dp 寫法（L2 floorX_Ξ = 0.079749）",
        },
        "12_full_chain_audit.md": {
            "-2.0695": "轉引 docs/02（M2 跨記錄）",
            "4.9538": "轉引 docs/02（M2 跨記錄）",
            "4.2257": "轉引 docs/02（inline 旁證）",
            "0.0026": "轉引 docs/05（文內算術）",
            "0.026": "轉引 docs/05（文內算術）",
            "0.079750": "轉引 M11 §14.6 值",
            "0.132744": "轉引 M11 §14.6 值",
            "0.0798": "轉引 M11 §14.6 之 4dp 寫法",
        },
        "01_problem_statement.md": {
            "0.7120": "M11 C2 表 em 錨（上游輸入）",
            "2.9590": "同上",
            "2.2600": "同上",
            "0.120": "同上之 σ",
            "0.203": "同上之 σ",
            "0.117": "同上之 σ",
        },
    }
    untraced, tot, acc = {}, 0, {}
    for name, txt in sorted(texts.items()):
        miss = set()
        for m in re.finditer(r"-?\d+\.\d{3,}", norm(txt)):
            s = m.group(0)
            dec = len(s.split(".")[1])
            if not traced(abs(float(s)), dec):
                if s in ACCEPTED.get(name, {}):
                    acc[name + "|" + s] = ACCEPTED[name][s]
                else:
                    miss.add(s)
        if miss:
            untraced[name] = sorted(miss)
            tot += len(miss)
    say("  文檔 ≥3 位小數之數值：未溯源 = %d 個（機器數值池 %d 個）" % (tot, len(blobs)))
    for k, v in untraced.items():
        say("     %-38s %d：%s" % (k, len(v), ", ".join(v[:14]) + (" …" if len(v) > 14 else "")))
    say("  已分類接受（跨記錄引用／文內算術／inline 旁證）= %d 個：" % len(acc))
    for k, v in sorted(acc.items()):
        say("     ▪ %-46s %s" % (k, v))
    if tot:
        NFAIL += 1
    RES["_trace"] = {"untraced_total": tot, "by_doc": untraced,
                     "accepted": acc, "blob_numbers": len(blobs)}

    hdr("§4 匯總")
    say("  FAIL = %d  MISSING = %d  未溯源 = %d" % (NFAIL, NMISS, tot))
    RES["_summary"] = {"claims": len(CLAIMS), "pass": npass, "fail": NFAIL,
                       "missing": NMISS, "untraced": tot}
    (DATA / "l3_doc_crosscheck.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    say("")
    say("寫出: %s" % (DATA / "l3_doc_crosscheck.json"))
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())