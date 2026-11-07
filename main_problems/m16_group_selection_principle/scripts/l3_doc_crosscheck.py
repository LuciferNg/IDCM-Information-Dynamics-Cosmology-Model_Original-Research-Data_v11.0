#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — M16 文檔交叉核（idcm-audit-verify §2.4）
================================================================
鐵律
  R1 唯讀：對 docs／scripts 只讀，只寫 data/l3_doc_crosscheck.json 與 data/run_l3.log
  R2 期望值只來自 L2（data/l2_independent_verify.json 之 _derived）
  R3 容差依文檔有效位：≥6 位 → 1e−4；≥4 位 → 5e−4；其餘 → 5e−3
  R4 負值以絕對值比
  R6 細整數（3／4／12／50／60）必須帶 anchor，只喺錨點 ±span 窗口抽數
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG, RES, NFAIL = [], {}, 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


# ══════════════════════════════════════════════════════════════
# ① CLAIMS
# ══════════════════════════════════════════════════════════════
CLAIMS = [
    # ── docs/01 算術歸約
    {"doc": "docs/01_r5_and_arithmetic_reduction.md", "key": "group_order", "label": "∏{3,4,5} = 60"},
    {"doc": "docs/01_r5_and_arithmetic_reduction.md", "key": "sum_deg", "label": "Σ = 12"},
    {"doc": "docs/01_r5_and_arithmetic_reduction.md", "key": "n_sets_sum12", "label": "Σ=12 之集合恰 8 個", "anchor": "恰好"},
    {"doc": "docs/01_r5_and_arithmetic_reduction.md", "key": "sum_deg2_dnt", "label": "Σd² = 50", "anchor": "唯一）"},
    {"doc": "docs/01_r5_and_arithmetic_reduction.md", "key": "sumd_plus_deep", "label": "Σd + N_d = 60"},
    # ── docs/02 掃描與修飾詞
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "n_groups", "label": "92,570 群", "anchor": "群"},
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "hits_prod", "label": "∏=|G| 單獨 17 群", "anchor": "17"},
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "hits_sum12", "label": "Σ=12 單獨 267 群", "anchor": "267"},
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "v1_no_distinct", "label": "去「相異」5 群 ∌A₅", "anchor": "5 群"},
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "v2_no_nontrivial", "label": "去「非平凡」586 群", "anchor": "586"},
    {"doc": "docs/02_global_scan_and_modifiers.md", "key": "family_frac", "label": "準則族 17/45 = 37.8%", "anchor": "37.8", "scale": 100},
    # ── docs/03 判準來源與判級
    {"doc": "docs/03_criterion_source_and_verdict.md", "key": "n_groups_order60", "label": "階 60 之 13 群", "anchor": "13"},
    {"doc": "docs/03_criterion_source_and_verdict.md", "key": "le_56", "label": "look-elsewhere 1.79%", "anchor": "1.79", "scale": 100},
    {"doc": "docs/03_criterion_source_and_verdict.md", "key": "deep_sector", "label": "深層 420", "anchor": "420"},
    # ── docs/04 判詞
    {"doc": "docs/04_verdict_r5_and_status.md", "key": "group_order", "label": "|G| = 60", "anchor": "|"},
    {"doc": "docs/04_verdict_r5_and_status.md", "key": "n_groups_order60", "label": "13 群", "anchor": "13 群（SmallGrp"},
    {"doc": "docs/04_verdict_r5_and_status.md", "key": "le_56", "label": "look-elsewhere 1.79%", "anchor": "1.79", "scale": 100},
    # ── docs/05 最低成本
    {"doc": "docs/05_minimum_cost_principle.md", "key": "dnt345_groups", "label": "DNT={3,4,5} 之 11 群", "anchor": "11"},
    {"doc": "docs/05_minimum_cost_principle.md", "key": "deep_sector", "label": "深層 420", "anchor": "420"},
    {"doc": "docs/05_minimum_cost_principle.md", "key": "dnt345_superset", "label": "含 {3,4,5} 之 24 群", "anchor": "24"},
    {"doc": "docs/05_minimum_cost_principle.md", "key": "dnt345_superset_min_order", "label": "最小階 60", "anchor": "min |G|"},
    # ── docs/06 變分公理（γ / Δ 精確表）
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_1", "label": "γ_1 = 0.15224093…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_2", "label": "γ_2 = 0.58578643…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_3", "label": "γ_3 = 1.23463313…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_5", "label": "γ_5 = 2.76536686…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_6", "label": "γ_6 = 3.41421356…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "gamma_7", "label": "γ_7 = 3.84775906…", "tol": 1e-6},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "delta_1", "label": "Δ_1 = +0.281305", "tol": 1e-5},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "delta_2", "label": "Δ_2 = +0.215301", "tol": 1e-5},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "delta_3", "label": "Δ_3 = +0.116520", "tol": 1e-5},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "N_d_val", "label": "中項帶 nR = N_d = 44", "anchor": "nR = N_d"},
    {"doc": "docs/06_variational_axiom_attack.md", "key": "cubic_B6", "label": "B=6 → 210", "anchor": "210"},
    # ── docs/07 MF2
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "omega_b_planck", "label": "Ω_b = 0.049302", "tol": 1e-5},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "sigma_omega_b", "label": "σ(Ω_b) = 0.000857", "tol": 1e-3},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "k_center", "label": "k = 1.0142", "tol": 1e-4},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "k_sigma", "label": "σ(k) = 0.0176"},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "z_k2", "label": "z(k=2) = −55.93", "tol": 1e-4},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "z_k3", "label": "z(k=3) = −112.67", "tol": 1e-4},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "z_k4", "label": "z(k=4) = −169.41", "tol": 1e-4},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "k_with_h_idcm", "label": "k(H₀=0.674071) = 1.0156", "tol": 1e-4},
    {"doc": "docs/07_mf2_closed_mf1prime_narrowed.md", "key": "le_56", "label": "1/56 → 1.79%（級聯）", "anchor": "1.79", "scale": 100},
    # ── docs/08 MF1′ 撤回
    {"doc": "docs/08_mf1prime_dissolution.md", "key": "n_sets_sum12", "label": "8 集", "anchor": "8 個"},
    {"doc": "docs/08_mf1prime_dissolution.md", "key": "n_surv_L2", "label": "生還 1 集", "anchor": "唯一"},
    {"doc": "docs/08_mf1prime_dissolution.md", "key": "le_8", "label": "1/8 = 12.5%", "anchor": "12.5", "scale": 100},
    {"doc": "docs/08_mf1prime_dissolution.md", "key": "le_56", "label": "1/56 = 1.79%", "anchor": "1.79", "scale": 100},
    {"doc": "docs/08_mf1prime_dissolution.md", "key": "gamma_3_plus_5", "label": "γ_3+γ_5 = 4", "anchor": "γ_3"},
]

TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")
COMMA = re.compile(r"(?<=\d),(?=\d{3})")


def tokens(text):
    out = []
    text = COMMA.sub("", text.replace("\u2212", "-"))
    for m in TOK.finditer(text):
        mant, exp = m.group(1), m.group(2)
        try:
            v = float(mant + ("e" + exp if exp else ""))
        except ValueError:
            continue
        sig = len(mant.replace(".", "").lstrip("0"))
        out.append((v, sig, m.start()))
    return out


def window(text, anchor, span):
    i = text.find(anchor)
    if i < 0:
        return None
    return text[max(0, i - span): i + len(anchor) + span]


def main():
    L2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    D = L2.get("_derived", {})
    chk("L2 存在且 FAIL=0", L2.get("_summary", {}).get("fail") == 0,
        f"FAIL={L2.get('_summary', {}).get('fail')}")

    hdr("§1 逐 claim 核對")
    per_doc = {}
    for c in CLAIMS:
        path = ROOT / c["doc"]
        txt = path.read_text(encoding="utf-8")
        key = c["key"]
        if key not in D:
            chk(f"{c['doc']} :: {c['label']}", False, f"L2 缺 key {key}")
            continue
        exp = D[key]["float"]
        if exp is None:
            exp = float(D[key]["str"].replace("/", " "))
        exp = exp * float(c.get("scale", 1))
        seg = txt
        if "anchor" in c:
            seg = window(txt, c["anchor"], c.get("span", 70))
            if seg is None:
                chk(f"{c['doc']} :: {c['label']}", False, f"錨點搵唔到：{c['anchor']}")
                continue
        best, bestrel, besttol = None, 9e9, 0.0
        for v, sig, _ in tokens(seg):
            if exp == 0:
                continue
            rel = abs(abs(v) - abs(exp)) / abs(exp)
            sig_eff = max(1, min(sig, 8))
            tol = c.get("tol") or (1e-4 if sig_eff >= 6 else (5e-4 if sig_eff >= 4 else 5e-3))
            if rel < bestrel:
                best, bestrel, besttol = v, rel, tol
        ok = best is not None and bestrel <= besttol
        chk(f"{c['doc'].split('/')[-1]} :: {c['label']}", ok,
            f"exp={exp:.10g} best={best} rel={bestrel:.3e} tol={besttol:.0e}")
        per_doc[c["doc"]] = per_doc.get(c["doc"], [0, 0])
        per_doc[c["doc"]][0] += 1
        per_doc[c["doc"]][1] += 1 if ok else 0

    hdr("§2 逐檔命中率（須 100%）")
    for d, (n, k) in sorted(per_doc.items()):
        chk(f"{d.split('/')[-1]} 命中 {k}/{n}", k == n, f"{100*k/n:.1f}%")

    hdr("§3 負控制（期望值 +1 必須 FAIL）")
    c = CLAIMS[0]
    exp = D[c["key"]]["float"]
    txt = (ROOT / c["doc"]).read_text(encoding="utf-8")
    bad = 0
    for v, sig, _ in tokens(txt):
        if abs(abs(v) - abs(exp + 1)) / abs(exp + 1) <= 1e-4:
            bad += 1
    chk("N1 期望值+1 在文中無匹配（負控制）", bad == 0, f"命中 {bad}")

    hdr("§4 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": nchk, "fail": NFAIL,
                       "docs": {d: {"claims": n, "hit": k} for d, (n, k) in per_doc.items()}}
    (DATA / "l3_doc_crosscheck.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2),
                                                 encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / 'l3_doc_crosscheck.json'}\n寫出: {DATA / 'run_l3.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())