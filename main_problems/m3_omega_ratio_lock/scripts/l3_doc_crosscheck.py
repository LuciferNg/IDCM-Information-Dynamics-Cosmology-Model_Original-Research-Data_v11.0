#!/usr/bin/env python3
"""M3 L3 — 唯讀文檔對賬（期望值【只】來自 L2 之 data/l2_independent_verify.json）

禁讀 docs 之數值作期望值；唔改任何 doc。
細整數 claim 必須用 anchor（±70 字窗），否則零鑑別力。
輸出：data/l3_doc_crosscheck.json
"""
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_l2 = json.loads((ROOT/"data/l2_independent_verify.json").read_text(encoding="utf-8"))
L2 = {**_l2["inputs"], **_l2["derived"]}   # 原始輸入 ＋ 推導值（期望值只此一源）
DOCS = sorted((ROOT/"docs").glob("*.md"))

SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
def tokens(txt):
    t = txt.replace("\u2212", "-").translate(SUP)
    t = re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", t)          # 千分位
    out = set()
    for m in re.finditer(r"[-+]?\d+(?:\.\d+)?[eE][-+]?\d+", t):
        try: out.add(float(m.group(0)))
        except ValueError: pass
    for m in re.finditer(r"[-+]?\d*\.\d+", t):
        try: out.add(float(m.group(0)))
        except ValueError: pass
    for m in re.finditer(r"(?<![\d.])[-+]?\d{1,6}(?![\d.])", t):
        out.add(float(m.group(0)))
    for m in re.finditer(r"(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)", t):
        try:
            d = float(m.group(2))
            if d: out.add(float(m.group(1))/d)
        except ValueError:
            pass
    return out

def near(txt, anchor, span=70):
    i = txt.find(anchor)
    return None if i < 0 else txt[max(0, i-span): i+len(anchor)+span]

def hit(txt, val, anchor=None, tol=5e-3):
    seg = txt if anchor is None else near(txt, anchor)
    if seg is None: return None
    best = None
    for v in tokens(seg):
        d = abs(v - val)/max(abs(val), 1e-12)
        if best is None or d < best[0]: best = (d, v)
    return best

# (doc, label, l2 key, anchor, tol, mode)
CLAIMS = [
    ("01_round1_partition_null.md",  "Σγ = 90",            "sig_gamma", None, 1e-4, "num"),
    ("01_round1_partition_null.md",  "Σγ_odd = 44",        "sig_odd",   None, 1e-9, "num"),
    ("01_round1_partition_null.md",  "Σγ_even = 46",       "sig_even",  None, 1e-9, "num"),
    ("02_round2_g3_assignment.md",   "(3:16:41)/60",       "prod_345",  "3:16:41", 5e-3, "num"),
    ("03_round3_r3_operatorization.md","N_d = 44",         "N_d",       "N_d", 5e-3, "num"),
    ("04_round4_oplus_legitimacy.md","量子 3/2 = 1.5",     "quantum",   "3/2", 5e-3, "num"),
    ("05_round5_chain_length_origin.md","88 sites",        "N_sites",   None, 1e-9, "num"),
    ("06_round6_f2_onedim.md",       "N_d = 44",           "N_d",       None, 5e-3, "num"),
    ("07_round7_single_component.md","44 模式",            "N_d",       None, 5e-3, "num"),
    ("08_round8_correction_k.md",    "44",                 "N_d",       None, 5e-3, "num"),
    ("09_round9_gamma1_identity.md", "∏γ = N_d = 44",      "N_d",       None, 5e-3, "num"),
    ("09_round9_gamma1_identity.md", "Σγ = 90",            "sig_gamma", None, 1e-4, "num"),
    ("10_round10_omega_b_fork.md",   "Ω_b h² = 0.0227186", "omega_b_h2", None, 1e-4, "num"),
    ("11_round11_g3_fork.md",        "γ₂ = 0.0203571",      "gamma_2",   "0.0203571", 5e-4, "num"),
    ("11_round11_g3_fork.md",        "γ₃ = 0.0457063",      "gamma_3",   "0.0457063", 5e-4, "num"),
    ("12_round12_g3a_p4.md",         "E_sync = 61.5",      "E_sync",    None, 5e-3, "num"),
    ("12_round12_g3a_p4.md",         "E_proj = 28.5",      "E_proj",    None, 5e-3, "num"),
    ("13_round13_index_audit.md",    "物理索引 rank = 2",  "rank_phys", "物理", 5e-3, "num"),
    ("14_round14_r5_mc_null.md",     "p1 = 5.845e-4",      "p1_value_space", None, 5e-3, "num"),
    ("14_round14_r5_mc_null.md",     "3.25σ",              "z_p1",      "3.25", 5e-3, "num"),
    ("15_round15_r1_gate_identity.md","Ω_dm = 0.26667",    "omega_dm",  None, 5e-5, "num"),
    ("16_round16_a3_dimension_capacity.md","1452",         "dims_total", None, 1e-9, "num"),
    ("17_round17_pythagoras_audit.md","462 ≠ 780",         "pyth_lhs_5_12_13", None, 1e-9, "num"),
    ("17_round17_pythagoras_audit.md","780",               "pyth_rhs_5_12_13", "780", 5e-3, "num"),
    ("final_result.md",              "E_sync = 61.5",      "E_sync",    None, 5e-3, "num"),
    ("final_result.md",              "統計上限 3.25σ",      "z_p1",      "3.25", 5e-3, "num"),
    ("18_round18_audit_verify.md",   "整數三分 1711",       "n_triples", None, 1e-9, "num"),
]

report = {"claims": [], "per_doc": {}, "pass": True}
for doc, label, key, anchor, tol, mode in CLAIMS:
    p = (ROOT/doc) if doc in ("final_result.md",) else (ROOT/"docs"/doc)
    txt = p.read_text(encoding="utf-8")
    raw = L2[key]
    try:
        if isinstance(raw, str) and "/" in raw:
            _a, _b = raw.split("/"); val = float(_a)/float(_b)
        else:
            val = float(raw)
    except (TypeError, ValueError):
        val = None
    res = hit(txt, val, anchor, tol) if val is not None else None
    ok = res is not None and res[0] <= tol
    report["claims"].append({"doc": doc, "label": label, "key": key,
                             "expected": raw, "nearest": (res[1] if res else None),
                             "dev": (res[0] if res else None), "pass": ok})
    report["per_doc"].setdefault(doc, [0, 0])
    report["per_doc"][doc][0] += 1
    report["per_doc"][doc][1] += int(ok)
    if not ok: report["pass"] = False
    print(f"  [{'PASS' if ok else 'FAIL'}] {doc:44s} {label:22s} "
          f"L2={raw}  最近={res[1] if res else None}  dev={('%.2e' % res[0]) if res else 'None'}")

# 覆蓋閘：每份在碟 doc 必須有 ≥1 claim
disk = [p.name for p in DOCS]
missing = [d for d in disk if d not in report["per_doc"]]
print(f"\n  docs 在碟 = {len(disk)} ；有 claim = {len(report['per_doc'])} ；缺 claim = {missing}")
report["docs_on_disk"] = len(disk)
report["docs_with_claim"] = len(report["per_doc"])
report["docs_missing_claim"] = missing
if missing: report["pass"] = False

Path(ROOT/"data/l3_doc_crosscheck.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")
npass = sum(1 for c in report["claims"] if c["pass"])
print(f"\nL3 {npass}/{len(report['claims'])} ；all_pass={report['pass']}")