#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 L3 -- 唯讀 doc 交叉核對（doc cross-check）。

硬規則（idcm-audit-verify §2）：
  * 唯讀：本腳本唔寫任何 doc
  * 期望值**只可**來自 L2（data/l2_independent_verify.json 之 derived 段）
  * 禁由 doc 內數字互檢（等於自我確認）

判據：每個 claim = (檔案, 文內出現之數值字串, L2 key)。容差依文內有效位：
  ≥6 位 1e−4 ／ ≥4 位 5e−4 ／ 其餘 5e−3；`tol` 可逐 claim 收緊。
  `excl` 類 = 已明文申報之「L2 未覆蓋」（pool 枚舉／上游產權），唔參與 PASS/FAIL，
  但會逐條列出，交付時必須報數。
輸出：data/l3_doc_crosscheck.json
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DOCS = os.path.join(ROOT, "docs")

with open(os.path.join(ROOT, "data", "l2_independent_verify.json"), encoding="utf-8") as fh:
    L2 = json.load(fh)["derived"]

# ---------------------------------------------------------------- claims
# (doc, 文內數值字串, L2 key, 額外容差 or None)
CLAIMS = [
    # ---- docs/37 三欄總表
    ("37_closeout_readiness.md", "784.96", "four_over_gamma1", None),
    ("37_closeout_readiness.md", "2.1498", "m_u_pred_MeV", None),
    ("37_closeout_readiness.md", "12.0055", "R_span", None),
    # ---- docs/11 FLAG 重判
    ("11_lattice_flag_span.md", "12.005533", "R_span", None),
    # σ_R 終局（2026-10-28）：canonical = 精確二階，由 L2 R_span_sigma_exact 生成
    ("11_lattice_flag_span.md", "0.039833", "R_span_sigma_exact", None),
    ("37_closeout_readiness.md", "0.039833", "R_span_sigma_exact", None),
    ("11_lattice_flag_span.md", "2.1498", "m_u_pred_MeV", None),
    # ---- docs/36 外部錨與判決
    ("36_r5_em_kaon_external_verdict.md", "1261.157", "dM2_pi", None),
    ("36_r5_em_kaon_external_verdict.md", "-3899.73", "dM2_K", None),
    ("36_r5_em_kaon_external_verdict.md", "2244.86", "dM2_K_gamma", None),
    ("36_r5_em_kaon_external_verdict.md", "2.2646", "em_lat_MeV", None),
    ("36_r5_em_kaon_external_verdict.md", "0.1209", "sigma_em", None),
    ("36_r5_em_kaon_external_verdict.md", "1.2722", "dashen_MeV", None),
    ("36_r5_em_kaon_external_verdict.md", "15.98", "z_single", None),
    ("36_r5_em_kaon_external_verdict.md", "14.43", "z_envelope", None),
    # ---- docs/29 單一形式
    ("29_operator_structure.md", "0.3327", "em_pred_single", None),
    ("29_operator_structure.md", "3.9340", "M_K0_minus_Kplus", None),
    # ---- docs/40（本輪自查更正）
    ("40_recheck_m7_m10_20261028.md", "2.421324", "O_K_backsolved", None),
    ("40_recheck_m7_m10_20261028.md", "3.524", "z_O_K", None),
    ("40_recheck_m7_m10_20261028.md", "2.131986", "n_s_required", None),
    ("40_recheck_m7_m10_20261028.md", "0.332667", "em_pred_single", None),
    ("40_recheck_m7_m10_20261028.md", "0.156667", "sigma_pred_X", None),
    ("40_recheck_m7_m10_20261028.md", "0.197869", "sigma_tot", None),
    ("40_recheck_m7_m10_20261028.md", "9.76", "z_with_X_budget", None),
    ("40_recheck_m7_m10_20261028.md", "2.560", "X_MeV", None),
    ("40_recheck_m7_m10_20261028.md", "0.094", "sigma_X_convention", None),
    # ---- docs/18 X 出處
    ("18_em_constraints.md", "2.56", "X_MeV", None),
]

# 已明文申報、L2 未覆蓋（唔參與 PASS/FAIL，但逐條報）
EXCLUDED = [
    ("docs/37 A4", "density 97.9% / 1.97% / 5.59% / 17 identity", "pool 枚舉（L1 scripts/03,09,16）；L2 未重建"),
    ("docs/37 C", "MC null p = 0.550 / 斜率 7.2× / 108 模型", "pool + MC（L1 scripts/31）；L2 未重建"),
    ("docs/37 A5", "O124 41.2σ / O85 R² = 0.310001 / a ≤ 0.9729 fm", "上游產權（O124/O85）— M4 只引用，不重建"),
    ("docs/37 C", "A = 1.64337 / 87γ₁⁴", "pool 枚舉（L1 scripts/09）"),
    ("docs/40 §4", "20.61% / 90 條", "M10 產權（docs/40 已明文：M4 未重建、不引用）"),
    ("docs/40 §4", "733.8× / 161.9%", "M10 W4 產權（M4 只獨立核比值，L2 未列 input）"),
    ("docs/36 §四", "1.78 / ε = 0.78(1)(+8/−11) / [0.67,0.86]", "外部錨（MILC18 → FLAG Eq.38），已作 L2 INPUT"),
]


def digits(s):
    s = s.lstrip("+-")
    return len(re.sub(r"[^0-9]", "", s.split("e")[0].split("E")[0]))


def tol_for(lit, override):
    if override is not None:
        return override
    d = digits(lit)
    if d >= 6:
        return 1e-4
    if d >= 4:
        return 5e-4
    return 5e-3


def norm(s):
    return s.replace("\u2212", "-").replace("\u00a0", " ")


results = []
missing = []
print("=" * 100)
print("M4 L3 doc 交叉核對（唯讀；期望值只來自 L2）")
print("=" * 100)
for doc, lit, key, ov in CLAIMS:
    path = os.path.join(DOCS, doc)
    txt = norm(open(path, encoding="utf-8").read())
    if lit not in txt:
        missing.append((doc, lit))
        print("MISS %-42s %-12s (字串不在檔內)" % (doc, lit))
        continue
    want = float(L2[key])
    got = float(lit.lstrip("+"))
    tol = tol_for(lit, ov)
    dev = abs(abs(got) - abs(want)) / abs(want) if want != 0 else abs(got)
    good = dev <= tol
    results.append({"doc": doc, "literal": lit, "l2_key": key, "l2_value": want,
                    "doc_value": got, "dev": dev, "tol": tol, "pass": good})
    print(("PASS " if good else "FAIL ") + "%-42s %-12s L2(%-22s)=%-20.10g dev=%.2e tol=%.0e"
          % (doc, lit, key, want, dev, tol))

npass = sum(1 for r in results if r["pass"])
nfail = [r for r in results if not r["pass"]]
print("\n" + "-" * 100)
print("L3 RESULT: %d/%d PASS（FAIL %d；MISS %d）" % (npass, len(results), len(nfail), len(missing)))
for r in nfail:
    print("  FAIL %s | %s vs L2 %s = %.10g | dev=%.2e > tol=%.0e"
          % (r["doc"], r["literal"], r["l2_key"], r["l2_value"], r["dev"], r["tol"]))
print("\n[EXCLUDED — L2 未覆蓋，已申報]")
for d, item, why in EXCLUDED:
    print("  %-14s %-52s %s" % (d, item, why))

out = os.path.join(ROOT, "data", "l3_doc_crosscheck.json")
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"summary": {"pass": npass, "total": len(results), "fail": len(nfail),
                           "miss": len(missing)},
               "claims": results, "missing": missing,
               "excluded": [{"doc": d, "item": i, "reason": w} for d, i, w in EXCLUDED]},
              fh, ensure_ascii=False, indent=2, sort_keys=True)
print("\nJSON -> data/l3_doc_crosscheck.json")
raise SystemExit(0 if (not nfail and not missing) else 1)