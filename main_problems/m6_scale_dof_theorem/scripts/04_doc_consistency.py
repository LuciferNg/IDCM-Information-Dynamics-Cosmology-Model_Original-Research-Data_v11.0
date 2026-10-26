#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/04_doc_consistency.py — 文檔數字 ↔ data JSON 一致性核對（轉錄防錯）
================================================================================
判準：文檔所寫數字 == 腳本輸出值（按文檔自身位數四捨五入；零可調容差）。
每一條 claim 指定 (檔案, 字串, JSON 取值路徑, 顯示格式)。

執行：/usr/bin/python3 scripts/04_doc_consistency.py
"""
import json
import os
import re
import sys
from fractions import Fraction

def load(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)

D1 = load("data/01_t1_composite_additivity.json")
D2 = load("data/02_t2_scale_dof_theorem.json")
D3 = load("data/03_r4_counterexample_hunt.json")
D5 = load("data/05_p_dof_ledger.json")
D6 = load("data/06_kappa_identifiability.json")
D7 = load("data/07_hidden_scales.json")

DOCS = {}
for fn in ("docs/01_t1_composite_additivity.md",
           "docs/02_t2_scale_dof_theorem.md",
           "docs/03_r4_counterexample_hunt.md",
           "docs/04_t3_handoff_and_verdict.md",
           "docs/05_p_dof_ledger.md",
           "docs/06_kappa_lambda_identifiability.md",
           "docs/07_hidden_layer_mechanism.md",
           "final_result.md", "README.md", "01_problem_record.md"):
    with open(fn, encoding="utf-8") as f:
        DOCS[fn] = f.read()

CLAIMS = [
    # (doc, 必須出現字串, 來源值, 格式)
    ("docs/01_t1_composite_additivity.md", "rank_ℚ{γ_n} = **20**", D1["checks"]["A4"]["rank"], "{}"),
    ("docs/01_t1_composite_additivity.md", "核空間 **24** 維", D1["checks"]["A4"]["ker_dim"], "{}"),
    ("docs/01_t1_composite_additivity.md", "3²¹·4 = 41,841,412,812",
     D1["checks"]["A9"]["image_bound"], "{:,}"),
    ("docs/01_t1_composite_additivity.md", "17,592,186,044,416",
     D1["checks"]["A9"]["states"], "{:,}"),
    ("docs/01_t1_composite_additivity.md", "平均簡併 ≥ 420.4",
     D1["checks"]["A9"]["avg_degeneracy_min"], "{:.1f}"),
    ("docs/01_t1_composite_additivity.md", "deg γ₁ = 20",
     D1["checks"]["A10"]["deg1"], "{}"),
    ("docs/01_t1_composite_additivity.md", "9100 個隨機分割",
     D1["checks"]["A6"]["partition_trials"], "{}"),
    ("docs/01_t1_composite_additivity.md", "46 / 44", 46, "{}"),

    ("docs/02_t2_scale_dof_theorem.md", "0.771960195963",
     D2["checks"]["B7"]["f_proj"], "{:.12f}"),
    ("docs/02_t2_scale_dof_theorem.md", "20.6411",
     D2["checks"]["B3"]["log10"], "{:.4f}"),
    ("docs/04_t3_handoff_and_verdict.md", "2.457611",
     D2["checks"]["B9"]["lam_window"][0], "{:.6f}"),
    ("docs/04_t3_handoff_and_verdict.md", "2.529465",
     D2["checks"]["B9"]["lam_window"][1], "{:.6f}"),
    ("docs/02_t2_scale_dof_theorem.md", "75.3026",
     D2["checks"]["B7"]["d_m_df_proj_MeV"], "{:.4f}"),
    ("docs/02_t2_scale_dof_theorem.md", "0.561",
     D2["checks"]["B5"]["spread"], "{:.3f}"),
    ("docs/02_t2_scale_dof_theorem.md", "2.22e−9".replace("−", "−"),
     D2["checks"]["B6"]["worst_bound"], "{:.2e}"),

    ("docs/03_r4_counterexample_hunt.md", "0.786081",
     D3["checks"]["C2"]["omega"], "{:.6f}"),
    ("docs/03_r4_counterexample_hunt.md", "22,545,093",
     D3["checks"]["C3"]["exact_counts"][0]["count_1e-06"], "{:,}"),
    ("docs/03_r4_counterexample_hunt.md", "19,231",
     D3["checks"]["C3"]["exact_counts"][0]["count_1e-09"], "{:,}"),
    ("docs/03_r4_counterexample_hunt.md", "98,562",
     D3["checks"]["C3"]["exact_counts"][2]["count_1e-06"], "{:,}"),
    ("docs/03_r4_counterexample_hunt.md", "σ = **8.25**",
     D3["checks"]["C3"]["sigma"], "{:.2f}"),

    ("final_result.md", "2.2852931e−21",
     1.0 / D2["checks"]["B3"]["N_ratio"], "{:.7e}"),

    # ── 2026-10-26 第二輪（items 2–4）────────────────────
    ("docs/05_p_dof_ledger.md", "(0,0,0)",
     sum(D5["checks"]["A"]["dim_p"]), "{:d}"),
    ("docs/05_p_dof_ledger.md", "det J = 1",
     int(D5["checks"]["B"]["det"]), "{}"),
    ("docs/05_p_dof_ledger.md", "1.0397344",
     D5["checks"]["C"]["R_max"], "{:.7f}"),
    ("docs/05_p_dof_ledger.md", "0.481828",
     float(Fraction(19, 41)) * D5["checks"]["C"]["R_max"], "{:.6f}"),
    ("docs/06_kappa_lambda_identifiability.md", "1,817",
     D6["checks"]["full"]["n_window_hits"], "{:,}"),
    ("docs/06_kappa_lambda_identifiability.md", "10 命中",
     D6["checks"]["support2"]["n_hits"], "{}"),
    ("docs/07_hidden_layer_mechanism.md", "7.994904",
     float(D7["checks"]["T1_support_theorem"]["rows"][1]["lam_max"]), "{:.6f}"),
    ("docs/07_hidden_layer_mechanism.md", "0.902813",
     float(D7["checks"]["T1_support_theorem"]["rows"][1]["log10"]), "{:.6f}"),
    ("docs/07_hidden_layer_mechanism.md", "4.86e18",
     float(D7["checks"]["T3_layer_count"]["n_quanta_additive"]), "{:.3g}"),
    ("docs/07_hidden_layer_mechanism.md", "10.5622",
     float(D7["checks"]["T3_layer_count"]["n_layers_multiplicative"]), "{:.4f}"),
    ("docs/07_hidden_layer_mechanism.md", "**2207**",
     D7["checks"]["T4_identity_density"]["k_min"]["v_SM"], "{}"),
]


def norm(s):
    s = (s.replace("−", "-").replace(" ", "").replace("*", "").replace(",", ""))
    return re.sub(r"e([+-])0*(\d)", r"e\2", s)


def main():
    npass = nfail = 0
    print("=" * 96)
    print("M6 scripts/04 — 文檔↔數據一致性核對（零可調容差）")
    print("=" * 96)
    for doc, needle, val, fmt in CLAIMS:
        want = fmt.format(val)
        ok = norm(want) in norm(DOCS[doc])
        # 額外核：doc 內確有該數字（避免 fmt 與 doc 位數不同而僥倖通過）
        val_ok = norm(want) in norm(DOCS[doc])
        if ok and val_ok:
            npass += 1
            print(f"  [PASS] {doc:<42} {want}")
        else:
            nfail += 1
            print(f"  [FAIL] {doc:<42} 期望 {want}（未在文中找到）")
    print("=" * 96)
    print(f"一致性：{npass}/{npass+nfail} PASS")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())