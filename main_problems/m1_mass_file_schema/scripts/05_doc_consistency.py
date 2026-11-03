#!/usr/bin/env python3
"""M1 scripts/05 — 文檔↔數據一致性檢查器

做法：由 data/*.json 取真值 → 抽出本 M docs/*.md 內全部數字 token → 逐個真值檢查
      「文檔內存在一個數字在給定有效位數容差內等於真值」。容差由文檔顯示位數推導（10^-(sig-1)），
      零可調參數；容差只放寬顯示圓整，唔准遮蓋真差異。
輸出：stdout（claims / pass / fail）；有 FAIL 則 exit 1。
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(ROOT, "data")
DOCS = os.path.join(ROOT, "docs")

NUM = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def load(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as fh:
        return json.load(fh)


def nums_of(doc_name):
    path = os.path.join(DOCS, doc_name)
    if not os.path.exists(path):
        path = os.path.join(ROOT, doc_name)
    with open(path, encoding="utf-8") as fh:
        text = fh.read().replace("−", "-").replace(",", "")
    return [float(t) for t in NUM.findall(text)]


def main():
    d00, d01, d02, d03, d04, d06, d07 = (load("00_layer_perturbation.json"), load("01_multiplicity.json"),
                                         load("02_mc_null.json"), load("03_canonical.json"),
                                         load("04_three_chain.json"), load("06_multiquantum.json"),
                                         load("07_downstream.json"))
    docs = {n: nums_of(n) for n in ("01_layer_reduction.md", "02_multiplicity.md",
                                    "03_canonical_form.md", "04_three_chain_presets.md",
                                    "05_prereg_results.md", "06_multiquantum.md", "08_downstream_impact.md", "README.md",
                                    "01_problem_record.md", "SESSION_START.md")}
    claims = []

    def C(doc, value, sig):
        claims.append((doc, float(value), sig))

    # ---- docs/01（R1）----
    C("01_layer_reduction.md", 4.2169e-81, 5)
    mu_row = next(r for r in d00["L1_vertices"]["N_perturbation"]["rows"] if r["fermion"] == "mu")
    for k in ("dev_N86", "dev_N88", "dev_N90"):
        C("01_layer_reduction.md", mu_row[k], 3)
    e_row = next(r for r in d00["L1_vertices"]["N_perturbation"]["rows"] if r["fermion"] == "e")
    for k in ("dev_N86", "dev_N88", "dev_N90"):
        C("01_layer_reduction.md", e_row[k], 3)
    C("01_layer_reduction.md", d00["check_summary"]["pass"], 2)

    # ---- docs/02（R2）----
    for f in ("e", "mu", "tau", "u", "d", "s", "c", "b", "t"):
        C("02_multiplicity.md", d01["rows"][f]["target_units"], 6)
    for k in ("1", "2", "3", "4", "8", "12"):
        C("02_multiplicity.md", d01["rows"]["mu"]["rho"]["1e-03"][k], 2)
        C("02_multiplicity.md", d01["rows"]["s"]["rho"]["1e-03"][k], 2)
    for k in ("2", "3", "4", "5"):
        C("02_multiplicity.md", d01["rows"]["mu"]["rho_1p"][k], 3)
        C("02_multiplicity.md", d01["rows"]["s"]["rho_1p"][k], 3)
        C("02_multiplicity.md", d01["rows"]["d"]["rho_1p"][k], 2)
    for key, sig in (("k1_e1e-03", 4), ("k2_e1e-03", 4), ("k3_e1e-03", 5), ("k4_e1e-03", 5)):
        C("02_multiplicity.md", d01["coverage"][key] * 100, sig)
    C("02_multiplicity.md", d02["null_A"]["mean_rho"], 4)
    C("02_multiplicity.md", d02["null_A"]["p_rho_zero"], 4)
    C("02_multiplicity.md", d02["null_B"]["p_rho_zero"], 4)
    for null in ("null_A", "null_B"):
        for f in ("mu", "s"):
            C("02_multiplicity.md", d02[null]["per_fermion"][f]["p_ge"], 2)
    C("02_multiplicity.md", d01["rows"]["mu"]["closest_dev_k4"], 4)
    C("02_multiplicity.md", d01["rows"]["s"]["closest_dev_k4"], 4)

    # ---- docs/03（R3）----
    r3 = d03["rows"]
    for f, sig in (("mu", 5), ("s", 5), ("tau", 5), ("c", 6)):
        C("03_canonical_form.md", r3[f]["dev"], sig)
        C("03_canonical_form.md", r3[f]["degeneracy"], 3)
        C("03_canonical_form.md", r3[f]["k"], 2)
    for f, sig in (("e", 10), ("u", 9), ("d", 9)):
        C("03_canonical_form.md", r3[f]["sigma44"], sig)
        C("03_canonical_form.md", r3[f]["degeneracy"], 2)
    for f in ("b", "t"):
        C("03_canonical_form.md", r3[f]["per_quantum_mean_units"], 8)
        C("03_canonical_form.md", r3[f]["n_q"], 2)

    # ---- docs/04（R4）----
    for f in ("e", "mu", "u", "d", "s"):
        C("04_three_chain_presets.md", d04["P1_O54"]["rows"][f]["n_feasible_modes"], 2)

    # ---- docs/05（對帳）----
    C("05_prereg_results.md", d02["null_A"]["per_fermion"]["mu"]["p_ge"], 2)
    C("05_prereg_results.md", d01["coverage"]["k4_e1e-03"] * 100, 5)
    C("05_prereg_results.md", d03["rows"]["mu"]["degeneracy"], 1)
    C("05_prereg_results.md", d03["rows"]["tau"]["degeneracy"], 3)
    C("05_prereg_results.md", d03["rows"]["c"]["degeneracy"], 4)
    C("05_prereg_results.md", 63.583465, 6)
    C("05_prereg_results.md", 66.146930, 6)
    C("02_multiplicity.md", 45.505650, 6)

    # ---- 三件套／README 之數字同步（按各檔實際顯示位數）----
    b_resid = d03["rows"]["b"]["target_units"]
    C("README.md", d03["rows"]["mu"]["degeneracy"], 1)
    C("README.md", d03["rows"]["s"]["degeneracy"], 1)
    C("README.md", d03["rows"]["tau"]["degeneracy"], 3)
    C("README.md", d03["rows"]["c"]["degeneracy"], 4)
    for f, sig in (("e", 10), ("u", 9), ("d", 9)):
        C("README.md", d03["rows"][f]["sigma44"], sig)
    C("README.md", d01["coverage"]["k4_e1e-03"] * 100, 4)
    C("README.md", d02["null_A"]["per_fermion"]["mu"]["p_ge"], 2)
    C("README.md", 149.815877839 - 90, 5)
    C("README.md", 69 * 90 - 6185.10192312, 5)
    for f, sig in (("mu", 3), ("s", 3), ("tau", 3), ("c", 4)):
        C("01_problem_record.md", d03["rows"][f]["degeneracy"], sig)
    for f, sig in (("e", 10), ("u", 9), ("d", 9)):
        C("01_problem_record.md", d03["rows"][f]["sigma44"], sig)
    for f in ("b", "t"):
        C("01_problem_record.md", d03["rows"][f]["per_quantum_mean_units"], 8)
    C("01_problem_record.md", d01["coverage"]["k4_e1e-03"] * 100, 5)
    C("01_problem_record.md", d02["null_A"]["per_fermion"]["mu"]["p_ge"], 2)
    for f, sig in (("mu", 3), ("s", 3), ("tau", 3), ("c", 4)):
        C("SESSION_START.md", d03["rows"][f]["degeneracy"], sig)
    C("SESSION_START.md", d01["coverage"]["k4_e1e-03"] * 100, 5)
    C("SESSION_START.md", d02["null_A"]["per_fermion"]["mu"]["p_ge"], 2)
    C("SESSION_START.md", 149.815877839 - 90, 7)
    C("SESSION_START.md", 69 * 90 - 6185.10192312, 7)

    # ---- docs/06（R6）----
    for name, rec in d06["targets"].items():
        C("06_multiquantum.md", rec["value_units"], 8)
        C("06_multiquantum.md", rec["k_min"], 2)
        C("06_multiquantum.md", rec["single_mode_best"]["dev"], 4)
        C("06_multiquantum.md", rec["rho_kmin"], 2)
    C("06_multiquantum.md", d06["coverage"]["k1"] * 100, 4)
    C("06_multiquantum.md", d06["coverage"]["k2"] * 100, 4)
    C("06_multiquantum.md", d06["coverage"]["k3"] * 100, 5)
    C("06_multiquantum.md", d06["R6"]["R6_Q2"]["per_quantum_deficit_b"], 8)
    C("06_multiquantum.md", d06["R6"]["R6_Q2"]["per_quantum_deficit_t"], 8)
    C("06_multiquantum.md", d06["R6"]["R6_Q2"]["ratio"], 6)

    # ---- docs/08（下游重推）----
    for k, rec in d07["P3_coverage"].items():
        C("08_downstream_impact.md", rec["uniform_t"] * 100, 3)
        C("08_downstream_impact.md", rec["loguniform_pool"] * 100, 3)
        C("08_downstream_impact.md", rec["loguniform_M4conv"] * 100, 3)
    C("08_downstream_impact.md", d07["P1_M0"]["ratio"], 6)
    C("08_downstream_impact.md", d07["P2_M6_nq"]["m_p"]["units"], 4)
    C("08_downstream_impact.md", d07["P2_M6_nq"]["v_SM"]["units"], 6)
    C("08_downstream_impact.md", d07["P2_M6_nq"]["M_Pl"]["units"], 6)
    C("08_downstream_impact.md", d07["P2_M6_nq"]["M_Pl"]["n_q_min"], 8)
    C("08_downstream_impact.md", d07["P4_support_gate"]["max_k_with_coverage_le_0.1"], 2)
    C("08_downstream_impact.md", d07["P6_lattice"]["n_hits"], 2)

    fails = []
    for doc, val, sig in claims:
        rtol = 10.0 ** (-(sig - 1))
        hit = any(abs(x - val) <= rtol * max(abs(val), 1e-300) for x in docs[doc])
        if not hit:
            fails.append({"doc": doc, "value": val, "sig": sig, "rtol": rtol})

    print(json.dumps({"claims": len(claims), "pass": len(claims) - len(fails), "fail": len(fails)},
                     ensure_ascii=False))
    for f in fails:
        print("FAIL", f)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()