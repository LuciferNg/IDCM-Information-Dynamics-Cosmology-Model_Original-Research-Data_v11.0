#!/bin/bash
# M15 全鏈 driver —— 依序跑 L1（01–05）→ L2 → L3 → audit（06）
# 用法：bash scripts/run_chain.sh
set -u
cd "$(dirname "$0")/.." || exit 1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/usr/bin/python3
RC=0

run() {  # run <script> <log>
  echo "== $1 -> data/$2"
  $PY "scripts/$1" > "data/$2" 2>&1 || { echo "   FAIL rc=$?"; RC=1; }
}

run 01_label_theorems.py            run_01_label_theorems.log
run 02_a5_harmonic_multiplicity.py  run_02_a5_harmonic.log
run 03_count_census.py              run_03_count_census.log
run 04_adic_map_nogo.py             run_04_adic_map.log
run 05_cy3_generation_count.py      run_05_cy3.log
run 07_h11_minus_M_candidate.py     run_07_h11_M.log
run 08_ngen_equals_B_minus_1.py     run_08_ngen_B1.log
run 09_tenc_c_downgrade_audit.py    run_09_downgrade.log
run 10_ngen_provenance_scan.py      run_10_provenance.log
run 11_gate_count_candidate.py      run_11_gate_count.log
run 12_gate_budget_identity.py      run_12_budget.log
run l2_independent_verify.py        run_l2_independent_verify.log
run l3_doc_crosscheck.py            run_l3_doc_crosscheck.log

echo "== 06_audit_verify.py（不入 L1；自行兩跑各 L1）"
$PY scripts/06_audit_verify.py > data/run_06_audit_verify.log 2>&1 || { echo "   FAIL rc=$?"; RC=1; }
$PY scripts/06_audit_verify.py > data/_audit_run2.log 2>&1
if cmp -s data/run_06_audit_verify.log data/_audit_run2.log; then
  echo "   審計自跑兩次 byte-identical"
else
  echo "   審計自跑兩次唔一致"; RC=1
fi
rm -f data/_audit_run2.log
tail -4 data/run_06_audit_verify.log

echo "== done rc=$RC"
exit $RC