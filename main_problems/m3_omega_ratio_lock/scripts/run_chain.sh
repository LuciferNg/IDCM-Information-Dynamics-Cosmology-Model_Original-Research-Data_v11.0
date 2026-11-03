#!/usr/bin/env bash
# M3 全鏈 driver：L1（17 支）→ L2 → L3 → audit 18
# 用法：bash scripts/run_chain.sh [--run]   （預設 dry-run）
set -u
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY=/usr/bin/python3
DO=${1:-}

run() {  # run <outfile> <script>
  if [ "$DO" = "--run" ]; then
    echo "== $2 -> data/$1"
    $PY "scripts/$2" > "data/$1" 2>&1 || { echo "FAIL $2"; exit 1; }
  else
    echo "run $1 <- $2"
  fi
}

run 01_round1.out 01_round1_spectral_identities.py
run 02_round2.out 02_round2_g3_assignment.py
run 03_round3.out 03_round3_r3_operatorization.py
run 04_round4.out 04_round4_oplus_legitimacy.py
run 05_round5.out 05_round5_chain_length_origin.py
run 06_round6.out 06_round6_f2_onedim.py
run 07_round7.out 07_round7_single_component.py
run 08_round8.out 08_round8_correction_k.py
run 09_round9.out 09_round9_gamma1_identity.py
run 10_round10.out 10_round10_omega_b_fork.py
run 11_round11.out 11_round11_g3_fork.py
run 12_round12.out 12_round12_g3a_p4.py
run 13_round13.out 13_round13_index_audit.py
run 14_round14.out 14_round14_r5_mc_null.py
run 15_round15.out 15_round15_r1_gate_identity.py
run 16_round16.out 16_round16_a3_scale_invariance.py
run 17_round17.out 17_round17_pythagoras_audit.py

if [ "$DO" = "--run" ]; then
  echo "== L2"; $PY scripts/l2_independent_verify.py > data/run_l2_independent_verify.log 2>&1 || exit 1
  echo "== L3"; $PY scripts/l3_doc_crosscheck.py > data/run_l3_doc_crosscheck.log 2>&1 || exit 1
  echo "== audit"; $PY scripts/18_audit_verify.py > data/run_18_audit_verify.log 2>&1 || exit 1
  tail -2 data/run_18_audit_verify.log
else
  echo "run l2 <- l2_independent_verify.py"
  echo "run l3 <- l3_doc_crosscheck.py"
  echo "run audit <- 18_audit_verify.py"
fi