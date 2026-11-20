#!/usr/bin/env bash
# M22 driver — 全鏈重跑（依序；任何一支非 0 即停）
set -u
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
PY=/usr/bin/python3
fail=0
for s in scripts/00_reproduction_gate.py scripts/01_a5_invariants.py \
         scripts/02_amplitude_nogo.py scripts/03_comb_field_stats.py \
         scripts/04_null_and_sensitivity.py scripts/05_chirality.py \
         scripts/07_mr0_carrier.py scripts/08_mr1_partition.py \
         scripts/09_mr2_chirality_test.py scripts/10_mr3_lcdm_baseline.py \
         scripts/11_mr0_r7.py scripts/12_mr1_equivariance.py; do
  echo "=== $s ==="
  $PY "$s" > "data/$(basename "${s%.py}").log" 2>&1
  rc=$?
  tail -1 "data/$(basename "${s%.py}").log"
  echo "rc=$rc"
  [ $rc -ne 0 ] && fail=1 && break
done
echo "CHAIN_EXIT=$fail"
exit $fail
