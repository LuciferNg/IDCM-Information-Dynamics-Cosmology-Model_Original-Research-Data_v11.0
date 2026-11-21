#!/usr/bin/env bash
# M19 / scripts/run_chain.sh — 全鏈 driver（scripts → L1 → L2 → L3 → L2b/coverage → audit）
# 用法：bash scripts/run_chain.sh [--dry-run]
# 任何一支 FAIL ⟹ 傳播 exit 1。
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
PY=/usr/bin/python3
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 1
DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1
rc_all=0

run () {                       # run <log> <script>
  local log="$1"; local script="$2"
  if [ "$DRY" = "1" ]; then echo "[dry] $script -> $log"; return 0; fi
  echo "── run $script"
  if ! $PY "$script" > "$log" 2>&1; then
    echo "   FAIL rc≠0 : $script  (see $log)"; rc_all=1; return 1
  fi
  grep -E '^RESULT ' "$log" | tail -1 || { echo "   FAIL 無 RESULT 行 : $script"; rc_all=1; return 1; }
}

# ── L1／核心鏈（順序固定）
run data/run_00_reproduction_gate.log        scripts/00_reproduction_gate.py
run data/run_01_frame_census.log             scripts/01_frame_census.py
run data/run_02_bit_decomposition.log        scripts/02_bit_decomposition.py
run data/run_03_admissibility_and_readings.log scripts/03_admissibility_and_readings.py
run data/run_04_cbt_scan_density.log         scripts/04_cbt_scan_density.py
run data/run_05_registry_build.log           scripts/05_registry_build.py
run data/run_07_slot_scope_theorem.log       scripts/07_slot_scope_theorem.py
run data/run_09_naming_operator_class.log    scripts/09_naming_operator_class.py
run data/run_10_assignment_prior.log         scripts/10_assignment_prior.py

# ── L2／L3／coverage／audit（依賴上面產物）
run data/run_08_reverse_numeric_coverage.log scripts/08_reverse_numeric_coverage.py
run data/run_l1_reproducibility.log          scripts/l1_reproducibility.py
run data/run_l2_independent_recompute.log    scripts/l2_independent_recompute.py
run data/run_l3_doc_crosscheck.log           scripts/l3_doc_crosscheck.py
run data/run_06_audit_verify.log             scripts/06_audit_verify.py

echo "── chain rc=$rc_all"
exit $rc_all
