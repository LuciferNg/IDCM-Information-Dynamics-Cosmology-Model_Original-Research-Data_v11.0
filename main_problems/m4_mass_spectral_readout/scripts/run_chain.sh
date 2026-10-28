#!/bin/bash
# M4 audit-verify: re-run the whole L1 chain twice, capture logs + timings + exit codes.
# Canonical: OMP_NUM_THREADS=1 (BLAS thread thrash in WSL otherwise).
set -u
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
OUT=data/audit_20261028
mkdir -p "$OUT"

ORDER="01_quark_scan 02_mc_mcs_judge 03_reading_strength 04_lepton_anchor 05_hadron_dm 06_hop_probe 07_hit_followup 08_span_12_judge 09_density_gate 10_gamma1_gamma6_bridge 11_lattice_flag_span 12_A_spectral_identity 13_mdmu_second_instance 14_double_probe 15_exponent8 16_scale_operator 17_em_decomposition 18_em_constraints 19_delivery_processing 20_identifiability 21_upstream_verification 22_o216_o86_verification 23_full_audit 24_audit_round2 25_import_verification 26_em_consistency_window 27_em_sign_theorem 28_matrix_element_resolution 29_operator_structure 30_prediction_envelope 31_flavour_insertion 31b_flavour_insertion_enum 33_handoff_intake_20261026 34_u_d_role_citation_verify 35_v9_table_audit 36_em_kaon_external_judgement 36b_em_kaon_external_judgement 38_g3_decay_constant_precheck 39_m_series_handoff_intake 40_recheck_m7_m10"

echo "chain start $(date -Iseconds)"
: > "$OUT/exit_codes.txt"
: > "$OUT/timings.txt"
for s in $ORDER; do
  f="scripts/${s}.py"
  [ -f "$f" ] || { echo "MISSING $f" | tee -a "$OUT/exit_codes.txt"; continue; }
  t0=$(date +%s.%N)
  /usr/bin/python3 "$f" > "$OUT/run_${s}.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  printf '%-42s rc=%d  %.1fs\n' "$s" "$rc" "$(echo "$t1 - $t0" | bc)" | tee -a "$OUT/exit_codes.txt"
  echo "$s $(echo "$t1 - $t0" | bc)" >> "$OUT/timings.txt"
done
echo "chain end $(date -Iseconds)"
echo "--- nonzero exits ---"
grep -v 'rc=0' "$OUT/exit_codes.txt" || echo "(none)"