#!/usr/bin/env bash
# M1 推導鏈 driver
#   用法： ./scripts/run_chain.sh            # dry-run（只列步驟）
#          ./scripts/run_chain.sh --run      # 逐支執行，寫 data/run_<name>.log
set -u
cd "$(dirname "$0")/.." || exit 1
export PYTHONDONTWRITEBYTECODE=1
PY=/usr/bin/python3
STEPS=(
  "00_layer_perturbation.py"
  "01_multiplicity_engine.py"
  "02_mc_null.py"
  "03_canonical_form.py"
  "04_three_chain_presets.py"
  "06_multiquantum.py"
  "07_downstream_propagation.py"
  "05_doc_consistency.py"
)
RUN=0
[ "${1:-}" = "--run" ] && RUN=1
mkdir -p data
for s in "${STEPS[@]}"; do
  if [ "$RUN" -eq 0 ]; then
    echo "[dry-run] $PY scripts/$s"
  else
    log="data/run_${s%.py}.log"
    echo "[run] scripts/$s -> $log"
    "$PY" "scripts/$s" > "$log" 2>&1
    rc=$?
    echo "      exit=$rc"
    [ $rc -ne 0 ] && { echo "CHAIN STOPPED at $s"; exit $rc; }
  fi
done
[ "$RUN" -eq 0 ] && echo "[dry-run] 完（用 --run 真跑）"
exit 0