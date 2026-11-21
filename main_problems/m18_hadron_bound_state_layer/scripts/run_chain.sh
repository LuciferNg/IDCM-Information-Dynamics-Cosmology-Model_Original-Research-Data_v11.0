#!/usr/bin/env bash
# M18 driver：全鏈 L1（01–18）→ 審計（08）→ 二層／三層審計（19）
# 用法：
#   ./scripts/run_chain.sh            # dry-run（列出將會跑嘅步驟）
#   ./scripts/run_chain.sh --run      # 真跑
#   ./scripts/run_chain.sh --run 2>/dev/null   # 靜默
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
cd "$(dirname "$0")/.." || exit 1
PY=/usr/bin/python3
MODE="${1:-}"

run() {   # run <log> <script>
  local log="$1" py="$2"
  if [ "$MODE" != "--run" ]; then
    echo "[dry] $PY scripts/$py > data/$log"
    return 0
  fi
  echo "[run] $py"
  $PY "scripts/$py" > "data/$log" 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "[FAIL] $py rc=$rc（見 data/$log）" >&2
    exit 1
  fi
  return 0
}

for n in 01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 21 22; do
  f=$(ls scripts/${n}_*.py | head -1 | sed 's|scripts/||')
  run "run_${f%.py}.log" "$f"
done

# 審計最後跑（次序固定：改腳本 → 重跑 → 才跑 audit）
run run_20_open_problems_census.log 20_open_problems_census.py
run run_08_full_audit.log      08_full_audit.py
run run_19_chain_audit.log      19_chain_audit.py

echo "[done] chain $([ "$MODE" = "--run" ] && echo executed || echo dry-run)"
