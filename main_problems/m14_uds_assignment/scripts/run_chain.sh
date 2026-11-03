#!/usr/bin/env bash
# M14 全鏈 driver（idcm-audit-verify §3.1 ⑥）
#   bash scripts/run_chain.sh          # dry-run
#   bash scripts/run_chain.sh --run    # 真跑 L1 → L2 → L3 → audit
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=/usr/bin/python3
MODE="${1:-dry}"

mapfile -t L1 < <(ls scripts/[0-9][0-9]_*.py 2>/dev/null | grep -v audit | sort)
mapfile -t L2 < <(ls scripts/l2_*.py 2>/dev/null | sort)
mapfile -t L3 < <(ls scripts/l3_*.py 2>/dev/null | sort)
mapfile -t AU < <(ls scripts/*audit*.py 2>/dev/null | sort)
CHAIN=("${L1[@]}" "${L2[@]}" "${L3[@]}" "${AU[@]}")
fail=0
for s in "${CHAIN[@]}"; do
  base=$(basename "$s" .py)
  log="data/run_${base}.log"
  if [ "$MODE" = "--run" ]; then
    "$PY" "$s" > "$log" 2>&1; rc=$?
    p=$(grep -c '^\[PASS\]' "$log" || true); f=$(grep -c '^\[FAIL\]' "$log" || true)
    printf '%-42s exit=%d PASS=%s FAIL=%s\n' "$base" "$rc" "${p:-0}" "${f:-0}"
    { [ "$rc" -ne 0 ] || [ "${f:-0}" -ne 0 ]; } && fail=1
  else
    printf '[dry] %s %s > %s\n' "$PY" "$s" "$log"
  fi
done
[ "$MODE" != "--run" ] && echo && echo "（dry-run：未執行。加 --run 真跑。）"
exit "$fail"
