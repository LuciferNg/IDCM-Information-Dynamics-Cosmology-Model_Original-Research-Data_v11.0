#!/usr/bin/env bash
# M8 推導鏈 driver（預設 dry-run；加 --run 才真跑）
# 用法:
#   bash scripts/run_chain.sh            # dry-run：只印將執行之命令
#   bash scripts/run_chain.sh --run      # 真跑全部，輸出寫 data/run_NN.log
#   bash scripts/run_chain.sh --run 05   # 只跑到 05 為止（前綴亦可）
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=/usr/bin/python3
export PYTHONDONTWRITEBYTECODE=1

MODE="${1:-dry}"
STOP="${2:-}"

# 順序即依賴序：L1（01–08）→ L2 獨立重推（10）→ L3 文檔交叉核（11）→ 全審計（09）
CHAIN=(
  01_fan_geometry
  02_two_chart_junction
  03_m0_absolute_usage_audit
  04_two_factor_decomposition
  05_hierarchy_reduction
  06_three_input_falsifiability
  07_target_shift_discriminating_power
  08_forcing_candidates
  10_independent_verifier
  11_doc_crosscheck
  09_full_audit
)

fail=0
for s in "${CHAIN[@]}"; do
  num="${s:0:2}"
  if [ "$MODE" = "--run" ]; then
    "$PY" "scripts/${s}.py" > "data/run_${num}.log" 2>&1
    rc=$?
    p=$(grep -c '^\[PASS\]' "data/run_${num}.log" || true)
    f=$(grep -c '^\[FAIL\]' "data/run_${num}.log" || true)
    printf '%-34s exit=%d  PASS=%s  FAIL=%s\n' "$s" "$rc" "$p" "$f"
    [ "$rc" -ne 0 ] && fail=1
  else
    printf '[dry] %s scripts/%s.py > data/run_%s.log\n' "$PY" "$s" "$num"
  fi
  [ -n "$STOP" ] && { [[ "$s" == "$STOP"* ]] && break; }
done

if [ "$MODE" != "--run" ]; then
  echo
  echo "（dry-run：未執行。加 --run 真跑。）"
fi
exit "$fail"