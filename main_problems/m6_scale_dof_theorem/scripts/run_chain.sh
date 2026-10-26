#!/usr/bin/env bash
# M6 全鏈重跑 driver（預設 dry-run；加 --run 才真跑）
# 用法：  scripts/run_chain.sh          # 只列將要跑的命令
#         scripts/run_chain.sh --run    # 全部重跑，輸出寫入 data/run_NN.log
set -u
cd "$(dirname "$0")/.." || exit 1
DO=0
[ "${1:-}" = "--run" ] && DO=1
for s in scripts/0[1-9]_*.py scripts/1[0-9]_*.py; do
  base=${s##*/}
  num=${base%%_*}
  echo "run data/run_${num}.log ${s}"
  if [ "$DO" -eq 1 ]; then
    PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 "$s" > "data/run_${num}.log" 2>&1
    rc=$?
    printf '  -> exit=%d  %s\n' "$rc" "$(grep -oE '[0-9]+/[0-9]+ PASS' "data/run_${num}.log" | tail -1)"
  fi
done