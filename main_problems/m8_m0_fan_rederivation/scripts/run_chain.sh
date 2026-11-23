#!/usr/bin/env bash
# M8 推導鏈 driver（預設 dry-run；加 --run 才真跑）
# 用法:
#   bash scripts/run_chain.sh            # dry-run：只印將執行之命令
#   bash scripts/run_chain.sh --run      # 真跑全部，輸出寫 data/run_NN.log
#   bash scripts/run_chain.sh --run 05   # 只跑到 05 為止（前綴亦可）
#
# 次序（**動態掃描**，2026-11-24 改：唔再硬編碼清單，新增腳本自動入鏈）：
#   ① 數字命名 `NN_*.py`，名稱含 audit 者剔除（＝L1，按號碼排序）
#   ② `l2_*.py`（獨立重推）
#   ③ `l3_*.py`（文檔交叉核）
#   ④ 名稱含 audit 者殿後（會重跑全部）
# 任一 FAIL>0 或 exit≠0 ⟹ 整體 exit 1。
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PY=/usr/bin/python3
export PYTHONDONTWRITEBYTECODE=1

MODE="${1:-dry}"
STOP="${2:-}"

mapfile -t AU < <(ls scripts/*full_audit*.py 2>/dev/null | sort)
L1=()
for f in $(ls scripts/[0-9][0-9]_*.py 2>/dev/null | sort); do
  b=$(basename "$f"); skip=0
  for a in ${AU[@]+"${AU[@]}"}; do
    [ "$(basename "$a")" = "$b" ] && { skip=1; break; }
  done
  [ "$skip" = 0 ] && L1+=("$f")
done
mapfile -t L2 < <(ls scripts/l2_*.py 2>/dev/null | sort)
mapfile -t L3 < <(ls scripts/l3_*.py 2>/dev/null | sort)
CHAIN=(${L1[@]+"${L1[@]}"} ${L2[@]+"${L2[@]}"} ${L3[@]+"${L3[@]}"} ${AU[@]+"${AU[@]}"})

fail=0
for s in "${CHAIN[@]}"; do
  base=$(basename "$s" .py)
  num="${base:0:2}"
  log="data/run_${num}.log"
  if [ "$MODE" = "--run" ]; then
    "$PY" "$s" > "$log" 2>&1
    rc=$?
    p=$(grep -c '^\[PASS\]' "$log" || true)
    f=$(grep -c '^\[FAIL\]' "$log" || true)
    printf '%-38s exit=%d  PASS=%s  FAIL=%s\n' "$base" "$rc" "$p" "$f"
    [ "$rc" -ne 0 ] && fail=1
    [ "${f:-0}" -ne 0 ] && fail=1
  else
    printf '[dry] %s %s > %s\n' "$PY" "$s" "$log"
  fi
  if [ -n "$STOP" ] && [[ "$base" == "$STOP"* ]]; then break; fi
done

if [ "$MODE" != "--run" ]; then
  echo
  echo "（dry-run：未執行。加 --run 真跑。）"
fi
exit "$fail"
