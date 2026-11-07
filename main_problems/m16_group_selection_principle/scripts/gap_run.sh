#!/bin/bash
# scripts/gap_run.sh — 本地 GAP 啟動器（無 root；二進位放記錄內持久位 data/gap/root）
# 用法：bash scripts/gap_run.sh <gap-script>
# 前置：bash scripts/gap_setup.sh（一次性；見該檔版本釘）
# 註：本機 GAP 未綁 ARGV，故 .g 腳本一律把輸出路徑寫死，不靠命令列參數。
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REC="$(dirname "$HERE")"
GAPROOT="$REC/data/gap/root"
GAPBIN="$GAPROOT/usr/lib/x86_64-linux-gnu/gap/gap"
GAPLIB="$GAPROOT/usr/share/gap"
if [ ! -x "$GAPBIN" ]; then
  echo "[gap_run] GAP 未部署：請先跑 bash scripts/gap_setup.sh" >&2
  exit 127
fi
export LD_LIBRARY_PATH="$GAPROOT/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd "$REC"
exec "$GAPBIN" -l "$GAPLIB" -q -b "$@"