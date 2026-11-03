#!/bin/bash
# scripts/gap_run.sh — 本地 GAP 啟動器（無 root：deb 解包至 scratch，僅供本 session 使用）
# 用法：bash scripts/gap_run.sh <gap-script> [args...]
GAPROOT=/home/wsl/.hermes/cache/scratch/gapdebs/root
export LD_LIBRARY_PATH="$GAPROOT/usr/lib/x86_64-linux-gnu"
exec "$GAPROOT/usr/lib/x86_64-linux-gnu/gap/gap" -l "$GAPROOT/usr/share/gap" -q -b "$@"