#!/bin/bash
# scripts/gap_setup.sh — 把 GAP 4.14.0 部署到本記錄內（持久位），供 00/02/02c/04 重跑
#
# 背景：舊版 gap_run.sh 指向 /home/wsl/.hermes/cache/scratch/gapdebs/root
#       ⟹ 該目錄會被 24h prune 清走 ⟹ 重跑變 FileNotFoundError（審計 §3.1③ 路徑衛生違規）。
# 修法：安裝到 <record>/data/gap/root（記錄內持久位，git 忽略二進位，只 commit 本腳本之版本釘）。
#
# 用法：bash scripts/gap_setup.sh            # 冪等：已存在則跳過
#       bash scripts/gap_setup.sh --force    # 重裝
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REC="$(dirname "$HERE")"
DEST="$REC/data/gap"
ROOT="$DEST/root"

# ── 版本釘（apt-cache policy 實測；改版須同步本檔與 docs/09 之審計記錄）
PKGS="gap gap-core gap-libs libgap9 gap-smallgrp gap-gapdoc gap-character-tables gap-atlasrep gap-primgrp gap-transgrp"
PIN="4.14.0-3"          # gap 主體
PIN_SG="1.5.4-2"        # SmallGrp

if [ "${1:-}" != "--force" ] && [ -x "$ROOT/usr/lib/x86_64-linux-gnu/gap/gap" ]; then
  echo "[gap_setup] already present: $ROOT"
  exit 0
fi

mkdir -p "$DEST/debs"
cd "$DEST/debs"
echo "[gap_setup] downloading (gap $PIN / smallgrp $PIN_SG) ..."
apt-get download $PKGS 2>&1 | tail -3

echo "[gap_setup] unpacking ..."
rm -rf "$ROOT"
mkdir -p "$ROOT"
for d in *.deb; do dpkg-deb -x "$d" "$ROOT"; done

GAPBIN="$ROOT/usr/lib/x86_64-linux-gnu/gap/gap"
GAPLIB="$ROOT/usr/share/gap"
[ -x "$GAPBIN" ] || { echo "[gap_setup] FAIL: no binary at $GAPBIN"; exit 1; }
[ -d "$GAPLIB" ] || { echo "[gap_setup] FAIL: no lib at $GAPLIB"; exit 1; }
export LD_LIBRARY_PATH="$ROOT/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
echo "[gap_setup] ok: $GAPBIN"
"$GAPBIN" -l "$GAPLIB" -q -b -c 'Print("GAP ", GAPInfo.Version, " | SmallGrp ", PackageInfo("smallgrp")[1].Version, "\n"); QUIT;'