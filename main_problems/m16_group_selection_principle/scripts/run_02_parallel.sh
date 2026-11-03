#!/bin/bash
# scripts/run_02_parallel.sh — 並行掃描驅動（10 片）＋ 合併
set -u
cd "$(dirname "$0")/.." || exit 1
NP=${1:-10}
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 scripts/02b_chunk_driver.py "$NP" || exit 1
rm -f data/run_02_parallel.log
pids=()
for f in data/chunks/scan_*.g; do
  j=$(basename "$f" .g)
  bash scripts/gap_run.sh "$f" > "data/${j}.log" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} GAP workers" | tee -a data/run_02_parallel.log
fail=0
for p in "${pids[@]}"; do
  wait "$p" || fail=1
done
echo "workers done fail=$fail" | tee -a data/run_02_parallel.log
cat data/scan_*.log 2>/dev/null >> data/run_02_parallel.log
# 合併（每片輸出檔內已按 n 排序，此處再全域排序）
cat data/chunk_*.txt > data/run_02_scan_all_500.txt
sort -k1,1n -k2,2n -o data/run_02_scan_all_500.txt data/run_02_scan_all_500.txt
wc -l data/run_02_scan_all_500.txt | tee -a data/run_02_parallel.log
exit $fail