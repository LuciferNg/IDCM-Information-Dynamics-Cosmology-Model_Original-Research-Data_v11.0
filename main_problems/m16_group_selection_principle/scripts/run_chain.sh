#!/bin/bash
# scripts/run_chain.sh — M16 全鏈驅動：L1 → L2 → L3 → audit
# 用法：bash scripts/run_chain.sh            # 全跑
#       bash scripts/run_chain.sh --dry-run  # 只列步驟
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

step() { echo "── $* ──"; }

run() { # run <log> <cmd...>
  local log="$1"; shift
  step "$* → $log"
  if [ "$DRY" = "1" ]; then return 0; fi
  ( "$@" ) > "$log" 2>&1 || { echo "FAIL: $*（見 $log）"; exit 1; }
}

# 全量掃描（重；~35 min）：scripts/02_gap_full_scan.g 經 11 分片驅動
# 只在 bash scripts/run_chain.sh --heavy 時跑；平時以 L1 之 md5 ≡ 入庫 作等價證明。
if [ "${1:-}" = "--heavy" ]; then
  run data/run_02_driver.log bash scripts/run_02_parallel.sh 11
fi

# 前置：GAP 部署（冪等）
run data/run_gap_setup.log bash scripts/gap_setup.sh

# L1 重現性（雙跑）
run data/run_l1_boot.log /usr/bin/python3 scripts/11_l1_reproducibility.py

# L2 獨立重推導
run data/run_l2_boot.log /usr/bin/python3 scripts/l2_independent_verify.py

# L3 文檔交叉核
run data/run_l3_boot.log /usr/bin/python3 scripts/l3_doc_crosscheck.py

# 審計
run data/run_audit_boot.log /usr/bin/python3 scripts/10_audit_verify.py

if [ "$DRY" = "0" ]; then
  echo "── 摘要 ──"
  /usr/bin/python3 - <<'EOF'
import json, pathlib
d = pathlib.Path("data")
l1 = json.loads((d/"run_l1.json").read_text())
l2 = json.loads((d/"l2_independent_verify.json").read_text())["_summary"]
l3 = json.loads((d/"l3_doc_crosscheck.json").read_text())["_summary"]
au = json.loads((d/"audit_verify.json").read_text())["_summary"]
print(f"L1 {l1['pass']}/{l1['total']} ｜ L2 {l2['checks']-l2['fail']}/{l2['checks']} ｜ "
      f"L3 {l3['checks']-l3['fail']}/{l3['checks']} ｜ audit {au['checks']-au['fail']}/{au['checks']}")
EOF
fi
echo "done (dry=$DRY)"