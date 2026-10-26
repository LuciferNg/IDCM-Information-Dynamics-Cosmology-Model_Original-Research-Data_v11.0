#!/usr/bin/env bash
# M7 推導鏈 driver（idcm-audit-verify §3.1 ⑥）
# 用法：scripts/run_chain.sh          # dry-run（預設，只印命令）
#       scripts/run_chain.sh --run    # 真跑，逐 script 覆寫 data/run_<name>.log
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
DRY=1; [ "${1:-}" = "--run" ] && DRY=0
CHAIN=(
  01_reproduction_gate
  02_r5_mode_content_nogo
  02b_r5_mode_content_nogo
  03_grep_evidence
  05_additivity_scope
  06_gap_closure
  07_audit_full_chain
  04_doc_consistency   # 必須最後：核對全部 producer 之 data/*.json
)
fail=0
for s in "${CHAIN[@]}"; do
  cmd="PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 scripts/${s}.py > data/run_${s}.log 2>&1"
  if [ "$DRY" -eq 1 ]; then
    echo "[dry] $cmd"
  else
    printf '[run] %-32s' "$s"
    if eval "$cmd"; then echo " ok"; else echo " FAIL"; fail=1; break; fi
  fi
done
exit "$fail"
