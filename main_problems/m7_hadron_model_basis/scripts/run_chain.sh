#!/usr/bin/env bash
# M7 推導鏈 driver（idcm-audit-verify §3.1 ⑥）
# 用法：scripts/run_chain.sh          # dry-run（預設，只印命令）
#       scripts/run_chain.sh --run    # 真跑，逐 script 覆寫 data/run_<name>.log
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
DRY=1; [ "${1:-}" = "--run" ] && DRY=0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # 防 BLAS 線程互搶（02 由 >300s 降至 7.6s）
CHAIN=(
  01_reproduction_gate
  02_r5_mode_content_nogo
  02b_r5_mode_content_nogo
  03_grep_evidence
  05_additivity_scope
  06_gap_closure
  10_watermark          # 水位重生（07 C2 之輸入）
  08_handoff_matrix     # register 收訖矩陣（唯讀）
  09_m4_receipt         # M4 交付收據核實
  11_o120_import_verify # O120 docs/06 引入驗證（29/29）
  12_o209_assignment    # O209 唯一分配重建（13/13）
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
