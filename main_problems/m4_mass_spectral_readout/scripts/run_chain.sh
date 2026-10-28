#!/bin/bash
# M4 audit-verify: re-run the whole L1 chain twice, capture logs + timings + exit codes.
# Canonical: OMP_NUM_THREADS=1 (BLAS thread thrash in WSL otherwise).
set -u
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
OUT=data/audit_20261028
mkdir -p "$OUT"

# ★ 動態枚舉（2026-10-28 審計修正）：原本 ORDER 係硬編碼白名單，止於 40 ⟹ 新增之 42/43
#   永遠唔會入 L1 全鏈 ⟹ 「L1 全鏈重跑」實為 40/43（覆蓋盲區）。改為自動枚舉，
#   只排除審計自身（41_chain_audit.py，否則自我重跑覆蓋 exit codes）。
EXCLUDE="41_chain_audit"
# ★ 排序規則（2026-10-28 第二修正）：`sort -V` 會把 `36b_…` 排在 `36_…` **之前**，
#   令已作廢之 v1 覆蓋 v2 之輸出（實測：data/36_em_kaon_external.json 被 v1 重寫）。
#   改用 (num, suffix) 排序：`36` < `36b` ⟹ v2 永遠最後寫入。
ORDER=$(cd scripts && /usr/bin/python3 -c "
import glob, os, re
fs = [f for f in glob.glob('[0-9][0-9]*.py') if not f.startswith('${EXCLUDE}')]
def key(f):
    m = re.match(r'(\d{2})([a-z]?)_', f)
    return (int(m.group(1)), m.group(2), f)
print(' '.join(os.path.splitext(f)[0] for f in sorted(fs, key=key)))
")

echo "chain start $(date -Iseconds)"
: > "$OUT/exit_codes.txt"
: > "$OUT/timings.txt"
for s in $ORDER; do
  f="scripts/${s}.py"
  [ -f "$f" ] || { echo "MISSING $f" | tee -a "$OUT/exit_codes.txt"; continue; }
  t0=$(date +%s.%N)
  /usr/bin/python3 "$f" > "$OUT/run_${s}.log" 2>&1
  rc=$?
  t1=$(date +%s.%N)
  dt=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}')   # bc 未必安裝，改用 awk
  printf '%-42s rc=%d  %ss\n' "$s" "$rc" "$dt" | tee -a "$OUT/exit_codes.txt"
  echo "$s $dt" >> "$OUT/timings.txt"
done
echo "chain end $(date -Iseconds)"
echo "--- nonzero exits ---"
grep -v 'rc=0' "$OUT/exit_codes.txt" || echo "(none)"