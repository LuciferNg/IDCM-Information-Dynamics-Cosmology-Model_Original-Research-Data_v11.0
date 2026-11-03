#!/usr/bin/env bash
# M11 run_chain.sh — 全鏈驅動（idcm-audit-verify §3.1 ⑥）
# 用法：
#   bash scripts/run_chain.sh            # dry-run（只印將執行嘅命令）
#   bash scripts/run_chain.sh --fast     # 略過 06（κ 掃描，~7 分鐘）與 11（審計自身）
#   bash scripts/run_chain.sh --full     # 全部（含 06）；L2→L3→audit 依序，任一 FAIL 即 exit 1
# ★ BLAS 必須鎖單線程（WSL 互搶會令 numpy 解算 timeout）
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
PY=/usr/bin/python3
cd "$(dirname "$0")/.." || exit 2
MODE="${1:---dry}"

run() {           # run <log> <script>
  local log="$1" scr="$2"
  if [ "$MODE" = "--dry" ]; then
    echo "[dry] $PY scripts/$scr > data/$log"
    return 0
  fi
  echo "── $scr"
  "$PY" "scripts/$scr" > "data/$log" 2>&1 || { echo "✗ FAIL: $scr（見 data/$log）"; return 1; }
}

rounds() {
  run run_01_r1_criteria.log            01_r1_criteria.py
  run run_02_r2_inventory_scan.log      02_r2_inventory_scan.py
  run run_03_r3_concave.log             03_r3_concave_identifiability.py
  run run_04_r4_necessity.log           04_r4_necessity.py
  run run_05_r5_k0_identity.log         05_r5_k0_identity.py
  if [ "$MODE" = "--full" ]; then
    run run_06_r6b_kappa_sweep.log      06_r6b_kappa_sweep.py || return 1
  else
    echo "── 06_r6b_kappa_sweep.py（--fast 略過；--full 才跑）"
  fi
  run run_07_em_anchor_scan.log         07_em_anchor_scan.py
  run run_08_g0_x_operator.log          08_g0_x_operator.py
  run run_09_route_c_total_dm.log       09_route_c_total_dm.py
  run run_10_route_c_degeneracy.log     10_route_c_degeneracy.py
  run run_13_k0_branch_sensitivity.log  13_k0_branch_sensitivity.py
}

echo "== L1（$MODE）=="
rounds || exit 1

if [ "$MODE" != "--dry" ]; then
  echo "== L2（獨立重推導）=="
  "$PY" scripts/l2_independent_verify.py > data/run_l2.log 2>&1 || { echo "✗ L2 FAIL"; exit 1; }
  echo "== L3（文檔交叉核）=="
  "$PY" scripts/l3_doc_crosscheck.py > data/run_l3.log 2>&1 || { echo "✗ L3 FAIL"; exit 1; }
  echo "== audit（全鏈）=="
  "$PY" scripts/11_audit_chain.py > data/run_11_audit_chain.log 2>&1 || { echo "✗ audit FAIL"; exit 1; }
  echo "✓ 全鏈綠燈（L1→L2→L3→audit）"
fi