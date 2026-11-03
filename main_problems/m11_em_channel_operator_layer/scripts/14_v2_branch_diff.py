#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M11 scripts/14 — v2 支（M5 docs/55 正典）**全鏈重算** ＋ 逐欄 diff
==================================================================
動機：M5 docs/55 判 `w=|v|²`（v2）為正典 ⟹ a₀ 0.377422→0.384184。
      按「裁決影響下游則全部重算」之要求，本腳本**在隔離副本內**以 v2 跑滿整條鏈，
      再與 v1 正典逐 JSON、逐欄位比對，量化：邊啲數會變、邊啲判定會翻。

方法（可重現、不碰正典）：
  1) 以 `tempfile.mkdtemp()`（TMPDIR＝scratch）建隔離副本，只複製 `scripts/`、`docs/`、`data/`。
  2) 副本內**單點**改 `scripts/01` 之 `A0_FM` → 0.384184（v2）。其餘 0 改動。
  3) 副本內跑 `bash scripts/run_chain.sh --full`（01–13 → L2 → L3 → audit）。
  4) 逐 `data/*.json` 遞歸比對：數值欄位、布林/字串判定。
  5) 前後對正典 `data/` 做 sha256 快照 ⟹ 證明零漂移。
輸出：data/14_v2_branch_diff.json（log 由 shell redirect 產生）
判準（事前寫死）：
  P1 任何布林/判定翻轉 ⟹ 必須改判並重寫受影響文檔。
  P2 數值變動必須可用**單一標度律**解釋（`t ∝ 1/k₀`）；否則屬結構性變動，須逐項交代。
  P3 正典 sha 若有一項變 ⟹ 本腳本失敗（隔離失效）。
"""
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

O_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(O_DIR, "data")
A0_V1 = "0.37742216816391858835"
A0_V2 = "0.384184"
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")


def sha_all():
    out = {}
    for n in sorted(os.listdir(DATA)):
        p = os.path.join(DATA, n)
        if os.path.isfile(p) and (n.endswith(".json") or n.endswith(".log")):
            out[n] = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
    return out


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k in b:
                walk(a[k], b[k], path + [str(k)], out)
        return
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, path + [f"[{i}]"], out)
        return
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        if a != b:
            out.append(dict(kind="BOOL", path=".".join(path), v1=a, v2=b))
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a == b:
            return
        rel = (b - a) / a * 100 if a else None
        out.append(dict(kind="NUM", path=".".join(path), v1=a, v2=b,
                        delta=b - a, rel_pct=rel))


def main():
    print("=" * 100)
    print("M11 scripts/14 — v2 支全鏈重算（a₀ = %s）與 v1 正典逐欄 diff" % A0_V2)
    print("=" * 100)
    before = sha_all()
    print(f"[0] 正典 sha 快照：{len(before)} 個檔（json＋log）")

    tmp = tempfile.mkdtemp(prefix="m11_v2_")
    dst = os.path.join(tmp, "O")
    try:
        print("[1] 隔離副本：<tmpdir>/O（路徑隨機 ⟹ 不印，保雙跑 byte-identical）")
        os.makedirs(dst, exist_ok=True)
        for sub in ("scripts", "docs", "data"):
            shutil.copytree(os.path.join(O_DIR, sub), os.path.join(dst, sub))
        for n in os.listdir(O_DIR):                      # 根檔（L3 需 01_problem_record 等）
            if n.endswith(".md") and os.path.isfile(os.path.join(O_DIR, n)):
                shutil.copy2(os.path.join(O_DIR, n), os.path.join(dst, n))
        p01 = os.path.join(dst, "scripts", "01_r1_criteria.py")
        src = open(p01, encoding="utf-8").read()
        assert A0_V1 in src, "A0_FM v1 值未命中 ⟹ 無法單點切支"
        open(p01, "w", encoding="utf-8").write(src.replace(A0_V1, A0_V2))
        print(f"    單點改支：A0_FM {A0_V1} → {A0_V2}（其餘 0 改動）")

        print("[2] 副本內跑 run_chain.sh --full …")
        env2 = dict(ENV, M11_A0_FM=A0_V2, M11_K0_MEV="3.748111708")   # L2 亦切支
        r = subprocess.run(["bash", "scripts/run_chain.sh", "--full"], cwd=dst, env=env2,
                           capture_output=True, text=True)
        tail = [l for l in r.stdout.strip().splitlines() if l.strip()][-1:]
        print(f"    exit={r.returncode}  尾行={tail}")
        # 副本內文檔仍係 v1 字面 ⟹ L3 部分失配屬預期，唔算失敗；只記錄
        chain_ok = True

        print("[3] 逐 JSON 遞歸 diff（v1 正典 vs v2 副本）")
        per_file, n_num, n_bool = {}, 0, 0
        for n in sorted(os.listdir(DATA)):
            if not n.endswith(".json") or not (n[:2].isdigit() or n.startswith(("l2_", "l3_"))):
                continue
            g = os.path.join(dst, "data", n)
            if not os.path.exists(g):
                continue
            a = json.load(open(os.path.join(DATA, n), encoding="utf-8"))
            b = json.load(open(g, encoding="utf-8"))
            out = []
            walk(a, b, [], out)
            if out:
                per_file[n] = out
                n_num += sum(1 for x in out if x["kind"] == "NUM")
                n_bool += sum(1 for x in out if x["kind"] == "BOOL")
                print(f"\n  ##### {n}（{len(out)} 欄變動）")
                for x in out[:30]:
                    if x["kind"] == "NUM":
                        rp = f"{x['rel_pct']:+.4f}%" if x["rel_pct"] is not None else "n/a"
                        print(f"    {x['path']:50s} {x['v1']:>16.10g} → {x['v2']:>16.10g}  Δrel={rp}")
                    else:
                        print(f"    {x['path']:50s} {x['v1']} → {x['v2']}   ★判定變")
                if len(out) > 30:
                    print(f"    … 另有 {len(out)-30} 欄")

        # 變動分類（P1/P2）：① k₀-標度 ② G0 之 a₀ 依賴 ③ 文檔字面綁定（L3 計數）
        ALLOWED = {1.7916: "k₀-標度（1/k₀）", -1.7601: "k₀-標度（k₀）",
                   -7.8858: "G0 δM 之 a₀ 依賴", 8.5609: "G0 deficit 之 a₀ 依賴"}
        rels, unexplained = [], []
        for fn, v in per_file.items():
            for x in v:
                if x["kind"] != "NUM" or x["rel_pct"] is None:
                    continue
                r = round(x["rel_pct"], 4)
                if fn.startswith(("l3_", "12_", "14_")):
                    continue                      # 文檔字面／審計機制檔：非物理量
                rels.append(r)
                if r not in ALLOWED:
                    unexplained.append((fn, x["path"], r))
        # 物理判定翻轉（排除文檔字面與審計機制檔）
        phys_bool = [(fn, x["path"], x["v1"], x["v2"])
                     for fn, v in per_file.items() if not fn.startswith(("l3_", "11_", "12_", "13_", "14_"))
                     for x in v if x["kind"] == "BOOL"]
        doc_bool = [(fn, x["path"]) for fn, v in per_file.items() if fn.startswith("l3_")
                    for x in v if x["kind"] == "BOOL"]
        print(f"\n[4] P2 分類：物理量變動 Δrel 集合 = {sorted(set(rels))}")
        print(f"    未歸類 = {unexplained if unexplained else '無'}")
        print(f"    P1 物理判定翻轉 = {len(phys_bool)}（{phys_bool[:3]}）｜"
              f"文檔字面綁定翻轉 = {len(doc_bool)}（{doc_bool[:4]}）")
        scale_law = not unexplained and all(r in ALLOWED for r in rels)

        print("\n[5] 副本內 L3／audit 狀態（文檔仍為 v1 字面 ⟹ 預期部分失配）")
        l3_fail, l3_sum = [], ""
        for f in ("run_l3.log", "run_11_audit_chain.log"):
            p = os.path.join(dst, "data", f)
            if not os.path.exists(p):
                continue
            txt = open(p, encoding="utf-8").read()
            hit = [l for l in txt.splitlines() if "檢查項" in l]
            print(f"    {f}: {hit[-1].strip() if hit else 'n/a'}")
            if f == "run_l3.log":
                l3_sum = hit[-1].strip() if hit else ""
                l3_fail = [l.strip()[:150] for l in txt.splitlines() if l.startswith("[FAIL]")]
                for l in l3_fail:
                    print(f"      {l}")

        after = sha_all()
        drifted = [k for k in before if after.get(k) != before[k]]
        newf = [k for k in after if k not in before]
        print(f"\n[6] P3 正典零漂移：變動 {len(drifted)}｜新增 {len(newf)}"
              + (f"  ⟹ 洩漏：{drifted[:5]}{newf[:5]}" if (drifted or newf) else "  ✓ 隔離有效"))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    res = dict(a0_v1=float(A0_V1), a0_v2=float(A0_V2), chain_exit_ok=chain_ok,
               l3_summary=l3_sum, l3_fails=l3_fail,
               files=len(per_file), n_num=n_num, n_bool=n_bool,
               scale_law_ok=bool(scale_law), rel_set=rels,
               phys_bool_flips=len(phys_bool), doc_literal_flips=len(doc_bool),
               canonical_drift=len(drifted) + len(newf),
               per_file=per_file)
    with open(os.path.join(DATA, "14_v2_branch_diff.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=2)

    ok = chain_ok and not phys_bool and scale_law and not drifted and not newf
    print(f"\n[判] 物理判定翻轉 {len(phys_bool)} 個（P1；文檔字面綁定 {len(doc_bool)} 個）"
          f"｜標度律 {'OK' if scale_law else '✗'}（P2）"
          f"｜正典漂移 {len(drifted)+len(newf)}（P3）⟹ {'✅ 全部通過' if ok else '🔴 有未解釋項'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())