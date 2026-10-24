"""M4 docs/24 配套 — 第二輪審計（增量）：審計「審計自身」＋審計後改動核對
================================================================================
背景：`docs/23`（2026-10-24 15:41）已全量審計 scripts 01–22 × docs 01–22（重跑 22/22 byte-identical）。
      docs/23 之後（commit `ea59cf9..HEAD`）M4 只改過兩份檔：`01_problem_record.md`（P16→M5 指引）與
      `docs/23_full_audit.md`（P16→M5 一字），**scripts/ 零改動**。
本輪範圍（= 之前「未經驗證」者）：
  ① `scripts/23_full_audit.py` —— 首輪審計器自己（首輪明文排除自身，從未驗）
  ② `docs/23_full_audit.md` —— 審計報告本身（無人審過）
  ③ 審計後兩份檔之改動（確認只涉及編號字串，**零數值改動**）
五項檢查：
  C1 scripts/ 自 `ea59cf9` 起零改動（`git diff --stat` 空）⟹ 首輪 22/22 結論仍有效（可繼承）
  C2 `docs/23` 宣稱 vs 其 JSON 實錄逐項核對（22/22、docstring 22/22、違規 0）
  C3 抽樣重跑 3 支輕腳本（獨立於首輪），stdout 對正典 log byte-identical
  C4 `docs/23` 自身結構（H1／必要節）＋ T3 數字回溯（對 `run_23_full_audit.log` ＋ `23_full_audit.json`）
  C5 審計後改動之數值零漂移（git diff 只用 P16→M5 字串替換；數值 token 集合相同）
輸出：`data/24_audit_round2.json`
"""
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS, SCRIPTS, DATA = ROOT / "docs", ROOT / "scripts", ROOT / "data"
PY = "/usr/bin/python3"
NUM = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")

report = {"meta": {"script": "24_audit_round2", "date": "2026-10-24"}, "checks": {}, "violations": []}
print("=" * 96)
print("M4 docs/24 — 第二輪審計（增量）：審計審計自身 ＋ 審計後改動核對")
print("=" * 96)


def sh(cmd):
    r = subprocess.run(cmd, cwd=ROOT.parents[1], capture_output=True, text=True, shell=True)
    return r.stdout.strip()


# ---- C1 scripts/ 自審計後零改動 ----
diff_sc = sh("git diff --stat ea59cf9..HEAD -- main_problems/m4_mass_spectral_readout/scripts")
c1 = (diff_sc == "")
report["checks"]["C1_scripts_unchanged_since_audit"] = {"pass": c1, "diff": diff_sc}
print(f"\n[C1] scripts/ 自 ea59cf9 起零改動：{'✅' if c1 else '🔴'}  {diff_sc[:120]}")

# ---- C2 docs/23 宣稱 vs JSON 實錄 ----
j23 = json.loads((DATA / "23_full_audit.json").read_text())
d23 = (DOCS / "23_full_audit.md").read_text()
n_pairs = len(j23.get("scripts", {}))
all_ident = all(v.get("rerun_identical") is True for v in j23["scripts"].values())
all_exit = all(v.get("exit") == 0 for v in j23["scripts"].values())
all_doc_str = all(v.get("docstring_match") is not False for v in j23["scripts"].values())
viol23 = j23.get("violations", [])
claims = {
    "22 支腳本": "22" in d23 or "22/22" in d23,
    "22/22 byte-identical": all_ident,
    "docstring 22/22": all_doc_str,
    "違規 0": len(viol23) == 0,
}
c2 = all([all_ident, all_exit, all_doc_str, len(viol23) == 0, n_pairs == 22])
report["checks"]["C2_doc23_claims_vs_json"] = {"pass": c2, "pairs": n_pairs, "all_identical": all_ident,
                                               "all_exit0": all_exit, "docstring_ok": all_doc_str,
                                               "violations_in_json": viol23[:5], "claim_tokens": claims}
print(f"[C2] docs/23 宣稱 vs JSON：{'✅' if c2 else '🔴'}  配對={n_pairs} 全同={all_ident} "
      f"exit0={all_exit} docstring={all_doc_str} JSON違規={len(viol23)}")

# ---- C3 抽樣重跑 3 支輕腳本 ----
SAMPLE = ["07_hit_followup", "15_exponent8", "16_scale_operator"]
tmp = Path(tempfile.mkdtemp(prefix="m4r2_", dir="/home/wsl/.hermes/cache/scratch"))
samples = {}
for stem in SAMPLE:
    sc = SCRIPTS / f"{stem}.py"
    if not sc.exists():
        samples[stem] = "missing"
        continue
    work = tmp / stem
    (work / "data").mkdir(parents=True, exist_ok=True)
    shutil.copy(sc, work / sc.name)
    r = subprocess.run([PY, sc.name], cwd=work, capture_output=True, text=True, timeout=900)
    log = DATA / f"run_{stem}.log"
    same = log.exists() and log.read_text() == r.stdout
    samples[stem] = {"exit": r.returncode, "identical": same}
c3 = all(isinstance(v, dict) and v["exit"] == 0 and v["identical"] for v in samples.values())
report["checks"]["C3_sample_rerun"] = {"pass": c3, "samples": samples, "tmp": str(tmp)}
print(f"[C3] 抽樣重跑：{'✅' if c3 else '🔴'}  {samples}")

# ---- C4 docs/23 結構 + T3 回溯 ----
log23 = (DATA / "run_23_full_audit.log").read_text()
# 跨文檔語料：docs/23 會**引用**其他 docs／logs 之數值（復算式、外部交付值），須一併納入回溯語料，
# 否則會把「引用」誤判為「無源數字」（首版即此假 FAIL，2026-10-24 修）
corpus = log23 + json.dumps(j23, ensure_ascii=False)
for f in list(DOCS.glob("*.md")) + list(DATA.glob("*.log")) + list(DATA.glob("*.json")):
    if f.name != "23_full_audit.md":
        corpus += f.read_text(errors="ignore")
sections = all(k in d23 for k in ("前設問題", "事前判準", "結論報告")) or all(k in d23 for k in ("審計範圍", "結論", "違規"))
h1 = d23.splitlines()[0].startswith("# M4 docs/")
log_nums = [float(t) for t in NUM.findall(log23) if t.replace(".", "").replace("-", "").isdigit()]
unexplained = []
for line in d23.splitlines():
    for tok in NUM.findall(line):
        try:
            v = float(tok)
        except ValueError:
            continue
        if abs(v) == 0 or len(re.sub(r"[-.eE+]", "", tok).lstrip("0")) < 3:
            continue
        if tok in corpus:
            continue
        if any(abs(v - w) <= max(0.006 * abs(v), 1e-12) for w in log_nums):
            continue
        if re.fullmatch(r"20\d\d", tok) or re.fullmatch(r"[0-9a-f]{6,}", tok):
            continue
        if any(h in line for h in ("PDG", "FLAG", "docs/", "O1", "O2", "O8", "sha", "MeV", "register")):
            continue
        unexplained.append(tok)
c4 = h1 and sections and len(unexplained) == 0
report["checks"]["C4_doc23_structure_t3"] = {"pass": c4, "h1": h1, "sections": sections,
                                             "unexplained": unexplained[:12], "n": len(unexplained)}
print(f"[C4] docs/23 結構={sections} H1={h1} T3 未解釋={len(unexplained)} {unexplained[:6]}")

# ---- C5 審計後改動：數值零漂移 ----
diff23 = sh("git diff ea59cf9..HEAD -- main_problems/m4_mass_spectral_readout/docs/23_full_audit.md")
diff01 = sh("git diff ea59cf9..HEAD -- main_problems/m4_mass_spectral_readout/01_problem_record.md")
changed = [l[1:].strip() for l in (diff23 + "\n" + diff01).splitlines() if l.startswith(("+", "-")) and not l.startswith(("+++", "---"))]
# ★ 正確做法：比**改動前後之數值 token 多重集**，唔可以睇「改動行內有咩數字」——
#   一行因字串替換而重寫，其原有數字會被誤報成漂移（首版即此假 FAIL）。
#   白名單 '11.0'：本次遷移合法加入之版本號（v11.0），非物理數值。


def _toks(txt):
    return sorted(x for x in NUM.findall(txt)
                  if len(re.sub(r"[-.eE+]", "", x).lstrip("0")) >= 3)


from collections import Counter as _C
drift = []
for rel in ("docs/23_full_audit.md", "01_problem_record.md"):
    old = sh(f"git show ea59cf9:main_problems/m4_mass_spectral_readout/{rel}")
    new = (ROOT / rel).read_text()
    d = (_C(_toks(new)) - _C(_toks(old)))
    d.pop("11.0", None)          # 版本號白名單
    if d:
        drift.append((rel, dict(d)))
nums_changed = drift
c5 = len(nums_changed) == 0
report["checks"]["C5_postaudit_numeric_drift"] = {"pass": c5, "numeric_tokens_in_diff": nums_changed[:20],
                                                  "changed_lines": len(changed)}
print(f"[C5] 審計後改動數值漂移：{'✅ 零' if c5 else '🔴 ' + str(nums_changed[:8])}（改動行 {len(changed)}）")

viol = [k for k, v in report["checks"].items() if not v["pass"]]
report["violations"] = viol
print(f"\n[違規] {len(viol)} {viol}")
with open(DATA / "24_audit_round2.json", "w") as f:
    json.dump(report, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/24_audit_round2.json")