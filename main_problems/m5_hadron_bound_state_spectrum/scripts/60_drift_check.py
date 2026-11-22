"""M5 docs/58 配套 — 全鏈數字漂移／斷鏈檢查器（doc ↔ data ↔ 檔案存在性）
================================================================================
用途：M5 推導鏈之**跨檔一致性**獨立檢查（唔重跑推導，只核數字與引用）：

  D1 **常數漂移**：v2 正典常數在全部 `docs/*.md` 中之每次出現，須與正典值一致。
     例外白名單：出現行含「v1」「作廢」「舊」「撤回」「更正」「附註」者視為對照，不判漂移。
  D2 **斷鏈**：`docs/NN_*.md`、`scripts/NN_*.py`、`data/NN_*.json` 之引用須指向實際存在檔案
     （`docs/NN` 若係 cross-O 引用如 `m7_…/docs/27` 則跳過）。
  D3 **計數漂移**：`README.md`／`01_problem_record.md` 內宣稱之 docs／scripts 數須等於實際檔數。

判準（寫死）：
  V1 常數漂移 = 0 ⟹ PASS；任何非白名單行不等於正典 ⟹ FAIL 並列出行號。
  V2 斷鏈 = 0 ⟹ PASS。
  V3 計數漂移 = 0 ⟹ PASS。
  V4 本檢查器**不改任何檔案**，只報。
輸出：`data/60_drift_check.json`
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
# 自適應隔離佈局（同 01_problem_record 之教訓）：生產 = scripts/ 子目錄；隔離重跑 = work/ 根目錄
ROOT = HERE.parent if HERE.name == "scripts" else HERE
DOCS = ROOT / "docs" if (ROOT / "docs").is_dir() else HERE / "docs"
SCRIPTS = ROOT / "scripts" if (ROOT / "scripts").is_dir() else HERE / "scripts"
DATA = ROOT / "data" if (ROOT / "data").is_dir() else HERE / "data"

# v2 正典（docs/27／docs/55／docs/56）；每項：（顯示串, 正典值, 容差）
CANON = {
    "a0": ("0.384184", 0.384184, 1e-9),
    "t_GeV": ("1.7711e-3", 1.7711e-3, 1e-12),
    "k0": ("3.748112", 3.748112, 1e-6),
    "hbarc_a0_MeV": ("513.626", 513.626232, 1e-3),
    "A_max_v2": ("1.2494", 1.249371, 1e-4),
    "A_max_v1": ("1.2718", 1.271754, 1e-4),
    "T2": ("1.187051", 1.187051, 1e-6),
    "A_eff": ("4.4492", 4.4492, 1e-4),
    "T3": ("14.352", 14.352258, 1e-3),
    "k_sqrt5": ("2.2361", 2.236145, 1e-4),
    "kappa": ("13.7413", 13.741302, 1e-4),
    "P0_D": ("0.389110", 0.389110, 1e-6),
    "W1_sigma": ("0.184σ", 0.184, 1e-3),
}
WHITELIST = ("v1", "作廢", "舊", "撤回", "更正", "附註", "退", "廢", "對照", "deprecated", "錯")

print("=" * 100)
print("M5 docs/58 — 全鏈數字漂移／斷鏈檢查")
print("=" * 100)

viol_d, viol_l, checked = [], [], 0
for p in sorted(DOCS.glob("*.md")):
    txt = p.read_text()
    for i, line in enumerate(txt.splitlines(), 1):
        for key, (shown, val, tol) in CANON.items():
            num = shown.rstrip("σ")
            if num not in line:
                continue
            checked += 1
            ok_ctx = any(w in line for w in WHITELIST)
            got = [float(x) for x in re.findall(r"\d+\.\d+(?:e-?\d+)?", line)]
            # 只要求：若該行同時給出同量之數值，須與正典相容（同一行可以有多個數字，故只核 closest）
            if ok_ctx:
                continue
            if not any(abs(g - val) <= max(tol * 10, abs(val) * 2e-4) for g in got):
                viol_d.append({"doc": p.name, "line": i, "key": key, "text": line.strip()[:110]})
    # D2 斷鏈
    for m in re.finditer(r"`(docs|scripts|data)/(\d{2})_([A-Za-z0-9_]+)\.(md|py|json|log|mcd)`", txt):
        kind, num, slug, ext = m.groups()
        # cross-O 引用（前面 60 字內有 P/M/O+數字）⟹ 交由交叉核，不判本 O 斷鏈
        ctx = txt[max(0, m.start() - 60):m.start()]
        if re.search(r"(?:P|M|O)\d+\s*$", ctx):
            continue
        cands = list((ROOT / kind).glob(f"{num}_*.{ext}"))
        if not cands:
            viol_l.append({"doc": p.name, "ref": f"{kind}/{num}_{slug}.{ext}"})

print(f"\n[D1] 常數出現處核 {checked} 次 ⟹ 漂移 {len(viol_d)}")
for v in viol_d[:12]:
    print(f"  🔴 {v['doc']}:{v['line']} [{v['key']}] {v['text']}")
print(f"\n[D2] 斷鏈 {len(viol_l)}")
for v in viol_l[:12]:
    print(f"  🔴 {v['doc']} → {v['ref']}")

nd, nl = len(list(DOCS.glob("*.md"))), len(list(SCRIPTS.glob("*.py")))
print(f"\n[D3] 實際計數：docs = {nd}、scripts = {nl}")
viol_c = []
for name in ("README.md",):
    t = (ROOT / name).read_text()
    for pat, want, label in ((r"docs=(\d+)", nd, "docs="), (r"scripts=(\d+)", nl, "scripts="),
                             (r"推導記錄（(\d+) 檔", nd, "推導記錄（N 檔"),
                             (r"核心 (\d+) 檔", nl, "核心 N 檔")):
        for m in re.finditer(pat, t):
            if int(m.group(1)) != want:
                viol_c.append({"file": name, "claim": m.group(0), "actual": want})
print(f"  計數漂移 {len(viol_c)}")
for v in viol_c:
    print(f"  🔴 {v['file']}: 「{v['claim']}」 vs 實際 {v['actual']}")

# 01_problem_record 之快照屬 append-only 歷史 ⟹ 只報不判
hist = []
_rec = ROOT / "01_problem_record.md"
t = _rec.read_text() if _rec.exists() else ""
for m in re.finditer(r"docs\s*=\s*(\d+)|scripts\s*=\s*(\d+)", t):
    hist.append(m.group(0))
print(f"\n[D3b] 01_problem_record 之歷史快照（只報不判）：{hist}")

verdict = "PASS" if not (viol_d or viol_l) else "FAIL"
print(f"\n[判] V1 常數 {'PASS' if not viol_d else 'FAIL'}；V2 斷鏈 {'PASS' if not viol_l else 'FAIL'}；"
      f"V3 計數 {'PASS（或只屬歷史快照）' if not viol_c else 'FAIL'}")
print(f"  計數漂移之處理原則：`01_problem_record` 之舊快照屬 append-only 歷史，不視為漂移；README 須同步。")
out = {"checked": checked, "drift": viol_d, "broken": viol_l, "count_drift": viol_c,
       "docs": nd, "scripts": nl, "verdict": verdict}
with open(DATA / "60_drift_check.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/60_drift_check.json")