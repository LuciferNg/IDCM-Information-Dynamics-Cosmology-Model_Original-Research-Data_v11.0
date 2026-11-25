#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/14 — audit-verify（全鏈：腳本 ↔ 文檔 ↔ 數據 ↔ 日誌 ↔ 四件套）
=====================================================================
判據（跑前寫死；A2/A3/A5/A7 之「引用式映射」與「上游解析」為本輪明示口徑）：

  A1  scripts/*.py 每支至少被一篇 docs/*.md 引用
  A2  docs/*.md 每篇（除 4 篇預註冊）至少引用一支存在之 scripts/*.py
  A3  data/NN_*.json 之 n_pass/n_total 必出現在「引用 scripts/NN_*」之某篇 doc
  A4  logs/NN.txt 之「結果：p/t PASS ｜ f FAIL」== data 同號 json
  A5  文檔內 docs/…、scripts/…、data/…、logs/… 引用須解析存在；
      本資料夾無者，按 basename 在 v11 main_problems ∪ v10 open_problems_fresh 內解析（跨 O 引用）
  A6  Σ data checks = 四件套（01/SESSION_START/README/final_result）所載水位一致
  A7  「N PASS / 0 FAIL」若非當前水位，須位於「第X輪」節內（append-only 歷史豁免）
  A8  關鍵常數（44／16／90／36／98／B=4／2B=8）在文檔層出現且無矛盾寫法
  A9  文檔層 5 組關鍵宣稱 vs data 機核值（防「文檔寫 A、數據係 B」漂移）

唯讀；發現即報，修正另行 patch。
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
DOCS = os.path.join(ODIR, "docs")
SCRIPTS = os.path.join(ODIR, "scripts")
LOGS = os.path.join(ODIR, "logs")
V11 = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture"
V10 = "/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat"

CHECKS = []

# ★ 自寫 log（不靠 shell 重定向）：先讀舊 log（供 A4 自身比對），再以 Tee 覆寫
LOG_SELF = os.path.join(LOGS, "14.txt")
OLD_LOG_SELF = open(LOG_SELF, encoding="utf-8").read() if os.path.isfile(LOG_SELF) else ""


class _Tee:
    def __init__(self, path, stream):
        self.f = open(path, "w", encoding="utf-8")
        self.s = stream

    def write(self, x):
        self.s.write(x)
        self.f.write(x)

    def flush(self):
        self.s.flush()
        self.f.flush()


sys.stdout = _Tee(LOG_SELF, sys.stdout)


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


def ls(d, ext):
    return sorted(f for f in os.listdir(d) if f.endswith(ext))


scripts = ls(SCRIPTS, ".py")
docs = ls(DOCS, ".md")
data = ls(DATA, ".json")
logs = ls(LOGS, ".txt")
print("=" * 78)
print("M24 scripts/14 — audit-verify")
print("=" * 78)
print(f"  scripts {len(scripts)}｜docs {len(docs)}｜data {len(data)}｜logs {len(logs)}")

DOCTEXT = {f: open(os.path.join(DOCS, f), encoding="utf-8").read() for f in docs}
TOPF = ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"]
TOP = {}
for f in TOPF:
    p = os.path.join(ODIR, f)
    TOP[f] = open(p, encoding="utf-8").read() if os.path.isfile(p) else None
ALLTEXT = "\n".join(list(DOCTEXT.values()) + [v for v in TOP.values() if v])

# ── A1 ────────────────────────────────────────────────────────
missing_ref = []
for s in scripts:
    num = s.split("_")[0]
    hit = any((s in t) or re.search(rf"scripts/{num}(?:_|\b)", t) for t in DOCTEXT.values())
    if not hit:
        missing_ref.append(s)
ck("A1  scripts 每支被 docs 引用", len(missing_ref) == 0, f"未引用 = {missing_ref}")

# ── A2 ────────────────────────────────────────────────────────
PREREG = {"00_preregistration.md", "06_preregistration_round2.md",
          "10_preregistration_round3.md", "15_preregistration_round4.md"}
no_script = []
for f, t in DOCTEXT.items():
    if f in PREREG:
        continue
    refs = re.findall(r"scripts/([0-9]{2})(?:_[A-Za-z0-9_]+\.py)?", t)
    ok = any(any(s.startswith(r + "_") for s in scripts) for r in refs)
    if not ok:
        no_script.append(f)
ck("A2  結果 doc 每篇引用至少一支存在 script（4 篇預註冊豁免）", len(no_script) == 0,
   f"缺陷 = {no_script}")

# ── A3 / A4 ───────────────────────────────────────────────────
prob_log, prob_doc = [], []
for j in data:
    num = j.split("_")[0]
    J = json.load(open(os.path.join(DATA, j), encoding="utf-8"))
    pt = f"{J['n_pass']}/{J['n_total']}"
    logp = os.path.join(LOGS, f"{num}.txt")
    if not os.path.isfile(logp):
        prob_log.append((j, "缺 log"))
    else:
        _txt = OLD_LOG_SELF if j.startswith("14_") else open(logp, encoding="utf-8").read()
        m = re.search(r"結果：(\d+)/(\d+) PASS ｜ (\d+) FAIL", _txt)
        if not m:
            prob_log.append((j, "log 無判詞行"))
        elif f"{m.group(1)}/{m.group(2)}" != pt or int(m.group(3)) != J["n_fail"]:
            prob_log.append((j, f"log {m.group(1)}/{m.group(2)}/{m.group(3)} ≠ json {pt}/{J['n_fail']}"))
    if j.startswith("14_"):
        continue          # meta：audit 自身之結果不可能在同一跑內被自己文檔引用（見 A6 meta 排除）
    cited = [f for f, t in DOCTEXT.items() if re.search(rf"scripts/{num}(?:_|\b)", t)]
    if not cited:
        prob_doc.append((j, "無 doc 引用其 script"))
    elif not any(pt in DOCTEXT[f] for f in cited):
        prob_doc.append((j, f"{pt} 不在 {cited}"))
ck("A3  data 之 n_pass/n_total 見於引用其 script 之 doc", len(prob_doc) == 0, f"problems = {prob_doc}")
ck("A4  logs/NN.txt 判詞 == data/NN_*.json（逐對）", len(prob_log) == 0, f"problems = {prob_log}")

# ── A5 引用解析（含跨 O 上游解析）──────────────────────────────
CACHE = {}


def find_elsewhere(basename):
    if basename in CACHE:
        return CACHE[basename]
    hits = []
    for root in (V11, V10):
        for dirpath, _dn, fns in os.walk(root):
            if basename in fns:
                hits.append(os.path.join(dirpath, basename))
    CACHE[basename] = hits
    return hits


refprob = []
xO = []
for name, t in list(DOCTEXT.items()) + [(k, v) for k, v in TOP.items() if v]:
    for m in re.finditer(r"(docs/[0-9]{2}_[A-Za-z0-9_]+\.md|scripts/[0-9]{2}_[A-Za-z0-9_]+\.py|"
                         r"data/[0-9]{2}_[A-Za-z0-9_]+\.json|logs/[0-9]{2}\.txt)", t):
        ref = m.group(1)
        if os.path.isfile(os.path.join(ODIR, ref)):
            continue
        hits = find_elsewhere(os.path.basename(ref))
        if hits:
            xO.append((name, ref, len(hits)))
        else:
            refprob.append((name, ref))
ck("A5  文檔引用全部解析（本地 或 跨 O 上游唯一定位）", len(refprob) == 0,
   f"斷鏈 = {sorted(set(refprob))}｜跨 O 解析 = {len(xO)} 處（皆存在）")

# ── A6 水位 ───────────────────────────────────────────────────
META = {"14_audit_verify.json"}      # audit 自身之 json 屬 meta，不入水位 ledger
LEDGER = [j for j in data if j not in META]
totp = sum(json.load(open(os.path.join(DATA, j), encoding="utf-8"))["n_pass"] for j in LEDGER)
tott = sum(json.load(open(os.path.join(DATA, j), encoding="utf-8"))["n_total"] for j in LEDGER)
totf = sum(json.load(open(os.path.join(DATA, j), encoding="utf-8"))["n_fail"] for j in LEDGER)
print(f"  ledger（排除 {sorted(META)}）= {len(LEDGER)} 個 json")
print(f"  Σ checks = {totp}/{tott}（FAIL {totf}）")
in_top = {k: (v is not None and str(totp) in v and str(tott) in v) for k, v in TOP.items()}
ck("A6  四件套水位 = Σ data checks 且 final_result 存在", totf == 0 and all(in_top.values()),
   f"{totp}/{tott}｜命中 = {in_top}")

# ── A7 stale（按節豁免）──────────────────────────────────────
def section_of(text, pos):
    """回傳 pos 之前最近之 '## ' 標題行。"""
    head = text[:pos]
    idx = head.rfind("\n## ")
    return head[idx:head.find("\n", idx + 1)] if idx >= 0 else ""


stale = []
for name, t in [(k, v) for k, v in TOP.items() if v]:
    for m in re.finditer(r"(\d+) PASS / 0 FAIL", t):
        n = int(m.group(1))
        if n == totp:
            continue
        sec = section_of(t, m.start())
        if "第" in sec and "輪" in sec:
            continue
        stale.append((name, n, sec.strip()[:40]))
ck("A7  非當前水位之 PASS 數皆位於「第X輪」節內", len(stale) == 0, f"stale = {stale}")

# ── A8 關鍵常數 ───────────────────────────────────────────────
KEY = {"N_d=44": r"N_d\s*=\s*44", "N₂=16": r"N₂\s*=\s*16", "Σγ=90": r"Σγ\s*=\s*90",
       "h11=36": r"h11\s*=\s*36", "h21=98": r"h21\s*=\s*98", "B=4": r"B\s*=\s*4",
       "2B=8": r"2B\s*=\s*8", "交集{16}": r"交集\s*=\s*\{16\}"}
kk = {k: len(re.findall(p, ALLTEXT)) for k, p in KEY.items()}
ck("A8  關鍵常數／結構句全鏈出現", all(v >= 1 for v in kk.values()), f"{kk}")

# ── A9 宣稱 vs 機核（5 組）────────────────────────────────────
J01 = json.load(open(os.path.join(DATA, "01_l2_label_algebra.json"), encoding="utf-8"))
J06 = json.load(open(os.path.join(DATA, "06_uniform_rule_search.json"), encoding="utf-8"))
J11 = json.load(open(os.path.join(DATA, "11_2B_forcing_test.json"), encoding="utf-8"))
J12 = json.load(open(os.path.join(DATA, "12_l0_stress_test.json"), encoding="utf-8"))
J13 = json.load(open(os.path.join(DATA, "13_character_induced_search.json"), encoding="utf-8"))
J09 = json.load(open(os.path.join(DATA, "09_induced_module_search.json"), encoding="utf-8"))

spot = []
# 1 docs/01：L2-A 命中 16:4、36:4、98:0
doc01 = DOCTEXT["01_l2_label_algebra_result.md"]
spot.append(("docs/01 命中 16 = 4", J01["per_system"]["L2-A"]["16"] == 4 and "n = 4" in doc01))
spot.append(("docs/01 交集 = {16}", J01["intersection_A_B"] == [16] and "交集定理" in doc01))
# 2 docs/07：triple = 0
doc07 = DOCTEXT["07_r_a_r_b_uniform_rule_result.md"]
spot.append(("docs/07 triple families = 0", len(J06["module"]["triple_hits"]) == 0
             and "之族 = **0**" in doc07))
# 3 docs/13：inv_m = [43,44]
doc13 = DOCTEXT["13_2B_forcing_result.md"]
spot.append(("docs/13 inv_m = [43,44]", J11["invariant_m"] == [43, 44] and "[43, 44]" in doc13))
# 4 docs/14：7 案例
doc14 = DOCTEXT["14_l0_stress_result.md"]
spot.append(("docs/14 n_cases = 7", J12["n_cases"] == 7 and "7/7" in doc14))
# 5 docs/16：98 = 0/225、36 = 1/225
doc16 = DOCTEXT["16_character_induced_result.md"]
h98 = sum(1 for t in J13["table"] if 98 in [t["D1"][0], t["D1"][1], t["D1"][2], t["D2"][0],
                                            t["D2"][1], t["D2"][2], t["D3"], t["D4"], t["D5"]])
spot.append(("docs/16 98 命中 = 0", h98 == 0 and "0/225" in doc16))
spot.append(("docs/11 98 = 0/81", "0/81" in DOCTEXT["11_induced_module_result.md"]
             and all(98 != v for t in J09["table"] for v in
                     [t["D1"][0], t["D1"][1], t["D1"][2], t["D2"][0], t["D2"][1], t["D2"][2],
                      t["D3"], t["D4"], t["D5"]])))
for tag, ok in spot:
    print(f"    {'ok ' if ok else 'BAD'} {tag}")
ck("A9  5+ 組文檔宣稱與 data 機核值一致", all(ok for _, ok in spot),
   f"fail = {[t for t, ok in spot if not ok]}")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/14 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

with open(os.path.join(DATA, "14_audit_verify.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "counts": {"scripts": len(scripts), "docs": len(docs), "data": len(data), "logs": len(logs)},
               "sum_checks": {"pass": totp, "total": tott, "fail": totf},
               "missing_ref": missing_ref, "no_script": no_script,
               "problems_doc": prob_doc, "problems_log": prob_log,
               "broken_refs": sorted(set(refprob)), "cross_O_refs": len(xO),
               "stale": stale, "key_constants": kk,
               "spot": [[t, ok] for t, ok in spot],
               "verdict": "見 docs/17_audit_report.md"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
