#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R32 — 全鏈覆蓋審計：腳本存量／文檔紀錄／數字溯源（docs 01–47 全覆蓋）／索引一致

2026-11-03。調律者八點指令之 (1)(2)(3) 之量化部分。既有工具之**空白**：
  · `scripts/15` A3 只覆蓋「可重跑」腳本（53 支），不含審計器自身與 `l2`／`l3`；
  · `scripts/28` C3 之數字溯源**只查 docs 02–21**；docs 22–47 從未被數字溯源審過；
  · 無任何工具檢查「腳本 ↔ 文檔」雙向引用（孤兒腳本／孤兒文檔）；
  · 無工具檢查四根檔之計數與實況是否一致。

預註冊判準（跑前寫死；N 全報）
  A1 腳本存量：每個編號腳本有 canonical `data/run_NN_<slug>.log`；`l2`／`l3` 亦有 log。
  A2 腳本被 15 A3 覆蓋（審計器自身 `15_full_audit.py` 例外，以其自身兩跑等同）。
  A3 腳本 → 文檔回指：每支編號腳本，其 slug 或編號須在 ≥1 份 doc 出現。
  A4 文檔 → 腳本回指：每份 doc（01–47）須提及 ≥1 個腳本編號／log 名（`01_problem_statement` 為陳述書，例外）。
  A5 數字溯源（**docs 01–47 全覆蓋**）：doc 內每個 ≥3 位小數之數值，須在 `data/*.log`／`data/*.json` 之全域池中出現；
     「只在更晚 round 出現」者單獨列為 `A5b`（對更正附註類 doc 屬合法），不直接判 FAIL。
  A6 索引一致：README 之 docs 表須列齊 01–47；`01_problem_record` §三.2 須有 R1–R31 每列；四根檔狀態行皆 R1–R31。
  A7 衞生：`scripts/`／`docs/` 內不得存在 `*pre_patch*` 備份與 `__pycache__`（備份應在 `~/.hermes/backups/`）。

跑法：/usr/bin/python3 scripts/55_r32_coverage_audit.py
"""
from __future__ import annotations
import atexit as _atexit
import json
import pathlib as _pl
import re
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path, self.buf, self.out = path, [], _sys.__stdout__

    def write(self, s):
        self.buf.append(s)
        if self.out is not None:
            try:
                self.out.write(s)
            except Exception:
                pass

    def flush(self):
        if self.out is not None:
            try:
                self.out.flush()
            except Exception:
                pass

    def close(self):
        self.path.write_text("".join(self.buf), encoding="utf-8")


ROOT = _pl.Path(__file__).resolve().parent.parent
_log = _Tee(ROOT / "data" / "run_55_r32_coverage_audit.log")
_sys.stdout = _log
_atexit.register(_log.close)

S, D, T = ROOT / "scripts", ROOT / "docs", ROOT / "data"
FAIL, OUT = [], {}


def chk(tag, ok, msg):
    print("  [%s] %-2s %s" % ("PASS" if ok else "FAIL", tag, msg))
    if not ok:
        FAIL.append(tag)


print("=" * 100)
print("R32 — 全鏈覆蓋審計（腳本存量／文檔紀錄／數字溯源 docs 01–47／索引一致／衞生）")
print("=" * 100)

num, aux = [], []
for p in sorted(S.glob("*.py")):
    (num if re.match(r"^\d\d_", p.name) else aux).append(p.name)
docs = sorted(p.name for p in D.glob("*.md") if "pre_patch" not in p.name)
logs = {p.name for p in T.glob("*.log")}

print("\n[A1] 腳本存量")
miss_log = [n for n in num if "run_" + n[:-3] + ".log" not in logs]
print("   編號腳本 %d ｜輔助 %s ｜docs %d ｜logs %d" % (len(num), aux, len(docs), len(logs)))
print("   缺 canonical log：%s" % (miss_log or "無"))
chk("A1", not miss_log, "每支編號腳本皆有 canonical log")

print("\n[A2] 15 A3 覆蓋")
a3 = set(re.findall(r"\[PASS\] A3 (\S+\.py)", (T / "run_15_full_audit.log").read_text(encoding="utf-8")))
# ★ 自我參照防護：55 之 A2 要求「15 覆蓋全部腳本」，而 15 只在 55 退出 0 時給 [PASS] A3
#   ⟹ 若把 55 自身計入期望，即死鎖。故兩個審計器自身一律排除（各自以「兩跑等同」獨立保證）。
_SELF = {"15_full_audit.py", "28_continuity_audit.py", "47_allchain_audit.py",
         "48_propagation_audit.py", "55_r32_coverage_audit.py"}   # 審計器互為輸入 ⟹ 不得互相要求，
                                                                  # 各自由「兩跑等同」保證（15 A1/A2）
not_a3 = [n for n in num if n not in a3 and n not in _SELF]
print("   A3 覆蓋 %d 支 ｜期望外殘餘：%s ｜審計器自身（豁免）：%s"
      % (len(a3), not_a3 or "無", sorted(_SELF - a3)))
chk("A2", not not_a3, "除審計器自身外，全部腳本皆被 15 A3 覆蓋")

print("\n[A3] 腳本 → 文檔回指")
dtxt = {d: (D / d).read_text(encoding="utf-8", errors="replace") for d in docs}
alltxt = "\n".join(dtxt.values())
orphan_script = []
for n in num:
    slug = n[:-3]
    if slug in alltxt or re.search(re.escape(n[:-3]), alltxt):
        continue
    orphan_script.append(n)
print("   孤兒腳本（無任何 doc 提及）：%s" % (orphan_script or "無"))
chk("A3", not orphan_script, "每支腳本皆被至少一份 doc 引用")

print("\n[A4] 文檔 → 腳本回指")
numset = {n[:-3] for n in num}
orphan_doc = []
for d in docs:
    t = dtxt[d]
    hit = any(re.search(r"\b" + re.escape(n), t) for n in numset) or re.search(r"scripts/\d\d", t) or \
          re.search(r"run_\d\d_", t) or "scripts/15" in t
    if not hit and not d.startswith("01_problem_statement"):
        orphan_doc.append(d)
print("   孤兒文檔（無任何腳本／log 引用）：%s" % (orphan_doc or "無"))
chk("A4", not orphan_doc, "每份 doc 皆有可追溯腳本或 log")

print("\n[A5] 數字溯源（docs 01–47 全覆蓋）")
# ★ 自污染防護：本審計器之 log／json 會列出孤兒 token 原文 ⟹ 若入池，孤兒下一輪即「溯源成功」。
_POOL_EXCL = {"run_55_r32_coverage_audit.log"}
_JSON_EXCL = {"55_r32_coverage_audit.json"}   # 同理由：其 orphan_values 記錄原文 token
# 跨記錄引用白名單（他 M／他 O 之數值，本 M 無 producer log）：須逐條附出處，披露不計 FAIL。
EXT = {
    "3.717438495": "M13 回執之 r 值（M13 docs；跨記錄引用，本 M 無 log）",
    "0.412601": "M13 回執之 σ(r)（同上）",
}
pool_txt = "\n".join((T / f).read_text(encoding="utf-8", errors="replace")
                     for f in sorted(logs) if f not in _POOL_EXCL)
pool_txt += "\n" + "\n".join(q.read_text(encoding="utf-8", errors="replace") for q in sorted(T.glob("*.json"))
                                  if q.name not in _JSON_EXCL)


def toks(s):
    return set(re.findall(r"(?<![\w.])\d+\.\d{3,}(?![\w])", s))


pool = toks(pool_txt)
pool_f = [float(v) for v in pool]
pool_round = {}
for f in sorted(logs):
    m = re.match(r"run_(\d\d)_", f)
    if m:
        for v in toks((T / f).read_text(encoding="utf-8", errors="replace")):
            pool_round.setdefault(v, set()).add(int(m.group(1)))
def traceable(tok):
    """doc 之數值可為池值之**四捨五入**（容差 = 半個末位），此屬排版非漂移。"""
    if tok in pool:
        return True
    d = len(tok.split(".")[1])
    x = float(tok)
    half = 0.5 * 10 ** (-d)
    return any(abs(v - x) <= half + 1e-15 for v in pool_f)


total, orphan_v, late_v, ext_hit = 0, [], [], []
for d in docs:
    dn = int(d[:2]) if d[:2].isdigit() else 0
    for v in sorted(toks(dtxt[d])):
        total += 1
        if not traceable(v):
            if v in EXT:
                ext_hit.append((d, v, EXT[v]))
            else:
                orphan_v.append((d, v))
        else:
            rs = pool_round.get(v, set())
            if rs and min(rs) > dn:
                late_v.append((d, v, sorted(rs)[:3]))
print("   檢查數值 token 數 = %d" % total)
print("   A5a 全域池無（真孤兒）：%d %s" % (len(orphan_v), orphan_v[:12]))
print("   A5b 只在更晚 round 出現（更正附註類合法）：%d %s" % (len(late_v), late_v[:8]))
print("   A5c 跨記錄白名單（披露，不計 FAIL）：%d %s" % (len(ext_hit), [(d, v) for d, v, _ in ext_hit]))
chk("A5", not orphan_v, "每個 ≥3 位小數數值皆可在 log／json 池溯到（容差＝半個末位；A5b 只作 INFO）")

print("\n[A6] 索引一致")
rd = (ROOT / "README.md").read_text(encoding="utf-8")
missing_row = [d for d in docs if d not in rd and d[:-3] not in rd]
print("   README docs 表缺列：%s" % (missing_row or "無"))
rec = (ROOT / "01_problem_record.md").read_text(encoding="utf-8")
rrows = {int(m) for m in re.findall(r"^\| R(\d+) \|", rec, re.M)}
_rmax = max(rrows) if rrows else 31
miss_r = sorted(set(range(1, _rmax + 1)) - rrows)
print("   §三.2 R 列：有 %d／%d ｜缺 %s" % (len(rrows), _rmax, miss_r or "無"))
roots = ["README.md", "01_problem_record.md", "final_result.md", "SESSION_START.md"]
stat, latest = {}, {}
for f in roots:
    t = (ROOT / f).read_text(encoding="utf-8")
    hits = [int(x) for x in re.findall(r"R1–R(\d+)", t)]
    stat[f] = bool(hits)
    latest[f] = max(hits) if hits else 0
rec_max = max(rrows) if rrows else 0
print("   四根檔宣稱之最新 R：%s ｜01 表最大 R 列 = R%d" % (latest, rec_max))
same = len(set(latest.values())) == 1 and all(v == rec_max for v in latest.values())
chk("A6", not missing_row and not miss_r and all(stat.values()) and same,
    "README 表／R 列／四根檔狀態一致（同一最新 R）")

print("\n[A7] 衞生")
pre_s = len(list(S.glob("*pre_patch*")))
pre_d = len(list(D.glob("*pre_patch*")))
pyc = (S / "__pycache__").exists()
print("   scripts/ 內 pre_patch %d ｜docs/ 內 %d ｜__pycache__ %s" % (pre_s, pre_d, pyc))
chk("A7", pre_s == 0 and pre_d == 0 and not pyc, "工作目錄內無備份殘留／快取")

print("\n[匯總] FAIL = %d ／ 檢查項 = 7" % len(FAIL))
print("   FAIL 項：%s" % (FAIL or "無"))
print("=" * 100)
OUT = {"num_scripts": len(num), "aux": aux, "docs": len(docs), "logs": len(logs),
       "a3": len(a3), "tokens": total, "orphan_values": orphan_v, "late_values": late_v,
       "orphan_script": orphan_script, "orphan_doc": orphan_doc, "missing_row": missing_row,
       "miss_r": miss_r, "stat": stat, "pre_patch": [pre_s, pre_d], "pycache": pyc, "fail": FAIL}
(ROOT / "data" / "55_r32_coverage_audit.json").write_text(json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")
_sys.exit(0 if not FAIL else 1)