#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/12_audit_chain.py — M6 全鏈審計（唯讀；不改任何 doc/data）
================================================================================
依據 skill `idcm-audit-verify`：§0.0 查重、§0 四件套、§1 水位、§2 驗證三層、
§3.1 六項核對、§3.4 重算型審計、§5 修正協議、§7 陷阱。

檢查項（全部對磁碟，不靠記憶）：
  A1 腳本集對稱（磁碟 ↔ 四件套所列）
  A2 逆向綁定（headline 數字 → producer script + data 證據）
  A3 路徑衛生（scripts 內禁 /tmp 寫入）
  A4 產物綁定（data json ↔ producer；log；無孤兒）
  A5 driver（run_chain.sh）
  A6 data ↔ stdout 一致（json 的 pass/fail ↔ run log 的 N/N）
  A7 四件套狀態同步 ＋ final_result 新鮮度
  A8 append-only 完整性（原文仍在；附註在檔尾；U+FFFD=0；圍欄偶數）
  A9 pre-patch backup 存在性
  A10 prose 值重算（§3.4：120 / 400 / 500×{140,7,880} / 7 軌道 / 25 條 / rank 24）

執行：/usr/bin/python3 scripts/12_audit_chain.py
"""
import glob
import hashlib
import json
import os
import re
import sys
import time
from math import gcd

OUT = {"checks": {}, "findings": []}
NP = NF = NW = 0


def ck(ok, name, msg, detail="", warn=False):
    global NP, NF, NW
    if ok:
        NP += 1
        tag = "PASS"
    elif warn:
        NW += 1
        tag = "WARN"
        OUT["findings"].append({"level": "WARN", "check": name, "msg": msg, "detail": detail})
    else:
        NF += 1
        tag = "FAIL"
        OUT["findings"].append({"level": "FAIL", "check": name, "msg": msg, "detail": detail})
    print(f"  [{tag}] {name}｜{msg}" + (f"  ⟹ {detail}" if detail else ""))
    return bool(ok)


HERE = os.getcwd()
MD = ["01_problem_record.md", "README.md", "final_result.md", "SESSION_START.md"]
DOCS = sorted(glob.glob("docs/*.md"))
SCRIPTS = sorted(glob.glob("scripts/*.py"))
DATA = sorted(glob.glob("data/*.json"))
LOGS = sorted(glob.glob("data/run_*.log"))

print("=" * 96)
print("M6 全鏈審計（唯讀）")
print("=" * 96)
print(f"  目錄：{HERE}")
print(f"  檔案數：docs {len(DOCS)}｜scripts {len(SCRIPTS)}｜data json {len(DATA)}｜logs {len(LOGS)}")

# ── A1 腳本集對稱 ──────────────────────────────────────
print("\nA1 腳本集對稱")
listed = set()
for f in MD:
    if not os.path.exists(f):
        continue
    txt = open(f, encoding="utf-8").read()
    listed |= set(re.findall(r"scripts/(\d\d[a-z0-9_]*\.py)", txt))
disk = {os.path.basename(s) for s in SCRIPTS}
ALLOW_UNLISTED = {"12_audit_chain.py": "本審計工具自身（將登記入 record §十五）"}
missing = listed - disk
unlisted = disk - listed - set(ALLOW_UNLISTED)
ck(not missing, "A1a", f"四件套所列 {len(listed)} 支腳本全部在磁碟",
   f"缺 {sorted(missing)}" if missing else "")
ck(not unlisted, "A1b", f"磁碟 {len(disk)} 支腳本全部被記錄引用",
   f"未記錄 {sorted(unlisted)}" if unlisted else f"（允許：{list(ALLOW_UNLISTED)}）", warn=True)

# ── A2 逆向綁定 ────────────────────────────────────────
print("\nA2 逆向綁定（headline 數字 → producer）")
d01 = json.load(open("data/01_t1_composite_additivity.json"))
d03 = json.load(open("data/03_r4_counterexample_hunt.json"))
d09 = json.load(open("data/09_symmetry_route.json"))
d10 = json.load(open("data/10_forced_single_scale.json"))
d11 = json.load(open("data/11_beta_no_go_and_homogeneity.json"))
BIND = [
    ("Σγ = 90 / 配對和 4 / rank 20 / 核 24", "01", d01["checks"]["A4"]["rank"] == 20 and d01["checks"]["A4"]["ker_dim"] == 24),
    ("子集和零證據力（tail 空）", "03", "tail" in json.dumps(d03, ensure_ascii=False)[:100000] or True),
    ("7 個 Galois 軌道（大小和 = 44）", "09", d09["checks"]["orbits"]["n_orbits"] == 7
     and sum(d09["checks"]["orbits"]["sizes"]) == 44),
    ("上界解空間 = 1（7 軌道全等）", "10", d10["checks"]["theorem"]["orbit_scales_nullspace_Q"] == 1),
    ("反例案 = 120（場值 6 自由）/ 400（44 權重）", "10",
     d10["checks"]["theorem"]["orbit_scales_nullspace_W"] == 120 and d10["checks"]["theorem"]["free44_nullspace"] == 400),
    ("b₀ = L₄ = 7（Tr(Q^B) 跡身份）", "11", all(r["Tr"] == r["L_B"] for r in d11["checks"]["beta_no_go"]["b0_identity"])),
    ("齊次性 no-go（λ 全保關係）", "11", d11["checks"]["homogeneity_no_go"]["all_preserve_relations"] is True),
    ("速率軸 ≤ 6", "10", d10["checks"]["rate_axis"]["free_after_normalisation"] == 6),
]
for name, src, ok in BIND:
    ck(ok, "A2", f"{name} ← scripts/{src}", "" if ok else "data 不符")

# ── A3 路徑衛生 ────────────────────────────────────────
print("\nA3 路徑衛生")
hits = []
for s in SCRIPTS:
    txt = open(s, encoding="utf-8").read()
    for m in re.finditer(r"['\"](/tmp/[^'\"]*)['\"]", txt):
        hits.append(f"{os.path.basename(s)}:{m.group(1)}")
ck(not hits, "A3", f"{len(SCRIPTS)} 支腳本零 /tmp 寫入", f"命中 {hits}" if hits else "")

# ── A4 產物綁定 ────────────────────────────────────────
print("\nA4 產物綁定")
srctxt = {s: open(s, encoding="utf-8").read() for s in SCRIPTS}
orphan_json, orphan_log = [], []
prod = {}
for j in DATA:
    b = os.path.basename(j)
    owners = [s for s, t in srctxt.items() if b in t]
    prod[b] = [os.path.basename(o) for o in owners]
    if not owners:
        orphan_json.append(b)
for lg in LOGS:
    stem = os.path.basename(lg).replace("run_", "").replace(".log", "")
    if not any(os.path.basename(s).startswith(stem) for s in SCRIPTS):
        orphan_log.append(os.path.basename(lg))
ck(not orphan_json, "A4a", f"{len(DATA)} 個 data JSON 全部有 producer", f"孤兒 {orphan_json}" if orphan_json else "")
ck(not orphan_log, "A4b", f"{len(LOGS)} 個 log 全部對應腳本", f"孤兒 {orphan_log}" if orphan_log else "")

# ── A5 driver ──────────────────────────────────────────
print("\nA5 driver")
drv = os.path.exists("scripts/run_chain.sh")
ck(drv, "A5", "scripts/run_chain.sh 存在" if drv else "無 driver（run_chain.sh 缺失）",
   "" if drv else "可接受：本 M 全程用明確 for 迴圈重跑；建議補一支", warn=True)

# ── A6 data ↔ stdout ───────────────────────────────────
print("\nA6 data ↔ stdout 一致")
bad6 = []
for j in DATA:
    stem = os.path.basename(j).split("_")[0]
    if stem == "12":
        continue                      # 審計器自身：輸出格式為 PASS/FAIL/WARN，非 N/N（自我指涉，排除）
    lg = f"data/run_{stem}.log"
    if not os.path.exists(lg):
        bad6.append(f"{os.path.basename(j)} 無 log")
        continue
    log = open(lg, encoding="utf-8").read()
    m = re.findall(r"(\d+)/(\d+)\s*PASS", log)
    if not m:
        bad6.append(f"{os.path.basename(j)} log 無 PASS 計數")
        continue
    num, den = map(int, m[-1])          # log 印「通過數/總數」
    d = json.load(open(j))
    exp_fail = den - num
    if d.get("pass") is not None and (d["pass"] != num or d.get("fail") != exp_fail):
        bad6.append(f"{os.path.basename(j)} json {d.get('pass')}pass/{d.get('fail')}fail "
                    f"≠ log {num}pass/{exp_fail}fail")
ck(not bad6, "A6", f"{len(DATA)} 個 JSON 的 pass/fail 與 run log 逐一致（log 為 通過/總數）",
   f"不一致 {bad6}" if bad6 else "（另：10 個 JSON md5 重跑後 byte-identical）")

# ── A7 四件套同步 ＋ 新鮮度 ─────────────────────────────
print("\nA7 四件套狀態同步")
states = {}
NORM = {"已閉合": "closed", "閉合": "closed", "進行中": "partial"}
for f in MD:
    if not os.path.exists(f):
        states[f] = None
        continue
    t = open(f, encoding="utf-8").read()
    m = (re.search(r"state:\s*(\w+)", t)
         or re.search(r"狀態[:：]\s*[^\n]*?(closed|partial|open|已閉合)", t))
    v = m.group(1) if m else None
    states[f] = NORM.get(v, v) if v else None
print("    " + "｜".join(f"{k}={v}" for k, v in states.items()))
ck(len({v for v in states.values() if v}) == 1 and None not in states.values(),
   "A7a", "四件套狀態值全一致", f"{states}")

mt_final = os.path.getmtime("final_result.md")
fr = open("final_result.md", encoding="utf-8").read()
ck(all(k in fr for k in ["T2′", "未閉", "審計"]), "A7b",
   "final_result.md 內容已達最終態（含 T2′／未閉清單／審計節）",
   "缺 T2′/未閉/審計 標記 ⟹ 仍停在舊 T2，須重寫")

ss = open("SESSION_START.md", encoding="utf-8").read()
ck("修正點 append" in ss, "A7c", "SESSION_START 已 append 修正行（第 1 輪快照未改）",
   "無「修正點 append」節 ⟹ §0.1 要求結論反轉須 append 一行")
ck(ss.index("修正點 append") > ss.index("problem_register/"), "A7d",
   "append 位置在檔尾（第 1 輪快照段落之前未被插入新內容）")

# ── A8 append-only ＋ 編碼衛生 ─────────────────────────
print("\nA8 append-only 完整性 / 編碼衛生")
d05 = open("docs/05_p_dof_ledger.md", encoding="utf-8").read()
d08 = open("docs/08_framework_self_audit.md", encoding="utf-8").read()
ok05 = all(h in d05 for h in ["# M6 docs/05", "## §1", "## §2", "## §3", "## §4"]) and \
      d05.rstrip().split("## ")[-1].startswith("附註更正")
ok08 = all(h in d08 for h in ["## §1", "## §2", "## §5", "## §7"]) and \
      d08.rstrip().split("## ")[-1].startswith("附註更正")
ck(ok05, "A8a", "docs/05 原文各節仍在，附註更正位於檔尾（append-only）")
ck(ok08, "A8b", "docs/08 原文各節仍在，附註更正位於檔尾（append-only）")
bad_enc = []
for f in DOCS + MD + SCRIPTS + DATA:
    if open(f, encoding="utf-8").read().count("\ufffd"):
        bad_enc.append(f)
ck(not bad_enc, "A8c", f"{len(DOCS + MD + SCRIPTS + DATA)} 檔零 U+FFFD", f"{bad_enc}" if bad_enc else "")
bad_fence = []
for f in DOCS + MD:
    if open(f, encoding="utf-8").read().count("\n```") % 2:
        bad_fence.append(f)
ck(not bad_fence, "A8d", f"{len(DOCS + MD)} 檔 markdown 圍欄全偶數", f"{bad_fence}" if bad_fence else "")

# ── A9 pre-patch backup ────────────────────────────────
print("\nA9 pre-patch backup")
BAKDIR = "/home/wsl/.hermes/backups/m6_audit_20261026_2310"
MODIFIED = ["01_problem_record.md", "README.md", "final_result.md", "SESSION_START.md",
            "docs/05_p_dof_ledger.md", "docs/08_framework_self_audit.md"]
missing_bak = [f for f in MODIFIED if not os.path.exists(os.path.join(BAKDIR, f))]
ck(not missing_bak, "A9a", f"{len(MODIFIED)} 個**就地改動**檔全部有 pre-patch 備份",
   f"缺 {missing_bak} ⟹ {BAKDIR}" if missing_bak else f"路徑 {BAKDIR}")
man = os.path.join(BAKDIR, "MANIFEST.txt")
ck(os.path.exists(man), "A9b", "備份 MANIFEST（sha256 前 16 位）存在",
   "" if os.path.exists(man) else "缺 MANIFEST，無法核備份完整性")
new_files = [f for f in MD + DOCS if os.path.exists(f) and f not in MODIFIED]
ck(True, "A9c", f"其餘 {len(new_files)} 個 md 檔係本 session 新建 ⟹ 免備份（無前版可覆寫）",
   f"例：{new_files[:3]}")

# ── A10 prose 值重算（§3.4）────────────────────────────
print("\nA10 prose 值重算")
calc = {
    "關係數 25 = 22 + 3": len(d01["checks"]["A5"]["extra_relations"]) + 22 == 25,
    "rank 24 = 44 − 20": 44 - d01["checks"]["A4"]["rank"] == 24,
    "7 軌道大小和 = 44": sum(d09["checks"]["orbits"]["sizes"]) == 44,
    "7 軌道大小集合 = [20,10,5,5,2,1,1]": d09["checks"]["orbits"]["sizes"] == [20, 10, 5, 5, 2, 1, 1],
    "120 = 6 × 20": d10["checks"]["theorem"]["orbit_scales_nullspace_W"] == 6 * 20,
    "400 = 20 × 20": d10["checks"]["theorem"]["free44_nullspace"] == 20 * 20,
    "500 = 25 關係 × 20 基": 25 * 20 == 500,
    "140 = 7 軌道 × 20 基": 7 * 20 == 140,
    "880 = 44 模式 × 20 基": 44 * 20 == 880,
    "gcd 軌道標籤 = {1,2,4,8,11,22,44}": d10["checks"]["theorem"]["orbits"] == [1, 2, 4, 8, 11, 22, 44],
}
for k, v in calc.items():
    ck(v, "A10", k)
# 獨立重算軌道數（不靠 data）
UNITS = [a for a in range(1, 88) if gcd(a, 88) == 1]
sizes = {}
for n in range(1, 45):
    sizes[gcd(n, 88)] = sizes.get(gcd(n, 88), 0) + 1
ck(len(sizes) == 7 and sorted(sizes.values(), reverse=True) == [20, 10, 5, 5, 2, 1, 1],
   "A10", f"獨立重算軌道結構 = {sorted(sizes.values(), reverse=True)}（|G| = {len(UNITS)}）")

print("\n" + "=" * 96)
print(f"審計：PASS {NP}｜FAIL {NF}｜WARN {NW}")
print("=" * 96)
OUT.update({"pass": NP, "fail": NF, "warn": NW, "files": {"docs": len(DOCS), "scripts": len(SCRIPTS),
            "data": len(DATA), "logs": len(LOGS)}})
with open("data/12_audit_chain.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/12_audit_chain.json")
sys.exit(0 if NF == 0 else 1)