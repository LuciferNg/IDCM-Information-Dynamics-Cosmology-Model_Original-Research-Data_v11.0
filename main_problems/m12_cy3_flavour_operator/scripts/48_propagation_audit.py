#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
48 - M12 傳播審計 (propagation audit of the all-chain audit, R26 級)

目的
----
全鏈審計 (scripts/47, A1-A10) 修好四件套之後，檢查「修正有沒有下游影響」。
本腳本唔重驗 M12 內部推導 (由 15/28/47/l2/l3 負責)，只驗 *跨記錄 / 跨檔案* 之
指針、索引、收據、慣例分離、內部取代指針五類。

判準（跑前寫死，不得事後增刪）
-----------------------------
P1 外部指針有效性：M1 docs/08 引用 M12 docs/02 之錨點，必須仍落在
   「非可分（束縛）族 … m = M0[Σσγ + λΣσ_iσ_j w_ij]」該行或其鄰行 (±3 行)。
P2 收件收據：register 內凡 target_pid=M12 且 kind=handed_off 之事件，
   其後必須存在一條 M12 之 note/progress 事件，note 內同時含「收據」與來源 pid。
P3 根索引同步：v11 根 README 之 M 系列表必須包含 register 內全部 v11.0/M* 之 pid，
   且狀態字串與 register 最新 state 相符（registered/partial/closed 三檔）。
P4 慣例分離：γ 恆等之名必須帶慣例。核實
   (a) /90 慣例配 45−n  ⟹ max|γ_n+γ_{45−n}−4| < 1e-25
   (b) /88 慣例配 44−n  ⟹ max|γ_n+γ_{44−n}−4| < 1e-25
   (c) 交叉配對 (/90 配 44−n ；/88 配 45−n) 必須 FAIL（> 1e-3）
   任一條不成立 ⟹ 恆等「帶慣例」之聲稱不成立，跨記錄引用會錯。
P5 內部取代指針：凡被後續輪次取代 / 撤回 / 判死之 doc，其檔頭須有
   標記 "⚠️ **跨文檔指針" 之行，列舉 (a) 取代者 (b) 被取代內容 (c) 以何者為準。
   指定集合 = {17, 24, 28, 33, 34}（逐項理由見 docs/40）。
P6 狀態快照：根 README "狀態快照" 段之 M 系列總數與三分組計數，
   須等於 register 實算值。
P7 記錄落碟：docs/40_propagation_audit.md 存在且非空（本輪產物必須落碟）。

不使用任何事後門檻；P1-P7 皆為布林，無 tolerance 可調（P4 除外：1e-25 / 1e-3
為 mpmath 60 dps 下之機器零與「明顯不同」分界，跑前寫死）。
"""
import json
import os
import re
import sys

import mpmath as mp

# ── canonical log 自寫（2026-10-29 修）────────────────────────────────────────
# 原缺陷：本腳本只 print 到 stdout，canonical log 靠 shell redirect 產生
#          ⟹ scripts/15 之 A3「跑前快照」法失效（腳本不覆寫 ⟹ 狀態一變即漂移）。
# 此處令 stdout 雙寫（終端 + canonical log），atexit 落碟。
import atexit as _atexit
import pathlib as _pl
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path = path
        self.buf = []
        self.out = _sys.__stdout__

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


_log = _Tee(_pl.Path(__file__).resolve().parent.parent / "data" / "run_48_propagation_audit.log")
_sys.stdout = _log
_atexit.register(_log.close)
# ────────────────────────────────────────────────────────────────────────────

mp.mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
M12 = os.path.dirname(HERE)
V11 = os.path.dirname(os.path.dirname(M12))
ORIG = os.path.dirname(V11)
REG = os.path.join(ORIG, "problem_register", "register.jsonl")

PASS = 0
FAIL = 0


def chk(name, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
        print("[PASS] %s%s" % (name, ("  " + detail) if detail else ""))
    else:
        FAIL += 1
        print("[FAIL] %s%s" % (name, ("  " + detail) if detail else ""))


def read(p):
    with open(p, encoding="utf-8") as f:
        return f.read()


def lines(p):
    return read(p).splitlines()


print("48 - M12 傳播審計（全鏈審計之後之下游影響）")
print("M12 = %s" % M12)
print("V11 = %s" % V11)

# ---------------------------------------------------------------- P1
print("\n[P1] 外部指針有效性：M1 docs/08 -> M12 docs/02")
m1 = os.path.join(V11, "main_problems", "m1_mass_file_schema",
                  "docs", "08_downstream_impact.md")
tgt = os.path.join(M12, "docs", "02_r1_domain_additivity.md")
m1_hits = []
for i, ln in enumerate(lines(m1), 1):
    if "docs/02_r1_domain_additivity.md" in ln and "M12" in ln:
        m1_hits.append(i)
cited = None
for i in m1_hits:
    for m in re.finditer(r"L(\d+)", lines(m1)[i - 1]):
        cited = int(m.group(1))
chk("P1a M1 確有引用 M12 docs/02", len(m1_hits) > 0, "hit lines=%s" % m1_hits)
body = lines(tgt)
carrier = [i for i, ln in enumerate(body, 1) if "Σσ_iσ_j w_ij" in ln or "m = M0[" in ln
           or "m = M₀[" in ln]
chk("P1b M12 docs/02 內載體行存在", len(carrier) > 0, "carrier lines=%s" % carrier)
if cited is not None and carrier:
    near = min(abs(cited - c) for c in carrier)
    chk("P1c 引用行號仍落載體 ±3 行內", near <= 3,
        "cited L%d, carrier L%s, dist=%d" % (cited, carrier, near))
else:
    chk("P1c 引用行號仍落載體 ±3 行內", False, "cited=%s carrier=%s" % (cited, carrier))

# ---------------------------------------------------------------- P2
print("\n[P2] 收件收據：target_pid=M12 之 handed_off 逐條")
ev = [json.loads(l) for l in read(REG).splitlines() if l.strip()]
incoming = [e for e in ev if e.get("target_pid") == "M12" and e.get("kind") == "handed_off"]
print("    incoming handed_off 條數 = %d" % len(incoming))
for e in incoming:
    src, ts = e.get("pid"), e.get("ts", "")
    # 來源 pid 必須係「完整 token」—— 否則 M1 會被 M12/M11 之文字誤配（假 PASS）
    src_re = re.compile(r"(?<![A-Za-z0-9])" + re.escape(src) + r"(?![0-9])")
    later = [x for x in ev
             if x.get("pid") == "M12" and x.get("ts", "") > ts
             and ("收據" in (x.get("note") or "") or "收件" in (x.get("note") or ""))
             and src_re.search(x.get("note") or "")]
    chk("P2 %s @%s" % (src, ts[:19]), len(later) > 0,
        "收據事件 %d 條" % len(later))

# ---------------------------------------------------------------- P3 / P6
print("\n[P3/P6] 根 README M 系列表 vs register 最新 state")
state = {}
for e in ev:
    if e.get("state") is not None and e.get("pid"):
        state[e["pid"]] = e["state"]
    if e.get("verdict") is not None and e.get("pid"):
        state[e["pid"]] = e["verdict"]
mids = sorted([p for p in state if re.match(r"^M\d+$", p)],
              key=lambda s: int(s[1:]))
print("    register M pid = %s" % mids)
readme = os.path.join(V11, "README.md")
rl = lines(readme)
table = {}
for ln in rl:
    m = re.match(r"^\|\s*M(\d+)\s*\|", ln)
    if m:
        table[int(m.group(1))] = ln
print("    根 README 表內 M = %s" % sorted(table))
for pid in mids:
    n = int(pid[1:])
    chk("P3 %s 在根 README M 表" % pid, n in table)
cn = {"registered": "已註冊", "partial": "部分閉合", "closed": "已閉合",
      "negative": "已閉合", "retracted": "撤回", "merged": "已合併"}
for pid in mids:
    n = int(pid[1:])
    if n not in table:
        continue
    want = cn.get(state[pid], state[pid])
    chk("P3b %s 狀態字串含「%s」" % (pid, want), want in table[n],
        "state=%s" % state[pid])

snap = {}
for ln in rl:
    if ln.startswith("M 系列："):
        m = re.search(r"(\d+)\s*個", ln)
        if m:
            snap["total"] = int(m.group(1))
    m = re.match(r"^\s+(registered|partial|closed)\s*:\s*(.+)$", ln)
    if m:
        items = []
        for s in m.group(2).split(","):
            mm = re.search(r"M(\d+)", s)
            if mm:
                items.append("M" + mm.group(1))
        snap[m.group(1)] = items
calc = {}
for pid in mids:
    k = {"registered": "registered", "partial": "partial", "closed": "closed",
         "negative": "closed"}.get(state[pid])
    calc.setdefault(k, []).append(pid)
print("    快照 = %s" % snap)
print("    實算 = %s" % calc)
chk("P6a 總數相符", snap.get("total") == len(mids),
    "快照 %s vs 實算 %d" % (snap.get("total"), len(mids)))
for k in ("registered", "partial", "closed"):
    got = sorted(snap.get(k, []), key=lambda s: int(s[1:]))
    want = calc.get(k, [])
    chk("P6b %s 分組相符" % k, got == want, "快照 %s vs 實算 %s" % (got, want))

# ---------------------------------------------------------------- P4
print("\n[P4] 慣例分離：同名 γ 恆等之配對索引")


def maxdev(N, pair):
    g = lambda n: 4 * mp.sin(n * mp.pi / N) ** 2
    return max(abs(g(n) + g(pair(n)) - 4) for n in range(1, 45))


d_a = maxdev(90, lambda n: 45 - n)
d_b = maxdev(88, lambda n: 44 - n)
d_c1 = maxdev(90, lambda n: 44 - n)
d_c2 = maxdev(88, lambda n: 45 - n)
print("    (a) /90 配 45-n : %s" % mp.nstr(d_a, 6))
print("    (b) /88 配 44-n : %s" % mp.nstr(d_b, 6))
print("    (c) 交叉 /90 配 44-n : %s ；/88 配 45-n : %s" % (mp.nstr(d_c1, 6), mp.nstr(d_c2, 6)))
chk("P4a /90 x 45-n 精確", d_a < mp.mpf("1e-25"))
chk("P4b /88 x 44-n 精確", d_b < mp.mpf("1e-25"))
chk("P4c 交叉配對必須失效", d_c1 > mp.mpf("1e-3") and d_c2 > mp.mpf("1e-3"),
    "cross dev 均 > 1e-3 ⟹ 配對不可互換")

# ---------------------------------------------------------------- P5
print("\n[P5] 內部取代指針（五份被取代 doc）")
SUPERSEDED = {
    17: "R8 遺留「帶號核仍存活但受阻」-> docs/25 已測死、docs/37 證不需帶號核",
    24: "G1「帶號核未測」-> docs/25 已測、docs/36 Z4 窮盡、docs/37/38 判死",
    28: "R15「真結論 = f-provenance」-> docs/30/31 已攻（R17/R18 全判死）",
    33: "R20 §E「(3,16,41) 一線之差」-> docs/34 改判死",
    34: "R21 §A / §B（41 為殘餘、60=|A5| 巧合）-> docs/35 撤回",
}
for n, why in sorted(SUPERSEDED.items()):
    p = os.path.join(M12, "docs", "%02d_%s.md" % (n, ""))
    cand = [f for f in os.listdir(os.path.join(M12, "docs"))
            if f.startswith("%02d_" % n)]
    ok = False
    det = "檔案缺失"
    for c in cand:
        txt = read(os.path.join(M12, "docs", c))
        ok = "⚠️ **跨文檔指針" in txt
        det = c
    chk("P5 docs/%02d" % n, ok, det)

# ---------------------------------------------------------------- P7
print("\n[P7] 本輪記錄落碟")
p40 = os.path.join(M12, "docs", "40_propagation_audit.md")
chk("P7 docs/40 存在且非空",
    os.path.exists(p40) and os.path.getsize(p40) > 1000,
    "%d bytes" % (os.path.getsize(p40) if os.path.exists(p40) else 0))

print("\nRESULT: %d/%d PASS" % (PASS, PASS + FAIL))
sys.exit(0 if FAIL == 0 else 1)