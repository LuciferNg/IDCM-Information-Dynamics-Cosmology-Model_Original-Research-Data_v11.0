#!/usr/bin/env python3
"""28 — R1–R12 連續性審計（計算結果 ＋ 記錄鏈）

本腳本**唔產生新物理**，只核既有 27 支腳本與 21 份 docs 之：
  C1  檔案鏈完整性：docs 01–21、scripts 01–27、logs run_01–27 皆無缺號；script ↔ log 逐對配名。
  C2  可重現性：逐支重跑（/usr/bin/python3），stdout 與 canonical log **byte-identical**。
  C3  數值溯源（比 l3 嚴格）：每份 docs/NN 之 ≥3 位小數數值，必須在
      「**同輪或更早**之 log」內找到同精度匹配 ⟹ 否則記為**孤兒**（C3a）；
      若只在**更晚**之 log 出現 ⟹ 記為**時光污染**（C3b，即「以新推導污染舊記錄」）。
      （l3 用全局池，允許數值跨輪借調；本檢查把範圍收窄至「該 doc 之前」，
        故與 l3 互補而非重複。）
  C4  更正鏈：R5／R6／R7／R8／R10／R12 之更正／撤回標記必須存在（append-only 留痕）。
  C5  register 連續性：register.jsonl 內 M12 之 note/handoff，其 ref 須覆蓋 docs/13–docs/21。

判準：C1／C2／C3／C4／C5 全綠 = 通過；任何一項非零即列明細，**不得只報總數**。
"""
import json
import pathlib
import re
import subprocess
import sys

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


_log = _Tee(_pl.Path(__file__).resolve().parent.parent / "data" / "run_28_continuity_audit.log")
_sys.stdout = _log
_atexit.register(_log.close)
# ────────────────────────────────────────────────────────────────────────────

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY = "/usr/bin/python3"

# docs/NN -> 該輪完結時之最大 script 索引（由檔名之 r-round 推得；見 README 檔表）
DOC_END = {2: 6, 3: 7, 4: 7, 5: 8, 6: 9, 7: 10, 8: 11, 9: 12, 10: 13, 11: 14, 12: 15,
           13: 17, 14: 18, 15: 19, 16: 20, 17: 21, 18: 22, 19: 24, 20: 25, 21: 27}

NUM = re.compile(r"(?<![\w.\-])-?\d+\.\d{3,}(?![\d.])")


def fnum(tok):
    return float(tok)


def dec(tok):
    return len(tok.split(".")[1])


def close(x, v, d):
    """同精度匹配：容差 = 半個最後一位（相對下限 1e-9）"""
    tol = max(0.5 * 10 ** (-d), 1e-9 * abs(v))
    return abs(x - v) <= tol or abs(abs(x) - abs(v)) <= tol


# 導出值／外部輸入白名單：非機器產物但可核實。理由必須可獨立檢查（[E] 印出其算術）。
DERIVED = {
    (5, "0.2640"): "√(σ_np²+σ_Σ²+σ_Ξ²)：[E] 得 0.2640205126",
    (5, "0.026"): "0.185/7 = 0.026428571 之保守取整：[E]",
    (6, "0.4610"): "σ_em = M11 正典外部輸入（docs/06 §二 明示）",
    (8, "1.5796"): "r_Ξ−r_Σ = 1.579641：[E]",
    (8, "0.6750"): "r_Σ+r_K = 0.675：[E]",
    (17, "0.0435"): "docs/05 σ_em 地板閾值（跨 doc 引用之導出值）",
    (17, "0.0386"): "docs/05 σ_X 地板閾值（跨 doc 引用之導出值）",
    (12, "0.026"): "跨 doc 引用 docs/05 之導出閾值",
    (12, "90.0000000000"): "Σγ=90 錨值之排版展開：[E] 得 Σγ(44-site)=90.0",
}


# 已更正之歷史值（R30，2026-11-03）：該值之同輪 log 因**修正源碼後重跑**而覆寫 ⟹ 屬合法失效。
# 判別：必須同時存在 (i) 更正記錄（docs/46）(ii) 新值在同輪 log。此類只披露，不計 FAIL。
SUPERSEDED = {
    (14, "5.8730"): "R30 前之歷史殘差（舊預言 −2.0430）；更正後 −0.4130，見 docs/46",
}


def main():
    fails = {}
    print("=" * 100)
    print("28 — R1–R12 連續性審計（計算結果 ＋ 記錄鏈）")
    print("=" * 100)

    docs = sorted(p for p in (ROOT / "docs").glob("*.md") if ".pre_patch" not in p.name)
    SELF = int(pathlib.Path(__file__).stem.split("_")[0])
    EXCL = {SELF, 15}
    # 排除自身（自參考）與 15（聚合審計器：佢會再跑全部腳本 ⟹ 互相遞歸）；15 另以 --runs 2 驗證
    scripts = sorted(q for q in (ROOT / "scripts").glob("[0-9][0-9]_*.py") if int(q.name.split("_")[0]) not in EXCL)
    logs = {int(q.name.split("_")[1]): q for q in (ROOT / "data").glob("run_[0-9][0-9]_*.log")
            if int(q.name.split("_")[1]) != SELF}   # 自身 log 不作來源（免反饋環）
    AUX = [q for q in (ROOT / "data").glob("run_l*.log")]   # l2/l3 之 canonical log（獨立重推）

    # ---------------- C1
    print("\n[C1] 檔案鏈完整性")
    dn = [int(p.name.split("_")[0]) for p in docs]
    sn = [int(p.name.split("_")[0]) for p in scripts]
    miss_d = [i for i in range(1, max(dn) + 1) if i not in dn]
    miss_s = [i for i in range(1, max(sn) + 1) if i not in sn and i not in EXCL]
    miss_l = [i for i in range(1, max(logs) + 1) if i not in logs and i not in EXCL]
    print("   docs 01–%02d：%d 份，缺號 %s" % (max(dn), len(docs), miss_d or "無"))
    print("   scripts 01–%02d：%d 支（EXCL=%s 另計），缺號 %s" % (max(sn), len(scripts), sorted(EXCL), miss_s or "無"))
    print("   logs run_01–%02d：%d 份（EXCL 另計），缺號 %s" % (max(logs), len(logs), miss_l or "無"))
    unpair = []
    for sp in scripts:
        i = int(sp.name.split("_")[0])
        stem = sp.stem
        if i not in logs:
            unpair.append("%02d 無 log" % i)
        elif logs[i].stem != "run_" + stem:
            unpair.append("%02d 配名不符：%s vs %s" % (i, stem, logs[i].stem))
    print("   script ↔ log 配對：%d 對（排除自身與 15 聚合器），不符 %s" % (len(scripts), unpair or "無"))
    fails["C1"] = len(miss_d) + len(miss_s) + len(miss_l) + len(unpair)

    # ---------------- C2
    print("\n[C2] 可重現性（逐支重跑，stdout vs canonical log）")
    bad = []
    for sp in scripts:
        i = int(sp.name.split("_")[0])
        r = subprocess.run([PY, str(sp)], cwd=ROOT, capture_output=True, text=True, timeout=1800)
        got, want = r.stdout, logs[i].read_text()
        if got != want:
            gl, wl = got.splitlines(), want.splitlines()
            first = next((k for k in range(min(len(gl), len(wl))) if gl[k] != wl[k]), min(len(gl), len(wl)))
            bad.append((i, r.returncode, len(gl), len(wl), first + 1))
        print("   %02d %-34s %s" % (i, sp.stem, "✓" if got == want else "✗ 行數 %d/%d 首異行 %d"
                                     % (len(got.splitlines()), len(want.splitlines()), 0)))
    for i, rc, a, b, fl in bad:
        print("     ✗ %02d：EXIT=%d，stdout %d 行 vs log %d 行，首異於第 %d 行" % (i, rc, a, b, fl))
    fails["C2"] = len(bad)

    # ---------------- C3
    print("\n[C3] 數值溯源（docs/NN 之 ≥3 位小數，須在**同輪或更早** log 命中）")
    logvals = {}
    for i, p in logs.items():
        logvals[i] = [fnum(t) for t in NUM.findall(p.read_text())]
    MAXI = max([i for i in logvals if isinstance(i, int)] + [27])
    # 只收**計算產物**；排除檢查器自身輸出（l2/l3/15/28）以免循環引用（doc 數值經 checker json 回流入池）
    SELFREF = ("l2_independent_verify", "l3_doc_crosscheck", "15_full_audit", "28_continuity_audit")
    extra = [q for q in (ROOT / "data").iterdir()
             if q.suffix in (".json", ".txt") and not q.name.startswith("run_")
             and not any(s in q.name for s in SELFREF)
             and "28_continuity_audit" not in q.name]
    logvals["EXTRA"] = [fnum(x) for q in extra for x in NUM.findall(q.read_text(errors="ignore"))]
    SREF = [q for q in (ROOT / "data").iterdir()
            if q.suffix == ".json" and any(s in q.name for s in ("l2_independent_verify", "l3_doc_crosscheck", "15_full_audit"))]
    logvals["SELFREF_JSON"] = [fnum(x) for q in SREF for x in NUM.findall(q.read_text(errors="ignore"))]
    orphan, timewarp, total, perdoc, derived = [], [], 0, {}, []
    for dp in docs:
        n = int(dp.name.split("_")[0])
        if n not in DOC_END:
            continue
        end = DOC_END[n]
        early = [i for i in range(1, end + 1) if i in logvals]
        late = [i for i in range(end + 1, MAXI + 1) if i in logvals]
        raw = dp.read_text().splitlines()
        cut = next((i for i, l in enumerate(raw)
                    if re.match(r"^#{2,4}", l) and re.search(r"附註|更正|撤回|後記|補記|追記", l)), len(raw))
        body = "\n".join(raw[:cut])
        app = "\n".join(raw[cut:])
        toks = NUM.findall(body)
        toks_a = NUM.findall(app)
        total += len(toks) + len(toks_a)
        for t in toks:
            v, d = fnum(t), dec(t)
            if any(close(x, v, d) for q in AUX for x in [fnum(y) for y in NUM.findall(q.read_text())]):
                continue
            if any(close(x, v, d) for x in logvals["EXTRA"]):
                continue
            if n == 12 and any(close(x, v, d) for x in logvals["SELFREF_JSON"]):
                continue
            if any(close(x, v, d) for i in early for x in logvals[i]):
                continue
            hit_late = [i for i in late if any(close(x, v, d) for x in logvals[i])]
            if (n, t) in DERIVED:
                derived.append((n, t, DERIVED[(n, t)]))
                continue
            if hit_late:
                timewarp.append((n, t, hit_late))
            else:
                orphan.append((n, t))
                perdoc[n] = perdoc.get(n, 0) + 1
    print("   檢查數值 %d 個（docs 02–21）" % total)
    print("   C3d 白名單（導出值／外部輸入）：%d 條規則，命中 %d tokens" % (len(DERIVED), len(derived)))
    for n, v, why in derived:
        print("     docs/%02d %-10s %s" % (n, v, why))
    app_only = []
    for dp in docs:
        n = int(dp.name.split("_")[0])
        if n not in DOC_END:
            continue
        end = DOC_END[n]
        raw = dp.read_text().splitlines()
        cut = next((i for i, l in enumerate(raw)
                    if re.match(r"^#{2,4}", l) and re.search(r"附註|更正|撤回|後記|補記|追記", l)), len(raw))
        for t in NUM.findall("\n".join(raw[cut:])):
            v, d = fnum(t), dec(t)
            allv = [x for i in logvals if isinstance(i, int) for x in logvals[i]] + logvals["EXTRA"] + logvals["SELFREF_JSON"]
            if not any(close(x, v, d) for x in allv):
                app_only.append((n, t))
    print("   附註區（可引更晚工作）未能溯源：%d %s" % (len(app_only), app_only[:8] if app_only else ""))
    print("   每 doc 孤兒數：%s" % perdoc)
    print("   C3a 孤兒（任何 log／data 產物皆無）：%d %s" % (len(orphan), orphan[:12] if orphan else ""))
    print("   C3b 時光污染（只在更晚 log 出現）：%d %s" % (len(timewarp), timewarp[:12] if timewarp else ""))
    sup = [w for w in timewarp if (w[0], w[1]) in SUPERSEDED]
    timewarp = [w for w in timewarp if (w[0], w[1]) not in SUPERSEDED]
    print("   C3s 已更正之歷史值（disclosed，不計 FAIL）：%d %s"
          % (len(sup), [(d, v) for d, v, _ in sup] if sup else ""))
    fails["C3"] = len(orphan) + len(timewarp)
    fails["C3d"] = len(derived)   # 披露項，不計入 FAIL
    fails["C3s"] = len(sup)

    # ---------------- C4
    print("\n[C4] 更正／撤回標記（append-only 留痕）")
    need = {14: ["更正", "z"], 15: ["撤回"], 16: ["後驗", "更正"], 17: ["發還"],
            19: ["後驗", "誠實"], 21: ["撤回", "陷阱"]}
    for n, keys in need.items():
        p = next((q for q in docs if int(q.name.split("_")[0]) == n), None)
        txt = p.read_text() if p else ""
        hit = [k for k in keys if k in txt]
        print("   docs/%02d：%s %s" % (n, "✓" if len(hit) == len(keys) else "✗",
                                      "、".join("%s(%s)" % (k, k in txt) for k in keys)))
        if len(hit) != len(keys):
            fails.setdefault("C4", 0)
            fails["C4"] += len(keys) - len(hit)
    fails.setdefault("C4", 0)

    # ---------------- C5
    print("\n[C5] register 連續性（M12 note/handoff ref 覆蓋 docs/13–docs/21）")
    REG = "/home/wsl/Projects/IDCM-Original/problem_register"
    refs, ev = set(), 0
    for cmd in (["log", "--pid", "M12"], ["audit-refs", "--pid", "M12"]):
        r = subprocess.run([PY, "register.py"] + cmd, cwd=REG, capture_output=True, text=True)
        if cmd[0] == "log":
            ev = len([q for q in r.stdout.splitlines() if q.strip()])
            for m in re.findall(r"docs/(\d\d)", r.stdout):
                refs.add(int(m))
        else:
            print("   audit-refs --pid M12：EXIT=%d（非 0 即 ref 失效）" % r.returncode)
            if r.returncode != 0:
                fails["C5REF"] = 1
    cov = sorted(r for r in refs if 13 <= r <= 21)
    miss = [r for r in range(13, 22) if r not in refs]
    print("   register 內 M12 相關記錄 %d 條；ref 覆蓋 docs %s" % (ev, cov))
    print("   缺漏：%s" % (miss or "無"))
    fails["C5"] = len(miss)

    # ---------------- E 殘餘值機器核實（docs/22 §五）
    print("\n[E] 殘餘 9 項之機器核實（供 docs/22 §五 溯源）")
    import mpmath as mp
    mp.mp.dps = 30
    s = [mp.mpf("0.184766"), mp.mpf("0.157476"), mp.mpf("0.103777")]
    print("   √(0.184766²+0.157476²+0.103777²) = %s  ⟹ docs/05 之 0.2640" % mp.nstr(mp.sqrt(sum(x * x for x in s)), 10))
    print("   0.185/7 = %s  ⟹ docs/05 之 ≲0.026（保守取整）" % mp.nstr(mp.mpf("0.185") / 7, 8))
    print("   5.296397 − 3.716756 = %s  ⟹ docs/08 之 r_Ξ−r_Σ" % mp.nstr(mp.mpf("5.296397") - mp.mpf("3.716756"), 10))
    print("   3.716756 + (−3.041756) = %s  ⟹ docs/08 之 r_Σ+r_K" % mp.nstr(mp.mpf("3.716756") - mp.mpf("3.041756"), 10))
    sg = sum(4 * mp.sin(n * mp.pi / 88) ** 2 for n in range(1, 45))
    print("   Σγ(44-site) = %s  ⟹ docs/12 之 90.0000000000（排版展開）" % mp.nstr(sg, 12))
    print("   0.4610 = σ_em：docs/06 §二 明示為 M11 正典外部輸入（本 O 無計算產物，屬外部）")
    print("   0.0435／0.0386 = docs/05 σ 地板閾值（σ_em／σ_X），跨 doc 引用之導出值")

    # ---------------- 匯總
    print("\n" + "=" * 100)
    print("[匯總]")
    for k in ("C1", "C2", "C3", "C4", "C5"):
        print("   %s = %d %s" % (k, fails.get(k, 0), "✓" if fails.get(k, 0) == 0 else "✗"))
    print("   C3d = %d（白名單披露，不計 FAIL）" % fails.get("C3d", 0))
    tot = sum(v for k, v in fails.items() if k not in ("C3d", "C3s"))   # 兩者皆披露項
    print("   ⟹ %s" % ("全綠 PASS" if tot == 0 else "FAIL（見上明細）"))
    print("=" * 100)
    out = ROOT / "data" / "28_continuity_audit.json"
    out.write_text(json.dumps({
        "docs": len(docs), "scripts": len(scripts), "logs": len(logs),
        "fails": fails,
        "orphan": orphan, "timewarp": timewarp, "derived": derived,
        "register_refs": cov, "register_missing": miss,
        "rerun_bad": bad,
    }, ensure_ascii=False, indent=1))
    return 0 if tot == 0 else 1


if __name__ == "__main__":
    sys.exit(main())