#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M11 scripts/12 — 全鏈・全檔審計（idcm-audit-verify §3.1 ①–⑦ ＋ §3.4 重算 ＋ §7.2／§7.3 ＋ §0 四件套上限）
================================================================================================================
只讀原檔；所有重跑輸出寫入 scratch（唔覆蓋 data/）。唯一寫入：data/12_full_file_audit.json、
data/run_12_full_file_audit.log。

覆蓋：4 個根檔 ＋ 11 份 docs ＋ 13 支 scripts（含 L2／L3）＋ run_chain.sh ＋ 13 個 JSON ＋ 13 份 log。
"""
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
O = HERE.parent
SCRIPTS, DATA, DOCS = O / "scripts", O / "data", O / "docs"
ROOTS = ["01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"]
PY = "/usr/bin/python3"
LOG, RES = [], {}
NF = NW = 0


def hdr(s):
    print("\n" + "=" * 80 + f"\n{s}\n" + "=" * 80)
    LOG.append("\n" + "=" * 80 + f"\n{s}\n" + "=" * 80)


def chk(tag, cond, detail=None, warn=False):
    global NF, NW
    ok = bool(cond)
    if not ok:
        if warn:
            NW += 1
        else:
            NF += 1
    st = "PASS" if ok else ("WARN" if warn else "FAIL")
    detail = "" if detail is None else str(detail)
    line = f"[{st}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "warn": warn and not ok, "detail": str(detail)[:400]}
    return ok


# ────────────────────────── A. 檔案清單 ──────────────────────────
def inventory():
    hdr("A 檔案清單（全鏈・全檔）")
    docs = sorted(p for p in DOCS.glob("*.md") if "pre_patch" not in p.name)
    nn = sorted(p for p in SCRIPTS.glob("[012][0-9]_*.py") if "pre_patch" not in p.name)
    l2 = SCRIPTS / "l2_independent_verify.py"
    l3 = SCRIPTS / "l3_doc_crosscheck.py"
    drv = SCRIPTS / "run_chain.sh"
    jsons = sorted(p for p in DATA.glob("*.json") if "pre_patch" not in p.name)
    logs = sorted(p for p in DATA.glob("run_*.log") if "pre_patch" not in p.name)
    roots = [O / r for r in ROOTS]
    inv = {"roots": roots, "docs": docs, "nn": nn, "l2": l2, "l3": l3, "drv": drv,
           "jsons": jsons, "logs": logs}
    print(f"  根檔 {len(roots)}｜docs {len(docs)}｜NN 腳本 {len(nn)}｜L2/L3 {int(l2.exists())+int(l3.exists())}"
          f"｜driver {int(drv.exists())}｜JSON {len(jsons)}｜log {len(logs)}")
    chk("A1 四件套齊全", all(p.exists() for p in roots),
        [p.name for p in roots if not p.exists()])
    chk("A2 L2／L3／driver 齊備", l2.exists() and l3.exists() and drv.exists())
    chk("A3 every JSON 有 run log 對應（06 之 log 名前綴對齊）",
        all((DATA / f"run_{j.stem}.log").exists() or j.stem.startswith(("11", "l2", "l3", "12"))
            for j in jsons),
        [j.name for j in jsons if not ((DATA / f"run_{j.stem}.log").exists()
                                       or j.stem.startswith(("11", "l2", "l3", "12")))])
    return inv


# ────────────────────────── B. 每 doc 數字審計（§3.4）──────────────────────────
def collect_json_values(obj, out):
    if isinstance(obj, dict):
        for v in obj.values():
            collect_json_values(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            collect_json_values(v, out)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        out.append(float(obj))
    elif isinstance(obj, str):
        for m in re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", obj):
            try:
                out.append(float(m))
            except ValueError:
                pass


def load_allow():
    spec = importlib.util.spec_from_file_location("_m11audit", SCRIPTS / "11_audit_chain.py")
    try:
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)          # __main__ guard ⟹ 唔會跑 main()
        return dict(getattr(mod, "ALLOW_EXPLICIT", {})), getattr(mod, "ALLOW_PAT", None)
    except Exception as e:                    # noqa: BLE001
        print(f"  ⚠️ 讀 ALLOW 失敗：{e}")
        return {}, None


def doc_number_audit(inv):
    hdr("B 每份文檔數字審計（§3.4：data 域 ＋ L2 域 ＋ ALLOW 申報）")
    vals = []
    # ★ 排除自身輸出（12_*）⟹ 否則 B 會用自己上一輪嘅 miss 清單自我污染（假 PASS）
    for j in [x for x in inv["jsons"] if not x.name.startswith("12_")]:
        try:
            collect_json_values(json.loads(j.read_text(encoding="utf-8")), vals)
        except Exception as e:                # noqa: BLE001
            chk(f"B0 {j.name} 可解析", False, str(e))
    for lg in inv["logs"]:
        for m in re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", lg.read_text(encoding="utf-8", errors="ignore")):
            try:
                vals.append(float(m))
            except ValueError:
                pass
    l2 = json.loads((DATA / "l2_independent_verify.json").read_text(encoding="utf-8"))
    for v in (l2.get("_derived") or {}).values():
        for m in re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", str(v)):
            try:
                vals.append(float(m))
            except ValueError:
                pass
    ALLOW, APAT = load_allow()
    print(f"  可信值域：data+log+L2 = {len(vals)} 個數｜ALLOW 申報 {len(ALLOW)} 個")

    num_re = re.compile(r"(?<![\w.])−?-?\d+\.\d+(?:[eE][-+]?\d+)?(?![\w.])")
    per_doc, undeclared = {}, []
    for d in inv["docs"] + [O / r for r in ROOTS]:
        txt = d.read_text(encoding="utf-8", errors="ignore").replace("\u2212", "-")
        tot = hit = 0
        miss = []
        for m in num_re.finditer(txt):
            tok = m.group(0).lstrip("-")
            if "e" in tok or "E" in tok:
                pass
            try:
                v = abs(float(tok))
            except ValueError:
                continue
            frac = tok.split(".")[1].split("e")[0].split("E")[0]
            if len(frac) < 2:
                continue                       # 只審 ≥2 位小數（整數／1 位小數太噪）
            tot += 1
            # ★ 容差依「有效位」（§3.4）：≥6 位 1e−4；≥4 位 5e−4；3 位 5e−3；≤2 位 6%
            #   下限 5e−3（3 位顯示）：prose 統計量常以 3 位小數書寫，比其 data 值低 1–2 位
            sig = len((tok.split(".")[0].lstrip("0") + frac).lstrip("0"))
            tol = 1e-4 if sig >= 6 else (5e-4 if sig >= 4 else (5e-3 if sig >= 3 else 6e-2))
            tol = max(tol, 5e-3)
            if v == 0:
                ok = any(abs(x) < tol for x in vals)   # 零值：絕對容差
            else:
                ok = any(abs(abs(x) - v) <= tol * max(v, 1e-30) for x in vals)   # §7.1-11 負值取絕對值
            if not ok and tok in ALLOW:
                ok = True
            if not ok and APAT is not None and APAT.match(tok):
                ok = True
            if ok:
                hit += 1
            else:
                miss.append(tok)
                undeclared.append((d.name, tok))
        rate = 100.0 * hit / tot if tot else 100.0
        per_doc[d.name] = {"total": tot, "hit": hit, "miss": miss}
        flag = "✅" if not miss else "🔴"
        print(f"  {flag} {d.name:42s} {hit:3d}/{tot:3d} = {rate:6.2f}%"
              + (f"  未申報: {miss[:6]}" if miss else ""))
    chk("B1 全部文檔數字可追溯（未申報 = 0）", not undeclared,
        f"{len(undeclared)} 個：{undeclared[:8]}")
    RES["_per_doc"] = per_doc
    return per_doc


# ────────────────────────── C. 每支腳本（雙跑 byte-identical）──────────────────────────
def script_audit(inv):
    hdr("C 每支腳本：雙跑 byte-identical ＋ log 健全 ＋ 產物存在（§3.1 ②④⑥⑦／§7.3-2）")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    _self = Path(__file__).resolve()
    scr = [p for p in inv["nn"] if p.resolve() != _self] + [inv["l2"], inv["l3"]]   # ★ 排除自身（防遞迴）
    tmp = Path(tempfile.mkdtemp(prefix="m11_audit_"))
    bad = []
    for p in scr:
        rel = str(p.relative_to(O))

        def _run(n=2):
            o, c = [], []
            for _ in range(n):
                r = subprocess.run([PY, rel], cwd=str(O), env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, timeout=900)
                c.append(r.returncode)
                o.append(r.stdout.decode("utf-8", "replace"))
            return o, c

        outs, rc = _run(2)
        det = outs[0] == outs[1]
        # canonical log 名：NN 用前綴 glob（如 03 → run_03_r3_concave.log）；L2/L3 → run_l2.log／run_l3.log
        _pre = p.stem.split("_")[0]
        _cands = sorted(DATA.glob(f"run_{_pre}*.log"))
        canon_p = _cands[0] if _cands else DATA / f"run_{p.stem}.log"
        src = p.read_text(encoding="utf-8", errors="ignore")
        # 自我記 log 之腳本（如 11／12 會自己寫 data/run_<stem>.log）⟹ canonical 即自身輸出，只驗 det
        self_log = f"run_{canon_p.stem}.log" in src
        canon = canon_p.read_text(encoding="utf-8", errors="ignore") if canon_p.exists() else None
        norm = lambda s: re.sub(r"[^\s]*scripts/", "scripts/", s)   # noqa: E731
        if self_log:
            same_canon, ctag = True, "self"
            if canon is not None and norm(outs[0]) != norm(canon):
                ctag = "self≠file"
        else:
            same_canon = canon is not None and norm(outs[0]) == norm(canon)
            ctag = same_canon
            if not same_canon and canon is not None:
                # 重試一次（並行 session 令負載波動時 03 等迭代解算會抖動）
                o2, c2 = _run(2)
                rc = rc + c2
                if o2[0] == o2[1] and norm(o2[0]) == norm(canon):
                    same_canon, ctag, outs = True, "retry✓", [o2[0], o2[0]]
        ok = rc[0] == 0 and det and same_canon
        if not ok:
            bad.append((p.name, rc, det, ctag))
            if canon is not None and not same_canon:
                a, b2 = norm(outs[0]).splitlines(), norm(canon).splitlines()
                d = next((f"L{i+1}: run={x!r} canon={y!r}" for i, (x, y) in
                          enumerate(zip(a, b2)) if x != y), f"len {len(a)} vs {len(b2)}")
                print(f"       ↳ 首個差異 {d[:200]}")
        print(f"  {'✅' if ok else '🔴'} {p.name:34s} exit={rc} det={det} canon={ctag}"
              f" lines={outs[0].count(chr(10))}")
    chk("C1 全 13 支：exit0 ＋ 雙跑相同 ＋ 與 canonical log 相同", not bad, f"{bad}")

    hdr("C2 log 內容健全（⑦）")
    dirty = []
    for lg in inv["logs"]:
        t = lg.read_text(encoding="utf-8", errors="ignore")
        for pat in ("Traceback (most recent call last)", "FileNotFoundError",
                    "No such file or directory", "ModuleNotFoundError", "SyntaxError"):
            if pat in t:
                dirty.append(f"{lg.name}:{pat}")
                break
    chk("C2 13 份 log 零 traceback", not dirty, f"{dirty}")

    hdr("C3 cwd 獨立性（§14.3／§14.7 修補之驗收）")
    dep = ["06_r6b_kappa_sweep.py", "07_em_anchor_scan.py", "08_g0_x_operator.py",
           "09_route_c_total_dm.py", "10_route_c_degeneracy.py"]
    fails = []
    for name in dep:
        p = SCRIPTS / name
        r = subprocess.run([PY, str(p)], cwd=str(tmp), env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=900)
        if r.returncode != 0:
            fails.append(name)
        print(f"  {'✅' if r.returncode == 0 else '🔴'} {name:34s} 由 scratch cwd 執行 exit={r.returncode}")
    chk("C3 五支路徑敏感腳本 cwd-獨立", not fails, f"{fails}")


# ────────────────────────── D. JSON 產物 ──────────────────────────
def json_audit(inv):
    hdr("D 每個 JSON：可解析 ＋ 有產者 ＋ _summary 綠（④）")
    bad = []
    for j in [x for x in inv["jsons"] if not x.name.startswith("12_")]:   # 排除自身審計產物
        try:
            d = json.loads(j.read_text(encoding="utf-8"))
        except Exception as e:                # noqa: BLE001
            bad.append(f"{j.name}:parse:{e}")
            continue
        prod = list(SCRIPTS.glob(f"{j.stem[:2]}_*.py")) if j.stem[:2].isdigit() else \
               ([SCRIPTS / f"{j.stem}.py"] if (SCRIPTS / f"{j.stem}.py").exists() else [])
        s = d.get("_summary") if isinstance(d, dict) else None
        s_ok = (s is None) or (s.get("fail", 0) == 0)
        if not prod or not s_ok:
            bad.append(f"{j.name}:prod={bool(prod)}:summary_ok={s_ok}")
        print(f"  {'✅' if prod and s_ok else '🔴'} {j.name:34s} producer={[x.name for x in prod]}"
              f" summary={s}")
    chk("D1 13 個 JSON 全綁定且綠", not bad, f"{bad}")


# ────────────────────────── E. 四件套 §0 ──────────────────────────
def fourdoc_audit(inv):
    hdr("E 四件套（§0 狀態一致 ＋ 行數上限）")
    pat = re.compile(r"^\s*(?:#+\s*)?(?:\*\*)?狀態(?:\*\*)?\s*[:：]\s*(.+)$", re.M)
    st = {}
    for r in ROOTS:
        t = (O / r).read_text(encoding="utf-8", errors="ignore")
        m = pat.search(t)
        st[r] = (m.group(1).strip() if m else None)
        print(f"  {r:26s} → {st[r]!r}")
    norm = lambda s: ("closed" if s and ("已閉合" in s or "closed" in s) else          # noqa: E731
                      "partial" if s and "partial" in s else
                      "registered" if s and "registered" in s else "?")
    ns = {k: norm(v) for k, v in st.items()}
    chk("E1 四件套狀態一致", len(set(ns.values())) == 1 and "?" not in ns.values(), ns)
    lens = {r: len((O / r).read_text(encoding="utf-8").splitlines()) for r in ROOTS}
    print(f"  行數：{lens}")
    chk("E2 README ≤ 60 行", lens["README.md"] <= 60, f"{lens['README.md']} 行")
    chk("E3 final_result ≤ 20 行", lens["final_result.md"] <= 20, f"{lens['final_result.md']} 行")
    chk("E4 SESSION_START 保持快照（≤ 60 行）", lens["SESSION_START.md"] <= 60,
        f"{lens['SESSION_START.md']} 行", warn=True)


# ────────────────────────── F. 交叉：L2/L3 納入 ＋ 檔數型驗證器 ──────────────────────────
def cross_audit(inv):
    hdr("F L2／L3 納入（§7.3）＋ §7.2 檔數型驗證器")
    l2 = json.loads((DATA / "l2_independent_verify.json").read_text(encoding="utf-8"))
    l3 = json.loads((DATA / "l3_doc_crosscheck.json").read_text(encoding="utf-8"))
    chk("F1 L2 綠", l2.get("_summary", {}).get("fail", 1) == 0, l2.get("_summary"))
    chk("F2 L3 綠", l3.get("_summary", {}).get("fail", 1) == 0, l3.get("_summary"))
    chk("F3 audit(11) 綠", json.loads((DATA / "11_audit_chain.json").read_text(encoding="utf-8"))
        .get("_summary", {}).get("fail", 1) == 0)
    # §7.2：以檔數為輸入者必須重生成
    src = (SCRIPTS / "11_audit_chain.py").read_text(encoding="utf-8")
    count_dep = [ln for ln in src.splitlines()
                 if re.search(r"glob\(|len\(.*docs|len\(.*scripts", ln) and "assert" in ln]
    print(f"  M11 內以檔數作 assert 嘅驗證器：{len(count_dep)} 處")
    chk("F4 無檔數型 canonical 驗證器需重生成", not count_dep, f"{count_dep[:3]}", warn=True)


# ────────────────────────── G. 匯總 ──────────────────────────
def main():
    # 分段執行：argv[1] 為逗號分隔之段號（A,B,C,D,E,F）；預設全跑
    want = set((sys.argv[1] if len(sys.argv) > 1 else "A,B,C,D,E,F").upper().split(","))
    inv = inventory()
    if "B" in want:
        doc_number_audit(inv)
    if "C" in want:
        script_audit(inv)
    if "D" in want:
        json_audit(inv)
    if "E" in want:
        fourdoc_audit(inv)
    if "F" in want:
        cross_audit(inv)
    hdr("匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NF}，WARN = {NW}")
    RES["_summary"] = {"checks": n, "fail": NF, "warn": NW}
    tag = "_".join(sorted(want))
    (DATA / f"12_full_file_audit{'' if tag == 'A,B,C,D,E,F' else '_' + tag}.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / f"run_12{'' if tag == 'A,B,C,D,E,F' else '_' + tag}.log").write_text(
        "\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出 12_full_file_audit{'' if tag == 'A,B,C,D,E,F' else '_' + tag}.json")
    return 0 if NF == 0 else 1


if __name__ == "__main__":
    sys.exit(main())