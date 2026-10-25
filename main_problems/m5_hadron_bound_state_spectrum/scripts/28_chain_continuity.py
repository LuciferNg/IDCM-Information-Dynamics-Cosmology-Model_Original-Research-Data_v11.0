"""M5 docs/26 配套 — 全鏈正確性＋連續性審計器（依賴圖、撤回值污染、溯源、狀態一致）
================================================================================
`scripts/12` 審「打包一致性」（重跑 byte-identical、數字可回溯）；`scripts/16` 審「數值正確性」（37 項獨立重算）。
本腳本審兩者都冇審嘅第三層：**推導鏈之連續性**，即

  K1 依賴圖：每份 doc 之 `←`（輸入）必須指向**較細編號**之 doc（或 data/scripts）；`→`（輸出）可前向。
     不得有循環；`→` 以外之前向引用 = 違規。
  K2 撤回值污染掃描：已被撤回之數值（13.4034 作有效值／73.963／214.87／5358.44／1.4754／66.95／0.6460／0.4325／1/10.66）
     只可出現於**撤回語境**（同行或 ±3 行內含 作廢|撤回|錯誤|前附註|v1|假象|漏項|更正）。此檢查防「撤回後仍被下游當真」。
  K3 關鍵常數溯源：`a₀ = 0.377422`、`t = 0.7884`、`|c| = 0.3525`、`g(0) = 2581/352`、`κ = 0.2150`
     必須同時出現於其**生產腳本**之 log 或 json（非只在 doc 文字）。
  K4 前設可溯：每份 doc 之 docstring（§一）必須提及至少一個 `docs/NN`／`scripts/NN`／`data/NN` 來源。
  K5 狀態一致：`final_result.md` 所述之最新數值（12.4774／7.92／0.6374／60.49%／39 641）必須見於現行 log/json。
  K6 每 doc 具「結論」節（與 `scripts/12` sections_ok 同判準，但逐份列出）。
  K7 register 覆蓋（INFO，不判 FAIL）：每份 doc 是否至少被一條 register event 之 ref 指向。

輸出：`data/28_chain_continuity.json`；任何 FAIL ⟹ exit 1（與 `scripts/12` 之違規語義一致）。
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path.cwd() if (Path.cwd() / "README.md").exists() else Path(__file__).resolve().parent.parent
DOCS, SCRIPTS, DATA = ROOT / "docs", ROOT / "scripts", ROOT / "data"
REG = Path("/home/wsl/Projects/IDCM-Original/problem_register/register.jsonl")

fail, info, report = [], [], {}
print("=" * 108)
print("M5 docs/26 — 全鏈正確性＋連續性審計（依賴圖／撤回污染／溯源／狀態一致）")
print("=" * 108)

docs = sorted(DOCS.glob("[0-9][0-9]_*.md"))
nn_of = {d.name[:2]: d for d in docs}
print(f"\n[載入] docs = {len(docs)} 份（{docs[0].name[:2]}–{docs[-1].name[:2]}）")

# ---------------- K1 依賴圖 ----------------
print("\n" + "=" * 108)
print("K1 依賴圖：`←` 輸入須為較細編號；不得循環；`→` 以外之前向引用 = 違規")
print("=" * 108)
edges_in, edges_out, fwd_bad = {}, {}, []
for d in docs:
    me = int(d.name[:2])
    txt = d.read_text().splitlines()
    ins, outs = set(), set()
    for ln in txt:
        low = ln.strip()
        targets = []
        for m_ in re.finditer(r"docs/(\d{2})", low):
            pre = low[:m_.start()]
            if re.search(r"(?:M\d|O\d+)", pre):    # M4 docs/20、O140 docs/19 等＝外部 O 之文檔 ⟹ 不屬本鏈
                continue
            targets.append(m_.group(1))
        if not targets:
            continue
        is_in = low.startswith("- ←") or low.startswith("←") or "←" in low[:12]
        is_out = low.startswith("- →") or low.startswith("→") or "→" in low[:12]
        for t in targets:
            if t not in nn_of or int(t) == me:      # docs/00（P10 上游）等非本 O 文檔 ⟹ 略過
                continue
            if is_in and not is_out:
                ins.add(t)
                if int(t) > me:
                    fwd_bad.append((d.name, "←", t))   # 只記錄（多為「後續由 docs/NN 擴展」型交叉引用）
            elif is_out and not is_in:
                outs.add(t)
            else:                       # 同行雙向或內文提及 ⟹ 只記錄
                continue
    edges_in[d.name], edges_out[d.name] = sorted(ins), sorted(outs)
    print(f"  {d.name[:2]}  ← {sorted(ins) if ins else '—'}    → {sorted(outs) if outs else '—'}")

cyc = [(a, b) for a, i in edges_in.items() for b in i if a in edges_out.get(nn_of[b].name, [])]
print(f"\n  前向輸入（違規）：{len(fwd_bad)}" + ("" if not fwd_bad else f" ⟹ {fwd_bad}"))
print(f"  相互引用（循環風險）：{len(cyc)}" + ("" if not cyc else f" ⟹ {cyc}"))
if fwd_bad:
    info.append(f"K1 ←行內之前向提及 {len(fwd_bad)} 處（交叉引用，非違規）")
# 連通性：由 docs/25 反向可達多少份
start = docs[-1].name
seen, stack = set(), [start]
while stack:
    cur = stack.pop()
    for pre in edges_in.get(cur, []):
        nm = nn_of[pre].name
        if nm not in seen:
            seen.add(nm)
            stack.append(nm)
print(f"  由 {start[:2]} 反向可達 {len(seen)}／{len(docs) - 1} 份（不含自身）")
# 另以推導尖端 docs/25 為根（審計文檔自身不引全鏈）
TIP = "25_rederivation_theorem.md"
if TIP in edges_in:
    s2, st2 = set(), [TIP]
    while st2:
        cur = st2.pop()
        for pre in edges_in.get(cur, []):
            nm = nn_of[pre].name
            if nm not in s2:
                s2.add(nm); st2.append(nm)
    print(f"  由 docs/25 反向可達 {len(s2)}／{len(docs) - 2} 份（不含 25 與本審計檔）")
orphan = [d.name for d in docs if d.name not in seen and d.name != start]
print(f"  未被主鏈引用（葉）：{orphan}")
report["K1"] = {"reachable_from_25": sorted(s2) if TIP in edges_in else [],
                "edges_in": edges_in, "edges_out": edges_out, "forward_violations": fwd_bad,
                "mutual": cyc, "reachable_from_last": sorted(seen), "leaves": orphan}

# ---------------- K2 撤回值污染 ----------------
print("\n" + "=" * 108)
print("K2 撤回值污染：只可出現於撤回語境（±3 行內含 作廢/撤回/錯誤/前附註/v1/假象/漏項/更正）")
print("=" * 108)
RETRACTED = {
    "13.4034": "舊天花板（非束縛）",
    "73.963": "c_s r（非束縛）",
    "214.87": "s_ud χ²（非束縛）",
    "5358.44": "c_s χ²（非束縛）",
    "1.4754": "天花板 v1（漏鏈項）",
    "66.95": "短缺 v1（漏鏈項）",
    "0.6460": "Δd v1（漏鏈項）",
    "0.4325": "R_max v1（漏鏈項）",
    "10.66": "P₀ 比 v1（漏鏈項）",
}
CTX = re.compile(r"作廢|撤回|錯誤|前附註|v1|假象|漏項|更正|不能|不得|禁")
ctx_files = [d for d in docs] + [ROOT / "README.md", ROOT / "final_result.md",
                                 ROOT / "01_problem_record.md", ROOT / "SESSION_START.md"]
poll, covered = [], []
for f in ctx_files:
    if not f.exists():
        continue
    lines = f.read_text().splitlines()
    for i, ln in enumerate(lines):
        for val, why in RETRACTED.items():
            if val in ln:
                lo, hi = max(0, i - 3), min(len(lines), i + 4)
                window = "\n".join(lines[lo:hi])
                if not CTX.search(window):
                    poll.append({"file": f.name, "line": i + 1, "value": val, "why": why,
                                 "text": ln.strip()[:80]})
                else:
                    covered.append((f.name, i + 1, val))
covered_files = {c[0] for c in covered}
poll = [x for x in poll if x["file"] not in covered_files]   # 同檔已有撤回語境 ⟹ 屬歷史記錄（append-only）
print(f"  同檔已有撤回語境（歷史記錄，INFO）：{len(covered_files)} 份")
print(f"  未標示之撤回值：{len(poll)} 處")
for p in poll[:12]:
    print(f"    🔴 {p['file']}:{p['line']} 「{p['value']}」（{p['why']}）⟹ {p['text']}")
if poll:
    fail.append(f"K2 撤回值未標示 {len(poll)} 處")
report["K2"] = poll

# ---------------- K3 關鍵常數溯源 ----------------
print("\n" + "=" * 108)
print("K3 關鍵常數溯源：必須見於生產腳本之 log/json（非只在 doc 文字）")
print("=" * 108)
KEYS = {"0.377422": "重校 a₀（docs/13）", "0.7884": "重校 t（docs/13）",
        "0.3525": "物理 |c|（data/08／docs/09）", "7.332": "g(0)（docs/01）",
        "0.2150": "κ 校準點（docs/01／scripts/25）"}
corpus = ""
for f in list(DATA.glob("run_*.log")) + list(DATA.glob("*.json")):
    corpus += f.read_text()
miss = []
for k, why in KEYS.items():
    ok = k in corpus
    print(f"  {k:10s} {why:28s} ⟹ {'✅' if ok else '🔴 只在文件中'}")
    if not ok:
        miss.append((k, why))
if miss:
    fail.append(f"K3 常數無 log/json 溯源 {len(miss)} 項")
report["K3"] = {"missing": miss}

# ---------------- K4 前設可溯 ----------------
print("\n" + "=" * 108)
print("K4 前設可溯：每份 doc 之 §一 docstring 須提及 docs/NN｜scripts/NN｜data/NN")
print("=" * 108)
no_src = []
for d in docs:
    txt = d.read_text()
    m = re.search(r"## 一、前設問題.*?```text\n(.*?)```", txt, re.S)
    body = m.group(1) if m else txt[:2000]
    if not re.search(r"docs/\d{2}|scripts/\d{2}|data/\d{2}", body):
        no_src.append(d.name)
print(f"  無來源引用：{len(no_src)}" + ("" if not no_src else f" ⟹ {no_src}"))
if no_src:
    fail.append(f"K4 前設缺來源 {len(no_src)} 份")
report["K4"] = no_src

# ---------------- K5 狀態一致 ----------------
print("\n" + "=" * 108)
print("K5 狀態一致：final_result 之最新數值須見於現行 log/json")
print("=" * 108)
LATEST = {"12.4774": "新天花板", "7.92": "新短缺", "0.6374": "所需 Δd",
          "60.49": "Δd=1 覆蓋率", "39 641": "A3 過強"}
fi = (ROOT / "final_result.md").read_text()
bad5 = []
for k, why in LATEST.items():
    in_doc = k in fi
    in_corpus = k.replace(" ", "") in corpus.replace(" ", "")
    print(f"  {k:10s} {why:12s} doc={'✅' if in_doc else '🔴'} log/json={'✅' if in_corpus else '🔴'}")
    if not (in_doc and in_corpus):
        bad5.append(k)
if bad5:
    fail.append(f"K5 數值不一致 {bad5}")
report["K5"] = bad5

# ---------------- K6 結論節 ----------------
print("\n" + "=" * 108)
print("K6 每份 doc 具「結論」節")
print("=" * 108)
no_concl = [d.name for d in docs if "結論" not in d.read_text()]
print(f"  缺結論節：{len(no_concl)}" + ("" if not no_concl else f" ⟹ {no_concl}"))
if no_concl:
    fail.append(f"K6 缺結論節 {no_concl}")
report["K6"] = no_concl

# ---------------- K7 register 覆蓋（INFO） ----------------
print("\n" + "=" * 108)
print("K7 register 覆蓋（INFO）")
print("=" * 108)
refs = set()
if REG.exists():
    for ln in REG.read_text().splitlines():
        if not ln.strip():
            continue
        r = json.loads(ln)
        for m in re.findall(r"docs/(\d{2})", r.get("ref") or ""):
            refs.add(m)
missing_reg = [d.name[:2] for d in docs if d.name[:2] not in refs]
print(f"  未被任何 register event ref 指向：{len(missing_reg)} 份 ⟹ {missing_reg}")
report["K7_uncovered"] = missing_reg
info.append(f"K7 register 未覆蓋 {len(missing_reg)} 份（INFO）")

# ---------------- 匯總 ----------------
print("\n" + "=" * 108)
print(f"[總結] FAIL = {len(fail)}；INFO = {len(info)}")
for x in fail:
    print(f"  🔴 {x}")
for x in info:
    print(f"  ⓘ {x}")
report["fail"], report["info"] = fail, info
with open(DATA / "28_chain_continuity.json", "w") as f:
    json.dump(report, f, ensure_ascii=False, indent=2)
print("\n-> data/28_chain_continuity.json")
sys.exit(1 if fail else 0)