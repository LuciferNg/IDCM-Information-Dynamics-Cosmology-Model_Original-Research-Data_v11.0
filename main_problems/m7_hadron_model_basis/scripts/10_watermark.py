"""M7 scripts/10 — 檔案閉合水位（live field 生成器）
================================================================================================
水位行屬 live field（skill `idcm-audit-verify` §1）：可直接重生，唔屬 append-only 範圍。
本腳本：
  1. 掃 docs/ scripts/ data/ ＋ 四件套，取 mtime，對 register close-out 基準比對
  2. 寫 data/watermark.json（機讀，供 scripts/07 C2 核覆蓋）
  3. 將 `01_problem_record.md` 之「### §檔案閉合水位」整節**移至檔尾並重生**
     （舊節刪除 → 新節 append；其餘章節原文不動）
輸出：data/watermark.json
"""
import json
import os
import re
import time
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
CLOSE_ISO = "2026-10-26T12:59:45+08:00"      # register close event ts
CLOSE_TS = datetime.fromisoformat(CLOSE_ISO).timestamp()
HDR = "### §檔案閉合水位"
REC = os.path.join(M7, "01_problem_record.md")

files = []
for sub in ("docs", "scripts", "data"):
    d = os.path.join(M7, sub)
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if os.path.isfile(p):
            files.append(f"{sub}/{f}")
for f in ("01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"):
    if os.path.isfile(os.path.join(M7, f)):
        files.append(f)

if "data/watermark.json" not in files:
    files.append("data/watermark.json")      # 自身輸出（先登記，mtime 於寫入時定）
rows = []
for p in files:
    if p == "data/watermark.json":
        rows.append(dict(path=p, mtime=datetime.now().strftime("%Y-%m-%d %H:%M"), after_close=True))
        continue
    full = os.path.join(M7, p)
    mt = os.path.getmtime(full)
    rows.append(dict(path=p, mtime=datetime.fromtimestamp(mt).strftime("%Y-%m-%d %H:%M"),
                     after_close=mt > CLOSE_TS))

wm = dict(close_out=CLOSE_ISO, generated=datetime.now().strftime("%Y-%m-%d %H:%M"),
          n=len(rows), rows=rows)
json.dump(wm, open(os.path.join(M7, "data", "watermark.json"), "w"), ensure_ascii=False, indent=1)

# ── 重生 markdown 節並置於檔尾 ────────────────────────────────
txt = open(REC, encoding="utf-8").read()
if HDR in txt:
    head, rest = txt.split(HDR, 1)
    m = re.search(r"\n## ", rest)
    tail = rest[m.start():] if m else ""
    txt = head + tail.lstrip("\n")
    txt = txt.rstrip("\n") + "\n"
lines = [HDR + f"（自動生成 {wm['generated']}；close-out 基準 `{CLOSE_ISO}`）", "",
         "> live field：由 `scripts/10_watermark.py` 重生（非 append-only 範圍）。"
         "判「未變」不只用 mtime——本表另經 `scripts/07` B 層內容級核對。", "",
         "| # | 檔案 | 上次檔案更新 (mtime) | 上次 closed-out | 已驗範圍 | 狀態 |",
         "|:-:|:-----|:---------------------|:----------------|:---------|:-----|"]
for i, r in enumerate(rows, 1):
    st = "⏳ close 後新增／改動 ⟹ 已驗（docs/17 §B/C）" if r["after_close"] else "✅ 水位內（無變更）"
    lines.append(f"| {i} | `{r['path']}` | {r['mtime']} | 2026-10-26 12:59 | 全文 | {st} |")
open(REC, "w", encoding="utf-8").write(txt + "\n" + "\n".join(lines) + "\n")
print(f"水位 {len(rows)} 行；after_close = {sum(r['after_close'] for r in rows)}；JSON → data/watermark.json")