"""M7 R1 配套 — 概念詞語料證據（grep 0 命中 ⟹ 定義自建；凍結檔 hash 前置）
================================================================================================
判準（寫於執行前）：
  E1 逐詞統計 v9.0／v10.0／v11.0（**排除本 M 目錄**）之命中數；
     「命中」= 檔案含該字串（files_with），另報總行數（lines）。
  E2 任何詞若只在 M4 docs/32／M7 自身出現 ⟹ 判「概念未形式化」（= 本 M 自建之依據）。
  E3 定義凍結：對 docs/11 之 §二 定義區塊原文算 SHA-256；之後任何修改 = 新 doc 明寫「移動球門」。
  E0 語料定義（v2，2026-10-26 修）：**只掃各版本 git HEAD 之 tracked files**（`git ls-files`），
     並排除本 M 目錄。理由：untracked 之並行 WIP（如 m6_*）會令計數漂移而不可重現；
     tracked-only ⟹ 數字綁定 commit hash，可精確復現。指紋寫入 JSON。

輸出：data/03_grep_evidence.json
"""
import hashlib
import json
import os
import re
import subprocess

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
M7 = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
VERS = ["IDCM-v9.0--W-Field-Lattice", "IDCM-v10.0--It-From-Beat", "IDCM-v11.0--Reality-Architecture"]

TERMS = [
    ("w_k 佔用", r"(?<![A-Za-z0-9_])w_\{?k\}?|w_k\(h\)"),
    ("佔有率／佔據數", r"佔有率|佔據數|occupation number|occupancy"),
    ("多夸克態", r"多夸克|multiquark|multi-quark"),
    ("Fock 展開", r"Fock"),
    ("投影子差", r"投影子差|projector difference|P_\{?k\}?\s*[-−]\s*P"),
    ("isovector", r"isovector|同位旋向量"),
    ("ūu − d̄d", r"ūu|ubar u|ūu|d̄d|dbar d"),
    ("夸克含量", r"夸克含量|quark content|q̄q|qqbar"),
    ("海夸克／介子雲", r"海夸克|sea quark|介子雲|meson cloud|pion cloud"),
    ("質量泛函非可加", r"非可加|non-additiv|質量泛函|mass functional"),
    ("⟨O⟩ 矩陣元", r"⟨O⟩|matrix element 矩陣元"),
]


def corpus(ver):
    """回傳 (files_rel, head, mode)：優先用 git tracked（可重現），否則全掃（fallback）。"""
    base = os.path.join(ROOT, ver)
    if not os.path.isdir(base):
        return [], None, "missing"
    try:
        head = subprocess.run(["git", "-C", base, "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        tracked = subprocess.run(["git", "-C", base, "ls-files"],
                                 capture_output=True, text=True, check=True).stdout.splitlines()
        rel = os.path.relpath(M7, base)
        files = [f for f in tracked
                 if f.endswith((".md", ".py", ".json", ".tex")) and not f.startswith(rel)]
        return files, head, "git-tracked"
    except (subprocess.CalledProcessError, OSError):
        out = []
        for dirpath, dirnames, filenames in os.walk(base):
            if M7 in os.path.abspath(dirpath) or "/.git" in dirpath:
                continue
            for fn in filenames:
                if fn.endswith((".md", ".py", ".json", ".tex")):
                    out.append(os.path.relpath(os.path.join(dirpath, fn), base))
        return sorted(out), None, "walk-fallback"


CORPUS = {}
for _v in VERS:
    CORPUS[_v] = corpus(_v)


def grep(term_re, ver):
    base = os.path.join(ROOT, ver)
    out = {"files": [], "lines": 0}
    if not os.path.isdir(base):
        return out
    for rel in CORPUS[ver][0]:
        fp = os.path.join(base, rel)
        if True:
            try:
                txt = open(fp, encoding="utf-8", errors="ignore").read()
            except OSError:
                continue
            hits = re.findall(term_re, txt)
            if hits:
                out["files"].append(rel)
                out["lines"] += len(hits)
    return out


def sha256_block(path, start_marker, end_marker):
    txt = open(path, encoding="utf-8").read()
    i = txt.find(start_marker)
    j = txt.find(end_marker, i + 1)
    if i < 0 or j < 0:
        return None, None
    block = txt[i:j]
    return hashlib.sha256(block.encode("utf-8")).hexdigest(), block


def main():
    res = {}
    print("=" * 100)
    print("E0 語料指紋（git HEAD tracked files，排除本 M）")
    for v in VERS:
        hh, nn, mm = CORPUS[v][1], len(CORPUS[v][0]), CORPUS[v][2]
        print(f"  {v:<38} HEAD={hh}  files={nn}  mode={mm}")
    print("\nE1 概念詞語料命中（排除本 M 目錄）")
    print(f"{'詞':22}" + "".join(f"{v.split('--')[0]:>18}" for v in VERS))
    for name, rx in TERMS:
        row = {}
        cells = []
        for v in VERS:
            r = grep(rx, v)
            row[v] = dict(nfiles=len(r["files"]), nhits=r["lines"], files=r["files"][:6])
            cells.append(f"{len(r['files'])}f/{r['lines']}h")
        res[name] = row
        print(f"  {name:20}" + "".join(f"{c:>18}" for c in cells))

    # E2 逐詞判定：是否僅見於 M4 引用鏈（docs/32 等）而非形式化定義
    print("\nE2 判定（v9.0 為原文權威；v9 命中 = 已形式化，v9 0 命中 = 概念未形式化）")
    verdicts = {}
    for name, _ in TERMS:
        n9 = res[name][VERS[0]]["nfiles"]
        verdicts[name] = "v9 已見（須逐條核是否同義）" if n9 else "v9 0 命中 ⟹ 未形式化"
        print(f"  {name:20} v9 = {n9} ⟹ {verdicts[name]}")

    # E3 定義凍結 hash
    print("\nE3 定義凍結 hash（docs/11 §二 定義區塊）")
    doc = os.path.join(M7, "docs", "11_definition_freeze.md")
    h, block = sha256_block(doc, "<!-- DEF-FREEZE-BEGIN -->", "<!-- DEF-FREEZE-END -->")
    if h:
        print(f"  SHA-256 = {h}")
        print(f"  區塊長度 = {len(block)} chars")
    else:
        print("  🔴 未找到凍結標記（docs/11 未寫或標記缺失）")

    out = dict(terms=res, verdicts=verdicts, def_sha256=h,
               exclusion="m7_hadron_model_basis 已排除",
               corpus_scope="git HEAD tracked files only（untracked WIP 不計）",
               fingerprint={v: dict(head=CORPUS[v][1], nfiles=len(CORPUS[v][0]), mode=CORPUS[v][2])
                            for v in VERS},
               versions=VERS)
    p = os.path.normpath(os.path.join(M7, "data", "03_grep_evidence.json"))
    with open(p, "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(f"\nJSON → {p}")


if __name__ == "__main__":
    main()