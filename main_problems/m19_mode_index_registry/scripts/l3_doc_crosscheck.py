#!/usr/bin/env python3
"""
M19 / l3 — 文檔逐 claim 交叉核（每條 doc 數字必須喺機器產物內出現）

判據（跑前凍結）:
  L3.1  本清單每條 claim 之 token 必須喺指定產物（log／JSON）內逐字出現
        （比較前只做 Unicode 正規化：− → -、空白壓縮；唔做任何數值寬容）
  L3.2  覆蓋 docs/01–06 全部核心數字；漏一條 = FAIL（唔可以「大意相同」）
  L3.3  docs/ 內不得出現「未經產物支持」之新數字 —— 以本清單為完備集抽查
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def norm(s):
    return " ".join(s.replace("\u2212", "-").split())


CLAIMS = [
    ("docs/01", "RESULT 11/11 PASS", "run_00_reproduction_gate.log"),
    ("docs/01", "0.005095770779492917275388486", "run_00_reproduction_gate.log"),
    ("docs/01", "|dev|=2.75388e-19", "run_00_reproduction_gate.log"),
    ("docs/01", "γ_44 − 4 = 0.0", "run_00_reproduction_gate.log"),
    ("docs/01", "max|dev|=6.22302e-61", "run_00_reproduction_gate.log"),
    ("docs/01", "min Δγ=0.005095770779", "run_00_reproduction_gate.log"),
    ("docs/02", "RESULT 8/8 PASS", "run_01_frame_census.log"),
    ("docs/02", "max rel resid = 0.0", "run_01_frame_census.log"),
    ("docs/02", "集外 5/9", "run_01_frame_census.log"),
    ("docs/02", "1214615.612988", "run_01_frame_census.log"),
    ("docs/03", "17,280 = 9!/21", "run_02_bit_decomposition.log"),
    ("docs/03", "結構 4.392317 bit", "run_02_bit_decomposition.log"),
    ("docs/03", "輸入 14.076816 bit", "run_02_bit_decomposition.log"),
    ("docs/03", "18.469133", "run_02_bit_decomposition.log"),
    ("docs/03", "4:2,880", "run_03_admissibility_and_readings.log"),
    ("docs/03", "5:720", "run_03_admissibility_and_readings.log"),
    ("docs/03", "7:2,880", "run_03_admissibility_and_readings.log"),
    ("docs/03", "4320", "run_03_admissibility_and_readings.log"),
    ("docs/03", "2880", "run_03_admissibility_and_readings.log"),
    ("docs/04", "RESULT 5/5 PASS", "run_04_cbt_scan_density.log"),
    ("docs/04", "p = 1.983e-03", "run_04_cbt_scan_density.log"),
    ("docs/04", "p = 1.189e-02", "run_04_cbt_scan_density.log"),
    ("docs/04", "理論 1.984e-03", "run_04_cbt_scan_density.log"),
    ("docs/04", "理論 1.190e-02", "run_04_cbt_scan_density.log"),
    ("docs/04", "T1 0.0%", "run_04_cbt_scan_density.log"),
    ("docs/05", "3 結構類別 ＋ 6 輸入", "run_05_registry_build.log"),
    ("docs/05", "非法：裸 k_q → BLOCK", "run_05_registry_build.log"),
    ("docs/05", "GX1", "run_05_registry_build.log"),
    ("docs/06", "全 T3（χ ≡ n mod 4 全施加）可行雙射 = 0", "run_07_slot_scope_theorem.log"),
    ("docs/06", "6 粒子 vs 4 槽", "run_07_slot_scope_theorem.log"),
]

print("=" * 92)
print("M19 l3 — 文檔逐 claim 交叉核")
print("=" * 92)

missing = []
for doc, token, art in CLAIMS:
    path = DATA / art
    txt = norm(path.read_text()) if path.exists() else ""
    if norm(token) not in txt:
        missing.append((doc, token, art))
chk("L3.1 全部 claim token 喺產物內逐字出現", not missing, f"{len(CLAIMS)} 條；缺 {missing}")

docs = sorted((ROOT / "docs").glob("*.md"))
covered = {d for d, _, _ in CLAIMS}
chk("L3.2 覆蓋 docs/01–06",
    {f"docs/{p.name[:2]}" for p in docs} >= covered and len(covered) >= 6,
    f"{len(docs)} 份 doc；claim 覆蓋 {sorted(covered)}")

# L3.3：docs 內所有「多位小數／科學記數」必須有產物支持（數值容差 5e-3，容許四捨五入展示）
import re


def floats_of(text):
    out = []
    for m in re.finditer(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", text.replace("\u2212", "-")):
        try:
            out.append(float(m.group(0)))
        except ValueError:
            pass
    return out


art_floats = []
for p in list(DATA.glob("*.log")) + list(DATA.glob("0*.json")):
    art_floats.extend(floats_of(p.read_text()))

unsupported = []
for p in docs:
    body = p.read_text()
    for m in re.finditer(r"\d+\.\d{3,}(?:[eE][-\u2212]?\d+)?|\d\.\d+[eE][-\u2212]\d+", body):
        tok = m.group(0)
        # 略過 § 章節號（如 §7.1）
        start = m.start()
        if start >= 1 and body[start - 1] in "§":
            continue
        try:
            v = float(tok.replace("\u2212", "-"))
        except ValueError:
            continue
        if not any(abs(f - v) <= 5e-3 * max(abs(v), 1e-12) for f in art_floats):
            unsupported.append((p.name, tok))
chk("L3.3 文檔內多位小數／科學記數全部有產物支持（容差 5e-3）",
    not unsupported, f"未支持 {len(unsupported)}：{unsupported[:8]}")

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
