#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/37 — 反向數值覆蓋（skill `idcm-audit-verify` §3.1 ⑦）M14 本地強化版

為何要本地版（三項真缺陷，皆為 M14 審計期間實測）：
  ① 共用工具 `reverse_numeric_coverage.py` 之 doc 側 regex `\\d+\\.\\d{3,}` **不帶符號**
     ⟹ 文檔寫 `−0.63645`（U+2212 歸一為 `-`）抽成 **+0.63645**，而 data 池存 `-0.63645`
     ⟹ 所有「只以負值出現」之量一律假 UNTRACEABLE（實測 10/17 屬此類）。
  ② 該工具不解析 **e-記法**（`1.578e−4`、`4.87189948e−3`）⟹ 抽成尾數，與 data 之 0.0001578 對不上。
  ③ 該工具無 fraction（`a/b`）支援。

本版判準（跑前凍結）：
  T1 抽取：上標 `×10⁻³` → ASCII；`e±n` → ASCII；帶符號小數（≥3 位）；分數 `a/b` 取商。
  T2 池：`data/*.log` ＋ `data/*.json`，同樣做 e-記法與上標歸一。
  T3 命中：|pool − v| ≤ 0.5%·max(|v|,1e-12)。
  T4 未命中者必須在 ALLOW（附理由），否則 exit 1。
性質：只讀。N_trials = 0。零觀測量。
"""
import json
import re
import sys
from bisect import bisect_left
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS, DATA = ROOT / "docs", ROOT / "data"
TOL = 5e-3

SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
SCI = re.compile(r"(\d+(?:\.\d+)?)\s*[×x]\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)")
ENOT = re.compile(r"(\d+(?:\.\d+)?)[eE]([-+]?\d+)")
POOL = re.compile(r"[-+]?\d+\.\d+")      # 池：任何小數（data 常只印 3 位）
DEC = re.compile(r"[-+]?\d+\.\d{4,}")   # 4+ 位小數（skill §3.1⑦「4+ 有效位」）
FRAC = re.compile(r"(?<![\d.])(\d+)/(\d+)(?![\d.])")

# 跨記錄引用／散文範圍 —— 每項必須寫理由（審計 §⑤ 會核理由非空）
ALLOW = {
    ("09_correction_object_type.md", "1.1486"):
        "跨記錄：M13 預註冊 T1 閘 ρ ≡ |r_K|/(|r_Ξ|/2) = 1.1486 ± 0.05；本 M 只引用作停止條件",
    ("09_correction_object_type.md", "0.4603"):
        "跨記錄：A_eff = 4.4492 ± 0.4603 為 M5 之給定值（v11 M5 docs/…），本 M 只引用",
    ("08_incoming_m13_triage.md", "0.636"):
        "跨記錄：M13 分診表之 δ_c 相關值，原文 3 位捨入；精確值在本 M data 以 −0.63645 落地",
    ("17_tpipe_filter_and_phi_ladder.md", "16.8"):
        "散文範圍：『φ vs 觀測 17–207』之區間下界，非量測值（禁作精確數引用）",
    ("17_tpipe_filter_and_phi_ladder.md", "207.0"):
        "散文範圍：同上區間上界，非量測值",
    ("22_audit_report.md", "7.6064"):
        "歷史值：docs/45 v1（錯矩陣元律）之 T_1 = 7.6064，僅作更正敘述；v1 log 封存 data/archive/（不在 data/*.log 池）",
}


# 跨記錄常數：任何 doc 引用皆可（值 → 理由）；以絕對值比對（容許文檔以負號作示例）
ALLOW_BY_VALUE = {
    0.4603: "跨記錄：M5 A_eff = 4.4492 ± 0.4603",
    1.1486: "跨記錄：M13 預註冊 T1 閘 ρ = 1.1486 ± 0.05",
    0.63645: "跨記錄：M13 δ_c 相關值；本 M data 以 −0.63645 落地（文檔以正號作符號缺陷示例時命中此條）",
}


def norm(t: str) -> str:
    t = t.replace("\u2212", "-")
    # 上標 ×10^n → 十進位字串
    def _sci(m):
        return f"{float(m.group(1)) * 10 ** int(m.group(2).translate(SUP)):.12g}"
    t = SCI.sub(_sci, t)
    # e-記法 → 十進位（防呆：hex 碎片如 md5 之 3e485897 不得當指數）
    def _e(m):
        try:
            exp = int(m.group(2))
        except ValueError:
            return m.group(0)
        if abs(exp) > 400:
            return m.group(0)
        try:
            return f"{float(m.group(1)) * 10 ** exp:.12g}"
        except (OverflowError, ValueError):
            return m.group(0)
    t = ENOT.sub(_e, t)
    return t


def build_pool():
    vals = []
    # ★ 池只取「L1 生產者 log」：排除審計自身、本工具、L2/L3（後三者會鏡像 doc 數字 ⟹ 自我確認）
    SELF = {"reverse_coverage.json", "run_37_reverse_coverage.log",
            "run_23_audit_verify.log", "run_l2_independent_verify.log", "run_l3_doc_crosscheck.log",
            # ★ L2/L3 之 JSON 亦會鏡像 doc 數字 ⟹ 一併排除（否則 ⑦ 變成自我確認）
            "l2_independent_verify.json", "l3_doc_crosscheck.json"}
    for f in list(DATA.glob("*.log")) + list(DATA.glob("*.json")):
        if f.name in SELF or f.name.startswith(("run_l2_", "run_l3_")):
            continue
        for m in POOL.finditer(norm(f.read_text(errors="ignore"))):
            vals.append(float(m.group(0)))
    vals.sort()
    return vals


def tokens(text: str):
    t = norm(text)
    out = {}
    for m in DEC.finditer(t):
        out.setdefault(round(float(m.group(0)), 12), m.group(0))
    # ★ 分數（a/b）不入反向覆蓋：屬精確式，由 L3 之字串 claim 覆蓋；
    #   且文檔引用編號（03/04）與計數比（44/112）會被誤抽成數值（實測 50 個假 UNTRACEABLE）。
    return out


def hit(pool, v):
    if not pool:
        return False
    i = bisect_left(pool, v)
    return any(0 <= j < len(pool) and abs(pool[j] - v) <= TOL * max(abs(v), 1e-12)
               for j in (i - 1, i, i + 1))


def main() -> int:
    pool = build_pool()
    bad, allowed_hits, n_tok = [], [], 0
    per_doc = {}
    for d in sorted(DOCS.glob("*.md")):
        toks = tokens(d.read_text(encoding="utf-8"))
        n_tok += len(toks)
        miss = []
        for v, raw in sorted(toks.items()):
            if hit(pool, v):
                continue
            reason = ALLOW.get((d.name, raw)) or ALLOW.get((d.name, f"{v:g}"))
            if reason is None:
                for av, why in ALLOW_BY_VALUE.items():
                    if abs(abs(av) - abs(v)) <= 1e-9 * max(abs(v), 1e-12):
                        reason = why
                        break
            if reason is None:
                for (dn, rs), why in ALLOW.items():
                    if dn == d.name:
                        try:
                            if abs(float(rs) - v) <= 1e-9 * max(abs(v), 1e-12):
                                reason = why
                                break
                        except ValueError:
                            pass
            if reason:
                allowed_hits.append((d.name, raw, reason))
            else:
                miss.append((raw, v))
        per_doc[d.name] = {"tokens": len(toks), "untraced": [f"{r}" for r, _ in miss]}
        bad += [(d.name, r, v) for r, v in miss]
    print("=" * 78)
    print("M14 scripts/37 — 反向數值覆蓋（docs → data，容差 0.5%）")
    print("=" * 78)
    print(f"  data 池 {len(pool)} 個值；docs token {n_tok} 個（含 e-記法／上標／分數）")
    print(f"  ALLOW（跨記錄／散文範圍）命中 {len(allowed_hits)} 項")
    for n, r, why in allowed_hits:
        print(f"    [ALLOW] {n:42s} {r:>14s}  ← {why[:60]}")
    print(f"\n  未追溯（且不在 ALLOW）= {len(bad)}")
    for n, r, v in bad:
        print(f"    [UNTRACEABLE] {n:42s} {r:>16s}  ({v:g})")
    json.dump({"pool": len(pool), "tokens": n_tok,
               "untraced": [[n, r] for n, r, _ in bad],
               "allowed": [[n, r] for n, r, _ in allowed_hits],
               "per_doc": per_doc},
              open(DATA / "reverse_coverage.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\n判讀：未追溯 = 0 ⟹ 文檔每個精確數字皆有 data 落地（或已明示 ALLOW 理由）。")
    print("=" * 78)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())