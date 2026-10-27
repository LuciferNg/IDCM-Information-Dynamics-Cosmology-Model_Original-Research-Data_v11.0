#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — 文檔交叉核（通用模板；idcm-audit-verify §2.4；所有 M/O 記錄必跑）
================================================================================
用法
  1. 先跑 L2（`scripts/l2_independent_verify.py`）產出 `data/l2_independent_verify.json`
  2. `cp` 到 `<記錄>/scripts/l3_doc_crosscheck.py`，只改 `CLAIMS` 區塊
  3. `/usr/bin/python3 scripts/l3_doc_crosscheck.py`
     → `data/l3_doc_crosscheck.json` ＋ `data/run_l3.log`；FAIL>0 時 exit 1

鐵律
  R1 **唯讀**：本檔對任何 doc／script 只 open 讀，零寫入（只寫 data/ 下兩個輸出）。
  R2 期望值**只可來自 L2**（`_derived[key]`）；禁抄 docs 內數字（否則等於自我確認）。
  R3 比對用**數值容差依文檔有效位**：≥6 位 → 1e−4；≥4 位 → 5e−4；其餘 → 5e−3。
     禁「字串相等」式比對 —— 文檔必然用不同精度／記法（純小數 vs e-記法）。
  R4 **負值以絕對值比**：token 抽取唔保留正負號（文檔用 U+2212），若不取 abs，
     任何負期望值（dex、Δ、−log…）都會永遠 FAIL（M8 實錄）。
  R5 容差要**仍分辨得出真錯**：M8 實測 78.55 vs 79.5538 = 1.26% ⟹ 必須 FAIL；
     5.223686e35 vs 2.3e35 = 56% ⟹ 必須 FAIL。若某個 claim 用預設容差捉唔到，
     改 `tol=` 顯式收緊，**唔准**放寬吞數。
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG = []
RES = {}
NFAIL = 0


# ============================================================================
# ① CLAIMS —— (doc, L2 key, 模式, 標籤)   ⚠️ 示例內容：換成本記錄嘅宣稱
#    模式："num"（數值容差，預設）｜"str"（L2 之 str 必須逐字出現）
#    可用 "tol": 顯式相對容差覆寫（更緊，用嚟捉二階換算錯）
# ============================================================================
CLAIMS = [
    {"doc": "docs/36_k_anomaly_significance.md", "key": "z_K", "label": "K z = 2.627σ"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "z_B", "label": "B z = 0.0279σ"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "c_model_K", "label": "c_model(K) = −0.229293"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "c_K", "label": "c_data(K) = −0.0135"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "sigma_c_K", "label": "σ_c(K) = 0.08216"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "delta_c_K", "label": "|Δc| = 0.215793"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "gain_K", "label": "精度改善 1.904×"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "dist_pct_K", "label": "鑑別力 279.1%"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "u", "label": "u = 18.4089393"},
    {"doc": "docs/36_k_anomaly_significance.md", "key": "v", "label": "v = 0.0634782072"},
    {"doc": "docs/31_scale_identifiability.md", "key": "u", "label": "docs/31 U = 18.408939"},
    {"doc": "docs/31_scale_identifiability.md", "key": "v", "label": "docs/31 V = 0.0634782"},
]

TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def tokens(text):
    """抽全部數值 token → [(float, 有效位數)]（U+2212 先歸一）"""
    out = []
    for m in TOK.finditer(text.replace("\u2212", "-")):
        mant, exp = m.group(1), m.group(2)
        try:
            v = float(mant) * (10.0 ** int(exp)) if exp else float(mant)
        except (ValueError, OverflowError):
            continue
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        out.append((v, sig))
    return out


def match(toks, target, tol_override=None):
    """回傳 (命中?, 最近值, 相對差)；負期望值以絕對值比（R4）"""
    t = abs(target)
    best, best_d = None, float("inf")
    for v, sig in toks:
        if v == 0:
            continue
        rel = abs(abs(v) - t) / t
        tol = tol_override if tol_override else (1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3))
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


def main():
    hdr("§0 讀 L2（期望值唯一來源）")
    if not L2PATH.exists():
        print(f"✗ 缺 {L2PATH} —— 先跑 scripts/l2_independent_verify.py")
        return 2
    l2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    der = l2.get("_derived", {})
    chk("E1 L2 有 _derived 且非空", bool(der), f"keys = {sorted(der)[:8]}...")

    hdr("§1 逐 claim 交叉核（唯讀）")
    cache, by_doc = {}, {}
    for c in CLAIMS:
        doc, key, label = c["doc"], c["key"], c.get("label", c["key"])
        mode = c.get("mode", "num")
        p = ROOT / doc
        if not p.exists():
            chk(f"C.{doc} 存在", False, "檔案不存在")
            continue
        if doc not in cache:
            cache[doc] = p.read_text(encoding="utf-8")
        if key not in der:
            chk(f"C.{doc}.{label}（L2 key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1

        if mode == "str":
            want = der[key]["str"]
            ok = want in cache[doc]
            if ok:
                by_doc[doc][0] += 1
            chk(f"C.{doc.split('/')[-1][:16]}.{label}（字串）", ok,
                f"L2 = {want!r}" if ok else f"✗ 文檔冇 {want!r}")
            continue

        exp = der[key].get("float")
        if exp is None:
            chk(f"C.{doc.split('/')[-1][:16]}.{label}", False, "L2 值超出 float 範圍（改用 str 模式）")
            continue
        ok, got, rel = match(tokens(cache[doc]), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{doc.split('/')[-1][:16]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok
            else f"✗ L2 {exp:.6g}，文檔最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率")
    for d, (h, t) in by_doc.items():
        print(f"    {d:46s} {h}/{t}")
        chk(f"F.{d} 全命中（{h}/{t}）", h == t)

    hdr("§3 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "by_doc": by_doc}

    (DATA / "l3_doc_crosscheck.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / 'l3_doc_crosscheck.json'}\n寫出: {DATA / 'run_l3.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
